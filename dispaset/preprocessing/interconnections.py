"""
Functions for handling interconnections between zones in Dispa-SET.

This module contains all interconnection-related functions that were previously in utils.py.
It handles the creation, formatting, and processing of interconnections between zones.
"""

import sys
import logging
import numpy as np
import pandas as pd


def interconnections(simulation_list, ntc_inter, historical_flows):
    """
    Function that checks for the possible interconnections of the zones included
    in the simulation. If the interconnections occurs between two of the zones
    defined by the user to perform the simulation with, it extracts the NTC between
    those two zones. If the interconnection occurs between one of the zones
    selected by the user and one country outside the simulation, it extracts the
    physical flows; it does so for each pair (country inside-country outside) and
    sums them together creating the interconnection of this country with the RoW.

    :param simulation_list:     List of simulated zones
    :param ntc_inter:           Day-ahead net transfer capacities (pd dataframe)
    :param historical_flows:    Historical flows (pd dataframe)
    :return:                    Tuple containing (internal interconnections, external interconnections, all interconnections)
    """
    index = ntc_inter.index.tz_localize(None).intersection(historical_flows.index.tz_localize(None))
    if len(index) == 0:
        logging.error('The two input dataframes (NTCs and Historical flows) must have the same index. '
                      'No common values have been found')
        sys.exit(1)
    elif len(index) < len(ntc_inter) or len(index) < len(historical_flows):
        diff = np.maximum(len(historical_flows), len(ntc_inter)) - len(index)
        logging.warning('The two input dataframes (NTCs and Historical flows) do not share the same index, '
                        'although some values are common. The intersection has been considered and ' + str(diff) +
                        ' data points have been lost')
    # Checking that all values are positive:
    if (ntc_inter.values < 0).any():
        pos = np.where(ntc_inter.values < 0)
        logging.warning('At least one NTC value is negative, for example in line ' + str(ntc_inter.columns[pos[1][0]]) +
                        ' and time step ' + str(ntc_inter.index[pos[0][0]]))
    all_connections = []
    simulation_connections = []
    # List all connections from the dataframe headers:
    con_list = historical_flows.columns.tolist() + [x for x in ntc_inter.columns.tolist() if
                                                 x not in historical_flows.columns.tolist()]
    for connection in con_list:
        # Ensure connection has the correct format with "->" separator
        if ' -> ' in connection:
            z = connection.split(' -> ')
        else:
            z = connection.split('->')
        if len(z) == 2:
            connection = z[0].strip() + ' -> ' + z[1].strip()
            z = [z[0].strip(), z[1].strip()]  # Update z with stripped values
        
        # Check if this connection is between simulated zones
        if z[0] in simulation_list and z[1] in simulation_list:
            # Internal connection between simulated zones
            all_connections.append(connection)
            simulation_connections.append(connection)
        elif z[0] in simulation_list or z[1] in simulation_list:
            # External connection with at least one simulated zone
            all_connections.append(connection)

    df_zones_simulated = pd.DataFrame(index=index)
    for interconnection in simulation_connections:
        # Handle both formats when looking up in ntc_inter
        if interconnection in ntc_inter.columns:
            df_zones_simulated[interconnection] = ntc_inter[interconnection]
            logging.info('Detected interconnection ' + interconnection +
                         '. The historical NTCs will be imposed as maximum flow value if the NTC method is used')
        else:
            # Try without spaces around arrow
            no_spaces = interconnection.replace(' -> ', '->')
            if no_spaces in ntc_inter.columns:
                df_zones_simulated[interconnection] = ntc_inter[no_spaces]
                logging.info('Detected interconnection ' + interconnection +
                            '. The historical NTCs will be imposed as maximum flow value if the NTC method is used')
    interconnections1 = df_zones_simulated.columns

    # Display a warning if a zone is isolated:
    for z in simulation_list:
        if not any([z in conn for conn in interconnections1]) and len(simulation_list) > 1:
            logging.warning('Zone ' + z + ' does not appear to be connected to any other zone in the NTC table. '
                            'It should be simulated in isolation')

    df_row_temp = pd.DataFrame(index=index)
    conn_names = []
    for interconnection in all_connections:
        # Handle both formats when looking up in historical_flows
        if interconnection in historical_flows.columns:
            if interconnection not in simulation_connections:
                df_row_temp[interconnection] = historical_flows[interconnection]
                conn_names.append(interconnection)
        else:
            # Try without spaces around arrow
            no_spaces = interconnection.replace(' -> ', '->')
            if no_spaces in historical_flows.columns and no_spaces not in simulation_connections:
                df_row_temp[interconnection] = historical_flows[no_spaces]
                conn_names.append(interconnection)

    # Process flows between simulated zones and RoW
    df_zones_row = pd.DataFrame(index=index)
    
    # First collect all the simulated zones involved in external flows
    simulated_zones_in_external_flows = set()
    for connection in conn_names:
        if ' -> ' in connection:
            z = connection.split(' -> ')
        else:
            z = connection.split('->')
        
        zone_from = z[0].strip()
        zone_to = z[1].strip()
        
        if zone_from in simulation_list:
            simulated_zones_in_external_flows.add(zone_from)
        if zone_to in simulation_list:
            simulated_zones_in_external_flows.add(zone_to)
    
    # Process each simulated zone with external connections
    # In the new approach, all flows for a zone with RoW will be merged into a single entry
    # where positive flow means export from zone to RoW, negative flow means import to zone from RoW
    for zone in simulated_zones_in_external_flows:
        net_flows = pd.Series(0, index=index)  # Initialize a series of zeros with the right index
        
        for connection in conn_names:
            if ' -> ' in connection:
                z = connection.split(' -> ')
            else:
                z = connection.split('->')
            
            zone_from = z[0].strip()
            zone_to = z[1].strip()
            
            # Zone -> External: Add as positive flow (export)
            if zone_from == zone and zone_to not in simulation_list:
                net_flows += df_row_temp[connection]
                logging.info(f'Detected interconnection {connection}, happening between simulated zone {zone} '
                             f'and the rest of the world. Adding as positive flow to {zone} -> RoW')
            
            # External -> Zone: Add as negative flow (import)
            elif zone_to == zone and zone_from not in simulation_list:
                net_flows -= df_row_temp[connection]  # Subtract for imports
                logging.info(f'Detected interconnection {connection}, happening between the rest of the world '
                             f'and simulated zone {zone}. Adding as negative flow to {zone} -> RoW')
        
        # Store the net flow for this zone
        if not net_flows.equals(pd.Series(0, index=index)):  # Only add if there are non-zero flows
            flow_name = f'{zone} -> RoW'
            df_zones_row[flow_name] = net_flows
            logging.info(f'Created merged bidirectional flow {flow_name} (positive=export, negative=import)')
            
    interconnections2 = df_zones_row.columns
    inter = list(interconnections1) + list(interconnections2)
    return df_zones_simulated, df_zones_row, inter


def merge_lines(interconnections_sim, interconnections_row, price_transmission_raw):
    """
    Merges bidirectional interconnections (e.g., 'A -> B' and 'B -> A') into a single representation 'A -> B'.

    - For NTCs (`interconnections_sim`), it averages the values. If the average difference between
      the two directions exceeds 20% at any point, a critical warning is logged.
    - For fixed flows (`interconnections_row`), it calculates the net flow: Flow(A->B) = Flow(A->B)_original - Flow(B->A)_original.
    - For prices (`price_transmission_raw`), it averages the values. If the average difference between
      the two directions exceeds 20% at any point, a critical warning is logged.

    :param interconnections_sim: DataFrame with NTCs for simulated interconnections.
    :param interconnections_row: DataFrame with historical/fixed flows for RoW interconnections.
    :param price_transmission_raw: DataFrame with raw transmission prices.
    :return: Tuple containing the modified (interconnections_sim, interconnections_row, price_transmission_raw).
    """
    # Work on copies to avoid modifying original inputs unexpectedly
    sim_df = interconnections_sim.copy()
    row_df = interconnections_row.copy()
    price_df = price_transmission_raw.copy()

    # Get the list of all unique interconnection lines
    all_connections = list(set(sim_df.columns.tolist() + row_df.columns.tolist()))
    processed_lines = set()
    final_interconnections_list = []

    for line in all_connections:
        if line in processed_lines:
            continue

        # Extract nodes, handling potential whitespace variations
        if ' -> ' in line:
            nodes = line.split(' -> ')
        else:
            nodes = line.split('->')

        if len(nodes) != 2:
            logging.warning(f"Skipping improperly formatted line: {line}")
            final_interconnections_list.append(line)
            processed_lines.add(line)
            continue

        node1 = nodes[0].strip()
        node2 = nodes[1].strip()
        # Define the canonical direction (lexicographically first node first)
        canonical_line = f"{min(node1, node2)} -> {max(node1, node2)}"
        original_line = line # Keep track of the line name we encountered first
        opposite_line = f"{node2} -> {node1}"

        # Determine which line name to keep (the canonical one)
        line_to_keep = canonical_line
        line_to_drop = opposite_line if original_line == line_to_keep else original_line

        # Check if the opposite line exists in the original list (and hasn't been processed)
        opposite_exists = opposite_line in all_connections and opposite_line not in processed_lines

        if opposite_exists:
            logging.warning(f"Merging symmetrical lines: {original_line} and {opposite_line} into {line_to_keep}")

            # --- Process interconnections_sim (NTC) ---
            df_name = "interconnections_sim"
            line_keep_exists_sim = line_to_keep in sim_df.columns
            line_drop_exists_sim = line_to_drop in sim_df.columns

            if line_keep_exists_sim and line_drop_exists_sim:
                try:
                    val_keep = sim_df[line_to_keep]
                    val_drop = sim_df[line_to_drop]
                    # Check for significant difference
                    with np.errstate(divide='ignore', invalid='ignore'): # Ignore division by zero or NaN issues for the check
                        mean_val = (val_keep + val_drop) / 2
                        diff_percent = np.abs(val_keep - val_drop) / mean_val * 100
                        diff_percent = diff_percent.fillna(0) # Treat NaN percentages (e.g., from 0/0) as no difference
                        max_diff = diff_percent.max() if isinstance(diff_percent, pd.Series) else diff_percent

                    if max_diff > 20:
                         logging.warning(f"CRITICAL WARNING: NTC values for {line_to_keep} and {line_to_drop} differ by more than 20% (max diff: {max_diff:.2f}%). Averaging them.", )

                    # Average the values
                    avg_col = mean_val # Use the already calculated mean
                    sim_df[line_to_keep] = avg_col
                    sim_df.drop(columns=[line_to_drop], inplace=True)

                except TypeError:
                     logging.warning(f"Columns {line_to_keep}, {line_to_drop} in {df_name} are not numeric. Keeping values from {line_to_keep} if available.")
                     if line_to_drop in sim_df.columns:
                          sim_df.drop(columns=[line_to_drop], inplace=True) # Drop the non-numeric opposite anyway

            elif line_keep_exists_sim:
                logging.debug(f"Symmetrical line {line_to_drop} not found in {df_name} for {line_to_keep}. Keeping {line_to_keep}.")
                # Value already exists under the canonical name
            elif line_drop_exists_sim:
                logging.warning(f"NTC data only found for {line_to_drop} in {df_name}. Renaming to {line_to_keep}.")
                sim_df.rename(columns={line_to_drop: line_to_keep}, inplace=True)
            # If neither exists, do nothing for sim_df

            # --- Process interconnections_row (Fixed Flows) ---
            # Calculate net flow: flow(A->B) = flow(A->B)_orig - flow(B->A)_orig
            df_name = "interconnections_row"
            line_keep_exists_row = line_to_keep in row_df.columns
            line_drop_exists_row = line_to_drop in row_df.columns
            
            # Only process if at least one of the lines exists in row_df
            if line_keep_exists_row or line_drop_exists_row:
                flow_keep_dir = row_df.get(line_to_keep, 0) # Default to 0 if column missing
                flow_drop_dir = row_df.get(line_to_drop, 0) # Default to 0 if column missing

                # Ensure Series alignment if flow_keep_dir/flow_drop_dir are Series
                if isinstance(flow_keep_dir, pd.Series) or isinstance(flow_drop_dir, pd.Series):
                     common_index = row_df.index # Use a common index
                     if not isinstance(flow_keep_dir, pd.Series):
                         flow_keep_dir = pd.Series(flow_keep_dir, index=common_index)
                     if not isinstance(flow_drop_dir, pd.Series):
                         flow_drop_dir = pd.Series(flow_drop_dir, index=common_index)
                     # Reindex to ensure alignment before subtraction
                     flow_keep_dir = flow_keep_dir.reindex(common_index, fill_value=0)
                     flow_drop_dir = flow_drop_dir.reindex(common_index, fill_value=0)

                net_flow = flow_keep_dir - flow_drop_dir
                row_df[line_to_keep] = net_flow # Assign calculated net flow
                if line_drop_exists_row:
                    row_df.drop(columns=[line_to_drop], inplace=True)


            # --- Process price_transmission_raw (Price) ---
            df_name = "price_transmission_raw"
            line_keep_exists_price = line_to_keep in price_df.columns
            line_drop_exists_price = line_to_drop in price_df.columns

            if line_keep_exists_price and line_drop_exists_price:
                try:
                    val_keep = price_df[line_to_keep]
                    val_drop = price_df[line_to_drop]

                    # Check for significant difference
                    with np.errstate(divide='ignore', invalid='ignore'):
                        mean_val = (val_keep + val_drop) / 2
                        diff_percent = np.abs(val_keep - val_drop) / mean_val * 100
                        diff_percent = diff_percent.fillna(0)
                        max_diff = diff_percent.max() if isinstance(diff_percent, pd.Series) else diff_percent

                    if max_diff > 20:
                        logging.warning(f"CRITICAL WARNING: Transmission prices for {line_to_keep} and {line_to_drop} differ by more than 20% (max diff: {max_diff:.2f}%). Averaging them.")

                    # Average the values
                    avg_col = mean_val
                    price_df[line_to_keep] = avg_col
                    price_df.drop(columns=[line_to_drop], inplace=True)

                except TypeError:
                     logging.warning(f"Columns {line_to_keep}, {line_to_drop} in {df_name} are not numeric. Keeping values from {line_to_keep} if available.")
                     if line_to_drop in price_df.columns:
                          price_df.drop(columns=[line_to_drop], inplace=True) # Drop the non-numeric opposite anyway

            elif line_keep_exists_price:
                 logging.debug(f"Symmetrical line {line_to_drop} not found in {df_name} for {line_to_keep}. Keeping {line_to_keep}.")
                 # Value already exists under the canonical name
            elif line_drop_exists_price:
                logging.warning(f"Price data only found for {line_to_drop} in {df_name}. Renaming to {line_to_keep}.")
                price_df.rename(columns={line_to_drop: line_to_keep}, inplace=True)
             # If neither exists, do nothing for price_df


            # Mark both lines as processed and add the canonical line to the final list
            final_interconnections_list.append(line_to_keep)
            processed_lines.add(original_line)
            processed_lines.add(opposite_line)

        elif line not in processed_lines: # Process lines without a symmetrical counterpart found in the list
             final_interconnections_list.append(line)
             processed_lines.add(line)


    return sim_df, row_df, price_df


def incidence_matrix(sets, set_used, parameters, param_used, nodes='n'):
    """
    This function generates the incidence matrix of the lines within the nodes
    A particular case is considered for the node "Rest Of the World", which is no explicitly defined in DispaSET

    :param sets:        Dictionary with the sets
    :param set_used:    The set to be used (e.g., 'l')
    :param parameters:  Dictionary with the parameters
    :param param_used:  The parameter to be used
    :param nodes:       The nodes to be considered (default: 'n')
    :return:            Updated parameter
    """
    for i, l in enumerate(sets[set_used]):
        # Handle both cases: with and without spaces around the arrow
        if ' -> ' in l:
            [from_node, to_node] = l.split(' -> ')
        else:
            [from_node, to_node] = l.split('->')
        if (from_node.strip() in sets[nodes]) and (to_node.strip() in sets[nodes]):
            parameters[param_used]['val'][i, sets[nodes].index(to_node.strip())] = 1
            parameters[param_used]['val'][i, sets[nodes].index(from_node.strip())] = -1
        elif (from_node.strip() in sets[nodes]) and (to_node.strip() == 'RoW'):
            parameters[param_used]['val'][i, sets[nodes].index(from_node.strip())] = -1
        elif (from_node.strip() == 'RoW') and (to_node.strip() in sets[nodes]):
            parameters[param_used]['val'][i, sets[nodes].index(to_node.strip())] = 1
        else:
            logging.error("The line " + str(
                l) + " contains unrecognized nodes (" + from_node.strip() + ' or ' + to_node.strip() + ")")

    return parameters[param_used]


def adjust_ntc(inputs, value=None, write_gdx=False, dest_path=''):
    """
    Function used to modify the net transfer capacities in the Dispa-SET generated input data
    The function update the Inputs.p file in the simulation directory at each call

    :param inputs:      Input data dictionary OR path to the simulation directory containing Inputs.p
    :param value:       Absolute value of the desired capacity (! Applied only if scaling != 1 !)
    :param write_gdx:   boolean defining if Inputs.gdx should be also overwritten with the new data
    :param dest_path:   Simulation environment path to write the new input data. If unspecified, no data is written!
    :return:            New SimData dictionary
    """
    import pickle
    import os
    import shutil
    from ..misc.gdx_handler import write_variables

    if isinstance(inputs, str):
        path = inputs
        inputfile = path + '/Inputs.p'
        if not os.path.exists(path):
            sys.exit('Path + "' + path + '" not found')
        with open(inputfile, 'rb') as f:
            SimData = pickle.load(f)
    elif isinstance(inputs, dict):
        SimData = inputs
        path = SimData['config']['SimulationDirectory']
    else:
        logging.error('The input data must be either a dictionary or string containing a valid directory')
        sys.exit(1)

    if value is not None:
        SimData['parameters']['FlowMaximum']['val']=SimData['parameters']['FlowMaximum']['val']*value
    else:
        pass

    if dest_path == '':
        logging.info('Not writing any input data to the disk')
    else:
        if not os.path.isdir(dest_path):
            shutil.copytree(path, dest_path)
            logging.info('Created simulation environment directory ' + dest_path)
        logging.info('Writing input files to ' + dest_path)
        with open(os.path.join(dest_path, 'Inputs.p'), 'wb') as pfile:
            pickle.dump(SimData, pfile, protocol=pickle.HIGHEST_PROTOCOL)
        if write_gdx:
            write_variables(SimData['config'], 'Inputs.gdx', [SimData['sets'], SimData['parameters']])
            shutil.copy('Inputs.gdx', dest_path + '/')
            os.remove('Inputs.gdx')
    return SimData


def create_price_transmission_data(interconnections_sim, flows, price_transmission_raw, default_price=0):
    """
    Create price transmission data dictionary for all interconnections.
    
    :param interconnections_sim:     DataFrame with simulated interconnections
    :param flows:                    DataFrame with flows data
    :param price_transmission_raw:   DataFrame with raw price transmission data
    :param default_price:            Default price to use when no value is provided
    :return:                         DataFrame with price transmission data for all lines
    """
    price_transmission_data = {}
    
    for l in (interconnections_sim.columns.tolist() + flows.columns.tolist()):
        if l in price_transmission_raw:
            price_transmission_data[l] = price_transmission_raw[l].reindex(interconnections_sim.index)
        else:
            price_transmission_data[l] = pd.Series(default_price, index=interconnections_sim.index)
            if default_price > 0:
                logging.warning(f'No detailed values were found the transmission prices of line {l}. '
                               f'Using default value {default_price}')
    
    # Create DataFrame all at once to avoid fragmentation
    return pd.DataFrame(price_transmission_data, index=interconnections_sim.index) 