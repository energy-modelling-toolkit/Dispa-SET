"""
This file gathers different functions used in the DispaSET pre-processing tools

@author: Sylvain Quoilin
"""

from __future__ import division

import logging
import os
import shutil
import sys
import itertools
import numpy as np
import pandas as pd

from ..common import commons
from ..misc.gdx_handler import write_variables
from ..misc.str_handler import clean_strings, shrink_to_64


def pd_timestep(hours):
    """
    Function that converts time steps in hours into pandas frequencies (e.g '1h', '15min', ...)
    """
    if not isinstance(hours, (int, float)):
        logging.critical('Time steps must be provided in hours (integer or float number')
        sys.exit(1)
    if hours == 1:
        return '1h'
    elif hours == 0.25:
        return '15min'
    elif hours == 24:
        return '24h'
    else:
        return ''


def EfficiencyTimeSeries(config, plants):
    """
    Function that calculates an efficiency time series for each unit
    In case of generation unit, the efficiency is constant in time (for now)

    :param config:          Dispa-SET config file
    :param plants:          Pandas dataframe with the original list of units

    :returns:               Dataframe with a time series of the efficiency for each unit
    """
    Efficiencies = pd.DataFrame(columns=plants.index, index=config['idx_long'])
    for u in plants.index:
        Efficiencies[u] = plants.loc[u, 'Efficiency']
    return Efficiencies


def BoundarySectorEfficiencyTimeSeries(config, plants, zones_bs):
    """
    Function that calculates boundary sector efficiency time series for each unit
    In case of generation unit, the efficiency is defined as the EfficiencySectorX
    In case of of power consumption units, the efficiency is defined as the ChargingEfficiencySectorX, (for now)

    :param config:          Dispa-SET config file
    :param plants:          Pandas dataframe with the original list of units
    :param zones_bs:        Boundary sector zones

    :returns:               Dictionary of Dataframes with a time series of the efficiency for each unit
    """
    # Get the sector columns (Sector1, Sector2, etc.)
    sector_cols = [col for col in plants if col.startswith('Sector') and not col.startswith('Efficiency')]
    Efficiencies = {}
    
    for n in zones_bs:
        Efficiencies[n] = pd.DataFrame(columns=plants.index, index=config['idx_long'])
        for s in sector_cols:
            for u in plants.index:
                if ((plants.loc[u, 'Technology'] in commons['tech_p2bs'] or plants.loc[u, 'Technology'] in commons['tech_boundary_sector']) and 
                    (plants.loc[u, s] == n)):
                    # For p2x units, use EfficiencySectorX for Power2XConversionMultiplier
                    # For x2p units, use EfficiencySectorX for X2PowerConversionMultiplier
                    eff_col = 'EfficiencySector' + s[6:]  # Convert 'Sector1' to 'EfficiencySector1'
                    Efficiencies[n][u] = plants.loc[u, eff_col]
        if n in Efficiencies:
            Efficiencies[n] = Efficiencies[n].fillna(0)
    
    return {'Efficiency': Efficiencies}


def select_units(units, config):
    """
    Function returning a new list of units by removing the ones that have unknown
    technology, zero capacity, or unknown zone

    :param units:       Pandas dataframe with the original list of units
    :param config:      Dispa-SET config dictionnary
    :return:            New list of units
    """
    for unit in units.index:
        if units.loc[unit, 'Technology'] == 'Other':
            logging.warning('Removed Unit ' + str(units.loc[unit, 'Unit']) + ' since its technology is unknown')
            units.drop(unit, inplace=True)
        elif (units.loc[unit, 'PowerCapacity'] == 0) and ((units.loc[unit, 'STOMaxChargingPower'] == 0) or (
                type(units.loc[unit, 'STOMaxChargingPower']) == np.float64)) and (units.loc[unit, 'Technology'] not in
                                                                                  commons['tech_p2bs']):
            logging.warning('Removed Unit ' + str(units.loc[unit, 'Unit']) + ' since it has a null capacity')
            units.drop(unit, inplace=True)
        elif units.loc[unit, 'Zone'] not in config['zones']:
            logging.warning('Removed Unit ' + str(units.loc[unit, 'Unit']) + ' since its zone (' + str(
                units.loc[unit, 'Zone']) + ') is not in the list of zones')
            units.drop(unit, inplace=True)
    units.index = range(len(units))
    return units


def incidence_matrix(sets, set_used, parameters, param_used, nodes='n'):
    """
    This function generates the incidence matrix of the lines within the nodes
    A particular case is considered for the node "Rest Of the World", which is no explicitely defined in DispaSET
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


def interconnections(Simulation_list, NTC_inter, Historical_flows):
    """
    Function that checks for the possible interconnections of the zones included
    in the simulation. If the interconnections occurs between two of the zones
    defined by the user to perform the simulation with, it extracts the NTC between
    those two zones. If the interconnection occurs between one of the zones
    selected by the user and one country outside the simulation, it extracts the
    physical flows; it does so for each pair (country inside-country outside) and
    sums them together creating the interconnection of this country with the RoW.

    :param Simulation_list:     List of simulated zones
    :param NTC_inter:                 Day-ahead net transfer capacities (pd dataframe)
    :param Historical_flows:    Historical flows (pd dataframe)
    """
    index = NTC_inter.index.tz_localize(None).intersection(Historical_flows.index.tz_localize(None))
    if len(index) == 0:
        logging.error('The two input dataframes (NTCs and Historical flows) must have the same index. '
                      'No common values have been found')
        sys.exit(1)
    elif len(index) < len(NTC_inter) or len(index) < len(Historical_flows):
        diff = np.maximum(len(Historical_flows), len(NTC_inter)) - len(index)
        logging.warning('The two input dataframes (NTCs and Historical flows) do not share the same index, '
                        'although some values are common. The intersection has been considered and ' + str(diff) +
                        ' data points have been lost')
    # Checking that all values are positive:
    if (NTC_inter.values < 0).any():
        pos = np.where(NTC_inter.values < 0)
        logging.warning('At least one NTC value is negative, for example in line ' + str(NTC_inter.columns[pos[1][0]]) +
                        ' and time step ' + str(NTC_inter.index[pos[0][0]]))
    all_connections = []
    simulation_connections = []
    # List all connections from the dataframe headers:
    ConList = Historical_flows.columns.tolist() + [x for x in NTC_inter.columns.tolist() if
                                                   x not in Historical_flows.columns.tolist()]
    for connection in ConList:
        # Ensure connection has the correct format with "->" separator
        if ' -> ' in connection:
            z = connection.split(' -> ')
        else:
            z = connection.split('->')
        if len(z) == 2:
            connection = z[0].strip() + ' -> ' + z[1].strip()
            z = [z[0].strip(), z[1].strip()]  # Update z with stripped values
        
        # Check if this connection is between simulated zones
        if z[0] in Simulation_list and z[1] in Simulation_list:
            # Internal connection between simulated zones
            all_connections.append(connection)
            simulation_connections.append(connection)
        elif z[0] in Simulation_list or z[1] in Simulation_list:
            # External connection with at least one simulated zone
            all_connections.append(connection)

    df_zones_simulated = pd.DataFrame(index=index)
    for interconnection in simulation_connections:
        # Handle both formats when looking up in NTC_inter
        if interconnection in NTC_inter.columns:
            df_zones_simulated[interconnection] = NTC_inter[interconnection]
            logging.info('Detected interconnection ' + interconnection +
                         '. The historical NTCs will be imposed as maximum flow value if the NTC method is used')
        else:
            # Try without spaces around arrow
            no_spaces = interconnection.replace(' -> ', '->')
            if no_spaces in NTC_inter.columns:
                df_zones_simulated[interconnection] = NTC_inter[no_spaces]
                logging.info('Detected interconnection ' + interconnection +
                            '. The historical NTCs will be imposed as maximum flow value if the NTC method is used')
    interconnections1 = df_zones_simulated.columns

    # Display a warning if a zone is isolated:
    for z in Simulation_list:
        if not any([z in conn for conn in interconnections1]) and len(Simulation_list) > 1:
            logging.warning('Zone ' + z + ' does not appear to be connected to any other zone in the NTC table. '
                            'It should be simulated in isolation')

    df_RoW_temp = pd.DataFrame(index=index)
    connNames = []
    for interconnection in all_connections:
        # Handle both formats when looking up in Historical_flows
        if interconnection in Historical_flows.columns:
            if interconnection not in simulation_connections:
                df_RoW_temp[interconnection] = Historical_flows[interconnection]
                connNames.append(interconnection)
        else:
            # Try without spaces around arrow
            no_spaces = interconnection.replace(' -> ', '->')
            if no_spaces in Historical_flows.columns and no_spaces not in simulation_connections:
                df_RoW_temp[interconnection] = Historical_flows[no_spaces]
                connNames.append(interconnection)

    # Process flows between simulated zones and RoW
    df_zones_RoW = pd.DataFrame(index=index)
    
    # First collect all the simulated zones involved in external flows
    simulated_zones_in_external_flows = set()
    for connection in connNames:
        if ' -> ' in connection:
            z = connection.split(' -> ')
        else:
            z = connection.split('->')
        
        zone_from = z[0].strip()
        zone_to = z[1].strip()
        
        if zone_from in Simulation_list:
            simulated_zones_in_external_flows.add(zone_from)
        if zone_to in Simulation_list:
            simulated_zones_in_external_flows.add(zone_to)
    
    # Process each simulated zone with external connections
    # In the new approach, all flows for a zone with RoW will be merged into a single entry
    # where positive flow means export from zone to RoW, negative flow means import to zone from RoW
    for zone in simulated_zones_in_external_flows:
        net_flows = pd.Series(0, index=index)  # Initialize a series of zeros with the right index
        
        for connection in connNames:
            if ' -> ' in connection:
                z = connection.split(' -> ')
            else:
                z = connection.split('->')
            
            zone_from = z[0].strip()
            zone_to = z[1].strip()
            
            # Zone -> External: Add as positive flow (export)
            if zone_from == zone and zone_to not in Simulation_list:
                net_flows += df_RoW_temp[connection]
                logging.info(f'Detected interconnection {connection}, happening between simulated zone {zone} '
                             f'and the rest of the world. Adding as positive flow to {zone} -> RoW')
            
            # External -> Zone: Add as negative flow (import)
            elif zone_to == zone and zone_from not in Simulation_list:
                net_flows -= df_RoW_temp[connection]  # Subtract for imports
                logging.info(f'Detected interconnection {connection}, happening between the rest of the world '
                             f'and simulated zone {zone}. Adding as negative flow to {zone} -> RoW')
        
        # Store the net flow for this zone
        if not net_flows.equals(pd.Series(0, index=index)):  # Only add if there are non-zero flows
            flow_name = f'{zone} -> RoW'
            df_zones_RoW[flow_name] = net_flows
            logging.info(f'Created merged bidirectional flow {flow_name} (positive=export, negative=import)')
            
    interconnections2 = df_zones_RoW.columns
    inter = list(interconnections1) + list(interconnections2)
    return df_zones_simulated, df_zones_RoW, inter


# Helpers

def _mylogspace(low, high, N):
    """
    Self-defined logspace function in which low and high are the first and last values of the space
    """
    # shifting all values so that low = 1
    space = np.logspace(0, np.log10(high + low + 1), N) - (low + 1)
    return space


def _find_nearest(array, value):
    """
    Self-defined function to find the index of the nearest value in a vector
    """
    idx = (np.abs(array - value)).argmin()
    return idx


def _reverse_dict(dict_):
    """
    Reverse Dictionary (Key, Value) to (Value, Key)
    :param dict_:     Dictionary to reverse

    """
    new_dic = {}
    for k, v in dict_.items():
        for x in v:
            new_dic[x] = k
    return new_dic


def _split_list(list_):
    """
    Split list elements into string with " - " seperator 
    :param list_:     List to split

    """
    res = str()
    # remove empty elements from the list:
    newlist = [l for l in list_ if (str(l) != 'nan') and (str(l) != '')]
    for l in newlist:
        if l != newlist[-1]:
            res += str(l) + " - "
        else:
            res += str(l)
    return res


def _list2dict(list_, agg): return {key: agg for key in list_}


def _flatten_list(l):
    """
    Function that unfolds nested lists
    Example: 
        [1, 3, ['aa','bb'],4] is turned into [1,3, 'aa', 'bb', 4]
    """
    flat_list = []
    for sublist in l:
        if isinstance(sublist, list):
            for item in sublist:
                flat_list.append(item)
        else:
            flat_list.append(sublist)
    return flat_list


def _merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy.
    Used for compatibility Python 2 and 3
    inspired by: https://stackoverflow.com/questions/38987/how-to-merge-two-dictionaries-in-a-single-expression 
    """
    z = x.copy()
    z.update(y)
    return z


def _get_index(df_, idx):
    """Helper function to get former indexes and units from a dataframe"""
    former_indexes = [_flatten_list(list(df_.loc[i]['FormerIndexes'].values)) for i in idx]
    former_units = [_flatten_list(list(df_.loc[i]['FormerUnits'].values)) for i in idx]
    return former_indexes, former_units


def _create_mapping(merged_df):
    mapping = {"NewIndex": {}, 'FormerIndexes': merged_df['FormerIndexes'].to_dict()}
    mapping['NewIndex'] = _reverse_dict(mapping['FormerIndexes'])
    return mapping


def _new_unit_names(df_merged, df_, string_keys):
    # if merged unit, create name -> else take old name for unit
    keys = ['FormerIndexes'] + string_keys
    create_unit_name = lambda x: str(x.FormerIndexes) + " - " + df_.iloc[x.FormerIndexes[0]]['Unit'] if len(
        x.FormerIndexes) == 1 else shrink_to_64(clean_strings(_split_list(list(x[keys].values))))
    df_merged['Unit'] = df_merged.apply(create_unit_name, axis=1)
    return df_merged.set_index('Unit', drop=False)


def _linearize_ramping(plants):
    '''
    Function that converts the integer constraints for the power plants into ramping rates
    '''
    #ramping_up = (lambda row: min(1/ (row["MinDownTime"]+1E-9) / 60, row["RampUpRate"]))
    #ramping_down = (lambda row: min(1/ (row["MinUpTime"]+1E-9) / 60, row["RampDownRate"]))
    #plants["RampUpRate"] = plants.apply(ramping_up, axis=1)
    #plants["RampDownRate"] = plants.apply(ramping_down, axis=1)
    plants["RampUpRate"] = plants['PartLoadMin'] * ( 1 / (np.maximum(1,plants['MinDownTime'])*60 + 1e-9)) + (1 - plants['PartLoadMin'])*np.minimum(plants['RampUpRate'],1/60)
    plants["RampDownRate"] = plants['PartLoadMin'] * ( 1 / (np.maximum(1,plants['MinUpTime'])*60 + 1e-9)) + (1 - plants['PartLoadMin'])*np.minimum(plants['RampDownRate'],1/60)


def group_plants(plants, method, df_grouped=False, group_list=None):
    """
    This function returns the final dataframe with the merged units and their characteristics

    :param plants:          Pandas dataframe with each power plant and their characteristics
                            (following the DispaSET format)
    :param method:          Select clustering method ('Standard'/'LP'/None)
    :param df_grouped:      Set to True if this plants dataframe has already been grouped and contains the column
                            "FormerIndexes"
    :param group_list:      List of columns whose values must be identical in order to group two units
    :return:                A list with the merged plants and the mapping between the original and merged units

    """
    # Definition of the merged power plants dataframe:
    # if (group_list is None) and ((plants['Zone_th'] != np.nan).all()) and ((plants['Zone_h2'] != np.nan).all()):
    #     group_list = ['Zone', 'Zone_th', 'Zone_h2', 'Technology', 'Fuel', 'CHPType']
    if (group_list is None) and ('FuelPricebyUnit' in plants.columns) and ((plants['Sector1'] != np.nan).all()):
        group_list = ['Zone','FuelPricebyUnit', 'Sector1', 'Technology', 'Fuel', 'CHPType']
    elif (group_list is None) and ((plants['Sector1'] != np.nan).all()):
        group_list = ['Zone', 'Sector1', 'Technology', 'Fuel', 'CHPType']
    elif (group_list is None) and ((plants['Sector1'] != np.nan).all()) and ((plants['Sector2'] != np.nan).all()):
        group_list = ['Zone', 'Sector1', 'Technology', 'Fuel', 'CHPType', 'Sector2']
    else:
        group_list = ['Zone', 'Technology', 'Fuel', 'CHPType']
    plants_merged = pd.DataFrame(columns=plants.columns)
    grouped = plants.groupby(group_list, as_index=False)
    agg_dict = create_agg_dict(plants, method=method)
    plants_merged = pd.concat([plants_merged, grouped.agg(agg_dict)])
    if method == "Integer clustering":
        plants_merged['StartUpCost'] = plants_merged['StartUpCost'] / plants_merged['Nunits']
        plants_merged['NoLoadCost'] = plants_merged['NoLoadCost'] / plants_merged['Nunits']
    idx = [list(i.values) for i in list(grouped.groups.values())]
    if not df_grouped:
        plants_merged['FormerIndexes'] = [list(plants.loc[i]['index'].values) for i in idx]
        plants_merged['FormerUnits'] = [list(plants.loc[i]['Unit'].values) for i in idx]
    else:
        # case in which the plants have already been clustered once => nested lists in FormerIndexes
        former_indexes, former_units = _get_index(plants, idx)
        plants_merged['FormerIndexes'] = list(former_indexes)
        plants_merged['FormerUnits'] = list(former_units)

    return plants_merged


def create_agg_dict(df_, method="Standard"):
    """
    This function returns a dictionnary with the proper aggregation method
    for each columns of the units table, depending on the clustering method

    Author: Matthias Zech
    """


    # lambda functions for other aggregations than standard aggregators like min/max,...
    wm_pcap = lambda x: np.average(x.astype(float),
                                   weights=df_.loc[x.index, "PowerCapacity"])  # weighted mean with weight=PowerCapacity
    wm_nunit = lambda x: np.average(x.astype(float),
                                    weights=df_.loc[x.index, "Nunits"])  # weighted mean with weight=NUnits
    get_ramping_cost = lambda x: wm_pcap(
        (1 - df_.loc[x.index, "PartLoadMin"]) * x + df_.loc[x.index, "StartUpCost"] / df_.loc[x.index, "PowerCapacity"])
    min_load = lambda x: np.min(x * df_.loc[x.index, "PowerCapacity"]) / df_.loc[x.index, "PowerCapacity"].sum()

    if method in ("Standard", "MILP"):
        sum_cols = ["PowerCapacity", "STOCapacity", "STOMaxChargingPower", "InitialPower", "CHPMaxHeat"]
        weighted_avg_cols = [
                            "RampUpRate",
                            "RampDownRate",
                            "MinUpTime",
                            "MinDownTime",
                            "StartUpCost",
                            "NoLoadCost",
                            "Efficiency",
                            "MinEfficiency",
                            "STOChargingEfficiency",
                            "CO2Intensity",
                            "STOSelfDischarge",
                            "CHPPowerToHeat",
                            "CHPPowerLossFactor",
                            'COP',
                            'TNominal',
                            'coef_COP_a',
                            'coef_COP_b',
                            'WaterConsumption',
                            'WaterWithdrawal',
                            'RampingCost'
                            ]
        min_cols = ["StartUpTime"]
        nunits = ["Nunits"]

        # Define aggregators
        agg_dict = _list2dict(sum_cols, 'sum')
        agg_dict = _merge_two_dicts(agg_dict, _list2dict(weighted_avg_cols, wm_pcap))
        agg_dict = _merge_two_dicts(agg_dict, _list2dict(min_cols, 'min'))
        agg_dict = _merge_two_dicts(agg_dict, _list2dict(['PartLoadMin'], min_load))
        # agg_dict = _merge_two_dicts(agg_dict, _list2dict(ramping_cost, get_ramping_cost))
        agg_dict = _merge_two_dicts(agg_dict, _list2dict(nunits, lambda x: 1))
        agg_dict = dict((k, v) for k, v in agg_dict.items() if k in df_.columns)  # remove unnecesary columns
        return agg_dict

    elif method == "LP clustered":
        sum_cols = ["PowerCapacity", "STOCapacity", "STOMaxChargingPower", "InitialPower", "CHPMaxHeat"]
        weighted_avg_cols = [
                            "RampUpRate",
                            "RampDownRate",
                            "MinUpTime",
                            "MinDownTime",
                            "StartUpCost",
                            "NoLoadCost",
                            "Efficiency",
                            "MinEfficiency",
                            'PartLoadMin',
                            "STOChargingEfficiency",
                            "CO2Intensity",
                            "STOSelfDischarge",
                            "CHPPowerToHeat",
                            "CHPPowerLossFactor",
                            'COP',
                            'TNominal',
                            'coef_COP_a',
                            'coef_COP_b',
                            'WaterConsumption',
                            'WaterWithdrawal',
                            'RampingCost'
                            ]
        min_cols = ["StartUpTime"]
        # ramping_cost = ["RampingCost"]
        nunits = ["Nunits"]

        # Define aggregators
        agg_dict = _list2dict(sum_cols, 'sum')
        agg_dict = _merge_two_dicts(agg_dict, _list2dict(weighted_avg_cols, wm_pcap))
        agg_dict = _merge_two_dicts(agg_dict, _list2dict(min_cols, 'min'))
        #agg_dict = _merge_two_dicts(agg_dict, _list2dict(['PartLoadMin'], lambda x: 0))
        # agg_dict = _merge_two_dicts(agg_dict, _list2dict(ramping_cost, get_ramping_cost))
        agg_dict = _merge_two_dicts(agg_dict, _list2dict(nunits, lambda x: 1))
        agg_dict = dict((k, v) for k, v in agg_dict.items() if k in df_.columns)  # remove unnecesary columns
        return agg_dict

    elif method == "Integer clustering":
        sum_cols = ["Nunits", "StartUpCost",'NoLoadCost',]

        weighted_avg_cols = ['PowerCapacity',
                             'RampUpRate',
                             'RampDownRate',
                             'MinUpTime',
                             'MinDownTime',
                             'Efficiency',
                             'MinEfficiency',
                             'STOChargingEfficiency',
                             'CO2Intensity',
                             'STOSelfDischarge',
                             'STOCapacity',
                             'STOMaxChargingPower',
                             'PartLoadMin',
                             'StartUpTime',
                             'RampingCost',
                             'CHPPowerToHeat',
                             'CHPPowerLossFactor',
                             'CHPMaxHeat',
                             'COP',
                             'TNominal',
                             'coef_COP_a',
                             'coef_COP_b',
                             'WaterConsumption',
                             'WaterWithdrawal'
                             ]

        # Define aggregators
        agg_dict = _list2dict(sum_cols, 'sum')
        agg_dict = _merge_two_dicts(agg_dict, _list2dict(weighted_avg_cols, wm_nunit))
        agg_dict = dict((k, v) for k, v in agg_dict.items() if k in df_.columns)  # remove unnecesary columns

        return agg_dict

    else:
        logging.critical('Clustering method not properly specified. Should be one of the following options:'
                         ' LP clustered, MILP, Standard, Integer clustering')
        sys.exit(1)


def clustering(plants_in, method="Standard", Nslices=20, PartLoadMax=0.1, Pmax=30):
    """
    Merge excessively disaggregated power Units.

    :param plants_in:       Pandas dataframe with each power plant and their characteristics
                            (following the DispaSET format)
    :param method:          Select clustering method ('Standard'/'LP'/None)
    :param Nslices:         Number of slices used to fingerprint each power plant characteristics.
                            Slices in the power plant data to categorize them
                            (fewer slices involves that the plants will be aggregated more easily)
    :param PartLoadMax:     Maximum part-load capability for the unit to be clustered
    :param Pmax:            Maximum power for the unit to be clustered
    :return:                A list with the merged plants and the mapping between the original and merged units

    @author: Matthias Zech
    """

    # do not alter the original plants table:
    plants = plants_in.copy()
    # Checking the the required columns are present in the input pandas dataframe:
    required_inputs = ['Unit', 'PowerCapacity', 'PartLoadMin', 'RampUpRate', 'RampDownRate', 'StartUpTime',
                       'MinUpTime', 'MinDownTime', 'NoLoadCost', 'StartUpCost', 'Efficiency']
    for input_value in required_inputs:
        if input_value not in plants.columns:
            logging.error("The plants dataframe requires a '" + input_value + "' column for clustering")
            sys.exit(1)
    if "Nunits" not in plants:
        plants["Nunits"] = 1
    plants['PowerCapacity'] = plants['PowerCapacity'].astype(float)
    plants.loc[plants['PowerCapacity'] == 0, 'PowerCapacity'] = 1e-9

    Nunits = len(plants)
    plants.index = range(Nunits)
    plants_merged = pd.DataFrame(columns=plants.columns)

    # Fill nan values:
    string_keys = ['Zone', 'Technology', 'Fuel', 'CHPType', 'Sector1']
    for key in string_keys:
        plants[key] = plants[key].fillna("")
    for key in ['PartLoadMin', 'StartUpTime', 'MinUpTime', 'MinDownTime', 'NoLoadCost', 'StartUpCost',
                'WaterWithdrawal', 'WaterConsumption']:
        plants[key] = plants[key].fillna(0)
    for key in ['RampUpRate', 'RampDownRate']:
        plants[key] = plants[key].fillna(1e9)

    # Checking the validity of the selected clustering method
    plants["index"] = plants.index

    OnlyOnes = (plants["Nunits"] == 1).all()
    if method in ["Standard", "MILP"]:
        if OnlyOnes:
            ####### Three cluster groups in the standard MILP formulation
            ###### 1) Highly flexible
            ###### 2) Low Pmin
            ###### 3) Similar characteristics --> similarity expressed via fingerprints
            # First, cluster by same string keys and flexible and low_pmax
            # Join grouped data with inflexible and no low_pmmax data
            # Group joined dataframe by string keys including same technical characteristics using fingerprints
            # The more Nslices, the more heterogenity between data, the less is merged
            # Definition of the fingerprint value of each power plant, i.e. the pattern of the slices number in
            # which each of its characteristics falls:

            # helper_cols = ['flex', 'low_pmin', 'low_pmax', 'fingerprints']
            highly_flexible = (
                    (plants["RampUpRate"] > 1 / 60)
                    & (plants["RampDownRate"] > 1 / 60)
                    & (plants["StartUpTime"] < 1)
                    & (plants["MinDownTime"] <= 1)
                    & (plants["MinUpTime"] <= 1)
            )

            low_pmax = plants["PowerCapacity"] <= Pmax
            plants["flex"] = highly_flexible
            plants["low_pmax"] = low_pmax
            plants["FormerIndexes"] = pd.Series(plants.index.values).apply(lambda x: [x])
            plants["FormerUnits"] = pd.Series(plants['Unit'].values).apply(lambda x: [x])

            condition = (plants["low_pmax"]) | (plants["flex"])
            first_cluster = plants[condition]  # all data without other clustering
            first_cluster = group_plants(first_cluster, method, False, string_keys)

            # first_cluster = first_cluster.append(plants[~condition], ignore_index=True)
            first_cluster = pd.concat([first_cluster, plants[~condition]], ignore_index=True)
            # Slicing:
            bounds = {
                "PartLoadMin": np.linspace(0, 1, Nslices),
                "RampUpRate": np.linspace(0, 1, Nslices),
                "RampDownRate": np.linspace(0, 1, Nslices),
                "StartUpTime": _mylogspace(0, 36, Nslices),
                "MinUpTime": _mylogspace(0, 168, Nslices),
                "MinDownTime": _mylogspace(0, 168, Nslices),
                "NoLoadCost": np.linspace(0, 50, Nslices),
                "StartUpCost": np.linspace(0, 500, Nslices),
                "Efficiency": np.linspace(0, 1, Nslices),
                "WaterWithdrawal": np.linspace(0, 200, 250),
                "WaterConsumption": np.linspace(0, 20, Nslices),
            }

            fingerprints = []
            for i in first_cluster.index:
                fingerprints.append(
                    [
                        _find_nearest(bounds["PartLoadMin"], first_cluster["PartLoadMin"][i]),
                        _find_nearest(bounds["RampUpRate"], first_cluster["RampUpRate"][i]),
                        _find_nearest(bounds["RampDownRate"], first_cluster["RampDownRate"][i]),
                        _find_nearest(bounds["StartUpTime"], first_cluster["StartUpTime"][i]),
                        _find_nearest(bounds["MinUpTime"], first_cluster["MinUpTime"][i]),
                        _find_nearest(bounds["MinDownTime"], first_cluster["MinDownTime"][i]),
                        _find_nearest(bounds["NoLoadCost"], first_cluster["NoLoadCost"][i]),
                        _find_nearest(bounds["StartUpCost"], first_cluster["StartUpCost"][i]),
                        _find_nearest(bounds["Efficiency"], first_cluster["Efficiency"][i]),
                        _find_nearest(bounds["WaterConsumption"], first_cluster["WaterConsumption"][i]),
                        _find_nearest(bounds["WaterWithdrawal"], first_cluster["WaterWithdrawal"][i]),
                    ]
                )

            first_cluster["fingerprints"] = fingerprints

            # the elements of the list are irrelevant for the clustering
            first_cluster["fingerprints"] = first_cluster["fingerprints"].astype(str)
            low_pmin = first_cluster["PartLoadMin"] <= PartLoadMax
            if not first_cluster[low_pmin].empty:
                second_cluster = group_plants(first_cluster[low_pmin], method, True, string_keys + ["fingerprints"])
                # plants_merged = second_cluster.append(first_cluster[~low_pmin], ignore_index=True)
                plants_merged = pd.concat([second_cluster, first_cluster[~low_pmin]], ignore_index=True)
            else:
                plants_merged = first_cluster[:]

            plants = plants.drop(["flex", "low_pmax", "FormerIndexes"], axis=1)
            plants_merged = plants_merged.drop(["index", "fingerprints", "flex", "low_pmax"], axis=1)
        else:  # not all only ones
            logging.warning("The standard (or MILP) clustering method is only applicable if all values of the "
                            "Nunits column in the power plant data are set to one. At least one different value has "
                            "been encountered. No clustering will be applied")
            plants_merged = plants.copy()
            plants_merged["FormerIndexes"] = plants["index"].apply(lambda x: [x])
            plants_merged["FormerUnits"] = plants["Unit"].apply(lambda x: [x])

    elif method == "LP clustered":
        if not OnlyOnes:
            logging.warning("The LP clustering method aggregates all the units of the same type. Individual units are "
                            "not considered")
            list_mult = [
                        "PowerCapacity",
                        "STOCapacity",
                        "STOMaxChargingPower",
                        "InitialPower",
                        "CHPMaxHeat",
                        ]
            # Restricting the list of values to multiply to those who are present in the plants table:
            list_mult = [x for x in list_mult if x in plants]
            # Modifying the table to remove multiple-units plants:
            plants[list_mult] = plants[list_mult].multiply(plants["Nunits"], axis="index")
            plants["Nunits"] = 1
            OnlyOnes = True

        plants_merged = group_plants(plants, method="LP clustered")

    elif method == "LP":
        if not OnlyOnes:
            logging.warning("The LP method aggregates all identical units by multiplying by the Nunits variable")
            list_mult = [
                        "PowerCapacity",
                        "STOCapacity",
                        "STOMaxChargingPower",
                        "InitialPower",
                        "CHPMaxHeat",
                        ]
            # Restricting the list of values to multiply to those who are present in the plants table:
            list_mult = [x for x in list_mult if x in plants]
            # Modifying the table to remove multiple-units plants:
            plants[list_mult] = plants[list_mult].multiply(plants["Nunits"], axis="index")
            plants["Nunits"] = 1
            OnlyOnes = True

        plants_merged = plants
        # formers indexes and units:        
        plants_merged["FormerIndexes"] = plants["index"].apply(lambda x: [x])
        plants_merged["FormerUnits"] = plants["Unit"].apply(lambda x: [x])

    elif method == "Integer clustering":
        plants_merged = group_plants(plants, method="Integer clustering")
        # Correcting the Nunits field of the clustered plants (must be integer):

    elif method == "No clustering":
        plants_merged = plants.copy()
        plants_merged["FormerIndexes"] = plants["index"].apply(lambda x: [x])
        plants_merged["FormerUnits"] = plants["Unit"].apply(lambda x: [x])

    else:
        logging.error('Method argument ("' + str(method) + '") not recognized in the clustering function')
        sys.exit(1)

    plants_merged = _new_unit_names(plants_merged, plants, string_keys)
    # Modify the Unit names with the original index number. In case of merged plants,
    # indicate all indexes + the plant type and fuel

    mapping = _create_mapping(plants_merged)
    if Nunits != len(plants_merged):
        logging.info("Clustered " + str(Nunits) + " original units into " + str(len(plants_merged)) + " new units")
    else:
        logging.warning("Did not cluster any unit")

    # indexes of units which were not clustered:
    idx_merged = [i for i in plants_merged.index if len(plants_merged.loc[i, 'FormerIndexes']) == 1]
    idx_orig = [plants_merged.loc[i, 'FormerIndexes'][0] for i in idx_merged]
    columns = plants_merged.columns.drop(['Unit', 'FormerIndexes', 'FormerUnits'])
    plants_merged.loc[idx_merged, columns] = plants.loc[idx_orig, columns].values
    if method in ['LP','LP clustered']:
        # Transforming the min up/down times into ramping rates
        _linearize_ramping(plants_merged)
        # Transforming the start-up cost into ramping for the plants that did not go through any clustering:
        ramping_lbd = (lambda row: row["StartUpCost"] / row["PowerCapacity"] if row.RampingCost == 0
                       else row.RampingCost)
        plants_merged["RampingCost"] = plants_merged.apply(ramping_lbd, axis=1)
    # reorder columns:
    new_columns = [key for key in plants.columns if key in plants_merged]
    plants_merged = plants_merged[new_columns + list(plants_merged.columns.drop(new_columns))]
    return plants_merged, mapping

# TODO: CHECK THIS FUNCTION WITH BACKWARD COMPATIBILITY
def PTDF_matrix(config, feeder):
    '''
    Function that calculate the ptdf matrix from the grid data table
    '''
    zones = config['zones']

    # Add Line number (1-based)
    feeder['Line'] = np.arange(len(feeder)) + 1

    # Split Code_Line into Code_Node_i and Code_Node_j
    feeder['Code_Line'] = feeder.index
    split_nodes = feeder['Code_Line'].str.split('->', n=1, expand=True)
    if split_nodes.shape[1] < 2:
        raise ValueError(f"Error splitting 'Code_Line'. Ensure all lines use '->' separator. Problematic lines:\n{feeder[split_nodes[1].isna()]}")
    feeder['Code_Node_i'] = split_nodes[0].str.strip()
    feeder['Code_Node_j'] = split_nodes[1].str.strip()

    # Initialize Node ID columns
    feeder['Node_i'] = pd.NA
    feeder['Node_j'] = pd.NA

    # Assign Node IDs based on zone list order (1-based index)
    for zone_id, zone in enumerate(zones, 1):
        conditions_i = (feeder['Code_Node_i'] == zone)
        conditions_j = (feeder['Code_Node_j'] == zone)

        feeder.loc[conditions_i, 'Node_i'] = zone_id
        feeder.loc[conditions_j, 'Node_j'] = zone_id

    # Check for nodes that were not mapped
    unmapped_i = feeder['Node_i'].isna()
    unmapped_j = feeder['Node_j'].isna()
    if unmapped_i.any() or unmapped_j.any():
        missing_i_codes = feeder.loc[unmapped_i, 'Code_Node_i'].unique()
        missing_j_codes = feeder.loc[unmapped_j, 'Code_Node_j'].unique()
        all_missing = set(missing_i_codes) | set(missing_j_codes)
        logging.error(f"Could not map the following nodes found in GridData Code_Line to integer IDs "
                      f"based on config zones {zones}: {all_missing}. Please check GridData and config.")
        # Raise an error to prevent proceeding with invalid data
        raise ValueError(f"Unmapped nodes found in GridData: {all_missing}")

    # Convert node IDs to integer type
    # Use Int64 to potentially handle NA if error wasn't raised, but the check above should prevent NAs.
    feeder['Node_i'] = feeder['Node_i'].astype(int)
    feeder['Node_j'] = feeder['Node_j'].astype(int)

    #ptdf calculation
    
    HM=feeder[feeder['HM']!=0].index.tolist() #Lines for maintenance 
    lines_under_maintenance=sum(feeder['HM']!=0) #Number of lines under maintenance
    bina=list(itertools.product(*[[1,0]]*lines_under_maintenance)) 
    tab_bin = pd.DataFrame(bina,columns=list(np.where(feeder['HM']!=0)[0])) #Posibles combinaciones
    pt_dic={}
    ld_dic={}

    len_datapt=[0]
    len_datadf=[0]

    nodes=len(set.union(set(feeder['Node_i']),set(feeder['Node_j'])))
    
    availability=[1]*len(tab_bin)
    for k in range(len(tab_bin)):
        li_ava=[]
        combi=list(tab_bin.iloc[k])
        for m in range(len(combi)):
            if combi[m]==0:
                li_ava+=[HM[m]]
            else:
                continue
        new_feeder=feeder.drop(li_ava)  
        l=len(new_feeder)
        n=len(set.union(set(new_feeder['Node_i']),set(new_feeder['Node_j'])))

        Bp=np.diag([1/(1j*feeder['X'][mq]) for mq in range(l)])#Primitive Susceptances Matrix
        A=np.zeros((n,l))
        for ldw in range(l):
            # linenodes = feeder['Code_Line'][ldw].split(' -> ')
            # From=linenodes[0]
            # To=linenodes[1]
            From=feeder['Node_i'][ldw]
            To=feeder['Node_j'][ldw]
            for m in range(n):
                if From==m+1:
                    A[m,ldw]=1
                elif To==m+1:
                    A[m,ldw]=-1
                    
        if 0 in np.sum(abs(A),axis=1) or len(A)!=nodes:
            print('the combination ', combi, 'leaves one node isolated from the system, therefore it will not be taken into account')
            availability[k]=0
            pt_dic[k]=np.zeros((len(A[0,:]),nodes))
            ld_dic[k]=np.zeros((len(A[0,:]),len(A[0,:])))
            len_datapt+=[len_datapt[-1]+len(new_feeder['Line'])+2]
            len_datadf+=[len_datadf[-1]+len(new_feeder['Line'])+2]
            continue
            
        #--Matrix Bbus--
        Bbus=A@Bp@A.T
        
        #--Building the X matrix selecting the slack node from the variable "Slack"--
        # Find slack bus automatically: node with the most connections (common heuristic)
        max_connections = 0
        slack_bus = -1 # Initialize with invalid value
        for node in range(1, nodes + 1):
            # Count lines connected to this node
            num_connections = feeder[(feeder['Node_i'] == node) | (feeder['Node_j'] == node)].shape[0]
            if num_connections > max_connections:
                max_connections = num_connections
                slack_bus = node

        if slack_bus == -1:
             # This should theoretically not happen if nodes > 0
             logging.error("Could not automatically determine slack bus.")
             raise ValueError("Slack bus determination failed.")
        else:
             logging.info(f"Automatically selected Node {slack_bus} as slack bus (most connections: {max_connections})")

        # Convert 1-based slack_bus to 0-based index for matrix operations
        slack_idx = slack_bus - 1

        # Modify Bbus using the 0-based slack index
        # Sets up to "0" the columns of the slack node
        Bbus[:, slack_idx] = 0
        # Sets up to "0" the rows of the slack node
        Bbus[slack_idx, :] = 0
        # Sets up to "1" the diagonal element of the slack node
        Bbus[slack_idx, slack_idx] = 1

        X=np.linalg.inv(Bbus)
        X=abs(X)
        X[slack_idx,slack_idx]=0
        
        Al=A.T #Incidence Matrix Line-Node
        Bd=abs(Bp)

        #-------------Line-Node Sensibility Factors Calculation ---------------
        lodf=np.zeros((l,l))
        ptdf_rr=np.zeros((l,l))
        
        ptdf_rn=Bd@Al@X#PTDF Line-Node
        ptdf_rr=Bd@Al@X@Al.T#PTDF Line-Line
 
        #--LODFS CALCULATION--
        for s in range(l):
            for kl in range(l):
                if ptdf_rr[kl,kl]==1:
                    ptdf_rr[kl,kl]=0
                lodf[s,kl]=ptdf_rr[s,kl]*(1/(1-ptdf_rr[kl,kl]))
        
        np.fill_diagonal(lodf,0)    
            
        pt_dic[k]=ptdf_rn
        ld_dic[k]=lodf
        pt_nf,pt_nc=ptdf_rn.shape
        df_nf,df_nc=lodf.shape
        
        fil_pf=[]
        row_pf=list(feeder['Line'])
        for k1 in row_pf:
            fil_pf+=list(np.repeat(k1,pt_nc))
   
        fil_df=[]
        row_df=list(feeder['Line'])
        for k2 in row_df:
            fil_df+=list(np.repeat(k2,df_nc))

    for k in pt_dic:
        df_ptdic=pd.DataFrame(pt_dic[k])
    
    df_ptdic.set_index(feeder['Code_Line'],inplace=True, drop=True)
    df_ptdic.columns = config['zones']
    #df_ptdic.to_csv('PTDF.csv')

    #Calculation complete 
    print('The PTDF Matrix was calculated without errors')             
    return df_ptdic



def adjust_unit_capacity(SimData, u_idx, scaling=1, value=None, singleunit=False):
    """
    Function used to modify the installed capacities in the Dispa-SET generated input data
    The function update the Inputs.p file in the simulation directory at each call
    :param SimData:     Input data dictionary
    :param u_idx:       names of the units to be scaled
    :param scaling:     Scaling factor to be applied to the installed capacity
    :param value:       Absolute value of the desired capacity (! Applied only if scaling != 1 !)
    :param singleunit:  Set to true if the technology should remain lumped in a single unit
    :return:            New SimData dictionary
    """
    # a few checks:
    if len(u_idx) ==0:
        logging.warning('adjust_unit_capacity : list of units to be scaled is empty')
        return SimData
    if scaling > 1E10:
        logging.warning('adjust_unit_capacity: scaling factor is too high (' + str(scaling) + ')')
        return SimData

    # find the units to be scaled:
    units = SimData['units'].loc[u_idx,:]
    cond = SimData['units'].index.isin(u_idx)
    idx = pd.Series(np.where(cond)[0], index=units.index)
    TotalCapacity = (units.PowerCapacity * units.Nunits).sum()
    if scaling != 1:
        RequiredCapacity = TotalCapacity * scaling
    elif value is not None:
        RequiredCapacity = value
    else:
        RequiredCapacity = TotalCapacity
    if singleunit:
        Nunits_new = pd.Series(1, index=units.index)
    else:
        Nunits_new = (units.Nunits * RequiredCapacity / TotalCapacity).astype('float').round()
    Nunits_new[Nunits_new < 1] = 1
    Cap_new = units.PowerCapacity * RequiredCapacity / (units.PowerCapacity * Nunits_new).sum()
    for u in units.index:
        logging.info('Unit ' + u + ':')
        logging.info('    PowerCapacity: ' + str(SimData['units'].PowerCapacity[u]) + ' --> ' + str(Cap_new[u]))
        logging.info('    Nunits: ' + str(SimData['units'].Nunits[u]) + ' --> ' + str(Nunits_new[u]))
        factor = Cap_new[u] / SimData['units'].PowerCapacity[u]
        SimData['parameters']['PowerCapacity']['val'][idx[u]] = Cap_new[u]
        SimData['parameters']['Nunits']['val'][idx[u]] = Nunits_new[u]
        SimData['units'].loc[u, 'PowerCapacity'] = Cap_new[u]
        SimData['units'].loc[u, 'Nunits'] = Nunits_new[u]
        for col in ['CostStartUp', 'NoLoadCost', 'StorageCapacity', 'StorageChargingCapacity']:
            SimData['units'].loc[u, col] = SimData['units'].loc[u, col] * factor
        for param in ['CostShutDown', 'CostStartUp', 'PowerInitial', 'RampDownMaximum', 'RampShutDownMaximum',
                      'RampStartUpMaximum', 'RampUpMaximum', 'StorageCapacity']:
            SimData['parameters'][param]['val'][idx[u]] = SimData['parameters'][param]['val'][idx[u]] * factor

        for param in ['StorageChargingCapacity', 'StorageInitial']:
            # find index, if any:
            idx_s = np.where(np.array(SimData['sets']['s']) == u)[0]
            if len(idx_s) == 1:
                idx_s = idx_s[0]
                SimData['parameters'][param]['val'][idx_s] = SimData['parameters'][param]['val'][idx_s] * factor
    return SimData



def adjust_capacity(inputs, tech_fuel, scaling=1, value=None, singleunit=False, sto_fp_time_range=None, write_gdx=False, dest_path='', temp_path='Inputs.gdx'):
    """
    Function used to modify the installed capacities in the Dispa-SET generated input data
    The function update the Inputs.p file in the simulation directory at each call
    :param inputs:            Input data dictionary OR path to the simulation directory containing Inputs.p
    :param tech_fuel:         tuple with the technology and fuel type for which the capacity should be modified
    :param scaling:           Scaling factor to be applied to the installed capacity
    :param value:             Absolute value of the desired capacity (! Applied only if scaling != 1 !)
    :param singleunit:        Set to true if the technology should remain lumped in a single unit
    :param sto_fp_time_range: Ignored for non storage units. Specifies the range of full power time contained in the storage
                              the units have to fit in. If none, infinite range is considered
    :param temp_path:         The path to the location where to temporary build the Inputs.gdx
    :param write_gdx:         boolean defining if Inputs.gdx should be also overwritten with the new data
    :param dest_path:         Simulation environment path to write the new input data. If unspecified, no data is written!
    :return:                  New SimData dictionary
    """
    import pickle

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

    if not isinstance(tech_fuel, tuple):
        sys.exit('tech_fuel must be a tuple')


    # find the units to be scaled:
    cond = (SimData['units']['Technology'] == tech_fuel[0]) & (SimData['units']['Fuel'] == tech_fuel[1])

    # add condition to filter storage units outside the range of full power time
    if tech_fuel[0] in ["BATS", "BEVS", "CAES", "P2GS", "THMS"]:
        if sto_fp_time_range is None:
            sto_fp_time_range = (-np.infty, np.infty)
        elif not isinstance(sto_fp_time_range, tuple):
            logging.error("The range of full power time is invalid (not a tuple)")
            sys.exit(1)

        fp_time = 3600 * SimData["units"]["PowerCapacity"] / SimData["units"]["STOCapacity"]
        cond &= (sto_fp_time_range[0] <= fp_time) & (fp_time < sto_fp_time_range[1])

    u_idx = SimData['units'][cond].index.tolist()

    SimData = adjust_unit_capacity(SimData, u_idx, scaling=scaling, value=value, singleunit=singleunit)
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
            write_variables(SimData['config'], temp_path, [SimData['sets'], SimData['parameters']])
            shutil.copy(temp_path, os.path.join(dest_path, 'Inputs.gdx'))
            os.remove(temp_path)
    return SimData

def ranges_from_tresholds(thresholds, only_thresholds=False):
    """
    Outputs a list of ranges as a list of tuples from a list of thresholds.
    :param tresholds:        A list of tresholds to separate ranges at
    :param only_thresholds:  If false, add a bottom 0 threshold and ceil +np.infty bound, else restrict to given values
    """
    all = []
    if not only_thresholds:
        if thresholds[0] != 0:
            all.append(0)

        all += thresholds

        if all[-1] != np.inf:
            all.append(np.inf)
    else:
        all = thresholds

    res = []
    for i in range(len(all)-1):
        res.append((all[i], all[i+1]))

    # for x, y in zip(all[:-1], all[1:]):
    #     res.append((x, y))

    return res


def adjust_flexibility(inputs, flex_units, slow_units, flex_ratio, singleunit=False, write_gdx=False, dest_path=''):
    """
    Function used to modify the share of the flexible capacity in the Dispa-SET input data
    The function update the Inputs.p file in the simulation directory at each call

    :param inputs:      Input data dictionary OR path to the simulation directory containing Inputs.p
    :param flex_units:  Dispa-SET units table filtered with only the flexible ones
    :param slow_units:  Dispa-SET units table filtered with only the slow ones
    :param flex_ratio:  Target flexibility ratio (single number for all zones)
    :param singleunit:  Set to true if the technology should remain lumped in a single unit
    :param write_gdx:   boolean defining if Inputs.gdx should be also overwritten with the new data
    :param dest_path:   Simulation environment path to write the new input data. If unspecified, no data is written!
    :return:            New SimData dictionary
    """
    import pickle

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

    # find the units to be scaled:
    units = SimData['units']
    
    # current situation for all zones:"
    current_flex_cap = units.PowerCapacity[flex_units].sum() 
    current_total_cap = current_flex_cap + units.PowerCapacity[slow_units].sum()
    current_flex_ratio = current_flex_cap / current_total_cap
    
    #make new dataframe with the current country flex,slow and total installed capacities:
    zones = units.loc[flex_units.tolist()+slow_units.tolist(),:].Zone.unique().tolist()
    current = pd.DataFrame(index=zones,columns=['flex','slow','total','ratio'])
    for z in zones:
        current.flex[z] = units.loc[flex_units,:][units.loc[flex_units,:].Zone==z].PowerCapacity.sum()
        current.slow[z] = units.loc[slow_units,:][units.loc[slow_units,:].Zone==z].PowerCapacity.sum()
        current.total[z] = current.flex[z] + current.slow[z]
        current.ratio[z] = current.flex[z] / current.total[z]
        
    # target flexible capacity for all zones:"
    target_flex_cap = flex_ratio * current_total_cap
    
    # flexibile capacity to be added (positive) or removed (negative)
    delta_flex_cap = target_flex_cap - current_flex_cap
    
    if delta_flex_cap >0:
        # sort the current dataframe, highest flexibility first:
        current.sort_values('ratio',ascending=False,inplace=True)
        current['cum_sum'] = current.total.cumsum()    # save the cumulative zone capacities in a column
        # variable containing the remaining flexible capacity to be assigned to the zones:
        remaining = delta_flex_cap
        # Recursively add flexible capacity in each zone:
        for z in current.index:
            #weight of the current zone compared to the total of remaining zones
            weight = current.total[z]/(current_total_cap - current.cum_sum[z] + current.total[z])
            # added flexible capacity in this zone is bounded in order not to exceed the total capacity:
            added_flex_cap = min(weight*remaining,current.total[z] - current.flex[z])
            current.loc[z,'new_flex_cap'] = current.flex[z] + added_flex_cap
            current.loc[z,'new_slow_cap'] = current.slow[z] - added_flex_cap
            remaining -= added_flex_cap
    elif delta_flex_cap < 0:
        # sort the current dataframe, highest flexibility first:
        current.sort_values('ratio',ascending=True,inplace=True)
        current['cum_sum'] = current.total.cumsum()    # save the cumulative zone capacities in a column
        # variable containing the remaining flexible capacity to be assigned to the zones:
        remaining = -delta_flex_cap
        # Recursively add flexible capacity in each zone:
        for z in current.index:
            #weight of the current zone compared to the total of remaining zones
            weight = current.total[z]/(current_total_cap - current.cum_sum[z] + current.total[z])
            # added flexible capacity in this zone is bounded in order not to exceed the total capacity:
            removed_flex_cap = min(weight*remaining,current.flex[z])
            current.loc[z,'new_flex_cap'] = current.flex[z] - removed_flex_cap
            current.loc[z,'new_slow_cap'] = current.slow[z] + removed_flex_cap
            remaining -= removed_flex_cap
    else:
        current.loc[z,'new_flex_cap'] = current.flex
        current.loc[z,'new_slow_cap'] = current.slow
    del current['cum_sum']
    print(current)
    
    # last loop where units are actually scaled in each country:
    for z in zones:
        u_idx = units.loc[flex_units,:][units.loc[flex_units,:].Zone==z].index.tolist()
        SimData = adjust_unit_capacity(SimData, u_idx, scaling=current.loc[z,'new_flex_cap']/current.loc[z,'flex'], singleunit=singleunit)
        u_idx = units.loc[slow_units,:][units.loc[slow_units,:].Zone==z].index.tolist()
        SimData = adjust_unit_capacity(SimData, u_idx, scaling=current.loc[z,'new_slow_cap']/current.loc[z,'slow'], singleunit=singleunit)
    
    # Checking
    units_new = SimData['units']
    
    # current situation for all zones:"
    new_flex_cap = units_new.PowerCapacity[flex_units].sum() 
    new_total_cap = new_flex_cap + units_new.PowerCapacity[slow_units].sum()
    new_flex_ratio = new_flex_cap / new_total_cap  
    
    if (new_flex_ratio - flex_ratio) > 0.01:
        logging.error('the new flexbility ratio (' + str(new_flex_ratio) + ') is not equal to the desired one: ' + str(flex_ratio))
    

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
