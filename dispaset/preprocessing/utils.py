"""
This file gathers different functions used in the DispaSET pre-processing tools

@author: Sylvain Quoilin (sylvain.quoilin@ec.europa.eu)
"""

from __future__ import division

import logging
import os
import shutil
import sys

import numpy as np
import pandas as pd

from ..misc.gdx_handler import write_variables
from ..misc.str_handler import clean_strings, shrink_to_64


def EfficiencyTimeSeries(config,plants,Temperatures):
    '''
    Function that calculates an efficiency time series for each unit
    In case of generation unit, the efficiency is constant in time (for now)
    In case of of p2h units, the efficicncy is defined as the COP, which can be
    temperature-dependent or not
    If it is temperature-dependent, the formula is:
        COP = COP_nom + coef_a * (T-T_nom) + coef_b * (T-T_nom)^2
    
    :param plants:          Pandas dataframe with the original list of units
    :param Temperatures:    Dataframe with the temperature for all relevant units
    
    :returns:               Dataframe with a time series of the efficiency for each unit
    '''
    Efficiencies = pd.DataFrame(columns = plants.index,index=config['idx_long'])
    for u in plants.index:
        z = plants.loc[u,'Zone']
        if plants.loc[u,'Technology'] == 'P2HT' and 'Tnominal' in plants:
            eff = plants.loc[u,'COP'] + plants.loc[u,'coef_COP_a'] * (Temperatures[z] - plants.loc[u,'Tnominal'])
            + plants.loc[u,'coef_COP_a'] * (Temperatures[z] - plants.loc[u,'Tnominal'])**2
        elif plants.loc[u,'Technology'] == 'P2HT':
            eff = plants.loc[u,'COP']
        else:
            eff = plants.loc[u,'Efficiency']
        Efficiencies[u] = eff
    return Efficiencies

def select_units(units,config):
    '''
    Function returning a new list of units by removing the ones that have unknown
    technology, zero capacity, or unknown zone
    
    :param units:       Pandas dataframe with the original list of units
    :param config:      Dispa-SET config dictionnary
    :return:            New list of units
    '''
    for unit in units.index:
        if units.loc[unit,'Technology'] == 'Other':
            logging.warning('Removed Unit ' + str(units.loc[unit,'Unit']) + ' since its technology is unknown')
            units.drop(unit,inplace=True)
        elif units.loc[unit,'PowerCapacity'] == 0:
            logging.warning('Removed Unit ' + str(units.loc[unit,'Unit']) + ' since it has a null capacity')
            units.drop(unit,inplace=True)
        elif units.loc[unit,'Zone'] not in config['zones']:
            logging.warning('Removed Unit ' + str(units.loc[unit,'Unit']) + ' since its zone (' + str(units.loc[unit,'Zone'])+ ') is not in the list of zones')    
            units.drop(unit,inplace=True)
    units.index = range(len(units))
    return units

def incidence_matrix(sets, set_used, parameters, param_used):
    """
    This function generates the incidence matrix of the lines within the nodes
    A particular case is considered for the node "Rest Of the World", which is no explicitely defined in DispaSET
    """
    for i,l in enumerate(sets[set_used]):
        [from_node, to_node] = l.split('->')
        if (from_node.strip() in sets['n']) and (to_node.strip() in sets['n']):
            parameters[param_used]['val'][i, sets['n'].index(to_node.strip())] = 1
            parameters[param_used]['val'][i, sets['n'].index(from_node.strip())] = -1
        elif (from_node.strip() in sets['n']) and (to_node.strip() == 'RoW'):
            parameters[param_used]['val'][i, sets['n'].index(from_node.strip())] = -1
        elif (from_node.strip() == 'RoW') and (to_node.strip() in sets['n']):
            parameters[param_used]['val'][i, sets['n'].index(to_node.strip())] = 1
        else:
            logging.error("The line " + str(l) + " contains unrecognized nodes (" + from_node.strip() + ' or ' + to_node.strip() + ")")

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
    :param NTC:                 Day-ahead net transfer capacities (pd dataframe)
    :param Historical_flows:    Historical flows (pd dataframe)
    """
    index = NTC_inter.index.tz_localize(None).intersection(Historical_flows.index.tz_localize(None))
    if len(index)==0:
        logging.error('The two input dataframes (NTCs and Historical flows) must have the same index. No common values have been found')
        sys.exit(1)
    elif len(index) < len(NTC_inter) or len(index) < len(Historical_flows):
        diff = np.maximum(len(Historical_flows),len(NTC_inter)) - len(index)
        logging.warning('The two input dataframes (NTCs and Historical flows) do not share the same index, although some values are common. The intersection has been considered and ' + str(diff) + ' data points have been lost')
    # Checking that all values are positive:
    if (NTC_inter.values < 0).any():
        pos = np.where(NTC_inter.values < 0)
        logging.warning('WARNING: At least NTC value is negative, for example in line ' + str(NTC_inter.columns[pos[1][0]]) + ' and time step ' + str(NTC_inter.index[pos[0][0]]))
    if (Historical_flows.values < 0).any():
        pos = np.where(Historical_flows.values < 0)
        logging.warning('WARNING: At least one historical flow is negative, for example in line ' + str(Historical_flows.columns[pos[1][0]]) + ' and time step ' + str(Historical_flows.index[pos[0][0]]))
    all_connections = []
    simulation_connections = []
    # List all connections from the dataframe headers:
    ConList = Historical_flows.columns.tolist() + [x for x in NTC_inter.columns.tolist() if x not in Historical_flows.columns.tolist()]
    for connection in ConList:
        z = connection.split(' -> ')
        if z[0] in Simulation_list:
            all_connections.append(connection)
            if z[1] in Simulation_list:
                simulation_connections.append(connection)
        elif z[1] in Simulation_list:
            all_connections.append(connection)

    df_zones_simulated = pd.DataFrame(index=index)
    for interconnection in simulation_connections:
        if interconnection in NTC_inter.columns:
            df_zones_simulated[interconnection] = NTC_inter[interconnection]
            logging.info('Detected interconnection ' + interconnection + '. The historical NTCs will be imposed as maximum flow value')
    interconnections1 = df_zones_simulated.columns

    # Display a warning if a zone is isolated:
    for z in Simulation_list:
        if not any([z in conn for conn in interconnections1]) and len(Simulation_list)>1:
            logging.warning('Zone ' + z + ' does not appear to be connected to any other zone in the NTC table. It should be simulated in isolation')

    df_RoW_temp = pd.DataFrame(index=index)
    connNames = []
    for interconnection in all_connections:
        if interconnection in Historical_flows.columns and interconnection not in simulation_connections:
            df_RoW_temp[interconnection] = Historical_flows[interconnection]
            connNames.append(interconnection)

    compare_set = set()
    for k in connNames:
        if not k[0:2] in compare_set and k[0:2] in Simulation_list:
            compare_set.add(k[0:2])

    df_zones_RoW = pd.DataFrame(index=index)
    while compare_set:
        nameToCompare = compare_set.pop()
        exports = []
        imports = []
        for name in connNames:
            if nameToCompare[0:2] in name[0:2]:
                exports.append(connNames.index(name))
                logging.info('Detected interconnection ' + name + ', happening between a simulated zone and the rest of the world. The historical flows will be imposed to the model')
            elif nameToCompare[0:2] in name[6:8]:
                imports.append(connNames.index(name))
                logging.info('Detected interconnection ' + name + ', happening between the rest of the world and a simulated zone. The historical flows will be imposed to the model')

        flows_out = pd.concat(df_RoW_temp[connNames[exports[i]]] for i in range(len(exports)))
        flows_out = flows_out.groupby(flows_out.index).sum()
        flows_out.name = nameToCompare + ' -> RoW'
        df_zones_RoW[nameToCompare + ' -> RoW'] = flows_out
        flows_in = pd.concat(df_RoW_temp[connNames[imports[j]]] for j in range(len(imports)))
        flows_in = flows_in.groupby(flows_in.index).sum()
        flows_in.name = 'RoW -> ' + nameToCompare
        df_zones_RoW['RoW -> ' + nameToCompare] = flows_in
    interconnections2 = df_zones_RoW.columns
    inter = list(interconnections1) + list(interconnections2)
    return (df_zones_simulated, df_zones_RoW, inter)



def clustering(plants, method='Standard', Nslices=20, PartLoadMax=0.1, Pmax=30):
    """
    Merge excessively disaggregated power Units.

    :param plants:          Pandas dataframe with each power plant and their characteristics (following the DispaSET format)
    :param method:          Select clustering method ('Standard'/'LP'/None)
    :param Nslices:         Number of slices used to fingerprint each power plant characteristics. slices in the power plant data to categorize them  (fewer slices involves that the plants will be aggregated more easily)
    :param PartLoadMax:     Maximum part-load capability for the unit to be clustered
    :param Pmax:            Maximum power for the unit to be clustered
    :return:                A list with the merged plants and the mapping between the original and merged units
    """

    # Checking the the required columns are present in the input pandas dataframe:
    required_inputs = ['Unit', 'PowerCapacity', 'PartLoadMin', 'RampUpRate', 'RampDownRate', 'StartUpTime',
                       'MinUpTime', 'MinDownTime', 'NoLoadCost', 'StartUpCost', 'Efficiency']
    for input_value in required_inputs:
        if input_value not in plants.columns:
            logging.error("The plants dataframe requires a '" + input_value + "' column for clustering")
            sys.exit(1)
    if not "Nunits" in plants:
        plants['Nunits'] = 1

    # Checking the validity of the selected clustering method
    OnlyOnes = (plants['Nunits'] == 1).all()
    if method in ['Standard','MILP']:
        if not OnlyOnes:
            logging.warning("The standard (or MILP) clustering method is only applicable if all values of the Nunits column in the power plant data are set to one. At least one different value has been encountered. No clustering will be applied")
    elif method == 'LP clustered':
        if not OnlyOnes:
            logging.warning("The LP clustering method aggregates all the units of the same type. Individual units are not considered")
            # Modifying the table to remove multiple-units plants:
            for key in ['PowerCapacity', 'STOCapacity', 'STOMaxChargingPower','InitialPower','CHPMaxHeat']:
                if key in plants:
                    plants.loc[:,key] = plants.loc[:,'Nunits'] * plants.loc[:,key]
            plants['Nunits'] = 1
            OnlyOnes = True
    elif method == 'LP':
        pass
    elif method == 'Integer clustering':
        pass
    elif method == 'No clustering':
        pass
    else:
        logging.error('Method argument ("' + str(method) + '") not recognized in the clustering function')
        sys.exit(1)

    # Number of units:
    Nunits = len(plants)
    plants.index = range(Nunits)

    # Definition of the mapping variable, from the old power plant list the new (merged) one:
    map_old_new = np.zeros(Nunits)
    map_plant_orig = []

    # Slicing:
    bounds = {'PartLoadMin': np.linspace(0, 1, Nslices), 'RampUpRate': np.linspace(0, 1, Nslices),
              'RampDownRate': np.linspace(0, 1, Nslices), 'StartUpTime': _mylogspace(0, 36, Nslices),
              'MinUpTime': _mylogspace(0, 168, Nslices), 'MinDownTime': _mylogspace(0, 168, Nslices),
              'NoLoadCost': np.linspace(0, 50, Nslices), 'StartUpCost': np.linspace(0, 500, Nslices),
              'Efficiency': np.linspace(0, 1, Nslices)}

    # Definition of the fingerprint value of each power plant, i.e. the pattern of the slices number in which each of
    # its characteristics falls:
    fingerprints = []
    fingerprints_merged = []
    for i in plants.index:
        fingerprints.append([_find_nearest(bounds['PartLoadMin'], plants['PartLoadMin'][i]),
                             _find_nearest(bounds['RampUpRate'], plants['RampUpRate'][i]),
                             _find_nearest(bounds['RampDownRate'], plants['RampDownRate'][i]),
                             _find_nearest(bounds['StartUpTime'], plants['StartUpTime'][i]),
                             _find_nearest(bounds['MinUpTime'], plants['MinUpTime'][i]),
                             _find_nearest(bounds['MinDownTime'], plants['MinDownTime'][i]),
                             _find_nearest(bounds['NoLoadCost'], plants['NoLoadCost'][i]),
                             _find_nearest(bounds['StartUpCost'], plants['StartUpCost'][i]),
                             _find_nearest(bounds['Efficiency'], plants['Efficiency'][i])])

    # Definition of the merged power plants dataframe:
    plants_merged = pd.DataFrame(columns=plants.columns)

    # Find the columns containing string values (in addition to "Unit")
    #    string_keys = []
    #    for i in range(len(plants.columns)):
    #        if plants.columns[i] != 'Unit' and plants.dtypes[i] == np.dtype('O'):
    #            string_keys.append(plants.columns[i])
    string_keys = ['Zone', 'Technology', 'Fuel','CHPType']
    # First, fill nan values:
    for key in string_keys:
        plants[key].fillna('',inplace=True)

    for i in plants.index:  # i is the plant to be added to the new list
        merged = False
        plants_string = plants[string_keys].iloc[i].fillna('')
        for j in plants_merged.index:  # j corresponds to the clustered plants
            same_type = all(plants_string == plants_merged[string_keys].iloc[j].fillna(''))
            same_fingerprint = (fingerprints[i] == fingerprints_merged[j])
            low_pmin = (plants['PartLoadMin'][i] <= PartLoadMax)
            low_pmax = (plants['PowerCapacity'][i] <= Pmax)
            highly_flexible = plants['RampUpRate'][i] > 1 / 60 and (plants['RampDownRate'][i] > 1 / 60) and (
            plants['StartUpTime'][i] < 1) and (plants['MinDownTime'][i] <= 1) and (plants['MinUpTime'][i] <= 1)
            cluster = OnlyOnes and same_type and ((same_fingerprint and low_pmin) or highly_flexible or low_pmax)
            if method in ('Standard','MILP') and cluster:  # merge the two plants in plants_merged:
                P_old = plants_merged['PowerCapacity'][j]  # Old power in plants_merged
                P_add = plants['PowerCapacity'][i]  # Additional power to be added
                for key in plants_merged:
                    if key in ['RampUpRate', 'RampDownRate', 'MinUpTime', 'MinDownTime', 'NoLoadCost', 'Efficiency',
                               'MinEfficiency', 'STOChargingEfficiency', 'CO2Intensity', 'STOSelfDischarge',
                               'CHPPowerToHeat','CHPPowerLossFactor','COP','TNominal','coef_COP_a','coef_COP_b']:
                        # Do a weighted average:
                        plants_merged.loc[j, key] = (plants_merged[key][j] * P_old + plants[key][i] * P_add) / (
                        P_add + P_old)
                    elif key in ['PowerCapacity', 'STOCapacity', 'STOMaxChargingPower','InitialPower','CHPMaxHeat']:
                        # Do a sum:
                        plants_merged.loc[j, key] = plants_merged[key][j] + plants[key][i]
                    elif key in ['PartLoadMin', 'StartUpTime']:
                        # Take the minimum
                        plants_merged.loc[j, key] = np.minimum(plants_merged[key][j] * P_old,
                                                               plants[key][i] * P_add) / (P_add + P_old)
                    elif key == 'RampingCost':
                        # The starting cost must be added to the ramping cost
                        Cost_to_fullload = P_add * (1 - plants['PartLoadMin'][i]) * plants['RampingCost'][i] + \
                                           plants['StartUpCost'][i]
                        plants_merged.loc[j, key] = (P_old * plants_merged[key][j] + Cost_to_fullload) / (P_old + P_add)
                    elif key == 'Nunits':
                        plants_merged.loc[j, key] = 1
                map_old_new[i] = j
                map_plant_orig[j].append(i)
                merged = True
                break
            elif method == 'LP clustered' and same_type and OnlyOnes:
                P_old = plants_merged['PowerCapacity'][j]  # Old power in plants_merged
                P_add = plants['PowerCapacity'][i]  # Additional power to be added
                for key in plants_merged:
                    if key in ['RampUpRate', 'RampDownRate', 'MinUpTime', 'MinDownTime', 'NoLoadCost', 'Efficiency',
                               'MinEfficiency', 'STOChargingEfficiency', 'CO2Intensity', 'STOSelfDischarge']:
                        # Do a weighted average:
                        plants_merged.loc[j, key] = (plants_merged[key][j] * P_old + plants[key][i] * P_add) / (
                        P_add + P_old)
                    elif key in ['PowerCapacity', 'STOCapacity', 'STOMaxChargingPower','InitialPower','CHPMaxHeat']:
                        # Do a sum:
                        plants_merged.loc[j, key] = plants_merged[key][j] + plants[key][i]
                    elif key in ['PartLoadMin', 'StartUpTime']:
                        # impose 0
                        plants_merged.loc[j, key] = 0
                    elif key == 'RampingCost':
                        # The starting cost must be added to the ramping cost
                        Cost_to_fullload = P_add * (1 - plants['PartLoadMin'][i]) * plants['RampingCost'][i] + \
                                           plants['StartUpCost'][i]
                        plants_merged.loc[j, key] = (P_old * plants_merged[key][j] + Cost_to_fullload) / (P_old + P_add)
                    elif key == 'Nunits':
                        plants_merged.loc[j, key] = 1
                map_old_new[i] = j
                map_plant_orig[j].append(i)
                merged = True
                break
            elif method == 'Integer clustering' and same_type:
                for key in plants_merged:
                    if key in ['PowerCapacity','RampUpRate', 'RampDownRate', 'MinUpTime', 'MinDownTime', 'NoLoadCost', 'Efficiency',
                               'MinEfficiency', 'STOChargingEfficiency', 'CO2Intensity', 'STOSelfDischarge',
                               'STOCapacity', 'STOMaxChargingPower','InitialPower','PartLoadMin', 'StartUpTime','RampingCost',
                               'CHPPowerToHeat','CHPPowerLossFactor','CHPMaxHeat']:
                        # Do a weighted average:
                        plants_merged.loc[j, key] = (plants_merged.loc[j,key] * plants_merged.loc[j,'Nunits'] + plants.loc[i,key] * plants.loc[i,'Nunits']) / (plants_merged.loc[j,'Nunits'] + plants.loc[i,'Nunits'])
                plants_merged.loc[j, 'Nunits'] = plants_merged.loc[j,'Nunits'] + plants.loc[i,'Nunits']
                        
                map_old_new[i] = j
                map_plant_orig[j].append(i)
                merged = True
                break

        if not merged:  # Add a new plant in plants_merged:
            plants_merged = plants_merged.append(plants.loc[i], ignore_index=True)
            plants_merged = plants_merged.copy()
            map_plant_orig.append([i])
            map_old_new[i] = len(map_plant_orig) - 1
            fingerprints_merged.append(fingerprints[i])

    Nunits_merged = len(plants_merged)
    mapping = {'NewIndex': {}, 'FormerIndexes': {}}
    #    mapping['NewIdx'] = map_plant_orig
    #    mapping['OldIdx'] = map_old_new
    # Modify the Unit names with the original index number. In case of merged plants, indicate all indexes + the plant type and fuel
    for j in range(Nunits_merged):
        if len(map_plant_orig[j]) == 1:  # The plant has not been merged
            NewName = str(map_plant_orig[j]) + ' - ' + plants_merged['Unit'][j]
            NewName = shrink_to_64(clean_strings(NewName))
            NewName = NewName.rstrip()                          # remove space at the end because it is not considered by gams
            plants_merged.loc[j, 'Unit'] = NewName
            mapping['FormerIndexes'][NewName] = [map_plant_orig[j][0]]
            mapping['NewIndex'][map_plant_orig[j][0]] = NewName
        else:
            all_stringkeys = ''
            for key in string_keys:
                all_stringkeys = all_stringkeys + ' - ' + plants_merged[key][j]
            NewName = str(map_plant_orig[j]) + all_stringkeys
            NewName = shrink_to_64(clean_strings(NewName))
            NewName = NewName.rstrip()                          # remove space at the end because it is not considered by gams
            plants_merged.loc[j, 'Unit'] = NewName
            list_oldplants = [x for x in map_plant_orig[j]]
            mapping['FormerIndexes'][NewName] = list_oldplants
            for oldplant in list_oldplants:
                mapping['NewIndex'][oldplant] = NewName

    # Transforming the start-up cost into ramping for the plants that did not go through any clustering:
    if method == 'LP clustered':
        for i in range(Nunits_merged):
            if plants_merged['RampingCost'][i] == 0:
                Power = plants_merged['PowerCapacity'][i]
                Start_up = plants_merged['StartUpCost'][i]
                plants_merged.loc[i, 'RampingCost'] = Start_up / Power
                
    # Correcting the Nunits field of the clustered plants (must be integer):
    elif method == 'Integer clustering':
        for idx in plants_merged.index:
            N = np.round(plants_merged.loc[idx,'Nunits'])
        for key in ['PowerCapacity', 'STOCapacity', 'STOMaxChargingPower','InitialPower','NoLoadCost']:
            if key in plants_merged.columns:
                plants_merged.loc[idx,key] = plants_merged.loc[idx,key] * N / plants_merged.loc[idx,'Nunits']
        plants_merged.loc[idx,'Nunits'] = N
                
    # Updating the index of the merged plants dataframe with the new unit names, after some cleaning:
    plants_merged.index = plants_merged['Unit']

    if Nunits != len(plants_merged):
        logging.info('Clustered ' + str(Nunits) + ' original units into ' + str(len(plants_merged)) + ' new units')
    else:
        logging.warning('Did not cluster any unit')
    return plants_merged, mapping

## Helpers

def _mylogspace(low, high, N):
    """
    Self-defined logspace function in which low and high are the first and last values of the space
    """
    # shifting all values so that low = 1
    space = np.logspace(0, np.log10(high + low + 1), N) - (low + 1)
    return (space)


def _find_nearest(array, value):
    """
    Self-defined function to find the index of the nearest value in a vector
    """
    idx = (np.abs(array - value)).argmin()
    return idx


def adjust_storage(inputs,tech_fuel,scaling=1,value=None,write_gdx=False,dest_path=''):
    '''
    Function used to modify the storage capacities in the Dispa-SET generated input data
    The function update the Inputs.p file in the simulation directory at each call

    :param inputs:      Input data dictionary OR path to the simulation directory containing Inputs.p
    :param tech_fuel:   tuple with the technology and fuel type for which the capacity should be modified
    :param scaling:     Scaling factor to be applied to the installed capacity
    :param value:       Absolute value of the desired capacity (! Applied only if scaling != 1 !)
    :param write_gdx:   boolean defining if Inputs.gdx should be also overwritten with the new data
    :param dest_path:   Simulation environment path to write the new input data. If unspecified, no data is written!
    :return:            New SimData dictionary
    '''
    import pickle

    if isinstance(inputs,str):
        path = inputs
        inputfile = path + '/Inputs.p'
        if not os.path.exists(path):
            sys.exit('Path + "' + path + '" not found')
        with open(inputfile, 'rb') as f:
            SimData = pickle.load(f)
    elif isinstance(inputs,dict):
        SimData = inputs
    else:
        logging.error('The input data must be either a dictionary or string containing a valid directory')
        sys.exit(1)

    if not isinstance(tech_fuel,tuple):
        sys.exit('tech_fuel must be a tuple')

    # find the units to be scaled:
    cond = (SimData['units']['Technology'] == tech_fuel[0]) & (SimData['units']['Fuel'] == tech_fuel[1]) & (SimData['units']['StorageCapacity'] > 0)
    units = SimData['units'][cond]
    idx = pd.Series(np.where(cond)[0],index=units.index)
    TotalCapacity = (units.StorageCapacity*units.Nunits).sum()
    if scaling != 1:
        RequiredCapacity = TotalCapacity*scaling
    elif value is not None:
        RequiredCapacity = value
    else:
        RequiredCapacity = TotalCapacity
    factor = RequiredCapacity/TotalCapacity
    for u in units.index:
        logging.info('Unit ' + u +':')
        logging.info('    StorageCapacity: ' + str(SimData['units'].StorageCapacity[u]) + ' --> ' + str(SimData['units'].StorageCapacity[u]*factor))
        SimData['units'].loc[u,'StorageCapacity'] = SimData['units'].loc[u,'StorageCapacity']*factor
        SimData['parameters']['StorageCapacity']['val'][idx[u]] = SimData['parameters']['StorageCapacity']['val'][idx[u]]*factor

    if dest_path == '':
        logging.info('Not writing any input data to the disk')
    else:
        if not os.path.isdir(dest_path):
            shutil.copytree(path, dest_path)
            logging.info('Created simulation environment directory ' + dest_path)
        logging.info('Writing input files to ' + dest_path)
        import cPickle
        with open(os.path.join(dest_path, 'Inputs.p'), 'wb') as pfile:
            cPickle.dump(SimData, pfile, protocol=cPickle.HIGHEST_PROTOCOL)
        if write_gdx:
            write_variables(SimData['config'], 'Inputs.gdx', [SimData['sets'], SimData['parameters']])
            shutil.copy('Inputs.gdx', dest_path + '/')
            os.remove('Inputs.gdx')
    return SimData


def adjust_capacity(inputs,tech_fuel,scaling=1,value=None,singleunit=False,write_gdx=False,dest_path=''):
    '''
    Function used to modify the installed capacities in the Dispa-SET generated input data
    The function update the Inputs.p file in the simulation directory at each call

    :param inputs:      Input data dictionary OR path to the simulation directory containing Inputs.p
    :param tech_fuel:   tuple with the technology and fuel type for which the capacity should be modified
    :param scaling:     Scaling factor to be applied to the installed capacity
    :param value:       Absolute value of the desired capacity (! Applied only if scaling != 1 !)
    :param singleunit:  Set to true if the technology should remain lumped in a single unit
    :param write_gdx:   boolean defining if Inputs.gdx should be also overwritten with the new data
    :param dest_path:   Simulation environment path to write the new input data. If unspecified, no data is written!
    :return:            New SimData dictionary
    '''
    import pickle

    if isinstance(inputs, str):
        path = inputs
        inputfile = path + '/Inputs.p'
        if not os.path.exists(path):
            sys.exit('Path + "' + path + '" not found')
        with open(inputfile, 'rb') as f:
            SimData = pickle.load(f)
    elif isinstance(inputs,dict):
        SimData = inputs
        path = SimData['config']['SimulationDirectory']
    else:
        logging.error('The input data must be either a dictionary or string containing a valid directory')
        sys.exit(1)

    if not isinstance(tech_fuel,tuple):
        sys.exit('tech_fuel must be a tuple')

    # find the units to be scaled:
    cond = (SimData['units']['Technology'] == tech_fuel[0]) & (SimData['units']['Fuel'] == tech_fuel[1])
    units = SimData['units'][cond]
    idx = pd.Series(np.where(cond)[0],index=units.index)
    TotalCapacity = (units.PowerCapacity*units.Nunits).sum()
    if scaling != 1:
        RequiredCapacity = TotalCapacity*scaling
    elif value is not None:
        RequiredCapacity = value
    else:
        RequiredCapacity = TotalCapacity
    if singleunit:
        Nunits_new = pd.Series(1,index=units.index)
    else:
        Nunits_new = (units.Nunits * RequiredCapacity/TotalCapacity).round()
    Nunits_new[Nunits_new < 1] = 1
    Cap_new = units.PowerCapacity * RequiredCapacity/(units.PowerCapacity*Nunits_new).sum()
    for u in units.index:
        logging.info('Unit ' + u +':')
        logging.info('    PowerCapacity: ' + str(SimData['units'].PowerCapacity[u]) + ' --> ' + str(Cap_new[u]))
        logging.info('    Nunits: ' + str(SimData['units'].Nunits[u]) + ' --> ' + str(Nunits_new[u]))
        factor = Cap_new[u]/SimData['units'].PowerCapacity[u]
        SimData['parameters']['PowerCapacity']['val'][idx[u]] = Cap_new[u]
        SimData['parameters']['Nunits']['val'][idx[u]] = Nunits_new[u]
        SimData['units'].loc[u,'PowerCapacity'] = Cap_new[u]
        SimData['units'].loc[u,'Nunits'] = Nunits_new[u]
        for col in ['CostStartUp', 'NoLoadCost','StorageCapacity','StorageChargingCapacity']:
            SimData['units'].loc[u,col] = SimData['units'].loc[u,col] * factor
        for param in ['CostShutDown','CostStartUp','PowerInitial','RampDownMaximum','RampShutDownMaximum','RampStartUpMaximum','RampUpMaximum','StorageCapacity']:
            SimData['parameters'][param]['val'][idx[u]] = SimData['parameters'][param]['val'][idx[u]]*factor
        for param in ['StorageChargingCapacity']:
            # find index, if any:
            idx_s = np.where(np.array(SimData['sets']['s']) == u)[0]
            if len(idx_s) == 1:
                idx_s = idx_s[0]
                SimData['parameters'][param]['val'][idx_s] = SimData['parameters'][param]['val'][idx_s]*factor
    if dest_path == '':
        logging.info('Not writing any input data to the disk')
    else:
        if not os.path.isdir(dest_path):
            shutil.copytree(path,dest_path)
            logging.info('Created simulation environment directory ' + dest_path)
        logging.info('Writing input files to ' + dest_path)
        with open(os.path.join(dest_path, 'Inputs.p'), 'wb') as pfile:
            pickle.dump(SimData, pfile, protocol=pickle.HIGHEST_PROTOCOL)
        if write_gdx:
            write_variables(SimData['config'], 'Inputs.gdx', [SimData['sets'], SimData['parameters']])
            shutil.copy('Inputs.gdx', dest_path + '/')
            os.remove('Inputs.gdx')
    return SimData