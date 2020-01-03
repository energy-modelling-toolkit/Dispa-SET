# -*- coding: utf-8 -*-
"""
This is the main file of the DispaSET pre-processing tool. It comprises a single function that generated the DispaSET simulation environment.

@author: S. Quoilin, edited by M. Zech
"""

#todo check indices

import datetime as dt
import logging
import os
import shutil
import sys

import numpy as np
import pandas as pd
try:
    from future.builtins import int
except ImportError:
    logging.warning("Couldn't import future package. Numeric operations may differ among different versions due to incompatible variable types")
    pass

from .data_check import check_units, check_chp, check_sto, check_heat_demand, check_df, isStorage, check_MinMaxFlows,check_AvailabilityFactors, check_clustering
from .utils import clustering, interconnections, incidence_matrix
from .data_handler import UnitBasedTable,NodeBasedTable,merge_series, \
        define_parameter, write_to_excel, load_csv, load_config_excel, load_config_yaml
from ..misc.gdx_handler import write_variables
from ..common import commons  # Load fuel types, technologies, timestep, etc:
from .data_loader import DataLoader


GMS_FOLDER = os.path.join(os.path.dirname(__file__), '..', 'GAMS')

def build_simulation(config):
    '''The main function for building the optimization problem'''

    _define_default_values(config)
    config['SimulationType'] == 'LP' or config['SimulationType'] == 'LP clustered'
    data = DataLoader(config)

    idx_utc, idx_utc_noloc, idx_utc_year_noloc = get_indices(config)
    enddate_long = idx_utc_noloc[-1] + dt.timedelta(days=config['LookAhead'])
    idx_long = pd.DatetimeIndex(pd.date_range(start=idx_utc_noloc[0], end=enddate_long, freq=commons['TimeStep']))
    sets = load_sets(
        Nhours_long = len(idx_long),
        look_ahead = data.config['LookAhead'],
        plants_index = data.Plants_merged.index.tolist(),
        plants_sto_index = data.plants_sto.index.tolist(),
        Plants_chp_index = data.plants_chp.index.tolist(),
        countries = data.config['countries'],
        Interconnections = data.Interconnections,
        plants_uc = data.plants_expanded.index.tolist()
    )

    sets_param = load_params()  # the parameters with their formal structure without data
    time_range = (idx_utc[-1] - idx_utc[0]).days
    sim = config['SimulationDirectory']
    dispa_version = str(get_git_revision_tag())
    SimData = {'sets': sets, 'parameters': sets_param, 'config': config, 'units': data.Plants_merged, 'version': dispa_version}

    parameters = build_model_parameters(config, sets, sets_param, data, idx_long)
    build_sim_dir(config, sim, sets, parameters, SimData)
    return SimData


def get_git_revision_tag():
    '''Get version of DispaSET used for this run. tag + commit hash'''
    from subprocess import check_output

    try:
        return check_output(["git", "describe", "--tags", "--always"]).strip()
    except:
        return 'NA'


def get_indices(config):
    '''Gets all relevant indices for the model simulation'''
    # Indexes of the simulation:
    idx_std = pd.DatetimeIndex(pd.date_range(start=pd.datetime(*config['StartDate']),
                                            end=pd.datetime(*config['StopDate']),
                                            freq=commons['TimeStep'])
                            ) #todo check brackets on master

    idx_utc_noloc = idx_std - dt.timedelta(hours=1)
    idx_utc = idx_utc_noloc.tz_localize('UTC')
    # Indexes for the whole year considered in StartDate
    idx_utc_year_noloc = pd.DatetimeIndex(pd.date_range(start=pd.datetime(*(config['StartDate'][0],1,1,0,0)),
                                                        end=pd.datetime(*(config['StartDate'][0],12,31,23,59,59)),
                                                        freq=commons['TimeStep'])
                                        )
    return idx_utc, idx_utc_noloc, idx_utc_year_noloc

def _define_default_values(config):
    if not isinstance(config['default']['CostLoadShedding'], (float, int)):
        config['default']['CostLoadShedding'] = 1000
    if not isinstance(config['default']['CostHeatSlack'], (float, int)):
        config['default']['CostHeatSlack'] = 50


def build_model_parameters(config, sets, sets_param, data, idx_long):
    '''Prepares the model parameters based on indices, and the data''' # todo better naming 
    parameters = dict()
    #plants = self.data.plants
    Plants_merged = data.Plants_merged
    Plants_sto = data.plants_sto
    Plants_chp = data.plants_chp
    ReservoirLevels = data.ReservoirLevels

    # Define all the parameters and set a default value of zero:
    for var in sets_param:
        parameters[var] = define_parameter(sets_param[var], sets, value=0)
    
    for var in ["Investment", "EconomicLifetime"]:
        parameters[var] = define_parameter(sets_param[var], sets, value=0)
        expanded_plants = data.plants_expanded.index
        parameters[var]["val"] =  Plants_merged.loc[expanded_plants, var].values
        
        # Plants_merged['FixedCost'] = pd.merge(Plants_merged, self.data.all_cost, how='left', on=['Fuel', 'Technology'])['FixedCost'].values
    for var in ["CostFixed"]:
        sets_param[var] = ['u']
        parameters[var] = define_parameter(sets_param[var], sets, value=0)
        parameters[var]["val"] = Plants_merged['FixedCost'].values

    Nunits = len(Plants_merged)

    # List of parameters whose default value is 1
    for var in ['AvailabilityFactor', 'Efficiency', 'Curtailment', 'StorageChargingEfficiency',
                'StorageDischargeEfficiency', 'Nunits']:
        parameters[var] = define_parameter(sets_param[var], sets, value=1)

    # List of parameters whose default value is very high
    for var in ['RampUpMaximum', 'RampDownMaximum', 'RampStartUpMaximum', 'RampShutDownMaximum',
                'EmissionMaximum']:
        parameters[var] = define_parameter(sets_param[var], sets, value=1e7)

    # Boolean parameters:
    for var in ['Technology', 'Fuel', 'Reserve', 'Location']:
        parameters[var] = define_parameter(sets_param[var], sets, value='bool')

    # List of parameters whose value is known, and provided in the dataframe Plants_merged.
    for var in ['Efficiency', 'PowerCapacity', 'PartLoadMin', 'TimeUpMinimum', 'TimeDownMinimum', 'CostStartUp',
                'CostRampUp','StorageCapacity', 'StorageSelfDischarge']:
        parameters[var]['val'] = Plants_merged[var].values

    # List of parameters whose value is not necessarily specified in the dataframe Plants_merged
    for var in ['Nunits']:
        if var in Plants_merged:
            parameters[var]['val'] = Plants_merged[var].values


    # List of parameters whose value is known, and provided in the dataframe Plants_sto.
    for var in ['StorageChargingCapacity', 'StorageChargingEfficiency']:
        parameters[var]['val'] = Plants_sto[var].values

    # The storage discharge efficiency is actually given by the unit efficiency:
    parameters['StorageDischargeEfficiency']['val'] = Plants_sto['Efficiency'].values

    # List of parameters whose value is known, and provided in the dataframe Plants_chp
    for var in ['CHPPowerToHeat','CHPPowerLossFactor', 'CHPMaxHeat']:
        parameters[var]['val'] = Plants_chp[var].values

    # Storage profile and initial state:
    for i, s in enumerate(sets['s']):
        if s in ReservoirLevels:
            # get the time
            parameters['StorageInitial']['val'][i] = ReservoirLevels[s][idx_long[0]] * \
                                                    Plants_sto['StorageCapacity'][s] * Plants_sto['Nunits'][s]
            parameters['StorageProfile']['val'][i, :] = ReservoirLevels[s][idx_long].values
            if any(ReservoirLevels[s] > 1):
                logging.warning(s + ': The reservoir level is sometimes higher than its capacity!')
        else:
            logging.warning('Could not find reservoir level data for storage plant ' + s + '. Assuming 50% of capacity')
            parameters['StorageInitial']['val'][i] = 0.5 * Plants_sto['StorageCapacity'][s]
            parameters['StorageProfile']['val'][i, :] = 0.5

    # Storage Inflows:
    for i, s in enumerate(sets['s']):
        if s in data.ReservoirScaledInflows:
            parameters['StorageInflow']['val'][i, :] = data.ReservoirScaledInflows[s][idx_long].values * \
                                                    Plants_sto['PowerCapacity'][s]
    # CHP time series:
    for i, u in enumerate(sets['chp']):
        if u in data.HeatDemand:
            parameters['HeatDemand']['val'][i, :] = data.HeatDemand[u][idx_long].values
            parameters['CostHeatSlack']['val'][i, :] = data.CostHeatSlack[u][idx_long].values

    # Ramping rates are reconstructed for the non dimensional value provided (start-up and normal ramping are not differentiated)
    parameters['RampUpMaximum']['val'] = Plants_merged['RampUpRate'].values * Plants_merged['PowerCapacity'].values * 60
    parameters['RampDownMaximum']['val'] = Plants_merged['RampDownRate'].values * Plants_merged[
        'PowerCapacity'].values * 60
    parameters['RampStartUpMaximum']['val'] = Plants_merged['RampUpRate'].values * Plants_merged[
        'PowerCapacity'].values * 60
    parameters['RampShutDownMaximum']['val'] = Plants_merged['RampDownRate'].values * Plants_merged[
        'PowerCapacity'].values * 60

    # If Curtailment is not allowed, set to 0:
    if config['AllowCurtailment'] == 0:
        parameters['Curtailment'] = define_parameter(sets_param['Curtailment'], sets, value=0)

    # Availability Factors
    if len(data.AF.columns) != 0:
        for i, u in enumerate(sets['u']):
            if u in data.AF.columns:
                parameters['AvailabilityFactor']['val'][i, :] = data.AF[u].values

    # Demand
    # Dayahead['NL'][1800:1896] = Dayahead['NL'][1632:1728]
    reserve_2U_tot = {i: (np.sqrt(10 * data.PeakLoad[i] + 150 ** 2) - 150) for i in data.Load.columns}
    reserve_2D_tot = {i: (0.5 * reserve_2U_tot[i]) for i in data.Load.columns}

    values = np.ndarray([len(sets['mk']), len(sets['n']), len(sets['h'])])
    for i in range(len(sets['n'])):
        values[0, i, :] = data.Load[sets['n'][i]]
        values[1, i, :] = reserve_2U_tot[sets['n'][i]]
        values[2, i, :] = reserve_2D_tot[sets['n'][i]]

    parameters['Demand'] = {'sets': sets_param['Demand'], 'val': values}
    # Emission Rate:
    parameters['EmissionRate']['val'][:, 0] = Plants_merged['EmissionRate'].values

    # Load Shedding:
    for i, c in enumerate(sets['n']):
        parameters['LoadShedding']['val'][i] = data.LoadShedding[c] * data.PeakLoad[c]
        parameters['CostLoadShedding']['val'][i] = data.CostLoadShedding[c]

    # %%#################################################################################################################################################################################################
    # Variable Cost
    # Equivalence dictionary between fuel types and price entries in the config sheet:
    FuelEntries = {'BIO':'PriceOfBiomass', 'GAS':'PriceOfGas', 'HRD':'PriceOfBlackCoal', 'LIG':'PriceOfLignite', 'NUC':'PriceOfNuclear', 'OIL':'PriceOfFuelOil', 'PEA':'PriceOfPeat'}
    for unit in range(Nunits):
        found = False
        for FuelEntry in FuelEntries:
            if Plants_merged['Fuel'][unit] == FuelEntry:
                parameters['CostVariable']['val'][unit, :] = data.FuelPrices[FuelEntries[FuelEntry]] / Plants_merged['Efficiency'][unit] + \
                                                            Plants_merged['EmissionRate'][unit] * data.FuelPrices['PriceOfCO2']
                found = True
        # Special case for biomass plants, which are not included in EU ETS:
        if Plants_merged['Fuel'][unit] == 'BIO':
            parameters['CostVariable']['val'][unit, :] = data.FuelPrices['PriceOfBiomass'] / Plants_merged['Efficiency'][
                unit]
            found = True
        if not found:
            logging.warning('No fuel price value has been found for fuel ' + Plants_merged['Fuel'][unit] + ' in unit ' + \
                Plants_merged['Unit'][unit] + '. A null variable cost has been assigned')

    # %%#################################################################################################################################################################################################

    # Maximum Line Capacity
    for i, l in enumerate(sets['l']):
        if l in data.NTCs.columns:
            parameters['FlowMaximum']['val'][i, :] = data.NTCs[l]
        if l in data.Inter_RoW.columns:
            parameters['FlowMaximum']['val'][i, :] = data.Inter_RoW[l]
            parameters['FlowMinimum']['val'][i, :] = data.Inter_RoW[l]
    # Check values:
    check_MinMaxFlows(parameters['FlowMinimum']['val'], parameters['FlowMaximum']['val'])
    parameters['LineNode'] = incidence_matrix(sets, 'l', parameters, 'LineNode')

    # Outage Factors
    if len(data.Outages.columns) != 0:
        for i, u in enumerate(sets['u']):
            if u in data.Outages.columns:
                parameters['OutageFactor']['val'][i, :] = data.Outages[u].values
            else:
                logging.warning('Outages factors not found for unit ' + u + '. Assuming no outages')

    # Participation to the reserve market
    values = np.array([s in config['ReserveParticipation'] for s in sets['t']], dtype='bool')
    parameters['Reserve'] = {'sets': sets_param['Reserve'], 'val': values}

    # Technologies
    for unit in range(Nunits):
        idx = sets['t'].index(Plants_merged['Technology'][unit])
        parameters['Technology']['val'][unit, idx] = True

    # Fuels
    for unit in range(Nunits):
        idx = sets['f'].index(Plants_merged['Fuel'][unit])
        parameters['Fuel']['val'][unit, idx] = True

    # Location
    for i in range(len(sets['n'])):
        parameters['Location']['val'][:, i] = (Plants_merged['Zone'] == config['countries'][i]).values

    # CHPType parameter:
    sets['chp_type'] = ['Extraction','Back-Pressure', 'P2H']
    parameters['CHPType'] = define_parameter(['chp','chp_type'],sets,value=0)
    for i,u in enumerate(sets['chp']):
        if u in Plants_chp.index:
            if Plants_chp.loc[u,'CHPType'].lower() == 'extraction':
                parameters['CHPType']['val'][i,0] = 1
            elif Plants_chp.loc[u,'CHPType'].lower() == 'back-pressure':
                parameters['CHPType']['val'][i,1] = 1
            elif Plants_chp.loc[u,'CHPType'].lower() == 'p2h':
                parameters['CHPType']['val'][i,2] = 1
            else:
                logging.error('CHPType not valid for plant ' + u)
                sys.exit(1)

    # Initial Power
    if 'InitialPower' in Plants_merged:
        parameters['PowerInitial']['val'] = Plants_merged['InitialPower'].values
    else:
        for i in range(Nunits):
            # Nuclear and Fossil Gas greater than 350 MW are up (assumption):
            if Plants_merged['Fuel'][i] in ['GAS', 'NUC'] and Plants_merged['PowerCapacity'][i] > 350:
                parameters['PowerInitial']['val'][i] = (Plants_merged['PartLoadMin'][i] + 1) / 2 * \
                                                    Plants_merged['PowerCapacity'][i]
            # Config variables:
    sets['x_config'] = ['FirstDay', 'LastDay', 'RollingHorizon Length', 'RollingHorizon LookAhead','ValueOfLostLoad','QuickStartShare','CostOfSpillage','WaterValue']
    sets['y_config'] = ['year', 'month', 'day', 'val']
    dd_begin = idx_long[4]
    dd_end = idx_long[-2]

#TODO: integrated the parameters (VOLL, Water value, etc) from the excel config file
    values = np.array([
        [dd_begin.year, dd_begin.month, dd_begin.day, 0],
        [dd_end.year, dd_end.month, dd_end.day, 0],
        [0, 0, config['HorizonLength'], 0],
        [0, 0, config['LookAhead'], 0],
        [0, 0, 0, 1e5],     # Value of lost load
        [0, 0, 0, 0.5],       # allowed Share of quick start units in reserve
        [0, 0, 0, 1],       # Cost of spillage (EUR/MWh)
        [0, 0, 0, 100],       # Value of water (for unsatisfied water reservoir levels, EUR/MWh)
    ])
    parameters['Config'] = {'sets': ['x_config', 'y_config'], 'val': values}
    return parameters

def build_sim_dir(config, sim, sets, parameters, SimData):

    LP = config['SimulationType'] == 'LP' or config['SimulationType'] == 'LP clustered'
    CEP = config['CEP'] == 1
    gdx_out = "Inputs.gdx"


    if config['WriteGDX']:
        write_variables(config['GAMS_folder'], gdx_out, [sets, parameters])

    # if the sim variable was not defined:
    if 'sim' not in locals():
        logging.error('Please provide a path where to store the DispaSET inputs (in the "sim" variable)')
        sys.exit(1)

    if not os.path.exists(sim):
        os.makedirs(sim)

    def replace_text_by_dict(text, dic):
        '''Replace dictionary items in text'''
        for i, j in dic.items():
            text = text.replace(i, j)
        return text

        
    gams_file_changes = {'LP':LP, 'CEP':CEP}
    changes_infile_string = {'LP': ('$setglobal LPFormulation 0','$setglobal LPFormulation 1'), 'CEP': ('$setglobal CEPFormulation 0', '$setglobal CEPFormulation 1')}
    gams_file_changes_list = {changes_infile_string[k][0]: changes_infile_string[k][1] for k,v in gams_file_changes.items() if v == True}  #filter based on selection
    if len(gams_file_changes_list)>0:
        fin = open(os.path.join(GMS_FOLDER, 'UCM_h.gms'))
        fout = open(os.path.join(sim,'UCM_h.gms'), "wt")
        for line in fin:
            fout.write(replace_text_by_dict(line, gams_file_changes_list))
        fin.close()
        fout.close()
    else:
        shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h.gms'),
                        os.path.join(sim, 'UCM_h.gms'))

    gmsfile = open(os.path.join(sim, 'UCM.gpr'), 'w')
    gmsfile.write(
        '[PROJECT] \n \n[RP:UCM_H] \n1= \n[OPENWINDOW_1] \nFILE0=UCM_h.gms \nFILE1=UCM_h.gms \nMAXIM=1 \nTOP=50 \nLEFT=50 \nHEIGHT=400 \nWIDTH=400')
    gmsfile.close()
    shutil.copyfile(os.path.join(GMS_FOLDER, 'writeresults.gms'),
                    os.path.join(sim, 'writeresults.gms'))
    # Create cplex option file
    cplex_options = {'epgap': 0.05, # TODO: For the moment hardcoded, it has to be moved to a config file
                    'numericalemphasis': 0,
                    'scaind': 1,
                    'lpmethod': 0,
                    'relaxfixedinfeas': 0,
                    'mipstart':1,
                    'epint':0}

    lines_to_write = ['{} {}'.format(k, v) for k, v in cplex_options.items()]
    with open(os.path.join(sim, 'cplex.opt'), 'w') as f:
        for line in lines_to_write:
            f.write(line + '\n')

    logging.debug('Using gams file from ' + GMS_FOLDER)
    if config['WriteGDX']:
        shutil.copy(gdx_out, sim + '/')
        os.remove(gdx_out)
    # Copy bat file to generate gdx file directly from excel:
    shutil.copy(os.path.join(GMS_FOLDER, 'makeGDX.bat'),
                os.path.join(sim, 'makeGDX.bat'))

    if config['WriteExcel']:
        write_to_excel(sim, [sets, parameters])

    if config['WritePickle']:
        try:
            import cPickle as pickle
        except ImportError:
            import pickle
        with open(os.path.join(sim, 'Inputs.p'), 'wb') as pfile:
            pickle.dump(SimData, pfile, protocol=pickle.HIGHEST_PROTOCOL)
    logging.info('Build finished')

    if os.path.isfile(commons['logfile']):
        shutil.copy(commons['logfile'], os.path.join(sim, 'warn_preprocessing.log'))


def load_sets(Nhours_long, look_ahead, plants_index, plants_sto_index, Plants_chp_index, countries, Interconnections, plants_uc=None):
    '''Build the sets/indices'''
    plants_ue = [plant for plant in plants_index if plant not in plants_uc]
    sets = {
        'h': [str(x + 1) for x in range(Nhours_long)],
        'z': [str(x + 1) for x in range(Nhours_long - look_ahead * 24)],
        'mk': ['DA', '2U', '2D'],
        'n': countries,
        'u': plants_index,
        'l': Interconnections,
        'f': commons['Fuels'],
        'p': ['CO2'],
        's': plants_sto_index,
        'chp': Plants_chp_index,
        't': commons['Technologies'],
        'tr': commons['tech_renewables'],
        'uc': plants_uc,
        'ue': plants_ue
    }

    return sets

def load_params():
    '''Load all parameters'''
    sets_param = {
        'AvailabilityFactor': ['u', 'h'],
        'CHPPowerToHeat': ['chp'],
        'CHPPowerLossFactor': ['chp'],
        'CHPMaxHeat': ['chp'],
        'CostFixed': ['u'],
        'CostHeatSlack': ['chp', 'h'],
        'CostLoadShedding': ['n', 'h'],
        'CostRampUp': ['u'],
        'CostRampDown': ['u'],
        'CostShutDown': ['u'],
        'CostStartUp': ['u'],
        'CostVariable': ['u', 'h'],
        'Curtailment': ['n'],
        'Demand': ['mk', 'n', 'h'],
        'Efficiency': ['u'],
        'EmissionMaximum': ['n', 'p'],
        'EmissionRate': ['u', 'p'],
        'FlowMaximum': ['l', 'h'],
        'FlowMinimum': ['l', 'h'],
        'Fuel': ['u', 'f'],
        'HeatDemand': ['chp', 'h'],
        'Investment': ['uc'],
        'EconomicLifetime': ['uc'],
        'LineNode': ['l', 'n'],
        'LoadShedding': ['n', 'h'],
        'Location': ['u', 'n'],
        'Markup': ['u', 'h'],
        'Nunits': ['u'],
        'OutageFactor': ['u', 'h'],
        'PartLoadMin': ['u'],
        'PowerCapacity': ['u'],
        'PowerInitial': ['u'],
        'PriceTransmission': ['l', 'h'],
        'RampUpMaximum': ['u'],
        'RampDownMaximum': ['u'],
        'RampStartUpMaximum': ['u'],
        'RampShutDownMaximum': ['u'],
        'Reserve': ['t'],
        'StorageCapacity': ['u'],
        'StorageChargingCapacity': ['s'],
        'StorageChargingEfficiency': ['s'],
        'StorageDischargeEfficiency': ['s'],
        'StorageSelfDischarge': ['u'],
        'StorageInflow': ['s', 'h'],
        'StorageInitial': ['s'],
        'StorageMinimum': ['s'],
        'StorageOutflow': ['s', 'h'],
        'StorageProfile': ['s', 'h'],
        'Technology': ['u', 't'],
        'TimeUpMinimum': ['u'],
        'TimeDownMinimum': ['u'],
    }

    return sets_param


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

    if isinstance(inputs,str) or isinstance(inputs,unicode):
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
            write_variables(SimData['config']['GAMS_folder'], 'Inputs.gdx', [SimData['sets'], SimData['parameters']])
            shutil.copy('Inputs.gdx', dest_path + '/')
            os.remove('Inputs.gdx')
    return SimData


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

    if isinstance(inputs,str) or isinstance(inputs,unicode):
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
            shutil.copytree(path,dest_path)
            logging.info('Created simulation environment directory ' + dest_path)
        logging.info('Writing input files to ' + dest_path)
        import cPickle
        with open(os.path.join(dest_path, 'Inputs.p'), 'wb') as pfile:
            cPickle.dump(SimData, pfile, protocol=cPickle.HIGHEST_PROTOCOL)
        if write_gdx:
            write_variables(SimData['config']['GAMS_folder'], 'Inputs.gdx', [SimData['sets'], SimData['parameters']])
            shutil.copy('Inputs.gdx', dest_path + '/')
            os.remove('Inputs.gdx')
    return SimData
