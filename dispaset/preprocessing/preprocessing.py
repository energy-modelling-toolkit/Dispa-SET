# -*- coding: utf-8 -*-
"""
This is the main file of the DispaSET pre-processing tool. It comprises a single function that generated the DispaSET simulation environment.

@author: S. Quoilin
"""

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

from .data_check import check_units, check_chp, check_sto, check_p2h, check_heat_demand, check_df, isStorage, check_MinMaxFlows,check_AvailabilityFactors, check_clustering, check_temperatures
from .utils import clustering, interconnections, incidence_matrix, select_units
from .data_handler import UnitBasedTable,NodeBasedTable,merge_series, define_parameter, load_time_series

from ..solve import solve_GAMS_simple

from ..misc.gdx_handler import write_variables, gdx_to_list, gdx_to_dataframe
from ..common import commons  # Load fuel types, technologies, timestep, etc:

GMS_FOLDER = os.path.join(os.path.dirname(__file__), '..', 'GAMS')

def get_git_revision_tag():
    """Get version of DispaSET used for this run. tag + commit hash"""
    from subprocess import check_output
    try:
        return check_output(["git", "describe", "--tags", "--always"]).strip()
    except:
        return 'NA'

def build_simulation(config, profiles=None):
    """
    This function reads the DispaSET config, loads the specified data,
    processes it when needed, and formats it in the proper DispaSET format.
    The output of the function is a directory with all inputs and simulation files required to run a DispaSET simulation

    :param config:        Dictionary with all the configuration fields loaded from the excel file. Output of the 'LoadConfig' function.
    :param plot_load:     Boolean used to display a plot of the demand curves in the different zones
    :param profiles:      Profiles from mid term scheduling simulations
    """
    dispa_version = str(get_git_revision_tag())
    logging.info('New build started. DispaSET version: ' + dispa_version)
    # %%################################################################################################################
    #####################################   Main Inputs    ############################################################
    ###################################################################################################################

    # Boolean variable to check wether it is milp or lp:
    LP = config['SimulationType'] == 'LP' or config['SimulationType'] == 'LP clustered'

    # Day/hour corresponding to the first and last days of the simulation:
    # Note that the first available data corresponds to 2015.01.31 (23.00) and the
    # last day with data is 2015.12.31 (22.00)
    __, m_start, d_start, __, __, __ = config['StartDate']
    y_end, m_end, d_end, _, _, _ = config['StopDate']
    config['StopDate'] = (y_end, m_end, d_end, 23, 59, 00)  # updating stopdate to the end of the day
    
    # Indexes of the simulation:
    config['idx'] = pd.DatetimeIndex(pd.date_range(start=pd.datetime(*config['StartDate']),
                                             end=pd.datetime(*config['StopDate']),
                                             freq=commons['TimeStep'])).tz_localize(None)


    # Indexes for the whole year considered in StartDate
    idx_year = pd.DatetimeIndex(pd.date_range(start=pd.datetime(*(config['StartDate'][0],1,1,0,0)),
                                                        end=pd.datetime(*(config['StartDate'][0],12,31,23,59,59)),
                                                        freq=commons['TimeStep'])
                                          )
    
    # Indexes including the last look-ahead period
    enddate_long = config['idx'][-1] + dt.timedelta(days=config['LookAhead'])
    idx_long = pd.DatetimeIndex(pd.date_range(start=config['idx'][0], end=enddate_long, freq=commons['TimeStep']))
    Nhours_long = len(idx_long)
    config['idx_long'] = idx_long

    # %%#################################################################################################################
    #####################################   Data Loading    ###########################################################
    ###################################################################################################################

    # Start and end of the simulation:
    delta = config['idx'][-1] - config['idx'][0]
    days_simulation = delta.days + 1

    # Defining missing configuration fields for backwards compatibility:
    if not isinstance(config['default']['CostLoadShedding'],(float,int)):
        config['default']['CostLoadShedding'] = 1000
    if not isinstance(config['default']['CostHeatSlack'],(float,int)):
        config['default']['CostHeatSlack'] = 50
    
    # Load :
    Load = NodeBasedTable('Demand',config)
    PeakLoad = Load.max()

    if config['modifiers']['Demand'] != 1:
        logging.info('Scaling load curve by a factor ' + str(config['modifiers']['Demand']))
        Load = Load * config['modifiers']['Demand']
        PeakLoad = PeakLoad * config['modifiers']['Demand']

    # Interconnections:
    if os.path.isfile(config['Interconnections']):
        flows = load_time_series(config,config['Interconnections']).fillna(0)
    else:
        logging.warning('No historical flows will be considered (no valid file provided)')
        flows = pd.DataFrame(index=config['idx_long'])
    if os.path.isfile(config['NTC']):
        NTC = load_time_series(config,config['NTC']).fillna(0)
    else:
        logging.warning('No NTC values will be considered (no valid file provided)')
        NTC = pd.DataFrame(index=config['idx_long'])

    # Load Shedding:
    LoadShedding = NodeBasedTable('LoadShedding',config,default=config['default']['LoadShedding'])
    CostLoadShedding = NodeBasedTable('CostLoadShedding',config,default=config['default']['CostLoadShedding'])

    # Power plants:
    plants = pd.DataFrame()
    if os.path.isfile(config['PowerPlantData']):
        plants = pd.read_csv(config['PowerPlantData'])
    elif '##' in config['PowerPlantData']:
        for z in config['zones']:
            path = config['PowerPlantData'].replace('##', str(z))
            tmp = pd.read_csv(path)
            plants = plants.append(tmp, ignore_index=True)
    # reomve invalide power plants:
    plants = select_units(plants,config)

    # Some columns can be in two format (absolute or per unit). If not specified, they are set to zero:
    for key in ['StartUpCost','NoLoadCost']:
        if key in plants:
            pass
        elif key+'_pu' in plants:
            plants[key] = plants[key+'_pu'] * plants['PowerCapacity']
        else:
            plants[key] = 0
    # check plant list:
    check_units(config, plants)
    # If not present, add the non-compulsory fields to the units table:
    for key in ['CHPPowerLossFactor','CHPPowerToHeat','CHPType','STOCapacity','STOSelfDischarge','STOMaxChargingPower','STOChargingEfficiency', 'CHPMaxHeat']:
        if key not in plants.columns:
            plants[key] = np.nan


    # Defining the hydro storages:
    plants_sto = plants[[u in commons['tech_storage'] for u in plants['Technology']]]
    # check storage plants:
    check_sto(config, plants_sto)
    
    # Defining the CHPs:
    plants_chp = plants[[str(x).lower() in commons['types_CHP'] for x in plants['CHPType']]]

    # Defining the P2H units:
    plants_p2h = plants[plants['Technology']=='P2HT']
    
    # All heating units:
    plants_heat = plants_chp.append(plants_p2h)

    Outages = UnitBasedTable(plants,'Outages',config,fallbacks=['Unit','Technology'])
    AF = UnitBasedTable(plants,'RenewablesAF',config,fallbacks=['Unit','Technology'],default=1,RestrictWarning=commons['tech_renewables'])

    # Reservoir levels
    """
    :profiles:      if turned on reservoir levels are overwritten by newly calculated ones 
                    from the mid term scheduling simulations
    """
    if profiles is None:
        ReservoirLevels = UnitBasedTable(plants_sto,'ReservoirLevels',config,fallbacks=['Unit','Technology','Zone'],default=0)
    else:
        ReservoirLevels = UnitBasedTable(plants_sto,'ReservoirLevels',config,fallbacks=['Unit','Technology','Zone'],default=0)
        MidTermSchedulingProfiles = profiles
        ReservoirLevels.update(MidTermSchedulingProfiles)

    ReservoirScaledInflows = UnitBasedTable(plants_sto,'ReservoirScaledInflows',config,fallbacks=['Unit','Technology','Zone'],default=0)
    HeatDemand = UnitBasedTable(plants_heat,'HeatDemand',config,fallbacks=['Unit'],default=0)
    CostHeatSlack = UnitBasedTable(plants_heat,'CostHeatSlack',config,fallbacks=['Unit','Zone'],default=config['default']['CostHeatSlack'])
    Temperatures = NodeBasedTable('Temperatures',config)

    # data checks:
    check_AvailabilityFactors(plants,AF)
    check_heat_demand(plants,HeatDemand)
    check_temperatures(plants,Temperatures)

    # Fuel prices:
    fuels = ['PriceOfNuclear', 'PriceOfBlackCoal', 'PriceOfGas', 'PriceOfFuelOil', 'PriceOfBiomass', 'PriceOfCO2', 'PriceOfLignite', 'PriceOfPeat']
    FuelPrices = {}
    for fuel in fuels:
        FuelPrices[fuel] = NodeBasedTable(fuel,config,default=config['default'][fuel])

    # Interconnections:
    [Interconnections_sim, Interconnections_RoW, Interconnections] = interconnections(config['zones'], NTC, flows)

    if len(Interconnections_sim.columns) > 0:
        NTCs = Interconnections_sim.reindex(config['idx_long'])
    else:
        NTCs = pd.DataFrame(index=config['idx_long'])
    Inter_RoW = Interconnections_RoW.reindex(config['idx_long'])

    # Clustering of the plants:
    Plants_merged, mapping = clustering(plants, method=config['SimulationType'])
    # Check clustering:
    check_clustering(plants,Plants_merged)

    # Renaming the columns to ease the production of parameters:
    Plants_merged.rename(columns={'StartUpCost': 'CostStartUp',
                                  'RampUpMax': 'RampUpMaximum',
                                  'RampDownMax': 'RampDownMaximum',
                                  'MinUpTime': 'TimeUpMinimum',
                                  'MinDownTime': 'TimeDownMinimum',
                                  'RampingCost': 'CostRampUp',
                                  'STOCapacity': 'StorageCapacity',
                                  'STOMaxChargingPower': 'StorageChargingCapacity',
                                  'STOChargingEfficiency': 'StorageChargingEfficiency',
                                  'STOSelfDischarge': 'StorageSelfDischarge',
                                  'CO2Intensity': 'EmissionRate'}, inplace=True)
    
    for key in ['TimeUpMinimum','TimeDownMinimum']:
        if any([not x.is_integer() for x in Plants_merged[key].fillna(0).values.astype('float')]):
            logging.warning(key + ' in the power plant data has been rounded to the nearest integer value')
            Plants_merged.loc[:,key] = Plants_merged[key].fillna(0).values.astype('int32')

    if not len(Plants_merged.index.unique()) == len(Plants_merged):
        # Very unlikely case:
        logging.error('plant indexes not unique!')
        sys.exit(1)

    # Apply scaling factors:
    if config['modifiers']['Solar'] != 1:
        logging.info('Scaling Solar Capacity by a factor ' + str(config['modifiers']['Solar']))
        for u in Plants_merged.index:
            if Plants_merged.Technology[u] == 'PHOT':
                Plants_merged.loc[u, 'PowerCapacity'] = Plants_merged.loc[u, 'PowerCapacity'] * config['modifiers']['Solar']
    if config['modifiers']['Wind'] != 1:
        logging.info('Scaling Wind Capacity by a factor ' + str(config['modifiers']['Wind']))
        for u in Plants_merged.index:
            if Plants_merged.Technology[u] == 'WTON' or Plants_merged.Technology[u] == 'WTOF':
                Plants_merged.loc[u, 'PowerCapacity'] = Plants_merged.loc[u, 'PowerCapacity'] * config['modifiers']['Wind']
    if config['modifiers']['Storage'] != 1:
        logging.info('Scaling Storage Power and Capacity by a factor ' + str(config['modifiers']['Storage']))
        for u in Plants_merged.index:
            if isStorage(Plants_merged.Technology[u]):
                Plants_merged.loc[u, 'PowerCapacity'] = Plants_merged.loc[u, 'PowerCapacity'] * config['modifiers']['Storage']
                Plants_merged.loc[u, 'StorageCapacity'] = Plants_merged.loc[u, 'StorageCapacity'] * config['modifiers']['Storage']
                Plants_merged.loc[u, 'StorageChargingCapacity'] = Plants_merged.loc[u, 'StorageChargingCapacity'] * config['modifiers']['Storage']

    # Defining the hydro storages:
    Plants_sto = Plants_merged[[u in commons['tech_storage'] for u in Plants_merged['Technology']]]
    # check storage plants:
    check_sto(config, Plants_sto,raw_data=False)
    # Defining the CHPs:
    Plants_chp = Plants_merged[[x.lower() in commons['types_CHP'] for x in Plants_merged['CHPType']]].copy()
    # check chp plants:
    check_chp(config, Plants_chp)
    # For all the chp plants correct the PowerCapacity, which is defined in cogeneration mode in the inputs and in power generation model in the optimization model
    for u in Plants_chp.index:
        PowerCapacity = Plants_chp.loc[u, 'PowerCapacity']

        if Plants_chp.loc[u,'CHPType'].lower() == 'p2h':
            PurePowerCapacity = PowerCapacity
        else:
            if pd.isnull(Plants_chp.loc[u,'CHPMaxHeat']):  # If maximum heat is not defined, then it is defined as the intersection between two lines
                MaxHeat = PowerCapacity / Plants_chp.loc[u,'CHPPowerToHeat']
                Plants_chp.loc[u, 'CHPMaxHeat'] = 'inf'
            else:
                MaxHeat = Plants_chp.loc[u, 'CHPMaxHeat']
            PurePowerCapacity = PowerCapacity + Plants_chp.loc[u,'CHPPowerLossFactor'] * MaxHeat
        Plants_merged.loc[u,'PartLoadMin'] = Plants_merged.loc[u,'PartLoadMin'] * PowerCapacity / PurePowerCapacity  # FIXME: Is this correct?
        Plants_merged.loc[u,'PowerCapacity'] = PurePowerCapacity
    
    #for p2h, assigning the cop to the efficiency in order not to multiply the variables:
    Plants_merged.loc[Plants_merged['Technology']=='P2HT','Efficiency']=Plants_merged[Plants_merged['Technology']=='P2HT']['COP']
    Plants_p2h = Plants_merged[Plants_merged['Technology']=='P2HT'].copy()
    # check chp plants:
    check_p2h(config, Plants_p2h)
    # All heating plants:
    Plants_heat = Plants_chp.append(Plants_p2h)
        
    
    # Get the hydro time series corresponding to the original plant list: #FIXME Unused variable ?
    #StorageFormerIndexes = [s for s in plants.index if
    #                        plants['Technology'][s] in commons['tech_storage']]

    # Same with the CHPs:
    # Get the heat demand time series corresponding to the original plant list:
#    CHPFormerIndexes = [s for s in plants.index if
#                            plants['CHPType'][s] in commons['types_CHP']]
#    for s in CHPFormerIndexes:  # for all the old plant indexes
#        # get the old plant name corresponding to s:
#        oldname = plants['Unit'][s]
#        # newname = mapping['NewIndex'][s] #FIXME Unused variable ?
#        if oldname not in HeatDemand:
#            logging.warning('No heat demand profile found for CHP plant "' + str(oldname) + '". Assuming zero')
#            HeatDemand[oldname] = 0
#        if oldname not in CostHeatSlack:
#            logging.warning('No heat cost profile found for CHP plant "' + str(oldname) + '". Assuming zero')
#            CostHeatSlack[oldname] = 0
# 

    # merge the outages:
    for i in plants.index:  # for all the old plant indexes
        # get the old plant name corresponding to s:
        oldname = plants['Unit'][i]
        newname = mapping['NewIndex'][i]

    # Merging the time series relative to the clustered power plants:
    ReservoirScaledInflows_merged = merge_series(plants, ReservoirScaledInflows, mapping, method='WeightedAverage', tablename='ScaledInflows')
    ReservoirLevels_merged = merge_series(plants, ReservoirLevels, mapping, tablename='ReservoirLevels')
    Outages_merged = merge_series(plants, Outages, mapping, tablename='Outages')
    HeatDemand_merged = merge_series(plants, HeatDemand, mapping, tablename='HeatDemand',method='Sum')
    AF_merged = merge_series(plants, AF, mapping, tablename='AvailabilityFactors')
    CostHeatSlack_merged = merge_series(plants, CostHeatSlack, mapping, tablename='CostHeatSlack')
    

    # %%
    # checking data
    check_df(Load, StartDate=config['idx'][0], StopDate=config['idx'][-1], name='Load')
    check_df(AF_merged, StartDate=config['idx'][0], StopDate=config['idx'][-1],
             name='AF_merged')
    check_df(Outages_merged, StartDate=config['idx'][0], StopDate=config['idx'][-1], name='Outages_merged')
    check_df(Inter_RoW, StartDate=config['idx'][0], StopDate=config['idx'][-1], name='Inter_RoW')
    for key in FuelPrices:
        check_df(FuelPrices[key], StartDate=config['idx'][0], StopDate=config['idx'][-1], name=key)
    check_df(NTCs, StartDate=config['idx'][0], StopDate=config['idx'][-1], name='NTCs')
    check_df(ReservoirLevels_merged, StartDate=config['idx'][0], StopDate=config['idx'][-1],
             name='ReservoirLevels_merged')
    check_df(ReservoirScaledInflows_merged, StartDate=config['idx'][0], StopDate=config['idx'][-1],
             name='ReservoirScaledInflows_merged')
    check_df(HeatDemand_merged, StartDate=config['idx'][0], StopDate=config['idx'][-1],
             name='HeatDemand_merged')
    check_df(CostHeatSlack_merged, StartDate=config['idx'][0], StopDate=config['idx'][-1],
             name='CostHeatSlack_merged')
    check_df(LoadShedding, StartDate=config['idx'][0], StopDate=config['idx'][-1],
             name='LoadShedding')
    check_df(CostLoadShedding, StartDate=config['idx'][0], StopDate=config['idx'][-1],
             name='CostLoadShedding')

#    for key in Renewables:
#        check_df(Renewables[key], StartDate=config['idx'][0], StopDate=config['idx'][-1],
#                 name='Renewables["' + key + '"]')


    # %%################################################################################################################
    ############################################   Sets    ############################################################
    ###################################################################################################################

    # The sets are defined within a dictionary:
    sets = {}
    sets['h'] = [str(x + 1) for x in range(Nhours_long)]
    sets['z'] = [str(x + 1) for x in range(Nhours_long - config['LookAhead'] * 24)]
    sets['mk'] = ['DA', '2U', '2D']
    sets['n'] = config['zones']
    sets['au'] = Plants_merged.index.tolist()
    sets['l'] = Interconnections
    sets['f'] = commons['Fuels']
    sets['p'] = ['CO2']
    sets['s'] = Plants_sto.index.tolist()
    sets['u'] = Plants_merged[Plants_merged['Technology']!='P2H'].index.tolist()
    sets['chp'] = Plants_chp.index.tolist()
    sets['p2h'] = Plants_p2h.index.tolist()
    sets['th'] = Plants_heat.index.tolist()
    sets['t'] = commons['Technologies']
    sets['tr'] = commons['tech_renewables']

    ###################################################################################################################
    ############################################   Parameters    ######################################################
    ###################################################################################################################

    Nunits = len(Plants_merged)
    parameters = {}

    # Each parameter is associated with certain sets, as defined in the following list:
    sets_param = {}
    sets_param['AvailabilityFactor'] = ['au', 'h']
    sets_param['CHPPowerToHeat'] = ['chp']
    sets_param['CHPPowerLossFactor'] = ['chp']
    sets_param['CHPMaxHeat'] = ['chp']
    sets_param['CostFixed'] = ['au']
    sets_param['CostHeatSlack'] = ['th','h']
    sets_param['CostLoadShedding'] = ['n','h']
    sets_param['CostRampUp'] = ['au']
    sets_param['CostRampDown'] = ['au']
    sets_param['CostShutDown'] = ['au']
    sets_param['CostStartUp'] = ['au']
    sets_param['CostVariable'] = ['au', 'h']
    sets_param['Curtailment'] = ['n']
    sets_param['Demand'] = ['mk', 'n', 'h']
    sets_param['Efficiency'] = ['au']
    sets_param['EmissionMaximum'] = ['n', 'p']
    sets_param['EmissionRate'] = ['au', 'p']
    sets_param['FlowMaximum'] = ['l', 'h']
    sets_param['FlowMinimum'] = ['l', 'h']
    sets_param['Fuel'] = ['au', 'f']
    sets_param['HeatDemand'] = ['th','h']
    sets_param['LineNode'] = ['l', 'n']
    sets_param['LoadShedding'] = ['n','h']
    sets_param['Location'] = ['au', 'n']
    sets_param['Markup'] = ['au', 'h']
    sets_param['Nunits'] = ['au']
    sets_param['OutageFactor'] = ['au', 'h']
    sets_param['PartLoadMin'] = ['au']
    sets_param['PowerCapacity'] = ['au']
    sets_param['PowerInitial'] = ['au']
    sets_param['PriceTransmission'] = ['l', 'h']
    sets_param['RampUpMaximum'] = ['au']
    sets_param['RampDownMaximum'] = ['au']
    sets_param['RampStartUpMaximum'] = ['au']
    sets_param['RampShutDownMaximum'] = ['au']
    sets_param['Reserve'] = ['t']
    sets_param['StorageCapacity'] = ['au']
    sets_param['StorageChargingCapacity'] = ['s']
    sets_param['StorageChargingEfficiency'] = ['s']
    sets_param['StorageDischargeEfficiency'] = ['s']
    sets_param['StorageSelfDischarge'] = ['au']
    sets_param['StorageInflow'] = ['s', 'h']
    sets_param['StorageInitial'] = ['s']
    sets_param['StorageMinimum'] = ['s']
    sets_param['StorageOutflow'] = ['s', 'h']
    sets_param['StorageProfile'] = ['s', 'h']
    sets_param['Technology'] = ['au', 't']
    sets_param['TimeUpMinimum'] = ['au']
    sets_param['TimeDownMinimum'] = ['au']

    # Define all the parameters and set a default value of zero:
    for var in sets_param:
        parameters[var] = define_parameter(sets_param[var], sets, value=0)

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

    # %%
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
        if profiles is not None:
            if (config['InitialFinalReservoirLevel'] == 0) or (config['InitialFinalReservoirLevel'] == ""):
                if s in ReservoirLevels_merged:
                    # get the time
                    parameters['StorageInitial']['val'][i] = ReservoirLevels_merged[s][idx_long[0]] * \
                                                             Plants_sto['StorageCapacity'][s] * Plants_sto['Nunits'][s]
                    parameters['StorageProfile']['val'][i, :] = ReservoirLevels_merged[s][idx_long].values
                    if any(ReservoirLevels_merged[s] > 1):
                        logging.warning(s + ': The reservoir level is sometimes higher than its capacity!')
                else:
                    logging.warning( 'Could not find reservoir level data for storage plant ' + s + '. Assuming 50% of capacity')
                    parameters['StorageInitial']['val'][i] = 0.5 * Plants_sto['StorageCapacity'][s]
                    parameters['StorageProfile']['val'][i, :] = 0.5
            else:
                if (config['default']['ReservoirLevelInitial'] > 1) or (config['default']['ReservoirLevelFinal'] > 1):
                    logging.warning(s + ': The initial or final reservoir levels are higher than its capacity!' )
                parameters['StorageInitial']['val'][i] = config['default']['ReservoirLevelInitial'] * Plants_sto['StorageCapacity'][s]
                parameters['StorageProfile']['val'][i, :] = config['default']['ReservoirLevelFinal']
        else:
            if s in ReservoirLevels_merged:
                # get the time
                parameters['StorageInitial']['val'][i] = ReservoirLevels_merged[s][idx_long[0]] * \
                                                         Plants_sto['StorageCapacity'][s] * Plants_sto['Nunits'][s]
                parameters['StorageProfile']['val'][i, :] = ReservoirLevels_merged[s][idx_long].values
                if any(ReservoirLevels_merged[s] > 1):
                    logging.warning(s + ': The reservoir level is sometimes higher than its capacity!')
            else:
                logging.warning( 'Could not find reservoir level data for storage plant ' + s + '. Assuming 50% of capacity')
                parameters['StorageInitial']['val'][i] = 0.5 * Plants_sto['StorageCapacity'][s]
                parameters['StorageProfile']['val'][i, :] = 0.5

    # Storage Inflows:
    for i, s in enumerate(sets['s']):
        if s in ReservoirScaledInflows_merged:
            parameters['StorageInflow']['val'][i, :] = ReservoirScaledInflows_merged[s][idx_long].values * \
                                                       Plants_sto['PowerCapacity'][s]
    # CHP time series:
    for i, u in enumerate(sets['th']):
        if u in HeatDemand_merged:
            parameters['HeatDemand']['val'][i, :] = HeatDemand_merged[u][idx_long].values
            parameters['CostHeatSlack']['val'][i, :] = CostHeatSlack_merged[u][idx_long].values

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
    if len(AF_merged.columns) != 0:
        for i, u in enumerate(sets['au']):
            if u in AF_merged.columns:
                parameters['AvailabilityFactor']['val'][i, :] = AF_merged[u].values


    # Demand
    # Dayahead['NL'][1800:1896] = Dayahead['NL'][1632:1728]
    reserve_2U_tot = {i: (np.sqrt(10 * PeakLoad[i] + 150 ** 2) - 150) for i in Load.columns}
    reserve_2D_tot = {i: (0.5 * reserve_2U_tot[i]) for i in Load.columns}

    values = np.ndarray([len(sets['mk']), len(sets['n']), len(sets['h'])])
    for i in range(len(sets['n'])):
        values[0, i, :] = Load[sets['n'][i]]
        values[1, i, :] = reserve_2U_tot[sets['n'][i]]
        values[2, i, :] = reserve_2D_tot[sets['n'][i]]
    parameters['Demand'] = {'sets': sets_param['Demand'], 'val': values}

    # Emission Rate:
    parameters['EmissionRate']['val'][:, 0] = Plants_merged['EmissionRate'].values

    # Load Shedding:
    for i, c in enumerate(sets['n']):
        parameters['LoadShedding']['val'][i] = LoadShedding[c] * PeakLoad[c]
        parameters['CostLoadShedding']['val'][i] = CostLoadShedding[c]

    # %%#################################################################################################################################################################################################
    # Variable Cost
    # Equivalence dictionary between fuel types and price entries in the config sheet:
    FuelEntries = {'BIO':'PriceOfBiomass', 'GAS':'PriceOfGas', 'HRD':'PriceOfBlackCoal', 'LIG':'PriceOfLignite', 'NUC':'PriceOfNuclear', 'OIL':'PriceOfFuelOil', 'PEA':'PriceOfPeat'}
    for unit in range(Nunits):
        c = Plants_merged['Zone'][unit]    # zone to which the unit belongs
        found = False
        for FuelEntry in FuelEntries:
            if Plants_merged['Fuel'][unit] == FuelEntry:
                parameters['CostVariable']['val'][unit, :] = FuelPrices[FuelEntries[FuelEntry]][c] / Plants_merged['Efficiency'][unit] + \
                                                             Plants_merged['EmissionRate'][unit] * FuelPrices['PriceOfCO2'][c]
                found = True
        # Special case for biomass plants, which are not included in EU ETS:
        if Plants_merged['Fuel'][unit] == 'BIO':
            parameters['CostVariable']['val'][unit, :] = FuelPrices['PriceOfBiomass'][c] / Plants_merged['Efficiency'][
                unit]
            found = True
        if not found:
            logging.warning('No fuel price value has been found for fuel ' + Plants_merged['Fuel'][unit] + ' in unit ' + \
                  Plants_merged['Unit'][unit] + '. A null variable cost has been assigned')

    # %%#################################################################################################################################################################################################

    # Maximum Line Capacity
    for i, l in enumerate(sets['l']):
        if l in NTCs.columns:
            parameters['FlowMaximum']['val'][i, :] = NTCs[l]
        if l in Inter_RoW.columns:
            parameters['FlowMaximum']['val'][i, :] = Inter_RoW[l]
            parameters['FlowMinimum']['val'][i, :] = Inter_RoW[l]
    # Check values:
    check_MinMaxFlows(parameters['FlowMinimum']['val'],parameters['FlowMaximum']['val'])
    
    parameters['LineNode'] = incidence_matrix(sets, 'l', parameters, 'LineNode')

    # Outage Factors
    if len(Outages_merged.columns) != 0:
        for i, u in enumerate(sets['au']):
            if u in Outages_merged.columns:
                parameters['OutageFactor']['val'][i, :] = Outages_merged[u].values
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
        parameters['Location']['val'][:, i] = (Plants_merged['Zone'] == config['zones'][i]).values

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

    # Check if values are specified in the config and overwrite the default ones
    values_lst = ['ValueOfLostLoad','ShareOfQuickStartUnits','PriceOfSpillage','WaterValue']
    values_val = [1e5,0.5,1,100]
    values = np.array([[dd_begin.year, dd_begin.month, dd_begin.day, 0],
        [dd_end.year, dd_end.month, dd_end.day, 0],
        [0, 0, config['HorizonLength'], 0],
        [0, 0, config['LookAhead'], 0]])
    for vl in values_lst:
        if vl in config['default']:
            if config['default'][vl] == "":
                if vl == 'ValueOfLostLoad':
                    tmp_value = np.array([[0,0,0,values_val[0]]])
                if vl == 'ShareOfQuickStartUnits':
                    tmp_value = np.array([[0,0,0,values_val[1]]])
                if vl == 'PriceOfSpillage':
                    tmp_value = np.array([[0,0,0,values_val[2]]])
                if vl == 'WaterValue':
                    tmp_value = np.array([[0,0,0,values_val[3]]])
            else:
                tmp_value = np.array([[0,0,0,config['default'][vl]]])
        else:
            if vl == 'ValueOfLostLoad':
                tmp_value = np.array([[0,0,0,values_val[0]]])
            if vl == 'ShareOfQuickStartUnits':
                tmp_value = np.array([[0,0,0,values_val[1]]])
            if vl == 'PriceOfSpillage':
                tmp_value = np.array([[0,0,0,values_val[2]]])
            if vl == 'WaterValue':
                tmp_value = np.array([[0,0,0,values_val[3]]])
        values = np.append(values, tmp_value, axis=0)
    parameters['Config'] = {'sets': ['x_config', 'y_config'], 'val': values}

    # %%#################################################################################################################
    ######################################   Simulation Environment     ################################################
    ####################################################################################################################

    # Output folder:
    sim = config['SimulationDirectory']

    # Simulation data:
    SimData = {'sets': sets, 'parameters': parameters, 'config': config, 'units': Plants_merged, 'version': dispa_version}

    # list_vars = []
    gdx_out = "Inputs.gdx"
    if config['WriteGDX']:
        write_variables(config['GAMS_folder'], gdx_out, [sets, parameters])

    # if the sim variable was not defined:
    if 'sim' not in locals():
        logging.error('Please provide a path where to store the DispaSET inputs (in the "sim" variable)')
        sys.exit(1)

    if not os.path.exists(sim):
        os.makedirs(sim)
    if LP:
        fin = open(os.path.join(GMS_FOLDER, 'UCM_h.gms'))
        fout = open(os.path.join(sim,'UCM_h.gms'), "wt")
        for line in fin:
            fout.write(line.replace('$setglobal LPFormulation 0', '$setglobal LPFormulation 1'))
        fin.close()
        fout.close()
        # additionally allso copy UCM_h_simple.gms
        shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'),
                        os.path.join(sim, 'UCM_h_simple.gms'))
    else:
        shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h.gms'),
                        os.path.join(sim, 'UCM_h.gms'))
        # additionally allso copy UCM_h_simple.gms
        shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'),
                        os.path.join(sim, 'UCM_h_simple.gms'))
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

def get_temp_sim_results(config, gams_dir=None, cache=False, temp_path='.pickle'):
    """
    This function reads the simulation environment folder once it has been solved and loads
    the input variables together with the results.

    :param path:                Relative path to the simulation environment folder (current path by default)
    :param cache:               If true, caches the simulation results in a pickle file for faster loading the next time
    :param temp_path:           Temporary path to store the cache file
    :returns inputs,results:    Two dictionaries with all the outputs 
    """

    resultfile = config['SimulationDirectory'] + '/Results_simple.gdx'
    results = gdx_to_dataframe(gdx_to_list(gams_dir, resultfile, varname='all', verbose=True), fixindex=True,
                               verbose=True)
    results['OutputStorageLevel'] = results['OutputStorageLevel'].reindex(list(range(results['OutputStorageLevel'].index.min(),
                                                         results['OutputStorageLevel'].index.max() + 1)), fill_value=0)
    return results

def mid_term_scheduling(config, zones, profiles=None):
    """
    This function reads the DispaSET config file, searches for active zones,
    loads data for each zone individually and solves model using UCM_h_simple.gms
    
    :config:                    Read config file
    """

    # Day/hour corresponding to the first and last days of the simulation:
    # Note that the first available data corresponds to 2015.01.31 (23.00) and the
    # last day with data is 2015.12.31 (22.00)
    __, m_start, d_start, __, __, __ = config['StartDate']
    y_end, m_end, d_end, _, _, _ = config['StopDate']
    config['StopDate'] = (y_end, m_end, d_end, 23, 59, 00)  # updating stopdate to the end of the day
    
    # Indexes of the simualtion:
    idx_std = pd.DatetimeIndex(start=pd.datetime(*config['StartDate']), 
                               end=pd.datetime(*config['StopDate']),
                               freq=commons['TimeStep'])
    idx_utc_noloc = idx_std - dt.timedelta(hours=1)
    idx = idx_utc_noloc
    
    # Checking which type of hydro scheduling simulation is specified in the config file:
    # Solving reservoir levels for each zone individually
    if config['HydroScheduling'] == 'Zonal':
        no_of_zones = len(zones)
        results = {}
        temp_results = {}
        i = 0
        for c in zones:
            i = i + 1  
            logging.info('(Curently simulating Zone): ' + str(i) + ' out of ' + str(no_of_zones))
            temp_config = dict(config)
            temp_config['zones'] = [c]                    # Override zone that needs to be simmulated
            _ = build_simulation(temp_config)       # Create temporary SimData   
            r = solve_GAMS_simple(temp_config['SimulationDirectory'], 
                                  temp_config['GAMS_folder'])
            temp_results[c] = get_temp_sim_results(config)
    #        print('Zones simulated: ' + str(i) + '/' + str(no_of_zones))
        temp = pd.DataFrame()
        for c in zones:
            if 'OutputStorageLevel' not in temp_results[c]:
                logging.critical('Storage levels in zone ' + c + ' were not computed, please check that storage units '
                                'are present in the ' + c + ' power plant database! If not, unselect ' + c + ' form the '
                                'zones in the MTS module')
                sys.exit(0)
            else:
                results[c] = dict(temp_results[c]['OutputStorageLevel'])
                r = pd.DataFrame.from_dict(results[c], orient='columns')
                results_t = pd.concat([temp, r], axis = 1)
                temp = results_t
        temp = temp.set_index(idx)
        temp = temp.rename(columns={col: col.split(' - ')[1] for col in temp.columns})
    # Solving reservoir levels for all regions simultaneously
    elif config['HydroScheduling'] == 'Regional':
        if zones == None:
            temp_config = dict(config)
        else:
            temp_config = dict(config)
            temp_config['zones'] = zones               # Override zones that need to be simmulated
        _ = build_simulation(temp_config)        # Create temporary SimData   
        r = solve_GAMS_simple(temp_config['SimulationDirectory'], temp_config['GAMS_folder'])
        temp_results = get_temp_sim_results(config)
        if 'OutputStorageLevel' not in temp_results:
            logging.critical('Storage levels in the selected region were not computed, please check that storage units '
                             'are present in the power plant database! If not, unselect zones with no storage units form '
                             'the zones in the MTS module')
            sys.exit(0)
        else:
            results = dict(temp_results['OutputStorageLevel'])
        temp = pd.DataFrame()
        r = pd.DataFrame.from_dict(results, orient='columns')
        results_t = pd.concat([temp, r], axis = 1)
        temp = results_t
        temp = temp.set_index(idx)
        temp = temp.rename(columns={col: col.split(' - ')[1] for col in temp.columns})        
    else:
        logging.info('Mid term scheduling turned off')
    return temp

def build_full_simulation(config, mts_plot=None):
    '''
    Dispa-SET function that builds different simulation environments based on the hydro scheduling option in the config file
    Hydro scheduling options:
        Off      -  Hydro scheduling turned off, normal call of BuildSimulation function
        Zonal    -  Zonal variation of hydro scheduling, if zones are not individually specified in a list (e.a. zones = ['AT','DE'])
                    hydro scheduling is imposed on all active zones from the Config file
        Regional -  Regional variation of hydro scheduling, if zones from a specific region are not individually specified in a list 
                    (e.a. zones = ['AT','DE']), hydro scheduling is imposed on all active zones from the Config file simultaneously
    
    :config:                    Read config file
    :zones_mts:                 List of zones where new reservoir levels should be calculated eg. ['AT','BE',...'UK']
    :mts_plot:                  If ms_plot = True indicative plot with temporary computed reservoir levels is displayed
    '''
    y_start, m_start, d_start, __, __, __ = config['StartDate']
    y_stop, m_stop, d_stop, __, __, __ = config['StopDate']
    # Check existance of hydro scheduling module in the config file
    # Hydro scheduling turned off, build_simulation performed without temporary computed reservoir levels
    if (config['HydroScheduling'] == 'Off') or (config['HydroScheduling'] == ""):
        logging.info('Simulation without mid therm scheduling')
        SimData = build_simulation(config)
    # Hydro scheduling per Zone    
    else:
        # Dates for the mid term scheduling
        if config['HydroSchedulingHorizon'] == 'Annual':
            config['StartDate'] = (y_start, 1, 1, 00, 00, 00)   # updating start date to the beginning of the year
            config['StopDate'] = (y_start, 12, 31, 23, 59, 00)  # updating stopdate to the end of the year
            logging.info('Hydro scheduling is performed for the period between 01.01.' + str(y_start) + ' and 12.31.' + str(y_start)) 
        else:
            logging.info('Hydro scheduling is performed between Start and Stop dates!') 
        # Mid term scheduling zone selection and new profile calculation
        if config['mts_zones'] == None:
            new_profiles = mid_term_scheduling(config, config['zones'])
            logging.info('Simulation with all zones selected')
        else:
            new_profiles = mid_term_scheduling(config, config['mts_zones'])
        # Plot new profiles
        if mts_plot == True:
            new_profiles.plot()
            logging.info('Simulation with specified zones selected')
        else:
            logging.info('No temporary profiles selected for display')
         # Build simulation data with new profiles
        config['StartDate'] = (y_start, m_start, d_start, 00, 00, 00)   # updating start date to the beginning of the year
        config['StopDate'] = (y_stop, m_stop, d_stop, 23, 59, 00)  # updating stopdate to the end of the year
        SimData = build_simulation(config, new_profiles)
    return SimData
