import datetime as dt
import logging
import os
import shutil
import sys

import numpy as np
import pandas as pd
from future.builtins import int

from .data_check import check_units, check_sto, check_AvailabilityFactors, check_heat_demand, \
    check_temperatures, check_clustering, isStorage, check_chp, check_p2h, check_h2, check_df, check_MinMaxFlows, \
    check_FlexibleDemand, check_reserves, check_PtLDemand, check_heat
from .data_handler import NodeBasedTable, load_time_series, UnitBasedTable, merge_series, define_parameter, \
    load_geo_data, GenericTable
from .utils import select_units, interconnections, clustering, EfficiencyTimeSeries, incidence_matrix, pd_timestep
from .reserves import percentage_reserve, probabilistic_reserve, generic_reserve

from .. import __version__
from ..common import commons
from ..misc.gdx_handler import write_variables

GMS_FOLDER = os.path.join(os.path.dirname(__file__), '..', 'GAMS')


def build_single_run(config, profiles=None, PtLDemand=None, MTS=0):
    """
    This function reads the DispaSET config, loads the specified data,
    processes it when needed, and formats it in the proper DispaSET format.
    The output of the function is a directory with all inputs and simulation files required to run a DispaSET simulation

    :param config:        Dictionary with all the configuration fields loaded from the excel file. Output of the
                          'LoadConfig' function.
    :param profiles:      Profiles from mid term scheduling simulations
    :param PtLDemand:     Profiles for PtL demand from mid ter scheduling simulations
    :param MTS:           Boolean, 1 if in MTS
    :PtLDemand:           PtLDemand from mid term scheduling simulations
    """
    dispa_version = __version__
    logging.info('New build started. DispaSET version: ' + dispa_version)
    # %%###############################################################################################################
    #####################################   Main Inputs    ############################################################
    ###################################################################################################################

    # Boolean variable to check wether it is milp or lp:
    LP = config['SimulationType'] == 'LP' or config['SimulationType'] == 'LP clustered'

    # check time steps:
    if config['DataTimeStep'] != 1:
        logging.critical('The data time step can only be 1 hour in this version of Dispa-SET. A value of ' + str(
            config['DataTimeStep']) + ' hours was provided')
        sys.exit(1)
    if config['SimulationTimeStep'] not in (1, 24):
        logging.critical(
            'The simulation time step can only be 1 or 24 hour in this version of Dispa-SET. A value of ' + str(
                config['DataTimeStep']) + ' hours was provided')
        sys.exit(1)
    # Day/hour corresponding to the first and last days of the simulation:
    __, m_start, d_start, __, __, __ = config['StartDate']
    y_end, m_end, d_end, _, _, _ = config['StopDate']
    config['StopDate'] = (y_end, m_end, d_end, 23, 59, 00)  # updating stopdate to the end of the day

    # Indexes of the simulation:
    config['idx'] = pd.date_range(start=dt.datetime(*config['StartDate']),
                                  end=dt.datetime(*config['StopDate']),
                                  freq=commons['TimeStep']).tz_localize(None)

    # Indexes for the whole year considered in StartDate
    idx_year = pd.date_range(start=dt.datetime(*(config['StartDate'][0], 1, 1, 0, 0)),
                             end=dt.datetime(*(config['StartDate'][0], 12, 31, 23, 59, 59)),
                             freq=commons['TimeStep'])

    # Indexes including the last look-ahead period
    enddate_long = config['idx'][-1] + dt.timedelta(days=config['LookAhead'])
    idx_long = pd.date_range(start=config['idx'][0], end=enddate_long, freq=commons['TimeStep'])
    Nhours_long = len(idx_long)
    config['idx_long'] = idx_long

    # Indexes of with the specified simulation time step:
    idx_sim = pd.date_range(start=config['idx'][0], end=enddate_long, freq=pd_timestep(config['SimulationTimeStep']))
    config['idx_sim'] = idx_sim
    Nsim = len(idx_sim)

    # %%###############################################################################################################
    #####################################   Data Loading    ###########################################################
    ###################################################################################################################

    # Start and end of the simulation:
    delta = config['idx'][-1] - config['idx'][0]
    days_simulation = delta.days + 1

    # Load :
    Load = NodeBasedTable('Demand', config)
    PeakLoad = Load.max()

    if config['modifiers']['Demand'] != 1:
        logging.info('Scaling load curve by a factor ' + str(config['modifiers']['Demand']))
        Load = Load * config['modifiers']['Demand']
        PeakLoad = PeakLoad * config['modifiers']['Demand']

    # Interconnections:
    if os.path.isfile(config['Interconnections']):
        flows = load_time_series(config, config['Interconnections']).fillna(0)
    else:
        logging.warning('No historical flows will be considered (no valid file provided)')
        flows = pd.DataFrame(index=config['idx_long'])
    if os.path.isfile(config['NTC']):
        NTC = load_time_series(config, config['NTC']).fillna(0)
    else:
        logging.warning('No NTC values will be considered (no valid file provided)')
        NTC = pd.DataFrame(index=config['idx_long'])
    if os.path.isfile(config['PriceTransmission']):
        PriceTransmission_raw = load_time_series(config, config['PriceTransmission'])
    else:
        PriceTransmission_raw = pd.DataFrame(index=config['idx_long'])

    # Geo data
    if 'GeoData' in config and os.path.isfile(config['GeoData']):
        geo = load_geo_data(config['GeoData'], header=0)
    else:
        logging.warning('No geo spatial data available')
        geo = None

    # Load Shedding:
    LoadShedding = NodeBasedTable('LoadShedding', config, default=config['default']['LoadShedding'])
    CostLoadShedding = NodeBasedTable('CostLoadShedding', config, default=config['default']['CostLoadShedding'])
    ShareOfFlexibleDemand = NodeBasedTable('ShareOfFlexibleDemand', config,
                                           default=config['default']['ShareOfFlexibleDemand'])
    Reserve2D = NodeBasedTable('Reserve2D', config, default=None)
    Reserve2U = NodeBasedTable('Reserve2U', config, default=None)

    # Curtailment:
    CostCurtailment = NodeBasedTable('CostCurtailment', config, default=config['default']['CostCurtailment'])

    # Power plants:
    plants = pd.DataFrame()
    if os.path.isfile(config['PowerPlantData']):
        plants = pd.read_csv(config['PowerPlantData'],
                             na_values=commons['na_values'],
                             keep_default_na=False)
    elif '##' in config['PowerPlantData']:
        for z in config['zones']:
            path = config['PowerPlantData'].replace('##', str(z))
            tmp = pd.read_csv(path, na_values=commons['na_values'],
                              keep_default_na=False)
            plants = plants.append(tmp, ignore_index=True, sort=False)
    # remove invalid power plants:
    plants = select_units(plants, config)

    # Some columns can be in two format (absolute or per unit). If not specified, they are set to zero:
    for key in ['StartUpCost', 'NoLoadCost']:
        if key in plants:
            pass
        elif key + '_pu' in plants:
            plants[key] = plants[key + '_pu'] * plants['PowerCapacity']
        else:
            plants[key] = 0

    # If not present, add the non-compulsory fields to the units table:
    for key in ['CHPPowerLossFactor', 'CHPPowerToHeat', 'CHPType', 'STOCapacity', 'STOSelfDischarge',
                'STOMaxChargingPower', 'STOChargingEfficiency', 'CHPMaxHeat', 'WaterWithdrawal',
                'WaterConsumption']:
        if key not in plants.columns:
            plants[key] = np.nan

    # If the thermal and h2 zones are not defined in the units table, define one individual zone per power plant:
    for key in ['Zone_th', 'Zone_h2']:
        if key in plants.columns:
            plants[key] = plants[key].fillna('')
        else:
            plants[key] = plants['Unit']
            logging.info('No "' + key + '" header was found in the units table. One individual zone is defined per '
                                        'power plant')

    # check plant list:
    check_units(config, plants)

    # Defining the hydro storages:
    plants_sto = plants[[u in commons['tech_storage'] for u in plants['Technology']]]
    check_sto(config, plants_sto)

    # Defining the thermal storages:
    plants_thms = plants[[u in commons['tech_thermal_storage'] for u in plants['Technology']]]
    check_sto(config, plants_thms)

    # Merging all MTS storage units
    plants_all_sto = plants_sto.append(plants_thms)

    # Defining the heat only units:
    plants_heat = plants[[u in commons['tech_heat'] for u in plants['Technology']]]
    check_heat(config, plants_heat)

    # Defining the CHPs:
    plants_chp = plants[[str(x).lower() in commons['types_CHP'] for x in plants['CHPType']]]
    check_chp(config, plants_chp)

    # Defining the P2H units:
    plants_p2h = plants[[u in commons['tech_p2ht'] for u in plants['Technology']]]
    check_p2h(config, plants_p2h)

    # All heating units:
    plants_heat = plants_heat.append(plants_chp)
    plants_heat = plants_heat.append(plants_p2h)

    # Defining the P2H units:
    plants_h2 = plants[plants['Technology'] == 'P2GS']
    check_h2(config, plants_h2)

    Outages = UnitBasedTable(plants, 'Outages', config, fallbacks=['Unit', 'Technology'])
    AF = UnitBasedTable(plants, 'RenewablesAF', config, fallbacks=['Unit', 'Technology'], default=1,
                        RestrictWarning=commons['tech_renewables'])
    AF = AF.apply(pd.to_numeric)

    ReservoirLevels = UnitBasedTable(plants_all_sto, 'ReservoirLevels', config,
                                     fallbacks=['Unit', 'Technology', 'Zone'],
                                     default=0)
    ReservoirScaledInflows = UnitBasedTable(plants_sto, 'ReservoirScaledInflows', config,
                                            fallbacks=['Unit', 'Technology', 'Zone'], default=0)
    # HeatDemand = UnitBasedTable(plants_heat, 'HeatDemand', config, fallbacks=['Unit'], default=0)
    # CostHeatSlack = UnitBasedTable(plants_heat, 'CostHeatSlack', config, fallbacks=['Unit', 'Zone'],
    #                                default=config['default']['CostHeatSlack'])
    Temperatures = NodeBasedTable('Temperatures', config)

    if plants_h2.empty is True:
        H2RigidDemand = pd.DataFrame(index=config['idx_long'])
        H2FlexibleDemand = pd.DataFrame(index=config['idx_long'])
        CostH2Slack = pd.DataFrame(index=config['idx_long'])
    else:
        H2RigidDemand = UnitBasedTable(plants_h2, 'H2RigidDemand', config, fallbacks=['Unit'], default=0)
        H2FlexibleDemand = UnitBasedTable(plants_h2, 'H2FlexibleDemand', config, fallbacks=['Unit'], default=0)
        CostH2Slack = UnitBasedTable(plants_h2, 'CostH2Slack', config, fallbacks=['Unit', 'Zone'],
                                     default=config['default']['CostH2Slack'])

    # Detecting thermal zones:
    zones_th = plants_heat['Zone_th'].unique().tolist()
    if '' in zones_th:
        zones_th.remove('')

    HeatDemand = GenericTable(zones_th, 'HeatDemand', config, default=0)
    CostHeatSlack = GenericTable(zones_th, 'CostHeatSlack', config, default=config['default']['CostHeatSlack'])

    # Update reservoir levels with newly computed ones from the mid-term scheduling
    if profiles is not None:
        plants_all_sto.set_index(plants_all_sto.loc[:, 'Unit'], inplace=True, drop=True)
        for key in profiles.columns:
            if key not in ReservoirLevels.columns:
                logging.warning('The reservoir profile "' + key + '" provided by the MTS is not found in the '
                                                                  'ReservoirLevels table')
            elif key in list(ReservoirLevels.loc[:, plants_all_sto['Technology'] == 'SCSP'].columns):
                ReservoirLevels[key] = config['default']['ReservoirLevelInitial']
                logging.info('The reservoir profile "' + key + '" can not be seleceted for MTS, instead, default value '
                                                               'of: ' + str(
                    config['default']['ReservoirLevelInitial']) + ' will be used')
            else:
                ReservoirLevels[key].update(profiles[key])
                logging.info(
                    'The reservoir profile "' + key + '" provided by the MTS is used as target reservoir level')
    # Update PtL demand (H2FlexibleDemand with demand from mid term scheduling)
    if PtLDemand is not None and any(H2FlexibleDemand) > 0:
        for key in PtLDemand.columns:
            if key not in H2FlexibleDemand.columns:
                logging.warning('The H2 flexible demand "' + key + '" provided by the MTS is not found in the '
                                                                   'H2FlexibleDemand table')
            else:
                H2FlexibleDemand[key].update(PtLDemand[key])

    # data checks:
    check_AvailabilityFactors(plants, AF)
    check_heat_demand(plants, HeatDemand, zones_th)
    check_temperatures(plants, Temperatures)
    check_FlexibleDemand(ShareOfFlexibleDemand)
    check_reserves(Reserve2D, Reserve2U, Load)

    # Fuel prices:
    fuels = ['PriceOfNuclear', 'PriceOfBlackCoal', 'PriceOfGas', 'PriceOfFuelOil', 'PriceOfBiomass', 'PriceOfCO2',
             'PriceOfLignite', 'PriceOfPeat']
    FuelPrices = {}
    for fuel in fuels:
        FuelPrices[fuel] = NodeBasedTable(fuel, config, default=config['default'][fuel])

    # Interconnections:
    [Interconnections_sim, Interconnections_RoW, Interconnections] = interconnections(config['zones'], NTC, flows)

    if len(Interconnections_sim.columns) > 0:
        NTCs = Interconnections_sim.reindex(config['idx_long'])
    else:
        NTCs = pd.DataFrame(index=config['idx_long'])
    Inter_RoW = Interconnections_RoW.reindex(config['idx_long'])
    PriceTransmission = pd.DataFrame(index=NTCs.index, columns=NTCs.columns)
    for l in (NTCs.columns.tolist() + Inter_RoW.columns.tolist()):
        if l in PriceTransmission_raw:
            PriceTransmission[l] = PriceTransmission_raw[l].reindex(PriceTransmission.index)
        else:
            PriceTransmission[l] = config['default']['PriceTransmission']
            if config['default']['PriceTransmission'] > 0:
                logging.warning('No detailed values were found the transmission prices of line ' + l +
                                '. Using default value ' + str(config['default']['PriceTransmission']))

    # Clustering of the plants:
    Plants_merged, mapping = clustering(plants, method=config['SimulationType'])
    # Check clustering:
    check_clustering(plants, Plants_merged)

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

    for key in ['TimeUpMinimum', 'TimeDownMinimum']:
        if any([not x.is_integer() for x in Plants_merged[key].fillna(0).values.astype('float')]):
            logging.warning(key + ' in the power plant data has been rounded to the nearest integer value')
            Plants_merged.loc[:, key] = Plants_merged[key].fillna(0).values.astype('int32')

    if not len(Plants_merged.index.unique()) == len(Plants_merged):
        # Very unlikely case:
        logging.error('plant indexes not unique!')
        sys.exit(1)

    # Apply scaling factors:
    if config['modifiers']['Solar'] != 1:
        logging.info('Scaling Solar Capacity by a factor ' + str(config['modifiers']['Solar']))
        for u in Plants_merged.index:
            if Plants_merged.Technology[u] == 'PHOT':
                Plants_merged.loc[u, 'PowerCapacity'] = Plants_merged.loc[u, 'PowerCapacity'] * config['modifiers'][
                    'Solar']
            if Plants_merged.Technology[u] == 'SOTH':
                Plants_merged.loc[u, 'PowerCapacity'] = Plants_merged.loc[u, 'PowerCapacity'] * config['modifiers'][
                    'Solar']
    if config['modifiers']['Wind'] != 1:
        logging.info('Scaling Wind Capacity by a factor ' + str(config['modifiers']['Wind']))
        for u in Plants_merged.index:
            if Plants_merged.Technology[u] == 'WTON' or Plants_merged.Technology[u] == 'WTOF':
                Plants_merged.loc[u, 'PowerCapacity'] = Plants_merged.loc[u, 'PowerCapacity'] * config['modifiers'][
                    'Wind']
    if config['modifiers']['Storage'] != 1:
        logging.info('Scaling Storage Power and Capacity by a factor ' + str(config['modifiers']['Storage']))
        for u in Plants_merged.index:
            if isStorage(Plants_merged.Technology[u]):
                Plants_merged.loc[u, 'PowerCapacity'] = Plants_merged.loc[u, 'PowerCapacity'] * config['modifiers'][
                    'Storage']
                Plants_merged.loc[u, 'StorageCapacity'] = Plants_merged.loc[u, 'StorageCapacity'] * config['modifiers'][
                    'Storage']
                Plants_merged.loc[u, 'StorageChargingCapacity'] = Plants_merged.loc[u, 'StorageChargingCapacity'] * \
                                                                  config['modifiers']['Storage']

    # Defining the hydro storages:
    Plants_sto = Plants_merged[[u in commons['tech_storage'] for u in Plants_merged['Technology']]]
    # Defining the thermal storages:
    Plants_thms = Plants_merged[[u in commons['tech_thermal_storage'] for u in Plants_merged['Technology']]]
    # check storage plants:
    check_sto(config, Plants_sto, raw_data=False)
    check_sto(config, Plants_thms, raw_data=False)
    # Defining the CHPs:
    Plants_chp = Plants_merged[[x.lower() in commons['types_CHP'] for x in Plants_merged['CHPType']]].copy()
    # check chp plants:
    check_chp(config, Plants_chp)
    # For all the chp plants correct the PowerCapacity, which is defined in cogeneration mode in the inputs and
    # in power generation model in the optimization model
    for u in Plants_chp.index:
        PowerCapacity = Plants_chp.loc[u, 'PowerCapacity']

        if Plants_chp.loc[u, 'CHPType'].lower() == 'p2h':
            PurePowerCapacity = PowerCapacity
        else:
            # If maximum heat is not defined, then it is defined as the intersection between two lines
            if pd.isnull(Plants_chp.loc[u, 'CHPMaxHeat']):
                MaxHeat = PowerCapacity / Plants_chp.loc[u, 'CHPPowerToHeat']
                Plants_chp.loc[u, 'CHPMaxHeat'] = 'inf'
            else:
                MaxHeat = Plants_chp.loc[u, 'CHPMaxHeat']
            PurePowerCapacity = PowerCapacity + Plants_chp.loc[u, 'CHPPowerLossFactor'] * MaxHeat
        # FIXME: Is this correct?
        Plants_merged.loc[u, 'PartLoadMin'] = Plants_merged.loc[u, 'PartLoadMin'] * PowerCapacity / PurePowerCapacity
        Plants_merged.loc[u, 'PowerCapacity'] = PurePowerCapacity

    # Filter power to heat units
    Plants_p2h = Plants_merged[[u in commons['tech_p2ht'] for u in Plants_merged['Technology']]].copy()
    # check power to heat plants:
    check_p2h(config, Plants_p2h)

    # Filter heat only plants
    Plants_heat_only = Plants_merged[[u in commons['tech_heat'] for u in Plants_merged['Technology']]].copy()
    # Check heat only units
    check_heat(config, Plants_heat_only)

    # All heating plants:
    Plants_heat = Plants_heat_only.copy()
    Plants_heat = Plants_heat.append(Plants_chp)
    Plants_heat = Plants_heat.append(Plants_p2h)
    Plants_heat = Plants_heat.append(Plants_thms)

    Plants_h2 = Plants_merged[Plants_merged['Technology'] == 'P2GS'].copy()
    # check chp plants:
    check_h2(config, Plants_h2)

    # Water storage
    Plants_wat = Plants_merged[(Plants_merged['Fuel'] == 'WAT') & (Plants_merged['Technology'] != 'HROR')].copy()

    # Calculating the efficiency time series for each unit:
    Efficiencies = EfficiencyTimeSeries(config, Plants_merged, Temperatures)

    # Reserve calculation
    reserve_2U_tot = pd.DataFrame(index=Load.index, columns=Load.columns)
    reserve_2D_tot = pd.DataFrame(index=Load.index, columns=Load.columns)
    for z in Load.columns:
        if config['ReserveCalculation'] == 'Exogenous':
            if z in Reserve2U and z in Reserve2D:
                reserve_2U_tot[z] = Reserve2U[z]
                reserve_2D_tot[z] = Reserve2D[z]
            else:
                logging.critical('Exogenous reserve requirements (2D and 2U) not found for zone ' + z)
                sys.exit(1)
        else:
            if z in Reserve2U and z in Reserve2D:
                logging.info('Using exogenous reserve data for zone ' + z)
                reserve_2U_tot[z] = Reserve2U[z]
                reserve_2D_tot[z] = Reserve2D[z]
            elif config['ReserveCalculation'] == 'Percentage':
                logging.info('Using percentage-based reserve sizing for zone ' + z)
                reserve_2U_tot[z], reserve_2D_tot[z] = percentage_reserve(config, plants, Load, AF, z)
            elif config['ReserveCalculation'] == 'Probabilistic':
                logging.info('Using probabilistic reserve sizing for zone ' + z)
                reserve_2U_tot[z], reserve_2D_tot[z] = probabilistic_reserve(config, plants, Load, AF, z)
            else:
                logging.info('Using generic reserve calculation for zone ' + z)
                reserve_2U_tot[z], reserve_2D_tot[z] = generic_reserve(Load[z])

                # %% Store all times series and format

    # Formatting all time series (merging, resempling) and store in the FinalTS dict
    finalTS = {'Load': Load, 'Reserve2D': reserve_2D_tot, 'Reserve2U': reserve_2U_tot,
               'Efficiencies': Efficiencies, 'NTCs': NTCs, 'Inter_RoW': Inter_RoW,
               'LoadShedding': LoadShedding, 'CostLoadShedding': CostLoadShedding, 'CostCurtailment': CostCurtailment,
               'ScaledInflows': ReservoirScaledInflows, 'ReservoirLevels': ReservoirLevels,
               'Outages': Outages, 'AvailabilityFactors': AF, 'CostHeatSlack': CostHeatSlack,
               'HeatDemand': HeatDemand, 'ShareOfFlexibleDemand': ShareOfFlexibleDemand,
               'PriceTransmission': PriceTransmission, 'CostH2Slack': CostH2Slack,
               'H2RigidDemand': H2RigidDemand, 'H2FlexibleDemand': H2FlexibleDemand}

    # Merge the following time series with weighted averages
    for key in ['ScaledInflows', 'Outages', 'AvailabilityFactors', 'CostH2Slack']:
        finalTS[key] = merge_series(Plants_merged, plants, finalTS[key], tablename=key)
    # Merge the following time series by summing
    for key in ['H2RigidDemand', 'H2FlexibleDemand']:
        finalTS[key] = merge_series(Plants_merged, plants, finalTS[key], tablename=key, method='Sum')
    # Merge the following time series by weighted average based on storage capacity
    for key in ['ReservoirLevels']:
        finalTS[key] = merge_series(Plants_merged, plants, finalTS[key], tablename=key, method='StorageWeightedAverage')

    # Check that all times series data is available with the specified data time step:
    for key in FuelPrices:
        check_df(FuelPrices[key], StartDate=config['idx'][0], StopDate=config['idx'][-1], name=key)
    for key in finalTS:
        check_df(finalTS[key], StartDate=config['idx'][0], StopDate=config['idx'][-1], name=key)

    # Resemple to the required time step
    if config['DataTimeStep'] != config['SimulationTimeStep']:
        for key in FuelPrices:
            if len(FuelPrices[key].columns) > 0:
                FuelPrices[key] = FuelPrices[key].resample(pd_timestep(config['SimulationTimeStep'])).mean()
        for key in finalTS:
            if len(finalTS[key].columns) > 0:
                finalTS[key] = finalTS[key].resample(pd_timestep(config['SimulationTimeStep'])).mean()

    # %%###############################################################################################################
    ############################################   Sets    ############################################################
    ###################################################################################################################

    # The sets are defined within a dictionary:
    sets = {}
    sets['h'] = [str(x + 1) for x in range(Nsim)]
    sets['z'] = [str(x + 1) for x in range(int(Nsim - config['LookAhead'] * 24 / config['SimulationTimeStep']))]
    sets['mk'] = ['DA', '2U', '2D', 'Flex']
    sets['n'] = config['zones']
    sets['n_th'] = zones_th
    sets['au'] = Plants_merged.index.tolist()
    sets['l'] = Interconnections
    sets['f'] = commons['Fuels']
    sets['p'] = ['CO2']
    sets['s'] = Plants_sto.index.tolist()
    sets['u'] = Plants_merged[[u in [x for x in commons['Technologies'] if x not in commons['tech_heat'] +
                                     commons['tech_p2ht'] + commons['tech_thermal_storage']]
                               for u in Plants_merged['Technology']]].index.tolist()
    sets['chp'] = Plants_chp.index.tolist()
    sets['p2h'] = Plants_p2h.index.tolist()
    sets['p2h2'] = Plants_h2.index.tolist()
    sets['th'] = Plants_heat.index.tolist()
    sets['thms'] = Plants_thms.index.tolist()
    sets['t'] = commons['Technologies']
    sets['tr'] = commons['tech_renewables']
    sets['wat'] = Plants_wat.index.tolist()
    sets['hu'] = Plants_heat_only.index.tolist()
    sets['asu'] = Plants_merged[[u in [x for x in commons['Technologies'] if x in commons['tech_storage'] +
                                       commons['tech_thermal_storage']] for u in
                                 Plants_merged['Technology']]].index.tolist()

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
    sets_param['CostHeatSlack'] = ['n_th', 'h']
    sets_param['CostH2Slack'] = ['p2h2', 'h']
    sets_param['CostLoadShedding'] = ['n', 'h']
    sets_param['CostRampUp'] = ['au']
    sets_param['CostRampDown'] = ['au']
    sets_param['CostShutDown'] = ['au']
    sets_param['CostStartUp'] = ['au']
    sets_param['CostVariable'] = ['au', 'h']
    sets_param['Curtailment'] = ['n']
    sets_param['CostCurtailment'] = ['n', 'h']
    sets_param['Demand'] = ['mk', 'n', 'h']
    sets_param['Efficiency'] = ['p2h', 'h']
    sets_param['EmissionMaximum'] = ['n', 'p']
    sets_param['EmissionRate'] = ['au', 'p']
    sets_param['FlowMaximum'] = ['l', 'h']
    sets_param['FlowMinimum'] = ['l', 'h']
    sets_param['Fuel'] = ['au', 'f']
    sets_param['HeatDemand'] = ['n_th', 'h']
    sets_param['LineNode'] = ['l', 'n']
    sets_param['LoadShedding'] = ['n', 'h']
    sets_param['Location'] = ['au', 'n']
    sets_param['Location_th'] = ['au', 'n_th']
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
    sets_param['Reserve'] = [
        'au']  # changed this also in the gams file(in the definition and in the equations satifying the reserve demand)
    sets_param['StorageCapacity'] = ['au']
    sets_param['StorageChargingCapacity'] = ['au']
    sets_param['StorageChargingEfficiency'] = ['au']
    sets_param['StorageDischargeEfficiency'] = ['asu']
    sets_param['StorageSelfDischarge'] = ['au']
    sets_param['StorageInflow'] = ['s', 'h']
    sets_param['StorageInitial'] = ['asu']
    sets_param['StorageMinimum'] = ['s']
    sets_param['StorageOutflow'] = ['s', 'h']
    sets_param['StorageProfile'] = ['asu', 'h']
    sets_param['Technology'] = ['au', 't']
    sets_param['TimeUpMinimum'] = ['au']
    sets_param['TimeDownMinimum'] = ['au']
    sets_param['PtLDemandInput'] = ['p2h2', 'h']
    sets_param['MaxCapacityPtL'] = ['p2h2']
    sets_param['H2Demand'] = ['s', 'h']

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
    for var in ['Technology', 'Fuel', 'Reserve', 'Location', 'Location_th']:
        parameters[var] = define_parameter(sets_param[var], sets, value='bool')

    # %%
    # List of parameters whose value is known, and provided in the dataframe Plants_merged.
    for var in ['PowerCapacity', 'PartLoadMin', 'TimeUpMinimum', 'TimeDownMinimum', 'CostStartUp',
                'CostRampUp', 'StorageCapacity', 'StorageSelfDischarge', 'StorageChargingCapacity']:
        parameters[var]['val'] = Plants_merged[var].values

    # List of parameters whose value is not necessarily specified in the dataframe Plants_merged
    for var in ['Nunits']:
        if var in Plants_merged:
            parameters[var]['val'] = Plants_merged[var].values

    # List of parameters whose value is known, and provided in the dataframe Plants_sto.
    for var in ['StorageChargingEfficiency']:
        # parameters[var]['val'] = Plants_sto[var].values
        parameters[var]['val'] = Plants_merged[var].values

    # The storage discharge efficiency is actually given by the unit efficiency:
    parameters['StorageDischargeEfficiency']['val'] = np.concatenate((Plants_sto['Efficiency'].values,
                                                                      Plants_thms['Efficiency'].values), axis=None)

    # List of parameters whose value is known, and provided in the dataframe Plants_chp
    for var in ['CHPPowerToHeat', 'CHPPowerLossFactor', 'CHPMaxHeat']:
        parameters[var]['val'] = Plants_chp[var].values

    # Particular treatment of MaxCapacityPtL that is not a time-series and
    # that is given separetly from the Power plant database 
    if 'H2FlexibleCapacity' in config and config['H2FlexibleCapacity'] != '':
        MaxCapacityPtL = pd.read_csv(config['H2FlexibleCapacity'], index_col=0, keep_default_na=False)
        for i, u in enumerate(sets['p2h2']):
            for unit in MaxCapacityPtL.index:
                if unit in u:
                    parameters['MaxCapacityPtL']['val'][i] = MaxCapacityPtL.loc[unit]

                    # Storage profile and initial state:
    for i, s in enumerate(sets['asu']):
        if s in finalTS['ReservoirLevels'] and any(finalTS['ReservoirLevels'][s] > 0) and all(
                finalTS['ReservoirLevels'][s] - 1 <= 1e-11):
            # get the time series
            parameters['StorageProfile']['val'][i, :] = finalTS['ReservoirLevels'][s][idx_sim].values
        elif s in finalTS['ReservoirLevels'] and any(finalTS['ReservoirLevels'][s] > 0) and any(
                finalTS['ReservoirLevels'][s] - 1 > 1e-11):
            logging.critical(s + ': The reservoir level is sometimes higher than its capacity (>1) !')
            sys.exit(1)
        else:
            logging.warning(
                'Could not find reservoir level data for storage plant ' + s + '. Using the provided default initial '
                                                                               'and final values')
            parameters['StorageProfile']['val'][i, :] = np.linspace(config['default']['ReservoirLevelInitial'],
                                                                    config['default']['ReservoirLevelFinal'],
                                                                    len(idx_sim))
        # The initial level is the same as the first value of the profile:
        if s in Plants_sto.index:
            parameters['StorageInitial']['val'][i] = parameters['StorageProfile']['val'][i, 0] * \
                                                     finalTS['AvailabilityFactors'][s][idx_sim[0]] * \
                                                     Plants_sto['StorageCapacity'][s] * Plants_sto['Nunits'][s]
        if s in Plants_thms.index:
            parameters['StorageInitial']['val'][i] = parameters['StorageProfile']['val'][i, 0] * \
                                                     finalTS['AvailabilityFactors'][s][idx_sim[0]] * \
                                                     Plants_thms['StorageCapacity'][s] * Plants_thms['Nunits'][s]

    # Storage Inflows:
    for i, s in enumerate(sets['s']):
        if s in finalTS['ScaledInflows']:
            parameters['StorageInflow']['val'][i, :] = finalTS['ScaledInflows'][s][idx_sim].values * \
                                                       Plants_sto['PowerCapacity'][s]

    # Heat demands:
    for i, u in enumerate(sets['n_th']):
        if u in finalTS['HeatDemand']:
            parameters['HeatDemand']['val'][i, :] = finalTS['HeatDemand'][u][idx_sim].values
            parameters['CostHeatSlack']['val'][i, :] = finalTS['CostHeatSlack'][u][idx_sim].values

    # H2 time series:
    for i, u in enumerate(sets['s']):
        if u in finalTS['H2RigidDemand']:
            parameters['H2Demand']['val'][i, :] = finalTS['H2RigidDemand'][u][idx_sim].values
    for i, u in enumerate(sets['p2h2']):
        if u in finalTS['H2RigidDemand']:
            parameters['CostH2Slack']['val'][i, :] = finalTS['CostH2Slack'][u][idx_sim].values
        if u in finalTS['H2FlexibleDemand']:
            parameters['PtLDemandInput']['val'][i, :] = finalTS['H2FlexibleDemand'][u][idx_sim].values
    if 'H2FlexibleCapacity' in config and config['H2FlexibleCapacity'] != '':
        check_PtLDemand(parameters, config)

    # Ramping rates are reconstructed for the non dimensional value provided
    # (start-up and normal ramping are not differentiated)
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
    if len(finalTS['AvailabilityFactors'].columns) != 0:
        for i, u in enumerate(sets['au']):
            if u in finalTS['AvailabilityFactors'].columns:
                parameters['AvailabilityFactor']['val'][i, :] = finalTS['AvailabilityFactors'][u].values

    # Efficiencies (currently limited to p2h units, but can be extended to all units):
    if len(finalTS['Efficiencies']) != 0:
        for i, u in enumerate(sets['p2h']):
            parameters['Efficiency']['val'][i, :] = finalTS['Efficiencies'][u].values

    values = np.ndarray([len(sets['mk']), len(sets['n']), len(sets['h'])])
    for i in range(len(sets['n'])):
        values[0, i, :] = finalTS['Load'][sets['n'][i]]
        values[1, i, :] = finalTS['Reserve2U'][sets['n'][i]]
        values[2, i, :] = finalTS['Reserve2D'][sets['n'][i]]
        values[3, i, :] = finalTS['Load'][sets['n'][i]] * finalTS['ShareOfFlexibleDemand'][sets['n'][i]]
    parameters['Demand'] = {'sets': sets_param['Demand'], 'val': values}

    # Emission Rate:
    parameters['EmissionRate']['val'][:, 0] = Plants_merged['EmissionRate'].values

    # Load Shedding:
    for i, c in enumerate(sets['n']):
        parameters['LoadShedding']['val'][i] = finalTS['LoadShedding'][c] * PeakLoad[c]
        parameters['CostLoadShedding']['val'][i] = finalTS['CostLoadShedding'][c]
        parameters['CostCurtailment']['val'][i] = finalTS['CostCurtailment'][c]

    # %%###############################################################################################################
    # Variable Cost
    # Equivalence dictionary between fuel types and price entries in the config sheet:
    FuelEntries = {'BIO': 'PriceOfBiomass', 'GAS': 'PriceOfGas', 'HRD': 'PriceOfBlackCoal', 'LIG': 'PriceOfLignite',
                   'NUC': 'PriceOfNuclear', 'OIL': 'PriceOfFuelOil', 'PEA': 'PriceOfPeat'}
    for unit in range(Nunits):
        c = Plants_merged['Zone'][unit]  # zone to which the unit belongs
        found = False
        for FuelEntry in FuelEntries:
            if Plants_merged['Fuel'][unit] == FuelEntry:
                if Plants_merged['Technology'][unit] == 'ABHP':
                    parameters['CostVariable']['val'][unit, :] = FuelPrices[FuelEntries[FuelEntry]][c] / \
                                                                 Plants_merged['Efficiency'][unit] + \
                                                                 Plants_merged['EmissionRate'][unit] * \
                                                                 FuelPrices['PriceOfCO2'][c]
                    found = True
                else:
                    parameters['CostVariable']['val'][unit, :] = FuelPrices[FuelEntries[FuelEntry]][c] / \
                                                                 Plants_merged['Efficiency'][unit] + \
                                                                 Plants_merged['EmissionRate'][unit] * \
                                                                 FuelPrices['PriceOfCO2'][c]
                    found = True
        # Special case for biomass plants, which are not included in EU ETS:
        if Plants_merged['Fuel'][unit] == 'BIO':
            parameters['CostVariable']['val'][unit, :] = FuelPrices['PriceOfBiomass'][c] / \
                                                         Plants_merged['Efficiency'][unit]
            found = True
        if not found:
            logging.warning('No fuel price value has been found for fuel ' + Plants_merged['Fuel'][unit] +
                            ' in unit ' + Plants_merged['Unit'][unit] + '. A null variable cost has been assigned')

    # %%###############################################################################################################

    # Maximum Line Capacity
    for i, l in enumerate(sets['l']):
        if l in NTCs.columns:
            parameters['FlowMaximum']['val'][i, :] = finalTS['NTCs'][l]
        if l in Inter_RoW.columns:
            parameters['FlowMaximum']['val'][i, :] = finalTS['Inter_RoW'][l]
            parameters['FlowMinimum']['val'][i, :] = finalTS['Inter_RoW'][l]
        parameters['PriceTransmission']['val'][i, :] = finalTS['PriceTransmission'][l]

    # Check values:
    check_MinMaxFlows(parameters['FlowMinimum']['val'], parameters['FlowMaximum']['val'])

    parameters['LineNode'] = incidence_matrix(sets, 'l', parameters, 'LineNode')

    # Outage Factors
    if len(finalTS['Outages'].columns) != 0:
        for i, u in enumerate(sets['au']):
            if u in finalTS['Outages'].columns:
                parameters['OutageFactor']['val'][i, :] = finalTS['Outages'][u].values
            else:
                logging.warning('Outages factors not found for unit ' + u + '. Assuming no outages')

    # Participation to the reserve market
    list_of_participating_units = []  # new list
    for unit in Plants_merged.index:
        tech = Plants_merged.loc[unit, 'Technology']
        if tech in config['ReserveParticipation'] and Plants_merged.loc[unit, 'CHPType'] == '':
            list_of_participating_units.append(
                unit)  # if unit same technology as allowed without CHP and unit is no CHP then add to list
        elif tech in config['ReserveParticipation_CHP'] and Plants_merged.loc[unit, 'CHPType'] != '':
            list_of_participating_units.append(
                unit)  # if unit same technology as allowed with CHP and unit is CHP then add to list

    values = np.array([s in list_of_participating_units for s in sets['au']],
                      dtype='bool')  # same as before but with new list
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
    for i in range(len(sets['n_th'])):
        parameters['Location_th']['val'][:, i] = (Plants_merged['Zone_th'] == zones_th[i]).values

    # CHPType parameter:
    sets['chp_type'] = ['Extraction', 'Back-Pressure', 'P2H']
    parameters['CHPType'] = define_parameter(['chp', 'chp_type'], sets, value=0)
    for i, u in enumerate(sets['chp']):
        if u in Plants_chp.index:
            if Plants_chp.loc[u, 'CHPType'].lower() == 'extraction':
                parameters['CHPType']['val'][i, 0] = 1
            elif Plants_chp.loc[u, 'CHPType'].lower() == 'back-pressure':
                parameters['CHPType']['val'][i, 1] = 1
            elif Plants_chp.loc[u, 'CHPType'].lower() == 'p2h':
                parameters['CHPType']['val'][i, 2] = 1
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
    sets['x_config'] = ['FirstDay', 'LastDay', 'RollingHorizon Length', 'RollingHorizon LookAhead',
                        'SimulationTimeStep', 'ValueOfLostLoad', 'QuickStartShare', 'CostOfSpillage', 'WaterValue',
                        'DemandFlexibility']
    sets['y_config'] = ['year', 'month', 'day', 'val']
    dd_begin = idx_sim[0]
    dd_end = idx_sim[-1]

    values = np.array(
        [[dd_begin.year, dd_begin.month, dd_begin.day, 0],
         [dd_end.year, dd_end.month, dd_end.day, 0],
         [1e-5, 0, config['HorizonLength'], 0],
         [1e-5, 0, config['LookAhead'], 0],
         [1e-5, 0, 0, config['SimulationTimeStep']],
         [1e-5, 0, 0, config['default']['ValueOfLostLoad']],
         [1e-5, 0, 0, config['default']['ShareOfQuickStartUnits']],
         [1e-5, 0, 0, config['default']['PriceOfSpillage']],
         [1e-5, 0, 0, config['default']['WaterValue']],
         [1e-5, 0, 0, config['default']['DemandFlexibility']]]
    )
    # the 1E-5 values are needed to be sure that all sets are written with the config parameter
    parameters['Config'] = {'sets': ['x_config', 'y_config'], 'val': values}

    # %%################################################################################################################
    ######################################   Simulation Environment     ################################################
    ####################################################################################################################

    # Output folder:
    sim = config['SimulationDirectory']

    # Simulation data:
    SimData = {'sets': sets, 'parameters': parameters, 'config': config, 'units': Plants_merged,
               'geo': geo, 'version': dispa_version}

    # list_vars = []
    gdx_out = "Inputs.gdx"
    if config['WriteGDX']:
        write_variables(config, gdx_out, [sets, parameters])

    # if the sim variable was not defined:
    if 'sim' not in locals():
        logging.error('Please provide a path where to store the DispaSET inputs (in the "sim" variable)')
        sys.exit(1)

    if not os.path.exists(sim):
        os.makedirs(sim)

    if MTS:
        if not LP:
            logging.error('Simulation in MTS must be LP')
            sys.exit(1)
        else:
            fin = open(os.path.join(GMS_FOLDER, 'UCM_h.gms'))
            fout = open(os.path.join(sim, 'UCM_h.gms'), "wt")
            for line in fin:
                line = line.replace('$setglobal LPFormulation 0', '$setglobal LPFormulation 1')
                line = line.replace('$setglobal MTS 0', '$setglobal MTS 1')
                fout.write(line)
            fin.close()
            fout.close()
            # additionally allso copy UCM_h_simple.gms
            shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'),
                            os.path.join(sim, 'UCM_h_simple.gms'))

    elif LP:
        fin = open(os.path.join(GMS_FOLDER, 'UCM_h.gms'))
        fout = open(os.path.join(sim, 'UCM_h.gms'), "wt")
        for line in fin:
            line = line.replace('$setglobal LPFormulation 0', '$setglobal LPFormulation 1')
            fout.write(line)
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
    cplex_options = {'epgap': 0.005,  # TODO: For the moment hardcoded, it has to be moved to a config file
                     'numericalemphasis': 0,
                     'scaind': 1,
                     'lpmethod': 0,
                     'relaxfixedinfeas': 0,
                     'mipstart': 1,
                     'epint': 0,
                     'heuristiceffort': 2,
                     'lbheur': 1,
                     'probe': 1}

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

    # Remove previously-created debug files:
    debugfile = os.path.join(sim, 'debug.gdx')
    if os.path.isfile(debugfile):
        try:
            os.remove(debugfile)
        except OSError:
            print('Could not erase previous debug file ' + debugfile)

    return SimData
