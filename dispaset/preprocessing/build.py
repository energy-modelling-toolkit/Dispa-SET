import datetime as dt
import logging
import os
import shutil
import sys

import numpy as np
import pandas as pd
from future.builtins import int

from .data_check import check_units, check_sto, check_AvailabilityFactors, check_temperatures, check_clustering, \
    isStorage, check_chp, check_p2h, check_df, check_MinMaxFlows, \
    check_FlexibleDemand, check_reserves, check_p2bs, check_boundary_sector, \
    check_BSFlexDemand, check_BSFlexSupply
from .data_handler import NodeBasedTable, load_time_series, UnitBasedTable, merge_series, define_parameter, \
    load_geo_data, GenericTable
from .reserves import percentage_reserve, probabilistic_reserve, generic_reserve
from .utils import select_units, interconnections, clustering, EfficiencyTimeSeries, \
    BoundarySectorEfficiencyTimeSeries, incidence_matrix, pd_timestep
from .. import __version__
from ..common import commons
from ..misc.gdx_handler import write_variables
from difflib import SequenceMatcher

GMS_FOLDER = os.path.join(os.path.dirname(__file__), '..', 'GAMS')


def build_single_run(config, profiles=None, PtLDemand=None, SectorXFlexDemand=None, SectorXFlexSupply=None, MTS=0,
                     profilesSectorX=None):
    """
    This function reads the DispaSET config, loads the specified data,
    processes it when needed, and formats it in the proper DispaSET format.
    The output of the function is a directory with all inputs and simulation files required to run a DispaSET simulation

    :param config:        Dictionary with all the configuration fields loaded from the excel file. Output of the
                          'LoadConfig' function.
    :param profiles:      Profiles from mid term scheduling simulations
    :param PtLDemand:     Profiles for PtL demand from mid term scheduling simulations
    :param SectorXFlexDemand:  Profiles for BS Flex demand from mid term scheduling simulations
    :param SectorXFlexDemand:  Profiles for BS Flex supply from mid term scheduling simulations
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
    
    # Boolean variable to check wether it is NTC or DC-POWERFLOW:
    grid_flag = config.get('TransmissionGridType', "")  # If key does not exist it returns ""

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
    
    # PTDF Matrix:
    if (grid_flag == "DC-Power Flow"):
        if os.path.isfile(config['PTDFMatrix']):
            PTDF = pd.read_csv(config['PTDFMatrix'],
                               na_values=commons['na_values'],
                               keep_default_na=False, index_col=0)
        else:
            logging.warning('No PTDF Matrix will be considered (no valid file provided)')
            PTDF = pd.DataFrame()

    # Boundary Sector Interconnections:
    if os.path.isfile(config['BoundarySectorInterconnections']):
        BS_flows = load_time_series(config, config['BoundarySectorInterconnections']).fillna(0)
    else:
        logging.warning('No historical boundary sector flows will be considered (no valid file provided)')
        BS_flows = pd.DataFrame(index=config['idx_long'])
    if os.path.isfile(config['BoundarySectorNTC']):
        BS_NTC = load_time_series(config, config['BoundarySectorNTC']).fillna(0)
    else:
        logging.warning('No boundary sector NTC values will be considered (no valid file provided)')
        BS_NTC = pd.DataFrame(index=config['idx_long'])

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

    # Boundary Sectors:
    BoundarySector = pd.DataFrame()
    if os.path.isfile(config['BoundarySectorData']):
        BoundarySector = pd.read_csv(config['BoundarySectorData'],
                                     na_values=commons['na_values'],
                                     keep_default_na=False, index_col='Sector')
    else:
        for key in ['SectorXStorageCapacity', 'SectorXStorageSelfDischarge']:
            BoundarySector[key] = np.nan

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
            # plants = plants.append(tmp, ignore_index=True, sort=False)
            plants = pd.concat([plants, tmp], ignore_index=True, sort=False)
    # remove invalid power plants:
    plants = select_units(plants, config)
    filter_col = [col for col in plants if col.startswith('Sector')]
    plants[filter_col] = plants[filter_col].astype(str)

    # Some columns can be in two format (absolute or per unit). If not specified, they are set to zero:
    for key in ['StartUpCost', 'NoLoadCost']:
        if key in plants:
            if key + '_pu' in plants:
                logging.warning(
                    'Column ' + key + '_pu found in the power plant table but not used since column ' + key + ' exists')
        elif key + '_pu' in plants:
            plants[key] = plants[key + '_pu'] * plants['PowerCapacity']
        else:
            plants[key] = 0

    # If not present, add the non-compulsory fields to the units table:
    for key in ['CHPPowerLossFactor', 'CHPPowerToHeat', 'CHPType', 'STOCapacity', 'STOSelfDischarge',
                'STOMaxChargingPower', 'STOChargingEfficiency', 'CHPMaxHeat', 'WaterWithdrawal',
                'WaterConsumption', 'Sector1', 'Sector2']:
        if key not in plants.columns:
            plants[key] = np.nan

    # If the thermal zones are not defined in the units table, define one individual zone per power plant:
    for key in ['Zone_th']:
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
    # plants_all_sto = plants_sto.append(plants_thms)
    plants_all_sto = pd.concat([plants_sto, plants_thms])

    # Defining the heat only units:
    # plants_heat = plants[[u in commons['tech_heat'] for u in plants['Technology']]]
    # check_heat(config, plants_heat)

    # Defining the boundary sector only units:
    plants_bs = plants[[u in commons['tech_boundary_sector'] for u in plants['Technology']]]
    check_boundary_sector(config, plants_bs)

    # Defining the CHPs:
    plants_chp = plants[[str(x).lower() in commons['types_CHP'] for x in plants['CHPType']]]
    check_chp(config, plants_chp)

    # Defining the P2H units:
    plants_p2h = plants[[u in commons['tech_p2ht'] for u in plants['Technology']]]
    check_p2h(config, plants_p2h)

    # All heating units:
    # plants_heat = pd.concat([plants_heat, plants_chp])
    # plants_heat = pd.concat([plants_heat, plants_p2h])

    # Defining the P2BS units:
    # TODO: Check if plants should be grouped by technology or by energy in one of the boundary sectors
    plants_p2bs = plants[[u in commons['tech_p2bs'] for u in plants['Technology']]]
    check_p2bs(config, plants_p2bs)

    # Define all Boundary Sector units:
    plants_all_bs = pd.concat([plants_p2bs, plants_bs])
    plants_all_bs = pd.concat([plants_all_bs, plants_p2h])
    plants_all_bs = pd.concat([plants_all_bs, plants_chp])
    plants_all_bs = pd.concat([plants_all_bs, plants_thms])

    Outages = UnitBasedTable(plants, 'Outages', config, fallbacks=['Unit', 'Technology'])
    AF = UnitBasedTable(plants, 'RenewablesAF', config, fallbacks=['Unit', 'Technology'], default=1,
                        RestrictWarning=commons['tech_renewables'])
    AF = AF.apply(pd.to_numeric)

    ReservoirLevels = UnitBasedTable(plants_all_sto, 'ReservoirLevels', config,
                                     fallbacks=['Unit', 'Technology', 'Zone'],
                                     default=0)
    ReservoirScaledInflows = UnitBasedTable(plants_all_sto, 'ReservoirScaledInflows', config,
                                            fallbacks=['Unit', 'Technology', 'Zone'], default=0)
    ReservoirScaledOutflows = UnitBasedTable(plants_all_sto, 'ReservoirScaledOutflows', config,
                                             fallbacks=['Unit', 'Technology', 'Zone'], default=0)
    Temperatures = NodeBasedTable('Temperatures', config)

    # Detecting boundary zones:
    filter_bs_cols = [col for col in plants_all_bs if col.startswith('Sector')]
    zones_bs = []
    for col in plants_all_bs[filter_bs_cols].columns:
        tmp = plants_all_bs[col].unique().tolist()
        zones_bs += tmp
    zones_bs = list(set(zones_bs))
    if '' in zones_bs:
        zones_bs.remove('')
    if np.nan in zones_bs:
        zones_bs = [x for x in zones_bs if pd.isnull(x) == False]
    if 'nan' in zones_bs:
        zones_bs.remove('nan')
    SectorXDemand = GenericTable(zones_bs, 'SectorXDemand', config, default=0)
    CostXNotServed = GenericTable(zones_bs, 'CostXNotServed', config,
                                  default=config['default']['CostXNotServed'])

    BoundarySector = BoundarySector.reindex(zones_bs)
    BoundarySector.fillna(0, inplace=True)
    SectorXReservoirLevels = GenericTable(zones_bs, 'SectorXReservoirLevels', config, default=0)

    # Boundary Sector Max Spillage
    if os.path.isfile(config['BoundarySectorMaxSpillage']):
        BS_spillage = load_time_series(config, config['BoundarySectorMaxSpillage']).fillna(0)
    else:
        BS_spillage = pd.DataFrame(index=idx_long)
        logging.warning('No maximum spillage capacity provided.')
        # sys.exit('No maximum spillage capacity provided. This parameter is necessary.')

    # Boundary Sector Forced Spillage
    if os.path.isfile(config['BoundarySectorMaxSpillage']):
        BS_forced_spillage = pd.DataFrame(0, index=BS_spillage.index, columns=BS_spillage.columns)
    else:
        BS_forced_spillage = pd.DataFrame(index=idx_long)

    # Read BS Flexible demand & supply
    SectorXFlexibleDemand = GenericTable(zones_bs, 'SectorXFlexibleDemand', config, default=0)
    SectorXFlexibleSupply = GenericTable(zones_bs, 'SectorXFlexibleSupply', config, default=0)

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

    if profilesSectorX is not None:
        for key in profilesSectorX.columns:
            if key not in SectorXReservoirLevels.columns:
                logging.warning('The reservoir profile "' + key + '" provided by the MTS is not found in the '
                                                                  'ReservoirLevels table')
            else:
                SectorXReservoirLevels[key].update(profilesSectorX[key])
                logging.info(
                    'The reservoir profile "' + key + '" provided by the MTS is used as target reservoir level')

    # Update SectorXFlexDemand (SectorXFlexibleDemand with demand from mid term scheduling)
    if SectorXFlexDemand is not None and any(SectorXFlexibleDemand) > 0:
        for key in SectorXFlexDemand.columns:
            if key not in SectorXFlexibleDemand.columns:
                logging.warning('The BS flexible demand "' + key + '" provided by the MTS is not found in the '
                                                                   'SectorXFlexibleDemand table')
            else:
                SectorXFlexibleDemand[key].update(SectorXFlexDemand[key])

    # Update SectorXFlexSupply (SectorXFlexibleSupply with demand from mid term scheduling)
    if SectorXFlexSupply is not None and any(SectorXFlexibleSupply) > 0:
        for key in SectorXFlexSupply.columns:
            if key not in SectorXFlexibleSupply.columns:
                logging.warning('The BS flexible demand "' + key + '" provided by the MTS is not found in the '
                                                                   'SectorXFlexibleSupply table')
            else:
                SectorXFlexibleSupply[key].update(SectorXFlexSupply[key])

    # data checks:
    check_AvailabilityFactors(plants, AF)
    # check_heat_demand(plants, HeatDemand, zones_th)
    # TODO: Check bounadry sector demands
    check_temperatures(plants, Temperatures)
    check_FlexibleDemand(ShareOfFlexibleDemand)
    check_reserves(Reserve2D, Reserve2U, Load)

    # Fuel prices:
    fuels = ['PriceOfNuclear', 'PriceOfBlackCoal', 'PriceOfGas', 'PriceOfFuelOil', 'PriceOfBiomass', 'PriceOfCO2',
             'PriceOfLignite', 'PriceOfPeat', 'PriceOfAmmonia']
    FuelPrices = {}
    for fuel in fuels:
        FuelPrices[fuel] = NodeBasedTable(fuel, config, default=config['default'][fuel])

    # Bidirectional Interconnections:
    if (grid_flag == "DC-Power Flow"):
        bidirectional_connections = NTC
        data_dict = {'connection1':[],'connection2':[],'ratio':[]}
        aux = bidirectional_connections
        for i in bidirectional_connections:
            aux = aux.drop([i], axis=1)
            for j in aux:
                ratio = SequenceMatcher(None, j, i).quick_ratio()
                data_dict['connection1'].append(i)
                data_dict['connection2'].append(j)
                data_dict['ratio'].append(ratio)
        df_ratio = pd.DataFrame(data_dict)
        eqratio = df_ratio['ratio'] == 1
        reptd_con = df_ratio[eqratio]
        reptd_con.reset_index(level =0, inplace = True)
        reptd_con = reptd_con.drop(['index','ratio'], axis=1)
        NTC_new = pd.DataFrame()
        for i in range(len(reptd_con)):
            ope = reptd_con.iloc[i]
            max_val = (NTC[ope].max()).max()
            NTC_aux = NTC[reptd_con.loc[i,'connection1']]
            NTC_aux.values[:] = max_val
            NTC_new = pd.concat([NTC_new,NTC_aux], axis=1)
            for j in ope:
                NTC = NTC.drop([j], axis=1)
            NTC = pd.concat([NTC,NTC_new], axis=1)
        NTC.index = NTC.index.astype('datetime64[ns]') 
               
        # Check the PTDF Matrix (TODO: check the order or sequence):   
        for key in NTC.columns:
            if key not in PTDF.index:
                logging.warning('The transmission line' + str(key) + ' is not in the provided PTDF')
                logging.error('The PTDF Matrix provided is not valid')
                sys.exit(1)
        for key in config['zones']:
            if key not in PTDF.columns:
                logging.warning('The node' + str(key) + ' is not in the provided zones')
                logging.error('The PTDF Matrix provided is not valid')
                sys.exit(1)
    
    # Interconnections:
    [Interconnections_sim, Interconnections_RoW, Interconnections] = interconnections(config['zones'], NTC, flows)

    # Boundary Sector Interconnections:
    [BSInterconnections_sim, BSInterconnections_RoW, BSInterconnections] = interconnections(zones_bs, BS_NTC, BS_flows)

    # Boundary Sector Spillage:
    [BSSpillage_sim, BSSpillage_RoW, BSSpillage] = interconnections(zones_bs, BS_spillage, BS_forced_spillage)

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

    if len(BSInterconnections_sim.columns) > 0:
        BS_NTCs = BSInterconnections_sim.reindex(config['idx_long'])
    else:
        BS_NTCs = pd.DataFrame(index=config['idx_long'])
    BS_Inter_RoW = BSInterconnections_RoW.reindex(config['idx_long'])

    if len(BSSpillage_sim.columns) > 0:
        BS_Spillages = BSSpillage_sim.reindex(config['idx_long'])
    else:
        BS_Spillages = pd.DataFrame(index=config['idx_long'])
    BS_Spillage_RoW = BSSpillage_RoW.reindex(config['idx_long'])

    # Clustering of the plants:
    Plants_merged, mapping = clustering(plants, method=config['SimulationType'])
    # Check clustering:
    check_clustering(plants, Plants_merged)

    # Renaming the columns to ease the production of parameters:
    BoundarySector.rename(columns={'STOCapacity': 'SectorXStorageCapacity',
                                   'STOSelfDischarge': 'SectorXStorageSelfDischarge',
                                   # 'STOMaxChargingPower': 'BoundarySectorStorageChargingCapacity',
                                   'STOMinSOC': 'SectorXStorageMinimum'}, inplace=True)

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

    filter_col = [col for col in Plants_merged if col.startswith('Sector')]
    Plants_merged[filter_col] = Plants_merged[filter_col].astype(str)

    # Defining the hydro storages:
    Plants_sto = Plants_merged[[u in commons['tech_storage'] for u in Plants_merged['Technology']]]
    # Defining the thermal storages:
    Plants_thms = Plants_merged[[u in commons['tech_thermal_storage'] for u in Plants_merged['Technology']]]
    # Defining all storage units:
    Plants_all_sto = Plants_merged[[u in [x for x in commons['Technologies'] if x in commons['tech_storage'] +
                                          commons['tech_thermal_storage']] for u in
                                    Plants_merged['Technology']]]
    # check storage plants:
    check_sto(config, Plants_sto, raw_data=False)
    check_sto(config, Plants_thms, raw_data=False)
    check_sto(config, Plants_all_sto, raw_data=False)
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

    # # Filter heat only plants
    # Plants_heat_only = Plants_merged[[u in commons['tech_heat'] for u in Plants_merged['Technology']]].copy()
    # # Check heat only units
    # check_heat(config, Plants_heat_only)

    # Filter boundary sector only plants
    Plants_boundary_sector_only = Plants_merged[
        [u in commons['tech_boundary_sector'] for u in Plants_merged['Technology']]].copy()
    # Check heat only units
    check_boundary_sector(config, Plants_boundary_sector_only)

    # All heating plants:
    # Plants_heat = Plants_heat_only.copy()
    Plants_heat = Plants_chp.copy()
    Plants_heat = pd.concat([Plants_heat, Plants_p2h])
    Plants_heat = pd.concat([Plants_heat, Plants_thms])

    # Filter power to boundary sector plants
    Plants_p2bs = Plants_merged[[u in commons['tech_p2bs'] for u in Plants_merged['Technology']]].copy()
    # Check heat only units
    check_p2bs(config, Plants_p2bs)

    # Water storage
    Plants_wat = Plants_merged[(Plants_merged['Fuel'] == 'WAT') & (Plants_merged['Technology'] != 'HROR')].copy()

    # Calculating the efficiency time series for each unit:
    Efficiencies = EfficiencyTimeSeries(config, Plants_merged, Temperatures)

    # Calculating boundary sector efficiencies
    BoundarySectorEfficiencies = BoundarySectorEfficiencyTimeSeries(config, Plants_merged, zones_bs)

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
               'Efficiencies': Efficiencies,
               'NTCs': NTCs, 'Inter_RoW': Inter_RoW,
               'BS_NTCs': BS_NTCs, 'BS_Inter_RoW': BS_Inter_RoW,
               'EfficienciesBoundarySector': BoundarySectorEfficiencies['Efficiency'],
               'ChargingEfficienciesBoundarySector': BoundarySectorEfficiencies['ChargingEfficiency'],
               'LoadShedding': LoadShedding, 'CostLoadShedding': CostLoadShedding, 'CostCurtailment': CostCurtailment,
               'ScaledInflows': ReservoirScaledInflows, 'ScaledOutflows': ReservoirScaledOutflows,
               'ReservoirLevels': ReservoirLevels, 'Outages': Outages, 'AvailabilityFactors': AF,
               # 'CostHeatSlack': CostHeatSlack, 'HeatDemand': HeatDemand,
               'ShareOfFlexibleDemand': ShareOfFlexibleDemand,
               'PriceTransmission': PriceTransmission,
               'SectorXDemand': SectorXDemand, 'CostXNotServed': CostXNotServed,
               'SectorXFlexibleDemand': SectorXFlexibleDemand, 'SectorXFlexibleSupply': SectorXFlexibleSupply,
               'BSMaxSpillage': BS_Spillages, 'SectorXReservoirLevels': SectorXReservoirLevels}

    # Merge the following time series with weighted averages
    for key in ['ScaledInflows', 'ScaledOutflows', 'Outages', 'AvailabilityFactors']:
        finalTS[key] = merge_series(Plants_merged, plants, finalTS[key], tablename=key)

    # Merge the following time series by weighted average based on storage capacity
    for key in ['ReservoirLevels']:
        finalTS[key] = merge_series(Plants_merged, plants, finalTS[key], tablename=key, method='StorageWeightedAverage')

    # Check that all times series data is available with the specified data time step:
    for key in FuelPrices:
        check_df(FuelPrices[key], StartDate=config['idx'][0], StopDate=config['idx'][-1], name=key)
    for key in finalTS:
        if key in ['EfficienciesBoundarySector','ChargingEfficienciesBoundarySector']:
            for ts in finalTS[key]:
                check_df(finalTS[key][ts], StartDate=config['idx'][0], StopDate=config['idx'][-1], name=key)
        else:
            check_df(finalTS[key], StartDate=config['idx'][0], StopDate=config['idx'][-1], name=key)

    # Resemple to the required time step
    if config['DataTimeStep'] != config['SimulationTimeStep']:
        for key in FuelPrices:
            if len(FuelPrices[key].columns) > 0:
                FuelPrices[key] = FuelPrices[key].resample(pd_timestep(config['SimulationTimeStep'])).mean()
        for key in finalTS:
            if key in ['EfficienciesBoundarySector','ChargingEfficienciesBoundarySector']:
                for ts in finalTS[key]:
                    if len(finalTS[key][ts].columns) > 0:
                        finalTS[key][ts] = finalTS[key][ts].resample(pd_timestep(config['SimulationTimeStep'])).mean()
            else:
                if len(finalTS[key].columns) > 0:
                    finalTS[key] = finalTS[key].resample(pd_timestep(config['SimulationTimeStep'])).mean()

    # Formatting PTDF matrix according to sets['l'] and sets['n'] order or sequence
    if (grid_flag == "DC-Power Flow"):
        PTDF = PTDF.T
        PTDF_l = pd.DataFrame()
        for i in Interconnections:
            PTDF_l = pd.concat([PTDF_l,PTDF[i]], axis=1)
        PTDF = PTDF_l.T
        PTDF_n = pd.DataFrame()
        for j in config['zones']:
            PTDF_n = pd.concat([PTDF_n,PTDF[j]], axis=1)
        PTDF = PTDF_n
    
    # %%###############################################################################################################
    ############################################   Sets    ############################################################
    ###################################################################################################################

    # The sets are defined within a dictionary:
    sets = {}
    sets['h'] = [str(x + 1) for x in range(Nsim)]
    sets['z'] = [str(x + 1) for x in range(int(Nsim - config['LookAhead'] * 24 / config['SimulationTimeStep']))]
    sets['mk'] = ['DA', '2U', '2D', 'Flex']
    sets['n'] = config['zones']
    sets['nx'] = zones_bs
    sets['au'] = Plants_merged.index.tolist()
    sets['l'] = Interconnections
    sets['lx'] = BSInterconnections
    sets['slx'] = BSSpillage
    sets['f'] = commons['Fuels']
    sets['p'] = ['CO2']
    sets['s'] = Plants_sto.index.tolist()
    sets['u'] = Plants_merged[[u in [x for x in commons['Technologies'] if x not in
                                     commons['tech_p2ht'] + commons['tech_thermal_storage'] +
                                     commons['tech_p2bs'] + commons['tech_boundary_sector']]
                               for u in Plants_merged['Technology']]].index.tolist()
    sets['chp'] = Plants_chp.index.tolist()
    sets['p2h'] = Plants_p2h.index.tolist()
    sets['p2x'] = Plants_p2bs.index.tolist()
    sets['th'] = Plants_heat.index.tolist()
    sets['thms'] = Plants_thms.index.tolist()
    sets['t'] = commons['Technologies']
    sets['tr'] = commons['tech_renewables']
    sets['tc'] = list(set(commons['Technologies']) - set(commons['tech_renewables']) - set(commons['tech_p2ht']) -
                      set(commons['tech_thermal_storage']))
    sets['wat'] = Plants_wat.index.tolist()
    sets['xu'] = Plants_boundary_sector_only.index.tolist()
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
    sets_param['SectorXStorageFinalMin'] = ['nx']
    sets_param['CostXNotServed'] = ['nx', 'h']
    sets_param['CostLoadShedding'] = ['n', 'h']
    sets_param['CostRampUp'] = ['au']
    sets_param['CostRampDown'] = ['au']
    sets_param['CostShutDown'] = ['au']
    sets_param['CostStartUp'] = ['au']
    sets_param['CostVariable'] = ['au', 'h']
    sets_param['Curtailment'] = ['n']
    sets_param['CostCurtailment'] = ['n', 'h']
    sets_param['Demand'] = ['mk', 'n', 'h']
    sets_param['Efficiency'] = ['au', 'h']
    sets_param['X2PowerConversionMultiplier'] = ['nx', 'au', 'h']
    sets_param['Power2XConversionMultiplier'] = ['nx', 'au', 'h']
    sets_param['EmissionMaximum'] = ['n', 'p']
    sets_param['EmissionRate'] = ['au', 'p']
    sets_param['FlowMaximum'] = ['l', 'h']
    sets_param['FlowMinimum'] = ['l', 'h']
    sets_param['FlowXMaximum'] = ['lx', 'h']
    sets_param['FlowXMinimum'] = ['lx', 'h']
    sets_param['Fuel'] = ['au', 'f']
    sets_param['SectorXDemand'] = ['nx', 'h']
    sets_param['LineNode'] = ['l', 'n']
    sets_param['LineXNode'] = ['lx', 'nx']
    sets_param['SectorXSpillageNode'] = ['slx', 'nx']
    sets_param['SectorXMaximumSpillage'] = ['slx', 'h']
    sets_param['LoadShedding'] = ['n', 'h']
    sets_param['Location'] = ['au', 'n']
    sets_param['LocationX'] = ['au', 'nx']
    sets_param['Markup'] = ['au', 'h']
    sets_param['Nunits'] = ['au']
    sets_param['OutageFactor'] = ['au', 'h']
    sets_param['PartLoadMin'] = ['au']
    sets_param['PowerCapacity'] = ['au']
    sets_param['PowerCapacityBoundarySector'] = ['au', 'nx']
    sets_param['PowerInitial'] = ['u']
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
    sets_param['StorageInflow'] = ['asu', 'h']
    sets_param['StorageInitial'] = ['asu']
    sets_param['StorageMinimum'] = ['asu']
    sets_param['StorageOutflow'] = ['asu', 'h']
    sets_param['StorageProfile'] = ['asu', 'h']
    sets_param['Technology'] = ['au', 't']
    sets_param['TimeUpMinimum'] = ['au']
    sets_param['TimeDownMinimum'] = ['au']
    sets_param['SectorXFlexDemandInput'] = ['nx', 'h']
    sets_param['SectorXFlexDemandInputInitial'] = ['nx']
    sets_param['SectorXFlexMaxCapacity'] = ['nx']
    sets_param['SectorXFlexSupplyInput'] = ['nx', 'h']
    sets_param['SectorXFlexSupplyInputInitial'] = ['nx']
    sets_param['SectorXFlexMaxSupply'] = ['nx']
    sets_param['SectorXStorageCapacity'] = ['nx']
    sets_param['SectorXStorageSelfDischarge'] = ['nx']
    sets_param['SectorXStorageMinimum'] = ['nx']
    sets_param['SectorXStorageInitial'] = ['nx']
    sets_param['SectorXStorageProfile'] = ['nx', 'h']
    if (grid_flag == "DC-Power Flow"):
        sets_param['PTDF'] = ['l', 'n']

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
    # for var in ['Technology', 'Fuel', 'Reserve', 'Location', 'Location_th', 'LocationX']:
    for var in ['Technology', 'Fuel', 'Reserve', 'Location', 'LocationX']:
        parameters[var] = define_parameter(sets_param[var], sets, value='bool')
    # for var in [col for col in plants if col.startswith('Sector')]:
    #     parameters[var] = define_parameter(sets_param[var], sets, value='bool')

    # %%
    # List of parameters whose value is known and provided in the dataframe BoundarySector
    for var in ['SectorXStorageCapacity', 'SectorXStorageSelfDischarge']:
        parameters[var]['val'] = BoundarySector[var].values

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

    # Particular treatment of SectorXFlexMaxCapacity that is not a time-series and that is given from the BS Inputs database
    if 'SectorXFlexibleDemand' in config and config['SectorXFlexibleDemand'] != '':
        for i, u in enumerate(sets['nx']):
            if u in BoundarySector.index:
                parameters['SectorXFlexMaxCapacity']['val'][i] = BoundarySector.loc[u, 'MaxFlexDemand']

    # Particular treatment of SectorXFlexMaxCapacity that is not a time-series and that is given from the BS Inputs database
    if 'SectorXFlexibleSupply' in config and config['SectorXFlexibleSupply'] != '':
        for i, u in enumerate(sets['nx']):
            if u in BoundarySector.index:
                parameters['SectorXFlexMaxSupply']['val'][i] = BoundarySector.loc[u, 'MaxFlexSupply']

    # Particular treatment of SectorXStorageMinimum:
    for i, u in enumerate(sets['nx']):
        if u in BoundarySector.index:
            parameters['SectorXStorageMinimum']['val'][i] = BoundarySector.loc[
                                                                u, 'SectorXStorageMinimum'] * \
                                                            BoundarySector.loc[
                                                                u, 'SectorXStorageCapacity']

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

    for i, nx in enumerate(sets['nx']):
        if nx in finalTS['SectorXReservoirLevels'] and any(finalTS['SectorXReservoirLevels'][nx] > 0) and all(
                finalTS['SectorXReservoirLevels'][nx] - 1 <= 1e-11):
            # get the time series
            parameters['SectorXStorageProfile']['val'][i, :] = finalTS['SectorXReservoirLevels'][nx][idx_sim].values
        elif nx in finalTS['SectorXReservoirLevels'] and any(finalTS['SectorXReservoirLevels'][nx] > 0) and any(
                finalTS['SectorXReservoirLevels'][nx] - 1 > 1e-11):
            logging.critical(nx + ': The reservoir level is sometimes higher than its capacity (>1) !')
            sys.exit(1)
        else:
            logging.warning(
                'Could not find reservoir level data for storage plant ' + nx + '. Using the provided default initial '
                                                                                'and final values')
            parameters['SectorXStorageProfile']['val'][i, :] = np.linspace(config['default']['ReservoirLevelInitial'],
                                                                           config['default']['ReservoirLevelFinal'],
                                                                           len(idx_sim))
        parameters['SectorXStorageInitial']['val'][i] = parameters['SectorXStorageProfile']['val'][i, 0] * \
                                                        BoundarySector['SectorXStorageCapacity'][nx]

    # Storage Inflows:
    for i, s in enumerate(sets['asu']):
        if s in finalTS['ScaledInflows']:
            parameters['StorageInflow']['val'][i, :] = finalTS['ScaledInflows'][s][idx_sim].values * \
                                                       Plants_all_sto['PowerCapacity'][s]
            if s in finalTS['ScaledOutflows']:
                parameters['StorageOutflow']['val'][i, :] = finalTS['ScaledOutflows'][s][idx_sim].values * \
                                                            Plants_all_sto['PowerCapacity'][s]

    # # Heat demands:
    # for i, u in enumerate(sets['n_th']):
    #     if u in finalTS['HeatDemand']:
    #         parameters['HeatDemand']['val'][i, :] = finalTS['HeatDemand'][u][idx_sim].values
    #         parameters['CostHeatSlack']['val'][i, :] = finalTS['CostHeatSlack'][u][idx_sim].values

    # Boundary sector demands:
    for i, u in enumerate(sets['nx']):
        if u in finalTS['SectorXDemand']:
            parameters['SectorXDemand']['val'][i, :] = finalTS['SectorXDemand'][u][idx_sim].values

        parameters['CostXNotServed']['val'][i, :] = finalTS['CostXNotServed'][u][idx_sim].values

    # Boundary Sector time series
    for i, u in enumerate(sets['nx']):
        if u in finalTS['SectorXFlexibleDemand']:
            parameters['SectorXFlexDemandInput']['val'][i, :] = finalTS['SectorXFlexibleDemand'][u][idx_sim].values
        if u in finalTS['SectorXFlexibleSupply']:
            parameters['SectorXFlexSupplyInput']['val'][i, :] = finalTS['SectorXFlexibleSupply'][u][idx_sim].values

    if 'SectorXFlexibleDemand' in config:
        check_BSFlexDemand(parameters, config)

    if 'SectorXFlexibleSupply' in config:
        check_BSFlexSupply(parameters, config)

    # Ramping rates are reconstructed for the non dimensional value provided
    # (start-up and normal ramping are not differentiated)
    parameters['RampUpMaximum']['val'] = Plants_merged['RampUpRate'].values * Plants_merged['PowerCapacity'].values * 60
    parameters['RampDownMaximum']['val'] = Plants_merged['RampDownRate'].values * Plants_merged[
        'PowerCapacity'].values * 60
    parameters['RampStartUpMaximum']['val'] = Plants_merged['RampUpRate'].values * Plants_merged[
        'PowerCapacity'].values * 60 + Plants_merged['PartLoadMin'].values * Plants_merged['PowerCapacity'].values
    parameters['RampShutDownMaximum']['val'] = Plants_merged['RampDownRate'].values * Plants_merged[
        'PowerCapacity'].values * 60 + Plants_merged['PartLoadMin'].values * Plants_merged['PowerCapacity'].values

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
        for i, u in enumerate(sets['au']):
            parameters['Efficiency']['val'][i, :] = finalTS['Efficiencies'][u].values

    # Assign charging and discharging efficiencies for boundary sectors
    values = np.ndarray([len(sets['nx']), len(sets['au']), len(sets['h'])])
    for i in range(len(sets['nx'])):
        for u in range(len(sets['au'])):
            values[i, u, :] = finalTS['EfficienciesBoundarySector'][sets['nx'][i]][sets['au'][u]]
    parameters['X2PowerConversionMultiplier'] = {'sets': sets_param['X2PowerConversionMultiplier'], 'val': values}

    values = np.ndarray([len(sets['nx']), len(sets['au']), len(sets['h'])])
    for i in range(len(sets['nx'])):
        for u in range(len(sets['au'])):
            values[i, u, :] = finalTS['ChargingEfficienciesBoundarySector'][sets['nx'][i]][sets['au'][u]]
    parameters['Power2XConversionMultiplier'] = {'sets': sets_param['Power2XConversionMultiplier'], 'val': values}

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
                   'NUC': 'PriceOfNuclear', 'OIL': 'PriceOfFuelOil', 'PEA': 'PriceOfPeat', 'AMO': 'PriceOfAmmonia'}
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
    if (grid_flag == "NTC"):
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
        
    elif (grid_flag == "DC-Power Flow"):
        for i, l in enumerate(sets['l']):
            if l in NTCs.columns:
                parameters['FlowMaximum']['val'][i, :] = finalTS['NTCs'][l]
                parameters['FlowMinimum']['val'][i, :] = (finalTS['NTCs'][l])*-1
                if l in Inter_RoW.columns:
                    parameters['FlowMaximum']['val'][i, :] = finalTS['Inter_RoW'][l]
                    parameters['FlowMinimum']['val'][i, :] = finalTS['Inter_RoW'][l]
                    parameters['PriceTransmission']['val'][i, :] = finalTS['PriceTransmission'][l]
    
        # Check values:
        check_MinMaxFlows(parameters['FlowMinimum']['val'], parameters['FlowMaximum']['val'])
    
        parameters['LineNode'] = incidence_matrix(sets, 'l', parameters, 'LineNode')        
        

    else:
        logging.error('Please provide a valid Transmission grid type')
        sys.exit(1)    

    # Maximum Boundary Sector Line Capacity
    for i, lx in enumerate(sets['lx']):
        if lx in BS_NTCs.columns:
            parameters['FlowXMaximum']['val'][i, :] = finalTS['BS_NTCs'][lx]
        if lx in BS_Inter_RoW.columns:
            parameters['FlowXMaximum']['val'][i, :] = finalTS['BS_Inter_RoW'][lx]
            parameters['FlowXMinimum']['val'][i, :] = finalTS['BS_Inter_RoW'][lx]
    for i, slx in enumerate(sets['slx']):
        if slx in BS_Spillages.columns:
            parameters['SectorXMaximumSpillage']['val'][i, :] = finalTS['BSMaxSpillage'][lx]

    # Check values:
    check_MinMaxFlows(parameters['FlowXMinimum']['val'], parameters['FlowXMaximum']['val'])

    parameters['LineXNode'] = incidence_matrix(sets, 'lx', parameters, 'LineXNode', nodes='nx')

    parameters['SectorXSpillageNode'] = incidence_matrix(sets, 'slx', parameters, 'SectorXSpillageNode', nodes='nx')

    # PTDF MATRIX 
    if (grid_flag == "DC-Power Flow"):
        if len(PTDF.columns) != 0:
            for i, u in enumerate(sets['n']):
                if u in PTDF.columns:
                    parameters['PTDF']['val'][:, i] = PTDF[u].values


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
    # for i in range(len(sets['n_th'])):
    #     parameters['Location_th']['val'][:, i] = (Plants_merged['Zone_th'] == zones_th[i]).values

    sectors = [col for col in plants if col.startswith('Sector')]
    for i in range(len(sets['nx'])):
        for s in sectors:
            parameters['LocationX']['val'][:, i] = np.logical_or(parameters['LocationX']['val'][:, i],
                                                                 (Plants_merged[s] == zones_bs[i]).values)

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
    if 'InitialPower' in Plants_merged.columns and Plants_merged['InitialPower'].notna().any():
        Plants_merged['InitialPower'].fillna(0, inplace=True)
    else:
        Plants_merged['InitialPower'] = finalTS['AvailabilityFactors'].iloc[0, :] * Plants_merged['PowerCapacity']

        for z in config['zones']:
            tmp_units = Plants_merged.loc[Plants_merged['Zone'] == z].copy()
            tmp_load = finalTS['Load'].iloc[0, :].loc[z]
            power_initial_res = tmp_units.loc[tmp_units['Fuel'].isin(['WAT', 'WIN', 'SUN']) &
                                              tmp_units['Technology'].isin(['HROR', 'WTON', 'WTOF', 'PHOT'])][
                'InitialPower'].sum()
            if tmp_load > power_initial_res:
                lack_of_power = tmp_load - power_initial_res
                for f in ['NUC', 'GAS', 'HRD', 'LIG', 'BIO', 'OIL']:
                    for t in ['COMC', 'STUR', 'GTUR', 'ICEN']:
                        n_units = tmp_units.loc[tmp_units['Fuel'].isin([f]) & tmp_units['Technology'].isin([t])].shape[
                            0]
                        power_initial_ft = tmp_units.loc[
                            tmp_units['Fuel'].isin([f]) & tmp_units['Technology'].isin([t]), 'InitialPower'].sum()
                        tmp_units.loc[:, 'Share'] = tmp_units.loc[
                                                        tmp_units['Fuel'].isin([f]) & tmp_units['Technology'].isin(
                                                            [t]), 'InitialPower'] / power_initial_ft
                        if lack_of_power - power_initial_ft > 0:
                            lack_of_power = lack_of_power - power_initial_ft
                        if lack_of_power - power_initial_ft < 0:
                            tmp_units.loc[tmp_units['Fuel'].isin([f]) & tmp_units['Technology'].isin(
                                [t]), 'InitialPower'] = lack_of_power * tmp_units.loc[
                                tmp_units['Fuel'].isin([f]) & tmp_units['Technology'].isin(
                                    [t]), 'Share']
                            lack_of_power = 0
                        if lack_of_power == 0:
                            tmp_units.loc[tmp_units['Fuel'].isin([f]) & tmp_units['Technology'].isin(
                                [t]), 'InitialPower'] = 0
            else:
                tmp_units.loc[~tmp_units['Fuel'].isin(['WAT', 'WIN', 'SUN']) &
                              ~tmp_units['Technology'].isin(['HROR', 'WTON', 'WTOF', 'PHOT']), 'InitialPower'] = 0
                lack_of_power = 0
            Plants_merged.update(tmp_units)

            if lack_of_power > 0:
                logging.error('In zone: ' + z + ' there is insufficient conventional + renewable ' +
                              'generation capacity of: ' + str(lack_of_power) +
                              '. If NTC + storage is not sufficient ShedLoad in ' + z +
                              ' is likely to occour. Check the inputs!')

        Plants_merged.loc[Plants_merged['Technology'].isin(commons['tech_renewables']), 'InitialPower'] = 0
        # Plants_merged.loc[Plants_merged['Technology'].isin(commons['tech_heat']), 'InitialPower'] = 0
        Plants_merged.loc[Plants_merged['Technology'].isin(commons['tech_p2ht']), 'InitialPower'] = 0
        Plants_merged.loc[Plants_merged['Technology'].isin(commons['tech_p2bs']), 'InitialPower'] = 0
        Plants_merged.loc[Plants_merged['Technology'].isin(commons['tech_thermal_storage']), 'InitialPower'] = 0

    if 'InitialPower' in Plants_merged:
        technologies = [x for x in commons['Technologies'] if x not in
                        commons['tech_p2ht'] + commons['tech_thermal_storage'] + commons['tech_p2bs'] +
                        commons['tech_boundary_sector']]
        parameters['PowerInitial']['val'] = Plants_merged.loc[Plants_merged['Technology'].isin(technologies)][
            'InitialPower'].values

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
    SimData = {'sets': sets, 'parameters': parameters, 'config': config, 'units_nonclustered': plants,
               'units': Plants_merged,
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

            if (grid_flag == "DC-Power Flow"):
                fin = open(os.path.join(GMS_FOLDER, 'UCM_h.gms'))
                fout = open(os.path.join(sim, 'UCM_h.gms'), "wt")
                logging.info('Simulation with DC-Power Flow')
                for line in fin:
                    line = line.replace('$setglobal LPFormulation 0', '$setglobal LPFormulation 1')
                    line = line.replace('$setglobal MTS 0', '$setglobal MTS 1')
                    line = line.replace('$setglobal TransmissionGrid 0', '$setglobal TransmissionGrid 1')
                    fout.write(line)
                fin.close()
                fout.close()
                
                # additionally allso copy UCM_h_simple.gms
                fin1 = open(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'))
                fout1 = open(os.path.join(sim, 'UCM_h_simple.gms'), "wt")
                logging.info('Simulation with DC-Power Flow')
                for line in fin1:
                    line = line.replace('$setglobal TransmissionGrid 0', '$setglobal TransmissionGrid 1')
                    fout1.write(line)
                fin1.close()
                fout1.close()
                # shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'),
                #                 os.path.join(sim, 'UCM_h_simple.gms'))
            elif (grid_flag == "NTC"):
                fin = open(os.path.join(GMS_FOLDER, 'UCM_h.gms'))
                fout = open(os.path.join(sim, 'UCM_h.gms'), "wt")
                logging.info('Simulation with NTC')
                for line in fin:
                    line = line.replace('$setglobal LPFormulation 0', '$setglobal LPFormulation 1')
                    line = line.replace('$setglobal MTS 0', '$setglobal MTS 1')
                    fout.write(line)
                fin.close()
                fout.close()
                # additionally allso copy UCM_h_simple.gms
                shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'),
                                os.path.join(sim, 'UCM_h_simple.gms'))
            else:
                logging.error('Please provide a valid Transmission grid type')
                sys.exit(1)

    elif LP:
            if (grid_flag == "DC-Power Flow"):
                fin = open(os.path.join(GMS_FOLDER, 'UCM_h.gms'))
                fout = open(os.path.join(sim, 'UCM_h.gms'), "wt")
                logging.info('Simulation with DC-Power Flow')
                for line in fin:
                    line = line.replace('$setglobal LPFormulation 0', '$setglobal LPFormulation 1')
                    line = line.replace('$setglobal TransmissionGrid 0', '$setglobal TransmissionGrid 1')
                    fout.write(line)
                fin.close()
                fout.close()
                
                # additionally allso copy UCM_h_simple.gms
                fin1 = open(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'))
                fout1 = open(os.path.join(sim, 'UCM_h_simple.gms'), "wt")
                logging.info('Simulation with DC-Power Flow')
                for line in fin1:
                    line = line.replace('$setglobal TransmissionGrid 0', '$setglobal TransmissionGrid 1')
                    fout1.write(line)
                fin1.close()
                fout1.close()
                # shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'),
                #                 os.path.join(sim, 'UCM_h_simple.gms'))
            elif (grid_flag == "NTC"):
                fin = open(os.path.join(GMS_FOLDER, 'UCM_h.gms'))
                fout = open(os.path.join(sim, 'UCM_h.gms'), "wt")
                logging.info('Simulation with NTC')
                for line in fin:
                    line = line.replace('$setglobal LPFormulation 0', '$setglobal LPFormulation 1')
                    fout.write(line)
                fin.close()
                fout.close()
                # additionally allso copy UCM_h_simple.gms
                shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'),
                                os.path.join(sim, 'UCM_h_simple.gms'))
            else:
                logging.error('Please provide a valid Transmission grid type')
                sys.exit(1)

    else:
            if (grid_flag == "DC-Power Flow"):
                fin = open(os.path.join(GMS_FOLDER, 'UCM_h.gms'))
                fout = open(os.path.join(sim, 'UCM_h.gms'), "wt")
                logging.info('Simulation with DC-Power Flow')
                for line in fin:
                    line = line.replace('$setglobal TransmissionGrid 0', '$setglobal TransmissionGrid 1')
                    fout.write(line)
                fin.close()
                fout.close()
                
                # additionally allso copy UCM_h_simple.gms
                fin1 = open(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'))
                fout1 = open(os.path.join(sim, 'UCM_h_simple.gms'), "wt")
                logging.info('Simulation with DC-Power Flow')
                for line in fin1:
                    line = line.replace('$setglobal TransmissionGrid 0', '$setglobal TransmissionGrid 1')
                    fout1.write(line)
                fin1.close()
                fout1.close()
                # shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'),
                #                 os.path.join(sim, 'UCM_h_simple.gms'))
            elif (grid_flag == "NTC"):
                fin = open(os.path.join(GMS_FOLDER, 'UCM_h.gms'))
                fout = open(os.path.join(sim, 'UCM_h.gms'), "wt")
                logging.info('Simulation with NTC')
                for line in fin:
                    line = line.replace('$setglobal LPFormulation 0', '$setglobal LPFormulation 1')
                    fout.write(line)
                fin.close()
                fout.close()
                # additionally allso copy UCM_h_simple.gms
                shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'),
                                os.path.join(sim, 'UCM_h_simple.gms'))
            else:
                logging.error('Please provide a valid Transmission grid type')
                sys.exit(1)
                
                
        # shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h.gms'),
        #                 os.path.join(sim, 'UCM_h.gms'))
        # # additionally allso copy UCM_h_simple.gms
        # shutil.copyfile(os.path.join(GMS_FOLDER, 'UCM_h_simple.gms'),
        #                 os.path.join(sim, 'UCM_h_simple.gms'))
    
    gmsfile = open(os.path.join(sim, 'UCM.gpr'), 'w')
    gmsfile.write(
        '[PROJECT] \n \n[RP:UCM_H] \n1= \n[OPENWINDOW_1] \nFILE0=UCM_h.gms \nFILE1=UCM_h.gms \nMAXIM=1 \nTOP=50 \nLEFT=50 \nHEIGHT=400 \nWIDTH=400')
    gmsfile.close()
    shutil.copyfile(os.path.join(GMS_FOLDER, 'writeresults.gms'),
                    os.path.join(sim, 'writeresults.gms'))
    # Create cplex option file
    if config['CplexSetting'] == '' and config['CplexAccuracy'] == '':
        cplex_options = {'epgap': 0.0005,
                         # TODO: For the moment hardcoded, it has to be moved to a config file
                         'numericalemphasis': 0,
                         'mipdisplay': 4,
                         'scaind': 1,
                         'lpmethod': 0,
                         'relaxfixedinfeas': 0,
                         'mipstart': 1,
                         'mircuts': 1,
                         'quality': True,
                         'bardisplay': 2,
                         'epint': 0,
                         'lbheur': 1,
                         }
        logging.info('Default Cplex setting used')
    elif config['CplexSetting'] == '' and config['CplexAccuracy'] != '':
        cplex_options = {'epgap': float(config['CplexAccuracy']),  # TODO: For the moment hardcoded, it has to be moved to a config file
                         'numericalemphasis': 0,
                         'mipdisplay': 4,
                         'scaind': 1,
                         'lpmethod': 0,
                         'relaxfixedinfeas': 0,
                         'mipstart': 1,
                         'mircuts': 1,
                         'quality': True,
                         'bardisplay': 2,
                         'epint': 0,
                         'lbheur': 1,
                         }
    elif config['CplexSetting'] == 'Default':
        cplex_options = {'epgap': float(config['CplexAccuracy']),  # TODO: For the moment hardcoded, it has to be moved to a config file
                         'numericalemphasis': 0,
                         'mipdisplay': 4,
                         'scaind': 1,
                         'lpmethod': 0,
                         'relaxfixedinfeas': 0,
                         'mipstart': 1,
                         'mircuts': 1,
                         'quality': True,
                         'bardisplay': 2,
                         'epint': 0,
                         'lbheur': 1,
                         }
        logging.info('Default Cplex setting used')
    elif config['CplexSetting'] == 'Agressive':
        cplex_options = {'epgap': config['CplexAccuracy'],  # TODO: For the moment hardcoded, it has to be moved to a config file
                         'numericalemphasis': 0,
                         'mipdisplay': 4,
                         'scaind': 1,
                         'lpmethod': 0,
                         'relaxfixedinfeas': 0,
                         'mipstart': 1,
                         'mircuts': 1,
                         'quality': True,
                         'bardisplay': 2,
                         'epint': 0,
                         'heuristiceffort': 2,
                         'lbheur': 1,
                         # Probing parameters
                         'probe': 1,
                         # Cut parameters
                         'cuts': 5,
                         'covers': 3,
                         'cliques': 3,
                         'disjcuts': 3,
                         'liftprojcuts': 3,
                         'localimplied': 3,
                         'flowcovers': 2,
                         'flowpaths': 2,
                         'fraccuts': 2,
                         }
        logging.info('Agressive Cplex setting used')
    else:
        logging.critical('Cplex setting must be specified Default or Agressive')
        sys.exit(1)

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
        if MTS:
            with open(os.path.join(sim, 'Inputs_MTS.p'), 'wb') as pfile:
                pickle.dump(SimData, pfile, protocol=pickle.HIGHEST_PROTOCOL)
        else:
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
