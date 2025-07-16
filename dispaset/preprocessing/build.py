import datetime as dt
import logging
import numbers
import os
import shutil
import sys

import numpy as np
import pandas as pd

from .data_check import check_units, check_sto, check_AvailabilityFactors, \
    check_clustering, isStorage, check_chp, check_df, check_MinMaxFlows, \
    check_FlexibleDemand, check_reserves, check_p2bs, check_boundary_sector, \
    check_BSFlexMaxCapacity, check_BSFlexMaxSupply, check_FFRLimit, check_PrimaryReserveLimit, check_CostXNotServed,\
    check_grid_data
from .data_handler import NodeBasedTable, load_time_series, UnitBasedTable, merge_series, define_parameter, \
    load_geo_data, GenericTable
from .reserves import percentage_reserve, probabilistic_reserve, generic_reserve
from .utils import select_units, interconnections, clustering, EfficiencyTimeSeries, \
    BoundarySectorEfficiencyTimeSeries, incidence_matrix, pd_timestep, PTDF_matrix, merge_lines
from .boundary_sector import zone_to_bs_mapping
from .. import __version__
from ..common import commons
from ..misc.gdx_handler import write_variables
from difflib import SequenceMatcher

GMS_FOLDER = os.path.join(os.path.dirname(__file__), '..', 'GAMS')


def build_single_run(config, profiles=None, PtLDemand=None, SectorXFlexDemand=None, SectorXFlexSupply=None, MTS=0,
                     profilesSectorX=None):
    """
    This function builds a simulation environment folder and populates it with
    the relevant input data.

    :param config:              ConfigDict instance or path to a yaml configuration file
    :param profiles:            Dictionary with the hydropower profiles resulting from the MTS
    :param PtLDemand:          Dictionary with the P2L profiles resulting from the MTS
    :param SectorXFlexDemand:   Dictionary with the SectorX flexible demand profiles resulting from the MTS
    :param SectorXFlexSupply:   Dictionary with the SectorX flexible supply profiles resulting from the MTS
    :param MTS:                 Boolean stating if we are running the MTS or not
    :param profilesSectorX:     Dictionary with the SectorX profiles resulting from the MTS
    :returns:                   Simulation data (dict)
    """
    dispa_version = __version__
    logging.info('Building simulation environment')

    # Checking the config first:
    if isinstance(config, str):
        config = load_config(config)

    # Boolean variable to check wether it is milp or lp:
    LP = config['SimulationType'] == 'LP' or config['SimulationType'] == 'LP clustered'

    # Boolean variable to check wether it is NTC or DC-POWERFLOW:
    grid_flag = config.get('TransmissionGridType', '')  # If key does not exist it returns ""

    # Remove SectorCoupling_flag declaration since it's always 'On' in the next version

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
    if (grid_flag == 'DC-Power Flow'):
        if os.path.isfile(config['GridData']):
            GridData = pd.read_csv(config['GridData'],
                                   na_values=commons['na_values'],
                                   keep_default_na=False,
                                   index_col=0)
            PTDF = PTDF_matrix(config, GridData).fillna(0)
        else:
            logging.error('No valid grid data file provided')
            GridData = pd.DataFrame()

    # Inertia Limit:
    if 'InertiaLimit' in config and os.path.isfile(config['InertiaLimit']):
        InertiaLimit = load_time_series(config, config['InertiaLimit']).fillna(0)
    else:
        logging.warning('No Inertia Limit will be considered (no valid file provided)')
        InertiaLimit = pd.DataFrame(index=config['idx_long'], data={'Value': 0})
        # logging.warning('No Inertia Limit will be considered (no valid file provided)')
        # InertiaLimit = pd.DataFrame(index=config['idx_long'])

    # Boundary Sector Interconnections:
    if 'BoundarySectorInterconnections' in config and os.path.isfile(config['BoundarySectorInterconnections']):
        BS_flows = load_time_series(config, config['BoundarySectorInterconnections']).fillna(0)
    else:
        logging.warning('No historical boundary sector flows will be considered (no valid file provided)')
        BS_flows = pd.DataFrame(index=config['idx_long'])
    if 'BoundarySectorNTC' in config and os.path.isfile(config['BoundarySectorNTC']):
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
    # aqui deberia inicializar los df de los resultados que no estan con timestamp

    # FFR Required:
    if 'FFRLimit' in config and os.path.isfile(config['FFRLimit']):
        FFRLimit = load_time_series(config, config['FFRLimit']).fillna(0)
    else:
        logging.warning('No FFR requirement will be considered (no valid file provided)')
        FFRLimit = pd.DataFrame(index=config['idx_long'], data={'Value': 0})

    # Primary Reserve Required:
    if 'PrimaryReserveLimit' in config and os.path.isfile(config['PrimaryReserveLimit']):
        PrimaryReserveLimit = load_time_series(config, config['PrimaryReserveLimit']).fillna(0)
    else:
        logging.warning('No Primary Reserve requirement will be considered (no valid file provided)')
        PrimaryReserveLimit = pd.DataFrame(index=config['idx_long'], data={'Value': 0})

    # Curtailment:
    CostCurtailment = NodeBasedTable('CostCurtailment', config, default=config['default']['CostCurtailment'])

    # Boundary Sectors:
    BoundarySector = pd.DataFrame()
    if 'BoundarySectorData' in config and os.path.isfile(config['BoundarySectorData']):
        BoundarySector = pd.read_csv(config['BoundarySectorData'],
                                     na_values=commons['na_values'],
                                     keep_default_na=False, index_col='Sector')
        for key in ['STOCapacity', 'STOSelfDischarge', 'STOMaxPower', 'STOMinSOC', 'STOHours']:
            if key not in BoundarySector.columns:
                BoundarySector[key] = np.nan
                logging.warning(f'{key} is not defined in the boundary sector data table')
    else:
        for key in ['STOCapacity', 'STOSelfDischarge', 'STOMaxPower', 'STOMinSOC', 'STOHours']:
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

    filter_col = [col for col in plants if col.startswith('Sector') and not col.startswith('SectorX')]
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
    # MARCO CHANGE: AGREGAR COLUMNAS QUE SON NECESARIAS PARA LAS NUEVAS FUNCIONES
    for key in ['CHPPowerLossFactor', 'CHPPowerToHeat', 'CHPType', 'STOCapacity', 'STOSelfDischarge',
                'STOMaxChargingPower', 'STOChargingEfficiency', 'CHPMaxHeat', 'WaterWithdrawal',
                'WaterConsumption', 'InertiaConstant', 'STOHours', 'Droop']:
        if key not in plants.columns:
            plants[key] = np.nan

    for key in ['Sector1', 'Sector2']:
        if key not in plants.columns:
            plants[key] = ''

    # check plant list:
    check_units(config, plants)
    plants.index = plants['Unit']

    # Defining the hydro storages:
    plants_sto = plants[[u in commons['tech_storage'] for u in plants['Technology']]]
    check_sto(config, plants_sto)
    check_sto(config, plants_sto.dropna(subset=['STOCapacity']))

    # Defining the boundary sector only units:
    plants_bs = plants[[u in commons['tech_boundary_sector'] for u in plants['Technology']]]
    check_boundary_sector(config, plants_bs, BoundarySector)

    # Defining the CHPs:
    plants_chp = plants[[str(x).lower() in commons['types_CHP'] for x in plants['CHPType']]]
    check_chp(config, plants_chp)
    
    # Defining the reserves:
    # TODO: add list of VRE to the Reserve participation in config and remove from here
    plants_res = plants[[u in config['ReserveParticipation'] or u in config['ReserveParticipation_CHP'] or u in commons['tech_renewables'] for u in plants['Technology']]]
 

    # Defining the P2BS units:
    plants_p2bs = plants[[u in commons['tech_p2bs'] for u in plants['Technology']]]
    check_p2bs(config, plants_p2bs)

    # Define all Boundary Sector units:
    #plants_all_bs = pd.concat([plants_p2bs, plants_bs])
    plants_all_bs = pd.concat([plants_bs, plants_chp])

    Outages = UnitBasedTable(plants, 'Outages', config, fallbacks=['Unit', 'Technology'])
    AF = UnitBasedTable(plants, 'RenewablesAF', config, fallbacks=['Unit', 'Technology'], default=1,
                        RestrictWarning=commons['tech_renewables'])
    AF = AF.apply(pd.to_numeric)

    ReservoirLevels = UnitBasedTable(plants_sto, 'ReservoirLevels', config,
                                     fallbacks=['Unit', 'Technology', 'Zone'],
                                     default=0)

    if 'StorageAlertLevels' in config and os.path.isfile(config['StorageAlertLevels']):
        StorageAlertLevels = UnitBasedTable(plants_sto, 'StorageAlertLevels', config,
                                            fallbacks=['Unit', 'Technology', 'Zone'],
                                            default=0)
    else:
        logging.warning('No Storage Alert Levels will be considered (no valid file provided)')
        StorageAlertLevels = pd.DataFrame(index=config['idx_long'])

    if 'StorageFloodControl' in config and os.path.isfile(config['StorageFloodControl']):
        StorageFloodControl = UnitBasedTable(plants_sto, 'StorageFloodControl', config,
                                             fallbacks=['Unit', 'Technology', 'Zone'],
                                             default=1)
    else:
        logging.warning('No Storage Flood Control will be considered (no valid file provided)')
        StorageFloodControl = pd.DataFrame(index=config['idx_long'])

    if 'ReservoirScaledInflows' in config and os.path.isfile(config['ReservoirScaledInflows']):
        ReservoirScaledInflows = UnitBasedTable(plants_sto, 'ReservoirScaledInflows', config,
                                                fallbacks=['Unit', 'Technology', 'Zone'], default=0)
    else:
        logging.warning('No historical Reservoir Scaled Inflows will be considered (no valid file provided)')
        ReservoirScaledInflows = pd.DataFrame(index=config['idx_long'])

    if 'ReservoirScaledOutflows' in config and os.path.isfile(config['ReservoirScaledOutflows']):
        ReservoirScaledOutflows = UnitBasedTable(plants_sto, 'ReservoirScaledOutflows', config,
                                                 fallbacks=['Unit', 'Technology', 'Zone'], default=0)
    else:
        logging.warning('No historical outflows will be considered (no valid file provided)')
        ReservoirScaledOutflows = pd.DataFrame(index=config['idx_long'])
    # TODO: Check if plants_sto esta bien para asignar el CostOfSpillage
    if 'CostOfSpillage' in config and os.path.isfile(config['CostOfSpillage']):
        CostOfSpillage = UnitBasedTable(plants_sto, 'CostOfSpillage', config,
                                        fallbacks=['Unit', 'Technology', 'Zone'],
                                        default=0)
    elif 'CostOfSpillage' in config:
        CostOfSpillage = UnitBasedTable(plants_sto, 'CostOfSpillage', config, default=config['default']['CostOfSpillage'])
    else:
        logging.warning('No CostOfSpillage will be considered (no valid file provided)')
        CostOfSpillage = pd.DataFrame(index=config['idx_long'])

    # Detecting thermal zones:
    if 'Zone_th' in plants_chp.columns:
        zones_th = plants_chp['Zone_th'].unique().tolist()
    else:
        zones_th = []
    if '' in zones_th:
        zones_th.remove('')

    # addding the thermal zones to the boundary sectors
    for z in zones_th:
        if z not in BoundarySector.index:
            BoundarySector.loc[z, 'STOCapacity'] = 0
            BoundarySector.loc[z, 'CostXNotServed'] = config['default']['CostHeatSlack']

    # Detecting boundary zones:
    if len(BoundarySector.index) >0:
        zones_bs = BoundarySector.index.to_list()
        if '' in zones_bs:
            zones_bs.remove('')
        if np.nan in zones_bs:
            zones_bs = [x for x in zones_bs if pd.isnull(x) == False]
        if 'nan' in zones_bs:
            zones_bs.remove('nan')

        zones_bss = np.unique(zones_bs + list(BoundarySector.index)).tolist()

        if 'SectorXDemand' in config and os.path.isfile(config['SectorXDemand']):
            SectorXDemand = GenericTable(zones_bs, 'SectorXDemand', config, default=0)
        else:
            logging.warning('No SectorXDemand will be considered (no valid file provided)')
            SectorXDemand = pd.DataFrame(index=config['idx_long'])

        if 'CostXNotServed' in config and os.path.isfile(config['CostXNotServed']):
            CostXNotServed = GenericTable(zones_bs, 'CostXNotServed', config)
            check_CostXNotServed(config, CostXNotServed, zones_bs)
        else:
            # Create a DataFrame with the CostXNotServed values from BoundarySector
            CostXNotServed = pd.DataFrame(index=config['idx_long'])
            for zone in zones_bs:
                if zone in BoundarySector.index and 'CostXNotServed' in BoundarySector.columns:
                    CostXNotServed[zone] = BoundarySector.loc[zone, 'CostXNotServed']
                else:
                    logging.critical(f'CostXNotServed is not defined for boundary sector {zone} in the boundary sector data table')
                    sys.exit(1)
            check_CostXNotServed(config, CostXNotServed, zones_bs)
        BoundarySector = BoundarySector[BoundarySector.index.isin(zones_bs)]
        BoundarySector.fillna(0, inplace=True)

        if 'SectorXReservoirLevels' in config and os.path.isfile(config['SectorXReservoirLevels']):
            SectorXReservoirLevels = GenericTable(zones_bs, 'SectorXReservoirLevels', config, default=0)
        else:
            logging.warning('No SectorXReservoirLevels will be considered (no valid file provided)')
            SectorXReservoirLevels = pd.DataFrame(index=config['idx_long'])

        if 'SectorXAlertLevel' in config and os.path.isfile(config['SectorXAlertLevel']):
            SectorXAlertLevel = GenericTable(zones_bs, 'SectorXAlertLevel', config, default=0)
        else:
            logging.warning('No SectorXAlertLevel will be considered (no valid file provided)')
            SectorXAlertLevel = pd.DataFrame(index=config['idx_long'])

        if 'SectorXFloodControl' in config and os.path.isfile(config['SectorXFloodControl']):
            SectorXFloodControl = GenericTable(zones_bs, 'SectorXFloodControl', config, default=0)
        else:
            logging.warning('No SectorXFloodControl will be considered (no valid file provided)')
            SectorXFloodControl = pd.DataFrame(index=config['idx_long'])

        # Boundary Sector Max Spillage
        if 'BoundarySectorMaxSpillage' in config and os.path.isfile(config['BoundarySectorMaxSpillage']):
            BS_spillage = load_time_series(config, config['BoundarySectorMaxSpillage']).fillna(0)
        else:
            BS_spillage = pd.DataFrame(index=config['idx_long'])
            logging.warning('No maximum spillage capacity provided.')
            # sys.exit('No maximum spillage capacity provided. This parameter is necessary.')

        # Boundary Sector Forced Spillage
        if 'BoundarySectorMaxSpillage' in config and os.path.isfile(config['BoundarySectorMaxSpillage']):
            BS_forced_spillage = pd.DataFrame(0, index=BS_spillage.index, columns=BS_spillage.columns)
        else:
            BS_forced_spillage = pd.DataFrame(index=config['idx_long'])

        if 'CostXSpillage' in config and config['CostXSpillage'] != '' and os.path.isfile(config['CostXSpillage']):
            CostXSpillage = load_time_series(config, config['CostXSpillage']).fillna(0)
        else:
            logging.warning('No CostXSpillage will be considered (no valid file provided)')
            CostXSpillage = pd.DataFrame(index=config['idx_long'])

        # Read BS Flexible demand & supply
        if 'SectorXFlexibleDemand' in config and os.path.isfile(config['SectorXFlexibleDemand']):
            SectorXFlexibleDemand = GenericTable(zones_bs, 'SectorXFlexibleDemand', config, default=0)
        else:
            logging.warning('No SectorXFlexibleDemand will be considered (no valid file provided)')
            SectorXFlexibleDemand = pd.DataFrame(index=config['idx_long'])

        if 'SectorXFlexibleSupply' in config and os.path.isfile(config['SectorXFlexibleSupply']):
            SectorXFlexibleSupply = GenericTable(zones_bs, 'SectorXFlexibleSupply', config, default=0)
        else:
            logging.warning('No SectorXFlexibleSupply will be considered (no valid file provided)')
            SectorXFlexibleSupply = pd.DataFrame(index=config['idx_long'])
    else:
        zones_bs = list()
        BSInterconnections = list()
        SectorXDemand = pd.DataFrame(index=config['idx_long'])
        CostXNotServed = pd.DataFrame(index=config['idx_long'])
        SectorXReservoirLevels = pd.DataFrame(index=config['idx_long'])
        SectorXAlertLevel = pd.DataFrame(index=config['idx_long'])
        SectorXFloodControl = pd.DataFrame(index=config['idx_long'])
        BS_spillage = pd.DataFrame(index=config['idx_long'])
        BS_forced_spillage = pd.DataFrame(index=config['idx_long'])
        CostXSpillage = pd.DataFrame(index=config['idx_long'])
        SectorXFlexibleDemand = pd.DataFrame(index=config['idx_long'])
        SectorXFlexibleSupply = pd.DataFrame(index=config['idx_long'])

    # Update reservoir levels with newly computed ones from the mid-term scheduling
    if profiles is not None:
        plants_sto.set_index(plants_sto.loc[:, 'Unit'], inplace=True, drop=True)
        for key in profiles.columns:
            if key not in ReservoirLevels.columns:
                logging.warning('The reservoir profile "' + key + '" provided by the MTS is not found in the '
                                                                  'ReservoirLevels table')
            elif key in list(ReservoirLevels.loc[:, plants_sto['Technology'] == 'SCSP'].columns):
                ReservoirLevels[key] = config['default']['ReservoirLevelInitial']
                logging.info('The reservoir profile "' + key + '" can not be seleceted for MTS, instead, default value '
                                                               'of: ' + str(
                                                                   config['default']['ReservoirLevelInitial']) + ' will be used')
            else:
                ReservoirLevels.loc[:, key] = profiles[key]
                logging.info(
                    'The reservoir profile "' + key + '" provided by the MTS is used as target reservoir level')

    if profilesSectorX is not None:
        for key in profilesSectorX.columns:
            if key not in SectorXReservoirLevels.columns:
                logging.warning('The reservoir profile "' + key + '" provided by the MTS is not found in the '
                                                                  'ReservoirLevels table')
            else:
                SectorXReservoirLevels.loc[:, key] = profilesSectorX[key]
                logging.info(
                    'The reservoir profile "' + key + '" provided by the MTS is used as target reservoir level')

    # Update SectorXFlexDemand (SectorXFlexibleDemand with demand from mid term scheduling)
    if SectorXFlexDemand is not None and any(SectorXFlexibleDemand) > 0:
        for key in SectorXFlexDemand.columns:
            if key not in SectorXFlexibleDemand.columns:
                logging.warning('The BS flexible demand "' + key + '" provided by the MTS is not found in the '
                                                                   'SectorXFlexibleDemand table')
            else:
                SectorXFlexibleDemand.loc[:, key] = SectorXFlexDemand[key]

    # Update SectorXFlexSupply (SectorXFlexibleSupply with demand from mid term scheduling)
    if SectorXFlexSupply is not None and any(SectorXFlexibleSupply) > 0:
        for key in SectorXFlexSupply.columns:
            if key not in SectorXFlexibleSupply.columns:
                logging.warning('The BS flexible demand "' + key + '" provided by the MTS is not found in the '
                                                                   'SectorXFlexibleSupply table')
            else:
                SectorXFlexibleSupply.loc[:, key] = SectorXFlexSupply[key]

    # data checks:
    check_AvailabilityFactors(plants, AF)
    check_FlexibleDemand(ShareOfFlexibleDemand)
    check_reserves(Reserve2D, Reserve2U, Load)
    check_FFRLimit(FFRLimit, Load)
    check_PrimaryReserveLimit(PrimaryReserveLimit, Load)

    # Fuel prices:
    fuels = ['PriceOfNuclear', 'PriceOfBlackCoal', 'PriceOfGas', 'PriceOfFuelOil', 'PriceOfBiomass', 'PriceOfCO2',
             'PriceOfLignite', 'PriceOfPeat', 'PriceOfAmmonia']
    FuelPrices = {}
    for fuel in fuels:
        FuelPrices[fuel] = NodeBasedTable(fuel, config, default=config['default'][fuel])

    # Interconnections:
    NTC.index = NTC.index.astype('datetime64[ns]')
    [Interconnections_sim, Interconnections_RoW, _ ] = interconnections(config['zones'], NTC, flows)
    
    # Create a dictionary to hold all price transmission data
    price_transmission_data = {}
    for l in (Interconnections_sim.columns.tolist() + flows.columns.tolist()):
        if l in PriceTransmission_raw:
            price_transmission_data[l] = PriceTransmission_raw[l].reindex(Interconnections_sim.index)
        else:
            price_transmission_data[l] = pd.Series(config['default']['PriceTransmission'], index=Interconnections_sim.index)
            if config['default']['PriceTransmission'] > 0:
                logging.warning('No detailed values were found the transmission prices of line ' + l +
                                '. Using default value ' + str(config['default']['PriceTransmission']))
    
    # Create DataFrame all at once to avoid fragmentation
    PriceTransmission = pd.DataFrame(price_transmission_data, index=Interconnections_sim.index)
    
    # Merge the symetrical lines between two zones
    [Interconnections_sim, Interconnections_RoW, PriceTransmission] = merge_lines(Interconnections_sim, Interconnections_RoW, PriceTransmission_raw)

    # Bidirectional Net transfer capacities:
    if (grid_flag == 'DC-Power Flow'):
        NTCs = pd.DataFrame(index=config['idx_long'])
        for l in GridData['Code_Line']:
            NTCs[l] = GridData.loc[l, 'Fmax']
        # Check the Gridata Matrix consistency
        check_grid_data(NTCs.columns, PTDF, config)
        # Formatting PTDF matrix according to sets['l'] and sets['n'] order or sequence
        # Reindex to ensure rows and columns are in the right order
        PTDF = PTDF.reindex(index=NTCs.columns, columns=config['zones'])
    elif len(Interconnections_sim.columns) > 0:
        NTCs = Interconnections_sim.reindex(config['idx_long'])
        # In case of NTC-based simulation, the PTDF matrix remains blank:
        PTDF = pd.DataFrame(index=NTCs.columns, columns=config['zones'])  
    else:
        NTCs = pd.DataFrame(index=config['idx_long'])
        # if there are not lines, the PTDF matrix is empty:
        PTDF = pd.DataFrame(columns=config['zones'])
    # Historical flows with the rest of the world:
    Inter_RoW = Interconnections_RoW.reindex(config['idx_long'])
    all_lines = list(NTCs.columns) + list(Inter_RoW.columns)

    # Boundary Sector Interconnections:
    [BSInterconnections_sim, BSInterconnections_RoW, BSInterconnections] = interconnections(zones_bs, BS_NTC, BS_flows)
    # Boundary Sector Spillage:
    [BSSpillage_sim, BSSpillage_RoW, BSSpillage] = interconnections(zones_bs, BS_spillage, BS_forced_spillage)

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
                                   'STOMaxPower': 'SectorXStoragePowerMax',
                                   'STOMaxChargingPower': 'BoundarySectorStorageChargingCapacity',
                                   'STOMinSOC': 'SectorXStorageMinimum', 'STOHours': 'SectorXStorageHours'}, inplace=True)

    Plants_merged.rename(columns={'StartUpCost': 'CostStartUp',
                                  'RampUpMax': 'RampUpMaximum',
                                  'RampDownMax': 'RampDownMaximum',
                                  'MinUpTime': 'TimeUpMinimum',
                                  'MinDownTime': 'TimeDownMinimum',
                                  'RampingCost': 'CostRampUp',
                                  'STOCapacity': 'StorageCapacity',
                                  'STOHours': 'StorageHours',
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
                Plants_merged.loc[u, 'StorageHours'] = Plants_merged.loc[u, 'StorageHours'] * config['modifiers'][
                    'Storage']
                Plants_merged.loc[u, 'StorageChargingCapacity'] = Plants_merged.loc[u, 'StorageChargingCapacity'] * \
                config['modifiers']['Storage']

    filter_col = [col for col in Plants_merged if col.startswith('Sector') and not col.startswith('SectorX')]
    Plants_merged[filter_col] = Plants_merged[filter_col].astype(str)

    # Defining the hydro storages:
    Plants_sto = Plants_merged[[u in commons['tech_storage'] for u in Plants_merged['Technology']]]
    # Defining all storage units:
    plants_sto = Plants_merged[[u in [x for x in commons['Technologies'] if x in commons['tech_storage']] for u in
                                Plants_merged['Technology']]]
    # check storage plants:
    check_sto(config, Plants_sto, raw_data=False)
    check_sto(config, plants_sto, raw_data=False)
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

    # Filter conventional units
    Plants_conventional = Plants_merged[[u in commons['tech_conventional'] for u in Plants_merged['Technology']]].copy()

    # Filter batteries units
    Plants_batteries = Plants_merged[[u in commons['tech_batteries'] for u in Plants_merged['Technology']]].copy()

    # Defining the reserves:
    # TODO: add list of VRE to the Reserve participation in config and remove from here
    Plants_res = Plants_merged[[u in config['ReserveParticipation'] or u in config['ReserveParticipation_CHP']  or u in commons['tech_renewables'] for u in Plants_merged['Technology']]]

    # Filter boundary sector only plants
    Plants_boundary_sector_only = Plants_merged[
        [u in commons['tech_boundary_sector'] for u in Plants_merged['Technology']]].copy()
    # Check boundary sector units
    check_boundary_sector(config, Plants_boundary_sector_only, BoundarySector)

    # Filter power to boundary sector plants
    Plants_p2bs = Plants_merged[[u in commons['tech_p2bs'] for u in Plants_merged['Technology']]].copy()
    # Check heat only units
    check_p2bs(config, Plants_p2bs)

    # Filter boundary sector to power plants
    Plants_x2p = Plants_merged[[u in commons['tech_bs2p'] for u in Plants_merged['Technology']]].copy()

    # Water storage
    Plants_wat = Plants_merged[(Plants_merged['Fuel'] == 'WAT') & (Plants_merged['Technology'] != 'HROR')].copy()

    # Calculating the efficiency time series for each unit:
    Efficiencies = EfficiencyTimeSeries(config, Plants_merged)

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
    
                
    # New Reserves calculation
    # Calculate the percentage that each zone represents of the total load in each hour
    Load_shares = Load.div(Load.sum(axis=1), axis=0)    
    # Prorratear el valor total de FFR por los porcentajes de cada zona
    FFRLimit_allocated = Load_shares.mul(FFRLimit.iloc[:, 0], axis=0)
    PrimaryReserveLimit_allocated = Load_shares.mul(PrimaryReserveLimit.iloc[:, 0], axis=0)

    # %% Store all times series and format

    # Formatting all time series (merging, resempling) and store in the FinalTS dict
    finalTS = {'Load': Load, 'Reserve2D': reserve_2D_tot, 'Reserve2U': reserve_2U_tot, 'Reserve3U': reserve_2U_tot,
               'Efficiencies': Efficiencies,
               'NTCs': NTCs, 'Inter_RoW': Inter_RoW,'PriceTransmission': PriceTransmission,
               'BS_NTCs': BS_NTCs, 'BS_Inter_RoW': BS_Inter_RoW,
               'EfficienciesBoundarySector': BoundarySectorEfficiencies['Efficiency'],
               'LoadShedding': LoadShedding, 'CostLoadShedding': CostLoadShedding, 'CostCurtailment': CostCurtailment,
               'ScaledInflows': ReservoirScaledInflows, 'ScaledOutflows': ReservoirScaledOutflows,
               'ReservoirLevels': ReservoirLevels, 'Outages': Outages, 'AvailabilityFactors': AF,
               'StorageAlertLevels': StorageAlertLevels, 'StorageFloodControl': StorageFloodControl,
               'ShareOfFlexibleDemand': ShareOfFlexibleDemand,
               'SectorXDemand': SectorXDemand, 'CostXNotServed': CostXNotServed,
               'SectorXFlexibleDemand': SectorXFlexibleDemand, 'SectorXFlexibleSupply': SectorXFlexibleSupply,
               'BSMaxSpillage': BS_Spillages, 'SectorXReservoirLevels': SectorXReservoirLevels, 'SectorXAlertLevel': SectorXAlertLevel,
               'SectorXFloodControl': SectorXFloodControl,
               'CostOfSpillage': CostOfSpillage, 'CostXSpillage': CostXSpillage,
               'InertiaLimit': InertiaLimit, 'FFRU': FFRLimit_allocated, 'FFRD': FFRLimit_allocated, 
               'PFRU': PrimaryReserveLimit_allocated, 'PFRD': PrimaryReserveLimit_allocated}

    # Merge the following time series with weighted averages
    for key in ['ScaledInflows', 'ScaledOutflows', 'Outages', 'AvailabilityFactors']:
        finalTS[key] = merge_series(Plants_merged, plants, finalTS[key], tablename=key)

    # Merge the following time series with weighted averages
    for key in ['CostOfSpillage']:
        finalTS[key] = merge_series(Plants_merged, plants, finalTS[key], tablename=key)

    # Merge the following time series by weighted average based on storage capacity
    for key in ['ReservoirLevels', 'StorageAlertLevels', 'StorageFloodControl']:
        finalTS[key] = merge_series(Plants_merged, plants, finalTS[key], tablename=key, method='StorageWeightedAverage')

    # Check that all times series data is available with the specified data time step:
    for key in FuelPrices:
        check_df(FuelPrices[key], StartDate=config['idx'][0], StopDate=config['idx'][-1], name=key)
    for key in finalTS:
        if key in ['EfficienciesBoundarySector', 'ChargingEfficienciesBoundarySector']:
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
            if key in ['EfficienciesBoundarySector', 'ChargingEfficienciesBoundarySector']:
                for ts in finalTS[key]:
                    if len(finalTS[key][ts].columns) > 0:
                        finalTS[key][ts] = finalTS[key][ts].resample(pd_timestep(config['SimulationTimeStep'])).mean()
            else:
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
    sets['nx'] = zones_bs
    sets['nx_CC'] = pd.Series(zones_bs)[pd.Series(zones_bs) != 'S_EC'].tolist()
    sets['au'] = Plants_merged.index.tolist()
    sets['l'] = all_lines 
    sets['l_int'] = list(NTCs.columns)          #lines between zones
    sets['l_RoW'] = list(Inter_RoW.columns)  #lines between zones and RoW
    sets['lx'] = BSInterconnections

    # TODO: CHECK THIS SECTION OF CODE IT WAS AN ADAPTATION
    if 'BoundarySectorMaxSpillage' in config and os.path.isfile(config['BoundarySectorMaxSpillage']):
        sets['slx'] = BSSpillage  # ALIZON: CON BS
    else:
        sets['slx'] = ReservoirScaledInflows.columns.tolist()  # MARCO: SIN BS

    sets['f'] = commons['Fuels']
    sets['p'] = ['CO2']
    sets['s'] = Plants_sto.index.tolist()
    sets['sx'] = Plants_wat[Plants_wat['Technology'] == 'HDAMC'].index.tolist()
    sets['u'] = Plants_merged[[u in [x for x in commons['Technologies'] if x not in commons['tech_p2bs'] + commons['tech_boundary_sector']]
                               for u in Plants_merged['Technology']]].index.tolist()
    sets['chp'] = Plants_chp.index.tolist()
    sets['p2x'] = Plants_p2bs.index.tolist()
    sets['x2p'] = Plants_x2p.index.tolist()
    sets['t'] = commons['Technologies']
    sets['tr'] = commons['tech_renewables']
    sets['tc'] = list(set(commons['Technologies']) - set(commons['tech_renewables']))
    sets['wat'] = Plants_wat.index.tolist()
    sets['xu'] = Plants_boundary_sector_only.index.tolist()
    sets['asu'] = Plants_merged[[u in [x for x in commons['Technologies'] if x in commons['tech_storage']] for u in
                                 Plants_merged['Technology']]].index.tolist()
    sets['cu'] = Plants_conventional.index.tolist()
    sets['ba'] = Plants_batteries.index.tolist()
    sets['res_U'] = ['FFRU', 'PFRU', '2U', '3U']
    sets['res_D'] = ['FFRD', 'PFRD', '2D']
    sets['res'] = sets['res_U'] + sets['res_D']
    
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
    sets_param['CostXStorageAlert'] = ['nx', 'h']
    sets_param['CostXFloodControl'] = ['nx', 'h']
    sets_param['CostXSpillage'] = ['slx', 'h']
    sets_param['CostOfSpillage'] = ['asu', 'h']
    sets_param['CostStorageAlert'] = ['au', 'h']
    sets_param['CostFloodControl'] = ['au', 'h']
    sets_param['Curtailment'] = ['n']
    sets_param['CostCurtailment'] = ['n', 'h']
    sets_param['Demand'] = ['mk', 'n', 'h']
    sets_param['ReserveDemand'] = ['res', 'n', 'h']      #new for reserves FCR
    sets_param['Efficiency'] = ['au', 'h']
    sets_param['X2PowerConversionMultiplier'] = ['nx', 'x2p', 'h']
    sets_param['Power2XConversionMultiplier'] = ['nx', 'p2x', 'h']
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
    sets_param['Nunits'] = ['au']
    sets_param['OutageFactor'] = ['au', 'h']
    sets_param['PartLoadMin'] = ['au']
    sets_param['PowerCapacity'] = ['au']
    sets_param['PowerCapacityBoundarySector'] = ['au', 'nx']
    sets_param['PowerInitial'] = ['au']
    sets_param['PriceTransmission'] = ['l', 'h']
    sets_param['RampUpMaximum'] = ['au']
    sets_param['RampDownMaximum'] = ['au']
    sets_param['RampStartUpMaximum'] = ['au']
    sets_param['RampShutDownMaximum'] = ['au']
    sets_param['ReserveParticipation'] = ['res', 'au', 'h']  # changed this also in the gams file(in the definition and in the equations satifying the reserve demand)
    sets_param['StorageCapacity'] = ['au']
    sets_param['StorageHours'] = ['au']
    sets_param['StorageChargingCapacity'] = ['au']
    sets_param['StorageChargingEfficiency'] = ['au']
    sets_param['StorageDischargeEfficiency'] = ['asu']
    sets_param['StorageSelfDischarge'] = ['au']
    sets_param['StorageInflow'] = ['asu', 'h']
    sets_param['StorageInitial'] = ['asu']
    sets_param['StorageMinimum'] = ['asu']
    sets_param['StorageAlertLevel'] = ['asu', 'h']
    sets_param['StorageFloodControl'] = ['asu', 'h']
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
    sets_param['SectorXStoragePowerMax'] = ['nx']
    sets_param['SectorXStorageHours'] = ['nx']
    sets_param['SectorXStorageSelfDischarge'] = ['nx']
    sets_param['SectorXStorageMinimum'] = ['nx']
    sets_param['SectorXStorageInitial'] = ['nx']
    sets_param['SectorXStorageProfile'] = ['nx', 'h']
    sets_param['SectorXAlertLevel'] = ['nx', 'h']
    sets_param['SectorXFloodControl'] = ['nx', 'h']
    # if MTS == 0:
    sets_param['InertiaConstant'] = ['au']
    sets_param['InertiaLimit'] = ['h']
    sets_param['Droop'] = ['au']
    sets_param['PrimaryReserveLimit'] = ['h']
    sets_param['FFRLimit'] = ['h']
    sets_param['PTDF'] = ['l_int', 'n']

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
    for var in ['Technology', 'Fuel', 'ReserveParticipation', 'Location', 'LocationX']:
        parameters[var] = define_parameter(sets_param[var], sets, value='bool')
    # %%
        # List of parameters whose value is known and provided in the dataframe BoundarySector
        for var in ['SectorXStorageCapacity', 'SectorXStorageSelfDischarge', 'SectorXStorageHours', 'SectorXStoragePowerMax']:
            parameters[var]['val'] = BoundarySector[var].values

    # List of parameters whose value is known, and provided in the dataframe Plants_merged.
    if MTS == 0:
        for var in ['PowerCapacity', 'PartLoadMin', 'TimeUpMinimum', 'TimeDownMinimum', 'CostStartUp',
                    'CostRampUp', 'StorageCapacity', 'StorageHours', 'StorageSelfDischarge', 'StorageChargingCapacity',
                    'InertiaConstant', 'Droop']:
            parameters[var]['val'] = Plants_merged[var].values
    else:
        for var in ['PowerCapacity', 'PartLoadMin', 'TimeUpMinimum', 'TimeDownMinimum', 'CostStartUp',
                    'CostRampUp', 'StorageCapacity', 'StorageHours', 'StorageSelfDischarge', 'StorageChargingCapacity']:
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
    parameters['StorageDischargeEfficiency']['val'] = Plants_sto['Efficiency'].values

    # List of parameters whose value is known, and provided in the dataframe Plants_chp
    for var in ['CHPPowerToHeat', 'CHPPowerLossFactor', 'CHPMaxHeat']:
        parameters[var]['val'] = Plants_chp[var].values
    if len(BoundarySector.index) >0:
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

            # Setting the storage alert levels
        if len(finalTS['StorageAlertLevels'].columns) != 0:
            if s in plants_sto.index:
                parameters['StorageAlertLevels']['val'][i, :] = finalTS['StorageAlertLevels'][s][idx_sim].values
            else:
                logging.warning('Storage Alert Levels not found for unit ' + s + '. Assuming no Storage Alert Levels')

        if len(finalTS['StorageFloodControl'].columns) != 0:
            if s in plants_sto.index:
                parameters['StorageFloodControl']['val'][i, :] = finalTS['StorageFloodControl'][s][idx_sim].values
            else:
                logging.warning('Storage Flood Control not found for unit ' + s + '. Assuming no Storage Flood Control')

        if len(finalTS['CostOfSpillage'].columns) != 0:
            if s in plants_sto.index:
                parameters['CostOfSpillage']['val'][i, :] = finalTS['CostOfSpillage'][s][idx_sim].values
            else:
                logging.warning('Cost Of Spillage not found for unit ' + s + '. Assuming no Cost Of Spillage')

        # if s in plants_sto.index:
        #     parameters['StorageAlertLevels']['val'][i, :] = finalTS['StorageAlertLevels'][s][idx_sim].values
        #     parameters['StorageFloodControl']['val'][i, :] = finalTS['StorageFloodControl'][s][idx_sim].values
        #     parameters['CostOfSpillage']['val'][i, :] = finalTS['CostOfSpillage'][s][idx_sim].values

        # The initial level is the same as the first value of the profile:
        if s in Plants_sto.index:
            parameters['StorageInitial']['val'][i] = parameters['StorageProfile']['val'][i, 0] * \
                finalTS['AvailabilityFactors'][s][idx_sim[0]] * \
                Plants_sto.loc[s, 'StorageCapacity'] * Plants_sto.loc[s, 'Nunits']

    if len(BoundarySector.index) >0:
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
                parameters['SectorXStorageProfile']['val'][i, :] = np.where(
                    BoundarySector.loc[nx, 'SectorXStorageCapacity'] == 0, 0,
                    np.linspace(config['default']['ReservoirLevelInitial'], config['default']['ReservoirLevelFinal'],
                                len(idx_sim)))

            # Setting Storage Alert
            if nx in finalTS['SectorXAlertLevel'] and any(finalTS['SectorXAlertLevel'][nx] > 0) and all(
                    finalTS['SectorXAlertLevel'][nx] - 1 <= 1e-11):
                parameters['SectorXAlertLevel']['val'][i, :] = finalTS['SectorXAlertLevel'][nx][idx_sim].values
            if nx in finalTS['SectorXFloodControl'] and any(finalTS['SectorXFloodControl'][nx] > 0) and all(
                    finalTS['SectorXFloodControl'][nx] - 1 <= 1e-11):
                parameters['SectorXFloodControl']['val'][i, :] = finalTS['SectorXFloodControl'][nx][idx_sim].values
            parameters['SectorXStorageInitial']['val'][i] = parameters['SectorXStorageProfile']['val'][i, 0] * \
            BoundarySector.loc[nx, 'SectorXStorageCapacity']
    # Storage Inflows:
    for i, s in enumerate(sets['asu']):
        if s in finalTS['ScaledInflows']:
            parameters['StorageInflow']['val'][i, :] = finalTS['ScaledInflows'][s][idx_sim].values * \
                plants_sto.loc[s, 'PowerCapacity']
            if s in finalTS['ScaledOutflows']:
                parameters['StorageOutflow']['val'][i, :] = finalTS['ScaledOutflows'][s][idx_sim].values * \
                    plants_sto.loc[s, 'PowerCapacity']
    if len(BoundarySector.index) >0:
    # Boundary sector demands:
        for i, u in enumerate(sets['nx']):
            if u in finalTS['SectorXDemand']:
                parameters['SectorXDemand']['val'][i, :] = finalTS['SectorXDemand'][u][idx_sim].values
            # MARCO CHANGE
            if u in finalTS['CostXNotServed']:
                parameters['CostXNotServed']['val'][i, :] = finalTS['CostXNotServed'][u][idx_sim].values

        # Boundary Sector time series
        for i, u in enumerate(sets['nx']):
            if u in finalTS['SectorXFlexibleDemand']:
                parameters['SectorXFlexDemandInput']['val'][i, :] = finalTS['SectorXFlexibleDemand'][u][idx_sim].values
            if u in finalTS['SectorXFlexibleSupply']:
                parameters['SectorXFlexSupplyInput']['val'][i, :] = finalTS['SectorXFlexibleSupply'][u][idx_sim].values

        if 'SectorXFlexibleDemand' in config:
            check_BSFlexMaxCapacity(parameters, config, sets)

        if 'SectorXFlexibleSupply' in config:
            check_BSFlexMaxSupply(parameters, config, sets)

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

    if len(BoundarySector.index) >0:
        # Initialize conversion multipliers with zeros
        values_p2x = np.zeros([len(sets['nx']), len(sets['p2x']), len(sets['h'])])
        values_x2p = np.zeros([len(sets['nx']), len(sets['x2p']), len(sets['h'])])

        # Handle p2x units
        for i in range(len(sets['nx'])):
            for u in range(len(sets['p2x'])):
                if Plants_merged.loc[sets['p2x'][u], 'Technology'] in commons['tech_p2bs']:
                    values_p2x[i, u, :] = finalTS['EfficienciesBoundarySector'][sets['nx'][i]][sets['p2x'][u]]
                    logging.info(f"Unit {sets['p2x'][u]} is a p2x unit, setting Power2XConversionMultiplier")

        # Handle x2p units
        for i in range(len(sets['nx'])):
            for u in range(len(sets['x2p'])):
                if Plants_merged.loc[sets['x2p'][u], 'Technology'] in commons['tech_boundary_sector']:
                    values_x2p[i, u, :] = finalTS['EfficienciesBoundarySector'][sets['nx'][i]][sets['x2p'][u]]
                    logging.info(f"Unit {sets['x2p'][u]} is an x2p unit, setting X2PowerConversionMultiplier")

        # Set the parameters with the appropriate values
        parameters['Power2XConversionMultiplier'] = {'sets': sets_param['Power2XConversionMultiplier'], 'val': values_p2x}
        parameters['X2PowerConversionMultiplier'] = {'sets': sets_param['X2PowerConversionMultiplier'], 'val': values_x2p}

    values = np.ndarray([len(sets['mk']), len(sets['n']), len(sets['h'])])
    for i in range(len(sets['n'])):
        values[0, i, :] = finalTS['Load'][sets['n'][i]] * (1 - finalTS['ShareOfFlexibleDemand'][sets['n'][i]])
        values[1, i, :] = finalTS['Reserve2U'][sets['n'][i]]
        values[2, i, :] = finalTS['Reserve2D'][sets['n'][i]]
        values[3, i, :] = finalTS['Load'][sets['n'][i]] * finalTS['ShareOfFlexibleDemand'][sets['n'][i]]
    parameters['Demand'] = {'sets': sets_param['Demand'], 'val': values}

    # Reserves Demand
    values = np.ndarray([len(sets['res']), len(sets['n']), len(sets['h'])])
    for i in range(len(sets['n'])):
        values[0, i, :] = finalTS['FFRU'][sets['n'][i]]
        values[1, i, :] = finalTS['PFRU'][sets['n'][i]]
        values[2, i, :] = finalTS['Reserve2U'][sets['n'][i]]
        values[3, i, :] = finalTS['Reserve3U'][sets['n'][i]]
        values[4, i, :] = finalTS['FFRD'][sets['n'][i]]
        values[5, i, :] = finalTS['PFRD'][sets['n'][i]]
        values[6, i, :] = finalTS['Reserve2D'][sets['n'][i]]
    parameters['ReserveDemand'] = {'sets': sets_param['ReserveDemand'], 'val': values}

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
    CostVariable = pd.DataFrame()
    for unit in range(Nunits):
        c = Plants_merged['Zone'][unit]  # zone to which the unit belongs
        found = False
        for FuelEntry in FuelEntries:
            CostVariable = pd.concat([
            CostVariable, FuelPrices[FuelEntries[FuelEntry]][c] / Plants_merged['Efficiency'][unit] + \
                  Plants_merged['EmissionRate'][unit] * FuelPrices['PriceOfCO2'][c]], axis=1)
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

    # Calculate MaxCostVariable for each zone
    MaxCostVariable = pd.DataFrame(index=sets['n'])
    for zone in sets['n']:
        # Get units in this zone
        zone_units_df = Plants_merged[Plants_merged['Zone'] == zone]
        if not zone_units_df.empty:
            # Get the integer positions of these units in the original dataframe
            zone_unit_positions = [Plants_merged.index.get_loc(idx) for idx in zone_units_df.index]
            max_cost = np.max(parameters['CostVariable']['val'][zone_unit_positions])
            MaxCostVariable.loc[zone, 'MaxCost'] = max_cost if not np.isnan(max_cost) else 0
        else:
            MaxCostVariable.loc[zone, 'MaxCost'] = 0

    # Storage alert level:
    for unit in range(Nunits):
        c = Plants_merged['Zone'].iloc[unit]  # zone to which the unit belongs
        if not np.isnan(MaxCostVariable.loc[c, 'MaxCost']) and MaxCostVariable.loc[c, 'MaxCost'] > 0:
            if Plants_merged['Unit'].iloc[unit] in plants_sto['Unit']:
                parameters['CostStorageAlert']['val'][unit, :] = MaxCostVariable.loc[c, 'MaxCost']
        else:
            logging.warning('No alert price has been found for ' + Plants_merged['Unit'].iloc[unit] +
                            '. A null alert price has been assigned')
            if Plants_merged['Unit'].iloc[unit] in plants_sto['Unit']:
                parameters['CostStorageAlert']['val'][unit, :] = 0

    # Storage flood control:
    for unit in range(Nunits):
        c = Plants_merged['Zone'].iloc[unit]  # zone to which the unit belongs
        if not np.isnan(MaxCostVariable.loc[c, 'MaxCost']) and MaxCostVariable.loc[c, 'MaxCost'] > 0:
            if Plants_merged['Unit'].iloc[unit] in plants_sto['Unit']:
                parameters['CostFloodControl']['val'][unit, :] = MaxCostVariable.loc[c, 'MaxCost']
        else:
            logging.warning('No flood price has been found for ' + Plants_merged['Unit'].iloc[unit] +
                            '. A null flood price has been assigned')
            if Plants_merged['Unit'].iloc[unit] in plants_sto['Unit']:
                parameters['CostFloodControl']['val'][unit, :] = 0

    if len(BoundarySector.index) >0:
        BoundarySector['Sector'] = BoundarySector.index
        BoundarySector['Zone'] = BoundarySector.index.map(zone_to_bs_mapping(plants_all_bs))

    zones = list(MaxCostVariable.index)

    if len(BoundarySector)>0:
        for unit in range(len(BoundarySector)):
            found = False
            if BoundarySector['Zone'][unit] in zones:
                parameters['CostXStorageAlert']['val'][unit, :] = MaxCostVariable.loc[BoundarySector['Zone'][unit], 'MaxCost']
                found = True
            if not found:
                parameters['CostXStorageAlert']['val'][unit, :] = 0

        # Assign storage flood control level costs to the unit with highest variable costs inside the zone
        for unit in range(len(BoundarySector)):
            parameters['CostXFloodControl']['val'][unit, :] = 0

    # %%###############################################################################################################
    # Inertia limit to parameters dict
    # parameters['InertiaLimit']['val'] = InertiaLimit['InertiaLimit'].values
    # if MTS == 0:
    # if len(finalTS['InertiaLimit'].columns) != 0:
    #     parameters['InertiaLimit']['val'] = finalTS['InertiaLimit'].iloc[:, 0].values
    # # FFR to parameters dict
    # if len(finalTS['FFRLimit'].columns) != 0:
    #     parameters['FFRLimit']['val'] = finalTS['FFRLimit'].iloc[:, 0].values
    # # Primary Reserve to parameters dict
    # if len(finalTS['PrimaryReserveLimit'].columns) != 0:
    #     parameters['PrimaryReserveLimit']['val'] = finalTS['PrimaryReserveLimit'].iloc[:, 0].values

    # Maximum Line Capacity
    for i, l in enumerate(sets['l_int']):
        parameters['FlowMaximum']['val'][i, :] = finalTS['NTCs'][l]
        parameters['FlowMinimum']['val'][i, :] = -(finalTS['NTCs'][l])
        if l in finalTS['PriceTransmission'].columns:
            parameters['PriceTransmission']['val'][i, :] = finalTS['PriceTransmission'][l]
    N = len(sets['l_int'])
    for i,l in enumerate(sets['l_RoW']):
        parameters['FlowMaximum']['val'][N+i, :] = finalTS['Inter_RoW'][l]  
        parameters['FlowMinimum']['val'][N+i, :] = finalTS['Inter_RoW'][l]
        if l in finalTS['PriceTransmission'].columns:
            parameters['PriceTransmission']['val'][N+i, :] = finalTS['PriceTransmission'][l]

    # Check values:
    check_MinMaxFlows(parameters['FlowMinimum']['val'], parameters['FlowMaximum']['val'])

    # Incidence matrix:
    parameters['LineNode'] = incidence_matrix(sets, 'l', parameters, 'LineNode')

    # the PTDF matrix has dimensions l_int x n, where l_int is the number of lines between nodes and with RoW.
    parameters['PTDF']['val'] = PTDF.values

    # Cost Spillage without BS
    if 'CostOfSpillage' in config and os.path.isfile(config['CostOfSpillage']):
        for i, wat in enumerate(sets['wat']):
            if wat in finalTS['ScaledInflows'].columns:  # MARCO: SIN BS
                parameters['CostOfSpillage']['val'][i, :] = finalTS['CostOfSpillage'][wat]

    # Maximum Boundary Sector Line Capacity
    for i, lx in enumerate(sets['lx']):
        if lx in BS_NTCs.columns:
            parameters['FlowXMaximum']['val'][i, :] = finalTS['BS_NTCs'][lx]
        if lx in BS_Inter_RoW.columns:
            parameters['FlowXMaximum']['val'][i, :] = finalTS['BS_Inter_RoW'][lx]
            parameters['FlowXMinimum']['val'][i, :] = finalTS['BS_Inter_RoW'][lx]
    for i, slx in enumerate(sets['slx']):
        if slx in BS_Spillages.columns:
            parameters['SectorXMaximumSpillage']['val'][i, :] = finalTS['BSMaxSpillage'][slx]
    # Cost of Spillage boundary sector
        if 'CostXSpillage' in config and os.path.isfile(config['CostXSpillage']):
            parameters['CostXSpillage']['val'][i, :] = finalTS['CostXSpillage'][slx]  # ALIZON: CON BS

    # Check values:
    check_MinMaxFlows(parameters['FlowXMinimum']['val'], parameters['FlowXMaximum']['val'])

    parameters['LineXNode'] = incidence_matrix(sets, 'lx', parameters, 'LineXNode', nodes='nx')

    if 'BoundarySectorMaxSpillage' in config and os.path.isfile(config['BoundarySectorMaxSpillage']):
        parameters['SectorXSpillageNode'] = incidence_matrix(sets, 'slx', parameters, 'SectorXSpillageNode', nodes='nx')  # ALIZON: CON BS

    # Outage Factors
    if len(finalTS['Outages'].columns) != 0:
        for i, u in enumerate(sets['au']):
            if u in finalTS['Outages'].columns:
                parameters['OutageFactor']['val'][i, :] = finalTS['Outages'][u].values
            else:
                logging.warning('Outages factors not found for unit ' + u + '. Assuming no outages')

    # Participation to the reserve market
    # list_of_participating_units = []  # new list
    # for unit in Plants_merged.index:
    #     tech = Plants_merged.loc[unit, 'Technology']
    #     if tech in config['ReserveParticipation'] and Plants_merged.loc[unit, 'CHPType'] == '':
    #         list_of_participating_units.append(
    #             unit)  # if unit same technology as allowed without CHP and unit is no CHP then add to list
    #     elif tech in config['ReserveParticipation_CHP'] and Plants_merged.loc[unit, 'CHPType'] != '':
    #         list_of_participating_units.append(
    #             unit)  # if unit same technology as allowed with CHP and unit is CHP then add to list

    # values = np.array([s in list_of_participating_units for s in sets['au']],
    #                   dtype='bool')  # same as before but with new list
    # parameters['Reserve'] = {'sets': sets_param['Reserve'], 'val': values}
    
    # Binary table of participation to the reserve market
    # TODO: suggestion is to create commons by type of reserve instead of technologies
    constants = {
        'SystemFrequency': 50,
        'DeltaFrequencyMax': 0.8
        }
    values = np.zeros((len(sets['res']), len(sets['au']), len(sets['h'])))
    
    res_index = {r: idx for idx, r in enumerate(sets['res'])}
    au_index = {u: idx for idx, u in enumerate(sets['au'])}
    
    for u in sets['au']:
        if u not in Plants_res.index:
            logging.error('Reserve not valid for plant ' + u)
            sys.exit(1)
    
        i = au_index[u]
        tech = Plants_res.loc[u, 'Technology']
        droop = Plants_merged.loc[u, 'Droop']
    
        for r in sets['res']:
            j = res_index[r]
    
            eligible = False
            if r in sets['res_U']:
                if tech in commons['tech_batteries'] and r in sets['res_U'][:2]:     # FFRU, PFRU
                    eligible = True
                elif tech in commons['tech_conventional'] and r in sets['res_U'][1:]:  # PFRU, 2U, RR
                    eligible = True
                elif tech in commons['tech_renewables'] and r in sets['res_U'][2:]:   # 2U, RR
                    eligible = True
    
            elif r in sets['res_D']:
                if tech in commons['tech_batteries'] and r in sets['res_D'][:2]:     # FFRD, PFRD
                    eligible = True
                elif tech in commons['tech_conventional'] and r in sets['res_D'][1:]:  # PFRD, 2D
                    eligible = True
                elif tech in commons['tech_renewables'] and r in sets['res_D'][2:]:   # 2D
                    eligible = True
    
            # Si es elegible, calcular participación
            if eligible:
                if r in ['PFRU', 'FFRU'] and droop > 0:
                    factor = (1 / (droop * constants['SystemFrequency'])) * constants['DeltaFrequencyMax']
                    values[j, i, :] = factor
                else:
                    values[j, i, :] = 1  # Participación binaria para otras reservas
        
    parameters['ReserveParticipation'] = {'sets': sets_param['ReserveParticipation'], 'val': values}
   
    # parameters['Reserve'] = define_parameter(['au', 'res'], sets, value=0)
    # res_index = {r: idx for idx, r in enumerate(sets['res'])}
    
    # for i, u in enumerate(sets['au']):
    #     if u not in Plants_res.index:
    #         logging.error('Reserve not valid for plant ' + u)
    #         sys.exit(1)
    
    #     tech = Plants_res.loc[u, 'Technology']
    
    #     for r in sets['res']:
    #         j = res_index[r]
            
    #         if r in sets['res_U']:
    #             if tech in commons['tech_batteries'] and r in sets['res_U'][:2]:  # FFRU, PFRU
    #                 parameters['Reserve']['val'][i, j] = 1
    #             elif tech in commons['tech_conventional'] and r in sets['res_U'][1:]:  # PFRU, 2U, RR
    #                 parameters['Reserve']['val'][i, j] = 1
    #             elif tech in commons['tech_renewables'] and r in sets['res_U'][2:]:  # 2U, RR
    #                 parameters['Reserve']['val'][i, j] = 1
    
    #         elif r in sets['res_D']:
    #             if tech in commons['tech_batteries'] and r in sets['res_D'][:2]:  # FFRD, PFRD
    #                 parameters['Reserve']['val'][i, j] = 1
    #             elif tech in commons['tech_conventional'] and r in sets['res_D'][1:]:  # PFRD, 2D
    #                 parameters['Reserve']['val'][i, j] = 1
    #             elif tech in commons['tech_renewables'] and r in sets['res_D'][2:]:  # 2D
    #                 parameters['Reserve']['val'][i, j] = 1         
                      
    # Technologies
    for unit in range(Nunits):
        idx = sets['t'].index(Plants_merged['Technology'].iloc[unit])
        parameters['Technology']['val'][unit, idx] = True

    # Fuels
    for unit in range(Nunits):
        idx = sets['f'].index(Plants_merged['Fuel'].iloc[unit])
        parameters['Fuel']['val'][unit, idx] = True

    # Location
    for i in range(len(sets['n'])):
        parameters['Location']['val'][:, i] = (Plants_merged['Zone'] == config['zones'][i]).values

    sectors = [col for col in plants if col.startswith('Sector') and not col.startswith('SectorX')]
    for i in range(len(sets['nx'])):
        for s in sectors:
            parameters['LocationX']['val'][:, i] = np.logical_or(parameters['LocationX']['val'][:, i],
                                                                    (Plants_merged[s] == zones_bs[i]).values)
    # Adding the CHP to their respective boundary sectors:
    for i in range(len(sets['nx'])):
        # Use logical OR to combine CHP values with existing boundary sector values
        parameters['LocationX']['val'][:, i] = np.logical_or(parameters['LocationX']['val'][:, i],
                                                             (Plants_merged['Zone_th'] == zones_bs[i]).values)

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
        Plants_merged['InitialPower'] = finalTS['AvailabilityFactors'].iloc[0].multiply(Plants_merged['PowerCapacity'])

        for z in config['zones']:
            tmp_units = Plants_merged.loc[Plants_merged['Zone'] == z].copy()
            tmp_load = finalTS['Load'].iloc[0].loc[z]
            power_initial_res = tmp_units.loc[tmp_units['Fuel'].isin(['WAT', 'WIN', 'SUN']) &
                                              tmp_units['Technology'].isin(['HROR', 'WTON', 'WTOF', 'PHOT'])]['InitialPower'].sum()
            if tmp_load > power_initial_res:
                lack_of_power = tmp_load - power_initial_res
                for f in ['NUC', 'GAS', 'HRD', 'LIG', 'BIO', 'OIL']:
                    for t in ['COMC', 'STUR', 'GTUR', 'ICEN']:
                        n_units = tmp_units.loc[tmp_units['Fuel'].isin([f]) & tmp_units['Technology'].isin([t])].shape[0]
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
    Plants_merged.loc[Plants_merged['Technology'].isin(commons['tech_p2bs']), 'InitialPower'] = 0
    
    if 'InitialPower' in Plants_merged:
        parameters['PowerInitial']['val'] = Plants_merged['InitialPower'].values

    # Config variables:
    sets['x_config'] = ['FirstDay', 'LastDay', 'RollingHorizon Length', 'RollingHorizon LookAhead',
                        'SimulationTimeStep', 'ValueOfLostLoad', 'QuickStartShare', 'WaterValue',
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

    # Determine GDX and Pickle filenames based on MTS flag
    if MTS:
        gdx_in_filename = "Inputs_MTS.gdx" # Use capitalized 'I'
        gdx_out_filename = "Results_MTS.gdx"
        pickle_filename = "Inputs_MTS.p"
    else:
        gdx_in_filename = "Inputs.gdx"
        gdx_out_filename = "Results.gdx"
        pickle_filename = "Inputs.p"

    # Write GDX file with the determined input name
    write_variables(config, gdx_in_filename, [sets, parameters]) # Use gdx_in_filename

    # if the sim variable was not defined:
    if 'sim' not in locals():
        logging.error('Please provide a path where to store the DispaSET inputs (in the "sim" variable)')
        sys.exit(1)

    if not os.path.exists(sim):
        os.makedirs(sim)

    # Determine target GAMS filename
    if MTS:
        target_gms_filename = 'UCM_MTS.gms'
    else:
        target_gms_filename = 'UCM_h.gms'
        
    # Prepare GAMS script modifications dictionary
    gms_modifications = {}
    set_input_prefix = '$set InputFileName'  # Use prefix for robustness
    exec_unload_line = 'EXECUTE_UNLOAD "Results.gdx"' # Use exact match

    # Add base modifications for input/output files (using the now potentially capitalized gdx_in_filename)
    gms_modifications[set_input_prefix] = f'{set_input_prefix} {gdx_in_filename}'
    gms_modifications[exec_unload_line] = f'EXECUTE_UNLOAD "{gdx_out_filename}"'
    
    # Add specific modifications based on MTS, LP, and grid_flag
    if MTS:
        # MTS requires LP
        if not LP:
            logging.error('Simulation in MTS must be LP')
            sys.exit(1)
        gms_modifications['$setglobal LPFormulation 0'] = '$setglobal LPFormulation 1'
        gms_modifications['$setglobal MTS 0'] = '$setglobal MTS 1'
    else: # Detailed dispatch run (not MTS)
        if LP:
            gms_modifications['$setglobal LPFormulation 0'] = '$setglobal LPFormulation 1'
        # For MILP (not LP), ensure MTS flag is off (it should be by default, but explicit)
        elif '$setglobal MTS 0' not in gms_modifications: # Avoid overriding if somehow set
             gms_modifications['$setglobal MTS 0'] = '$setglobal MTS 0' # Explicitly keep default
             
    # Process and write the target GAMS file
    try:
        source_gms_path = os.path.join(GMS_FOLDER, 'UCM.gms')
        target_gms_path = os.path.join(sim, target_gms_filename)
        with open(source_gms_path, 'r') as fin, open(target_gms_path, 'wt') as fout:
            logging.info(f'Processing {source_gms_path} into {target_gms_path}')
            for line in fin:
                stripped_line = line.strip()
                modified = False
                
                # Handle '$set InputFileName' using startswith for robustness
                if stripped_line.startswith(set_input_prefix):
                    fout.write(gms_modifications[set_input_prefix] + '\n')
                    modified = True
                # Handle 'EXECUTE_UNLOAD' using exact match
                elif stripped_line == exec_unload_line:
                    fout.write(gms_modifications[exec_unload_line] + '\n')
                    modified = True
                else:
                    # Handle other modifications (like $setglobal) using exact match
                    for key, value in gms_modifications.items():
                        # Avoid reapplying the ones handled above
                        if key not in [set_input_prefix, exec_unload_line]:
                            if stripped_line == key:
                                fout.write(value + '\n')
                                modified = True
                                break # Apply only the first matching modification
                
                # If no modification was applied, write the original line
                if not modified:
                    fout.write(line)
    except FileNotFoundError:
        logging.error(f"Source GAMS file not found at {source_gms_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error processing GAMS file: {e}")
        sys.exit(1)

    # Create cplex option file
    if config['OptimalityGap'] == '':
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
    else:
        cplex_options = {'epgap': float(config['OptimalityGap']),  # TODO: For the moment hardcoded, it has to be moved to a config file
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

    lines_to_write = ['{} {}'.format(k, v) for k, v in cplex_options.items()]
    with open(os.path.join(sim, 'cplex.opt'), 'w') as f:
        for line in lines_to_write:
            f.write(line + '\n')

    logging.debug('Using gams file from ' + GMS_FOLDER)
    # Copy the correct input GDX file and remove the temporary one created by write_variables
    shutil.copy(gdx_in_filename, os.path.join(sim, gdx_in_filename)) # Use gdx_in_filename for copy destination
    os.remove(gdx_in_filename) # Use gdx_in_filename for removal

    try:
        import cPickle as pickle
    except ImportError:
        import pickle
    # Save pickle file with the correct name
    with open(os.path.join(sim, pickle_filename), 'wb') as pfile: # Use conditional pickle_filename
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
