"""
This files gathers different functions used in the DispaSET to check the input
data

__author__ = 'Sylvain Quoilin (sylvain.quoilin@ec.europa.eu)'
"""

import logging
import os
import sys

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

from ..common import commons  # Load fuel types, technologies, timestep, etc


def isVRE(tech):
    """
    Function that returns true the technology is a variable renewable energy technology
    """
    return tech in commons['tech_renewables']


def isStorage(tech):
    """
    Function that returns true the technology is a storage technology
    """
    return tech in commons['tech_storage']


def check_AvailabilityFactors(plants, AF):
    """
    Function that checks the validity of the provided availability factors and warns
    if a default value of 100% is used.
    """
    RES = commons['tech_renewables']
    for i, v in plants.iterrows():
        u = v['Unit']
        t = v['Technology']
        if t in RES and u not in AF:
            logging.error('Unit ' + str(u) + ' (technology ' + t + ') does not appear in the availbilityFactors table. '
                          'Please provide')
            raise ValueError('Please provide RES AF timeseries for ' + str(u))
        if u in AF:
            if pd.isna(AF[u]).any():
                Nna = pd.isna(AF[u]).count()
                logging.warning('The Availability factor of unit {} for technology {} contains {} '
                                'empty values.'.format(str(u), t, Nna))
            df_af = AF[u].dropna()
            if (df_af == 1).all(axis=None):
                logging.debug('The availability factor of unit ' + str(u) + ' + for technology ' + t +
                              ' is always 100%!')
            if ((df_af < 0) | (df_af > 1)).any(axis=None):
                Nup = df_af[df_af > 1].count()
                Ndo = df_af[df_af < 0].count()
                logging.error('The Availability factor of unit {} for technology {} should be between 0 and 1. '
                              'There are {} values above 1.0 and {} below 0.0'.format(str(u), t, Nup, Ndo))
        else:
            logging.error('Unit ' + str(u) + ' (technology ' + t + ') does not appear in the availbilityFactors table. '
                          'Its values will be set to 100%!')


def check_FlexibleDemand(flex):
    """
    Function that checks the validity of the provided flexibility demand time series
    """
    if (flex.dropna().values < 0).any():
        logging.error('Some flexibility demand values are negative. They must be comprised between 0 and 1')
        sys.exit(1)
    if (flex.dropna().values > 1).any():
        logging.error('Some flexibility demand values are more than 1. They must be comprised between 0 and 1')
        sys.exit(1)


def check_clustering(plants, plants_merged):
    """
    Function that checks that the installed capacities are still equal after the clustering process

    :param plants:  Non-clustered list of units
    :param plants_merged:  clustered list of units
    """
    # First, list all pairs of technology - fuel
    techs = pd.DataFrame([[plants.Technology[idx], plants.Fuel[idx]] for idx in plants.index])
    techs.drop_duplicates(inplace=True)
    for i in techs.index:
        tech = (techs.loc[i, 0], techs.loc[i, 1])
        units_old = plants[(plants.Technology == tech[0]) & (plants.Fuel == tech[1])]
        units_new = plants_merged[(plants_merged.Technology == tech[0]) & (plants_merged.Fuel == tech[1])]
        P_old = (units_old.PowerCapacity * units_old.Nunits).sum()
        P_new = (units_new.PowerCapacity * units_new.Nunits).sum()
        if np.abs(P_old - P_new) / (P_old + 0.0001) > 0.01:
            logging.error('The installed capacity for technology "' + tech[0] + '" and fuel "' + tech[1] +
                          '" is not equal between the original units table (P = ' + str(P_old) +
                          ') and the clustered table (P = ' + str(P_new) + ')')
            sys.exit(1)
    # Check the overall installed storage capacity:
    List_tech_storage = commons['tech_storage']
    isstorage = pd.Series(index=plants.index, dtype='bool')
    for u in isstorage.index:
        isstorage[u] = plants.Technology[u] in List_tech_storage
    isstorage_merged = pd.Series(index=plants_merged.index, dtype='bool')
    for u in isstorage_merged.index:
        isstorage_merged[u] = plants_merged.Technology[u] in List_tech_storage
    TotalStorage = (plants.STOCapacity[isstorage] * plants.Nunits[isstorage]).sum()
    TotalStorage_merged = (plants_merged.STOCapacity[isstorage_merged] * plants_merged.Nunits[isstorage_merged]).sum()
    if np.abs(TotalStorage - TotalStorage_merged) / (TotalStorage + 0.0001) > 0.01:
        logging.error('The total installed storage capacity is not equal between the original units table (' +
                      str(TotalStorage) + ') and the clustered table (' + str(TotalStorage_merged) + ')')
        # sys.exit(1)
    return True


def check_MinMaxFlows(df_min, df_max):
    """
    Function that checks that there is no incompatibility between the minimum and maximum flows
    """
    if (df_min > df_max).any():
        pos = np.where(df_min > df_max)
        logging.critical('At least one minimum flow is higher than the maximum flow, for example in line number ' +
                         str(pos[0][0]) + ' and time step ' + str(pos[1][0]))
        sys.exit(1)

    if (df_max < 0).any():
        pos = np.where(df_max < 0)
        logging.critical('At least one maximum flow is negative, for example in line number ' + str(pos[0][0]) +
                         ' and time step ' + str(pos[1][0]))
        sys.exit(1)

    return True


def check_NonNaNKeys(plants, NonNaNKeys):
    """
    Checking if keys are of type NonNaN

    :param plants:      plants dataframe
    :param NonNaNKeys:  list of NonNaN keys
    """
    for key in NonNaNKeys:
        for u in plants.index:
            if 'Unit' in plants:
                unitname = plants.loc[u, 'Unit']
            else:
                unitname = str(u)
            if isinstance(plants.loc[u, key], str):
                logging.critical('A non numeric value was detected in the power plants inputs for parameter "' + key +
                                 '"')
                sys.exit(1)
            if np.isnan(plants.loc[u, key]):
                logging.critical('The power plants data is missing for unit ' + unitname + ' and parameter "' + key +
                                 '"')
                sys.exit(1)


def check_StrKeys(plants, StrKeys):
    """
    Checking if keys are of type Str

    :param plants:      plants dataframe
    :param StrKeys:     list of Str keys
    """
    for key in StrKeys:
        for u in plants.index:
            if 'Unit' in plants:
                unitname = plants.loc[u, 'Unit']
            else:
                unitname = str(u)
            if not isinstance(plants.loc[u, key], str):
                logging.critical('A numeric value was detected in the power plants inputs for parameter "' + key +
                                 '". This column should contain strings only.')
                sys.exit(1)
            elif plants.loc[u, key] == '':
                logging.critical('An empty value was detected in the power plants inputs for unit "' + unitname +
                                 '" and parameter "' + key + '"')
                sys.exit(1)


def check_keys(plants, keys, unit):
    """
    Checking mandatory keys

    :param plants:      plants dataframe
    :param keys:        list of keys
    :param unit:        string denoting type of units being checked
    """
    for key in keys:
        if key not in plants:
            logging.critical('The power plants data does not contain the field "' + key +
                             '", which is mandatory for ' + unit + ' units')
            sys.exit(1)


def check_sto(config, plants, raw_data=True):
    """
    Function that checks the storage plant characteristics
    """
    if raw_data:
        keys = ['STOCapacity', 'STOSelfDischarge', 'STOMaxChargingPower', 'STOChargingEfficiency']
        NonNaNKeys = ['STOCapacity']
    else:
        keys = ['StorageCapacity', 'StorageSelfDischarge', 'StorageChargingCapacity', 'StorageChargingEfficiency']
        NonNaNKeys = ['StorageCapacity']

    if 'StorageInitial' in plants:
        logging.warning('The "StorageInitial" column is present in the power plant table, although it is deprecated '
                        '(it should now be defined in the ReservoirLevel data table). It will not be considered.')

    check_keys(plants, keys, 'Storage')
    check_NonNaNKeys(plants,NonNaNKeys)

    if raw_data:
        key = 'STOCapacity'
        P_charging = 'STOMaxChargingPower'
    else:
        key = 'StorageCapacity'
        P_charging = 'StorageChargingCapacity'
    for u in plants.index:
        maxpower = max(plants.loc[u, 'PowerCapacity'], plants.loc[u, P_charging])
        if plants.loc[u, key] > maxpower * 8760:
            logging.error('The Storage capacity for unit ' + plants.loc[u, 'Unit'] + ' is prohibitively high. '
                          'More than one year at full power is required to discharge the reservoir')
        elif plants.loc[u, key] > maxpower * 3000:
            logging.warning('The Storage capacity for unit ' + plants.loc[u, 'Unit'] + ' is very high.')
        elif (plants.loc[u, key] > maxpower * 24 * config['HorizonLength'] / config['SimulationTimeStep']) and (
                config['HydroScheduling'] not in ['Zonal', 'Regional']):
            logging.warning('The Storage capacity for unit ' + plants.loc[u, 'Unit'] +
                            ' is high. Make sure to provide a proper storage level profile')

    return True


def check_p2bs(config, plants):
    """
    Function that checks the p2bs unit characteristics
    """
    keys = [col for col in plants if col.startswith('ChargingEfficiencySector') or col.startswith('EfficiencySector')]
    NonNaNKeys = []
    StrKeys = [col for col in plants if col.startswith('Sector') and not col.startswith('SectorX')]

    if len(plants) == 0:  # If there are no P2HT units, exit the check
        return True

    for t in commons['tech_p2bs']:
        check_keys(plants, keys, t)
    check_NonNaNKeys(plants, NonNaNKeys)
    check_StrKeys(plants, StrKeys)

    # Check the COP values:
    #TODO: check variables that need to be checked (Eff, sector, power capacity etc.)

    # for u in plants.index:
    #     if plants.loc[u, 'Efficiency1'] < 0 or plants.loc[u, 'Efficiency2'] < 0:
    #         logging.critical('The Efficiency value of p2bs units must be >= 0. '
    #                          'The provided value for Efficiency1 of unit ' + u + ' is "' +
    #                          str(plants.loc[u, 'Efficiency1'] + '"')
    #                          )
    #         sys.exit(1)

    return True


def check_boundary_sector(config, plants, BoundarySector=None):
    """
    Check boundary sector units characteristics
    :param config: Config dictionary
    :param plants: DataFrame with plants data
    :param BoundarySector: Optional DataFrame with boundary sector data for consistency check
    """
    keys = ['PowerCapacity', 'Efficiency']
    NonNaNKeys = []
    StrKeys = ['Sector1','Sector2']

    if len(plants) == 0:  # If there are no boundary sector units, exit the check
        return True

    # Check basic characteristics
    for t in commons['tech_boundary_sector']:
        check_keys(plants, keys, t)
    check_NonNaNKeys(plants, NonNaNKeys)
    check_StrKeys(plants, StrKeys)

    # Check boundary sector name consistency
    if BoundarySector is not None:
        # Get valid sector names from BoundarySector table
        valid_sectors = BoundarySector.index.tolist()
        
        # Check each sector column
        for col in ['Sector1', 'Sector2']:
            if col in plants.columns:
                # Get unique non-empty sector values, excluding both pandas nan and string "nan"
                sectors = plants[col].dropna().unique()
                sectors = [s for s in sectors if str(s).strip() != '' and str(s).lower() != 'nan']
                
                # Check each sector value
                for sector in sectors:
                    if sector not in valid_sectors:
                        logging.critical('Boundary sector "{}" found in {} column of plants table is not defined in the BoundarySector table. Valid sectors are: {}'.format(
                            sector, col, valid_sectors))
                        sys.exit(1)

    return True


def check_chp(config, plants):
    """
    Function that checks the CHP plant characteristics
    """
    keys = ['CHPType', 'CHPPowerToHeat', 'CHPPowerLossFactor']
    NonNaNKeys = ['CHPPowerToHeat', 'CHPPowerLossFactor']
    StrKeys = ['CHPType']

    check_keys(plants, keys, 'CHP')
    check_NonNaNKeys(plants, NonNaNKeys)
    check_StrKeys(plants, StrKeys)

    # Check the efficiency values:
    unitname = str()
    for u in plants.index:
        if 'Unit' in plants:
            unitname = plants.loc[u, 'Unit']
        else:
            unitname = str(u)
        plant_PowerCapacity = plants.loc[u, 'PowerCapacity']
        plant_MaxHeat = plants.loc[u, 'CHPMaxHeat']
        plant_powertoheat = plants.loc[u, 'CHPPowerToHeat']
        plant_powerlossfactor = plants.loc[u, 'CHPPowerLossFactor']

        if plants.loc[u, 'CHPType'].lower() not in ['extraction', 'back-pressure', 'p2h']:
            logging.critical('The value of CHPType should be "extraction", "back-pressure" or "p2h". '
                             'The type of unit ' + u + ' is "' + str(plants.loc[u, 'CHPType'] + '"'))
            sys.exit(1)
        if 0 > plant_powertoheat > 10:
            logging.critical('The value of CHPPowerToHeat should be higher or equal to zero and lower than 10. '
                             'Unit ' + u + ' has a value of ' + str(plant_powertoheat))
            sys.exit(1)
        if 0 > plant_powerlossfactor > 1 and plants.loc[u, 'CHPType'].lower() != 'p2h':
            logging.critical('The value of CHPPowerLossFactor should be higher or equal to zero and lower than 1. '
                             'Unit ' + u + ' has a value of ' + str(plant_powerlossfactor))
            sys.exit(1)
        if plants.loc[u, 'CHPType'].lower() == 'back-pressure' and plant_powerlossfactor != 0:
            logging.critical('The value of CHPPowerLossFactor must be zero if the CHP types is "back-pressure". '
                             'Unit ' + u + ' has a value of ' + str(plant_powerlossfactor))
            sys.exit(1)
        if plants.loc[u, 'CHPType'].lower() == 'extraction':
            intersection_MaxHeat = plant_PowerCapacity / plant_powertoheat
            if not pd.isnull(plant_MaxHeat):
                if intersection_MaxHeat < plant_MaxHeat:
                    logging.warning('Given Maximum heat CHPMaxHeat ({}) is higher than the intersection point of the '
                                    'two other constraints ({}) (power loss factor and backpressure line) therefore it '
                                    'will not be ignored'.format(plant_MaxHeat, intersection_MaxHeat))
                    plant_MaxHeat = intersection_MaxHeat
            else:
                plant_MaxHeat = intersection_MaxHeat

        # Calculating the nominal total efficiency at the highest point:
        if plants.loc[u, 'CHPType'].lower() != 'p2h':
            Fuel = (plant_PowerCapacity + plant_powerlossfactor * plant_MaxHeat) / plants.loc[
                u, 'Efficiency']  # F = (P + C_v * Q)/eta_condensation
            TotalEfficiency = (plant_PowerCapacity + plant_MaxHeat) / Fuel  # eta_tot = (P + Q) / F
            logging.debug('Highest overall efficiency of CHP plant {} is {:.2f}'.format(u, TotalEfficiency))
            if TotalEfficiency < 0 or TotalEfficiency > 1.14:
                logging.critical('The calculated value of the total CHP efficiency for unit ' + unitname + ' is ' +
                                 str(TotalEfficiency) + ', which is unrealistic!')
                sys.exit(1)
            if TotalEfficiency > 0.95:
                logging.warning('The calculated value of the total CHP efficiency for unit ' + unitname + ' is ' +
                                str(TotalEfficiency) + ', which is very high!')

    # Check the optional MaxHeatCapacity parameter. While it adds another realistic boundary it is not a required
    # parameter for the definition of the CHP's operational envelope.:
    if 'CHPMaxHeat' in plants:
        for u in plants.index:
            plant_MaxHeat = plants.loc[u, 'CHPMaxHeat']
            if plant_MaxHeat <= 0:
                logging.warning('CHPMaxHeat for plant {} is {} which shuts down any heat '
                                'production.'.format(u, plant_MaxHeat))
    # Check the optional heat storage values:
    if 'STOCapacity' in plants:
        for u in plants.index:
            Qdot = plants.loc[u, 'PowerCapacity'] / plants.loc[u, 'CHPPowerToHeat']
            if plants.loc[u, 'STOCapacity'] < Qdot * 0.5:
                logging.warning('Unit ' + unitname + ': The value of the thermal storage capacity (' +
                                str(plants.loc[u, 'STOCapacity']) + 'MWh) seems very low compared to its '
                                'thermal power (' + str(Qdot) + 'MW).')
            elif plants.loc[u, 'STOCapacity'] > Qdot * 24:
                logging.warning('Unit ' + unitname + ': The value of the thermal storage capacity (' + str(
                                plants.loc[u, 'STOCapacity']) + 'MWh) seems very high compared to its thermal power (' +
                                str(Qdot) + 'MW).')

    if 'STOSelfDischarge' in plants:
        for u in plants.index:
            if plants.loc[u, 'STOSelfDischarge'] < 0:
                logging.error('Unit ' + unitname + ': The value of the thermal storage self-discharge (' +
                              str(plants.loc[u, 'STOSelfDischarge'] * 100) + '%/day) cannot be negative')
                sys.exit(1)
            elif plants.loc[u, 'STOSelfDischarge'] > 1:
                logging.warning('Unit ' + unitname + ': The value of the thermal storage self-discharge (' +
                                str(plants.loc[u, 'STOSelfDischarge'] * 100) + '%/day) seems very high')
            elif plants.loc[u, 'STOSelfDischarge'] > 24:
                logging.error('Unit ' + unitname + ': The value of the thermal storage self-discharge (' +
                              str(plants.loc[u, 'STOSelfDischarge'] * 100) + '%/day) is too high')
                sys.exit(1)

    return True


def check_units(config, plants):
    """
    Function that checks the power plant characteristics
    """

    keys = ['Unit', 'Fuel', 'Zone', 'Technology', 'PowerCapacity', 'PartLoadMin', 'RampUpRate', 'RampDownRate',
            'StartUpTime', 'MinUpTime', 'MinDownTime', 'NoLoadCost', 'StartUpCost', 'Efficiency', 'CO2Intensity',
            'WaterWithdrawal', 'WaterConsumption','InertiaConstant']
    NonNaNKeys = ['PowerCapacity', 'PartLoadMin', 'RampUpRate', 'RampDownRate', 'Efficiency', 'RampingCost',
                  'CO2Intensity']
    StrKeys = ['Unit', 'Zone', 'Fuel', 'Technology']

    # Special treatment for the Optional key Nunits:
    if 'Nunits' in plants:
        keys.append('Nunits')
        NonNaNKeys.append('Nunits')
        if any([not float(x).is_integer() for x in plants['Nunits']]):
            logging.error('Some values are not integers in the "Nunits" column of the plant database')
            sys.exit(1)
    else:
        logging.info('The columns "Nunits" is not present in the power plant database. '
                     'A value of one will be assumed by default')

    check_keys(plants, keys, 'all')
    check_NonNaNKeys(plants, NonNaNKeys)
    check_StrKeys(plants, StrKeys)

    lower = {'PowerCapacity': 0, 'PartLoadMin': 0, 'StartUpTime': 0, 'MinUpTime': 0, 'MinDownTime': 0, 'NoLoadCost': 0,
             'StartUpCost': 0, 'WaterWithdrawal': 0, 'WaterConsumption': 0}
    strictly_lower = {'RampUpRate': 0, 'RampDownRate': 0, 'Efficiency': 0}
    higher = {'PartLoadMin': 1, 'Efficiency': 1}
    higher_time = {'MinUpTime': 0, 'MinDownTime': 0}  # 'StartUpTime':0,

    # Special treatment for the Optional key Nunits:
    if 'Nunits' in plants:
        strictly_lower['Nunits'] = 0

    if len(plants['Unit'].unique()) != len(plants['Unit']):
        duplicates = plants['Unit'][plants['Unit'].duplicated()].tolist()
        logging.error('The names of the power plants are not unique. The following names are duplicates: ' +
                      str(duplicates) + '. "' + str(duplicates[0] + '" appears for example in the following zones: ' +
                      str(plants.Zone[plants['Unit'] == duplicates[0]].tolist())))
        sys.exit(1)

    for key in lower:
        if any(plants[key] < lower[key]):
            plantlist = plants[plants[key] < lower[key]]
            plantlist = plantlist['Unit'].tolist()
            logging.critical('The value of ' + key + ' should be higher or equal to zero. A negative value has been '
                             'found for units ' + str(plantlist))
            sys.exit(1)

    for key in strictly_lower:
        if any(plants[key] <= strictly_lower[key]):
            plantlist = plants[plants[key] <= strictly_lower[key]]
            plantlist = plantlist['Unit'].tolist()
            logging.critical('The value of ' + key + ' should be strictly higher than zero. '
                             'A null or negative value has been found for units ' + str(plantlist))
            sys.exit(1)

    for key in higher:
        if any(plants[key] > higher[key]):
            plantlist = plants[plants[key] > higher[key]]
            plantlist = plantlist[~plantlist['Technology'].str.contains("ABHP")]
            if not plantlist.empty:
                plantlist = plantlist['Unit'].tolist()
                logging.critical('The value of ' + key + ' should be lower or equal to one. '
                                 'A higher value has been found for units ' + str(plantlist))
                sys.exit(1)

    for key in higher_time:
        if any(plants[key] >= config['HorizonLength'] * 24):
            plantlist = plants[plants[key] >= config['HorizonLength'] * 24]
            plantlist = plantlist['Unit'].tolist()
            logging.critical('The value of ' + key + ' should be lower than the horizon length (' +
                             str(config['HorizonLength'] * 24) + ' hours). A higher value has been found for units ' +
                             str(plantlist))
            sys.exit(1)

    # Checking che compatibility between the selected simulation time and the power plant constraints:
    if config['SimulationType'] in ('LP', 'LP clustered'):
        for key in ['NoLoadCost', 'PartLoadMin', 'MinEfficiency', 'StartUpTime']:
            if (plants[key] > 0).any():
                logging.error('Non-null value(s) have been found for key ' + key + ' in the power plant list. '
                              'This cannot be modelled with the ' + config['SimulationType'] + ' formulation and '
                              'will therefore not be considered.')
    return True


def check_heat_demand(plants, data, zones_th):
    """
    Function that checks the validity of the heat demand profiles
    :param     plants:  List of plants
    :param     data: Dataframe with the heat demand time series
    :param     zones_th: list with the heating zones
    """
    plants.index = plants['Unit']
    plants_heating = plants[
        [str(plants['CHPType'][u]).lower() in commons['types_CHP'] or plants.loc[u, 'Technology'] == 'P2HT' for u in
         plants.index]]
    plants_chp = plants[[str(plants['CHPType'][u]).lower() in commons['types_CHP'] for u in plants.index]]

    for z in data:  # for each heating zone in the heating demand data
        if z not in zones_th:
            logging.error('The heat demand profile with header "' + str(
                z) + '" does not correspond to any heating zone. It will be ignored.')
        elif (data[z] == 0).all():
            logging.error('Heat demand data for zone "' + z + '" is either no found or always equal to zero')
        if z in plants_chp.index:  # special case in which the heating zone corresponds to a single CHP unit
            u = z
            Nunits = plants.loc[u, 'Nunits']
            plant_CHP_type = plants.loc[u, 'CHPType'].lower()
            if pd.isnull(plants.loc[u, 'CHPMaxHeat']):
                plant_Qmax = +np.inf
            else:
                plant_Qmax = plants.loc[u, 'CHPMaxHeat']
            if plant_CHP_type == 'extraction':
                Qmin = 0
                Qmax = min(plants.loc[u, 'PowerCapacity'] / plants.loc[u, 'CHPPowerToHeat'], plant_Qmax) * Nunits
            elif plant_CHP_type == 'back-pressure':
                Qmin = plants.loc[u, 'PowerCapacity'] * plants.loc[u, 'PartLoadMin'] / plants.loc[u, 'CHPPowerToHeat']
                Qmax = min(plants.loc[u, 'PowerCapacity'] / plants.loc[u, 'CHPPowerToHeat'], plant_Qmax) * Nunits
            elif plant_CHP_type == 'p2h':
                Qmin = 0
                Qmax = plant_Qmax * Nunits
            else:
                logging.error('The CHP type for unit ' + u + ' is not valid.')
            if np.isnan(Qmax) and plant_CHP_type != 'p2h':
                logging.error(
                    'CHPPowerToHeat is not defined for unit ' + str(u) + ' appearing in the heat demand profiles')
                sys.exit(1)
            elif data[u].max() > Qmax:
                logging.warning('The maximum thermal demand for unit ' + str(u) + ' (' + str(
                    data[u].max()) + ') is higher than its thermal capacity (' + str(
                    Qmax) + '). Slack heat will be used to cover that.')
            if data[u].min() < Qmin:
                logging.warning('The minimum thermal demand for unit ' + str(u) + ' (' + str(
                    data[u].min()) + ') is lower than its minimum thermal generation (' + str(Qmin) + ' MWth)')

    # check that a heating demand has been provided for all heating zones
    for z in zones_th:
        if z not in data:
            logging.critical('No heat demand data was found for thermal zone ' + z)
            sys.exit(1)

    return True


def check_reserves(Reserve2D, Reserve2U, Load):
    """
    Function that checks the validity of the reserve requirement time series
    :param Reserve2D:   DataFrame of reserves 2D
    :param Reserve2U:   DataFrame of reserves 2U
    :param Load:        DataFrame of Loads
    """
    for z in Load.columns:
        if z in Reserve2U:
            if (Reserve2U[z] < 0).any():
                logging.critical('The reserve 2U table contains negative values for zone ' + z)
                sys.exit(1)
            if (Load[z] - Reserve2U[z] < 0).any():
                logging.critical('The reserve 2U table contains negative values higher than demand for zone ' + z)
                sys.exit(1)
        else:
            logging.warning('No 2U reserve requirement data has been found for zone ' + z +
                            '. Using the standard formula')
        if z in Reserve2D:
            if (Reserve2D[z] < 0).any():
                logging.critical('The reserve 2D table contains negative values for zone ' + z)
                sys.exit(1)
            if (Load[z] - Reserve2D[z] < 0).any():
                logging.critical('The reserve 2D table contains values higher than demand for zone ' + z)
                sys.exit(1)
        else:
            logging.warning('No 2D reserve requirement data has been found for zone ' + z +
                            '. Using the standard formula')

def check_FFRLimit(FFRLimit, Load):
    """
    Function that checks the validity of the reserve requirement time series
    :param FFR:   DataFrame of FFR Limit
    :param Load:        DataFrame of Loads
    """
    if (FFRLimit.sum(axis=1) < 0).any():
        logging.critical('The FFR Limit table contains negative values')
        sys.exit(1)
    if (Load.sum(axis=1) - FFRLimit.sum(axis=1) < 0).any():
        logging.critical('The FFR Limit table contains values higher than demand')
        sys.exit(1)
    else:
        logging.warning('No FFR Limit requirement data has been found')
        
def check_PrimaryReserveLimit(PrimaryReserveLimit, Load):
    """
    Function that checks the validity of the reserve requirement time series
    :param PrimaryReserve:   DataFrame of Primary Reserve Limit
    :param Load:        DataFrame of Loads
    """
    if (PrimaryReserveLimit.sum(axis=1) < 0).any():
        logging.critical('The Primary Reserve Limit table contains negative values')
        sys.exit(1)
    if (Load.sum(axis=1) - PrimaryReserveLimit.sum(axis=1) < 0).any():
        logging.critical('The Primary Reserve Limit table contains values higher than demand')
        sys.exit(1)
    else:
        logging.warning('No Primary Reserve Limit requirement data has been found')


def check_df(df, StartDate=None, StopDate=None, name=''):
    """
    Function that check the time series provided as inputs
    """

    if isinstance(df.index, pd.DatetimeIndex):
        if StartDate not in df.index:
            logging.warning('The start date ' + str(StartDate) + ' is not in the index of the provided dataframe')
        if StopDate not in df.index:
            logging.warning('The stop date ' + str(StopDate) + ' is not in the index of the provided dataframe')
    if any(np.isnan(df)):
        for key in df:
            missing = np.sum(np.isnan(df[key]))
            # pos = np.where(np.isnan(df.sum(axis=1)))
            # idx_pos = [df.index[i] for i in pos]
            if missing > 1:
                logging.warning('There are ' + str(missing) + ' missing entries in the column ' + key +
                                ' of the dataframe ' + name)
    if not df.columns.is_unique:
        logging.error('The column headers of table "' + name + '" are not unique!. '
                      'The following headers are duplicated: ' + str(
                df.columns.get_duplicates()))
        sys.exit(1)
    return True


def check_BSFlexMaxCapacity(parameters, config, sets):
    """
    Function that checks if SectorXFlexMaxCapacity is at least as high as the average flexible demand input
    """
    for i, nx in enumerate(sets['nx']):
        # Calculate average flexible demand
        avg_demand = parameters['SectorXFlexDemandInput']['val'][i, :].mean()
        # Check if max capacity is sufficient
        if parameters['SectorXFlexMaxCapacity']['val'][i] < avg_demand:
            logging.critical(f'Boundary sector {nx}: SectorXFlexMaxCapacity ({parameters["SectorXFlexMaxCapacity"]["val"][i]}) is lower than average flexible demand ({avg_demand})')
            sys.exit(1)


def check_BSFlexMaxSupply(parameters, config, sets):
    """
    Function that checks if SectorXFlexMaxSupply is at least as high as the average flexible supply input
    """
    for i, nx in enumerate(sets['nx']):
        # Calculate average flexible supply
        avg_supply = parameters['SectorXFlexSupplyInput']['val'][i, :].mean()
        # Check if max supply is sufficient
        if parameters['SectorXFlexMaxSupply']['val'][i] < avg_supply:
            logging.critical(f'Boundary sector {nx}: SectorXFlexMaxSupply ({parameters["SectorXFlexMaxSupply"]["val"][i]}) is lower than average flexible supply ({avg_supply})')
            sys.exit(1)

def check_simulation_environment(SimulationPath, store_type='pickle', firstline=7):
    """
    Function to test the validity of disapset inputs
    :param SimulationPath:          Path to the simulation folder
    :param store_type:              choose between: "list", "excel", "pickle"
    :param firstline:               Number of the first line in the data (only if type=='excel')
    """

    import cPickle

    # minimum list of variable required for dispaSET:
    list_sets = [
        'h',
        'd',
        'mk',
        'n',
        'c',
        'p',
        'l',
        'f',
        's',
        't',
        'tr',
        'u']

    list_param = [
        'AvailabilityFactor',
        'CostFixed',
        'CostShutDown',
        'Curtailment',
        'Demand',
        'Efficiency',
        'Fuel',
        'CostVariable',
        'FuelPrice',
        'Markup',
        'CostStartUp',
        'EmissionMaximum',
        'EmissionRate',
        'FlowMaximum',
        'FlowMinimum',
        'LineNode',
        'Location',
        'LoadShedding',
        'OutageFactor',
        'PermitPrice',
        'PriceTransmission',
        'PowerCapacity',
        'PartLoadMin',
        'RampUpMaximum',
        'RampDownMaximum',
        'RampStartUpMaximum',
        'RampShutDownMaximum',
        'Reserve',
        'StorageDischargeEfficiency',
        'StorageCapacity',
        'StorageInflow',
        'StorageOutflow',
        'StorageInitial',
        'StorageMinimum',
        'StorageChargingEfficiency',
        'StorageChargingCapacity',
        'Technology',
        'TimeDownMinimum',
        'TimeUpMinimum',
        'TimeDownInitial',
        'TimeUpInitial',
        'PowerInitial']

    if store_type == 'list':
        if isinstance(SimulationPath, list):
            # The list of sets and parameters has been passed directly to the function, checking that all are present:
            SimulationPath_vars = [SimulationPath[i]['name'] for i in range(len(SimulationPath))]
            for var in list_sets + list_param:
                if var not in SimulationPath_vars:
                    logging.critical('The variable "' + var + '" has not been found in the list of input variables')
                    sys.exit(1)
        else:
            logging.critical('The argument must a list. Please correct or change the "type" argument')
            sys.exit(1)

    elif store_type == 'pickle':
        if os.path.exists(SimulationPath):
            if os.path.isfile(os.path.join(SimulationPath, 'Inputs.p')):
                variables = cPickle.load(open(os.path.join(SimulationPath, 'Inputs.p'), 'rb'))
                arg_vars = [variables[i]['name'] for i in range(len(variables))]
                for var in list_sets + list_param:
                    if var not in arg_vars:
                        logging.critical('Found Pickle file but does not contain valid DispaSET input (' + var +
                                         ' missing)')
                        sys.exit(1)
            else:
                logging.critical('Could not find the Inputs.p file in the specified directory')
                sys.exit(1)
        else:
            logging.critical('The function argument is not a valid directory')
            sys.exit(1)

    elif store_type == 'excel':
        if os.path.exists(SimulationPath):
            if not os.path.isfile(os.path.join(SimulationPath, 'InputDispa-SET - Sets.xlsx')):
                logging.critical("Could not find the file 'InputDispa-SET - Sets.xlsx'")
                sys.exit(1)
            for var in list_param:
                if os.path.isfile(os.path.join(SimulationPath, 'InputDispa-SET - ' + var + '.xlsx')):
                    a = 1
                else:
                    logging.critical("Could not find the file 'InputDispa-SET - " + var + ".xlsx'")
                    sys.exit(1)

        else:
            logging.critical('The function argument is not a valid directory')
            sys.exit(1)

    else:
        logging.critical('The "type" parameter must be one of the following : "list", "excel", "pickle"')
        sys.exit(1)

def check_CostXNotServed(config, CostXNotServed, zones_bs):
    """
    Check that CostXNotServed is properly defined for all boundary sectors
    
    :param config:            Dictionary with all the configuration fields loaded from the excel file
    :param CostXNotServed:    DataFrame with CostXNotServed values
    :param zones_bs:          List of boundary sector zones
    """
    if CostXNotServed.empty:
        logging.critical('CostXNotServed is not defined for any boundary sector')
        sys.exit(1)
    
    for zone in zones_bs:
        if zone not in CostXNotServed.columns:
            logging.critical('CostXNotServed is not defined for boundary sector ' + zone)
            sys.exit(1)
        if CostXNotServed[zone].isna().all():
            logging.critical('CostXNotServed is not defined for boundary sector ' + zone)
            sys.exit(1)
        if (CostXNotServed[zone] == 0).all():
            logging.warning('CostXNotServed is zero for boundary sector ' + zone + '. This may lead to unrealistic results.')
