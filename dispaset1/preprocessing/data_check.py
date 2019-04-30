"""
This files gathers different functions used in the DispaSET to check the input
data

__author__ = 'Sylvain Quoilin (sylvain.quoilin@ec.europa.eu)'
"""

import os
import sys
import numpy as np
import pandas as pd
import logging

def check_AvailabilityFactors(plants,AF):
    '''
    Function that checks the validity of the provided availability factors and warns
    if a default value of 100% is used.
    '''
    if (AF.values < 0).any():
        logging.error('Some Availaibility factors are negative')
        sys.exit(1)
    if (AF.values > 1).any():
        logging.critical('Some Availability factors are higher than one. They must be carefully checked')
    for t in ['WTON', 'WTOF', 'PHOT', 'HROR']:
        for i in plants[plants['Technology']==t].index:
            u = plants.loc[i,'Unit']
            if u in AF:
                if (AF[u].values == 1).all():
                    logging.critical('The availability factor of unit ' + str(u) + ' + for technology ' + t + ' is always 100%!')
            else:
                logging.critical('Unit ' + str(u) + ' (technology ' + t + ') does not appear in the availbilityFactors table. Its values will be set to 100%!')

def check_MinMaxFlows(df_min,df_max):
    '''
    Function that checks that there is no incompatibility between the minimum and maximum flows
    '''
    if (df_min > df_max).any():
        pos = np.where(df_min > df_max)
        logging.critical('ERROR: At least one minimum flow is higher than the maximum flow, for example in line number ' + str(pos[0][0]) + ' and time step ' + str(pos[1][0]))
        sys.exit(1)
        
    if (df_max < 0).any():
        pos = np.where(df_max < 0)
        logging.critical('ERROR: At least one maximum flow is negative, for example in line number ' + str(pos[0][0]) + ' and time step ' + str(pos[1][0]))
        sys.exit(1)

    return True


def check_chp(config, plants):
    """
    Function that checks the CHP plant characteristics
    """   
    keys = ['CHPType','CHPPowerToHeat','CHPPowerLossFactor']
    NonNaNKeys = ['CHPPowerToHeat','CHPPowerLossFactor']
    StrKeys = ['CHPType']
    
    for key in keys:
        if key not in plants:
            logging.critical('The power plants data does not contain the field "' + key + '", which is mandatory for CHP units')
            sys.exit(1)

    for key in NonNaNKeys:
        for u in plants.index:
            if type(plants.loc[u, key]) == str:
                logging.critical('A non numeric value was detected in the power plants inputs for parameter "' + key + '"')
                sys.exit(1)
            if np.isnan(plants.loc[u, key]):
                logging.critical('The power plants data is missing for unit number ' + str(u) + ' and parameter "' + key + '"')
                sys.exit(1)

    for key in StrKeys:
        for u in plants.index:
            if not type(plants.loc[u, key]) == str:
                logging.critical(
                    'A numeric value was detected in the power plants inputs for parameter "' + key + '". This column should contain strings only.')
                sys.exit(1)
            elif plants.loc[u, key] == '':
                logging.critical('An empty value was detected in the power plants inputs for unit "' + str(u) + '" and parameter "' + key + '"')
                sys.exit(1) 
    
    # Check the efficiency values:
    for u in plants.index:
        plant_PowerCapacity = plants.loc[u,'PowerCapacity']
        plant_MaxHeat = plants.loc[u, 'CHPMaxHeat']
        plant_powertoheat =  plants.loc[u,'CHPPowerToHeat']
        plant_powerlossfactor = plants.loc[u,'CHPPowerLossFactor']

        if plants.loc[u,'CHPType'].lower() not in ['extraction','back-pressure', 'p2h']:
            logging.critical('The value of CHPType should be "extraction", "back-pressure" or "p2h". The type of unit ' + u + ' is "' + str(plants.loc[u,'CHPType'] + '"'))
            sys.exit(1)              
        if 0 > plant_powertoheat > 10:
            logging.critical('The value of CHPPowerToHeat should be higher or equal to zero and lower than 10. Unit ' + u + ' has a value of ' + str(plant_powertoheat))
            sys.exit(1)         
        if 0 > plant_powerlossfactor > 1 and plants.loc[u,'CHPType'].lower() != 'p2h':
            logging.critical('The value of CHPPowerLossFactor should be higher or equal to zero and lower than 1. Unit ' + u + ' has a value of ' + str(plant_powerlossfactor))
            sys.exit(1)   
        if plants.loc[u,'CHPType'].lower() == 'back-pressure' and plant_powerlossfactor != 0:
            logging.critical('The value of CHPPowerLossFactor must be zero if the CHP types is "back-pressure". Unit ' + u + ' has a value of ' + str(plant_powerlossfactor))
            sys.exit(1)
        if plants.loc[u, 'CHPType'].lower() == 'extraction':
            intersection_MaxHeat = plant_PowerCapacity / plant_powertoheat
            if not pd.isnull(plant_MaxHeat):
                if intersection_MaxHeat < plant_MaxHeat:
                    logging.warning('Given Maximum heat CHPMaxHeat ({}) is higher than the intersection point of the two other constraints ({}) '
                                    '(power loss factor and backpressure line) therefore it will not be ignored'.format(plant_MaxHeat, intersection_MaxHeat) )
                    plant_MaxHeat = intersection_MaxHeat
            else:
                plant_MaxHeat = intersection_MaxHeat

        # Calculating the nominal total efficiency at the highest point:
        if plants.loc[u,'CHPType'].lower() != 'p2h':
            Fuel = (plant_PowerCapacity + plant_powerlossfactor * plant_MaxHeat)/plants.loc[u,'Efficiency'] # F = (P + C_v * Q)/eta_condensation
            TotalEfficiency = (plant_PowerCapacity + plant_MaxHeat) / Fuel             # eta_tot = (P + Q) / F
            logging.debug('Highest overall efficiency of CHP plant {} is {:.2f}'.format(u,TotalEfficiency))
            if TotalEfficiency < 0 or TotalEfficiency > 1.14:
                logging.critical('The calculated value of the total CHP efficiency for unit ' + u + ' is ' + str(TotalEfficiency) + ', which is unrealistic!')
                sys.exit(1)
            if TotalEfficiency > 0.95:
                logging.warn('The calculated value of the total CHP efficiency for unit ' + u + ' is ' + str(TotalEfficiency) + ', which is very high!')

    # Check the optional MaxHeatCapacity parameter. While it adds another realistic boundary it is not a required parameter for the definition of the CHP's operational envelope.:
    if 'CHPMaxHeat' in plants:
        for u in plants.index:
            plant_MaxHeat = plants.loc[u, 'CHPMaxHeat']
            if plant_MaxHeat <=0:
                logging.warn('CHPMaxHeat for plant {} is {} which shuts down any heat production.'.format(u, plant_MaxHeat))
    # Check the optional heat storage values:
    if 'STOCapacity' in plants:
        for u in plants.index:     
            Qdot = plants.loc[u,'PowerCapacity']/plants.loc[u,'CHPPowerToHeat']
            if plants.loc[u,'STOCapacity'] < Qdot * 0.5 :
                logging.warn('Unit ' + u + ': The value of the thermal storage capacity (' + str(plants.loc[u,'STOCapacity']) + 'MWh) seems very low compared to its thermal power (' + str(Qdot) + 'MW).')
            elif plants.loc[u,'STOCapacity'] > Qdot * 24:
                logging.warn('Unit ' + u + ': The value of the thermal storage capacity (' + str(plants.loc[u,'STOCapacity']) + 'MWh) seems very high compared to its thermal power (' + str(Qdot) + 'MW).')

    if 'STOSelfDischarge' in plants:
        for u in plants.index:     
            if plants.loc[u,'STOSelfDischarge'] < 0 :
                logging.error('Unit ' + u + ': The value of the thermal storage self-discharge (' + str(plants.loc[u,'STOSelfDischarge']*100) + '%/day) cannot be negative')
                sys.exit(1)
            elif plants.loc[u,'STOSelfDischarge'] > 1:
                logging.warn('Unit ' + u + ': The value of the thermal storage self-discharge (' + str(plants.loc[u,'STOSelfDischarge']*100) + '%/day) seems very high')
            elif plants.loc[u,'STOSelfDischarge'] > 24:
                logging.error('Unit ' + u + ': The value of the thermal storage self-discharge (' + str(plants.loc[u,'STOSelfDischarge']*100) + '%/day) is too high')
                sys.exit(1)                           

    return True

def check_units(config, plants):
    """
    Function that checks the power plant characteristics
    """

    keys = ['Unit', 'Fuel', 'Zone', 'Technology', 'PowerCapacity', 'PartLoadMin', 'RampUpRate', 'RampDownRate',
            'StartUpTime', 'MinUpTime', 'MinDownTime', 'NoLoadCost', 'StartUpCost', 'Efficiency', 'CO2Intensity']
    NonNaNKeys = ['PowerCapacity', 'PartLoadMin', 'RampUpRate', 'RampDownRate', 'Efficiency', 'RampingCost',
                  'StartUpTime', 'CO2Intensity']
    StrKeys = ['Unit', 'Zone', 'Fuel', 'Technology']

    # Special treatment for the Optional key Nunits:
    if 'Nunits' in plants:
        keys.append('Nunits')
        NonNaNKeys.append('Nunits')
        if any([not float(x).is_integer() for x in plants['Nunits']]):
            logging.error('Some values are not integers in the "Nunits" column of the plant database')
            sys.exit(1)
    else:
        logging.info('The columns "Nunits" is not present in the power plant database. A value of one will be assumed by default')

    for key in keys:
        if key not in plants:
            logging.critical('The power plants data does not contain the field "' + key + '", which is mandatory')
            sys.exit(1)

    for key in NonNaNKeys:
        for u in plants.index:
            if type(plants.loc[u, key]) == str:
                logging.critical('A non numeric value was detected in the power plants inputs for parameter "' + key + '"')
                sys.exit(1)
            if np.isnan(plants.loc[u, key]):
                logging.critical('The power plants data is missing for unit number ' + str(u) + ' and parameter "' + key + '"')
                sys.exit(1)

    for key in StrKeys:
        for u in plants.index:
            if not type(plants.loc[u, key]) == str:
                logging.critical(
                    'A numeric value was detected in the power plants inputs for parameter "' + key + '". This column should contain strings only.')
                sys.exit(1)
            elif plants.loc[u, key] == '':
                logging.critical('An empty value was detected in the power plants inputs for unit "' + str(
                    u) + '" and parameter "' + key + '"')
                sys.exit(1)

    lower = {'PowerCapacity': 0, 'PartLoadMin': 0, 'StartUpTime': 0, 'MinUpTime': 0, 'MinDownTime': 0, 'NoLoadCost': 0,
             'StartUpCost': 0}
    lower_hard = {'RampUpRate': 0, 'RampDownRate': 0, 'Efficiency': 0}
    higher = {'PartLoadMin': 1, 'Efficiency': 1}
    higher_time = {'MinUpTime': 0, 'MinDownTime': 0}  # 'StartUpTime':0,

    # Special treatment for the Optional key Nunits:
    if 'Nunits' in plants:
        lower_hard['Nunits'] = 0

    for key in lower:
        if any(plants[key] < lower[key]):
            plantlist = plants[plants[key] < lower[key]]
            plantlist = plantlist['Unit'].tolist()
            logging.critical(
                'The value of ' + key + ' should be higher or equal to zero. A negative value has been found for units ' + str(
                    plantlist))
            sys.exit(1)

    for key in lower_hard:
        if any(plants[key] <= lower_hard[key]):
            plantlist = plants[plants[key] <= lower_hard[key]]
            plantlist = plantlist['Unit'].tolist()
            logging.critical(
                'The value of ' + key + ' should be strictly higher than zero. A null or negative value has been found for units ' + str(
                    plantlist))
            sys.exit(1)

    for key in higher:
        if any(plants[key] > higher[key]):
            plantlist = plants[plants[key] > higher[key]]
            plantlist = plantlist['Unit'].tolist()
            logging.critical(
                'The value of ' + key + ' should be lower or equal to one. A higher value has been found for units ' + str(
                    plantlist))
            sys.exit(1)


    for key in higher_time:
        if any(plants[key] >= config['HorizonLength'] * 24):
            plantlist = plants[plants[key] >= config['HorizonLength'] * 24]
            plantlist = plantlist['Unit'].tolist()
            logging.critical('The value of ' + key + ' should be lower than the horizon length (' + str(
                config['HorizonLength'] * 24) + ' hours). A higher value has been found for units ' + str(plantlist))
            sys.exit(1)
            
    return True


def check_heat_demand(plants,data):
    '''
    Function that checks the validity of the heat demand profiles
    '''
    plants.index = plants['Unit']
    for u in data:
        if u in plants.index:
            plant_CHP_type = plants.loc[u,'CHPType'].lower()
            if pd.isnull(plants.loc[u, 'CHPMaxHeat']):
                plant_Qmax = +np.inf
            else:
                plant_Qmax = plants.loc[u,'CHPMaxHeat']
            if plant_CHP_type == 'extraction':
                Qmin = 0
                Qmax = min(plants.loc[u, 'PowerCapacity'] / plants.loc[u, 'CHPPowerToHeat'], plant_Qmax)
            elif plant_CHP_type == 'back-pressure':
                Qmin = plants.loc[u,'PowerCapacity'] * plants.loc[u,'PartLoadMin'] /plants.loc[u,'CHPPowerToHeat']
                Qmax = min(plants.loc[u, 'PowerCapacity'] / plants.loc[u, 'CHPPowerToHeat'], plant_Qmax)
            elif plant_CHP_type == 'p2h':
                Qmin = 0
                Qmax = plant_Qmax
            else:
                logging.error('The CHP type for unit ' + u + ' is not valid.')
            if np.isnan(Qmax) and plant_CHP_type!='p2h':
                logging.error('CHPPowerToHeat is not defined for unit ' + str(u) + ' appearing in the heat demand profiles')
                sys.exit(1)
            elif data[u].max() > Qmax:
                logging.warn('The maximum thermal demand for unit ' + str(u) + ' (' + str(data[u].max()) + ') is higher than its thermal capacity (' + str(Qmax) + '). Slack heat will be used to cover that.')
            if data[u].min() < Qmin:
                logging.warn('The minimum thermal demand for unit ' + str(u) + ' (' + str(data[u].min()) + ') is lower than its minimum thermal generation (' + str(Qmin) + ' MWth)')
        else:
            logging.warn('The heat demand profile with header "' + str(u) + '" does not correspond to any CHP plant. It will be ignored.')
    return True


def check_df(df, StartDate=None, StopDate=None, name=''):
    """
    Function that check the time series provided as inputs
    """
    ####### Edits #######
    for i in range(0, len(df.columns)):
        df.iloc[:,i] = pd.to_numeric(df.iloc[:,i], errors='ignore')
    
    if isinstance(df.index, pd.DatetimeIndex):
        if not StartDate in df.index:
            if missing != 0:
                if name == 'FuelPrices':
                    if len(key) == 2:
                        logging.warning(
                            'There are ' + str(missing) + ' missing entries in the column ' + key[1] + ' of the zone ' +
                            key[0] + 'of the dataframe ' + name)
                    else:
                        logging.warning('There are ' + str(
                            missing) + ' missing entries in the column ' + key + ' of the dataframe ' + name)
                else:
                    logging.warning('There are ' + str(
                        missing) + ' missing entries in the column ' + key + ' of the dataframe ' + name)
    return True
'''
            logging.warn('The start date ' + str(StartDate) + ' is not in the index of the provided dataframe')
        if not StopDate in df.index:
            logging.warn('The stop date ' + str(StopDate) + ' is not in the index of the provided dataframe')
    if any(np.isnan(df)):
        for key in df:
            missing = np.sum(np.isnan(df[key]))
            # pos = np.where(np.isnan(df.sum(axis=1)))
            # idx_pos = [df.index[i] for i in pos]
            if missing != 0:
                if name=='FuelPrices':
                    logging.warn('There are ' + str(missing) + ' missing entries in the column ' + key[1] + ' of the zone ' + key[0] + 'of the dataframe ' + name)
                else:
                    logging.warn('There are ' + str(missing) + ' missing entries in the column ' + key + ' of the dataframe ' + name)
'''
    ####### Edits #######

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
            vars = SimulationPath  # FIXME: vars not used
        else:
            logging.critical('The argument must a list. Please correct or change the "type" argument')
            sys.exit(1)

    elif store_type == 'pickle':
        if os.path.exists(SimulationPath):
            if os.path.isfile(os.path.join(SimulationPath, 'Inputs.p')):
                vars = cPickle.load(open(os.path.join(SimulationPath, 'Inputs.p'), 'rb'))
                arg_vars = [vars[i]['name'] for i in range(len(vars))]
                for var in list_sets + list_param:
                    if var not in arg_vars:
                        logging.critical('Found Pickle file but does not contain valid DispaSET input (' + var + ' missing)')
                        sys.exit(1)
            else:
                logging.critical('Could not find the Inputs.p file in the specified directory')
                sys.exit(1)
        else:
            logging.critical('The function argument is not a valid directory')
            sys.exit(1)

    elif store_type == 'excel':
        if os.path.exists(SimulationPath):
            if os.path.isfile(os.path.join(SimulationPath, 'InputDispa-SET - Sets.xlsx')):
                a = 1
            else:
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


def isVRE(tech):
    '''
    Function that returns true the technology is a variable renewable energy technology
    '''
    return tech in ['HROR','PHOT','WTON','WTOF']

def isStorage(tech):
    '''
    Function that returns true the technology is a storage technology
    '''
    return tech in ['HDAM','HPHS','CAES','BATS','BEVS','THMS','P2GS']

