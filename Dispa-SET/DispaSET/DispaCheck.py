'''
This files gathers different functions used in the DispaSET to check the input
data

__author__ = 'Sylvain Quoilin (sylvain.quoilin@ec.europa.eu)'
'''

import os
import sys
import numpy as np
import pandas as pd 
from warnings import warn

def check_units(config,plants):
    '''
    Function that checks the power plant characteristics
    '''
    
    keys = ['Unit','Fuel','Zone','Technology','PowerCapacity','PartLoadMin','RampUpRate','RampDownRate','StartUpTime','MinUpTime','MinDownTime','NoLoadCost','StartUpCost','Efficiency','CO2Intensity']
    NonNaNKeys = ['PowerCapacity','PartLoadMin','RampUpRate','RampDownRate','Efficiency','RampingCost','StartUpTime','CO2Intensity']
    StrKeys = ['Unit','Zone','Fuel','Technology',]    
    
    for key in keys:
        if key not in plants:
            sys.exit('The power plants data does not contain the field "' + key + '", which is mandatory')

    for key in NonNaNKeys:
        for u in plants.index:
            if type(plants.loc[u,key]) == str:
                sys.exit('A non numeric value was detected in the power plants inputs for parameter "' + key + '"')
            if np.isnan(plants.loc[u,key]):
                sys.exit('The power plants data is missing for unit number ' + str(u) + ' and parameter "' + key + '"')

    for key in StrKeys:
        for u in plants.index:
            if not type(plants.loc[u,key]) == str:
                sys.exit('A numeric value was detected in the power plants inputs for parameter "' + key + '". This column should contain strings only.')
            elif plants.loc[u,key] == '':
                sys.exit('An empty value was detected in the power plants inputs for unit "' + str(u) + '" and parameter "' + key + '"')

    lower = {'PowerCapacity':0,'PartLoadMin':0,'StartUpTime':0,'MinUpTime':0,'MinDownTime':0,'NoLoadCost':0,'StartUpCost':0}
    lower_hard = {'RampUpRate':0,'RampDownRate':0,'Efficiency':0}
    higher = {'PartLoadMin':1,'Efficiency':1}
    higher_time = {'MinUpTime':0,'MinDownTime':0}     #'StartUpTime':0,
    
    for key in lower:
        if any(plants[key] < lower[key]):
            plantlist = plants[plants[key] < lower[key]]
            plantlist = plantlist['Unit'].tolist()
            sys.exit('The value of ' + key + ' should be higher or equal to zero. A negative value has been found for units ' + str(plantlist))
 
    for key in lower_hard:
        if any(plants[key] <= lower_hard[key]):
            plantlist = plants[plants[key] <= lower_hard[key]]
            plantlist = plantlist['Unit'].tolist()
            sys.exit('The value of ' + key + ' should be strictly higher than zero. A null or negative value has been found for units ' + str(plantlist))

    for key in higher:
        if any(plants[key] > higher[key]):
            plantlist = plants[plants[key] > higher[key]]
            plantlist = plantlist['Unit'].tolist()
            sys.exit('The value of ' + key + ' should be lower or equal to one. A higher value has been found for units ' + str(plantlist))

    for key in higher_time:
        if any(plants[key] >= config['HorizonLength']*24):
            plantlist = plants[plants[key] >= config['HorizonLength']*24]
            plantlist = plantlist['Unit'].tolist()
            sys.exit('The value of ' + key + ' should be lower than the horizon length (' + str(config['HorizonLength']*24) + ' hours). A higher value has been found for units ' + str(plantlist))




    return True

def check_df(df,StartDate=None,StopDate=None,name=''):
    '''
    Function that check the time series provided as inputs
    '''
    
    if type(df.index) is pd.tseries.index.DatetimeIndex:
        if not StartDate in df.index:
            warn('The start date ' + str(StartDate) + ' is not in the index of the provided dataframe')
        if not StopDate in df.index:
            warn('The stop date ' + str(StopDate) + ' is not in the index of the provided dataframe')
    if any(np.isnan(df)):
        for key in df:
            missing =  np.sum(np.isnan(df[key]))
            #pos = np.where(np.isnan(df.sum(axis=1)))
            #idx_pos = [df.index[i] for i in pos]
            if missing != 0:
                print('WARNING: There are ' + str(missing) + ' missing entries in the column ' + key + ' of the dataframe ' + name)
    return True
    
def check_simulation_environment(SimulationPath,type='pickle',firstline=7):
    '''
    Function to test the validity of disapset inputs
    :param SimulationPath:          Path to the simulation folder
    :param type:                    choose between: "list", "excel", "pickle"
    :param firstline:               Number of the first line in the data (only if type=='excel')
    '''
    
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
    
    if type == 'list':
         if isinstance(SimulationPath,list):
             # The list of sets and parameters has been passed directly to the function, checking that all are present:
             SimulationPath_vars = [SimulationPath[i]['name'] for i in range(len(SimulationPath))]
             for var in list_sets + list_param:
                 if var not in SimulationPath_vars:
                     sys.exit('The variable "' + var + '" has not been found in the list of input variables')
             vars=SimulationPath
         else:
             sys.exsit('The argument must a list. Please correct or change the "type" argument')
            
    elif type == 'pickle':
        if os.path.exists(SimulationPath):
            if os.path.isfile(os.path.join(SimulationPath,'Inputs.p')):
                vars = cPickle.load(open(os.path.join(SimulationPath,'Inputs.p'), 'rb'))
                arg_vars = [vars[i]['name'] for i in range(len(vars))]
                for var in list_sets + list_param:
                     if var not in arg_vars:
                         sys.exit('Found Pickle file but does not contain valid DispaSET input (' + var + ' missing)')
            else:
                sys.exit('Could not find the Inputs.p file in the specified directory')
        else:
            sys.exit('The function argument is not a valid directory')
    
    elif type == 'excel':
        if os.path.exists(SimulationPath):
            if os.path.isfile(os.path.join(SimulationPath,'InputDispa-SET - Sets.xlsx')):
                a = 1
            else:
                sys.exit("Could not find the file 'InputDispa-SET - Sets.xlsx'")
            for var in list_param:
                if os.path.isfile(os.path.join(SimulationPath,'InputDispa-SET - ' + var +'.xlsx')):
                    a = 1
                else:
                    sys.exit("Could not find the file 'InputDispa-SET - " + var +".xlsx'")
                
        else:
            sys.exit('The function argument is not a valid directory')
    
    else:
        sys.exit('The "type" parameter must be one of the following : "list", "excel", "pickle"')
    
