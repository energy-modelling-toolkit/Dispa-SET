# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 22:53:40 2015

@author: sylvain


Function to test the validity of disapset inputs
"""
import sys
import os
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
 
arg = 'simulation'
type = 'pickle'
firstline = 7


if type == 'list':
     if isinstance(arg,list):
         # The list of sets and parameters has been passed directly to the function, checking that all are present:
         arg_vars = [arg[i]['name'] for i in range(len(arg))]
         for var in list_sets + list_param:
             if var not in arg_vars:
                 sys.exit('The variable "' + var + '" has not been found in the list of input variables')
         vars=arg
     else:
         sys.exsit('The argument must a list. Please correct or change the "type" argument')
        
elif type == 'pickle':
    if os.path.exists(arg):
        if os.path.isfile(os.path.join(arg,'Inputs.p')):
            vars = cPickle.load(open(os.path.join(arg,'Inputs.p'), 'rb'))
            arg_vars = [vars[i]['name'] for i in range(len(vars))]
            for var in list_sets + list_param:
                 if var not in arg_vars:
                     sys.exit('Found Pickle file but does not contain valid DispaSET input (' + var + ' missing)')
        else:
            sys.exit('Could not find the Inputs.p file in the specified directory')
    else:
        sys.exit('The function argument is not a valid directory')

elif type == 'excel':
    if os.path.exists(arg):
        if os.path.isfile(os.path.join(arg,'InputDispa-SET - Sets.xlsx')):
            a = 1
        else:
            sys.exit("Could not find the file 'InputDispa-SET - Sets.xlsx'")
        for var in list_param:
            if os.path.isfile(os.path.join(arg,'InputDispa-SET - ' + var +'.xlsx')):
                a = 1
            else:
                sys.exit("Could not find the file 'InputDispa-SET - " + var +".xlsx'")
            
    else:
        sys.exit('The function argument is not a valid directory')

else:
    sys.exit('The "type" parameter must be one of the following : "list", "excel", "pickle"')

