# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 12:02:10 2020

@author: Chiara Magni
"""

#Reserve Requirements : (3+5) rule

#dynamic method for the evaluation of reserve requirements based on forecasted 
#load and forecasted available wind and solar power.

import numpy as np
import pandas as pd
import sys,os
sys.path.append(os.path.abspath('../..')) 

WRITE_CSV_FILE = True
filename = 'reserves_3_5_rule'


# %% Inputs

# Folder destinations
input_folder ='../../Database'  # Standard input folder
output_folder ='../../Database/reserves/'  # Standard output folder

# %% import data
# countries 

countries=['AT', 'BE', 'BG', 'CH', 'CY', 'CZ', 'DE', 'DK',
 'EE', 'EL', 'ES', 'FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LT',
 'LU', 'LV', 'MT', 'NL', 'NO', 'PL', 'PT', 'RO', 'SE', 'SI', 'SK',]
countries.remove('CY')
countries.remove('MT')

# import load

inputfile_load = input_folder + '/TotalLoadValue/%s/1h/2016.csv'

load_dict = {}
for c in countries: 
    load_dict[c] = pd.read_csv(inputfile_load %c, header = None, index_col = 0)
    load_dict[c] = load_dict[c].iloc[:,0]

load = pd.DataFrame.from_dict(load_dict)


# import units

inputfile_units = input_folder + '/PowerPlants/%s/clustered.csv'

allunits_dict = {}
for c in countries: 
    allunits_dict[c] = pd.read_csv(inputfile_units %c, index_col = 0)


# import availability

inputfile_availability = input_folder + '/AvailabilityFactors/%s/1h/2016.csv'

availability_dict = {}
for c in countries: 
    availability_dict[c] = pd.read_csv(inputfile_availability %c, index_col = 0)

# %% 1 hour
date_str = '1/1/2016'
start_1h = pd.to_datetime(date_str)
hourly_periods = len(load)
drange_1h = pd.date_range(start_1h, periods=hourly_periods, freq='H')
hour = pd.DataFrame(drange_1h)


# %% reserve calculation

reserve = pd.DataFrame()

for c in countries:    
    
    # 3% of the load forecast
    load_fraction = load[c]*0.03
    load_fraction.reset_index(drop=True, inplace=True)

    
    # 5% of the windon forecast
    if allunits_dict[c]['Technology'].str.contains('WTON').any():
        tmp = allunits_dict[c].index[allunits_dict[c]['Technology']=='WTON']
        wton_cap = allunits_dict[c]['PowerCapacity'][tmp]*allunits_dict[c]['Nunits'][tmp]
        wton_fraction = (availability_dict[c]['WTON']*wton_cap[0])*0.05
        wton_fraction.reset_index(drop=True, inplace=True)
        
    else:
        wton_fraction = pd.Series(0, index=load_fraction.index)
    
    # 5% of the windoff forecast
        
    if allunits_dict[c]['Technology'].str.contains('WTOF').any():
        tmp = allunits_dict[c].index[allunits_dict[c]['Technology']=='WTOF']
        wtof_cap = allunits_dict[c]['PowerCapacity'][tmp]*allunits_dict[c]['Nunits'][tmp]
        wtof_fraction = (availability_dict[c]['WTOF']*wtof_cap[0])*0.05
        wtof_fraction.reset_index(drop=True, inplace=True)

    else:
        wtof_fraction = pd.Series(0, index=load_fraction.index)
    
    # 5% of solar forecast
    
    if allunits_dict[c]['Technology'].str.contains('PHOT').any():
        tmp = allunits_dict[c].index[allunits_dict[c]['Technology']=='PHOT']
        phot_cap = allunits_dict[c]['PowerCapacity'][tmp]*allunits_dict[c]['Nunits'][tmp]
        phot_fraction = (availability_dict[c]['PHOT']*phot_cap[0])*0.05
        phot_fraction.reset_index(drop=True, inplace=True)

    else:
        phot_fraction = pd.Series(0, index=load_fraction.index)
    
    # total reservedemand
    
    rr = load_fraction + wton_fraction + wtof_fraction + phot_fraction
    rr = pd.Series.to_frame(rr)
    rr = pd.DataFrame.rename(rr,columns={0: c})
    reserve = pd.concat([reserve,rr], axis=1)

reserve = pd.concat([hour,reserve],axis=1)
reserve = pd.DataFrame.set_index(reserve,0)


# %% Write csv file:
def write_csv_files(filename, reserve, write_csv=None):
    """
    Function that generates .csv files in the Output/Database/PowerPlants/ folder
    :filename:      clustered for example (has to be a string)
    :units:                     allunits for example
    """
    filename = filename + '.csv'
    if write_csv is True:
        reserve.to_csv(output_folder + filename)
    else:
        print('[WARNING ]: ' + 'WRITE_CSV_FILES = False, unable to write .csv files')


write_csv_files(filename, reserve, WRITE_CSV_FILE)


