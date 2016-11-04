# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 15:30:51 2016

@author: Sylvain Quoilin, JRC
"""

from __future__ import division
import DispaSET as ds
import argparse
import sys
import os 
import numpy as np
import pandas as pd

#%%

parser = argparse.ArgumentParser(description='Reads the dispa-set results from the simulation environment folder')
parser.add_argument("-r", "--read", help="Path to the simulation directory to be read")
args = parser.parse_args()

if args.read:
    path = args.read
    print "Using directory " + path + " to read the simulation results"  
else:
    print "No folder specified to build the simulation environment"
    try:
        import easygui
        path = easygui.diropenbox('Please select the config file')
        print "Using directory " + path + " to read the simulation results"
    except:
        path = ds.InputDir('Specify the path to the simulation directory (e.g. ../Simulation/): ')
        print "Using directory " + path + " to read the simulation results"
    build = True
    simulate = False

if not os.path.isdir(path):
    sys.exit('Could not find directory ' + path + ' \n Exiting')
    
    
if not ds.gdxcc_ok:
        print 'WARNING: the gdxcc library could not be loaded.'   

inputs,results = ds.GetResults(path=path,cache=True)

# Temporary fix (to be removed):
inputs['units'].index = ds.shrink_to_64(ds.clean_strings(inputs['units']['Unit'].tolist()))


#%%

# Formatting the input data: 
datain = ds.ds_to_df(inputs)

# Select the time period for plotting
rng = pd.DatetimeIndex(start='2015-01-2 00:00:00',end='2015-12-31 04:01:00',freq='h')

# plotting the detailed analysis for on of the zones, randomly
Nzones = len(inputs['sets']['n'])
c = inputs['sets']['n'][np.random.randint(Nzones)]

ds.CountryPlots(inputs,results,c,rng=rng)

r = ds.ResultAnalysis(inputs,results)

if 'OutputStorageLevel' in results:
    results['OutputStorageLevel'].plot(figsize=(15,10))
    

#%%

# bar plots

PPindicators = ds.PerPowerPlantIndicators(inputs,results)
ax = ds.EnergyBarPlot(datain,results,PPindicators)

#ax.barh(0.7,left=ax.get_xticks()-0.4, width=0.8,height=0.1,linewidth=2,color='k')


'''
bar plots with max demand, installed capacity in each country
line capacities wrong for some countries
'''
