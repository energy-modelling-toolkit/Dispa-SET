# -*- coding: utf-8 -*-
"""
Minimalist example file showing how to access the Dispa-SET api to read a configuration file, 
create a simulation environment folder and run the simulation in GAMS

@author: Sylvain Quoilin
"""

# Add the root folder of Dispa-SET to the path so that the library can be loaded:
import sys,os
sys.path.append(os.path.abspath('..'))

# Import Dispa-SET
import dispaset as ds

# Load the configuration file
config = ds.load_config('../ConfigFiles/ConfigBE.xlsx')

config['Heat']['RunHeat'] = 1
config['Heat']['NumberHP']=1000000
config['Heat']['NumberWH']=1000000
config['StartDate']=(2015,1,1,0,0,0)
config['StopDate']=(2015,12,31,23,59,59)
config['modifiers']['Demand']=1
config['modifiers']['Solar']=2
config['modifiers']['Wind']=2
# config['Heat']['NumberSimulatedHP']=1
# config['Heat']['NumberSimulatedWH']=1


# Build the simulation environment:
SimData = ds.build_simulation(config)

# Solve using GAMS:
#r = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])