# -*- coding: utf-8 -*-
"""
Minimalist example file showing how to access the Dispa-SET api to read a configuration file, 
create a simulation environment folder and run the simulation in GAMS or PYOMO

@author: Sylvain Quoilin
"""

# Add the root folder of Dispa-SET to the path so that the library can be loaded:
import sys,os
sys.path.append(os.path.abspath('..'))

# Import Dispa-SET
import dispaset as ds

# Load the configuration file
config = ds.load_config_excel('../ConfigFiles/ConfigTest.xlsx')

# Build the simulation environment:
SimData = ds.build_simulation(config)

# Solve using PYOMO/GAMS:
#r = ds.solve_pyomo(config['SimulationDirectory'])
r = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])