# -*- coding: utf-8 -*-
"""
Minimalist example file showing how to access the Dispa-SET api to read a configuration file, 
create a simulation environment folder and run the simulation in GAMS or PYOMO

The script can be run in ipython, but for some reasons, it does not work more than once with the gams API.

@author: Sylvain Quoilin
"""

# Change directory to the root folder of Dispa-SET:
import os
os.chdir('..')

# Import Dispa-SET
import DispaSET as ds

# Load the configuration file
config = ds.load_config_excel('ConfigFiles/ConfigEU.xlsx')

# Build the simulation environment:
SimData = ds.build_simulation(config)

# Solve using PYOMO/GAMS:
#r = ds.solve_pyomo(config['SimulationDirectory'])
#r = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])