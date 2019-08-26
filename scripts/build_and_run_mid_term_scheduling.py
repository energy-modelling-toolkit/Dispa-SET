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
import dispaset as ds

# Load the configuration file
config = ds.load_config_excel('ConfigFiles/Config_Mid_Term_Scheduling.xlsx')

# locate simmulation folder with temporary solutions
pathname = config['SimulationDirectory']
path, file = os.path.split(pathname)
head, tail = os.path.split(os.path.split(pathname)[0])
path = tail + '/' + file 

# List of countries where new reservoir levels should be calculated eg. ['AT','BE',...'UK']
countries = ['AT', 'CH']

# Calculate and plot new profiles
new_profiles = ds.mid_term_scheduling(config, countries, path)
new_profiles.plot()

# Build simulation data with new profiles
SimData = ds.build_simulation(config, new_profiles)

##Solve using GAMS:
r = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])

