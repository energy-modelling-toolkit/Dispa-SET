# -*- coding: utf-8 -*-
"""
Minimalist example file showing how to access the Dispa-SET api to read a configuration file, 
create a simulation environment folder and run the simulation in GAMS or PYOMO

The script can be run in ipython, but for some reasons, it does not work more than once with the gams API.

@author: Sylvain Quoilin
"""

# Change directory to the root folder of Dispa-SET:
import sys,os
sys.path.append(os.path.abspath('..'))

# Import Dispa-SET
import dispaset as ds

# Load the configuration file
config = ds.load_config('../ConfigFiles/Config_Mid_Term_Scheduling.xlsx')

## Build simulation data with new profiles
SimData = ds.build_simulation(config,mts_plot=True)

#Solve using GAMS:
r = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])

# Load the simulation results:
inputs,results = ds.get_sim_results(config['SimulationDirectory'],cache=False)

# Generate country-specific plots
ds.plot_zone(inputs,results,z='DE')

# Bar plot with the installed capacities in all countries:
cap = ds.plot_zone_capacities(inputs)

# Bar plot with the energy balances in all countries:
ds.plot_energy_zone_fuel(inputs,results,ds.get_indicators_powerplant(inputs,results))

# Analyse the results for each country and provide quantitative indicators:
r = ds.get_result_analysis(inputs,results)