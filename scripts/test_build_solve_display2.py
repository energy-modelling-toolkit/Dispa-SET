# -*- coding: utf-8 -*-
"""
Minimalist example file showing how to access the Dispa-SET api to
-  read a configuration file, 
- create a simulation environment folder and 
- run the simulation using either GAMS or Linopy
- load the results and display them

This script must be run from the root folder of Dispa-SET.

@author: Sylvain Quoilin
"""

# Add the root folder of Dispa-SET to the path so that the library can be loaded:
import sys,os
import pandas as pd
sys.path.append(os.path.abspath('..'))

# Import Dispa-SET
import dispaset as ds

# Load the configuration file
config = ds.load_config('ConfigFiles/ConfigTest.yml')

# Build the simulation environment:
SimData = ds.build_simulation(config)

r = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])


# Load the inputs and the results of the simulation
inputs,results = ds.get_sim_results(path=config['SimulationDirectory'],cache=False)

# if needed, define the plotting range for the dispatch plot:
rng = pd.date_range(start='2016-01-01',end='2016-12-31',freq='h')

# Generate country-specific plots
ds.plot_zone(inputs,results)

# Test the new boundary sector dispatch plot
ds.plot_dispatchX(inputs, results, z='Z1_h2')

# Bar plot with the installed capacities in all countries:
cap = ds.plot_zone_capacities(inputs,results)

# Bar plot with the energy balances in all countries:
ds.plot_energy_zone_fuel(inputs,results,ds.get_indicators_powerplant(inputs,results))

# Analyse the results for each country and provide quantitative indicators:
r = ds.get_result_analysis(inputs,results)

# Plot the reservoir levels
#ds.storage_levels(inputs,results)
#ds.plot_storage_levels(inputs,results,'NO')

#ds.plot_power_flow_tracing_matrix(inputs, results)

