# -*- coding: utf-8 -*-
"""
This script runs the Dispa-SET EU model with the 2016 data. The main steps are:
    - Load Dispa-SET
    - Load the config file for the EU model
    - build the mode
    - run the model
    - display and analyse the results

@author: Sylvain Quoilin
"""

# Add the root folder of Dispa-SET to the path so that the library can be loaded:
import sys,os
sys.path.append(os.path.abspath('..'))

# Import Dispa-SET
import dispaset as ds

# Load the configuration file
config = ds.load_config_excel('../ConfigFiles/ConfigEU.xlsx')

# Limit the simulation period (for testing purposes, comment the line to run the whole year)
config['StartDate'] = (2016, 1, 1, 0, 0, 0)
config['StopDate'] = (2016, 1, 7, 0, 0, 0)

# Build the simulation environment:
SimData = ds.build_simulation(config)

# Solve using GAMS:
_ = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])

# Load the simulation results:
inputs,results = ds.get_sim_results(config['SimulationDirectory'],cache=False)

# Generate country-specific plots
ds.plot_country(inputs,results)

# Bar plot with the installed capacities in all countries:
cap = ds.plot_country_capacities(inputs)

# Bar plot with the energy balances in all countries:
ds.plot_energy_country_fuel(inputs,results,ds.get_indicators_powerplant(inputs,results))

# Analyse the results for each country and provide quantitative indicators:
r = ds.get_result_analysis(inputs,results)

# Plot the reservoir levels
ds.storage_levels(inputs,results)