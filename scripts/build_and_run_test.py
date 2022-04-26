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
import os
import sys

sys.path.append(os.path.abspath('..'))

# Import Dispa-SET
import dispaset as ds

# Load the configuration file
config = ds.load_config('../ConfigFiles/ConfigTest.xlsx')

# Limit the simulation period (for testing purposes, comment the line to run the whole year)
# config['StartDate'] = (2016, 1, 1, 0, 0, 0)
# config['StopDate'] = (2016, 1, 7, 0, 0, 0)

# Build the simulation environment:
SimData = ds.build_simulation(config)

# Solve using GAMS:
_ = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])

# Load the simulation results:
inputs, results = ds.get_sim_results(config['SimulationDirectory'], cache=False)

# import pandas as pd
import pandas as pd

rng = pd.date_range('2016-1-01', '2016-1-07', freq='H')
# Generate country-specific plots
ds.plot_zone(inputs, results, z='Z1', rng=rng)

# Generate country-specific plots
ds.plot_zone(inputs, results)

# Bar plot with the installed capacities in all countries:
cap = ds.plot_zone_capacities(inputs, results)

# Bar plot with installed storage capacity
sto = ds.plot_tech_cap(inputs)

# Violin plot for CO2 emissions
ds.plot_co2(inputs, results, figsize=(9, 6), width=0.9)

# Bar plot with the energy balances in all countries:
ds.plot_energy_zone_fuel(inputs, results, ds.get_indicators_powerplant(inputs, results))

# Analyse the results for each country and provide quantitative indicators:
r = ds.get_result_analysis(inputs, results)

# Analyze power flow tracing
pft, pft_prct = ds.plot_power_flow_tracing_matrix(inputs, results, cmap="magma_r", figsize=(15, 10))

# Plot net flows on a map
ds.plot_net_flows_map(inputs, results, terrain=True, margin=3, bublesize=5000, figsize=(8, 7))

# Plot congestion in the interconnection lines on a map
ds.plot_line_congestion_map(inputs, results, terrain=True, margin=3, figsize=(9, 7), edge_width=3.5, bublesize=100)
