# -*- coding: utf-8 -*-
"""
Minimalist example file showing how to read the results of a Dispa-SET run

@author: Sylvain Quoilin
"""

# Add the root folder of Dispa-SET to the path so that the library can be loaded:
import sys,os
sys.path.append(os.path.abspath('..'))

# Import Dispa-SET
import dispaset as ds
import pandas as pd

# Load the inputs and the results of the simulation
inputs,results = ds.get_sim_results(path='../Simulations/simulation_test',cache=False)

# Get the actual date range from the results
if 'OutputPower' in results:
    rng = results['OutputPower'].index
else:
    # Fallback to a default range if PowerDemand is not available
    rng = pd.date_range(start='2016-01-01',end='2016-12-31',freq='h')

# Generate country-specific plots
ds.plot_zone(inputs,results,rng=rng)

# Bar plot with the installed capacities in all countries:
cap = ds.plot_zone_capacities(inputs,results)

# Bar plot with the energy balances in all countries:
ds.plot_energy_zone_fuel(inputs,results,ds.get_indicators_powerplant(inputs,results))

# Analyse the results for each country and provide quantitative indicators:
r = ds.get_result_analysis(inputs,results)

# Plot the reservoir levels
#ds.storage_levels(inputs,results)
#ds.plot_storage_levels(inputs,results,'NO')

ds.plot_power_flow_tracing_matrix(inputs, results)

ds.plot_net_flows_map(inputs,results)
ds.plot_line_congestion_map(inputs,results)
