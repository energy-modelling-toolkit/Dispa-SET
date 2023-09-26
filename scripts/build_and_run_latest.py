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
import pandas as pd
# Import Dispa-SET
import dispaset as ds

sys.path.append(os.path.abspath('..'))

# Load the configuration file
config = ds.load_config('../ConfigFiles/ConfigTestBoundarySector.xlsx')

# Limit the simulation period (for testing purposes, comment the line to run the whole year)
config['StartDate'] = (2015, 1, 1, 0, 0, 0)
config['StopDate'] = (2015, 1, 15, 0, 0, 0)

# Build the simulation environment:
SimData = ds.build_simulation(config)

# Solve using GAMS:
_ = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])

# Load the simulation results:
inputs, results = ds.get_sim_results(config['SimulationDirectory'], cache=False)
inputs_MTS, results_MTS = ds.get_sim_results(config['SimulationDirectory'], cache=False, inputs_file='Inputs_MTS.p',
                                             results_file='Results_MTS.gdx')

colors = {'LIG': '#af4b9180', 'PEA': '#af4b9199', 'HRD': 'darkviolet', 'OIL': 'magenta', 'GAS': '#d7642dff',
          'NUC': '#466eb4ff', 'SUN': '#e6a532ff', 'WIN': '#41afaaff', 'WAT': '#00a0e1ff', 'HYD': '#A0522D',
          'BIO': '#7daf4bff', 'AMO': '#ffff00ff', 'GEO': '#7daf4bbf', 'Storage': '#b93c46ff', 'FlowIn': 'red',
          'FlowOut': 'green', 'OTH': '#57D53B', 'WST': '#b9c337ff', 'HDAM': '#00a0e1ff', 'HDAMC': '#00a0e1ff',
          'HPHS': 'blue', 'THMS': '#C04000ff', 'BATS': '#41A317ff', 'BEVS': '#b9c33799', 'SCSP': '#e6a532ff',
          'P2GS': '#A0522D', 'ShedLoad': '#ffffffff', 'AIR': '#aed6f1ff', 'WHT': '#a93226ff', 'ELE': '#2C75FFff',
          'THE': '#c70509ff',
          'Z2_h2': 'cyan', 'Z2_w2': 'magenta', 'curtailment': 'red', 'Z1_h2': 'cyan',
          'Z1_th': 'yellow'}

rng = pd.date_range('2015-1-01', '2015-1-07', freq='H')
# Generate country-specific plots
ds.plot_zone(inputs, results, z='Z2', rng=rng, colors=colors)

# Generate country-specific plots
ds.plot_zone(inputs, results)

# Bar plot with the installed capacities in all countries:
cap = ds.plot_zone_capacities(inputs, results)

# Bar plot with installed storage capacity
sto = ds.plot_tech_cap(inputs)

# # Violin plot for CO2 emissions
# ds.plot_co2(inputs, results, figsize=(9, 6), width=0.9)
#
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
