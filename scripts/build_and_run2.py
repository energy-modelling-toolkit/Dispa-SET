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
config = ds.load_config_excel('../ConfigFiles/ConfigEU.xlsx')
config['StartDate'] = (2015, 1, 1, 0, 0, 0)
config['StopDate'] = (2015, 1, 7, 0, 0, 0)
#config['modifiers']['demand'] = 0.15
#config['default']['ShareOfFlexibleDemand'] = 0.0
#config['default']['WaterValue'] = 1000

# Build the simulation environment:
SimData = ds.build_simulation(config)

# Solve using GAMS:
r = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])


# Load the inputs and the results of the simulation
inputs,results = ds.get_sim_results(path=config['SimulationDirectory'],cache=False)

# if needed, define the plotting range for the dispatch plot:
import pandas as pd
rng = pd.date_range(start='2016-01-01',end='2016-12-31',freq='h')

# Generate country-specific plots
ds.plot_zone(inputs,results,rng=rng)

# Bar plot with the installed capacities in all countries:
cap = ds.plot_zone_capacities(inputs)

# Bar plot with the energy balances in all countries:
ds.plot_energy_zone_fuel(inputs,results,ds.get_indicators_powerplant(inputs,results))

# Analyse the results for each country and provide quantitative indicators:
r = ds.get_result_analysis(inputs,results)