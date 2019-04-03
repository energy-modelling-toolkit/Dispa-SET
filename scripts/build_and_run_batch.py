# -*- coding: utf-8 -*-
"""
Eexample file showing how to access the Dispa-SET api to read a configuration file, 
and build batch runs in GAMS or PYOMO

@author: K Kavvadias
"""

from itertools import product
import cPickle as pickle
import logging

import os
os.chdir('..')

# Import Dispa-SET
import dispaset as ds

# Define path to save results. Each run will be a pickled tuple: (inputs, results)
path_to_save = r'./scenario_runs'

# Load the configuration file
config = ds.load_config_excel('ConfigFiles/ConfigCY.xlsx')

# Define your different input files as a list. The number of total runs will be the product of the length of the defined lists.
heat_demand_scen = ['./Database/Heat_demand/CY/Vassilikos_CCP2.csv',
                    './Database/Heat_demand/CY/Vassilikos_CCP2.csv_80']

power_scen = ['./Database/PowerPlants/##/2015.csv',
              './Database/PowerPlants/##/2015_heat.csv',
              ]

cost_heat_slack_scen = [20.0, 51.0, 100.0]
res_mod_scen = [0.8,1.0,1.2]

all_combinations = list(product(heat_demand_scen, power_scen, cost_heat_slack_scen, res_mod_scen))
logging.info(str(len(all_combinations)) + ' runs will be made')  # 36 for this case

i = 0

for iheat, ipower, icost_heat, iresmod in all_combinations:
    i += 1
    scenario_name = "run_1_" + str(i)

    # Load run specific inputs
    config['HeatDemand'] = iheat
    config['PowerPlantData'] = ipower
    config['default']['CostHeatSlack'] = icost_heat
    config['modifiers']['Solar'] = config['modifiers']['wind'] = iresmod

    # Build simulation
    SimData = ds.build_simulation(config)
    # Solve using GAMS:
    r = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])

    # Load inputs and results to memory
    inputs, results = ds.get_sim_results(path=config['SimulationDirectory'], cache=True)
    # Save inputs and results to a file
    filename = os.path.join(path_to_save, scenario_name+".p")
    with open(filename, "wb") as f:
        pickle.dump((inputs, results), f, protocol=pickle.HIGHEST_PROTOCOL)
    logging.info('Saved {} to pickle file'.format(scenario_name))
