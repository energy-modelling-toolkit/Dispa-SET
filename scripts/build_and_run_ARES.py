# -*- coding: utf-8 -*-
"""
Minimalist example file showing how to access the Dispa-SET api to read a configuration file, 
create a simulation environment folder and run the simulation in GAMS

The script can be run in ipython, but for some reasons, it does not work more than once with the gams API.

@author: Sylvain Quoilin
"""

# Change directory to the root folder of Dispa-SET:
import sys, os

sys.path.append(os.path.abspath('..'))

# Import Dispa-SET
import dispaset as ds
import pickle

# Load the configuration file
config = ds.load_config('../ConfigFiles/Config_ARES_ALLPP.xlsx')

scenarios = ['Baseline_', 'Baseline_NTC_', 'TEMBA_Reference_', 'TEMBA_1.5deg_', 'TEMBA_2.0deg_']

scenario = 'Baseline_'

# define scenario name
def run_scenarios(scenario):
    config['SimulationDirectory'] = config['SimulationDirectory'][:-8] + scenario + str('2015')
    if scenario == 'Baseline_':
        config['NTC'] = config['NTC'][:-8] + '2018.csv'
    elif (scenario == 'Baseline_NTC_') or (scenario == 'TEMBA_Reference_') or (scenario == 'TEMBA_1.5deg_') or \
            (scenario == 'TEMBA_2.0deg_'):
        config['NTC'] = config['NTC'][:-8] + '2025.csv'

    if (scenario == 'Baseline_') or (scenario == 'Baseline_NTC_'):
        for year in range(1980, 1982, 1):
            # adjust simmulation year
            start_date = [year, 1, 1, 0, 0, 0]
            stop_date = [year, 12, 31, 23, 59, 0]
            config['StartDate'] = tuple(start_date)
            config['StopDate'] = tuple(stop_date)
            # rename simulation directory
            config['SimulationDirectory'] = config['SimulationDirectory'][:-4] + str(year)
            # load different weather years
            config['RenewablesAF'] = config['RenewablesAF'][:-8] + str(year) + '.csv'
            config['ReservoirScaledInflows'] = config['ReservoirScaledInflows'][:-8] + str(year) + '.csv'

            # Build simulation data with new profiles
            SimData = ds.build_simulation(config, mts_plot=True, MTSTimeStep=24)
            # Solve using GAMS:
            r = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])
    else:
        config['Demand'] = config['Demand'][:-12] + 'TEMBA_' + str(2015) + '.csv'
        config['Outages'] = config['Outages'][:-12] + scenario + str(2015) + '.csv'
        config['PowerPlantData'] = config['PowerPlantData'][:-12] + scenario + str(2015) + '.csv'
        for year in range(2035, 2046, 10):
            # adjust simmulation year
            start_date = [year, 1, 1, 0, 0, 0]
            stop_date = [year, 12, 31, 23, 59, 0]
            config['StartDate'] = tuple(start_date)
            config['StopDate'] = tuple(stop_date)
            # adjust simulation directory
            config['SimulationDirectory'] = config['SimulationDirectory'][:-4] + str(year)
            # use appropriate NTC's
            if (year >= 2030) and (year < 2040):
                config['NTC'] = config['NTC'][:-8] + '2030.csv'
            if (year >= 2040) and (year < 2050):
                config['NTC'] = config['NTC'][:-8] + '2040.csv'
            if year >= 2050:
                config['NTC'] = config['NTC'][:-8] + '2050.csv'
            # Remove transmission price
            config['default']['PriceTransmission'] = 0.0
            # use typical weather year
            config['RenewablesAF'] = config['RenewablesAF'][:-8] + str(2005) + '.csv'
            config['ReservoirScaledInflows'] = config['ReservoirScaledInflows'][:-8] + str(2001) + '.csv'
            # choose TEMBA demand projection for selected year:
            config['Demand'] = config['Demand'][:-8] + str(year) + '.csv'
            # chose TEMBA Outage factors:
            config['Outages'] = config['Outages'][:-8] + str(year) + '.csv'
            # chose power plants
            config['PowerPlantData'] = config['PowerPlantData'][:-8] + str(year) + '.csv'
            # adjust carbon price:
            if scenario == 'TEMBA_2.0deg_':
                config['default']['PriceOfCO2'] = 0.0043948452413487*year**3 - 26.4853074773464*year**2 + \
                                                  53203.1614853172*year - 35623844.7886163
            if scenario == 'TEMBA_1.5deg_':
                config['default']['PriceOfCO2'] = 0.0309047657613287*year**3 - 187.147451332543*year**2 + \
                                                  377765.2901937630*year - 254179108.818786
            if config['default']['PriceOfCO2'] < 0:
                config['default']['PriceOfCO2'] = 0.0
            config['default']['WaterValue'] = 401 + config['default']['PriceOfCO2']
            config['default']['CostLoadShedding'] = 400 + config['default']['PriceOfCO2']

            # Build simulation data with new profiles
            SimData = ds.build_simulation(config, mts_plot=True, MTSTimeStep=24)
            # Solve using GAMS:
            r = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])

run_scenarios(scenario)

# # Load the simulation results:
# def load_all_results(scenarios):
#     all_inputs = {}
#     all_results = {}
#     all_r = {}
#     all_c = {}
#     all_operation = {}
#     old_path = config['SimulationDirectory']
#     for scenario in scenarios:
#         config['SimulationDirectory'] = old_path[:-8] + scenario + str('2015')
#         for year in range(2025, 2046, 10):
#             config['SimulationDirectory'] = config['SimulationDirectory'][:-4] + str(year)
#             inputs, results = ds.get_sim_results(config['SimulationDirectory'],cache=False)
#             all_inputs[scenario + str(year)], all_results[scenario + str(year)] = inputs, results
#             all_r[scenario + str(year)] = ds.get_result_analysis(inputs, results)
#             all_c[scenario + str(year)] = ds.CostExPost(inputs, results)
#             all_operation[scenario + str(year)] = ds.get_units_operation_cost(inputs, results)
#     with open('E:\Projects\Github\DispaSET-SideTools\Inputs\ARES_Africa/TEMBA/' + 'TEMBA_results' + '.p', 'wb') as handle:
#         pickle.dump(all_inputs, handle, protocol=pickle.HIGHEST_PROTOCOL)
#         pickle.dump(all_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
#         pickle.dump(all_r, handle, protocol=pickle.HIGHEST_PROTOCOL)
#         pickle.dump(all_c, handle, protocol=pickle.HIGHEST_PROTOCOL)
#         pickle.dump(all_operation, handle, protocol=pickle.HIGHEST_PROTOCOL)
#     return
#
# # load_all_results(['TEMBA_Reference_', 'TEMBA_2.0deg_', 'TEMBA_1.5deg_'])
#
# config['SimulationDirectory'] = config['SimulationDirectory'][:-8] + scenario + str(2035)
# inputs, results = ds.get_sim_results(config['SimulationDirectory'],cache=False)
# import pandas as pd
# rng = pd.date_range('2035-1-1','2035-1-15',freq='H')
# # Generate country-specific plots
# ds.plot_zone(inputs,results,z='DZ',rng=rng)
# # ds.plot_zone(inputs,results,z='TZ')
#
# # Bar plot with the installed capacities in all countries:
# cap = ds.plot_zone_capacities(inputs)
#
# # # Bar plot with the energy balances in all countries:
# aa = ds.plot_energy_zone_fuel(inputs,results,ds.get_indicators_powerplant(inputs,results))
#
#
# #
# # Analyse the results for each country and provide quantitative indicators:
# r = ds.get_result_analysis(inputs,results)
#
# c = ds.CostExPost(inputs,results)
#
# aa['FlowIn'] = aa['FlowIn'] - aa['FlowOut']
# aa.drop('FlowOut',axis=1,inplace=True)
# aa.to_csv('aa.csv')
# r['UnitData'].to_csv('UnitData.csv')
# inputs['units'].to_csv('UnitInputs.csv')
#
# operation = ds.get_units_operation_cost(inputs,results)

