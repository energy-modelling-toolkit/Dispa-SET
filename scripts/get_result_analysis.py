# -*- coding: utf-8 -*-
"""
Not part of Dispa-SET sidetools module
@author: Andrea Mangipinto
"""

# %% Imports

# Add the root folder of Dispa-SET to the path so that the library can be loaded:
import sys, os

sys.path.append(os.path.abspath(r'../..'))

# Import Dispa-SET
import dispaset as ds

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

pd.options.display.float_format = '{:.4g}'.format

import matplotlib.ticker as mtick
import pickle
# import enlopy as el
#
# import seaborn as sns

plt.style.use('default')

# Path where saving the plot figures
output_folder = r"E:\Projects\ARES-Africa\Results\Figures"

#%% Load the inputs and the results/result analysis of the simulation
#  (Time consuming, after first time load pickle files at the bottom)

inputs = {}
results = {}

inputs['R_H'],results['R_H'] = ds.get_sim_results(path=r"E:\Projects\Github\Dispa-SET\Simulations\ARES_ALLPP_2015",cache=False)
inputs['R_W'],results['R_W'] = ds.get_sim_results(path=r"E:\Projects\Github\Dispa-SET\Simulations\ARES_ALLPP_2015_Wet",cache=False)
inputs['R_D'],results['R_D'] = ds.get_sim_results(path=r"E:\Projects\Github\Dispa-SET\Simulations\ARES_ALLPP_2015_Dry",cache=False)
inputs['RN_H'],results['RN_H'] = ds.get_sim_results(path=r"E:\Projects\Github\Dispa-SET\Simulations\ARES_ALLPP_2015_NTC",cache=False)
inputs['RN_W'],results['RN_W'] = ds.get_sim_results(path=r"E:\Projects\Github\Dispa-SET\Simulations\ARES_ALLPP_2015_NTC_Wet",cache=False)
inputs['RN_D'],results['RN_D'] = ds.get_sim_results(path=r"E:\Projects\Github\Dispa-SET\Simulations\ARES_ALLPP_2015_NTC_Dry",cache=False)

inputs['T_H'],results['T_H'] = ds.get_sim_results(path=r"E:\Projects\Github\Dispa-SET\Simulations\ARES_ALLPP_TEMBA",cache=False)
inputs['T_W'],results['T_W'] = ds.get_sim_results(path=r"E:\Projects\Github\Dispa-SET\Simulations\ARES_ALLPP_TEMBA_Wet",cache=False)
inputs['T_D'],results['T_D'] = ds.get_sim_results(path=r"E:\Projects\Github\Dispa-SET\Simulations\ARES_ALLPP_TEMBA_Dry",cache=False)
inputs['TN_H'],results['TN_H'] = ds.get_sim_results(path=r"E:\Projects\Github\Dispa-SET\Simulations\ARES_ALLPP_TEMBA_NTC",cache=False)
inputs['TN_W'],results['TN_W'] = ds.get_sim_results(path=r"E:\Projects\Github\Dispa-SET\Simulations\ARES_ALLPP_TEMBA_NTC_Wet",cache=False)
inputs['TN_D'],results['TN_D'] = ds.get_sim_results(path=r"E:\Projects\Github\Dispa-SET\Simulations\ARES_ALLPP_TEMBA_NTC_Dry",cache=False)

scenarios = list(results.keys())
# Analyse the results for each country and provide quantitative indicators:
r = {}
costs = {}
operation = {}
for s in scenarios:
    r[s] = ds.get_result_analysis(inputs[s],results[s])
    costs[s] = ds.CostExPost(inputs[s],results[s])
    operation[s] = ds.get_units_operation_cost(inputs[s], results[s])


with open('ARES_APP_Results.p', 'wb') as handle:
    pickle.dump(inputs, handle, protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(r, handle, protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(costs, handle, protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(operation, handle, protocol=pickle.HIGHEST_PROTOCOL)

aa = stop
with open('ARES_APP_Results.p', 'rb') as handle:
    inputs, results = pickle.load(handle)

# %% Set and list definition
# tech_fuel_raw = list(inputs_R_H['units'].iloc[:, 0].astype(str).str[3:].unique())

tech_fuel = ['HDAM_WAT', 'HROR_WAT', 'HPHS_WAT', 'WTON_WIN', 'WTOF_WIN', 'PHOT_SUN', 'STUR_SUN', 'COMC_BIO', 'STUR_BIO',
             'STUR_NUC', 'STUR_GEO', 'STUR_GAS', 'GTUR_GAS', 'COMC_GAS', 'STUR_OIL', 'COMC_OIL', 'STUR_HRD', 'STUR_LIG',
             'STUR_BIO_CHP', 'STUR_GAS_CHP', 'STUR_HRD_CHP', 'STUR_OIL_CHP', 'GTUR_GAS_CHP', 'STUR_LIG_CHP',
             'COMC_GAS_CHP', 'P2HT_OTH', 'Heat Slack']

tech_fuel_heat = ['STUR_BIO_CHP', 'STUR_GAS_CHP', 'STUR_HRD_CHP', 'STUR_OIL_CHP', 'GTUR_GAS_CHP', 'STUR_LIG_CHP',
                  'COMC_GAS_CHP', 'P2HT_OTH', 'Heat Slack', 'Storage Losses']
tech_fuel_storage = ['HDAM_WAT', 'HROR_WAT', 'HPHS_WAT', 'BEVS_OTH']
colors_tech = {}
for tech in tech_fuel:
    if 'Slack' not in tech:
        if 'CHP' not in tech:
            colors_tech[tech] = ds.commons['colors'][tech[-3:]]
        if 'CHP' in tech:
            colors_tech[tech] = ds.commons['colors'][tech[5:-4]]

colors_tech['HDAM_WAT'] = 'darkblue'
colors_tech['HROR_WAT'] = 'blue'
colors_tech['WTON_WIN'] = 'aquamarine'
colors_tech['STUR_SUN'] = 'yellow'
colors_tech['COMC_BIO'] = 'lime'
colors_tech['GTUR_GAS'] = 'brown'
colors_tech['STUR_GAS'] = 'maroon'
colors_tech['STUR_OIL'] = 'magenta'
colors_tech['COMC_GAS_CHP'] = 'tan'
colors_tech['GTUR_GAS_CHP'] = 'tomato'
colors_tech['STUR_GAS_CHP'] = 'lightcoral'
colors_tech['STUR_BIO_CHP'] = 'lawngreen'
colors_tech['STUR_HRD_CHP'] = 'darkviolet'
colors_tech['STUR_OIL_CHP'] = 'pink'
colors_tech['STUR_LIG_CHP'] = 'violet'
colors_tech['P2HT_ELE'] = 'gold'
colors_tech['Backup Heater'] = 'darkorange'
colors_tech['Shed Load'] = 'orangered'
colors_tech['Lost Load'] = 'darkred'
colors_tech['Curtailment'] = 'salmon'
colors_tech['LostLoad'] = colors_tech['Lost Load']
colors_tech['CostLoadShedding'] = colors_tech['Shed Load']

# %% ########################### Analysis and charts #################################
# Helper functions
def filter_by_zone(inputs,zones):
    tmp_z = inputs['units']['Zone'].isin(zones)
    return tmp_z

def filter_by_tech(inputs,techs):
    tmp_t = inputs['units']['Technology'].isin(techs)
    return tmp_t

def filter_by_fuel(inputs,fuels):
    tmp_f = inputs['units']['Fuel'].isin(fuels)
    return tmp_f

def get_demand(inputs, results, zones = None):
    demand = {}
    if zones is None:
        zones = inputs['sets']['n']
    else:
        zones = zones
    for z in zones:
        if 'OutputPowerConsumption' in results:
            loc = inputs['units']['Zone']
            demand_p2h = results['OutputPowerConsumption'].loc[:,
                         [u for u in results['OutputPowerConsumption'].columns if loc[u] == z]]
            demand_p2h = demand_p2h.sum(axis=1)
        else:
            demand_p2h = pd.Series(0, index=results['OutputPower'].index)
        demand_da = inputs['param_df']['Demand'][('DA', z)]
        demand[z] = pd.DataFrame(demand_da + demand_p2h, columns=[('DA', z)])
    demand = pd.concat(demand, axis=1)
    demand.columns = demand.columns.droplevel(-1)
    return demand

# %% Curtailment
def get_stats(inputs,results,variable):
    scenarios = list(results.keys())
    tmp = pd.DataFrame(index=['Total','Mean', 'Median', 'Max', 'Hours'], columns=scenarios)
    for s in scenarios:
        tmp.loc['Total',s] = results[s][variable].sum().sum()
        tmp.loc['Max',s] = results[s][variable].sum(axis=1).max()
        tmp.loc['Mean',s] = results[s][variable].sum(axis=1).mean()
        tmp.loc['Median', s] = results[s][variable].sum(axis=1).median()
        tmp.loc['Hours', s] = (results[s][variable].sum(axis=1) != 0).sum()
    return tmp

curtailment = get_stats(inputs, results, 'OutputCurtailedPower')
shed_load = get_stats(inputs, results, 'OutputShedLoad')

# Curtailment % with regard to V-RES generation
demand = pd.DataFrame(index=['DA'], columns=scenarios)
en_vres = pd.DataFrame(index=['OutputPower'], columns=scenarios)
zones = inputs['R_H']['sets']['n']
tech_res = ['HROR','WTON','WTOF','PHOT']
for s in scenarios:
    demand.loc['DA',s] = get_demand(inputs[s],results[s],zones).sum().sum()
    demand.loc['DA_max',s] = get_demand(inputs[s],results[s],zones).sum(axis=1).max()
    en_vres.loc['OutputPower',s] = results[s]['OutputPower'].loc[:,(filter_by_zone(inputs[s],zones)) & (filter_by_tech(inputs[s],tech_res))].sum().sum()
    curtailment.loc['Percentage total VRES',s] = curtailment.loc['Total',s] / (
                en_vres.loc['OutputPower',s] + curtailment.loc['Total',s]) * 100
    curtailment.loc['Percentage peak VRES', s] = curtailment.loc['Max', s] / (
                en_vres.loc['OutputPower', s] + curtailment.loc['Max', s]) * 100
    shed_load.loc['Percentage total Load',s] = shed_load.loc['Total',s] / (
                demand.loc['DA',s] + shed_load.loc['Total',s]) * 100
    shed_load.loc['Percentage peak Load', s] = shed_load.loc['Max', s] / (
                demand.loc['DA_max', s] + shed_load.loc['Max', s]) * 100

# %% Curtailment plot

ax = curtailment.loc['Percentage total VRES', :].T.plot(kind='bar', color='salmon', width=0.3, rot=0, position=1,
                                                         legend=True, figsize=(14, 5), fontsize=15)
ax2 = curtailment.loc['Percentage peak VRES', :].T.plot(kind='bar', color='firebrick', secondary_y=True, rot=0, ax=ax,
                                                         width=0.3, position=0, legend=True, mark_right=False,
                                                         fontsize=15)

ax.set_ylabel('Percentage Total VRES', fontsize=15)
ax2.set_ylabel('Percentage Peak VRES', fontsize=15)
ax.set_xlabel('Scenario', fontsize=15)
ax2.set_title("Curtailment", fontsize=15)

ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
ax2.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))

fig = ax2.get_figure()
fig.savefig(output_folder + "\Curtailment.png", bbox_inches='tight')
fig.clf()

# %% Shed Load plot

ax = shed_load.loc['Percentage total Load', :].T.plot(kind='bar', color='orangered', width=0.3, rot=0, figsize=(14, 5),
                                                       position=1, legend=True, fontsize=15)
ax2 = shed_load.loc['Percentage peak Load', :].T.plot(kind='bar', color='firebrick', secondary_y=True, rot=0, ax=ax,
                                                       width=0.3, position=0, legend=True, mark_right=False,
                                                       fontsize=15)

ax.set_ylabel('Percentage Total Load', fontsize=15)
ax2.set_ylabel('Percentage Peak Load', fontsize=15)
ax.set_xlabel('Scenario', fontsize=15)
ax2.set_title("Shed Load", fontsize=15)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=3))
ax2.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))

fig = ax2.get_figure()
fig.savefig(output_folder + "\Shed Load.png", bbox_inches='tight')
fig.clf()

# %% CO2 emissions
co2 = pd.DataFrame(columns=scenarios)
fossil = ['GAS','HRD','LIG','OIL','OTH','BIO','WST']
techs = ['GTUR','STUR','COMC','ICEN']
for s in scenarios:
    for f in fossil:
        for t in techs:
            idx = t + '_' + f
            co2.loc[idx ,s] = r[s]['UnitData']['CO2 [t]'].loc[(filter_by_zone(inputs[s],zones)) &
                                              (filter_by_fuel(inputs[s],[f]) & filter_by_tech(inputs[s], [t]))].sum()
co2 = co2[(co2.T != 0).any()] / 1e6


# %% co2 emissions plot
ax = co2.T.plot(kind='bar', stacked=True, rot=0, color=[colors_tech.get(x) for x in co2.index], legend=False,
                   figsize=(14, 5), fontsize=15)

# for container, hatch in zip(ax.containers,
#                             ("", '', '', '', '', '', '', '///', '///', '///', '///', '///', '///', '///', '///')):
#     for patch in container.patches:
#         patch.set_hatch(hatch)

ax.set_ylabel('[$CO_{2}$ [Mton]', fontsize=15)
ax.set_xlabel('Scenario', fontsize=15)
ax.set_title("$CO_{2}$ emissions", fontsize=15)
handles, labels = ax.get_legend_handles_labels()
ax.legend(reversed(handles), reversed(labels), loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fontsize=15)
ax.set_ylim(0, co2.sum().max()*1.05)

fig = ax.get_figure()
fig.savefig(output_folder + "\CO2 emissions.png", bbox_inches='tight')
fig.clf()

# %% Costs breakdown



costs_raw_ALL = pd.DataFrame(costs_ALL[0].sum(axis=0).values, columns=['ALLFLEX'], index=costs_ALL[0].columns).T
costs_raw_EV = pd.DataFrame(costs_EV[0].sum(axis=0).values, columns=['EVFLEX'], index=costs_EV[0].columns).T
costs_raw_TH = pd.DataFrame(costs_TH[0].sum(axis=0).values, columns=['THFLEX'], index=costs_TH[0].columns).T
costs_raw_NO = pd.DataFrame(costs_NO[0].sum(axis=0).values, columns=['NOFLEX'], index=costs_NO[0].columns).T
costs_raw_HY = pd.DataFrame(costs_HY[0].sum(axis=0).values, columns=['HYFLEX'], index=costs_HY[0].columns).T

costs = pd.DataFrame(index=scenarios, columns=costs_ALL[0].columns)

costs.loc['ALLFLEX', :] = costs_raw_ALL.values / 1e9
costs.loc['EVFLEX', :] = costs_raw_EV.values / 1e9
costs.loc['THFLEX', :] = costs_raw_TH.values / 1e9
costs.loc['NOFLEX', :] = costs_raw_NO.values / 1e9
costs.loc['HYFLEX', :] = costs_raw_HY.values / 1e9

for c in costs.columns:
    if costs[c].values.sum() == 0:
        costs.drop(columns=[c], inplace=True)

costs_var_raw_ALL = pd.DataFrame(
    (results_ALL['OutputPower'] * inputs_ALL['param_df']['CostVariable']).fillna(0).sum(axis=0).values,
    columns=['ALLFLEX'],
    index=(results_ALL['OutputPower'] * inputs_ALL['param_df']['CostVariable']).fillna(0).sum(axis=0).index).T
costs_var_raw_EV = pd.DataFrame(
    (results_EV['OutputPower'] * inputs_EV['param_df']['CostVariable']).fillna(0).sum(axis=0).values,
    columns=['EVFLEX'],
    index=(results_EV['OutputPower'] * inputs_EV['param_df']['CostVariable']).fillna(0).sum(axis=0).index).T
costs_var_raw_TH = pd.DataFrame(
    (results_TH['OutputPower'] * inputs_TH['param_df']['CostVariable']).fillna(0).sum(axis=0).values,
    columns=['THFLEX'],
    index=(results_TH['OutputPower'] * inputs_TH['param_df']['CostVariable']).fillna(0).sum(axis=0).index).T
costs_var_raw_NO = pd.DataFrame(
    (results_NO['OutputPower'] * inputs_NO['param_df']['CostVariable']).fillna(0).sum(axis=0).values,
    columns=['NOFLEX'],
    index=(results_NO['OutputPower'] * inputs_NO['param_df']['CostVariable']).fillna(0).sum(axis=0).index).T
costs_var_raw_HY = pd.DataFrame(
    (results_HY['OutputPower'] * inputs_HY['param_df']['CostVariable']).fillna(0).sum(axis=0).values,
    columns=['HYFLEX'],
    index=(results_HY['OutputPower'] * inputs_HY['param_df']['CostVariable']).fillna(0).sum(axis=0).index).T

costs_var = pd.DataFrame(index=scenarios, columns=tech_fuel)

for t in tech_fuel:
    costs_var.loc['ALLFLEX', t] = costs_var_raw_ALL.loc[:, costs_var_raw_ALL.columns.str.endswith(t)].values.sum() / 1e9
    costs_var.loc['EVFLEX', t] = costs_var_raw_EV.loc[:, costs_var_raw_EV.columns.str.endswith(t)].values.sum() / 1e9
    costs_var.loc['THFLEX', t] = costs_var_raw_TH.loc[:, costs_var_raw_TH.columns.str.endswith(t)].values.sum() / 1e9
    costs_var.loc['NOFLEX', t] = costs_var_raw_NO.loc[:, costs_var_raw_NO.columns.str.endswith(t)].values.sum() / 1e9
    costs_var.loc['HYFLEX', t] = costs_var_raw_HY.loc[:, costs_var_raw_HY.columns.str.endswith(t)].values.sum() / 1e9

for c in costs_var.columns:
    if costs_var[c].values.sum() == 0:
        costs_var.drop(columns=[c], inplace=True)

costs_breakdown = pd.concat([costs, costs_var], axis=1, sort=False)

costs_breakdown.drop(columns=['CostVariable'], inplace=True)
costs_breakdown.drop(columns=['CostHeat'], inplace=True)

cols = ['CostHeatSlack', 'STUR_BIO_CHP', 'STUR_GAS_CHP', 'STUR_HRD_CHP', 'STUR_OIL_CHP', 'GTUR_GAS_CHP', 'STUR_LIG_CHP',
        'COMC_GAS_CHP',
        'COMC_BIO', 'STUR_BIO', 'STUR_NUC', 'STUR_GAS', 'GTUR_GAS', 'COMC_GAS', 'STUR_OIL', 'COMC_OIL', 'STUR_HRD',
        'STUR_LIG',
        'CostStartUp', 'CostRampUp', 'Spillage', 'CostLoadShedding', 'LostLoad']

costs_breakdown = costs_breakdown[cols]
costs_breakdown.rename(columns={"CostHeatSlack": "CostBackupHeater"}, inplace=True)
costs_breakdown.rename(columns={"CostLoadShedding": "CostShedLoad"}, inplace=True)

costs_breakdown_plot = costs_breakdown
costs_breakdown_plot.drop(columns=['Spillage'], inplace=True)
costs_breakdown_plot.drop(columns=['LostLoad'], inplace=True)

costs_chp = costs_breakdown.loc[:, costs_breakdown.columns.str.endswith('CHP')].sum(axis=1)

tot_syst_cost = pd.DataFrame(index=['Total system Cost'], columns=scenarios)
tot_syst_cost.loc[:, 'ALLFLEX'] = costs_raw_ALL.sum(axis=1).values / 1e9
tot_syst_cost.loc[:, 'EVFLEX'] = costs_raw_EV.sum(axis=1).values / 1e9
tot_syst_cost.loc[:, 'NOFLEX'] = costs_raw_NO.sum(axis=1).values / 1e9
tot_syst_cost.loc[:, 'THFLEX'] = costs_raw_TH.sum(axis=1).values / 1e9
tot_syst_cost.loc[:, 'HYFLEX'] = costs_raw_HY.sum(axis=1).values / 1e9

# %% Costs breakdown plot
plt.style.use('default')
ax = costs_breakdown_plot.plot(kind='bar', stacked=True, rot=0,
                               color=[colors_tech.get(x) for x in costs_breakdown_plot.columns], legend=False,
                               figsize=(10, 5), fontsize=15)

for container, hatch in zip(ax.containers, ('///', '///', '///', '///', '///', '///', '///', '///',)):
    for patch in container.patches:
        patch.set_hatch(hatch)

ax.set_ylabel('Costs [Billion €]', fontsize=15)
ax.set_xlabel('Scenario', fontsize=15)
ax.set_title("Total System Costs", fontsize=15)
handles, labels = ax.get_legend_handles_labels()
ax.legend(reversed(handles), reversed(labels), loc='center left', bbox_to_anchor=(1, 0.5), ncol=2, fontsize=15)

fig = ax.get_figure()
fig.savefig(output_folder + "\Costs Breakdown.png", bbox_inches='tight')
fig.clf()

# # %% Shifted Load
# shift_load = pd.DataFrame(columns = scenarios)
# for s in scenarios:
#     shift_load.loc['diff'.s] = results[s]['OutputStorageLevel']
#
# sto_cap_ALL = pd.DataFrame(inputs_ALL['units']['StorageCapacity'][inputs_ALL['units']['StorageCapacity'] > 0].values,
#                            columns=["ALLFLEX"], index=list(
#         inputs_ALL['units']['StorageCapacity'][inputs_ALL['units']['StorageCapacity'] > 0].index)).T
# sto_cap_ALL_tot = sto_cap_ALL * inputs_ALL['units'].loc[sto_cap_ALL.columns, 'Nunits']
# sto_cap_EV = pd.DataFrame(inputs_EV['units']['StorageCapacity'][inputs_EV['units']['StorageCapacity'] > 0].values,
#                           columns=["EVFLEX"], index=list(
#         inputs_EV['units']['StorageCapacity'][inputs_EV['units']['StorageCapacity'] > 0].index)).T
# sto_cap_TH = pd.DataFrame(inputs_TH['units']['StorageCapacity'][inputs_TH['units']['StorageCapacity'] > 0].values,
#                           columns=["THFLEX"], index=list(
#         inputs_TH['units']['StorageCapacity'][inputs_TH['units']['StorageCapacity'] > 0].index)).T
# sto_cap_TH_tot = sto_cap_TH * inputs_TH['units'].loc[sto_cap_TH.columns, 'Nunits']
# sto_cap_NO = pd.DataFrame(inputs_NO['units']['StorageCapacity'][inputs_NO['units']['StorageCapacity'] > 0].values,
#                           columns=["NOFLEX"], index=list(
#         inputs_NO['units']['StorageCapacity'][inputs_NO['units']['StorageCapacity'] > 0].index)).T
# sto_cap_HY = pd.DataFrame(inputs_HY['units']['StorageCapacity'][inputs_HY['units']['StorageCapacity'] > 0].values,
#                           columns=["HYFLEX"], index=list(
#         inputs_HY['units']['StorageCapacity'][inputs_HY['units']['StorageCapacity'] > 0].index)).T
#
# delta_ALL = results_ALL['OutputStorageLevel'].diff()
# shift_load_raw_ALL = pd.DataFrame(delta_ALL[delta_ALL > 0].sum(), columns=['ALLFLEX'],
#                                   index=delta_ALL[delta_ALL > 0].columns).T
# stor_usage_raw_ALL = pd.DataFrame(-delta_ALL[delta_ALL < 0].sum(), columns=['ALLFLEX'],
#                                   index=delta_ALL[delta_ALL < 0].columns).T
# delta_EV = results_EV['OutputStorageLevel'].diff()
# shift_load_raw_EV = pd.DataFrame(delta_EV[delta_EV > 0].sum(), columns=['EVFLEX'],
#                                  index=delta_EV[delta_EV > 0].columns).T
# stor_usage_raw_EV = pd.DataFrame(-delta_EV[delta_EV < 0].sum(), columns=['EVFLEX'],
#                                  index=delta_EV[delta_EV < 0].columns).T
# delta_TH = results_TH['OutputStorageLevel'].diff()
# shift_load_raw_TH = pd.DataFrame(delta_TH[delta_TH > 0].sum(), columns=['THFLEX'],
#                                  index=delta_TH[delta_TH > 0].columns).T
# stor_usage_raw_TH = pd.DataFrame(-delta_TH[delta_TH < 0].sum(), columns=['THFLEX'],
#                                  index=delta_TH[delta_TH < 0].columns).T
# delta_NO = results_NO['OutputStorageLevel'].diff()
# shift_load_raw_NO = pd.DataFrame(delta_NO[delta_NO > 0].sum(), columns=['NOFLEX'],
#                                  index=delta_NO[delta_NO > 0].columns).T
# stor_usage_raw_NO = pd.DataFrame(-delta_NO[delta_NO < 0].sum(), columns=['NOFLEX'],
#                                  index=delta_NO[delta_NO < 0].columns).T
# delta_HY = results_HY['OutputStorageLevel'].diff()
# shift_load_raw_HY = pd.DataFrame(delta_HY[delta_HY > 0].sum(), columns=['HYFLEX'],
#                                  index=delta_HY[delta_HY > 0].columns).T
# stor_usage_raw_HY = pd.DataFrame(-delta_HY[delta_HY < 0].sum(), columns=['HYFLEX'],
#                                  index=delta_HY[delta_HY < 0].columns).T
#
# shift_load = pd.DataFrame(index=scenarios, columns=tech_fuel)
#
# for t in tech_fuel:
#     shift_load.loc['EVFLEX', t] = shift_load_raw_EV.loc[:, shift_load_raw_EV.columns.str.endswith(t)].values.sum() / 1e6
#     shift_load.loc['HYFLEX', t] = shift_load_raw_HY.loc[:, shift_load_raw_HY.columns.str.endswith(t)].values.sum() / 1e6
#
# for t in tech_fuel_heat:
#     shift_load.loc['ALLFLEX', t] = shift_load_raw_ALL.loc[:,
#                                    shift_load_raw_ALL.columns.str.endswith(t)].values.sum() / 1e6
#     shift_load.loc['THFLEX', t] = shift_load_raw_TH.loc[:, shift_load_raw_TH.columns.str.endswith(t)].values.sum() / 1e6
#     shift_load.loc['NOFLEX', t] = shift_load_raw_NO.loc[:, shift_load_raw_NO.columns.str.endswith(t)].values.sum() / 1e6
#
# for t in tech_fuel_storage:
#     shift_load.loc['ALLFLEX', t] = shift_load_raw_ALL.loc[:,
#                                    shift_load_raw_ALL.columns.str.endswith(t)].values.sum() / 1e6
#
# shift_load.fillna(0, inplace=True)
#
# for c in shift_load.columns:
#     if shift_load[c].values.sum() == 0:
#         shift_load.drop(columns=[c], inplace=True)
#
# shift_load.rename(columns={"P2HT_OTH": "P2HT_ELE"}, inplace=True)
#
# sto_cycle_ALL = stor_usage_raw_ALL / sto_cap_ALL_tot
# sto_cycle_EV_old = en_prod_raw_EV[sto_cap_EV.columns] / sto_cap_EV
# sto_cycle_EV = stor_usage_raw_EV / sto_cap_EV
# sto_cycle_TH = stor_usage_raw_TH / sto_cap_TH_tot
# sto_cycle_NO = stor_usage_raw_NO / sto_cap_NO
# sto_cycle_HY_old = en_prod_raw_HY[sto_cap_HY.columns] / sto_cap_HY
# sto_cycle_HY = stor_usage_raw_HY / sto_cap_HY
#
# # %% Shifted load plot
#
# ax = shift_load.plot(kind='bar', stacked=True, rot=0, color=[colors_tech.get(x) for x in shift_load.columns],
#                      legend=False, figsize=(7, 5), fontsize=15)
#
# for container, hatch in zip(ax.containers, ('', '', '', "///", '///', '///', '///', '///', '///', '///', '///')):
#     for patch in container.patches:
#         patch.set_hatch(hatch)
#
# ax.set_ylabel('Energy [TWh]', fontsize=15)
# ax.set_xlabel('Scenario', fontsize=15)
# ax.set_title("Shifted Load", fontsize=15)
# handles, labels = ax.get_legend_handles_labels()
#
# ax.legend(reversed(handles), reversed(labels), loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fontsize=15)
#
# fig = ax.get_figure()
# fig.savefig(output_folder + "\Load Shift.png", bbox_inches='tight')

# # %% Energy output
#
# def get_en_prod(results, scenario='None'):
#     # Extract raw energy prodction from results
#     en_prod_raw = pd.DataFrame(results['OutputPower'].sum(axis=0).values, columns=[scenario],
#                                index=results['OutputPower'].columns).T
#     # Define df with grouped energy production
#     en_prod = pd.DataFrame(columns=tech_fuel, index=[scenario])
#
#     # Fill the df
#     if isinstance(results['LostLoad_WaterSlack'], int):
#         en_prod.loc[scenario, 'Lost Load'] = -(
#                     results['LostLoad_2D'].sum().sum() + results['LostLoad_2U'].sum().sum() + results[
#                 'LostLoad_3U'].sum().sum() + results['LostLoad_MaxPower'].sum().sum() + results[
#                         'LostLoad_MinPower'].sum().sum() + results['LostLoad_RampDown'].sum().sum() + results[
#                         'LostLoad_RampUp'].sum().sum() + results['LostLoad_WaterSlack']) / 10 ** 6
#     else:
#         en_prod.loc[scenario, 'Lost Load'] = -(
#                     results['LostLoad_2D'].sum().sum() + results['LostLoad_2U'].sum().sum() + results[
#                 'LostLoad_3U'].sum().sum() + results['LostLoad_MaxPower'].sum().sum() + results[
#                         'LostLoad_MinPower'].sum().sum() + results['LostLoad_RampDown'].sum().sum() + results[
#                         'LostLoad_RampUp'].sum().sum() + results['LostLoad_WaterSlack'].sum()) / 10 ** 6
#     en_prod.loc[scenario, 'Shed Load'] = -results['OutputShedLoad'].sum().sum() / 10 ** 6
#     en_prod.loc[scenario, 'Curtailment'] = -results['OutputCurtailedPower'].sum().sum() / 10 ** 6
#
#     for t in tech_fuel:
#         en_prod.loc[scenario, t] = en_prod_raw.loc[:, en_prod_raw.columns.str.endswith(t)].values.sum() / 1e6
#
#     for c in en_prod.columns:
#         if en_prod[c].values.sum() == 0:
#             en_prod.drop(columns=[c], inplace=True)
#
#     # Reordering the order of columns to be coherent in the plot
#     en_prod_list = list(en_prod.columns)
#     en_prod_list = en_prod_list[-3:] + en_prod_list[:-3]
#     en_prod = en_prod[en_prod_list]
#
#     return en_prod
#
#
# # Example to create a df with en_prod for all 5 scenarios
#
# en_prod = pd.DataFrame(index=scenarios, columns=get_en_prod(results['R_H'], 'R_H').columns)
#
# for scenario in scenarios:
#     en_prod.loc[scenario, :] = get_en_prod(results[scenario], scenario).iloc[0, :]
#
# en_prod.fillna(0, inplace=True)
#
# # % Increase in RES absorption
# en_prod_sun = en_prod.loc[:, en_prod.columns.str.contains('_SUN')].sum(axis=1)
# en_prod_win = en_prod.loc[:, en_prod.columns.str.contains('_WIN')].sum(axis=1)
# en_prod_bio = en_prod.loc[:, en_prod.columns.str.contains('_BIO')].sum(axis=1)
#
# # Delta in ALLFLEX
# delta_en_prod_sun_ALL = (en_prod_sun.loc['ALLFLEX'] - en_prod_sun.loc['NOFLEX']) / en_prod_sun.loc['NOFLEX']
# delta_en_prod_win_ALL = (en_prod_win.loc['ALLFLEX'] - en_prod_win.loc['NOFLEX']) / en_prod_win.loc['NOFLEX']
# delta_en_prod_bio_ALL = (en_prod_bio.loc['ALLFLEX'] - en_prod_bio.loc['NOFLEX']) / en_prod_bio.loc['NOFLEX']
#
# # Delta in EVFLEX
# delta_en_prod_sun_EV = (en_prod_sun.loc['EVFLEX'] - en_prod_sun.loc['NOFLEX']) / en_prod_sun.loc['NOFLEX']
# delta_en_prod_win_EV = (en_prod_win.loc['EVFLEX'] - en_prod_win.loc['NOFLEX']) / en_prod_win.loc['NOFLEX']
# delta_en_prod_bio_EV = (en_prod_bio.loc['EVFLEX'] - en_prod_bio.loc['NOFLEX']) / en_prod_bio.loc['NOFLEX']
#
# # Delta in HYFLEX
# delta_en_prod_sun_HY = (en_prod_sun.loc['HYFLEX'] - en_prod_sun.loc['NOFLEX']) / en_prod_sun.loc['NOFLEX']
# delta_en_prod_win_HY = (en_prod_win.loc['HYFLEX'] - en_prod_win.loc['NOFLEX']) / en_prod_win.loc['NOFLEX']
# delta_en_prod_bio_HY = (en_prod_bio.loc['HYFLEX'] - en_prod_bio.loc['NOFLEX']) / en_prod_bio.loc['NOFLEX']
#
# # Energy Mix
#
# demand = inputs_ALL['param_df']['Demand']['DA']
# demand_p2h_raw_ALL = results_ALL['OutputPowerConsumption']
# demand_p2h_raw_EV = results_EV['OutputPowerConsumption']
# demand_p2h_raw_TH = results_TH['OutputPowerConsumption']
# demand_p2h_raw_NO = results_NO['OutputPowerConsumption']
# demand_p2h_raw_HY = results_HY['OutputPowerConsumption']
#
# sto_input_raw_ALL = results_ALL['OutputStorageInput']
# sto_input_raw_EV = results_EV['OutputStorageInput']
# sto_input_raw_TH = results_TH['OutputStorageInput']
# sto_input_raw_NO = results_NO['OutputStorageInput']
# sto_input_raw_HY = results_HY['OutputStorageInput']
#
# shed_load_ALL = pd.DataFrame(
#     results_ALL['OutputShedLoad'].loc[:, ~ results_ALL['OutputShedLoad'].columns.str.contains('x|y')],
#     index=demand.index, columns=inputs_ALL['param_df']['sets']['n'])
# shed_load_ALL.fillna(0, inplace=True)
# shed_load_EV = pd.DataFrame(
#     results_EV['OutputShedLoad'].loc[:, ~ results_EV['OutputShedLoad'].columns.str.contains('x|y')], index=demand.index,
#     columns=inputs_EV['param_df']['sets']['n'])
# shed_load_EV.fillna(0, inplace=True)
# shed_load_TH = pd.DataFrame(results_TH['OutputShedLoad'], index=demand.index,
#                             columns=inputs_TH['param_df']['sets']['n'])
# shed_load_TH.fillna(0, inplace=True)
# shed_load_NO = pd.DataFrame(results_NO['OutputShedLoad'], index=demand.index,
#                             columns=inputs_NO['param_df']['sets']['n'])
# shed_load_NO.fillna(0, inplace=True)
# shed_load_HY = pd.DataFrame(results_HY['OutputShedLoad'], index=demand.index,
#                             columns=inputs_HY['param_df']['sets']['n'])
# shed_load_HY.fillna(0, inplace=True)
#
# ll_max_ALL = pd.DataFrame(results_ALL['LostLoad_MaxPower'], index=demand.index,
#                           columns=inputs_ALL['param_df']['sets']['n'])
# ll_max_ALL.fillna(0, inplace=True)
# ll_max_EV = pd.DataFrame(results_EV['LostLoad_MaxPower'], index=demand.index,
#                          columns=inputs_EV['param_df']['sets']['n'])
# ll_max_EV.fillna(0, inplace=True)
# ll_max_TH = pd.DataFrame(results_TH['LostLoad_MaxPower'], index=demand.index,
#                          columns=inputs_TH['param_df']['sets']['n'])
# ll_max_TH.fillna(0, inplace=True)
# ll_max_NO = pd.DataFrame(results_NO['LostLoad_MaxPower'], index=demand.index,
#                          columns=inputs_NO['param_df']['sets']['n'])
# ll_max_NO.fillna(0, inplace=True)
# ll_max_HY = pd.DataFrame(results_HY['LostLoad_MaxPower'], index=demand.index,
#                          columns=inputs_HY['param_df']['sets']['n'])
# ll_max_HY.fillna(0, inplace=True)
#
# ll_min_ALL = pd.DataFrame(results_ALL['LostLoad_MinPower'], index=demand.index,
#                           columns=inputs_ALL['param_df']['sets']['n'])
# ll_min_ALL.fillna(0, inplace=True)
# ll_min_EV = pd.DataFrame(results_EV['LostLoad_MinPower'], index=demand.index,
#                          columns=inputs_EV['param_df']['sets']['n'])
# ll_min_EV.fillna(0, inplace=True)
# ll_min_TH = pd.DataFrame(results_TH['LostLoad_MinPower'], index=demand.index,
#                          columns=inputs_TH['param_df']['sets']['n'])
# ll_min_TH.fillna(0, inplace=True)
# ll_min_NO = pd.DataFrame(results_NO['LostLoad_MinPower'], index=demand.index,
#                          columns=inputs_NO['param_df']['sets']['n'])
# ll_min_NO.fillna(0, inplace=True)
# ll_min_HY = pd.DataFrame(results_HY['LostLoad_MinPower'], index=demand.index,
#                          columns=inputs_HY['param_df']['sets']['n'])
# ll_min_HY.fillna(0, inplace=True)
#
# demand_p2h_ALL = pd.DataFrame(index=demand_p2h_raw_ALL.index, columns=inputs_ALL['param_df']['sets']['n'])
# demand_p2h_EV = pd.DataFrame(index=demand_p2h_raw_EV.index, columns=inputs_ALL['param_df']['sets']['n'])
# demand_p2h_TH = pd.DataFrame(index=demand_p2h_raw_TH.index, columns=inputs_ALL['param_df']['sets']['n'])
# demand_p2h_NO = pd.DataFrame(index=demand_p2h_raw_NO.index, columns=inputs_ALL['param_df']['sets']['n'])
# demand_p2h_HY = pd.DataFrame(index=demand_p2h_raw_HY.index, columns=inputs_ALL['param_df']['sets']['n'])
#
# sto_input_ALL = pd.DataFrame(index=sto_input_raw_ALL.index, columns=inputs_ALL['param_df']['sets']['n'])
# sto_input_EV = pd.DataFrame(index=sto_input_raw_EV.index, columns=inputs_EV['param_df']['sets']['n'])
# sto_input_TH = pd.DataFrame(index=sto_input_raw_TH.index, columns=inputs_TH['param_df']['sets']['n'])
# sto_input_NO = pd.DataFrame(index=sto_input_raw_NO.index, columns=inputs_NO['param_df']['sets']['n'])
# sto_input_HY = pd.DataFrame(index=sto_input_raw_HY.index, columns=inputs_HY['param_df']['sets']['n'])
#
# for c in inputs_ALL['param_df']['sets']['n']:
#     demand_p2h_ALL.loc[:, c] = demand_p2h_raw_ALL.loc[:, demand_p2h_raw_ALL.columns.str.contains(' ' + c + '_')].sum(
#         axis=1)
#     demand_p2h_EV.loc[:, c] = demand_p2h_raw_EV.loc[:, demand_p2h_raw_EV.columns.str.contains(' ' + c + '_')].sum(
#         axis=1)
#     demand_p2h_TH.loc[:, c] = demand_p2h_raw_TH.loc[:, demand_p2h_raw_TH.columns.str.contains(' ' + c + '_')].sum(
#         axis=1)
#     demand_p2h_NO.loc[:, c] = demand_p2h_raw_NO.loc[:, demand_p2h_raw_NO.columns.str.contains(' ' + c + '_')].sum(
#         axis=1)
#     demand_p2h_HY.loc[:, c] = demand_p2h_raw_HY.loc[:, demand_p2h_raw_HY.columns.str.contains(' ' + c + '_')].sum(
#         axis=1)
#
# for c in inputs_ALL['param_df']['sets']['n']:
#     sto_input_ALL.loc[:, c] = sto_input_raw_ALL.loc[:, sto_input_raw_ALL.columns.str.contains(' ' + c + '_')].sum(
#         axis=1)
#     sto_input_EV.loc[:, c] = sto_input_raw_EV.loc[:, sto_input_raw_EV.columns.str.contains(' ' + c + '_')].sum(axis=1)
#     sto_input_TH.loc[:, c] = sto_input_raw_TH.loc[:, sto_input_raw_TH.columns.str.contains(' ' + c + '_')].sum(axis=1)
#     sto_input_NO.loc[:, c] = sto_input_raw_NO.loc[:, sto_input_raw_NO.columns.str.contains(' ' + c + '_')].sum(axis=1)
#     sto_input_HY.loc[:, c] = sto_input_raw_HY.loc[:, sto_input_raw_HY.columns.str.contains(' ' + c + '_')].sum(axis=1)
#
# energy_mix_den = pd.DataFrame(index=['0'], columns=scenarios)
# energy_mix_den.loc[:, 'ALLFLEX'] = (
#                                            demand + demand_p2h_ALL + sto_input_ALL - shed_load_ALL - ll_max_ALL + ll_min_ALL).sum().sum() / 1e6
# energy_mix_den.loc[:, 'EVFLEX'] = (
#                                               demand + demand_p2h_EV + sto_input_EV - shed_load_EV - ll_max_EV + ll_min_EV).sum().sum() / 1e6
# energy_mix_den.loc[:, 'THFLEX'] = (
#                                               demand + demand_p2h_TH + sto_input_TH - shed_load_TH - ll_max_TH + ll_min_TH).sum().sum() / 1e6
# energy_mix_den.loc[:, 'NOFLEX'] = (
#                                               demand + demand_p2h_NO + sto_input_NO - shed_load_NO - ll_max_NO + ll_min_NO).sum().sum() / 1e6
# energy_mix_den.loc[:, 'HYFLEX'] = (
#                                               demand + demand_p2h_HY + sto_input_HY - shed_load_HY - ll_max_HY + ll_min_HY).sum().sum() / 1e6
#
# energy_mix = pd.DataFrame(index=scenarios, columns=en_prod.columns)
#
# energy_mix.loc['ALLFLEX'] = en_prod.loc['ALLFLEX'] / energy_mix_den['ALLFLEX'].values
# energy_mix.loc['EVFLEX'] = en_prod.loc['EVFLEX'] / energy_mix_den['EVFLEX'].values
# energy_mix.loc['THFLEX'] = en_prod.loc['THFLEX'] / energy_mix_den['THFLEX'].values
# energy_mix.loc['NOFLEX'] = en_prod.loc['NOFLEX'] / energy_mix_den['NOFLEX'].values
# energy_mix.loc['HYFLEX'] = en_prod.loc['HYFLEX'] / energy_mix_den['HYFLEX'].values
#
# energy_mix = energy_mix.abs()
# # energy_mix.drop(index = ['Lost Load', 'Shed Load', 'Curtailment'], inplace = True)
#
#
# # %% Energy output plot
#
# # Use this to change what columns are not plotted from the original df
# en_prod_plot = en_prod
# en_prod_plot.drop(columns=['BEVS_OTH'], inplace=True)
#
# ax = en_prod_plot.plot(kind='bar', stacked=True, rot=0, color=[colors_tech.get(x) for x in en_prod_plot.columns],
#                        legend=False, figsize=(10, 5), fontsize=15)
#
# for container, hatch in zip(ax.containers, (
# '', '', "", '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '///', '///', '///', '///', '///',
# '///', '///')):
#     for patch in container.patches:
#         patch.set_hatch(hatch)
#
# ax.axhline(y=0, xmin=-1, xmax=1, color='black', linestyle='-', lw=1)
# ax.set_ylabel('Energy [TWh]', fontsize=15)
# ax.set_xlabel('Scenario', fontsize=15)
# ax.set_title("Energy Output per Fuel / Technology", fontsize=15)
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(reversed(handles), reversed(labels), loc='center left', bbox_to_anchor=(1, 0.5), ncol=2, fontsize=15)
#
# fig = ax.get_figure()
# fig.savefig(output_folder + "\Energy output.png", bbox_inches='tight')
#
# # %% Heating output
#
# heat_prod_raw_ALL = pd.DataFrame(results_ALL['OutputHeat'].sum(axis=0).values, columns=['ALLFLEX'],
#                                  index=results_ALL['OutputHeat'].columns).T
# heat_prod_raw_ALL['Heat Slack'] = results_ALL['OutputHeatSlack'].sum().sum()
# heat_prod_raw_EV = pd.DataFrame(results_EV['OutputHeat'].sum(axis=0).values, columns=['EVFLEX'],
#                                 index=results_EV['OutputHeat'].columns).T
# heat_prod_raw_EV['Heat Slack'] = results_EV['OutputHeatSlack'].sum().sum()
# heat_prod_raw_TH = pd.DataFrame(results_TH['OutputHeat'].sum(axis=0).values, columns=['THLFLEX'],
#                                 index=results_TH['OutputHeat'].columns).T
# heat_prod_raw_TH['Heat Slack'] = results_TH['OutputHeatSlack'].sum().sum()
# heat_prod_raw_NO = pd.DataFrame(results_NO['OutputHeat'].sum(axis=0).values, columns=['NOFLEX'],
#                                 index=results_NO['OutputHeat'].columns).T
# heat_prod_raw_NO['Heat Slack'] = results_NO['OutputHeatSlack'].sum().sum()
# heat_prod_raw_HY = pd.DataFrame(results_HY['OutputHeat'].sum(axis=0).values, columns=['HYFLEX'],
#                                 index=results_HY['OutputHeat'].columns).T
# heat_prod_raw_HY['Heat Slack'] = results_HY['OutputHeatSlack'].sum().sum()
#
# # Heat storage losses
#
# heat_sto_col_ALL = inputs_ALL['units']['StorageSelfDischarge'].loc[
#     inputs_ALL['units']['StorageSelfDischarge'] > 0].index
# sto_self_dis = 0.03 / 24
# sto_loss_ALL = results_ALL['OutputStorageLevel'].loc[:,
#                results_ALL['OutputStorageLevel'].columns.str.contains('_CHP|P2HT')] * sto_self_dis
# sto_loss_EV = results_EV['OutputStorageLevel'].loc[:,
#               results_EV['OutputStorageLevel'].columns.str.contains('_CHP|P2HT')] * sto_self_dis
# sto_loss_TH = results_TH['OutputStorageLevel'].loc[:,
#               results_TH['OutputStorageLevel'].columns.str.contains('_CHP|P2HT')] * sto_self_dis
# sto_loss_NO = results_NO['OutputStorageLevel'].loc[:,
#               results_NO['OutputStorageLevel'].columns.str.contains('_CHP|P2HT')] * sto_self_dis
# sto_loss_HY = results_HY['OutputStorageLevel'].loc[:,
#               results_HY['OutputStorageLevel'].columns.str.contains('_CHP|P2HT')] * sto_self_dis
#
# heat_prod = pd.DataFrame(index=scenarios, columns=tech_fuel_heat)
#
# for t in tech_fuel_heat:
#     heat_prod.loc['ALLFLEX', t] = heat_prod_raw_ALL.loc[:, heat_prod_raw_ALL.columns.str.endswith(t)].values.sum() / 1e6
#     heat_prod.loc['EVFLEX', t] = heat_prod_raw_EV.loc[:, heat_prod_raw_EV.columns.str.endswith(t)].values.sum() / 1e6
#     heat_prod.loc['THFLEX', t] = heat_prod_raw_TH.loc[:, heat_prod_raw_TH.columns.str.endswith(t)].values.sum() / 1e6
#     heat_prod.loc['NOFLEX', t] = heat_prod_raw_NO.loc[:, heat_prod_raw_NO.columns.str.endswith(t)].values.sum() / 1e6
#     heat_prod.loc['HYFLEX', t] = heat_prod_raw_HY.loc[:, heat_prod_raw_HY.columns.str.endswith(t)].values.sum() / 1e6
#
# heat_prod.loc['ALLFLEX', 'Storage Losses'] = sto_loss_ALL.sum().sum() / 1e6
# heat_prod.loc['EVFLEX', 'Storage Losses'] = sto_loss_EV.sum().sum() / 1e6
# heat_prod.loc['THFLEX', 'Storage Losses'] = sto_loss_TH.sum().sum() / 1e6
# heat_prod.loc['NOFLEX', 'Storage Losses'] = sto_loss_NO.sum().sum() / 1e6
# heat_prod.loc['HYFLEX', 'Storage Losses'] = sto_loss_HY.sum().sum() / 1e6
#
# for c in heat_prod.columns:
#     if heat_prod[c].values.sum() == 0:
#         heat_prod.drop(columns=[c], inplace=True)
#
# heat_prod.rename(columns={"P2HT_OTH": "P2HT_ELE"}, inplace=True)
# heat_prod.rename(columns={"Heat Slack": "Backup Heater"}, inplace=True)
#
# # Energy Heat mix
#
# heat_demand_raw_ALL = inputs_ALL['param_df']['HeatDemand']
# heat_demand_raw_EV = inputs_EV['param_df']['HeatDemand']
# heat_demand_raw_TH = inputs_TH['param_df']['HeatDemand']
# heat_demand_raw_NO = inputs_NO['param_df']['HeatDemand']
# heat_demand_raw_HY = inputs_HY['param_df']['HeatDemand']
#
# heat_demand_ALL = pd.DataFrame(index=heat_demand_raw_ALL.index, columns=inputs_ALL['param_df']['sets']['n'])
# heat_demand_EV = pd.DataFrame(index=heat_demand_raw_EV.index, columns=inputs_ALL['param_df']['sets']['n'])
# heat_demand_TH = pd.DataFrame(index=heat_demand_raw_TH.index, columns=inputs_ALL['param_df']['sets']['n'])
# heat_demand_NO = pd.DataFrame(index=heat_demand_raw_NO.index, columns=inputs_ALL['param_df']['sets']['n'])
# heat_demand_HY = pd.DataFrame(index=heat_demand_raw_HY.index, columns=inputs_ALL['param_df']['sets']['n'])
#
# heat_slack_ALL = pd.DataFrame(index=heat_demand_raw_ALL.index, columns=inputs_ALL['param_df']['sets']['n'])
# heat_slack_EV = pd.DataFrame(index=heat_demand_raw_EV.index, columns=inputs_ALL['param_df']['sets']['n'])
# heat_slack_TH = pd.DataFrame(index=heat_demand_raw_TH.index, columns=inputs_ALL['param_df']['sets']['n'])
# heat_slack_NO = pd.DataFrame(index=heat_demand_raw_NO.index, columns=inputs_ALL['param_df']['sets']['n'])
# heat_slack_HY = pd.DataFrame(index=heat_demand_raw_HY.index, columns=inputs_ALL['param_df']['sets']['n'])
#
# for c in inputs_ALL['param_df']['sets']['n']:
#     heat_demand_ALL.loc[:, c] = heat_demand_raw_ALL.loc[:, heat_demand_raw_ALL.columns.str.contains(' ' + c + '_')].sum(
#         axis=1)
#     heat_demand_EV.loc[:, c] = heat_demand_raw_EV.loc[:, heat_demand_raw_EV.columns.str.contains(' ' + c + '_')].sum(
#         axis=1)
#     heat_demand_TH.loc[:, c] = heat_demand_raw_TH.loc[:, heat_demand_raw_TH.columns.str.contains(' ' + c + '_')].sum(
#         axis=1)
#     heat_demand_NO.loc[:, c] = heat_demand_raw_NO.loc[:, heat_demand_raw_NO.columns.str.contains(' ' + c + '_')].sum(
#         axis=1)
#     heat_demand_HY.loc[:, c] = heat_demand_raw_HY.loc[:, heat_demand_raw_HY.columns.str.contains(' ' + c + '_')].sum(
#         axis=1)
#
# for c in inputs_ALL['param_df']['sets']['n']:
#     heat_slack_ALL.loc[:, c] = results_ALL['OutputHeatSlack'].loc[:,
#                                results_ALL['OutputHeatSlack'].columns.str.contains(' ' + c + '_')].sum(axis=1)
#     heat_slack_EV.loc[:, c] = results_EV['OutputHeatSlack'].loc[:,
#                               results_EV['OutputHeatSlack'].columns.str.contains(' ' + c + '_')].sum(axis=1)
#     heat_slack_TH.loc[:, c] = results_TH['OutputHeatSlack'].loc[:,
#                               results_TH['OutputHeatSlack'].columns.str.contains(' ' + c + '_')].sum(axis=1)
#     heat_slack_NO.loc[:, c] = results_NO['OutputHeatSlack'].loc[:,
#                               results_NO['OutputHeatSlack'].columns.str.contains(' ' + c + '_')].sum(axis=1)
#     heat_slack_HY.loc[:, c] = results_HY['OutputHeatSlack'].loc[:,
#                               results_HY['OutputHeatSlack'].columns.str.contains(' ' + c + '_')].sum(axis=1)
#
# heat_mix = pd.DataFrame(index=scenarios, columns=heat_prod.columns)
#
# heat_mix.loc['ALLFLEX', :] = heat_prod.loc['ALLFLEX'] / ((heat_demand_ALL - heat_slack_ALL).sum(axis=1).sum() / 1e6)
# heat_mix.loc['EVFLEX', :] = heat_prod.loc['EVFLEX'] / ((heat_demand_EV - heat_slack_EV).sum(axis=1).sum() / 1e6)
# heat_mix.loc['THFLEX', :] = heat_prod.loc['THFLEX'] / ((heat_demand_TH - heat_slack_TH).sum(axis=1).sum() / 1e6)
# heat_mix.loc['NOFLEX', :] = heat_prod.loc['NOFLEX'] / ((heat_demand_NO - heat_slack_NO).sum(axis=1).sum() / 1e6)
# heat_mix.loc['HYFLEX', :] = heat_prod.loc['HYFLEX'] / ((heat_demand_HY - heat_slack_HY).sum(axis=1).sum() / 1e6)
#
# heat_mix.drop(columns=['Backup Heater'], inplace=True)
#
# # Delta CHP Production
#
# delta_chp = (heat_prod.loc[:, heat_prod.columns.str.contains('CHP')].sum(axis=1)['ALLFLEX'] -
#              heat_prod.loc[:, heat_prod.columns.str.contains('CHP')].sum(axis=1)['NOFLEX']) / \
#             heat_prod.loc[:, heat_prod.columns.str.contains('CHP')].sum(axis=1)['NOFLEX']
#
# # %% Heat output plot
# ax = heat_prod.plot(kind='bar', stacked=True, rot=0, color=[colors_tech.get(x) for x in heat_prod.columns],
#                     legend=False, figsize=(10, 5), fontsize=15)
#
# for container, hatch in zip(ax.containers, ("///", '///', '///', '///', '///', '///', '///', '///', '///', '///')):
#     for patch in container.patches:
#         patch.set_hatch(hatch)
# # ax.axhline(y=1528, xmin=-1, xmax=1, color = 'black', linestyle='-', lw = 1)
#
# ax.set_ylabel('Heat [TWh]', fontsize=15)
# ax.set_xlabel('Scenario', fontsize=15)
# ax.set_title("Heat Output", fontsize=15)
# handles, labels = ax.get_legend_handles_labels()
#
# ax.legend(reversed(handles), reversed(labels), loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fontsize=15)
#
# fig = ax.get_figure()
# fig.savefig(output_folder + "\Heat output.png", bbox_inches='tight')



# %% Shadow Prices plots
# %% Shadow price ALLFLEX - NOFLEX
import matplotlib

shad_price_ALL = results_ALL['ShadowPrice']

for c in shad_price_ALL.columns:
    shad_price_ALL.loc[shad_price_ALL[c] >= 1e5, :] = 1e5

shad_price_ALL.sort_index(axis=1, inplace=True)

shad_price_NO = results_NO['ShadowPrice']

for c in shad_price_NO.columns:
    shad_price_NO.loc[shad_price_NO[c] >= 1e5, :] = 1e5

shad_price_NO.sort_index(axis=1, inplace=True)

fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(20, 14))

# sns.set(font_scale=1.17)
# fig, ax = plt.subplots(figsize=(10,5))
xticks_sh = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
ax1 = sns.heatmap(shad_price_NO.T, ax=axes[0], vmax=400, cmap='summer', yticklabels=True, cbar=False)
sns.heatmap(shad_price_NO.T, ax=axes[0], mask=shad_price_NO.T < 1e5, cbar=False, cmap='bwr')
ax2 = sns.heatmap(shad_price_ALL.T, ax=axes[1], cmap='summer', vmax=400, yticklabels=True, cbar=False)
sns.heatmap(shad_price_ALL.T, ax=axes[1], mask=shad_price_ALL.T < 1e5, cbar=False, cmap='bwr')

# PuOr
ax2.set_title('Shadow price for each zone - ALLFLEX', fontsize=18)
ax2.set_yticklabels(shad_price_ALL.columns, rotation=0)
ax2.set_xticks(range(1, 8784, 798))
ax2.set_xticklabels(xticks_sh, rotation=0)
ax2.tick_params(labelsize=15)

ax1.set_title('Shadow price for each zone - NOFLEX', fontsize=18)
ax1.set_yticklabels(shad_price_ALL.columns, rotation=0)
ax1.set_xticks(range(1, 8784, 798))
ax1.set_xticklabels(xticks_sh, rotation=0)
ax1.tick_params(labelsize=15)

mappable = ax1.get_children()[0]
cbar = plt.colorbar(mappable, ax=[axes[0], axes[1]], orientation='vertical')
cbar.ax.tick_params(labelsize=15)

plt.show()

fig = ax2.get_figure()
fig.savefig(output_folder + "\Shadow price ALLFLEX.png", bbox_inches='tight')

# %% Shadow price duration curve

# TODO If needed do the 28 plots where all the 5 scenarios are shown at the same time.

shad_price_ALL = results_ALL['ShadowPrice']

el.plot_LDC(shad_price_ALL, stacked=False, x_norm=False, legend=True)
plt.ylabel('Daily Generation Cost Electricity [M€]')
plt.xlabel('Hour')
plt.legend()

# %% Shadow price ALLFLEX

dates = {}

dates[0] = ['2016-01-18', '2016-01-24']
dates[1] = ['2016-02-22', '2016-02-28']
dates[2] = ['2016-06-05', '2016-06-11']
dates[3] = ['2016-07-10', '2016-07-16']
dates[4] = ['2016-08-07', '2016-08-13']

xticks_day = {}

xticks_day[0] = ['18\nJan', '19', '20', '21', '22', '23', '24']
xticks_day[1] = ['22\nFeb', '23', '24', '25', '26', '27', '28']
xticks_day[2] = ['5\nJun', '6', '7', '8', '9', '10', '11']
xticks_day[3] = ['10\nJul', '11', '12', '13', '14', '15', '16']
xticks_day[4] = ['7\nAug', '8', '9', '10', '11', '12', '13']

i = 0
dates = dates[i]
xticks_day = xticks_day[i]

shad_price_ALL = results_ALL['ShadowPrice']
shad_price_ALL = shad_price_ALL.loc[dates[0]:dates[1], :]

shad_price_ALL.sort_index(axis=1, inplace=True)

sns.set(font_scale=1.27)
# xticks_sh = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ax = sns.heatmap(shad_price_ALL.T, vmax=1000, cmap='coolwarm')
sns.set(rc={'figure.figsize': (8, 8)})
ax.set_title('Shadow price for each zone - ALLFLEX', fontsize=15)
ax.set_yticklabels(shad_price_ALL.columns)
# ax.set_xticks(range(1,8784,798))
# ax.set_xticklabels(xticks_sh, rotation = 45)
ax.set_xticks(range(12, 157, 24))
ax.set_xticklabels(xticks_day, rotation=45)

fig = ax.get_figure()
fig.savefig(output_folder + "\Shadow price ALLFLEX " + str(dates) + " .png", bbox_inches='tight')

# %% Shadow price EVFLEX

scenario = 'EVFLEX'
i = 2

dates = {}

dates[0] = ['2016-01-18', '2016-01-24']
dates[1] = ['2016-02-22', '2016-02-28']
dates[2] = ['2016-06-05', '2016-06-11']
dates[3] = ['2016-07-10', '2016-07-16']
dates[4] = ['2016-08-07', '2016-08-13']

xticks_day = {}

xticks_day[0] = ['18\nJan', '19', '20', '21', '22', '23', '24']
xticks_day[1] = ['22\nFeb', '23', '24', '25', '26', '27', '28']
xticks_day[2] = ['5\nJun', '6', '7', '8', '9', '10', '11']
xticks_day[3] = ['10\nJul', '11', '12', '13', '14', '15', '16']
xticks_day[4] = ['7\nAug', '8', '9', '10', '11', '12', '13']

dates = dates[i]
xticks_day = xticks_day[i]

shad_price = {}

shad_price[scenario] = results[scenario]['ShadowPrice']
shad_price[scenario] = shad_price[scenario].loc[dates[0]:dates[1], :]

shad_price[scenario].sort_index(axis=1, inplace=True)

sns.set(font_scale=1.27)
# xticks_sh = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ax = sns.heatmap(shad_price[scenario].T, vmax=1000, cmap='coolwarm')
sns.set(rc={'figure.figsize': (8, 8)})
ax.set_title('Shadow price for each zone - ' + scenario, fontsize=15)
ax.set_yticklabels(shad_price[scenario].columns)
# ax.set_xticks(range(1,8784,798))
# ax.set_xticklabels(xticks_sh, rotation = 45)
ax.set_xticks(range(12, 157, 24))
ax.set_xticklabels(xticks_day, rotation=45)

fig = ax.get_figure()
fig.savefig(output_folder + "\Shadow price " + scenario + " " + str(dates) + " .png", bbox_inches='tight')

# %% Shadow price THFLEX

scenario = 'THFLEX'
i = 4

dates = {}

dates[0] = ['2016-01-18', '2016-01-24']
dates[1] = ['2016-02-22', '2016-02-28']
dates[2] = ['2016-06-05', '2016-06-11']
dates[3] = ['2016-07-10', '2016-07-16']
dates[4] = ['2016-08-07', '2016-08-13']

xticks_day = {}

xticks_day[0] = ['18\nJan', '19', '20', '21', '22', '23', '24']
xticks_day[1] = ['22\nFeb', '23', '24', '25', '26', '27', '28']
xticks_day[2] = ['5\nJun', '6', '7', '8', '9', '10', '11']
xticks_day[3] = ['10\nJul', '11', '12', '13', '14', '15', '16']
xticks_day[4] = ['7\nAug', '8', '9', '10', '11', '12', '13']

dates = dates[i]
xticks_day = xticks_day[i]

shad_price = {}

shad_price[scenario] = results[scenario]['ShadowPrice']
shad_price[scenario] = shad_price[scenario].loc[dates[0]:dates[1], :]

shad_price[scenario].sort_index(axis=1, inplace=True)

sns.set(font_scale=1.27)
# xticks_sh = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ax = sns.heatmap(shad_price[scenario].T, vmax=1000, cmap='coolwarm')
sns.set(rc={'figure.figsize': (8, 8)})
ax.set_title('Shadow price for each zone - ' + scenario, fontsize=15)
ax.set_yticklabels(shad_price[scenario].columns)
# ax.set_xticks(range(1,8784,798))
# ax.set_xticklabels(xticks_sh, rotation = 45)
ax.set_xticks(range(12, 157, 24))
ax.set_xticklabels(xticks_day, rotation=45)

fig = ax.get_figure()
fig.savefig(output_folder + "\Shadow price " + scenario + " " + str(dates) + " .png", bbox_inches='tight')

# %% Shadow price NOFLEX

scenario = 'NOFLEX'
i = 4

dates = {}

dates[0] = ['2016-01-18', '2016-01-24']
dates[1] = ['2016-02-22', '2016-02-28']
dates[2] = ['2016-06-05', '2016-06-11']
dates[3] = ['2016-07-10', '2016-07-16']
dates[4] = ['2016-08-07', '2016-08-13']

xticks_day = {}

xticks_day[0] = ['18\nJan', '19', '20', '21', '22', '23', '24']
xticks_day[1] = ['22\nFeb', '23', '24', '25', '26', '27', '28']
xticks_day[2] = ['5\nJun', '6', '7', '8', '9', '10', '11']
xticks_day[3] = ['10\nJul', '11', '12', '13', '14', '15', '16']
xticks_day[4] = ['7\nAug', '8', '9', '10', '11', '12', '13']

dates = dates[i]
xticks_day = xticks_day[i]

shad_price = {}

shad_price[scenario] = results[scenario]['ShadowPrice']
shad_price[scenario] = shad_price[scenario].loc[dates[0]:dates[1], :]

shad_price[scenario].sort_index(axis=1, inplace=True)

sns.set(font_scale=1.27)
# xticks_sh = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ax = sns.heatmap(shad_price[scenario].T, vmax=1000, cmap='coolwarm')
sns.set(rc={'figure.figsize': (8, 8)})
ax.set_title('Shadow price for each zone - ' + scenario, fontsize=15)
ax.set_yticklabels(shad_price[scenario].columns)
# ax.set_xticks(range(1,8784,798))
# ax.set_xticklabels(xticks_sh, rotation = 45)
ax.set_xticks(range(12, 157, 24))
ax.set_xticklabels(xticks_day, rotation=45)

fig = ax.get_figure()
fig.savefig(output_folder + "\Shadow price " + scenario + " " + str(dates) + " .png", bbox_inches='tight')

# %% Shadow price HYFLEX

scenario = 'HYFLEX'
i = 4

dates = {}

dates[0] = ['2016-01-18', '2016-01-24']
dates[1] = ['2016-02-22', '2016-02-28']
dates[2] = ['2016-06-05', '2016-06-11']
dates[3] = ['2016-07-10', '2016-07-16']
dates[4] = ['2016-08-07', '2016-08-13']

xticks_day = {}

xticks_day[0] = ['18\nJan', '19', '20', '21', '22', '23', '24']
xticks_day[1] = ['22\nFeb', '23', '24', '25', '26', '27', '28']
xticks_day[2] = ['5\nJun', '6', '7', '8', '9', '10', '11']
xticks_day[3] = ['10\nJul', '11', '12', '13', '14', '15', '16']
xticks_day[4] = ['7\nAug', '8', '9', '10', '11', '12', '13']

dates = dates[i]
xticks_day = xticks_day[i]

shad_price = {}

shad_price[scenario] = results[scenario]['ShadowPrice']
shad_price[scenario] = shad_price[scenario].loc[dates[0]:dates[1], :]

shad_price[scenario].sort_index(axis=1, inplace=True)

sns.set(font_scale=1.27)
# xticks_sh = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ax = sns.heatmap(shad_price[scenario].T, vmax=1000, cmap='coolwarm')
sns.set(rc={'figure.figsize': (8, 8)})
ax.set_title('Shadow price for each zone - ' + scenario, fontsize=15)
ax.set_yticklabels(shad_price[scenario].columns)
# ax.set_xticks(range(1,8784,798))
# ax.set_xticklabels(xticks_sh, rotation = 45)
ax.set_xticks(range(12, 157, 24))
ax.set_xticklabels(xticks_day, rotation=45)

fig = ax.get_figure()
fig.savefig(output_folder + "\Shadow price " + scenario + " " + str(dates) + " .png", bbox_inches='tight')

# %% Shadow price EU Average - 5 scenarios

plt.style.use('default')

z = 'BE'  # 'ALL' to do the EU average

dates = {}

dates[0] = ['2016-01-18', '2016-01-24']
dates[1] = ['2016-02-22', '2016-02-28']
dates[2] = ['2016-06-05', '2016-06-11']
dates[3] = ['2016-07-10', '2016-07-16']
dates[4] = ['2016-08-07', '2016-08-13']

# dates[0] = ['2016-01-01','2016-12-31']
# dates[1] = ['2016-02-01','2016-03-01']
# dates[2] = ['2016-06-01','2016-07-01']
# dates[3] = ['2016-07-01','2016-08-01']
# dates[4] = ['2016-08-01','2016-09-01']

xticks_day = {}

xticks_day[0] = ['18\nJan', '18', '20', '21', '22', '23', '24']
xticks_day[1] = ['22\nFeb', '23', '24', '25', '26', '27', '28']
xticks_day[2] = ['5\nJun', '6', '7', '8', '9', '10', '11']
xticks_day[3] = ['10\nJul', '11', '12', '13', '14', '15', '16']
xticks_day[4] = ['7\nAug', '8', '9', '10', '11', '12', '13']

if z == 'ALL':
    for i in list(dates.keys()):

        shad_price_av = pd.DataFrame(index=results['ALLFLEX']['ShadowPrice'].index, columns=scenarios)

        for scenario in scenarios:
            shad_price_av.loc[:, scenario] = results[scenario]['ShadowPrice'].mean(axis=1)

        shad_price_av = shad_price_av.loc[dates[i][0]:dates[i][1], :]

        shad_price_av_plot = pd.DataFrame(shad_price_av.values, columns=shad_price_av.columns)

        ax = shad_price_av_plot.plot(kind='line', rot=0, fontsize=15, ylim=(0, 400), legend=True)
        ax.set_xticks(range(12, 157, 24))
        ax.set_xticklabels(xticks_day[i])
        ax.set_title("Shadow Price EU Average", fontsize=15)

        fig = ax.get_figure()
        fig.savefig(output_folder + "\Shadow price Average " + str(dates[i]) + " .png", bbox_inches='tight')

else:
    for i in list(dates.keys()):

        shad_price = pd.DataFrame(index=results['ALLFLEX']['ShadowPrice'].index, columns=scenarios)

        for scenario in scenarios:
            shad_price.loc[:, scenario] = results[scenario]['ShadowPrice'][z]

        shad_price = shad_price.loc[dates[i][0]:dates[i][1], :]

        shad_price_plot = pd.DataFrame(shad_price.values, columns=shad_price.columns)

        ax = shad_price_plot.plot(kind='line', rot=0, fontsize=15, ylim=(0, 400), legend=True)
        ax.set_xticks(range(12, 157, 24))
        ax.set_xticklabels(xticks_day[i])
        ax.set_title("Shadow Price " + z, fontsize=15)

        fig = ax.get_figure()
        fig.savefig(output_folder + "\Shadow price " + z + ' ' + str(dates[i]) + ".png", bbox_inches='tight')

# %% Heat Shadow Price EU Average

plt.style.use('default')

dates = {}

dates[0] = ['2016-01-18', '2016-01-24']
dates[1] = ['2016-02-22', '2016-02-28']
dates[2] = ['2016-06-05', '2016-06-11']
dates[3] = ['2016-07-10', '2016-07-16']
dates[4] = ['2016-08-07', '2016-08-13']

# dates[0] = ['2016-01-01','2016-12-31']
# dates[1] = ['2016-02-01','2016-03-01']
# dates[2] = ['2016-06-01','2016-07-01']
# dates[3] = ['2016-07-01','2016-08-01']
# dates[4] = ['2016-08-01','2016-09-01']

xticks_day = {}

xticks_day[0] = ['18\nJan', '18', '20', '21', '22', '23', '24']
xticks_day[1] = ['22\nFeb', '23', '24', '25', '26', '27', '28']
xticks_day[2] = ['5\nJun', '6', '7', '8', '9', '10', '11']
xticks_day[3] = ['10\nJul', '11', '12', '13', '14', '15', '16']
xticks_day[4] = ['7\nAug', '8', '9', '10', '11', '12', '13']

yr = 2016
hour = pd.date_range(start=str(yr) + '-01-01', end=str(yr) + '-12-31 23:00', freq='H')

for i in list(dates.keys()):

    heat_shad_price_av = pd.DataFrame(index=results['ALLFLEX']['HeatShadowPrice'].index, columns=scenarios)
    heat_shad_price_av.set_index(hour, inplace=True)

    for scenario in scenarios:
        for c in results[scenario]['HeatShadowPrice'].columns:
            if results[scenario]['HeatShadowPrice'].loc[:, c].sum() > 1e+100:
                for i in results[scenario]['HeatShadowPrice'].index:
                    if results[scenario]['HeatShadowPrice'].loc[i, c] > 1e+100:
                        results[scenario]['HeatShadowPrice'].loc[i, c] = 0

        heat_shad_price_av.loc[:, scenario] = results[scenario]['HeatShadowPrice'].mean(axis=1).values

    heat_shad_price_av = heat_shad_price_av.loc[dates[i][0]:dates[i][1], :]

    heat_shad_price_av_plot = pd.DataFrame(heat_shad_price_av.values, columns=heat_shad_price_av.columns)

    ax = heat_shad_price_av_plot.plot(kind='line', rot=0, fontsize=15, legend=True)
    # ax.set_xticks(range(12,157,24))
    # ax.set_xticklabels(xticks_day[i])
    ax.set_title("Heat Shadow Price EU Average", fontsize=15)

    fig = ax.get_figure()
    fig.savefig(output_folder + "\Heat Shadow price Average " + str(dates[i]) + " .png", bbox_inches='tight')

# %% Shadow price % difference

shad_price_ALL_av = results_ALL['ShadowPrice'].sum(axis=1).mean()
shad_price_EV_av = results_EV['ShadowPrice'].sum(axis=1).mean()
shad_price_TH_av = results_TH['ShadowPrice'].sum(axis=1).mean()
shad_price_NO_av = results_NO['ShadowPrice'].sum(axis=1).mean()
shad_price_HY_av = results_HY['ShadowPrice'].sum(axis=1).mean()

diff_per_HY = (shad_price_HY_av - shad_price_NO_av) / shad_price_NO_av
diff_per_EV = (shad_price_EV_av - shad_price_NO_av) / shad_price_NO_av
diff_per_TH = (shad_price_TH_av - shad_price_NO_av) / shad_price_NO_av
diff_per_ALL = (shad_price_ALL_av - shad_price_NO_av) / shad_price_NO_av

diff_per_EV_hy = (diff_per_EV - diff_per_HY) / diff_per_HY
diff_per_TH_hy = (diff_per_TH - diff_per_HY) / diff_per_HY



# %% Lost Load power

ll_list = ['LostLoad_2D', 'LostLoad_2U', 'LostLoad_3U', 'LostLoad_MaxPower', 'LostLoad_MinPower', 'LostLoad_RampDown',
           'LostLoad_RampUp', 'LostLoad_WaterSlack']
ll_short = ['2D', '2U', '3U', 'MaxPower', 'MinPower', 'RampDown', 'RampUp', 'WaterSlack']

lost_load_raw_ALL = pd.DataFrame(columns=['ALLFLEX'], index=ll_list).T
lost_load_raw_EV = pd.DataFrame(columns=['EVFLEX'], index=ll_list).T
lost_load_raw_TH = pd.DataFrame(columns=['THFLEX'], index=ll_list).T
lost_load_raw_NO = pd.DataFrame(columns=['NOFLEX'], index=ll_list).T
lost_load_raw_HY = pd.DataFrame(columns=['HYFLEX'], index=ll_list).T

for ll in ll_list:
    if 'Water' not in ll:
        lost_load_raw_ALL[ll] = results_ALL[ll].sum().sum()
        lost_load_raw_EV[ll] = results_EV[ll].sum().sum()
        lost_load_raw_TH[ll] = results_TH[ll].sum().sum()
        lost_load_raw_NO[ll] = results_NO[ll].sum().sum()
        lost_load_raw_HY[ll] = results_HY[ll].sum().sum()
    else:
        lost_load_raw_ALL[ll] = results_ALL[ll].sum()
        lost_load_raw_EV[ll] = results_EV[ll]
        lost_load_raw_TH[ll] = results_TH[ll]
        lost_load_raw_NO[ll] = results_NO[ll]
        lost_load_raw_HY[ll] = results_HY[ll]

lost_load = pd.DataFrame(index=scenarios, columns=lost_load_raw_ALL.columns)

lost_load.loc['ALLFLEX', :] = lost_load_raw_ALL.values / 1e6
lost_load.loc['EVFLEX', :] = lost_load_raw_EV.values / 1e6
lost_load.loc['THFLEX', :] = lost_load_raw_TH.values / 1e6
lost_load.loc['NOFLEX', :] = lost_load_raw_NO.values / 1e6
lost_load.loc['HYFLEX', :] = lost_load_raw_HY.values / 1e6

d_rename = dict(zip(ll_list, ll_short))
lost_load.rename(columns=d_rename, inplace=True)

ll_sum_ALL = results_ALL['LostLoad_2D'].add(results_ALL['LostLoad_2U'], fill_value=0).add(results_ALL['LostLoad_3U'],
                                                                                          fill_value=0).add(
    results_ALL['LostLoad_MaxPower'], fill_value=0).add(results_ALL['LostLoad_MinPower'], fill_value=0).add(
    results_ALL['LostLoad_RampDown'], fill_value=0).add(results_ALL['LostLoad_RampUp'], fill_value=0)
ll_max_ALL = ll_sum_ALL.sum(axis=1).max()

ll_sum_EV = results_EV['LostLoad_2D'].add(results_EV['LostLoad_2U'], fill_value=0).add(results_EV['LostLoad_3U'],
                                                                                       fill_value=0).add(
    results_EV['LostLoad_MaxPower'], fill_value=0).add(results_EV['LostLoad_MinPower'], fill_value=0).add(
    results_EV['LostLoad_RampDown'], fill_value=0).add(results_EV['LostLoad_RampUp'], fill_value=0)
ll_max_EV = ll_sum_EV.sum(axis=1).max()

ll_sum_TH = results_TH['LostLoad_2D'].add(results_TH['LostLoad_2U'], fill_value=0).add(results_TH['LostLoad_3U'],
                                                                                       fill_value=0).add(
    results_TH['LostLoad_MaxPower'], fill_value=0).add(results_TH['LostLoad_MinPower'], fill_value=0).add(
    results_TH['LostLoad_RampDown'], fill_value=0).add(results_TH['LostLoad_RampUp'], fill_value=0)
ll_max_TH = ll_sum_TH.sum(axis=1).max()

ll_sum_NO = results_NO['LostLoad_2D'].add(results_NO['LostLoad_2U'], fill_value=0).add(results_NO['LostLoad_3U'],
                                                                                       fill_value=0).add(
    results_NO['LostLoad_MaxPower'], fill_value=0).add(results_NO['LostLoad_MinPower'], fill_value=0).add(
    results_NO['LostLoad_RampDown'], fill_value=0).add(results_NO['LostLoad_RampUp'], fill_value=0)
ll_max_NO = ll_sum_NO.sum(axis=1).max()

ll_sum_HY = results_HY['LostLoad_2D'].add(results_HY['LostLoad_2U'], fill_value=0).add(results_HY['LostLoad_3U'],
                                                                                       fill_value=0).add(
    results_HY['LostLoad_MaxPower'], fill_value=0).add(results_HY['LostLoad_MinPower'], fill_value=0).add(
    results_HY['LostLoad_RampDown'], fill_value=0).add(results_HY['LostLoad_RampUp'], fill_value=0)
ll_max_HY = ll_sum_HY.sum(axis=1).max()

lost_load_max = pd.DataFrame(index=scenarios, columns=['Lost Load Max'])

lost_load_max.loc['ALLFLEX', :] = ll_max_ALL / 1e3
lost_load_max.loc['EVFLEX', :] = ll_max_EV / 1e3
lost_load_max.loc['THFLEX', :] = ll_max_TH / 1e3
lost_load_max.loc['NOFLEX', :] = ll_max_NO / 1e3
lost_load_max.loc['HYFLEX', :] = ll_max_HY / 1e3

# %% Lost Load power plot

ax = lost_load.plot(kind='bar', stacked=True, color=[colors_tech.get(x) for x in lost_load.columns], width=0.3, rot=0,
                    figsize=(7, 5), position=1, legend=True, fontsize=15)

# ax = lost_load.plot(kind='bar', stacked=True, rot = 0,color=[colors_tech.get(x) for x in lost_load.columns], legend = False,  figsize = (10,5), fontsize=15)
ax2 = lost_load_max.plot(kind='bar', secondary_y=True, rot=0, ax=ax, position=0, color='peachpuff', width=0.3,
                         legend=True, mark_right=False, fontsize=15)

ax.set_ylabel('Lost Load [TWh]', fontsize=15)
ax2.set_ylabel('Max Lost Load [GW]', fontsize=15)
ax.set_xlabel('Scenario', fontsize=15)
ax2.set_title("Lost Load", fontsize=15)

handles, labels = ax.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()
handles.append(handles2[0])
labels.append(labels2[0])
ax.legend(reversed(handles), reversed(labels))

fig = ax2.get_figure()
fig.savefig(output_folder + "\Lost Load.png", bbox_inches='tight')

# %% Lost Load hours

ll_hours = pd.DataFrame(index=scenarios, columns=['2D', '2U', '3U', 'MaxPower', 'MinPower', 'RampDown', 'RampUp'])


def ll_hours_fill(results, scenario, ll_type):
    ll_hours.loc[scenario, ll_type] = results['LostLoad_' + ll_type].sum(axis=1)[
        results['LostLoad_' + ll_type].sum(axis=1) != 0].count()
    return True


for l in ll_hours.columns:
    ll_hours_fill(results_ALL, 'ALLFLEX', l)
    ll_hours_fill(results_EV, 'EVFLEX', l)
    ll_hours_fill(results_TH, 'THFLEX', l)
    ll_hours_fill(results_NO, 'NOFLEX', l)
    ll_hours_fill(results_HY, 'HYFLEX', l)

# ll_2U_ALL_f = results_ALL['LostLoad_2U']
#
# for c in ll_2U_ALL_f.columns:
#    if ll_2U_ALL_f[c].values.sum() == 0:
#        ll_2U_ALL_f.drop(columns = [c], inplace = True)


# %% Lost load hours plot

# ax = shload_tot.T.plot(kind='bar', color='orangered', width=0.3, rot = 0, figsize = (7,5), position=1, legend  = True, fontsize=15)
# ax2 = shload_max.T.plot(kind='bar', color='firebrick', secondary_y=True, rot= 0, ax=ax, width=0.3, position=0, legend  = True, mark_right = False, fontsize=15)
fig = plt.figure()  # Create matplotlib figure

# ax = fig.add_subplot(111)
# width = 0.1
#
# pos = np.linspace(0,2.5, 4)
ax = ll_hours.plot(kind='bar', rot=0, width=1.2, legend=True, figsize=(10, 5), fontsize=15)

# ax1 = ll_hours.iloc[:,0].plot(kind='bar', rot = 0, width = width, color = 'navy', legend = True,  position=pos[3], figsize = (10,5), fontsize=15)
# ax2 = ll_hours.iloc[:,1].plot(kind='bar', rot = 0,  width = width,color = 'dodgerblue', legend = True,  position=pos[2], figsize = (10,5), fontsize=15)
# ax3 = ll_hours.iloc[:,2].plot(kind='bar', rot = 0,  width = width,color = 'aqua', legend = True,  position=pos[1], figsize = (10,5), fontsize=15)
# ax4 = ll_hours.iloc[:,3].plot(kind='bar', rot = 0, width = width,color = 'springgreen', legend = True,  position=pos[0], figsize = (10,5), fontsize=15)

# for container, hatch in zip(ax.containers, ( '///','///','///','///','///','///','///','///')):
#    for patch in container.patches:
#        patch.set_hatch(hatch)

ax.set_ylabel('Lost Load [h]', fontsize=15)
ax.set_xlabel('Scenario', fontsize=15)
ax.set_title("Lost Load hours", fontsize=15)
# handles, labels = ax.get_legend_handles_labels()
# ax1.legend(reversed(handles), reversed(labels))
# , loc='center left', bbox_to_anchor=(1, 0.5), ncol=2, fontsize=15)

fig = ax.get_figure()
fig.savefig(output_folder + "\Lost Load hours.png", bbox_inches='tight')


# %% Average cost of el without lost load

def ll_hours_index(results):
    if results['LostLoad_MaxPower'].values.sum() > 0:
        ll_hours_index = results['LostLoad_MaxPower'].index
    else:
        ll_hours_index = []
    return ll_hours_index


ll_hours_index_ALL = ll_hours_index(results_ALL)
ll_hours_index_EV = ll_hours_index(results_EV)
ll_hours_index_TH = ll_hours_index(results_TH)
ll_hours_index_NO = ll_hours_index(results_NO)
ll_hours_index_HY = ll_hours_index(results_HY)


def get_demand(inputs, results):
    demand = {}
    for z in inputs['sets']['n']:
        if 'OutputPowerConsumption' in results:
            demand_p2h = ds.filter_by_zone(results['OutputPowerConsumption'], inputs, z)
            demand_p2h = demand_p2h.sum(axis=1)
        else:
            demand_p2h = pd.Series(0, index=results['OutputPower'].index)
        demand_da = inputs['param_df']['Demand'][('DA', z)]
        demand[z] = pd.DataFrame(demand_da + demand_p2h, columns=[('DA', z)])
    demand = pd.concat(demand, axis=1)
    demand.columns = demand.columns.droplevel(-1)

    return demand


demand_f_ALL = get_demand(inputs_ALL, results_ALL).drop(list(ll_hours_index_ALL))
demand_f_EV = get_demand(inputs_EV, results_EV).drop(list(ll_hours_index_EV))
demand_f_TH = get_demand(inputs_TH, results_TH).drop(list(ll_hours_index_TH))
demand_f_NO = get_demand(inputs_NO, results_NO).drop(list(ll_hours_index_NO))
demand_f_HY = get_demand(inputs_HY, results_HY).drop(list(ll_hours_index_HY))

total_load_f_ALL = demand_f_ALL.sum().sum()
total_load_f_EV = demand_f_EV.sum().sum()
total_load_f_TH = demand_f_TH.sum().sum()
total_load_f_NO = demand_f_NO.sum().sum()
total_load_f_HY = demand_f_HY.sum().sum()

cost_f_ALL = results_ALL['OutputSystemCost'].drop(list(ll_hours_index_ALL)).sum()
cost_f_EV = results_EV['OutputSystemCost'].drop(list(ll_hours_index_EV)).sum()
cost_f_TH = results_TH['OutputSystemCost'].drop(list(ll_hours_index_TH)).sum()
cost_f_NO = results_NO['OutputSystemCost'].drop(list(ll_hours_index_NO)).sum()
cost_f_HY = results_HY['OutputSystemCost'].drop(list(ll_hours_index_HY)).sum()

cost_av_f_ALL = cost_f_ALL / total_load_f_ALL
cost_av_f_EV = cost_f_EV / total_load_f_EV
cost_av_f_TH = cost_f_NO / total_load_f_TH
cost_av_f_NO = cost_f_NO / total_load_f_NO
cost_av_f_HY = cost_f_HY / total_load_f_HY

# %%
''' Demand Breakdown '''

# %% EV demand calculation

#### Input demands
import pandas as pd

inputfile_el2050 = r"C:\Users\Andrea\OneDrive - Politecnico di Milano\Università\Tesi (OneDrive)\Article PROres1 coupling\Raw Data\JRC-EU-TIMES PROres1\PROres raw\Electricity generation 2050.csv"
inputfile_el2020 = r"C:\Users\Andrea\OneDrive - Politecnico di Milano\Università\Tesi (OneDrive)\Article PROres1 coupling\Raw Data\JRC-EU-TIMES PROres1\PROres raw\Electricity generation 2020.csv"
inputfile_heat_diff = r"C:\Users\Andrea\OneDrive - Politecnico di Milano\Università\Tesi (OneDrive)\Article PROres1 coupling\Raw Data\JRC-EU-TIMES PROres1\PROres Postprocessed\delta el_heat 2050-20.csv"

# PJ
el_2050 = pd.read_csv(inputfile_el2050, header=0, index_col=0, skiprows=1)
el_2020 = pd.read_csv(inputfile_el2020, header=0, index_col=0, skiprows=1)
heat_delta = pd.read_csv(inputfile_heat_diff, header=0, index_col=0, skiprows=1)
heat_delta = heat_delta.sum(axis=1)

el_2050 = el_2050.iloc[:, 0] / 3.6
el_2020 = el_2020.iloc[:, 0] / 3.6
heat_delta = heat_delta / 3.6

country_coeff = (el_2050 - heat_delta) / el_2020
country_coeff.drop(index=['CY', 'Mt'], inplace=True)

##### Import the power curves from dispaset
inputfolder = r"C:\Users\Andrea\OneDrive - Politecnico di Milano\Università\Tesi (OneDrive)\Article PROres1 coupling\Database\TotalLoadValue"
countries = list(country_coeff.index)

# import power curves from dispaset database
p_curve_dict = {}
for c in countries:
    inputfile = inputfolder + "/" + c + r"\1h\2016.csv"
    p_curve_dict[c] = pd.read_csv(inputfile, header=None, index_col=0)
    p_curve_dict[c] = p_curve_dict[c].iloc[:, 0]

# create dataframe from
p_curve = pd.DataFrame.from_dict(p_curve_dict)

p_curve_scaled = p_curve * country_coeff.T

delta_el_curve = p_curve_scaled - p_curve

#### Add electric demand

inputfile = r"C:\Users\Andrea\OneDrive - Politecnico di Milano\Università\Tesi (OneDrive)\Database Improvements\3. EV profiles\Other approaches\PyPSA model\EV charging demand\Pypsa EV charging demand.csv"
ev_demand = pd.read_csv(inputfile, header=0, index_col=0)

ev_demand_ad = ev_demand / (ev_demand.sum(axis=0) / 10 ** 6)
ev_demand_ad['UK'] = ev_demand_ad.pop('GB')

ev_demand_ad['CY'] = ev_demand_ad['EL']
ev_demand_ad['MT'] = ev_demand_ad['EL']

ev_demand_ad.drop(columns=['BA', 'RS'], inplace=True)

inputfile = r"C:\Users\Andrea\OneDrive - Politecnico di Milano\Università\Tesi (OneDrive)\Article PROres1 coupling\Raw Data\JRC-EU-TIMES PROres1\PROres Postprocessed\EV Demand.xlsx"
ev_demand_times = pd.read_excel(inputfile, header=0, index_col=0)  # TWh

ev_demand_scaled = ev_demand_ad * ev_demand_times.iloc[:, 0]

for h in ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17',
          '18', '19', '20', '21', '22', '23']:
    ev_demand_scaled.loc['2011-02-29 ' + h + ':00:00', :] = (ev_demand_scaled.loc['2011-02-28 ' + h + ':00:00',
                                                             :] + ev_demand_scaled.loc['2011-03-01 ' + h + ':00:00',
                                                                  :]) / 2
    ev_demand_scaled.sort_index(inplace=True)

# Shift 2011 to 2016:
# 01-01-2011 is saturday - 01-01-2011 is friday
# so 31-12-2011 (friday) must become the first day of the year, and then I can reset the index to 2016

temp_friday = ev_demand_scaled.loc['2011-12-31 00:00:00':'2011-12-31 23:00:00', :]
ev_demand_scaled = ev_demand_scaled.loc[:'2011-12-31', :]
ev_demand_scaled = pd.concat([temp_friday, ev_demand_scaled], axis=0, join='inner')

yr = 2016
hour = pd.date_range(start=str(yr) + '-01-01', end=str(yr) + '-12-31 23:00', freq='H')
ev_demand_scaled.set_index(hour, inplace=True)

p_curve_scaled_ev = p_curve_scaled + ev_demand_scaled
p_curve_scaled_ev.drop(columns=['CY', 'MT'], inplace=True)

###### Dispatch plot of loads data #######

z = 'DE'

rng = pd.date_range(start='2016-02-15', end='2016-02-22', freq='h')

demand_p2h = results_ALL['OutputPowerConsumption'].loc[:, results_ALL['OutputPowerConsumption'].columns.str.contains(z)]
ev_demand = ev_demand_scaled[z]
heat_demand_raw = inputs_ALL['param_df']['HeatDemand']

ev_demand_plot = ev_demand.loc[rng]
ev_demand_plot = ev_demand_plot / ev_demand_plot.max()

heat_demand = pd.DataFrame(index=heat_demand_raw.index, columns=inputs_ALL['param_df']['sets']['n'])

for c in inputs_ALL['param_df']['sets']['n']:
    heat_demand.loc[:, c] = heat_demand_raw.loc[:, heat_demand_raw.columns.str.contains(' ' + c + '_')].sum(axis=1)

heat_demand_f = pd.DataFrame(index=heat_demand.index, columns=['Heat Demand'])
heat_demand_f['Heat Demand'] = heat_demand.loc[:, z].values / 1e3
heat_demand_f = heat_demand_f.loc[rng, :]

plotdata = pd.DataFrame(index=heat_demand.index, columns=['Demand 2016', 'TIMES Increase', 'EV Demand'])

plotdata.loc[:, 'Demand 2016'] = p_curve[z].values / 1e3
plotdata.loc[:, 'TIMES Increase'] = delta_el_curve[z].values / 1e3
# plotdata.loc[:,'Power Consumption'] = demand_p2h.values
plotdata.loc[:, 'EV Demand'] = ev_demand.values / 1e3
# plotdata.loc[:,'Heat Demand'] = -heat_demand[z].values

plotdata = plotdata.loc[rng, :]

p_curve_scaled_ev_f = pd.DataFrame(index=heat_demand.index, columns=['Demand 2050'])
p_curve_scaled_ev_f.loc[:, 'Demand 2050'] = p_curve_scaled_ev[z].values / 1e3
p_curve_scaled_ev_f = p_curve_scaled_ev_f.loc[rng, :]

######## Reservoir levels comparison

sto_cap_HY = pd.DataFrame(inputs_HY['units']['StorageCapacity'][inputs_HY['units']['StorageCapacity'] > 0].values,
                          columns=["HYFLEX"], index=list(
        inputs_HY['units']['StorageCapacity'][inputs_HY['units']['StorageCapacity'] > 0].index)).T

mts_wat = inputs_HY['param_df']['StorageProfile'].loc[:,
          inputs_HY['param_df']['StorageProfile'].columns.str.contains('HDAM')]

year = 2016
dates = pd.DataFrame(pd.date_range(start=str(year) + '-01-01', end=str(year) + '-12-31 23:00', freq='H'),
                     columns=['dates'])

profile_old = {}
profile_old_wat_dict = {}

for c in inputs_ALL['param_df']['sets']['n']:
    profile_old[c] = pd.read_csv(
        r"C:\Users\Andrea\OneDrive - Politecnico di Milano\Università\Tesi (OneDrive)\Article PROres1 coupling\Database\HydroData\ScaledLevels/" + c + r"\1h\ScaledLevels_ProRes1_2050.csv",
        index_col=0)
    profile_old[c].set_index(dates['dates'], inplace=True)
    profile_old_wat_dict[c] = profile_old[c].iloc[:, 0]

profile_old_wat = pd.DataFrame.from_dict(profile_old_wat_dict)

for c in profile_old_wat.columns:
    if profile_old_wat[c].values.sum() == 0:
        profile_old_wat.drop(columns=[c], inplace=True)

columns_mts_list = list(mts_wat.columns)
columns_countries_list = list(profile_old_wat.columns)

d_rename = dict(zip(columns_mts_list, columns_countries_list))
mts_wat.rename(columns=d_rename, inplace=True)

sto_cap_HY_tot_ext = pd.DataFrame(1, index=dates.iloc[:, 0], columns=sto_cap_HY.columns)
sto_cap_HY_tot_ext = sto_cap_HY.values * sto_cap_HY_tot_ext

sto_cap_HY_tot_ext.rename(columns=d_rename, inplace=True)

profile_old_wat_extensive = sto_cap_HY_tot_ext * profile_old_wat
profile_old_wat_extensive = profile_old_wat_extensive.sum(axis=1)

mts_wat_extensive = sto_cap_HY_tot_ext * mts_wat
mts_wat_extensive = mts_wat_extensive.sum(axis=1)

levels_comparison = pd.DataFrame(index=mts_wat.index, columns=['Original', 'MTS'])

levels_comparison['Original'] = profile_old_wat_extensive / 1e6
levels_comparison['MTS'] = mts_wat_extensive / 1e6

# mts_wat_extensive.to_csv(r"C:\Users\Andrea\OneDrive - Politecnico di Milano\Università\Tesi (OneDrive)\Article PROres1 coupling\mts_levels_extensive.csv")
# profile_old_wat_extensive.to_csv(r"C:\Users\Andrea\OneDrive - Politecnico di Milano\Università\Tesi (OneDrive)\Article PROres1 coupling\old_levels_extensive.csv")

# %% Demand breakdown plot

xticks_sh = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
xticks_day = ['15\nFeb', '16', '17', '18', '19', '20', '21']

levels_comparison_plot = pd.DataFrame(levels_comparison.values, columns=['Historical Level', 'MTS Level'])
p_curve_scaled_ev_f_plot = pd.DataFrame(p_curve_scaled_ev_f.values, columns=p_curve_scaled_ev_f.columns)
plotdata_plot = pd.DataFrame(plotdata.values, columns=plotdata.columns)
heat_demand_f_plot = pd.DataFrame(heat_demand_f.values, columns=heat_demand_f.columns)
ev_demand_plot = pd.DataFrame(plotdata['EV Demand'].values, columns=['EV Demand'])

fig, axes = plt.subplots(nrows=2, ncols=2, sharex=False, figsize=(10, 5))

ax = plotdata_plot.plot(kind='area', ax=axes[0, 0], rot=0, fontsize=15, legend=False)
ax = p_curve_scaled_ev_f_plot.plot(kind='line', color='black', ax=axes[0, 0], rot=0, fontsize=15, legend=False)
ax2 = heat_demand_f_plot.plot(kind='line', color='darkred', ax=axes[1, 0], rot=0, fontsize=15, legend=False)
ax3 = ev_demand_plot.plot(kind='line', color='green', ax=axes[0, 1], rot=0, fontsize=15, legend=False)
ax4 = levels_comparison_plot.plot(kind='line', ylim=(0, 1.3 * 1e2), color=['deepskyblue', 'gold'], rot=0, ax=axes[1, 1],
                                  fontsize=15, legend=False)

# ylim = (0, 1.2*1e8),
# ax.set_ylabel('Power', fontsize=15)
# ax2.set_xlabel('Time', fontsize=15)
# ax2.set_ylabel('Heat', fontsize=15)
# ax4.set_xlabel('Time', fontsize=15)
ax2.set_xticks(range(12, 157, 24))
ax2.set_xticklabels(xticks_day)
ax3.set_xticks(range(12, 157, 24))
ax3.set_xticklabels(xticks_day)
ax.set_xticks(range(12, 157, 24))
ax.set_xticklabels(xticks_day)

ax4.set_xticks(range(360, 8424, 733))
ax4.set_xticklabels(xticks_sh, rotation=90)
ax4.yaxis.set_label_position("right")
# ax4.set_ylabel('Reservoir Level', fontsize=15)
ax.set_title("Electric Demand [$GW_{el}$]", fontsize=15)
ax2.set_title("Heat Demand [$GW_{th}$]", fontsize=15)
ax3.set_title("EV Demand [$GW_{el}$]", fontsize=15)
ax4.set_title("Reservoir Levels comparison [TW]", fontsize=15)
# $n_{10}$
handles, labels = ax.get_legend_handles_labels()
myorder = [0, 3, 2, 1]
labels = [labels[i] for i in myorder]
handles = [handles[i] for i in myorder]
handles2, labels2 = ax2.get_legend_handles_labels()
handles3, labels3 = ax3.get_legend_handles_labels()
handles4, labels4 = ax4.get_legend_handles_labels()

handles.append(handles2[0])
handles.append(handles3[0])
handles.append(handles4[0])
handles.append(handles4[1])

labels.append(labels2[0])
labels.append(labels3[0])
labels.append(labels4[0])
labels.append(labels4[1])

fig.tight_layout()

ax2.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.02, -0.35), ncol=4, fontsize=15)

fig.savefig(output_folder + "\Demand Breakdown.png", bbox_inches='tight')

# %% Plot of EV demand

xticks_lab = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun']
ev_demand_plot_tick = pd.DataFrame(ev_demand_plot.values)

ax = ev_demand_plot_tick.plot(kind='line', rot=0, figsize=(10, 5), fontsize=15, legend=False)

# ax.set_xticks(np.interp(xticks, np.linspace(1, 169, 169), np.arange(ev_demand_plot.size)))
#
# ax.set_xticklabels(xticks)
# ax.set_ylabel('Adimensional EV Demand', fontsize=15)
ax.set_xlabel('Day of the Week', fontsize=15)
ax.set_title("Adimensional EV Demand", fontsize=15)

ax.set_xticks(range(12, 157, 24))
ax.set_xticklabels(xticks_lab)

fig = ax.get_figure()
fig.savefig(output_folder + "\EV Demand.png", bbox_inches='tight')

# %% Heat demand shares

heat_demand_raw = inputs_ALL['param_df']['HeatDemand']

heat_demand = pd.DataFrame(index=heat_demand_raw.index, columns=inputs_ALL['param_df']['sets']['n'])

for c in inputs_ALL['param_df']['sets']['n']:
    heat_demand.loc[:, c] = heat_demand_raw.loc[:, heat_demand_raw.columns.str.contains(' ' + c + '_')].sum(axis=1)

heat_demand_dh = pd.DataFrame(index=heat_demand_raw.index, columns=heat_demand_raw.columns)
heat_demand_p2h = pd.DataFrame(index=heat_demand_raw.index, columns=heat_demand_raw.columns)

heat_demand_dh = heat_demand_raw.loc[:, heat_demand_raw.columns.str.contains('CHP')]
heat_demand_p2h = heat_demand_raw.loc[:, heat_demand_raw.columns.str.contains('P2HT')]

heat_demand_dh_tot = heat_demand_dh.sum(axis=1)
heat_demand_p2h_tot = heat_demand_p2h.sum(axis=1)

dh_perc = (heat_demand_dh_tot.sum()) / (heat_demand.sum(axis=1).sum())
p2h_perc = (heat_demand_p2h_tot.sum()) / (heat_demand.sum(axis=1).sum())

demand = inputs_ALL['param_df']['Demand']['DA']

demand_ratio = heat_demand.sum(axis=1).sum() / demand.sum().sum()

# for t in heat_demand_raw.columns:
#    if 'P2HT' in t:
#        heat_demand_p2h.loc[:,t] = heat_demand_raw.loc[:,t].sum(axis = 1)
#    else:
#        heat_demand_dh.loc[:,t] = heat_demand_raw.loc[:,heat_demand_raw.columns.str.contains(t)].sum(axis = 1)

# %% COP plot

cop = pd.read_excel(
    r"C:\Users\Andrea\OneDrive - Politecnico di Milano\Università\Tesi (OneDrive)\Article PROres1 coupling\Simulations and results\COP values for plot.xlsx",
    index_col=0)
cop = cop.drop(index=['CY', 'IS'])

ax = cop.T.plot(kind='line', rot=0, figsize=(10, 5), fontsize=15, legend=True)

ax.grid('on', alpha=0.5)
ax.set_ylabel('COP', fontsize=15)
ax.set_xlabel('Temperature [°C]', fontsize=15)
ax.set_title("Overall combined COP for heat pumps and electrical heaters", fontsize=15)
ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.45), ncol=2, fontsize=15)

fig = ax.get_figure()
fig.savefig(output_folder + "\COP curves.png", bbox_inches='tight')

# %% NTC Congestion

# libraries
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import matplotlib as mpl

countries = ['AT', 'BE', 'BG', 'CH', 'CZ', 'DE', 'DK', 'EE', 'EL', 'ES', 'FI', 'FR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU',
             'LV', 'NL', 'NO', 'PL', 'PT', 'RO', 'SE', 'SI', 'SK', 'UK']


# Build a dataframe with your connections

def build_cong_df(r_scenario):
    cong_scenario = pd.DataFrame(columns=['from', 'to', 'value'])

    for i in range(0, len(list(r_scenario['Congestion'].keys()))):
        cong_scenario.loc[i, 'from'] = list(r_scenario['Congestion'].keys())[i][:2]
        cong_scenario.loc[i, 'to'] = list(r_scenario['Congestion'].keys())[i][-2:]
        cong_scenario.loc[i, 'value'] = list(r_scenario['Congestion'].values())[i]

    cong_scenario["value"] = pd.to_numeric(cong_scenario["value"])

    return cong_scenario


cong_ALL = build_cong_df(r_ALL)
cong_EV = build_cong_df(r_EV)
cong_TH = build_cong_df(r_TH)
cong_NO = build_cong_df(r_NO)
cong_HY = build_cong_df(r_HY)

# Delta % in congestion

cong_delta_per = pd.DataFrame(index=['Cong Delta %'], columns=scenarios)
cong_delta_per.loc[:, 'ALLFLEX'] = (cong_ALL['value'].sum() - cong_NO['value'].sum()) / (cong_NO['value'].sum())
cong_delta_per.loc[:, 'EVFLEX'] = (cong_EV['value'].sum() - cong_NO['value'].sum()) / (cong_NO['value'].sum())
cong_delta_per.loc[:, 'THFLEX'] = (cong_TH['value'].sum() - cong_NO['value'].sum()) / (cong_NO['value'].sum())
cong_delta_per.loc[:, 'HYFLEX'] = (cong_HY['value'].sum() - cong_NO['value'].sum()) / (cong_NO['value'].sum())
cong_delta_per.fillna(0, inplace=True)
cong_delta_per = cong_delta_per * 100


# %% Draw map function

def draw_map_cong(cong_scenario, scenario, path=False, namefile=None):
    plt.figure(figsize=(20, 10))
    #    plt.subplot(111)

    #    #Custom adjust of the subplots
    #    plt.subplots_adjust(left=0.05,right=0.9,top=0.90,bottom=0.1,wspace=0.15,hspace=0.15)
    #    plt.subplots_adjust(left=0.05,right=0.95,top=0.90,bottom=0.05,wspace=0.15,hspace=0.05)

    #    fig = plt.figure(figsize=(7,10))
    y1_shift = 4.5
    x1_shift = 6.5

    x1 = -12
    x2 = 28
    y1 = 36
    y2 = 65
    #    y1 = 32.
    #    x1 = -18.
    #    x2 = 32.

    #    lat_ts=(x1+x2)/2
    m = Basemap(projection='merc', resolution='l', lat_0=(y1 + y2) / 2, llcrnrlat=y1, urcrnrlat=y2, llcrnrlon=x1,
                urcrnrlon=x2)
    m.drawcountries(linewidth=0.2)
    m.fillcontinents(color='white', lake_color='white')
    m.drawcoastlines(linewidth=0.2)

    position = {'AT': (3431171.010925041 - x1_shift * 1e5, 2196945.658932812 - y1_shift * 1e5),
                'BE': (2409122.63439698 - x1_shift * 1e5, 2782458.279736769 - y1_shift * 1e5),
                'BG': (4708739.694503188 - x1_shift * 1e5, 1523426.2665741867 - y1_shift * 1e5),
                'CH': (2847144.9315600675 - x1_shift * 1e5, 2143261.1721419306 - y1_shift * 1e5),
                'CZ': (3668436.7387408563 - x1_shift * 1e5, 2596767.239598376 - y1_shift * 1e5),
                'DE': (2956650.5058508394 - x1_shift * 1e5, 2811413.149362125 - y1_shift * 1e5),
                'DK': (3066156.080141611 - x1_shift * 1e5, 3733264.972073233 - y1_shift * 1e5),
                'EE': (4818245.26879396 - x1_shift * 1e5, 4345101.374055384 - y1_shift * 1e5),
                'EL': (4380222.971630873 - x1_shift * 1e5, 942744.7517043282 - y1_shift * 1e5),
                'ES': (1533078.0400708055 - x1_shift * 1e5, 1084664.6122998544 - y1_shift * 1e5),
                'FI': (4818245.26879396 - x1_shift * 1e5, 5495424.022817817 - y1_shift * 1e5),
                'FR': (2190111.4858154366 - x1_shift * 1e5, 1984171.5761665502 - y1_shift * 1e5),
                'HR': (3768436.7387408563 - x1_shift * 1e5, 1623426.2665741867 - y1_shift * 1e5),
                'HU': (4161211.823049329 - x1_shift * 1e5, 2143261.1721419306 - y1_shift * 1e5),
                'IE': (1095055.7429077183 - x1_shift * 1e5, 3167223.489534808 - y1_shift * 1e5),
                'IT': (3176418.223779655 - x1_shift * 1e5, 1498500.042589245 - y1_shift * 1e5),
                'LT': (4599234.120212416 - x1_shift * 1e5, 3533264.972073233 - y1_shift * 1e5),
                'LU': (2646388.3622127953 - x1_shift * 1e5, 2596767.239598376 - y1_shift * 1e5),
                'LV': (4708739.694503188 - x1_shift * 1e5, 3931681.3009275147 - y1_shift * 1e5),
                'NL': (2600757.3894058308 - x1_shift * 1e5, 3076765.6906770044 - y1_shift * 1e5),
                'NO': (3066156.080141611 - x1_shift * 1e5, 5012797.17202571 - y1_shift * 1e5),
                'PL': (4161211.823049329 - x1_shift * 1e5, 2987330.8735358263 - y1_shift * 1e5),
                'PT': (1095055.7429077183 - x1_shift * 1e5, 1013449.4429983455 - y1_shift * 1e5),
                'RO': (4708739.694503188 - x1_shift * 1e5, 1984171.5761665502 - y1_shift * 1e5),
                'SE': (3613683.9515954703 - x1_shift * 1e5, 5012797.17202571 - y1_shift * 1e5),
                'SI': (3613683.9515954703 - x1_shift * 1e5, 1984171.5761665502 - y1_shift * 1e5),
                'SK': (4106459.0359039432 - x1_shift * 1e5, 2415178.327890961 - y1_shift * 1e5),
                'UK': (1752089.1886523492 - x1_shift * 1e5, 3351332.1756189675 - y1_shift * 1e5)}

    # Add color bar to the right
    colors = range(8784)
    #    cmap = plt.cm.hot_r
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["green", "orange", "red"])
    vmin = min(colors)
    vmax = max(colors)

    G = nx.from_pandas_edgelist(cong_scenario, 'from', 'to', edge_attr='value', create_using=nx.DiGraph(directed=True))
    colors = [i['value'] for i in dict(G.edges).values()]

    nx.draw(G, position, with_labels=True, arrows=True, connectionstyle='arc3, rad = 0.1', edge_color=colors,
            node_color='lightblue', node_size=600, font_size=15, edge_cmap=cmap)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.ax.tick_params(labelsize=15)
    plt.title("Map of Congestion - " + scenario + " Scenario", fontsize=15)

    if path is not False:
        fig.savefig(path + '/' + namefile, bbox_inches='tight')
        plt.show()
    else:
        plt.show()


# Example

# path = output_folder
# draw_map_cong(cong_ALL, scenario = 'ALLFLEX', path = path, namefile = 'NTC ALLFLEX')

# %% Draw the NTC map

path = output_folder

draw_map_cong(cong_ALL, scenario='ALLFLEX', path=path, namefile='NTC ALLFLEX')
draw_map_cong(cong_EV, scenario='EVFLEX', path=path, namefile='NTC EVFLEX')
draw_map_cong(cong_TH, scenario='THFLEX', path=path, namefile='NTC THFLEX')
draw_map_cong(cong_NO, scenario='NOFLEX', path=path, namefile='NTC NOFLEX')
draw_map_cong(cong_HY, scenario='HYFLEX', path=path, namefile='NTC HYFLEX')

# %% Dis
# patch plot

# %% ALLFLEX
z = 'AT'

rng = pd.date_range(start='2016-02-20', end='2016-02-27', freq='h')

ds.plot_zone(inputs_ALL, results_ALL, rng=rng, z=z)

# %% EVFLEX
z = 'AT'

rng = pd.date_range(start='2016-02-20', end='2016-02-27', freq='h')

ds.plot_zone(inputs_EV, results_EV, rng=rng, z=z)

# %% THFLEX
z = 'AT'

rng = pd.date_range(start='2016-02-20', end='2016-02-27', freq='h')
ds.plot_zone(inputs_TH, results_TH, rng=rng, z=z)

# %% HYFLEX
z = 'AT'

rng = pd.date_range(start='2016-02-20', end='2016-02-27', freq='h')
ds.plot_zone(inputs_HY, results_HY, rng=rng, z=z)

# %% NOFLEX
z = 'AT'

rng = pd.date_range(start='2016-02-20', end='2016-02-27', freq='h')
ds.plot_zone(inputs_NO, results_NO, rng=rng, z=z)

# %% Levels comparison plot

# xticks_sh = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
# levels_comparison_plot = pd.DataFrame(levels_comparison.values, columns = levels_comparison.columns)
#
# ax = levels_comparison_plot.plot(kind='line', rot = 0, figsize = (7,5), fontsize=15)
#
# ax.set_ylabel('Reservoir Level', fontsize=15)
# ax.set_xlabel('Time', fontsize=15)
# ax.set_xticks(range(1,8784,799))
# ax.set_xticklabels(xticks_sh)
#
# handles, labels = ax.get_legend_handles_labels()
#
# ax.legend(handles, labels, loc = 4, fontsize=15)
# ax.set_title("Reservoir Levels comparison", fontsize=15)
#
# fig = ax.get_figure()
# fig.savefig(r"C:\Users\Andrea\OneDrive - Politecnico di Milano\Università\Tesi (OneDrive)\Article PROres1 coupling\Simulations and results\Charts\Levels comparison.png", bbox_inches='tight')

# %% Cost load curves plots

# TODO: Not used plot, if needed needs to be fixed

shad_price_country = pd.DataFrame(index=costs_ALL[0]['CostHeat'].index, columns=('ALL', 'EV', 'HY', 'NO', 'TH'))
cost_el_plot['ALL'] = (costs_ALL[1]['Spillage'] - (costs_ALL[0]['CostHeat'] + costs_ALL[0]['CostHeatSlack'])) / 10 ** 6
cost_el_plot['EV'] = (costs_EV[1]['Spillage'] - (costs_EV[0]['CostHeat'] + costs_EV[0]['CostHeatSlack'])) / 10 ** 6
cost_el_plot['TH'] = (costs_TH[1]['Spillage'] - (costs_TH[0]['CostHeat'] + costs_TH[0]['CostHeatSlack'])) / 10 ** 6
cost_el_plot['NO'] = (costs_NO[1]['Spillage'] - (costs_NO[0]['CostHeat'] + costs_NO[0]['CostHeatSlack'])) / 10 ** 6
cost_el_plot['HY'] = (costs_HY[1]['Spillage'] - (costs_HY[0]['CostHeat'] + costs_HY[0]['CostHeatSlack'])) / 10 ** 6

el.plot_LDC(cost_el_plot, x_norm=False, legend=True)
plt.ylabel('Daily Generation Cost Electricity [M€]')
plt.xlabel('Hour')
plt.legend()

cost_heat_plot = pd.DataFrame(index=costs_ALL[0]['CostHeat'].index, columns=('ALL', 'EV', 'HY', 'NO', 'TH'))
cost_heat_plot['ALL'] = (costs_ALL[0]['CostHeat'] + costs_ALL[0]['CostHeatSlack']) / 10 ** 6
cost_heat_plot['EV'] = (costs_EV[0]['CostHeat'] + costs_EV[0]['CostHeatSlack']) / 10 ** 6
cost_heat_plot['TH'] = (costs_TH[0]['CostHeat'] + costs_TH[0]['CostHeatSlack']) / 10 ** 6
cost_heat_plot['NO'] = (costs_NO[0]['CostHeat'] + costs_NO[0]['CostHeatSlack']) / 10 ** 6
cost_heat_plot['HY'] = (costs_HY[0]['CostHeat'] + costs_HY[0]['CostHeatSlack']) / 10 ** 6

el.plot_LDC(cost_heat_plot, x_norm=False, legend=True)
plt.ylabel('Daily Generation Cost for Heat [M€]')
plt.xlabel('Hour')

# %% Additional functions from the read results dispaset file

# if needed, define the plotting range for the dispatch plot:
import pandas as pd

rng = pd.date_range(start='2016-01-01', end='2016-12-31', freq='h')

# Generate country-specific plots
ds.plot_zone(inputs, results, rng=rng, z=z)

# Bar plot with the installed capacities in all countries:
cap = ds.plot_zone_capacities(inputs)

# Bar plot with the energy balances in all countries:
ds.plot_energy_zone_fuel(inputs, results, ds.get_indicators_powerplant(inputs, results))

# Calculate costs per type
costs = ds.CostExPost(inputs, results)

# Plot the reservoir levels
ds.storage_levels(inputs, results)

def drop_unit_data(r,scenario):
    r['UnitData'].to_csv(scenario + '.csv')
    # inputs['units'].to_csv('UnitInputs.csv')

drop_unit_data(r_R_H,'R_H')
drop_unit_data(r_R_W,'R_W')
drop_unit_data(r_R_D,'R_D')
drop_unit_data(r_RN_H,'RN_H')
drop_unit_data(r_RN_W,'RN_W')
drop_unit_data(r_RN_D,'RN_D')
drop_unit_data(r_T_H,'T_H')
drop_unit_data(r_T_W,'T_W')
drop_unit_data(r_T_D,'T_D')
drop_unit_data(r_TN_H,'TN_H')
drop_unit_data(r_TN_W,'TN_W')
drop_unit_data(r_TN_D,'TN_D')

def drop_generation(inputs, results, scenario):
    aa = ds.plot_energy_zone_fuel(inputs, results, ds.get_indicators_powerplant(inputs, results))
    aa.to_csv(scenario + '.csv')

drop_generation(inputs_R_H, results_R_H, 'R_H')
drop_generation(inputs_R_W, results_R_W, 'R_W')
drop_generation(inputs_R_D, results_R_D, 'R_D')
drop_generation(inputs_RN_H, results_RN_H, 'RN_H')
drop_generation(inputs_RN_W, results_RN_W, 'RN_W')
drop_generation(inputs_RN_D, results_RN_D, 'RN_D')
drop_generation(inputs_T_H, results_T_H, 'T_H')
drop_generation(inputs_T_W, results_T_W, 'T_W')
drop_generation(inputs_T_D, results_T_D, 'T_D')
drop_generation(inputs_TN_H, results_TN_H, 'TN_H')
drop_generation(inputs_TN_W, results_TN_W, 'TN_W')
drop_generation(inputs_TN_D, results_TN_D, 'TN_D')

