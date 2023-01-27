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
config = ds.load_config('../ConfigFiles/Config_BOLIVIA_POWERFLOWDC_TEST.xlsx')

# # Limit the simulation period (for testing purposes, comment the line to run the whole year)
# config['StartDate'] = (2026, 1, 1, 0, 0, 0)
# config['StopDate'] = (2026, 7, 1, 0, 0, 0)

# Build the simulation environment:
SimData = ds.build_simulation(config)

# # Solve using GAMS:
_ = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])

# Load the simulation results:
inputs,results = ds.get_sim_results(path='../Simulations/BOLIVIA_POWERFLOWDC_TEST',cache=False)
# inputs,results = ds.get_sim_results(path='../Simulations/BOLIVIA_POWERFLOWDC_TEST',cache=False, inputs_file='Inputs_MTS.p',results_file='Results_MTS.gdx')
# inputs, results = ds.get_sim_results(config, cache=False)
inputs_MTS, results_MTS = ds.get_sim_results(path='../Simulations/BOLIVIA_POWERFLOWDC_TEST', cache=False, inputs_file='Inputs_MTS.p',results_file='Results_MTS.gdx')

# # import pandas as pd
import pandas as pd

# Generate country-specific plots
#ds.plot_zone(inputs, results, rng=rng)
rng = pd.date_range('2026-01-01', '2026-01-07', freq='H')
for i in config['zones']:
    ds.plot_zone(inputs, results, z=i, rng=rng)
    
rng = pd.date_range('2026-01-01', '2026-01-07', freq='H')
for i in config['zones']:
    ds.plot_zone(inputs, results, z=i, rng=rng)
    
rng1 = pd.date_range('2026-01-01', '2026-01-07', freq='D')
for j in config['zones']:
    ds.plot_zone(inputs_MTS, results_MTS, z=j, rng=rng1)

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


# #%% ########RESULTADOS EN CSV

# # results['CapacityMargin'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/CapacityMargin.csv', header=True, index=True)
# # #results['EESH-ExpectedEnergyNotServed'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/EESH-ExpectedEnergyNotServed.csv', header=True, index=True)
# # results['ENSH-EnergyNotServedHourly'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/ENSH-EnergyNotServedHourly.csv', header=True, index=True)
# # results['ENSR-EnergyNotServedRamping'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/ENSR-EnergyNotServedRamping.csv', header=True, index=True)
# # results['H2ShadowPrice'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/H2ShadowPrice.csv', header=True, index=True)
# # results['HeatShadowPrice'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/HeatShadowPrice.csv', header=True, index=True)
# # # results['LOLF-LosOfLoadFrequency'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/LOLF-LosOfLoadFrequency.csv', header=True, index=True)
# # # results['LOLH-LosOfLoadHours'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/LOLH-LosOfLoadHours.csv', header=True, index=True)
# # # results['LOLP-LosOfLoadProbability'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/LOLP-LosOfLoadProbabilityy.csv', header=True, index=True)
# # results['LostLoad_2D'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/LostLoad_2D.csv', header=True, index=True)
# # results['LostLoad_2U'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/LostLoad_2U.csv', header=True, index=True)
# # results['LostLoad_3U'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/LostLoad_3U.csv', header=True, index=True)
# # results['LostLoad_MaxPower'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/LostLoad_MaxPower.csv', header=True, index=True)
# # results['LostLoad_MinPower'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/LostLoad_MinPower.csv', header=True, index=True)
# # results['LostLoad_RampDown'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/LostLoad_RampDown.csv', header=True, index=True)
# # results['LostLoad_RampUp'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/LostLoad_RampUp.csv', header=True, index=True)
# # # results['LostLoad_WaterSlack'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/LostLoad_WaterSlack.csv', header=True, index=True)
# # results['NodalPowerConsumption'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/NodalPowerConsumption.csv', header=True, index=True)
# # results['OutputCommitted'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputCommitted.csv', header=True, index=True)
# # results['OutputCostRampUpH'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputCostRampUpH.csv', header=True, index=True)
# # results['OutputCostStartUpH'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputCostStartUpH.csv', header=True, index=True)
# # results['OutputCurtailedHeat'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputCurtailedHeat.csv', header=True, index=True)
# # results['OutputCurtailedPower'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputCurtailedPower.csv', header=True, index=True)
# # results['OutputCurtailmentReserve_2U'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputCurtailmentReserve_2U.csv', header=True, index=True)
# # results['OutputCurtailmentReserve_3U'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputCurtailmentReserve_3U.csv', header=True, index=True)
# # results['OutputDemand_2D'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputDemand_2D.csv', header=True, index=True)
# # results['OutputDemand_2U'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputDemand_2U.csv', header=True, index=True)
# # results['OutputDemand_3U'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputDemand_3U.csv', header=True, index=True)
# # results['OutputDemandModulation'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputDemandModulation.csv', header=True, index=True)
# # results['OutputEmissions'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputEmissions.csv', header=True, index=True)
# # results['OutputFlow'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputFlow.csv', header=True, index=True)
# # results['OutputH2Output'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputH2Output.csv', header=True, index=True)
# # results['OutputH2Slack'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputH2Slack.csv', header=True, index=True)
# # results['OutputHeat'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputHeat.csv', header=True, index=True)
# # results['OutputHeatSlack'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputHeatSlack.csv', header=True, index=True)
# # results['OutputMaxOutageDown'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputMaxOutageDown.csv', header=True, index=True)
# # results['OutputMaxOutageUp'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputMaxOutageUp.csv', header=True, index=True)
# # # results['OutputOptimalityGap'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputOptimalityGap.csv', header=True, index=True)
# # # results['OutputOptimizationCheck'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputOptimizationCheck.csv', header=True, index=True)
# # # results['OutputOptimizationError'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputOptimizationError.csv', header=True, index=True)
# # results['OutputPower'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputPower.csv', header=True, index=True)
# # results['OutputPowerConsumption'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputPowerConsumption.csv', header=True, index=True)
# # results['OutputPowerMustRun'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputPowerMustRun.csv', header=True, index=True)
# # results['OutputPtLDemand'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputPtLDemand.csv', header=True, index=True)
# # results['OutputRampRate'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputRampRate.csv', header=True, index=True)
# # results['OutputReserve_2D'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputReserve_2D.csv', header=True, index=True)
# # results['OutputReserve_2U'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputReserve_2U.csv', header=True, index=True)
# # results['OutputReserve_3U'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputReserve_3U.csv', header=True, index=True)
# # results['OutputShedLoad'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputShedLoad.csv', header=True, index=True)
# # results['OutputShutDown'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputShutDown.csv', header=True, index=True)
# # results['OutputSpillage'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputSpillage.csv', header=True, index=True)
# # results['OutputStartUp'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputStartUp.csv', header=True, index=True)
# # results['OutputStorageInput'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputStorageInput.csv', header=True, index=True)
# # results['OutputStorageLevel'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputStorageLevel.csv', header=True, index=True)
# # results['OutputStorageSlack'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputStorageSlack.csv', header=True, index=True)
# # # results['OutputSystemCost'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputSystemCost.csv', header=True, index=True)
# # # results['OutputSystemCostD'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/OutputSystemCostD.csv', header=True, index=True)
# # results['ShadowPrice'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/ShadowPrice.csv', header=True, index=True)
# # results['ShadowPrice_2D'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/ShadowPrice_2D.csv', header=True, index=True)
# # results['ShadowPrice_2U'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/ShadowPrice_2U.csv', header=True, index=True)
# # results['ShadowPrice_3U'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/ShadowPrice_3U.csv', header=True, index=True)
# # results['ShadowPrice_RampDown_TC'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/ShadowPrice_RampDown_TC.csv', header=True, index=True)
# # results['ShadowPrice_RampUp_TC'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/ShadowPrice_RampUp_TC.csv', header=True, index=True)
# # # results['SMML-SystemMinusesMaximalLoad'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/SMML-SystemMinusesMaximalLoad.csv', header=True, index=True)
# # # results['SMNL-SystemMinusesNominalLoad'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/SMNL-SystemMinusesNominalLoad.csv', header=True, index=True)
# # results['status'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/status.csv', header=True, index=True)
# # results['StorageShadowPrice'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/StorageShadowPrice.csv', header=True, index=True)
# # results['TotalDemand'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/TotalDemand.csv', header=True, index=True)
# # results['UnitHourly2DRevenue'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/UnitHourly2DRevenue.csv', header=True, index=True)
# # results['UnitHourly2URevenue'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/UnitHourly2URevenue.csv', header=True, index=True)
# # results['UnitHourly3URevenue'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/UnitHourly3URevenue.csv', header=True, index=True)
# # results['UnitHourlyPowerRevenue'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/UnitHourlyPowerRevenue.csv', header=True, index=True)
# # results['UnitHourlyProductionCost'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/UnitHourlyProductionCost.csv', header=True, index=True)
# # results['UnitHourlyProfit'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/UnitHourlyProfit.csv', header=True, index=True)
# # results['UnitHourlyRampingCost'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/UnitHourlyRampingCost.csv', header=True, index=True)
# # results['UnitHourlyRevenue'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/UnitHourlyRevenue.csv', header=True, index=True)
# # results['UnitHourlyStartUpCost'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/UnitHourlyStartUpCost.csv', header=True, index=True)
# # results['UnitHourlyVariableCost'].to_csv('../Results/FINAL_AREAS_NOCLUSTERING/UnitHourlyVariableCost.csv', header=True, index=True)

# #%% CONSTRUYENDO TABLAS

# #TABLA 1

# tabla1 = pd.DataFrame(columns = ['ZONA', 'GENERACION HIDRO [TWh]', 'EMISIONES CO2 HIDRO [MT]', 'GENERACION TERMO [TWh]', 'EMISIONES CO2 TERMO [MT]','GENERACION VRES [TWh]', 'EMISIONES CO2 VRES [MT]'], dtype=float )
# hidro = r['UnitData'][(r['UnitData'].Fuel) == 'WAT']
# hidroce = hidro[(hidro.Zone) == 'CE']
# hidrono = hidro[(hidro.Zone) == 'NO']
# hidroor = hidro[(hidro.Zone) == 'OR']
# hidrosu = hidro[(hidro.Zone) == 'SU']

# gas = r['UnitData'][(r['UnitData'].Fuel) == 'GAS']
# oil = r['UnitData'][(r['UnitData'].Fuel) == 'OIL']
# win = r['UnitData'][(r['UnitData'].Fuel) == 'WIN']
# sun = r['UnitData'][(r['UnitData'].Fuel) == 'SUN']
# bio = r['UnitData'][(r['UnitData'].Fuel) == 'BIO']
# termo = pd.concat([gas, oil, bio], axis=0)
# termoce = termo[(termo.Zone) == 'CE']
# termono = termo[(termo.Zone) == 'NO']
# termoor = termo[(termo.Zone) == 'OR']
# termosu = termo[(termo.Zone) == 'SU']

# vres = pd.concat([win, sun], axis=0)
# vresce = vres[(vres.Zone) == 'CE']
# vresno = vres[(vres.Zone) == 'NO']
# vresor = vres[(vres.Zone) == 'OR']
# vressu = vres[(vres.Zone) == 'SU']

# ghidroce  = hidroce['Generation [TWh]'].sum()
# ghidrono  = hidrono['Generation [TWh]'].sum()
# ghidroor  = hidroor['Generation [TWh]'].sum()
# ghidrosu  = hidrosu['Generation [TWh]'].sum()

# gtermoce  = termoce['Generation [TWh]'].sum()
# gtermono  = termono['Generation [TWh]'].sum()
# gtermoor  = termoor['Generation [TWh]'].sum()
# gtermosu  = termosu['Generation [TWh]'].sum()

# gvresce  = vresce['Generation [TWh]'].sum()
# gvresno  = vresno['Generation [TWh]'].sum()
# gvresor  = vresor['Generation [TWh]'].sum()
# gvressu  = vressu['Generation [TWh]'].sum()

# ehidroce  = hidroce['CO2 [t]'].sum()
# ehidrono  = hidrono['CO2 [t]'].sum()
# ehidroor  = hidroor['CO2 [t]'].sum()
# ehidrosu  = hidrosu['CO2 [t]'].sum()

# etermoce  = termoce['CO2 [t]'].sum()
# etermono  = termono['CO2 [t]'].sum()
# etermoor  = termoor['CO2 [t]'].sum()
# etermosu  = termosu['CO2 [t]'].sum()

# evresce  = vresce['CO2 [t]'].sum()
# evresno  = vresno['CO2 [t]'].sum()
# evresor  = vresor['CO2 [t]'].sum()
# evressu  = vressu['CO2 [t]'].sum()

# datat1 = data=[['CE',ghidroce,ehidroce,gtermoce,etermoce,gvresce,evresce],
#                ['NO',ghidrono,ehidrono,gtermono,etermono,gvresno,evresno],
#                ['OR',ghidroor,ehidroor,gtermoor,etermoor,gvresor,evresor],
#                ['SU',ghidrosu,ehidrosu,gtermosu,etermosu,gvressu,evressu]]

# tabla1=pd.DataFrame(datat1,columns = ['Zona','Generacion HIDRO [TWh]', 'Emisiones CO2 HIDRO [MT]', 'Generacion TERMO [TWh]', 'Emisiones CO2 TERMO [MT]','Generacion VRES [TWh]', 'Emisiones CO2 VRES [MT]'])

# ### SHEDLOAD por zona

# outshed = results['OutputShedLoad'].sum()
# outshed['OR'] = 0 
# outshed['SU'] = 0 
# outshed = outshed.to_frame(name='OutputShedLoad')


# ### CURTAILMENT por zona

# outcurt = results['OutputCurtailedPower'].sum()
# outcurt['NO'] = 0 
# outcurt = outcurt.to_frame(name='OutputCurtailedPower')
# outcurt.reset_index(level =0, inplace = True)
# outcurt['newindex'] = [0,2,3,1]
# outcurt.set_index('newindex',inplace=True, drop=True)	
# outcurt = outcurt.sort_index()
# outcurt.set_index('index',inplace=True, drop=True)	

# ### LOSTLOADS por zona

# outlost = pd.DataFrame()
# outlost['Zona'] = ['CE','NO','OR','SU'] 
# outlost['LostLoads'] = [0,0,0,0] 
# outlost.set_index('Zona',inplace=True, drop=True)	

# ### SHADOWPRICES por zona

# outshad = results['ShadowPrice'].mean()
# outshad = outshad.to_frame(name='ShadowPrice')

# ###Spillage por zona
# import numpy as np
# outspi = results['OutputSpillage']
# outspi1= outspi.transpose()
# outspi1.reset_index(level =0, inplace = True)

# buslist = r['UnitData']
# # buslist = buslist.drop(['level_0'], axis=1)
# buslist.reset_index(level =0, inplace = True)
# buslist =  buslist[['index','Zone']]
# buslist = buslist.rename({'index':'Busname'}, axis=1)	

# busname = buslist.Busname
# busname.to_dict()
# busname = np.asarray(busname)
# zone = buslist.Zone
# zone.to_dict()
# zone = np.asarray(zone)
# for i, j in zip(busname, zone):
#     outspi1['index'] = outspi1['index'].str.replace(i,j, regex=False)

# outspi2= outspi1.transpose()
# outspi2.columns = outspi2.iloc[0]
# outspi2 = outspi2[1:]
# outspi3 = outspi2.groupby(level=0, axis=1).sum()

# outspice  = outspi3['CE'].sum()
# outspino  = outspi3['NO'].sum()
# outspisu  = outspi3['SU'].sum()
# outspior  = 0

# datat2 = data2=[['CE',ghidroce,ehidroce,gtermoce,etermoce,gvresce,evresce,outspice],
#                ['NO',ghidrono,ehidrono,gtermono,etermono,gvresno,evresno,outspino],
#                ['OR',ghidroor,ehidroor,gtermoor,etermoor,gvresor,evresor,outspior],
#                ['SU',ghidrosu,ehidrosu,gtermosu,etermosu,gvressu,evressu,outspisu]]

# tabla2=pd.DataFrame(datat2,columns = ['Zona','Generacion HIDRO [TWh]', 'Emisiones CO2 HIDRO [MT]', 
#                                       'Generacion TERMO [TWh]', 'Emisiones CO2 TERMO [MT]','Generacion VRES [TWh]', 
#                                       'Emisiones CO2 VRES [MT]','OutputSpillage'])
# tabla2.set_index('Zona',inplace=True, drop=True)
# frames = [tabla2, outshed, outcurt, outlost, outshad]
# tabla2 = pd.concat(frames,axis=1)

# tabla2.to_csv('../Results/ENDE_SCENARIO1/1.resumen.csv', header=True, index=True)


# #%% PLOT GENERATION BY ZONES

# import numpy as np
# outpow = results['OutputPower']
# outpow1= outpow.transpose()
# outpow1.reset_index(inplace = True)

# buslist = r['UnitData']
# # buslist = buslist.drop(['level_0'], axis=1)
# # buslist.reset_index(inplace = True)
# buslist =  buslist[['index','Fuel','Zone']]
# buslist['FuelZone'] = buslist['Fuel'] + '-' + buslist['Zone']
# buslist = buslist.rename({'index':'Busname'}, axis=1)	

# busname = buslist.Busname
# busname.to_dict()
# busname = np.asarray(busname)
# fuelzone = buslist.FuelZone
# fuelzone.to_dict()
# fuelzone = np.asarray(fuelzone)
# for i, j in zip(busname, fuelzone):
#     outpow1['index'] = outpow1['index'].str.replace(i,j, regex=False)

# outpow2= outpow1.transpose()
# outpow2.columns = outpow2.iloc[0]
# outpow2 = outpow2[1:]
# outpow3 = outpow2.groupby(level=0, axis=1).sum()

# df = pd.DataFrame()


# df = outpow3
# df['Datetime'] = pd.date_range(start='2025-12-31 23:00:00+00:00', end='2026-12-31 22:00:00+00:00', freq='H')
# df['Datetime'] = pd.to_datetime(df['Datetime'])
# df = df.drop(['level_0','index'], axis=1)
# df = df.rename({'Datetime':'TIMESTAMP'}, axis=1)	


# df['Total Hidro'] = df['WAT-CE'] + df['WAT-NO'] + df['WAT-SU']
# df['Total Term'] = df['GAS-CE'] + df['GAS-OR'] + df['GAS-NO'] + df['GAS-SU'] + df['BIO-NO'] + df['BIO-OR']
# df['Total Renov'] = df['SUN-CE'] + df['SUN-SU'] + df['WIN-CE'] + df['WIN-OR'] + df['WIN-SU']


# mask1 = (df['TIMESTAMP'] > '2026-06-01 22:00:00+00:00') & (df['TIMESTAMP'] <= '2026-06-07 23:00:00+00:00')
# df = df.loc[mask1]

# #%% STACKED AREA CHART
# import numpy as np
# import pandas as pd
# from plotly.offline import plot
# import plotly.graph_objs as go


# #CENTRAL
# x = df['TIMESTAMP']

# gasce = dict(
#     x = x,
#     y = df['GAS-CE'],
#     hovertemplate = '<b>%{x}'+
#                    '<br>Gas CE: <b>%{y}',
#     mode = 'lines',
#     line = dict(width = 0.5,
#                 color = 'darkorange'),
#     stackgroup = 'one'
#     )

# watce = dict(
#     x = x,
#     y = df['WAT-CE'],
#     hovertemplate = '<b>%{x}'+
#                 '<br>Hidro CE: <b>%{y}',
#     mode = 'lines',
#     line = dict(width = 0.5,
#                 color = 'steelblue'),
#     stackgroup = 'one'
#     )

# wince = dict(
#     x = x,
#     y = df['WIN-CE'],
#     hovertemplate = '<b>%{x}'+
#                 '<br>Wind CE: <b>%{y}',
#     mode = 'lines',
#     line = dict(width = 0.5,
#                 color = 'forestgreen'),
#     stackgroup = 'one'
#     )

# sunce = dict(
#     x = x,
#     y = df['SUN-CE'],
#     hovertemplate = '<b>%{x}'+
#                 '<br>Sun CE: <b>%{y}',
#     mode = 'lines',
#     line = dict(width = 0.5,
#                 color = 'yellow'),
#     stackgroup = 'one'
#     )

# data = [gasce, watce, wince, sunce]

# fig = dict(data = data)
# hovermode = "x unified"

# plot(fig, filename = 'stacked-area-plot-hover', validate = False)





# # #%%
# # #filtrar valores no nulos de OutputSpillage
# # OutputSpillage = results['OutputSpillage']
# # OutputSpillage_mask = OutputSpillage.sum(axis = 1) != 0 
# # OutputSpillage_filtered = OutputSpillage[OutputSpillage_mask]  #ver esto
# # OutputSpillage_filtered.reset_index(level =0, inplace = True)
# # values = OutputSpillage_filtered['index']
# # OutputSpillage_filtered.to_csv('../Results/FINAL_AREAS_NOCLUSTERING/2.OutputSpillage_filtered.csv', header=True, index=True)

# # #filtrar fechas con spillage en OutputFlow
# # OutputFlow = results['OutputFlow']
# # OutputFlow_filtered = OutputFlow[OutputFlow['index'].isin(values)] #ver esto
# # OutputFlow_filtered.to_csv('../Results/FINAL_AREAS_NOCLUSTERING/3.OutputFlow_filtered.csv', header=True, index=True)

# # #filtrar fechas con spillage en OutputPower
# # units_spillage = OutputSpillage.columns
# # OutputPower = results['OutputPower']
# # OutputPower = OutputPower[units_spillage]
# # OutputPower.reset_index(level =0, inplace = True)
# # OutputPower_filtered = OutputPower[OutputFlow['index'].isin(values)] #ver esto
# # OutputPower_filtered['index'] = OutputFlow_filtered['index']
# # OutputPower_filtered.to_csv('../Results/FINAL_AREAS_NOCLUSTERING/4.OutputPower_filtered.csv', header=True, index=True)

# # #promedio de shadowprice
# # ShadowPrice = results['ShadowPrice']
# # ShadowPriceCE = ShadowPrice['CE'].mean()
# # ShadowPriceNO = ShadowPrice['NO'].mean()
# # ShadowPriceOR = ShadowPrice['OR'].mean()
# # ShadowPriceSU = ShadowPrice['SU'].mean()
# # #%% OTRAS TABLAS
# # # TABLA 5 DATOS DEL SISTEMA POR ZONAS

# # tabla5 = pd.concat([results['CapacityMargin'].sum(), results['EENS-ExpectedEnergyNotServed'],
# #                     results['ENSH-EnergyNotServedHourly'].sum(),
# #                     results['ENSR-EnergyNotServedRamping'].sum(),
# #                     results['LOLF-LosOfLoadFrequency'],
# #                     results['LOLH-LosOfLoadHours'],
# #                     results['LOLP-LosOfLoadProbability'],
# #                     results['LostLoad_2D'].sum(),results['LostLoad_2U'].sum(),
# #                     results['LostLoad_3U'].sum(),results['LostLoad_MaxPower'].sum(),
# #                     results['LostLoad_MinPower'].sum(),results['LostLoad_RampDown'].transpose (),
# #                     results['OutputCurtailmentReserve_2U'].sum(),
# #                     results['OutputCurtailmentReserve_3U'].sum(),results['OutputDemand_2D'].sum(),
# #                     results['OutputDemand_2U'].sum(),results['OutputDemand_3U'].sum(),
# #                     results['OutputMaxOutageDown'].sum(),results['OutputMaxOutageUp'].sum(),
# #                     results['OutputShedLoad'].sum(),results['SMML-SystemMinusesMaximalLoad'],
# #                     results['SMNL-SystemMinusesNominalLoad']], axis=1)

# # tabla5.rename(columns = {0:'CapacityMargin',1:'EENS-ExpectedEnergyNotServed',
# #                          2:'ENSH-EnergyNotServedHourly',3:'ENSR-EnergyNotServedRamping',
# #                          4:'LOLF-LosOfLoadFrequency',5:'LOLH-LosOfLoadHours',
# #                          6:'LOLP-LosOfLoadProbability',7:'LostLoad_2D',8:'LostLoad_2U',
# #                          9:'LostLoad_3U',10:'LostLoad_MaxPower',11:'LostLoad_MinPower',
# #                          '2026-01-01 00:00:00':'LostLoad_RampDown',
# #                          12:'OutputCurtailmentReserve_2U',13:'OutputCurtailmentReserve_3U',
# #                          14:'OutputDemand_2D',15:'OutputDemand_2U',16:'OutputDemand_3U',
# #                          17:'OutputMaxOutageDown',18:'OutputMaxOutageUp',19:'OutputShedLoad',
# #                          20:'SMML-SystemMinusesMaximalLoad',21:'SMNL-SystemMinusesNominalLoad'}, 
# #               inplace=True)

# # ########PONER ESTO EN RESERVA ROTANTE Y FRIA
# # ShadowPrice = pd.DataFrame(columns=['ShadowPrice'])
# # ShadowPrice['ShadowPrice'] = results['ShadowPrice']

# # ShadowPrice_2D = pd.DataFrame(columns=['ShadowPrice_2D'])
# # ShadowPrice_2D['ShadowPrice_2D'] = results['ShadowPrice_2D']

# # ShadowPrice_2U = pd.DataFrame(columns=['ShadowPrice_2U'])
# # ShadowPrice_2U['ShadowPrice_2U'] = results['ShadowPrice_2U']

# # ShadowPrice_3U = pd.DataFrame(columns=['ShadowPrice_3U'])
# # ShadowPrice_3U['ShadowPrice_3U'] = results['ShadowPrice_3U']

# # StorageShadowPrice = pd.DataFrame(columns=['StorageShadowPrice'])
# # StorageShadowPrice['StorageShadowPrice'] = results['StorageShadowPrice']
# # #########


# # ########ANALIZAR DONDE PONER ESTO????

# # OutputStorageLevel = pd.DataFrame(columns=['OutputStorageLevel'])
# # OutputStorageLevel['OutputStorageLevel'] = results['OutputStorageLevel'].sum()

# # OutputPower = pd.DataFrame(columns=['OutputPower'])
# # OutputPower['OutputPower'] = results['OutputPower'].sum()

# # OutputStorageSlack = pd.DataFrame(columns=['OutputStorageSlack'])
# # OutputStorageSlack['OutputStorageSlack'] = results['OutputStorageSlack'].sum()

# # #########



# # ########PONER ESTO EN EMISIONES
# # OutputEmissions = pd.DataFrame(columns=['OutputEmissions'])
# # OutputEmissions['OutputEmissions'] = results['OutputEmissions'].sum()
# # #########




# # ########PONER ESTO EN CONGESTION
# # OutputFlow = pd.DataFrame(columns=['OutputFlow'])
# # OutputFlow['OutputFlow'] = results['OutputFlow'].sum()
# # ############


# # ########PONER ESTO EN POR UNIDADES DE GENERACION
# # OutputPower = pd.DataFrame(columns=['OutputPower'])
# # OutputPower['OutputPower'] = results['OutputPower'].sum()

# # OutputPowerMustRun = pd.DataFrame(columns=['OutputPowerMustRun'])
# # OutputPowerMustRun['OutputPowerMustRun'] = results['OutputPowerMustRun'].sum()

# # OutputRampRate = pd.DataFrame(columns=['OutputRampRate'])
# # OutputRampRate['OutputRampRate'] = results['OutputRampRate'].sum()

# # OutputReserve_2D = pd.DataFrame(columns=['OutputReserve_2D'])
# # OutputReserve_2D['OutputReserve_2D'] = results['OutputReserve_2D'].sum()

# # OutputReserve_2U = pd.DataFrame(columns=['OutputReserve_2U'])
# # OutputReserve_2U['OutputReserve_2U'] = results['OutputReserve_2U'].sum()

# # OutputReserve_3U = pd.DataFrame(columns=['OutputReserve_3U'])
# # OutputReserve_3U['OutputReserve_3U'] = results['OutputReserve_3U'].sum()

# # OutputShutDown = pd.DataFrame(columns=['OutputShutDown'])
# # OutputShutDown['OutputShutDown'] = results['OutputShutDown'].sum()

# # OutputStartUp = pd.DataFrame(columns=['OutputStartUp'])
# # OutputStartUp['OutputStartUp'] = results['OutputStartUp'].sum()

# # OutputStorageLevel = pd.DataFrame(columns=['OutputStorageLevel'])
# # OutputStorageLevel['OutputStorageLevel'] = results['OutputStorageLevel'].sum()

# # OutputStorageSlack = pd.DataFrame(columns=['OutputStorageSlack'])
# # OutputStorageSlack['OutputStorageSlack'] = results['OutputStorageSlack'].sum()



# # ############COSTO DEL SISTEMA

# # OutputSystemCost = pd.DataFrame(columns=['OutputSystemCost'])
# # OutputSystemCost['OutputSystemCost'] = results['OutputSystemCost'].sum()





