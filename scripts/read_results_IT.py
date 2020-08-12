# -*- coding: utf-8 -*-
"""
Minimalist example file showing how to read the results of a Dispa-SET run

@author: Sylvain Quoilin
"""

import pandas as pd
import numpy as np
import sys,os
sys.path.append(os.path.abspath('..'))
import dispaset as ds
from matplotlib import pyplot as plt

#%% options

LoadData = True
CalculateVariables = True
cost = False
plot = False
# %%Load the inputs and the results of the simulations 
#all scenarios
#scenarios = ['s1','s2_A','s2_B','s3_A','s3_B','s4_A','s4_B','s5','s6_A','s6_B','s7_A','s7_B','s8_A','s8_B','s9','s10_A','s10_B','s11_A','s11_B','s12_A','s12_B','s13','s14','s15','s16','s20','s30','s100','s110']
#what_if_scenarios 
#scenarios = ['s1','s2_A','s2_B','s3_A','s3_B','s4_A','s4_B','s5','s6_A','s6_B','s7_A','s7_B','s8_A','s8_B','s9','s10_A','s10_B','s11_A','s11_B','s12_A','s12_B']
#prores1 scenarios
scenarios = ['s13']#,'s14','s15','s16']
#technical potential scenarios
#scenarios = ['s20','s30','s40','s60','s70','s80','s100','s110','s120']
#no storage constraint_scenarios 
#scenarios = ['s2_C','s3_C','s4_C']


if LoadData == True:
    scenarios_dic = {}
    for SCEN in scenarios:
        print(SCEN)
        PATH = '/Users/tomfranssens/Dispa-SET.git/Simulations/simulation_##'
        PATH = PATH.replace('##', SCEN)
        INPUTS = 'inputs_##'
        INPUTS = INPUTS.replace('##', SCEN)
        RESULTS = 'results_##'
        RESULTS = RESULTS.replace('##', SCEN)
        scenarios_dic[INPUTS],scenarios_dic[RESULTS] = ds.get_sim_results(path=PATH,cache=False)
        print(SCEN)

#%% index, unit and zone
if 's1' in scenarios:
    index = scenarios_dic['inputs_s1']['config']['idx']
if 's13' in scenarios:
    index = scenarios_dic['inputs_s13']['config']['idx']
    COMC_GAS='[26] - IT_COMC_GAS_CHP'
    STUR_BIO='[25] - IT_STUR_BIO_CHP'
    ICEN_GAS='[27] - IT_ICEN_GAS_CHP'
    
if 's13' in scenarios:
    index = scenarios_dic['inputs_s13']['config']['idx']
CHP_GAS = '[10] - IT_COMC_GAS_CHP'
CHP_HRD = '[13] - IT_STUR_HRD_CHP'
CHP_OIL = '[14] - IT_STUR_OIL_CHP'
P2H_unit = '[21] - IT_P2HT_OTH'


unit_list = [P2H_unit,CHP_GAS,CHP_HRD,CHP_OIL]
if 's13' in scenarios:
    unit_list = [COMC_GAS,STUR_BIO,ICEN_GAS]

zone = 'IT'

#%% tests

#%% calculations
if CalculateVariables == True:

#costs, emissions
    calculations = pd.DataFrame(0,index = scenarios, columns = ['postsystemcostIT','gams_costs','cost_savings'
                                                                'co2_power','co2_heatslack','co2_total'])
# reserve demand
    reservedemand = pd.DataFrame(0,index= scenarios,columns=['upwards','downwards'])
# total availability
    yearly_availability = {}
    yearly_availability['Up'] = pd.DataFrame(0,index= unit_list,columns=scenarios)
    yearly_availability['Down'] = pd.DataFrame(0,index= unit_list,columns=scenarios)
# hourly availability
    hourly_availability = {}
# mean availability
    availability_mean = {}
    availability_mean['2U'] = pd.DataFrame(0,index= [P2H_unit,CHP_GAS,CHP_HRD,CHP_OIL,'total'],columns=scenarios)
    availability_mean['3U'] = pd.DataFrame(0,index= [P2H_unit,CHP_GAS,CHP_HRD,CHP_OIL,'total'],columns=scenarios)
    availability_mean['Up'] = pd.DataFrame(0,index= [P2H_unit,CHP_GAS,CHP_HRD,CHP_OIL,'total'],columns=scenarios)
    availability_mean['Down'] = pd.DataFrame(0,index= [P2H_unit,CHP_GAS,CHP_HRD,CHP_OIL,'total'],columns=scenarios)
    
#availability all units    
    availability_mean_all={}
#power output all units
    poweroutput = {}
    
# total availablility
    availability_total = {}
    availability_total['2U'] = pd.DataFrame(0,index= unit_list,columns=scenarios)
    availability_total['3U'] = pd.DataFrame(0,index= unit_list,columns=scenarios)
    availability_total['Up'] = pd.DataFrame(0,index= unit_list,columns=scenarios)
    availability_total['Down'] = pd.DataFrame(0,index= unit_list,columns=scenarios)
# load shedding sum, max, amount
    load_shedding = pd.DataFrame(0,index= scenarios,columns=['sum','max','amount'])
# curtailment sum, max, amount
    curtailment = pd.DataFrame(0,index= scenarios,columns=['sum','max','amount'])
# heat provision
    heat = pd.DataFrame(0,index= [P2H_unit,CHP_GAS,CHP_HRD,CHP_OIL,'total'],columns=scenarios)
# shadow prices
    shadowprices_mean = pd.DataFrame(0,index= ['DA','2U','3U','2D'],columns=scenarios)
    shadowprices_percentages = pd.DataFrame(0,index= ['DA','2U','3U','2D'],columns=scenarios)

    shadowprices_hourly = {}
    shadowprices_hourly['DA'] = pd.DataFrame(0,index= index,columns=scenarios)
    shadowprices_hourly['2U'] = pd.DataFrame(0,index= index,columns=scenarios)
    shadowprices_hourly['3U'] = pd.DataFrame(0,index= index,columns=scenarios)
    shadowprices_hourly['2D'] = pd.DataFrame(0,index= index,columns=scenarios)
    
    hours = pd.DataFrame(0,index= scenarios,columns=['DA','2U','3U','2D'])
    
    marketprice = pd.DataFrame(0.0,index= scenarios,columns=['DA','2U','3U','2D'])

# cash flows
    cash_flow = {}
    for u in unit_list:
        cash_flow[u] = pd.DataFrame(0,index= scenarios,columns=['costs','heat','DA','2U','2D','3U','total'])
# lost load reserves
    lostloads = pd.DataFrame(0,index= scenarios,columns=['2U','3U','2D','rampup','rampdown','minpower','maxpower'])
    
# calculations
    
    for SCEN in scenarios:
        INPUTS = 'inputs_##'
        INPUTS = INPUTS.replace('##', SCEN)
        RESULTS = 'results_##'
        RESULTS = RESULTS.replace('##', SCEN)
        inputs = scenarios_dic[INPUTS]
        results = scenarios_dic[RESULTS]
        print(SCEN)
        
# Italian units
        IT_units = []
        for u in inputs['units'].index:
            if inputs['units'].at[u,'Zone'] == zone:
                IT_units.append(u)

# lost load and hours of lost load 
        if zone in results['LostLoad_2D'] and results['LostLoad_2D'][zone].sum() != 0:
            print('2D lostload in Italy')
            
        if zone in results['LostLoad_2U'] and results['LostLoad_2U'][zone].sum() != 0:
            print('2U lostload in Italy')
        if zone in results['LostLoad_3U'] and results['LostLoad_3U'][zone].sum() != 0:
            print('3U lostload in Italy')
        if zone in results['LostLoad_MaxPower'] and results['LostLoad_MaxPower'][zone].sum() != 0:
            print('LostLoad_MaxPower in Italy')
        if zone in results['LostLoad_MinPower'] and results['LostLoad_MinPower'][zone].sum() != 0:
            print('LostLoad_MinPower in Italy')
        if zone in results['LostLoad_RampDown'] and results['LostLoad_RampDown'][zone].sum() != 0:
            print('LostLoad_RampDown in Italy')
        if zone in results['LostLoad_RampUp'] and results['LostLoad_RampUp'][zone].sum() != 0:
            print('LostLoad_RampUp in Italy')
        
        if zone in results['LostLoad_2D']:
            for i in results['LostLoad_2D'].index:
                if results['LostLoad_2D'].at[i,zone] != 0:
                    lostloads.loc[SCEN,'2D']=lostloads.loc[SCEN,'2D']+1
        if zone in results['LostLoad_2U']:
            for i in results['LostLoad_2U'].index:
                if results['LostLoad_2U'].at[i,zone] != 0:
                    lostloads.loc[SCEN,'2U']=lostloads.loc[SCEN,'2U']+1
        if zone in results['LostLoad_3U']:
            for i in results['LostLoad_3U'].index:
                if results['LostLoad_3U'].at[i,zone] != 0:
                    lostloads.loc[SCEN,'3U']=lostloads.loc[SCEN,'3U']+1
        if zone in results['LostLoad_MaxPower']:
            for i in results['LostLoad_MaxPower'].index:
                if results['LostLoad_MaxPower'].at[i,zone] != 0:
                    lostloads.loc[SCEN,'maxpower']=lostloads.loc[SCEN,'maxpower']+1
        if zone in results['LostLoad_MinPower']:
            for i in results['LostLoad_MinPower'].index:
                if results['LostLoad_MinPower'].at[i,zone] != 0:
                    lostloads.loc[SCEN,'minpower']=lostloads.loc[SCEN,'minpower']+1
        if zone in results['LostLoad_RampUp']:
            for i in results['LostLoad_RampUp'].index:
                if results['LostLoad_RampUp'].at[i,zone] != 0:
                    lostloads.loc[SCEN,'rampup']=lostloads.loc[SCEN,'rampup']+1
        if zone in results['LostLoad_RampDown']:
            for i in results['LostLoad_RampDown'].index:
                if results['LostLoad_RampDown'].at[i,zone] != 0:
                    lostloads.loc[SCEN,'rampdown']=lostloads.loc[SCEN,'rampdown']+1
                    
# power output
        IT_units=[]
        for u in inputs['units'].index:
            if inputs['units'].at[u,'Zone']==zone:
                IT_units.append(u)
        poweroutput[SCEN]=pd.DataFrame(0,index= IT_units,columns=['mean'])
        for u in results['OutputPower'].columns:
            if u in IT_units:
                poweroutput[SCEN].loc[u,'mean']=results['OutputPower'][u].mean()
        
        
        
        
# reserve demand and availability (unit specific)
        TMP1,TMP2,TMP3 = ds.reserve_availability_demand(inputs,results)

        reservedemand.loc[SCEN, 'upwards']=TMP3.at[zone,'upwards']
        reservedemand.loc[SCEN, 'downwards']=TMP3.at[zone,'downwards']
        
        for unit in unit_list:
            if unit in TMP2['2U'].index:
                availability_mean['2U'].loc[unit,SCEN] = TMP2['2U'].at[unit,'mean']
                availability_total['2U'].loc[unit,SCEN] = TMP2['2U'].at[unit,'total']
            if unit in TMP2['3U'].index:
                availability_mean['3U'].loc[unit,SCEN] = TMP2['3U'].at[unit,'mean']
                availability_total['3U'].loc[unit,SCEN] = TMP2['3U'].at[unit,'total']
            if unit in TMP2['Up'].index:
                availability_mean['Up'].loc[unit,SCEN] = TMP2['Up'].at[unit,'mean']
                availability_total['Up'].loc[unit,SCEN] = TMP2['Up'].at[unit,'total']
            if unit in TMP2['Down'].index:
                availability_mean['Down'].loc[unit,SCEN] = TMP2['Down'].at[unit,'mean']
                availability_total['Down'].loc[unit,SCEN] = TMP2['Down'].at[unit,'total']
            
        availability_mean['2U'].loc['total',SCEN] = availability_mean['2U'][SCEN].sum()
        availability_mean['3U'].loc['total',SCEN] = availability_mean['3U'][SCEN].sum()
        availability_mean['Down'].loc['total',SCEN] = availability_mean['Down'][SCEN].sum()

        hourly_availability[SCEN] = TMP1
        
        availability_mean_all[SCEN]=TMP2
        for u in availability_mean_all[SCEN]['2U'].index:
            if u not in IT_units:
                availability_mean_all[SCEN]['2U'].drop(u,inplace=True)

        for u in availability_mean_all[SCEN]['3U'].index:
            if u not in IT_units:
                availability_mean_all[SCEN]['3U'].drop(u,inplace=True)

        for u in availability_mean_all[SCEN]['Down'].index:
            if u not in IT_units:
                availability_mean_all[SCEN]['Down'].drop(u,inplace=True)

        for u in availability_mean_all[SCEN]['Up'].index:
            if u not in IT_units:
                availability_mean_all[SCEN]['Up'].drop(u,inplace=True)

        
# yearly availability (up is only the secondary reserves)
        for unit in unit_list:
            if unit in results['OutputReserve_2U']:
                yearly_availability['Up'].loc[unit,SCEN] = results['OutputReserve_2U'][unit].sum()
            if unit in results['OutputReserve_2D']:
                yearly_availability['Down'].loc[unit,SCEN] = results['OutputReserve_2D'][unit].sum()

# IT system cost with heatslack
        if cost == True:
            TMP1,TMP2,TMP3 = ds.CostExPost(inputs,results,False)
            calculations.loc[SCEN,'postsystemcostIT'] = TMP3[zone].sum()
            if SCEN in ['s1','s5','s9']:
                calculations.loc[SCEN,'postsystemcostIT'] = calculations.loc[SCEN,'postsystemcostIT']+(97790400.75*50)
                
# gams system cost without heatslack
        calculations.loc[SCEN,'gams_costs'] = results['OutputSystemCost'].sum()

# co2 power system in italy
        TMP= ds.emissions(inputs,results)
        calculations.loc[SCEN,'co2_power'] = TMP.at[zone,'emissions']
                
# co2 of the heatslack in Italy (not unit specific)
        if SCEN in ['s1', 's5', 's9']:
            emissionrate_gasboiler = 0.23320 # tCO2 / MWh
            calculations.loc[SCEN,'co2_heatslack'] = 97790400.75 * emissionrate_gasboiler # MWh * tCO2 / MWh = tCO2
        for unit in unit_list:
            if unit in results['OutputHeatSlack'].columns:
                emissionrate_gasboiler = 0.23320 # tCO2 / MWh
                calculations.loc[SCEN,'co2_heatslack'] = calculations.loc[SCEN,'co2_heatslack']+(results['OutputHeatSlack'][unit].sum() * emissionrate_gasboiler) # MWh * tCO2 / MWh = tCO2
        if 's13' in scenarios:
            for unit in ['[25] - IT_STUR_BIO_CHP','[26] - IT_COMC_GAS_CHP','[27] - IT_ICEN_GAS_CHP','[24] - IT_P2HT_OTH']:
                if unit in results['OutputHeatSlack'].columns:
                    emissionrate_gasboiler = 0.23320 # tCO2 / MWh
                    calculations.loc[SCEN,'co2_heatslack'] = calculations.loc[SCEN,'co2_heatslack']+(results['OutputHeatSlack'][unit].sum() * emissionrate_gasboiler) # MWh * tCO2 / MWh = tCO2
            
# co2 total emissions for Italy
        calculations.loc[SCEN,'co2_total'] = calculations.loc[SCEN,'co2_power']+calculations.loc[SCEN,'co2_heatslack']
        
# amount of curtailment
        TMP1 = ds.curtailment(inputs,results)
        curtailment.loc[SCEN,'sum']=TMP1.at['sum',zone]
        curtailment.loc[SCEN,'max']=TMP1.at['max',zone]
        curtailment.loc[SCEN,'amount']=TMP1.at['amount',zone]
        
# amount of load shedding
        TMP1 = ds.load_shedding(inputs,results)
        load_shedding.loc[SCEN,'sum']=TMP1.at['sum',zone]
        load_shedding.loc[SCEN,'max']=TMP1.at['max',zone]
        load_shedding.loc[SCEN,'amount']=TMP1.at['amount',zone]

# heat provisioned (unit specific + total)
        for unit in unit_list:
            if unit in results['OutputHeat'].columns:
                heat.loc[unit,SCEN] = results['OutputHeat'][unit].sum()/inputs['param_df']['HeatDemand'][unit].sum()
                heat.loc['total',SCEN] = heat.loc['total',SCEN] + results['OutputHeat'][unit].sum()
        heat.loc['total',SCEN] = heat.loc['total',SCEN]/97790400.75 # total heat/ total heat demand

# shadowprices
        TMP = ds.shadowprices(results, zone)
        shadowprices_hourly['DA'][SCEN]=TMP['DA']
        shadowprices_hourly['2U'][SCEN]=TMP['2U']
        shadowprices_hourly['3U'][SCEN]=TMP['3U']
        shadowprices_hourly['2D'][SCEN]=TMP['2D']
        TMP1=[]
        if TMP['DA'].max() > 10000:
            print('remove lostload index')
            for i in TMP.index:
                if TMP.loc[i,'DA']>10000:
                    TMP.loc[i,'DA']=TMP.at[i-1,'DA']
            shadowprices_mean.loc['DA',SCEN]=TMP['DA'].mean()
            for i in TMP.index:
                if TMP.loc[i,'DA']>0:
                    TMP1.append(TMP.at[i,'DA'])
                    marketprice.at[SCEN,'DA']=sum(TMP1)/len(TMP1)
        else:
            shadowprices_mean.loc['DA',SCEN]=TMP['DA'].mean()
            for i in TMP.index:
                if TMP.loc[i,'DA']>0:
                    TMP1.append(TMP.at[i,'DA'])
                    marketprice.at[SCEN,'DA']=sum(TMP1)/len(TMP1)

        TMP1=[]
        if TMP['2U'].max() > 10000:
            print('remove lostload index')
            for i in TMP.index:
                if TMP.loc[i,'2U']>10000:
                    TMP.loc[i,'2U']=TMP.at[i-1,'2U']
            shadowprices_mean.loc['2U',SCEN]=TMP['2U'].mean()
            for i in TMP.index:
                if TMP.loc[i,'2U']>0:
                    TMP1.append(TMP.at[i,'2U'])
                    marketprice.at[SCEN,'2U']=sum(TMP1)/len(TMP1)
        else:
            shadowprices_mean.loc['2U',SCEN]=TMP['2U'].mean()
            for i in TMP.index:
                if TMP.loc[i,'2U']>0:
                    TMP1.append(TMP.at[i,'2U'])
                    marketprice.at[SCEN,'2U']=sum(TMP1)/len(TMP1)

        TMP1=[]
        if TMP['3U'].max() > 10000:
            print('remove lostload index')
            for i in TMP.index:
                if TMP.loc[i,'3U']>10000:
                    TMP.loc[i,'3U']=TMP.at[i-1,'3U']
            shadowprices_mean.loc['3U',SCEN]=TMP['3U'].mean()
            for i in TMP.index:
                if TMP.loc[i,'3U']>0:
                    TMP1.append(TMP.at[i,'3U'])
                    marketprice.at[SCEN,'3U']=sum(TMP1)/len(TMP1)
        else:
            shadowprices_mean.loc['3U',SCEN]=TMP['3U'].mean()
            for i in TMP.index:
                if TMP.loc[i,'3U']>0:
                    TMP1.append(TMP.at[i,'3U'])
                    marketprice.at[SCEN,'3U']=sum(TMP1)/len(TMP1)

        TMP1=[]
        if TMP['2D'].max() > 10000:
            print('remove lostload index')
            for i in TMP.index:
                if TMP.loc[i,'2D']>10000:
                    TMP.loc[i,'2D']=TMP.at[i-1,'2D']
            shadowprices_mean.loc['2D',SCEN]=TMP['2D'].mean()
            for i in TMP.index:
                if TMP.loc[i,'2D']>0:
                    TMP1.append(TMP.at[i,'2D'])
                    marketprice.at[SCEN,'2D']=sum(TMP1)/len(TMP1)
        else:
            shadowprices_mean.loc['2D',SCEN]=TMP['2D'].mean()
            for i in TMP.index:
                if TMP.loc[i,'2D']>0:
                    TMP1.append(TMP.at[i,'2D'])
                    marketprice.at[SCEN,'2D']=sum(TMP1)/len(TMP1)
        '''
        if SCEN in ['s1','s2_A','s3_A','s4_A','s5','s6_A','s7_A','s8_A','s9','s10_A','s11_A','s12_A']:
            shadowprices_percentages.loc['DA',SCEN] = shadowprices_mean.loc['DA',SCEN]
            shadowprices_percentages.loc['2U',SCEN] = shadowprices_mean.loc['2U',SCEN]
            shadowprices_percentages.loc['3U',SCEN] = shadowprices_mean.loc['3U',SCEN]
            shadowprices_percentages.loc['2D',SCEN] = shadowprices_mean.loc['2D',SCEN]
        elif SCEN == 's2_B':
            shadowprices_percentages.loc['DA',SCEN] = (shadowprices_mean.loc['DA','s2_B']-shadowprices_mean.loc['DA','s2_A'])/shadowprices_mean.loc['DA','s2_A']*100
            shadowprices_percentages.loc['2U',SCEN] = (shadowprices_mean.loc['2U','s2_B']-shadowprices_mean.loc['2U','s2_A'])/shadowprices_mean.loc['2U','s2_A']*100
            shadowprices_percentages.loc['3U',SCEN] = (shadowprices_mean.loc['3U','s2_B']-shadowprices_mean.loc['3U','s2_A'])/shadowprices_mean.loc['3U','s2_A']*100
            shadowprices_percentages.loc['2D',SCEN] = (shadowprices_mean.loc['2D','s2_B']-shadowprices_mean.loc['2D','s2_A'])/shadowprices_mean.loc['2D','s2_A']*100
        elif SCEN == 's3_B':
            shadowprices_percentages.loc['DA',SCEN] = (shadowprices_mean.loc['DA','s3_B']-shadowprices_mean.loc['DA','s3_A'])/shadowprices_mean.loc['DA','s3_A']*100
            shadowprices_percentages.loc['2U',SCEN] = (shadowprices_mean.loc['2U','s3_B']-shadowprices_mean.loc['2U','s3_A'])/shadowprices_mean.loc['2U','s3_A']*100
            shadowprices_percentages.loc['3U',SCEN] = (shadowprices_mean.loc['3U','s3_B']-shadowprices_mean.loc['3U','s3_A'])/shadowprices_mean.loc['3U','s3_A']*100
            shadowprices_percentages.loc['2D',SCEN] = (shadowprices_mean.loc['2D','s3_B']-shadowprices_mean.loc['2D','s3_A'])/shadowprices_mean.loc['2D','s3_A']*100
        elif SCEN == 's4_B':
            shadowprices_percentages.loc['DA',SCEN] = (shadowprices_mean.loc['DA','s4_B']-shadowprices_mean.loc['DA','s4_A'])/shadowprices_mean.loc['DA','s4_A']*100
            shadowprices_percentages.loc['2U',SCEN] = (shadowprices_mean.loc['2U','s4_B']-shadowprices_mean.loc['2U','s4_A'])/shadowprices_mean.loc['2U','s4_A']*100
            shadowprices_percentages.loc['3U',SCEN] = (shadowprices_mean.loc['3U','s4_B']-shadowprices_mean.loc['3U','s4_A'])/shadowprices_mean.loc['3U','s4_A']*100
            shadowprices_percentages.loc['2D',SCEN] = (shadowprices_mean.loc['2D','s4_B']-shadowprices_mean.loc['2D','s4_A'])/shadowprices_mean.loc['2D','s4_A']*100
        elif SCEN == 's6_B':
            shadowprices_percentages.loc['DA',SCEN] = (shadowprices_mean.loc['DA','s6_B']-shadowprices_mean.loc['DA','s6_A'])/shadowprices_mean.loc['DA','s6_A']*100
            shadowprices_percentages.loc['2U',SCEN] = (shadowprices_mean.loc['2U','s6_B']-shadowprices_mean.loc['2U','s6_A'])/shadowprices_mean.loc['2U','s6_A']*100
            shadowprices_percentages.loc['3U',SCEN] = (shadowprices_mean.loc['3U','s6_B']-shadowprices_mean.loc['3U','s6_A'])/shadowprices_mean.loc['3U','s6_A']*100
            shadowprices_percentages.loc['2D',SCEN] = (shadowprices_mean.loc['2D','s6_B']-shadowprices_mean.loc['2D','s6_A'])/shadowprices_mean.loc['2D','s6_A']*100
        elif SCEN == 's7_B':
            shadowprices_percentages.loc['DA',SCEN] = (shadowprices_mean.loc['DA','s7_B']-shadowprices_mean.loc['DA','s7_A'])/shadowprices_mean.loc['DA','s7_A']*100
            shadowprices_percentages.loc['2U',SCEN] = (shadowprices_mean.loc['2U','s7_B']-shadowprices_mean.loc['2U','s7_A'])/shadowprices_mean.loc['2U','s7_A']*100
            shadowprices_percentages.loc['3U',SCEN] = (shadowprices_mean.loc['3U','s7_B']-shadowprices_mean.loc['3U','s7_A'])/shadowprices_mean.loc['3U','s7_A']*100
            shadowprices_percentages.loc['2D',SCEN] = (shadowprices_mean.loc['2D','s7_B']-shadowprices_mean.loc['2D','s7_A'])/shadowprices_mean.loc['2D','s7_A']*100
        elif SCEN == 's8_B':
            shadowprices_percentages.loc['DA',SCEN] = (shadowprices_mean.loc['DA','s8_B']-shadowprices_mean.loc['DA','s8_A'])/shadowprices_mean.loc['DA','s8_A']*100
            shadowprices_percentages.loc['2U',SCEN] = (shadowprices_mean.loc['2U','s8_B']-shadowprices_mean.loc['2U','s8_A'])/shadowprices_mean.loc['2U','s8_A']*100
            shadowprices_percentages.loc['3U',SCEN] = (shadowprices_mean.loc['3U','s8_B']-shadowprices_mean.loc['3U','s8_A'])/shadowprices_mean.loc['3U','s8_A']*100
            shadowprices_percentages.loc['2D',SCEN] = (shadowprices_mean.loc['2D','s8_B']-shadowprices_mean.loc['2D','s8_A'])/shadowprices_mean.loc['2D','s8_A']*100
        elif SCEN == 's10_B':
            shadowprices_percentages.loc['DA',SCEN] = (shadowprices_mean.loc['DA','s10_B']-shadowprices_mean.loc['DA','s10_A'])/shadowprices_mean.loc['DA','s10_A']*100
            shadowprices_percentages.loc['2U',SCEN] = (shadowprices_mean.loc['2U','s10_B']-shadowprices_mean.loc['2U','s10_A'])/shadowprices_mean.loc['2U','s10_A']*100
            shadowprices_percentages.loc['3U',SCEN] = (shadowprices_mean.loc['3U','s10_B']-shadowprices_mean.loc['3U','s10_A'])/shadowprices_mean.loc['3U','s10_A']*100
            shadowprices_percentages.loc['2D',SCEN] = (shadowprices_mean.loc['2D','s10_B']-shadowprices_mean.loc['2D','s10_A'])/shadowprices_mean.loc['2D','s10_A']*100
        elif SCEN == 's11_B':
            shadowprices_percentages.loc['DA',SCEN] = (shadowprices_mean.loc['DA','s11_B']-shadowprices_mean.loc['DA','s11_A'])/shadowprices_mean.loc['DA','s11_A']*100
            shadowprices_percentages.loc['2U',SCEN] = (shadowprices_mean.loc['2U','s11_B']-shadowprices_mean.loc['2U','s11_A'])/shadowprices_mean.loc['2U','s11_A']*100
            shadowprices_percentages.loc['3U',SCEN] = (shadowprices_mean.loc['3U','s11_B']-shadowprices_mean.loc['3U','s11_A'])/shadowprices_mean.loc['3U','s11_A']*100
            shadowprices_percentages.loc['2D',SCEN] = (shadowprices_mean.loc['2D','s11_B']-shadowprices_mean.loc['2D','s11_A'])/shadowprices_mean.loc['2D','s11_A']*100
        elif SCEN == 's12_B':
            shadowprices_percentages.loc['DA',SCEN] = (shadowprices_mean.loc['DA','s12_B']-shadowprices_mean.loc['DA','s12_A'])/shadowprices_mean.loc['DA','s12_A']*100
            shadowprices_percentages.loc['2U',SCEN] = (shadowprices_mean.loc['2U','s12_B']-shadowprices_mean.loc['2U','s12_A'])/shadowprices_mean.loc['2U','s12_A']*100
            shadowprices_percentages.loc['3U',SCEN] = (shadowprices_mean.loc['3U','s12_B']-shadowprices_mean.loc['3U','s12_A'])/shadowprices_mean.loc['3U','s12_A']*100
            shadowprices_percentages.loc['2D',SCEN] = (shadowprices_mean.loc['2D','s12_B']-shadowprices_mean.loc['2D','s12_A'])/shadowprices_mean.loc['2D','s12_A']*100
        '''

        #hours of prices > 0 + marketprices of those hours
        for i in shadowprices_hourly['2D'].index:
            if shadowprices_hourly['2D'].at[i,SCEN] > 0:
                hours.loc[SCEN,'2D']=hours.loc[SCEN,'2D']+1
        for i in shadowprices_hourly['2U'].index:
            if shadowprices_hourly['2U'].at[i,SCEN] > 0:
                hours.loc[SCEN,'2U']=hours.loc[SCEN,'2U']+1
        for i in shadowprices_hourly['3U'].index:
            if shadowprices_hourly['3U'].at[i,SCEN] > 0:
                hours.loc[SCEN,'3U']=hours.loc[SCEN,'3U']+1
        for i in shadowprices_hourly['DA'].index:
            if shadowprices_hourly['DA'].at[i,SCEN] > 0:
                hours.loc[SCEN,'DA']=hours.loc[SCEN,'DA']+1     
        

        # storage level
        '''
        scenarios_dic['results_s2_B']['OutputStorageLevel'][CHP_GAS].mean()
        scenarios_dic['results_s2_B']['OutputStorageLevel'][CHP_GAS].plot()
        scenarios_dic['results_s2_B']['OutputStorageLevel'][CHP_HRD].mean()
        scenarios_dic['results_s2_B']['OutputStorageLevel'][CHP_HRD].plot()
        scenarios_dic['results_s2_B']['OutputStorageLevel'][CHP_OIL].mean()
        scenarios_dic['results_s2_B']['OutputStorageLevel'][CHP_OIL].plot()
       '''
       
# cashflow costs and total
        for u in unit_list:
            if u in inputs['units'].index:
                TMP = ds.Cashflows(inputs, results, u)
                cash_flow[u].at[SCEN,'costs'] = TMP['costs'].sum()
                cash_flow[u].at[SCEN,'heat'] = TMP['heat'].sum()
                cash_flow[u].at[SCEN,'DA'] = TMP['DA'].sum()
                cash_flow[u].at[SCEN,'2U'] = TMP['2U'].sum()
                cash_flow[u].at[SCEN,'2D'] = TMP['2D'].sum()
                cash_flow[u].at[SCEN,'3U'] = TMP['3U'].sum()
        

    for u in cash_flow:
        cash_flow[u]['total'] = cash_flow[u].sum(axis=1)

'''
calculating the costs difference


test,TMP2,test2 = ds.CostExPost(scenarios_dic['inputs_s12_A'],scenarios_dic['results_s12_A'],False)
test1 = (test.sum(axis=1) - scenarios_dic['results_s12_A']['OutputSystemCost']).abs()

'''
#technical potential
'''
CHP_pot=(-63/((234-63)/(0.081-0.054)))+0.054
P2H_pot=(-2/((721-2)/(0.081-0.054)))+0.054
CHP_P2H_pot=(-1/((188-1)/(0.081-0.054)))+0.054

tech_pot=pd.DataFrame(0,index= [0.054,0.081,0.112],columns=['CHP','P2H','CHP+P2H'])
tech_pot.loc[0.054,'CHP']=63
tech_pot.loc[0.081,'CHP']=234
tech_pot.loc[0.112,'CHP']=546
tech_pot.loc[0.054,'P2H']=2
tech_pot.loc[0.081,'P2H']=721
tech_pot.loc[0.112,'P2H']=3041
tech_pot.loc[0.054,'CHP+P2H']=1
tech_pot.loc[0.081,'CHP+P2H']=188
tech_pot.loc[0.112,'CHP+P2H']=1008
tech_pot.loc[CHP_pot,'CHP']=0
tech_pot.loc[P2H_pot,'P2H']=0
tech_pot.loc[CHP_P2H_pot,'CHP+P2H']=0
'''



    
# %% function to plot data in specified time interval (rng)
def plot_week(data, rng=index, color='k'):
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(13, 6), frameon=True)
    axes.plot(rng, data[rng], color='k')
    plt.show()
    return True

def plot_clustered_stacked(dfall, labels=None,  H="/"):
    '''
    Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
    labels is a list of the names of the dataframe, used for the legend
    title is a string for the title of the plot
    H is the hatch used for identification of the different dataframe
    '''

    n_df = len(dfall)
    n_col = len(dfall[0].columns) 
    n_ind = len(dfall[0].index)
    axe = plt.subplot(111)

    for df in dfall : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ax=axe,
                      legend=False,
                      grid=False,
                      )  # make bar plots

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                rect.set_hatch(H * int(i / n_col)) #edited part     
                rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.index, rotation = 0)
    #axe.set_title(title)

    # Add invisible data to add another legend
    n=[]        
    for i in range(n_df):
        n.append(axe.bar(0, 0, color='gray', hatch=H * i))

    l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0.5])
    if labels is not None:
        l2 = plt.legend(n, labels, loc=[1.01, 0.1]) 
    axe.add_artist(l1)
    return axe


    
# %%  defining indexes
'''
january = index[0:744]
february = index[744:1416]
march = index[1416:2160]
april = index[2160:2880]
may = index[2880:3624]
june = index[3624:4344]
july = index[4344:5088]
august = index[5088:5832]
september = index[5832:6552]
october = index[6552:7296]
november = index[7296:8016]
decemcer = index[8016:8760]

april1=index[2620:2697]
april2 = index[2596:2676]
april3 = index[2112:2208]
april4 = index[2112:2160]
october1 = index[6636:6708]
february1 = index[792:864]
july1 = index[4646:4806]
'''

# %% graphs don't forget ProRes1
if plot == True:
    
    #%% reserve demand
    plot_reservedemand = pd.DataFrame(0,index=['RES min','RES mid','RES max'],columns=['Reserve Demand'])
    plot_reservedemand.loc['RES min','Reserve Demand'] = reservedemand.at['s1','upwards']
    plot_reservedemand.loc['RES mid','Reserve Demand'] = reservedemand.at['s5','upwards']
    plot_reservedemand.loc['RES max','Reserve Demand'] = reservedemand.at['s9','upwards']
    ax = plot_reservedemand.plot.bar(figsize=(18, 6), alpha=0.4, legend = False)
    ax.set_ylabel('Reserve Demand [%]')
    ax.title.set_text('reserve demand')
    
    #%% total cost between scenarios add participation difference
    
    totalcosts = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['none','CHP','P2H','CHP+P2H'] )
    totalcosts.loc['RES min','none'] = calculations.at['s1','postsystemcostIT']
    totalcosts.loc['RES mid','none'] = calculations.at['s5','postsystemcostIT']
    totalcosts.loc['RES max','none'] = calculations.at['s9','postsystemcostIT']
    totalcosts.loc['RES min','CHP'] = calculations.at['s2_A','postsystemcostIT']
    totalcosts.loc['RES mid','CHP'] = calculations.at['s6_A','postsystemcostIT']
    totalcosts.loc['RES max','CHP'] = calculations.at['s10_A','postsystemcostIT']
    totalcosts.loc['RES min','P2H'] = calculations.at['s3_A','postsystemcostIT']
    totalcosts.loc['RES mid','P2H'] = calculations.at['s7_A','postsystemcostIT']
    totalcosts.loc['RES max','P2H'] = calculations.at['s11_A','postsystemcostIT']
    totalcosts.loc['RES min','CHP+P2H'] = calculations.at['s4_A','postsystemcostIT']
    totalcosts.loc['RES mid','CHP+P2H'] = calculations.at['s8_A','postsystemcostIT']
    totalcosts.loc['RES max','CHP+P2H'] = calculations.at['s12_A','postsystemcostIT']
    ax = totalcosts.plot.bar(figsize=(18, 6), alpha=0.2)
    ax.set_ylabel('total Italian system costs [euro]')
    ax.title.set_text('Italian system cost')#with heat slack and without DHS participation
    
    
    CHP_totalcosts = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    CHP_totalcosts.loc['RES min','no coupling'] = calculations.at['s1','postsystemcostIT']
    CHP_totalcosts.loc['RES min','no participation'] = calculations.at['s2_A','postsystemcostIT']
    CHP_totalcosts.loc['RES min','participation'] = calculations.at['s2_B','postsystemcostIT']
    CHP_totalcosts.loc['RES mid','no coupling'] = calculations.at['s5','postsystemcostIT']
    CHP_totalcosts.loc['RES mid','no participation'] = calculations.at['s6_A','postsystemcostIT']
    CHP_totalcosts.loc['RES mid','participation'] = calculations.at['s6_B','postsystemcostIT']
    CHP_totalcosts.loc['RES max','no coupling'] = calculations.at['s9','postsystemcostIT']
    CHP_totalcosts.loc['RES max','no participation'] = calculations.at['s10_A','postsystemcostIT']
    CHP_totalcosts.loc['RES max','participation'] = calculations.at['s10_B','postsystemcostIT']

    P2H_totalcosts = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    P2H_totalcosts.loc['RES min','no coupling'] = calculations.at['s1','postsystemcostIT']
    P2H_totalcosts.loc['RES min','no participation'] = calculations.at['s3_A','postsystemcostIT']
    P2H_totalcosts.loc['RES min','participation'] = calculations.at['s3_B','postsystemcostIT']
    P2H_totalcosts.loc['RES mid','no coupling'] = calculations.at['s5','postsystemcostIT']
    P2H_totalcosts.loc['RES mid','no participation'] = calculations.at['s7_A','postsystemcostIT']
    P2H_totalcosts.loc['RES mid','participation'] = calculations.at['s7_B','postsystemcostIT']
    P2H_totalcosts.loc['RES max','no coupling'] = calculations.at['s9','postsystemcostIT']
    P2H_totalcosts.loc['RES max','no participation'] = calculations.at['s11_A','postsystemcostIT']
    P2H_totalcosts.loc['RES max','participation'] = calculations.at['s11_B','postsystemcostIT']

    CHP_P2H_totalcosts = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    CHP_P2H_totalcosts.loc['RES min','no coupling'] = calculations.at['s1','postsystemcostIT']
    CHP_P2H_totalcosts.loc['RES min','no participation'] = calculations.at['s4_A','postsystemcostIT']
    CHP_P2H_totalcosts.loc['RES min','participation'] = calculations.at['s4_B','postsystemcostIT']
    CHP_P2H_totalcosts.loc['RES mid','no coupling'] = calculations.at['s5','postsystemcostIT']
    CHP_P2H_totalcosts.loc['RES mid','no participation'] = calculations.at['s8_A','postsystemcostIT']
    CHP_P2H_totalcosts.loc['RES mid','participation'] = calculations.at['s8_B','postsystemcostIT']
    CHP_P2H_totalcosts.loc['RES max','no coupling'] = calculations.at['s9','postsystemcostIT']
    CHP_P2H_totalcosts.loc['RES max','no participation'] = calculations.at['s12_A','postsystemcostIT']
    CHP_P2H_totalcosts.loc['RES max','participation'] = calculations.at['s12_B','postsystemcostIT']
    
    MAX = max([CHP_totalcosts.max().max(),P2H_totalcosts.max().max(),CHP_P2H_totalcosts.max().max()])*1.05
    
    ax = CHP_totalcosts.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('Italian system cost including heat slack for the CHP scenarios')
    
    ax = P2H_totalcosts.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('Italian system cost including heat slack for the P2H scenarios')
    
    ax = CHP_P2H_totalcosts.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('Italian system cost including heat slack for the CHP+P2H scenarios')
    
    # cost savings true participation       
    costsdifference = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['CHP_participation','P2H_participation','CHP+P2H_participation'] )
    costsdifference.loc['RES min','CHP_participation'] = (calculations.at['s2_B','postsystemcostIT']-calculations.at['s2_A','postsystemcostIT'])/calculations.at['s2_A','postsystemcostIT']*100
    costsdifference.loc['RES min','P2H_participation'] = (calculations.at['s3_B','postsystemcostIT']-calculations.at['s3_A','postsystemcostIT'])/calculations.at['s3_A','postsystemcostIT']*100
    costsdifference.loc['RES min','CHP+P2H_participation'] = (calculations.at['s4_B','postsystemcostIT']-calculations.at['s4_A','postsystemcostIT'])/calculations.at['s4_A','postsystemcostIT']*100
    costsdifference.loc['RES mid','CHP_participation'] = (calculations.at['s6_B','postsystemcostIT']-calculations.at['s6_A','postsystemcostIT'])/calculations.at['s6_A','postsystemcostIT']*100
    costsdifference.loc['RES mid','P2H_participation'] = (calculations.at['s7_B','postsystemcostIT']-calculations.at['s7_A','postsystemcostIT'])/calculations.at['s7_A','postsystemcostIT']*100
    costsdifference.loc['RES mid','CHP+P2H_participation'] = (calculations.at['s8_B','postsystemcostIT']-calculations.at['s8_A','postsystemcostIT'])/calculations.at['s8_A','postsystemcostIT']*100
    costsdifference.loc['RES max','CHP_participation'] = (calculations.at['s10_B','postsystemcostIT']-calculations.at['s10_A','postsystemcostIT'])/calculations.at['s10_A','postsystemcostIT']*100
    costsdifference.loc['RES max','P2H_participation'] = (calculations.at['s11_B','postsystemcostIT']-calculations.at['s11_A','postsystemcostIT'])/calculations.at['s11_A','postsystemcostIT']*100
    costsdifference.loc['RES max','CHP+P2H_participation'] = (calculations.at['s12_B','postsystemcostIT']-calculations.at['s12_A','postsystemcostIT'])/calculations.at['s12_A','postsystemcostIT']*100
    
    MAX = costsdifference.max().max()*1.05
    MIN = costsdifference.min().min()*1.05
    
    ax = costsdifference.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (MIN,MAX), legend = False)
    ax.set_ylabel('[%]')
    ax.legend(['CHP_participation','P2H_participation','CHP+P2H_participation'])
    ax.title.set_text('cost savings through participation (negative is a cost reduction)')
 
    
    #%% total emissions without participation
    '''
    none_emissions = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['Power','Heat Slack'])
    none_emissions.loc['RES min','Power'] = calculations.at['s1','co2_power']
    none_emissions.loc['RES min','Heat Slack'] = calculations.at['s1','co2_heatslack']
    none_emissions.loc['RES mid','Power'] = calculations.at['s5','co2_power']
    none_emissions.loc['RES mid','Heat Slack'] = calculations.at['s5','co2_heatslack']
    none_emissions.loc['RES max','Power'] = calculations.at['s9','co2_power']
    none_emissions.loc['RES max','Heat Slack'] = calculations.at['s9','co2_heatslack']
    
    CHP_emissions = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['Power','Heat Slack'])
    CHP_emissions.loc['RES min','Power'] = calculations.at['s2_A','co2_power']
    CHP_emissions.loc['RES min','Heat Slack'] = calculations.at['s2_A','co2_heatslack']
    CHP_emissions.loc['RES mid','Power'] = calculations.at['s6_A','co2_power']
    CHP_emissions.loc['RES mid','Heat Slack'] = calculations.at['s6_A','co2_heatslack']
    CHP_emissions.loc['RES max','Power'] = calculations.at['s10_A','co2_power']
    CHP_emissions.loc['RES max','Heat Slack'] = calculations.at['s10_A','co2_heatslack']
    
    P2H_emissions = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['Power','Heat Slack'])
    P2H_emissions.loc['RES min','Power'] = calculations.at['s3_A','co2_power']
    P2H_emissions.loc['RES min','Heat Slack'] = calculations.at['s3_A','co2_heatslack']
    P2H_emissions.loc['RES mid','Power'] = calculations.at['s7_A','co2_power']
    P2H_emissions.loc['RES mid','Heat Slack'] = calculations.at['s7_A','co2_heatslack']
    P2H_emissions.loc['RES max','Power'] = calculations.at['s11_A','co2_power']
    P2H_emissions.loc['RES max','Heat Slack'] = calculations.at['s11_A','co2_heatslack']
    
    CHP_P2H_emissions = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['Power','Heat Slack'])
    CHP_P2H_emissions.loc['RES min','Power'] = 0
    CHP_P2H_emissions.loc['RES min','Heat Slack'] = 0
    CHP_P2H_emissions.loc['RES mid','Power'] = 0
    CHP_P2H_emissions.loc['RES mid','Heat Slack'] = 0
    CHP_P2H_emissions.loc['RES max','Power'] = 0
    CHP_P2H_emissions.loc['RES max','Heat Slack'] = 0
        
    plot_clustered_stacked([none_emissions, CHP_emissions, P2H_emissions,CHP_P2H_emissions],['none', 'CHP', 'P2H','CHP+P2H'])
    '''
    
    plot_emissions = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['none power system','none heat slack','CHP power system','CHP heat slack','P2H power system','P2H heat slack','CHP+P2H power system','CHP+P2H heat slack'] )
    plot_emissions.loc['RES min','none power system'] = calculations.at['s1','co2_power']
    plot_emissions.loc['RES min','none heat slack'] = calculations.at['s1','co2_heatslack']
    plot_emissions.loc['RES min','CHP power system'] = calculations.at['s2_A','co2_power']
    plot_emissions.loc['RES min','CHP heat slack'] = calculations.at['s2_A','co2_heatslack']
    plot_emissions.loc['RES min','P2H power system'] = calculations.at['s3_A','co2_power']
    plot_emissions.loc['RES min','P2H heat slack'] = calculations.at['s3_A','co2_heatslack']
    plot_emissions.loc['RES min','CHP+P2H power system'] = calculations.at['s4_A','co2_power']
    plot_emissions.loc['RES min','CHP+P2H heat slack'] = calculations.at['s4_A','co2_heatslack']

    plot_emissions.loc['RES mid','none power system'] = calculations.at['s5','co2_power']
    plot_emissions.loc['RES mid','none heat slack'] = calculations.at['s5','co2_heatslack']
    plot_emissions.loc['RES mid','CHP power system'] = calculations.at['s6_A','co2_power']
    plot_emissions.loc['RES mid','CHP heat slack'] = calculations.at['s6_A','co2_heatslack']
    plot_emissions.loc['RES mid','P2H power system'] = calculations.at['s7_A','co2_power']
    plot_emissions.loc['RES mid','P2H heat slack'] = calculations.at['s7_A','co2_heatslack']
    plot_emissions.loc['RES mid','CHP+P2H power system'] = calculations.at['s8_A','co2_power']
    plot_emissions.loc['RES mid','CHP+P2H heat slack'] = calculations.at['s8_A','co2_heatslack']

    plot_emissions.loc['RES max','none power system'] = calculations.at['s9','co2_power']
    plot_emissions.loc['RES max','none heat slack'] = calculations.at['s9','co2_heatslack']
    plot_emissions.loc['RES max','CHP power system'] = calculations.at['s10_A','co2_power']
    plot_emissions.loc['RES max','CHP heat slack'] = calculations.at['s10_A','co2_heatslack']
    plot_emissions.loc['RES max','P2H power system'] = calculations.at['s11_A','co2_power']
    plot_emissions.loc['RES max','P2H heat slack'] = calculations.at['s11_A','co2_heatslack']
    plot_emissions.loc['RES max','CHP+P2H power system'] = calculations.at['s12_A','co2_power']
    plot_emissions.loc['RES max','CHP+P2H heat slack'] = calculations.at['s12_A','co2_heatslack']

    MAX = (plot_emissions.at['RES min','CHP power system']+plot_emissions.at['RES min','CHP heat slack'])*1.05
    fig, ax = plt.subplots()
    ax = plot_emissions[['none power system','none heat slack']].plot.bar(figsize=(18, 6),stacked = True, alpha=0.8, width=0.1, position=2,colormap="rainbow",ax=ax,ylim = (0,MAX), legend = True)
    ax = plot_emissions[['CHP power system','CHP heat slack']].plot.bar(figsize=(18, 6),stacked = True,alpha=0.6, width=0.1, position=1,colormap="rainbow",ax=ax,ylim = (0,MAX), legend = True)
    ax = plot_emissions[['P2H power system','P2H heat slack']].plot.bar(figsize=(18, 6),stacked = True,alpha=0.4, width=0.1, position=0,colormap="rainbow",ax=ax,ylim = (0,MAX), legend = True)
    ax = plot_emissions[['CHP+P2H power system','CHP+P2H heat slack']].plot.bar(figsize=(18, 6),stacked = True,alpha=0.1, width=0.1, position=-1,colormap="rainbow",ax=ax,ylim = (0,MAX), legend = True)
    ax.set_ylabel('emissions [ton co2]')
    ax.title.set_text('emissions produced in Italy')
    
    CHP_emissions = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    CHP_emissions.loc['RES min','no coupling'] = calculations.at['s1','co2_total']
    CHP_emissions.loc['RES min','no participation'] = calculations.at['s2_A','co2_total']
    CHP_emissions.loc['RES min','participation'] = calculations.at['s2_B','co2_total']
    CHP_emissions.loc['RES mid','no coupling'] = calculations.at['s5','co2_total']
    CHP_emissions.loc['RES mid','no participation'] = calculations.at['s6_A','co2_total']
    CHP_emissions.loc['RES mid','participation'] = calculations.at['s6_B','co2_total']
    CHP_emissions.loc['RES max','no coupling'] = calculations.at['s9','co2_total']
    CHP_emissions.loc['RES max','no participation'] = calculations.at['s10_A','co2_total']
    CHP_emissions.loc['RES max','participation'] = calculations.at['s10_B','co2_total']

    P2H_emissions = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    P2H_emissions.loc['RES min','no coupling'] = calculations.at['s1','co2_total']
    P2H_emissions.loc['RES min','no participation'] = calculations.at['s3_A','co2_total']
    P2H_emissions.loc['RES min','participation'] = calculations.at['s3_B','co2_total']
    P2H_emissions.loc['RES mid','no coupling'] = calculations.at['s5','co2_total']
    P2H_emissions.loc['RES mid','no participation'] = calculations.at['s7_A','co2_total']
    P2H_emissions.loc['RES mid','participation'] = calculations.at['s7_B','co2_total']
    P2H_emissions.loc['RES max','no coupling'] = calculations.at['s9','co2_total']
    P2H_emissions.loc['RES max','no participation'] = calculations.at['s11_A','co2_total']
    P2H_emissions.loc['RES max','participation'] = calculations.at['s11_B','co2_total']

    CHP_P2H_emissions = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    CHP_P2H_emissions.loc['RES min','no coupling'] = calculations.at['s1','co2_total']
    CHP_P2H_emissions.loc['RES min','no participation'] = calculations.at['s4_A','co2_total']
    CHP_P2H_emissions.loc['RES min','participation'] = calculations.at['s4_B','co2_total']
    CHP_P2H_emissions.loc['RES mid','no coupling'] = calculations.at['s5','co2_total']
    CHP_P2H_emissions.loc['RES mid','no participation'] = calculations.at['s8_A','co2_total']
    CHP_P2H_emissions.loc['RES mid','participation'] = calculations.at['s8_B','co2_total']
    CHP_P2H_emissions.loc['RES max','no coupling'] = calculations.at['s9','co2_total']
    CHP_P2H_emissions.loc['RES max','no participation'] = calculations.at['s12_A','co2_total']
    CHP_P2H_emissions.loc['RES max','participation'] = calculations.at['s12_B','co2_total']
    
    MAX = max([CHP_emissions.max().max(),P2H_emissions.max().max(),CHP_P2H_emissions.max().max()])*1.05
    
    ax = CHP_emissions.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('Italian system emissions including heat slack for the CHP scenarios')
    
    ax = P2H_emissions.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('Italian system emissions including heat slack for the P2H scenarios')
    
    ax = CHP_P2H_emissions.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('Italian system emissions including heat slack for the CHP+P2H scenarios')

    emission_difference = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['CHP_participation','P2H_participation','CHP+P2H_participation'] )
    emission_difference.loc['RES min','CHP_participation'] = (calculations.at['s2_B','co2_total']-calculations.at['s2_A','co2_total'])/calculations.at['s2_A','co2_total']*100
    emission_difference.loc['RES min','P2H_participation'] = (calculations.at['s3_B','co2_total']-calculations.at['s3_A','co2_total'])/calculations.at['s3_A','co2_total']*100
    emission_difference.loc['RES min','CHP+P2H_participation'] = (calculations.at['s4_B','co2_total']-calculations.at['s4_A','co2_total'])/calculations.at['s4_A','co2_total']*100
    emission_difference.loc['RES mid','CHP_participation'] = (calculations.at['s6_B','co2_total']-calculations.at['s6_A','co2_total'])/calculations.at['s6_A','co2_total']*100
    emission_difference.loc['RES mid','P2H_participation'] = (calculations.at['s7_B','co2_total']-calculations.at['s7_A','co2_total'])/calculations.at['s7_A','co2_total']*100
    emission_difference.loc['RES mid','CHP+P2H_participation'] = (calculations.at['s8_B','co2_total']-calculations.at['s8_A','co2_total'])/calculations.at['s8_A','co2_total']*100
    emission_difference.loc['RES max','CHP_participation'] = (calculations.at['s10_B','co2_total']-calculations.at['s10_A','co2_total'])/calculations.at['s10_A','co2_total']*100
    emission_difference.loc['RES max','P2H_participation'] = (calculations.at['s11_B','co2_total']-calculations.at['s11_A','co2_total'])/calculations.at['s11_A','co2_total']*100
    emission_difference.loc['RES max','CHP+P2H_participation'] = (calculations.at['s12_B','co2_total']-calculations.at['s12_A','co2_total'])/calculations.at['s12_A','co2_total']*100

    MAX = emission_difference.max().max()*1.05
    MIN = emission_difference.min().min()*1.05
    
    ax = emission_difference.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (MIN,MAX), legend = False)
    ax.set_ylabel('[%]')
    ax.legend(['CHP_participation','P2H_participation','CHP+P2H_participation'])
    ax.title.set_text('emissions saving through participation (negative is a decreasement)')

    
    
    #%% curtailment @summed up

    CHP_curtailment = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    CHP_curtailment.loc['RES min','no coupling'] = curtailment.at['s1','sum']
    CHP_curtailment.loc['RES min','no participation'] = curtailment.at['s2_A','sum']
    CHP_curtailment.loc['RES min','participation'] = curtailment.at['s2_B','sum']
    CHP_curtailment.loc['RES mid','no coupling'] = curtailment.at['s5','sum']
    CHP_curtailment.loc['RES mid','no participation'] = curtailment.at['s6_A','sum']
    CHP_curtailment.loc['RES mid','participation'] = curtailment.at['s6_B','sum']
    CHP_curtailment.loc['RES max','no coupling'] = curtailment.at['s9','sum']
    CHP_curtailment.loc['RES max','no participation'] = curtailment.at['s10_A','sum']
    CHP_curtailment.loc['RES max','participation'] = curtailment.at['s10_B','sum']

    P2H_curtailment = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    P2H_curtailment.loc['RES min','no coupling'] = curtailment.at['s1','sum']
    P2H_curtailment.loc['RES min','no participation'] = curtailment.at['s3_A','sum']
    P2H_curtailment.loc['RES min','participation'] = curtailment.at['s3_B','sum']
    P2H_curtailment.loc['RES mid','no coupling'] = curtailment.at['s5','sum']
    P2H_curtailment.loc['RES mid','no participation'] = curtailment.at['s7_A','sum']
    P2H_curtailment.loc['RES mid','participation'] = curtailment.at['s7_B','sum']
    P2H_curtailment.loc['RES max','no coupling'] = curtailment.at['s9','sum']
    P2H_curtailment.loc['RES max','no participation'] = curtailment.at['s11_A','sum']
    P2H_curtailment.loc['RES max','participation'] = curtailment.at['s11_B','sum']

    CHP_P2H_curtailment = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    CHP_P2H_curtailment.loc['RES min','no coupling'] = curtailment.at['s1','sum']
    CHP_P2H_curtailment.loc['RES min','no participation'] = curtailment.at['s4_A','sum']
    CHP_P2H_curtailment.loc['RES min','participation'] = curtailment.at['s4_B','sum']
    CHP_P2H_curtailment.loc['RES mid','no coupling'] = curtailment.at['s5','sum']
    CHP_P2H_curtailment.loc['RES mid','no participation'] = curtailment.at['s8_A','sum']
    CHP_P2H_curtailment.loc['RES mid','participation'] = curtailment.at['s8_B','sum']
    CHP_P2H_curtailment.loc['RES max','no coupling'] = curtailment.at['s9','sum']
    CHP_P2H_curtailment.loc['RES max','no participation'] = curtailment.at['s12_A','sum']
    CHP_P2H_curtailment.loc['RES max','participation'] = curtailment.at['s12_B','sum']
    
    MAX = max([CHP_curtailment.max().max(),P2H_curtailment.max().max(),CHP_P2H_curtailment.max().max()])*1.05
    
    ax = CHP_curtailment.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('curtailment for CHP scenarios')
    
    ax = P2H_curtailment.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('curtailment for P2H scenarios')
    
    ax = CHP_P2H_curtailment.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('curtailment for CHP+P2H scenarios')

    
    
    
    #%% load shedding (total)
    '''
    CHP_loadshedding = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    CHP_loadshedding.loc['RES min','no coupling'] = load_shedding.at['s1','sum']
    CHP_loadshedding.loc['RES min','no participation'] = load_shedding.at['s2_A','sum']
    CHP_loadshedding.loc['RES min','participation'] = load_shedding.at['s2_B','sum']
    CHP_loadshedding.loc['RES mid','no coupling'] = load_shedding.at['s5','sum']
    CHP_loadshedding.loc['RES mid','no participation'] = load_shedding.at['s6_A','sum']
    CHP_loadshedding.loc['RES mid','participation'] = load_shedding.at['s6_B','sum']
    CHP_loadshedding.loc['RES max','no coupling'] = load_shedding.at['s9','sum']
    CHP_loadshedding.loc['RES max','no participation'] = load_shedding.at['s10_A','sum']
    CHP_loadshedding.loc['RES max','participation'] = load_shedding.at['s10_B','sum']

    P2H_loadshedding = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    P2H_loadshedding.loc['RES min','no coupling'] = load_shedding.at['s1','sum']
    P2H_loadshedding.loc['RES min','no participation'] = load_shedding.at['s3_A','sum']
    P2H_loadshedding.loc['RES min','participation'] = load_shedding.at['s3_B','sum']
    P2H_loadshedding.loc['RES mid','no coupling'] = load_shedding.at['s5','sum']
    P2H_loadshedding.loc['RES mid','no participation'] = load_shedding.at['s7_A','sum']
    P2H_loadshedding.loc['RES mid','participation'] = load_shedding.at['s7_B','sum']
    P2H_loadshedding.loc['RES max','no coupling'] = load_shedding.at['s9','sum']
    P2H_loadshedding.loc['RES max','no participation'] = load_shedding.at['s11_A','sum']
    P2H_loadshedding.loc['RES max','participation'] = load_shedding.at['s11_B','sum']

    CHP_P2H_loadshedding = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    CHP_P2H_loadshedding.loc['RES min','no coupling'] = load_shedding.at['s1','sum']
    CHP_P2H_loadshedding.loc['RES min','no participation'] = load_shedding.at['s4_A','sum']
    CHP_P2H_loadshedding.loc['RES min','participation'] = load_shedding.at['s4_B','sum']
    CHP_P2H_loadshedding.loc['RES mid','no coupling'] = load_shedding.at['s5','sum']
    CHP_P2H_loadshedding.loc['RES mid','no participation'] = load_shedding.at['s8_A','sum']
    CHP_P2H_loadshedding.loc['RES mid','participation'] = load_shedding.at['s8_B','sum']
    CHP_P2H_loadshedding.loc['RES max','no coupling'] = load_shedding.at['s9','sum']
    CHP_P2H_loadshedding.loc['RES max','no participation'] = load_shedding.at['s12_A','sum']
    CHP_P2H_loadshedding.loc['RES max','participation'] = load_shedding.at['s12_B','sum']
    
    MAX = max([CHP_loadshedding.max().max(),P2H_loadshedding.max().max(),CHP_P2H_loadshedding.max().max()])*1.05
    
    ax = CHP_loadshedding.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('loadshedding [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('total load shedding in MWh for CHP scenarios')
    
    ax = P2H_loadshedding.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('loadshedding [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('total load shedding in MWh for P2H scenarios')
    
    ax = CHP_P2H_loadshedding.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('loadshedding [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('total load shedding in MWh for CHP+P2H scenarios')
    
    #%% load shedding (amount)
    
    CHP_loadshedding = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    CHP_loadshedding.loc['RES min','no coupling'] = load_shedding.at['s1','amount']
    CHP_loadshedding.loc['RES min','no participation'] = load_shedding.at['s2_A','amount']
    CHP_loadshedding.loc['RES min','participation'] = load_shedding.at['s2_B','amount']
    CHP_loadshedding.loc['RES mid','no coupling'] = load_shedding.at['s5','amount']
    CHP_loadshedding.loc['RES mid','no participation'] = load_shedding.at['s6_A','amount']
    CHP_loadshedding.loc['RES mid','participation'] = load_shedding.at['s6_B','amount']
    CHP_loadshedding.loc['RES max','no coupling'] = load_shedding.at['s9','amount']
    CHP_loadshedding.loc['RES max','no participation'] = load_shedding.at['s10_A','amount']
    CHP_loadshedding.loc['RES max','participation'] = load_shedding.at['s10_B','amount']

    P2H_loadshedding = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    P2H_loadshedding.loc['RES min','no coupling'] = load_shedding.at['s1','amount']
    P2H_loadshedding.loc['RES min','no participation'] = load_shedding.at['s3_A','amount']
    P2H_loadshedding.loc['RES min','participation'] = load_shedding.at['s3_B','amount']
    P2H_loadshedding.loc['RES mid','no coupling'] = load_shedding.at['s5','amount']
    P2H_loadshedding.loc['RES mid','no participation'] = load_shedding.at['s7_A','amount']
    P2H_loadshedding.loc['RES mid','participation'] = load_shedding.at['s7_B','amount']
    P2H_loadshedding.loc['RES max','no coupling'] = load_shedding.at['s9','amount']
    P2H_loadshedding.loc['RES max','no participation'] = load_shedding.at['s11_A','amount']
    P2H_loadshedding.loc['RES max','participation'] = load_shedding.at['s11_B','amount']

    CHP_P2H_loadshedding = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['no coupling','no participation','participation'] )
    CHP_P2H_loadshedding.loc['RES min','no coupling'] = load_shedding.at['s1','amount']
    CHP_P2H_loadshedding.loc['RES min','no participation'] = load_shedding.at['s4_A','amount']
    CHP_P2H_loadshedding.loc['RES min','participation'] = load_shedding.at['s4_B','amount']
    CHP_P2H_loadshedding.loc['RES mid','no coupling'] = load_shedding.at['s5','amount']
    CHP_P2H_loadshedding.loc['RES mid','no participation'] = load_shedding.at['s8_A','amount']
    CHP_P2H_loadshedding.loc['RES mid','participation'] = load_shedding.at['s8_B','amount']
    CHP_P2H_loadshedding.loc['RES max','no coupling'] = load_shedding.at['s9','amount']
    CHP_P2H_loadshedding.loc['RES max','no participation'] = load_shedding.at['s12_A','amount']
    CHP_P2H_loadshedding.loc['RES max','participation'] = load_shedding.at['s12_B','amount']
    
    MAX = max([CHP_loadshedding.max().max(),P2H_loadshedding.max().max(),CHP_P2H_loadshedding.max().max()])*1.05
    
    ax = CHP_loadshedding.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('loadshedding [h]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('total amount of load shedding in h for CHP scenarios')
   
    ax = P2H_loadshedding.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('loadshedding [h]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('total amount of load shedding in h for P2H scenarios')
    
    ax = CHP_P2H_loadshedding.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('loadshedding [h]')
    ax.legend(['no coupling','no participation','participation'])
    ax.title.set_text('total amount of load shedding in h for CHP+P2H scenarios')
    '''
#%% availability
    '''
    # upwards
    MAX = 50
    plot_availability = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP_GAS','CHP_HRD','CHP_OIL','P2H','CHP+P2H'] )
    plot_availability.loc['RES min','CHP_GAS'] = availability_total['Up'].at[CHP_GAS,'s2_B']
    plot_availability.loc['RES min','CHP_HRD'] = availability_total['Up'].at[CHP_HRD,'s2_B']
    plot_availability.loc['RES min','CHP_OIL'] = availability_total['Up'].at[CHP_OIL,'s2_B']
    plot_availability.loc['RES min','P2H'] = availability_total['Up'].at[P2H_unit,'s3_B']
    plot_availability.loc['RES min','CHP+P2H_P2H'] = availability_total['Up'].at[P2H_unit,'s4_B']
    plot_availability.loc['RES min','CHP+P2H_GAS'] = availability_total['Up'].at[CHP_GAS,'s4_B']
    plot_availability.loc['RES min','CHP+P2H_HRD'] = availability_total['Up'].at[CHP_HRD,'s4_B']
    plot_availability.loc['RES min','CHP+P2H_OIL'] = availability_total['Up'].at[CHP_OIL,'s4_B']
    plot_availability.loc['RES mid','CHP_GAS'] = availability_total['Up'].at[CHP_GAS,'s6_B']
    plot_availability.loc['RES mid','CHP_HRD'] = availability_total['Up'].at[CHP_HRD,'s6_B']
    plot_availability.loc['RES mid','CHP_OIL'] = availability_total['Up'].at[CHP_OIL,'s6_B']
    plot_availability.loc['RES mid','P2H'] = availability_total['Up'].at[P2H_unit,'s7_B']
    plot_availability.loc['RES mid','CHP+P2H_P2H'] = availability_total['Up'].at[P2H_unit,'s8_B']
    plot_availability.loc['RES mid','CHP+P2H_GAS'] = availability_total['Up'].at[CHP_GAS,'s8_B']
    plot_availability.loc['RES mid','CHP+P2H_HRD'] = availability_total['Up'].at[CHP_HRD,'s8_B']
    plot_availability.loc['RES mid','CHP+P2H_OIL'] = availability_total['Up'].at[CHP_OIL,'s8_B']
    plot_availability.loc['RES max','CHP_GAS'] = availability_total['Up'].at[CHP_GAS,'s10_B']
    plot_availability.loc['RES max','CHP_HRD'] = availability_total['Up'].at[CHP_HRD,'s10_B']
    plot_availability.loc['RES max','CHP_OIL'] = availability_total['Up'].at[CHP_OIL,'s10_B']
    plot_availability.loc['RES max','P2H'] = availability_total['Up'].at[P2H_unit,'s11_B']
    plot_availability.loc['RES max','CHP+P2H_P2H'] = availability_total['Up'].at[P2H_unit,'s12_B']
    plot_availability.loc['RES max','CHP+P2H_GAS'] = availability_total['Up'].at[CHP_GAS,'s12_B']
    plot_availability.loc['RES max','CHP+P2H_HRD'] = availability_total['Up'].at[CHP_HRD,'s12_B']
    plot_availability.loc['RES max','CHP+P2H_OIL'] = availability_total['Up'].at[CHP_OIL,'s12_B']

    ax = plot_availability[['CHP_GAS','CHP_HRD','CHP_OIL']].plot.bar(figsize=(18, 6),stacked = True, alpha=0.2, width=0.1, position=1.5,colormap="rainbow",ylim = (0,MAX), legend = True)
    ax = plot_availability[['CHP+P2H_P2H','CHP+P2H_GAS','CHP+P2H_HRD','CHP+P2H_OIL']].plot.bar(figsize=(18, 6),stacked = True,alpha=0.2, width=0.1, position=-0.5,colormap="rainbow",ylim = (0,MAX), legend = True)
    ax = plot_availability['P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=0.5,ylim = (0,MAX), legend = True)
    ax.set_ylabel('availability [%]')
    ax.title.set_text('upwards availability (2U+3U) [%]')
        
    #downwards
    MAX = 100
    plot_availability = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP_GAS','CHP_HRD','CHP_OIL','P2H','CHP+P2H'] )
    plot_availability.loc['RES min','CHP_GAS'] = availability_total['Down'].at[CHP_GAS,'s2_B']
    plot_availability.loc['RES min','CHP_HRD'] = availability_total['Down'].at[CHP_HRD,'s2_B']
    plot_availability.loc['RES min','CHP_OIL'] = availability_total['Down'].at[CHP_OIL,'s2_B']
    plot_availability.loc['RES min','P2H'] = availability_total['Down'].at[P2H_unit,'s3_B']
    plot_availability.loc['RES min','CHP+P2H_P2H'] = availability_total['Down'].at[P2H_unit,'s4_B']
    plot_availability.loc['RES min','CHP+P2H_GAS'] = availability_total['Down'].at[CHP_GAS,'s4_B']
    plot_availability.loc['RES min','CHP+P2H_HRD'] = availability_total['Down'].at[CHP_HRD,'s4_B']
    plot_availability.loc['RES min','CHP+P2H_OIL'] = availability_total['Down'].at[CHP_OIL,'s4_B']
    plot_availability.loc['RES mid','CHP_GAS'] = availability_total['Down'].at[CHP_GAS,'s6_B']
    plot_availability.loc['RES mid','CHP_HRD'] = availability_total['Down'].at[CHP_HRD,'s6_B']
    plot_availability.loc['RES mid','CHP_OIL'] = availability_total['Down'].at[CHP_OIL,'s6_B']
    plot_availability.loc['RES mid','P2H'] = availability_total['Down'].at[P2H_unit,'s7_B']
    plot_availability.loc['RES mid','CHP+P2H_P2H'] = availability_total['Down'].at[P2H_unit,'s8_B']
    plot_availability.loc['RES mid','CHP+P2H_GAS'] = availability_total['Down'].at[CHP_GAS,'s8_B']
    plot_availability.loc['RES mid','CHP+P2H_HRD'] = availability_total['Down'].at[CHP_HRD,'s8_B']
    plot_availability.loc['RES mid','CHP+P2H_OIL'] = availability_total['Down'].at[CHP_OIL,'s8_B']
    plot_availability.loc['RES max','CHP_GAS'] = availability_total['Down'].at[CHP_GAS,'s10_B']
    plot_availability.loc['RES max','CHP_HRD'] = availability_total['Down'].at[CHP_HRD,'s10_B']
    plot_availability.loc['RES max','CHP_OIL'] = availability_total['Down'].at[CHP_OIL,'s10_B']
    plot_availability.loc['RES max','P2H'] = availability_total['Down'].at[P2H_unit,'s11_B']
    plot_availability.loc['RES max','CHP+P2H_P2H'] = availability_total['Down'].at[P2H_unit,'s12_B']
    plot_availability.loc['RES max','CHP+P2H_GAS'] = availability_total['Down'].at[CHP_GAS,'s12_B']
    plot_availability.loc['RES max','CHP+P2H_HRD'] = availability_total['Down'].at[CHP_HRD,'s12_B']
    plot_availability.loc['RES max','CHP+P2H_OIL'] = availability_total['Down'].at[CHP_OIL,'s12_B']

    ax = plot_availability[['CHP_GAS','CHP_HRD','CHP_OIL'],].plot.bar(figsize=(18, 6),stacked = True, alpha=0.2, width=0.1, position=1.5,colormap="rainbow",ylim = (0,MAX), legend = True)
    ax = plot_availability['P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=0.5,ylim = (0,MAX), legend = True)
    ax = plot_availability['CHP+P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=-0.5,ylim = (0,MAX), legend = False)
    ax.set_ylabel('availability [%]')
    ax.title.set_text('downwards availability [%]')
    
    # spinning upwards reserves
    MAX = 300*1.05
    plot_availability = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP_GAS','CHP_HRD','CHP_OIL','P2H','CHP+P2H'] )
    plot_availability.loc['RES min','CHP_GAS'] = availability_total['2U'].at[CHP_GAS,'s2_B']
    plot_availability.loc['RES min','CHP_HRD'] = availability_total['2U'].at[CHP_HRD,'s2_B']
    plot_availability.loc['RES min','CHP_OIL'] = availability_total['2U'].at[CHP_OIL,'s2_B']
    plot_availability.loc['RES min','P2H'] = availability_total['2U'].at[P2H_unit,'s3_B']
    plot_availability.loc['RES min','CHP+P2H'] = 0
    plot_availability.loc['RES mid','CHP_GAS'] = availability_total['2U'].at[CHP_GAS,'s6_B']
    plot_availability.loc['RES mid','CHP_HRD'] = availability_total['2U'].at[CHP_HRD,'s6_B']
    plot_availability.loc['RES mid','CHP_OIL'] = availability_total['2U'].at[CHP_OIL,'s6_B']
    plot_availability.loc['RES mid','P2H'] = availability_total['2U'].at[P2H_unit,'s7_B']
    plot_availability.loc['RES mid','CHP+P2H'] = 0
    plot_availability.loc['RES max','CHP_GAS'] = availability_total['2U'].at[CHP_GAS,'s10_B']
    plot_availability.loc['RES max','CHP_HRD'] = availability_total['2U'].at[CHP_HRD,'s10_B']
    plot_availability.loc['RES max','CHP_OIL'] = availability_total['2U'].at[CHP_OIL,'s10_B']
    plot_availability.loc['RES max','P2H'] = availability_total['2U'].at[P2H_unit,'s11_B']
    plot_availability.loc['RES max','CHP+P2H'] = 0

    ax = plot_availability[['CHP_GAS','CHP_HRD','CHP_OIL']].plot.bar(figsize=(18, 6),stacked = True, alpha=0.2, width=0.1, position=1.5,colormap="rainbow",ylim = (0,MAX), legend = True)
    ax = plot_availability['P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=0.5,ylim = (0,MAX), legend = True)
    ax = plot_availability['CHP+P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=-0.5,ylim = (0,MAX), legend = False)
    ax.set_ylabel('availability [%]')
    
    # non spinning reserves
    MAX = 100*1.05
    plot_availability = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP_GAS','CHP_HRD','CHP_OIL','P2H','CHP+P2H'] )
    plot_availability.loc['RES min','CHP_GAS'] = availability_total['3U'].at[CHP_GAS,'s2_B']
    plot_availability.loc['RES min','CHP_HRD'] = availability_total['3U'].at[CHP_HRD,'s2_B']
    plot_availability.loc['RES min','CHP_OIL'] = availability_total['3U'].at[CHP_OIL,'s2_B']
    plot_availability.loc['RES min','P2H'] = availability_total['3U'].at[P2H_unit,'s3_B']
    plot_availability.loc['RES min','CHP+P2H'] = 0
    plot_availability.loc['RES mid','CHP_GAS'] = availability_total['3U'].at[CHP_GAS,'s6_B']
    plot_availability.loc['RES mid','CHP_HRD'] = availability_total['3U'].at[CHP_HRD,'s6_B']
    plot_availability.loc['RES mid','CHP_OIL'] = availability_total['3U'].at[CHP_OIL,'s6_B']
    plot_availability.loc['RES mid','P2H'] = availability_total['3U'].at[P2H_unit,'s7_B']
    plot_availability.loc['RES mid','CHP+P2H'] = 0
    plot_availability.loc['RES max','CHP_GAS'] = availability_total['3U'].at[CHP_GAS,'s10_B']
    plot_availability.loc['RES max','CHP_HRD'] = availability_total['3U'].at[CHP_HRD,'s10_B']
    plot_availability.loc['RES max','CHP_OIL'] = availability_total['3U'].at[CHP_OIL,'s10_B']
    plot_availability.loc['RES max','P2H'] = availability_total['3U'].at[P2H_unit,'s11_B']
    plot_availability.loc['RES max','CHP+P2H'] = 0

    ax = plot_availability[['CHP_GAS','CHP_HRD','CHP_OIL']].plot.bar(figsize=(18, 6),stacked = True, alpha=0.2, width=0.1, position=1.5,colormap="rainbow",ylim = (0,MAX), legend = True)
    ax = plot_availability['P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=0.5,ylim = (0,MAX), legend = True)
    ax = plot_availability['CHP+P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=-0.5,ylim = (0,MAX), legend = False)
    ax.set_ylabel('availability [%]')

    # yearly availability up
    MAX = 30000000
    plot_yearly_availability = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP_GAS','CHP_HRD','CHP_OIL','P2H'] )
    plot_yearly_availability.loc['RES min','CHP_GAS'] = yearly_availability['Up'].at[CHP_GAS,'s2_B']
    plot_yearly_availability.loc['RES min','CHP_HRD'] = yearly_availability['Up'].at[CHP_HRD,'s2_B']
    plot_yearly_availability.loc['RES min','CHP_OIL'] = yearly_availability['Up'].at[CHP_OIL,'s2_B']
    plot_yearly_availability.loc['RES min','P2H'] = yearly_availability['Up'].at[P2H_unit,'s3_B']
    plot_yearly_availability.loc['RES mid','CHP_GAS'] = yearly_availability['Up'].at[CHP_GAS,'s6_B']
    plot_yearly_availability.loc['RES mid','CHP_HRD'] = yearly_availability['Up'].at[CHP_HRD,'s6_B']
    plot_yearly_availability.loc['RES mid','CHP_OIL'] = yearly_availability['Up'].at[CHP_OIL,'s6_B']
    plot_yearly_availability.loc['RES mid','P2H'] = yearly_availability['Up'].at[P2H_unit,'s7_B']
    plot_yearly_availability.loc['RES max','CHP_GAS'] = yearly_availability['Up'].at[CHP_GAS,'s10_B']
    plot_yearly_availability.loc['RES max','CHP_HRD'] = yearly_availability['Up'].at[CHP_HRD,'s10_B']
    plot_yearly_availability.loc['RES max','CHP_OIL'] = yearly_availability['Up'].at[CHP_OIL,'s10_B']
    plot_yearly_availability.loc['RES max','P2H'] = yearly_availability['Up'].at[P2H_unit,'s11_B']
    
    ax = plot_yearly_availability[['CHP_GAS','CHP_HRD','CHP_OIL']].plot.bar(figsize=(18, 6),stacked = True, alpha=0.2, width=0.1, position=-0.5,colormap="rainbow", legend = True)
    ax = plot_yearly_availability['P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=0.5, ylim = (0,MAX),legend = True)
    ax.title.set_text('total sum 2U availability')

    # yearly availability down
    MAX = 3000000
    plot_yearly_availability = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP','P2H','CHP+P2H','P2H+CHP'] )
    plot_yearly_availability.loc['RES min','CHP'] = yearly_availability['Down'].at[CHP_GAS,'s2_B']+yearly_availability['Down'].at[CHP_HRD,'s2_B']+yearly_availability['Down'].at[CHP_OIL,'s2_B']
    plot_yearly_availability.loc['RES min','P2H'] = yearly_availability['Down'].at[P2H_unit,'s3_B']
    plot_yearly_availability.loc['RES min','CHP+P2H'] = yearly_availability['Down'].at[CHP_GAS,'s4_B']+yearly_availability['Down'].at[CHP_HRD,'s4_B']+yearly_availability['Down'].at[CHP_OIL,'s4_B']
    plot_yearly_availability.loc['RES min','P2H+CHP'] = yearly_availability['Down'].at[P2H_unit,'s4_B']
    plot_yearly_availability.loc['RES mid','CHP'] = yearly_availability['Down'].at[CHP_GAS,'s6_B']+yearly_availability['Down'].at[CHP_HRD,'s6_B']+yearly_availability['Down'].at[CHP_OIL,'s6_B']
    plot_yearly_availability.loc['RES mid','P2H'] = yearly_availability['Down'].at[P2H_unit,'s7_B']
    plot_yearly_availability.loc['RES mid','CHP+P2H'] = yearly_availability['Down'].at[CHP_GAS,'s8_B']+yearly_availability['Down'].at[CHP_HRD,'s8_B']+yearly_availability['Down'].at[CHP_OIL,'s8_B']
    plot_yearly_availability.loc['RES mid','P2H+CHP'] = yearly_availability['Down'].at[P2H_unit,'s8_B']
    plot_yearly_availability.loc['RES max','CHP'] = yearly_availability['Down'].at[CHP_GAS,'s10_B']+yearly_availability['Down'].at[CHP_HRD,'s10_B']+yearly_availability['Down'].at[CHP_OIL,'s10_B']
    plot_yearly_availability.loc['RES max','P2H'] = yearly_availability['Down'].at[P2H_unit,'s11_B']
    plot_yearly_availability.loc['RES max','CHP+P2H'] = yearly_availability['Down'].at[CHP_GAS,'s12_B']+yearly_availability['Down'].at[CHP_HRD,'s12_B']+yearly_availability['Down'].at[CHP_OIL,'s12_B']
    plot_yearly_availability.loc['RES max','P2H+CHP'] = yearly_availability['Down'].at[P2H_unit,'s12_B']
    
    ax = plot_yearly_availability[['CHP+P2H','P2H+CHP']].plot.bar(figsize=(18, 6),stacked = True,alpha=0.8, width=0.1, position=-0.5,color=['gray','#b93c46ff'],ylim = (0,MAX), legend = False)
    ax = plot_yearly_availability['CHP'].plot.bar(figsize=(18, 6),stacked = True, alpha=0.8, width=0.1, position=1.5,color=['#b93c46ff'],ylim = (0,MAX), legend = True)
    ax = plot_yearly_availability['P2H'].plot.bar(figsize=(18, 6),stacked = True,alpha=0.8, width=0.1, position=0.5,color=['gray'],ylim = (0,MAX), legend = True)
    ax.title.set_text('total sum availability downwards')
    
    #'#b93c46ff','#b93c46b2','#b93c4666'


    #downwards
    MAX = 35*1.05
    plot_availability = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP','P2H','CHP+P2H','P2H+CHP'] )
    plot_availability.loc['RES min','CHP'] = availability_total['Down'].at[CHP_GAS,'s2_B']+availability_total['Down'].at[CHP_HRD,'s2_B']+availability_total['Down'].at[CHP_OIL,'s2_B']
    plot_availability.loc['RES min','P2H'] = availability_total['Down'].at[P2H_unit,'s3_B']
    plot_availability.loc['RES min','CHP+P2H'] = availability_total['Down'].at[CHP_GAS,'s4_B']+availability_total['Down'].at[CHP_HRD,'s4_B']+availability_total['Down'].at[CHP_OIL,'s4_B']
    plot_availability.loc['RES min','P2H+CHP'] = availability_total['Down'].at[P2H_unit,'s4_B']
    plot_availability.loc['RES mid','CHP'] = availability_total['Down'].at[CHP_GAS,'s6_B']+availability_total['Down'].at[CHP_HRD,'s6_B']+availability_total['Down'].at[CHP_OIL,'s6_B']
    plot_availability.loc['RES mid','P2H'] = availability_total['Down'].at[P2H_unit,'s7_B']
    plot_availability.loc['RES mid','CHP+P2H'] = availability_total['Down'].at[CHP_GAS,'s8_B']+availability_total['Down'].at[CHP_HRD,'s8_B']+availability_total['Down'].at[CHP_OIL,'s8_B']
    plot_availability.loc['RES mid','P2H+CHP'] = availability_total['Down'].at[P2H_unit,'s8_B']
    plot_availability.loc['RES max','CHP'] = availability_total['Down'].at[CHP_GAS,'s10_B']+availability_total['Down'].at[CHP_HRD,'s10_B']+availability_total['Down'].at[CHP_OIL,'s10_B']
    plot_availability.loc['RES max','P2H'] = availability_total['Down'].at[P2H_unit,'s11_B']
    plot_availability.loc['RES max','CHP+P2H'] = availability_total['Down'].at[CHP_GAS,'s12_B']+availability_total['Down'].at[CHP_HRD,'s12_B']+availability_total['Down'].at[CHP_OIL,'s12_B']
    plot_availability.loc['RES max','P2H+CHP'] = availability_total['Down'].at[P2H_unit,'s12_B']

    ax = plot_availability[['CHP+P2H','P2H+CHP']].plot.bar(figsize=(18, 6),stacked = True,alpha=0.8, width=0.1, position=-0.5,color=['gray','#b93c46ff'],ylim = (0,MAX), legend = False)
    ax = plot_availability['CHP'].plot.bar(figsize=(18, 6),stacked = True, alpha=0.8, width=0.1, position=1.5,color=['#b93c46ff'],ylim = (0,MAX), legend = True)
    ax = plot_availability['P2H'].plot.bar(figsize=(18, 6),stacked = True,alpha=0.8, width=0.1, position=0.5,color=['gray'],ylim = (0,MAX), legend = True)
    ax.set_ylabel('availability [%]')
    ax.title.set_text('downwards availability total [%]')


    '''
    
    #downwards
    MAX = 35*1.05
    plot_availability = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP','P2H','CHP+P2H','P2H+CHP'] )
    plot_availability.loc['RES min','CHP'] = availability_mean['Down'].at[CHP_GAS,'s2_B']+availability_mean['Down'].at[CHP_HRD,'s2_B']+availability_mean['Down'].at[CHP_OIL,'s2_B']
    plot_availability.loc['RES min','P2H'] = availability_mean['Down'].at[P2H_unit,'s3_B']
    plot_availability.loc['RES min','CHP+P2H'] = availability_mean['Down'].at[CHP_GAS,'s4_B']+availability_mean['Down'].at[CHP_HRD,'s4_B']+availability_mean['Down'].at[CHP_OIL,'s4_B']
    plot_availability.loc['RES min','P2H+CHP'] = availability_mean['Down'].at[P2H_unit,'s4_B']
    plot_availability.loc['RES mid','CHP'] = availability_mean['Down'].at[CHP_GAS,'s6_B']+availability_mean['Down'].at[CHP_HRD,'s6_B']+availability_mean['Down'].at[CHP_OIL,'s6_B']
    plot_availability.loc['RES mid','P2H'] = availability_mean['Down'].at[P2H_unit,'s7_B']
    plot_availability.loc['RES mid','CHP+P2H'] = availability_mean['Down'].at[CHP_GAS,'s8_B']+availability_mean['Down'].at[CHP_HRD,'s8_B']+availability_mean['Down'].at[CHP_OIL,'s8_B']
    plot_availability.loc['RES mid','P2H+CHP'] = availability_mean['Down'].at[P2H_unit,'s8_B']
    plot_availability.loc['RES max','CHP'] = availability_mean['Down'].at[CHP_GAS,'s10_B']+availability_mean['Down'].at[CHP_HRD,'s10_B']+availability_mean['Down'].at[CHP_OIL,'s10_B']
    plot_availability.loc['RES max','P2H'] = availability_mean['Down'].at[P2H_unit,'s11_B']
    plot_availability.loc['RES max','CHP+P2H'] = availability_mean['Down'].at[CHP_GAS,'s12_B']+availability_mean['Down'].at[CHP_HRD,'s12_B']+availability_mean['Down'].at[CHP_OIL,'s12_B']
    plot_availability.loc['RES max','P2H+CHP'] = availability_mean['Down'].at[P2H_unit,'s12_B']

    ax = plot_availability[['CHP+P2H','P2H+CHP']].plot.bar(figsize=(18, 6),stacked = True,alpha=0.8, width=0.1, position=-0.5,color=['#b93c46ff','gray'],ylim = (0,MAX), legend = False)
    ax = plot_availability['CHP'].plot.bar(figsize=(18, 6),stacked = True, alpha=0.8, width=0.1, position=1.5,color=['#b93c46ff'],ylim = (0,MAX), legend = True)
    ax = plot_availability['P2H'].plot.bar(figsize=(18, 6),stacked = True,alpha=0.8, width=0.1, position=0.5,color=['gray'],ylim = (0,MAX), legend = True)
    ax.set_ylabel('availability [%]')
    ax.title.set_text('downwards availability mean [%]')
    

    # spinning upwards reserves
    MAX = 35*1.05
    plot_availability = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP','P2H','CHP+P2H','P2H+CHP'] )
    plot_availability.loc['RES min','CHP'] = availability_mean['2U'].at[CHP_GAS,'s2_B']+availability_mean['2U'].at[CHP_HRD,'s2_B']+availability_mean['2U'].at[CHP_OIL,'s2_B']
    plot_availability.loc['RES min','P2H'] = availability_mean['2U'].at[P2H_unit,'s3_B']
    plot_availability.loc['RES min','CHP+P2H'] = availability_mean['2U'].at[CHP_GAS,'s4_B']+availability_mean['2U'].at[CHP_HRD,'s4_B']+availability_mean['2U'].at[CHP_OIL,'s4_B']
    plot_availability.loc['RES min','P2H+CHP'] = availability_mean['2U'].at[P2H_unit,'s4_B']
    plot_availability.loc['RES mid','CHP'] = availability_mean['2U'].at[CHP_GAS,'s6_B']+availability_mean['2U'].at[CHP_HRD,'s6_B']+availability_mean['2U'].at[CHP_OIL,'s6_B']
    plot_availability.loc['RES mid','P2H'] = availability_mean['2U'].at[P2H_unit,'s7_B']
    plot_availability.loc['RES mid','CHP+P2H'] = availability_mean['2U'].at[CHP_GAS,'s8_B']+availability_mean['2U'].at[CHP_HRD,'s8_B']+availability_mean['2U'].at[CHP_OIL,'s8_B']
    plot_availability.loc['RES mid','P2H+CHP'] = availability_mean['2U'].at[P2H_unit,'s8_B']
    plot_availability.loc['RES max','CHP'] = availability_mean['2U'].at[CHP_GAS,'s10_B']+availability_mean['2U'].at[CHP_HRD,'s10_B']+availability_mean['2U'].at[CHP_OIL,'s10_B']
    plot_availability.loc['RES max','P2H'] = availability_mean['2U'].at[P2H_unit,'s11_B']
    plot_availability.loc['RES max','CHP+P2H'] = availability_mean['2U'].at[CHP_GAS,'s12_B']+availability_mean['2U'].at[CHP_HRD,'s12_B']+availability_mean['2U'].at[CHP_OIL,'s12_B']
    plot_availability.loc['RES max','P2H+CHP'] = availability_mean['2U'].at[P2H_unit,'s12_B']

    ax = plot_availability[['CHP+P2H','P2H+CHP']].plot.bar(figsize=(18, 6),stacked = True,alpha=0.8, width=0.1, position=-0.5,color=['#b93c46ff','gray'],ylim = (0,MAX), legend = False)
    ax = plot_availability['CHP'].plot.bar(figsize=(18, 6),stacked = True, alpha=0.8, width=0.1, position=1.5,color=['#b93c46ff'],ylim = (0,MAX), legend = True)
    ax = plot_availability['P2H'].plot.bar(figsize=(18, 6),stacked = True,alpha=0.8, width=0.1, position=0.5,color=['gray'],ylim = (0,MAX), legend = True)
    ax.set_ylabel('availability [%]')
    ax.title.set_text('2U  availability [%]')
#%% technical potential
    ax = tech_pot.plot(figsize=(18, 6),alpha=0.8,color=['#b93c46ff','gray','b'],xlim = (0.05,0.116),ylim = (0,3000*1.05), legend = True)
    ax.set_ylabel('lost load [h]')
    ax.title.set_text('Technical potential [%]')
    
    
    
#%% shadow prices
    
    
    
    
    
    
    
    

#%% profit
    '''
    Plot_Profit = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP_GAS','CHP_HRD','CHP_OIL','P2H'] )
    Plot_Profit.loc['RES min','CHP'] = 0
    Plot_Profit.loc['RES min','P2H'] = availability['s3_B']['2U'].at['[21] - IT_P2HT_OTH','sum']
    Plot_Profit.loc['RES min','CHP+P2H'] = 0

    Plot_Profit.loc['RES mid','CHP'] = 0
    Plot_Profit.loc['RES mid','P2H'] = availability['s7_B']['2U'].at['[21] - IT_P2HT_OTH','sum']
    Plot_Profit.loc['RES mid','CHP+P2H'] = 0

    Plot_Profit.loc['RES max','CHP'] = 0
    Plot_Profit.loc['RES max','P2H'] = availability['s11_B']['2U'].at['[21] - IT_P2HT_OTH','sum']
    Plot_Profit.loc['RES max','CHP+P2H'] = 0
    
    df1 = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP', 'P2H','CHP+P2H'] )

    fig, axes = plt.subplots(nrows=1, ncols=2,figsize=(18, 6))
    plot_availability.plot(ax=axes[0],kind='bar')
    df1.plot(ax=axes[1], kind='bar');
    '''
#%% analyzing graph

def analyzinggraph(scenarios_dic,index,scenarios):
    z='IT'
    unit = '[21] - IT_P2HT_OTH'


    #dictionaries
    GenerationOutput = pd.DataFrame()
    demand = pd.DataFrame()
    shedload = pd.DataFrame()
    powerconsumption = pd.DataFrame()
    storagelevel = pd.DataFrame()
    curtailment = pd.DataFrame()

    for SCEN in scenarios:
        INPUTS = 'inputs_##'
        INPUTS = INPUTS.replace('##', SCEN)
        RESULTS = 'results_##'
        RESULTS = RESULTS.replace('##', SCEN)
        inputs = scenarios_dic[INPUTS]
        results = scenarios_dic[RESULTS]
        
        # GenerationOutput from powerplants and storage (what storage?)
        Generation=ds.get_plot_data(inputs, results, z)
        GenerationOutput[SCEN] = Generation.sum(axis=1)
        
        # demand
        if 'OutputPowerConsumption' in results:
            demand_p2h = ds.filter_by_zone(results['OutputPowerConsumption'], inputs, z)
            demand_p2h = demand_p2h.sum(axis=1)
        else:
            demand_p2h = pd.Series(0, index=results['OutputPower'].index)
        if ('Flex', z) in inputs['param_df']['Demand']:
            demand_flex = inputs['param_df']['Demand'][('Flex', z)]
        else:
            demand_flex = pd.Series(0, index=results['OutputPower'].index)
        if 'OutputDemandModulation' in results and z in results['OutputDemandModulation']:
            shifted_load = -results['OutputDemandModulation'][z]
            shifted_load = pd.Series(shifted_load, index=results['OutputPower'].index).fillna(0)
        else:
            shifted_load = pd.Series(0, index=results['OutputPower'].index)
    
        
        demand_da = inputs['param_df']['Demand'][('DA', z)]
        demand[SCEN] = demand_da + demand_p2h + demand_flex+shifted_load
        
        #shed load
        if z in results['OutputShedLoad']:
            shedload[SCEN] = results['OutputShedLoad'][z]
        else:
            shedload[SCEN] = 0
        
        #curtailment
        if z in results['OutputCurtailedPower']:
            curtailment[SCEN] = results['OutputCurtailedPower'][z]
        else:
            curtailment[SCEN] = 0

        #powerconsumption
        if unit in results['OutputPowerConsumption']:
            powerconsumption[SCEN] = results['OutputPowerConsumption'][unit]
        else:
            powerconsumption[SCEN] = 0
            
        #storage level
        if unit in results['OutputStorageLevel']:
            storagelevel[SCEN] = results['OutputStorageLevel'][unit]
        else:
            storagelevel[SCEN] = 0
            
    pd.DataFrame.fillna(curtailment,0,inplace=True)
    pd.DataFrame.fillna(shedload,0,inplace=True)
    pd.DataFrame.fillna(powerconsumption,0,inplace=True)
    pd.DataFrame.fillna(storagelevel,0,inplace=True)

    fig, axes = plt.subplots(nrows=5, ncols=2, sharex=True,sharey=False, figsize=(26, 24), frameon=True)
    fig.tight_layout()
    
    axes[0,0].plot(index, GenerationOutput[scenarios[0]][index],color='k',label='Power output s9')
    axes[0,0].plot(index, GenerationOutput[scenarios[1]][index],color='b',label='Power output s11_A')
    axes[0,0].plot(index, GenerationOutput[scenarios[2]][index],color='m',label='Power output s11_B')
    
    axes[1,0].plot(index, demand[scenarios[0]][index],color='k',label='demand s9')
    axes[1,0].plot(index, demand[scenarios[1]][index],color='b',label='demand s11_A')
    axes[1,0].plot(index, demand[scenarios[2]][index],color='m',label='demand s11_B')
        
    axes[0,1].plot(index, shedload[scenarios[0]][index],color='k',label='shed load s9')
    axes[0,1].plot(index, shedload[scenarios[1]][index],color='b',label='shed load s11_A')
    axes[0,1].plot(index, shedload[scenarios[2]][index],color='m',label='shed load s11_B')

    axes[1,1].plot(index, curtailment[scenarios[0]][index],color='k',label='curtailment s9')
    axes[1,1].plot(index, curtailment[scenarios[1]][index],color='b',label='curtailment s11_A')
    axes[1,1].plot(index, curtailment[scenarios[2]][index],color='m',label='curtailment s11_B')
    
    axes[3,1].plot(index, storagelevel[scenarios[0]][index],color='k',label='thermal storage s9')
    axes[3,1].plot(index, storagelevel[scenarios[1]][index],color='b',label='thermal storage s11_A')
    axes[3,1].plot(index, storagelevel[scenarios[2]][index],color='m',label='thermal storage s11_B')
    
    axes[2,1].plot(index, powerconsumption[scenarios[0]][index],color='k',label='powerconsumption s9')
    axes[2,1].plot(index, powerconsumption[scenarios[1]][index],color='b',label='powerconsumption s11_A')
    axes[2,1].plot(index, powerconsumption[scenarios[2]][index],color='m',label='powerconsumption s11_B')

    axes[2,0].plot(index, demand[scenarios[0]][index],color='k',label='demand s9')
    axes[2,0].plot(index, GenerationOutput[scenarios[0]][index],color='b',label='Power output s9')
    
    axes[3,0].plot(index, demand[scenarios[1]][index],color='k',label='demand s11_A')
    axes[3,0].plot(index, GenerationOutput[scenarios[1]][index],color='b',label='Power output s11_A')
    
    axes[4,0].plot(index, demand[scenarios[2]][index],color='k',label='demand s11_B')
    axes[4,0].plot(index, GenerationOutput[scenarios[2]][index],color='b',label='Power output s11_B')


    axes[0,0].legend()
    axes[1,0].legend()
    axes[2,0].legend()
    axes[3,0].legend()
    axes[4,0].legend()
    axes[0,1].legend()
    axes[1,1].legend()
    axes[2,1].legend()
    axes[3,1].legend()
    axes[4,1].legend()
    
    return True

#%% loading graphs
#analyzinggraph(scenarios_dic,april1,['s1','s3_A','s3_B'])
#analyzinggraph(scenarios_dic,april2,['s5','s7_A','s7_B'])
#analyzinggraph(scenarios_dic,july1,['s9','s11_A','s11_B'])

#ds.plot_zone(scenarios_dic['inputs_s1'],scenarios_dic['results_s1'],'IT',rng=october1)
#ds.plot_zone(scenarios_dic['inputs_s3_A'],scenarios_dic['results_s3_A'],'IT',rng=october1)

#%% 2050 scenario
    


    #%% availability
    
    
    
    
    
    #%% market prices
    
    
    
    
    
    #%% profit
    
    
    
    
    #%% total system cost
    system_costs_2050 = pd.DataFrame(0, index = ['none','DHS','RES','DHS+RES'], columns =['participation'])
    system_costs_2050.loc['none','participation'] = calculations.at['s13','postsystemcostIT']
    system_costs_2050.loc['DHS','participation'] = calculations.at['s14','postsystemcostIT']
    system_costs_2050.loc['RES','participation'] = calculations.at['s15','postsystemcostIT']
    system_costs_2050.loc['DHS+RES','participation'] = calculations.at['s16','postsystemcostIT']
    
    MAX = system_costs_2050.max().max()*1.05
    
    ax = system_costs_2050.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('total Italian system costs [euro]')
    ax.title.set_text('Italian system cost for the 2050 scenarios')
    
    
    costsdifference = pd.DataFrame(0, index = ['none','DHS','RES','DHS+RES'], columns =['participation'] )
    costsdifference.loc['DHS','participation'] = (calculations.at['s14','postsystemcostIT']-calculations.at['s13','postsystemcostIT'])/calculations.at['s13','postsystemcostIT']*100
    costsdifference.loc['RES','participation'] = (calculations.at['s15','postsystemcostIT']-calculations.at['s13','postsystemcostIT'])/calculations.at['s13','postsystemcostIT']*100
    costsdifference.loc['DHS+RES','participation'] = (calculations.at['s16','postsystemcostIT']-calculations.at['s13','postsystemcostIT'])/calculations.at['s13','postsystemcostIT']*100
    
    MAX = costsdifference.max().max()*1.05
    MIN = costsdifference.min().min()*1.05
    
    ax = costsdifference.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (MIN,MAX), legend = False)
    ax.set_ylabel('[%]')
    ax.legend(['CHP_participation','P2H_participation','CHP+P2H_participation'])
    ax.title.set_text('cost savings through participation (negative is a cost reduction)')


    #%% curtailment
    curtailment_2050 = pd.DataFrame(0, index = ['none','DHS','RES','DHS+RES'], columns =['participation'])
    curtailment_2050.loc['none','participation'] = curtailment.at['s13','sum']
    curtailment_2050.loc['DHS','participation'] = curtailment.at['s14','sum']
    curtailment_2050.loc['RES','participation'] = curtailment.at['s15','sum']
    curtailment_2050.loc['DHS+RES','participation'] = curtailment.at['s16','sum']
    
    MAX = curtailment_2050.max().max()*1.05
    
    ax = curtailment_2050.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.title.set_text('curtailment for the 2050 scenarios')

    #%% total emissions
    
    plot_emissions = pd.DataFrame(0, index = ['none','DHS','RES','DHS+RES'], columns = [])
    plot_emissions.loc['none','none power system'] = calculations.at['s13','co2_power']
    plot_emissions.loc['none','none heat slack'] = calculations.at['s13','co2_heatslack']
    plot_emissions.loc['DHS','DHS power system'] = calculations.at['s14','co2_power']
    plot_emissions.loc['DHS','DHS heat slack'] = calculations.at['s14','co2_heatslack']
    plot_emissions.loc['RES','RES power system'] = calculations.at['s15','co2_power']
    plot_emissions.loc['RES','RES heat slack'] = calculations.at['s15','co2_heatslack']
    plot_emissions.loc['DHS+RES','DHS+RES power system'] = calculations.at['s16','co2_power']
    plot_emissions.loc['DHS+RES','DHS+RES heat slack'] = calculations.at['s16','co2_heatslack']


    MAX = 80000000#(plot_emissions.at['RES min','CHP power system']+plot_emissions.at['RES min','CHP heat slack'])*1.05
    fig, ax = plt.subplots()
    ax = plot_emissions[['none power system','none heat slack']].plot.bar(figsize=(18, 6),stacked = True, alpha=0.6, width=0.8,colormap="rainbow",ax=ax,ylim = (0,MAX), legend = True)
    ax = plot_emissions[['DHS power system','DHS heat slack']].plot.bar(figsize=(18, 6),stacked = True,alpha=0.6, width=0.8,colormap="rainbow",ax=ax,ylim = (0,MAX), legend = True)
    ax = plot_emissions[['RES power system','RES heat slack']].plot.bar(figsize=(18, 6),stacked = True,alpha=0.6, width=0.8,colormap="rainbow",ax=ax,ylim = (0,MAX), legend = True)
    ax = plot_emissions[['DHS+RES power system','DHS+RES heat slack']].plot.bar(figsize=(18, 6),stacked = True,alpha=0.6, width=0.8,colormap="rainbow",ax=ax,ylim = (0,MAX), legend = True)
    ax.set_ylabel('emissions [ton co2]')
    ax.title.set_text('emissions produced in Italy')

    
    emission_difference = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns =['none','DHS','RES','DHS+RES'] )
    emission_difference.loc['RES min','CHP_participation'] = (calculations.at['s2_B','co2_total']-calculations.at['s2_A','co2_total'])/calculations.at['s2_A','co2_total']*100
    emission_difference.loc['RES min','P2H_participation'] = (calculations.at['s3_B','co2_total']-calculations.at['s3_A','co2_total'])/calculations.at['s3_A','co2_total']*100
    emission_difference.loc['RES min','CHP+P2H_participation'] = (calculations.at['s4_B','co2_total']-calculations.at['s4_A','co2_total'])/calculations.at['s4_A','co2_total']*100
    emission_difference.loc['RES mid','CHP_participation'] = (calculations.at['s6_B','co2_total']-calculations.at['s6_A','co2_total'])/calculations.at['s6_A','co2_total']*100
    emission_difference.loc['RES mid','P2H_participation'] = (calculations.at['s7_B','co2_total']-calculations.at['s7_A','co2_total'])/calculations.at['s7_A','co2_total']*100
    emission_difference.loc['RES mid','CHP+P2H_participation'] = (calculations.at['s8_B','co2_total']-calculations.at['s8_A','co2_total'])/calculations.at['s8_A','co2_total']*100
    emission_difference.loc['RES max','CHP_participation'] = (calculations.at['s10_B','co2_total']-calculations.at['s10_A','co2_total'])/calculations.at['s10_A','co2_total']*100
    emission_difference.loc['RES max','P2H_participation'] = (calculations.at['s11_B','co2_total']-calculations.at['s11_A','co2_total'])/calculations.at['s11_A','co2_total']*100
    emission_difference.loc['RES max','CHP+P2H_participation'] = (calculations.at['s12_B','co2_total']-calculations.at['s12_A','co2_total'])/calculations.at['s12_A','co2_total']*100

    MAX = emission_difference.max().max()*1.05
    MIN = emission_difference.min().min()*1.05
    
    ax = emission_difference.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (MIN,MAX), legend = False)
    ax.set_ylabel('[%]')
    ax.legend(['CHP_participation','P2H_participation','CHP+P2H_participation'])
    ax.title.set_text('emissions saving through participation (negative is a decreasement)')

    #%% 2U availability
    #CHP unit
    testgas=(scenarios_dic['results_s2_B']['OutputPower'][CHP_GAS]/511)
    for x in [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]:
        for i in testgas.index:
            if testgas.loc[i]>=x and testgas.loc[i]<(x+1):
                testgas.loc[i]=testgas.at[i]-x        
    testgas.mean()
        

    testhrd=(scenarios_dic['results_s4_B']['OutputPower'][CHP_HRD]/1078)
    for x in [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]:
        for i in testhrd.index:
            if testhrd.loc[i]>=x and testhrd.loc[i]<(x+1):
                testhrd.loc[i]=testhrd.at[i]-x
    testhrd.mean()

    testoil=(scenarios_dic['results_s4_B']['OutputPower'][CHP_OIL]/476)
    for x in [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]:
        for i in testoil.index:
            if testoil.loc[i]>=x and testoil.loc[i]<(x+1):
                testoil.loc[i]=testoil.at[i]-x
    testoil.mean()


    
