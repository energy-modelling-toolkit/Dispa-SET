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

LoadData = False
CalculateVariables = False
cost = True
plot = True
# %%Load the inputs and the results of the simulations
scenarios = ['s1','s2_A','s2_B','s3_A','s3_B','s5','s6_A','s6_B','s7_A','s7_B','s9','s10_A','s10_B','s11_A','s11_B']

if LoadData == True:
    scenarios_dic = {}
    for SCEN in scenarios:
        PATH = '/Users/tomfranssens/Dispa-SET.git/Simulations/simulation_##'
        PATH = PATH.replace('##', SCEN)
        INPUTS = 'inputs_##'
        INPUTS = INPUTS.replace('##', SCEN)
        RESULTS = 'results_##'
        RESULTS = RESULTS.replace('##', SCEN)
        scenarios_dic[INPUTS],scenarios_dic[RESULTS] = ds.get_sim_results(path=PATH,cache=False)
        

#%% index, unit and zone
index = scenarios_dic['inputs_s1']['config']['idx']
P2H_unit = '[21] - IT_P2HT_OTH'
CHP_GAS = '[10] - IT_COMC_GAS_CHP'
CHP_HRD = '[13] - IT_STUR_HRD_CHP'
CHP_OIL = '[14] - IT_STUR_OIL_CHP'
unit_list = [P2H_unit,CHP_GAS,CHP_HRD,CHP_OIL]
zone = 'IT'

#%% calculations
if CalculateVariables == True:

#costs, emissions
    calculations = pd.DataFrame(0,index = scenarios, columns = ['postsystemcostIT','gams_costs',
                                                                'co2_power','co2_heatslack','co2_total'])
#reserve demand
    reservedemand = pd.DataFrame(0,index= scenarios,columns=['upwards','downwards'])
#hourly availability
    hourly_availability = {}
# mean availability
    availability_mean = {}
    availability_mean['2U'] = pd.DataFrame(0,index= unit_list,columns=scenarios)
    availability_mean['3U'] = pd.DataFrame(0,index= unit_list,columns=scenarios)
    availability_mean['Up'] = pd.DataFrame(0,index= unit_list,columns=scenarios)
    availability_mean['Down'] = pd.DataFrame(0,index= unit_list,columns=scenarios)
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


#    shadowprices = {}
#    cashflow = {}
#    cashflows_sum = {}
#    profit = {}


    
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

# lost load
        for lostload in [results['LostLoad_2D'],results['LostLoad_2U'],results['LostLoad_3U'],
                         results['LostLoad_MaxPower'],results['LostLoad_MinPower'],
                         results['LostLoad_RampDown'],results['LostLoad_RampUp']]:
            if zone in lostload and lostload[zone].sum() != 0:
                print(lostload)
        
# reserve demand and availability (unit specific) include the total !!!!
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

        hourly_availability[SCEN] = TMP1
        
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
        #shadowprices[SCEN] = ds.shadowprices(results, zone)

# profit one unit (cost negative and graf with profit not cashflow)
        #if unit in inputs['units'].index:
        #    cashflow[SCEN], cashflows_sum[SCEN], profit[SCEN] = ds.Cashflows(inputs, results, unit)
    
        
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
    plot_reservedemand = pd.DataFrame(0,index=['RES min','RES mid','RES max','ProRes1'],columns=['Reserve Demand'])
    plot_reservedemand.loc['RES min','Reserve Demand'] = reservedemand.at['s1','upwards']
    plot_reservedemand.loc['RES mid','Reserve Demand'] = reservedemand.at['s5','upwards']
    plot_reservedemand.loc['RES max','Reserve Demand'] = reservedemand.at['s9','upwards']
    plot_reservedemand.loc['ProRes1','Reserve Demand'] = 0
    ax = plot_reservedemand.plot.bar(figsize=(18, 6), alpha=0.2, legend = False)
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
    totalcosts.loc['RES min','CHP+P2H'] = 0
    totalcosts.loc['RES mid','CHP+P2H'] = 0
    totalcosts.loc['RES max','CHP+P2H'] = 0
    ax = totalcosts.plot.bar(figsize=(18, 6), alpha=0.2)
    ax.set_ylabel('total Italian system costs [euro]')
    ax.title.set_text('Italian system cost including heat slack')
    
    #%% total emissions between scenarios with mid RES for now
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
    
    
    #%% curtailment

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
    CHP_P2H_curtailment.loc['RES min','no coupling'] = 0#curtailment.at['s1','sum']
    CHP_P2H_curtailment.loc['RES min','no participation'] = 0#curtailment.at['s4_A','sum']
    CHP_P2H_curtailment.loc['RES min','participation'] = 0#curtailment.at['s4_B','sum']
    CHP_P2H_curtailment.loc['RES mid','no coupling'] = 0#curtailment.at['s5','sum']
    CHP_P2H_curtailment.loc['RES mid','no participation'] = 0#curtailment.at['s8_A','sum']
    CHP_P2H_curtailment.loc['RES mid','participation'] = 0#curtailment.at['s8_B','sum']
    CHP_P2H_curtailment.loc['RES max','no coupling'] = 0#curtailment.at['s9','sum']
    CHP_P2H_curtailment.loc['RES max','no participation'] = 0#curtailment.at['s12_A','sum']
    CHP_P2H_curtailment.loc['RES max','participation'] = 0#curtailment.at['s12_B','sum']
    
    MAX = max([CHP_curtailment.max().max(),P2H_curtailment.max().max(),CHP_P2H_curtailment.max().max()])*1.05
    
    ax = CHP_curtailment.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    
    ax = P2H_curtailment.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    
    ax = CHP_P2H_curtailment.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('curtailment [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    
    
    
    
    #%% load shedding
    
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
    CHP_P2H_loadshedding.loc['RES min','no participation'] = 0#load_shedding.at['s4_A','sum']
    CHP_P2H_loadshedding.loc['RES min','participation'] = 0#load_shedding.at['s4_B','sum']
    CHP_P2H_loadshedding.loc['RES mid','no coupling'] = load_shedding.at['s5','sum']
    CHP_P2H_loadshedding.loc['RES mid','no participation'] = 0#load_shedding.at['s8_A','sum']
    CHP_P2H_loadshedding.loc['RES mid','participation'] = 0#load_shedding.at['s8_B','sum']
    CHP_P2H_loadshedding.loc['RES max','no coupling'] = load_shedding.at['s9','sum']
    CHP_P2H_loadshedding.loc['RES max','no participation'] = 0#load_shedding.at['s12_A','sum']
    CHP_P2H_loadshedding.loc['RES max','participation'] = 0#load_shedding.at['s12_B','sum']
    
    MAX = max([CHP_loadshedding.max().max(),P2H_loadshedding.max().max(),CHP_P2H_loadshedding.max().max()])*1.05
    
    ax = CHP_loadshedding.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('loadshedding [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    
    ax = P2H_loadshedding.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('loadshedding [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    
    ax = CHP_P2H_loadshedding.plot.bar(figsize=(18, 6), alpha=0.2,ylim = (0,MAX), legend = False)
    ax.set_ylabel('loadshedding [MWh]')
    ax.legend(['no coupling','no participation','participation'])
    
#%% availability
    # upwards
    MAX = 300*1.05
    plot_availability = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP_GAS','CHP_HRD','CHP_OIL','P2H','CHP+P2H'] )
    plot_availability.loc['RES min','CHP_GAS'] = availability_total['Up'].at[CHP_GAS,'s2_B']
    plot_availability.loc['RES min','CHP_HRD'] = availability_total['Up'].at[CHP_HRD,'s2_B']
    plot_availability.loc['RES min','CHP_OIL'] = availability_total['Up'].at[CHP_OIL,'s2_B']
    plot_availability.loc['RES min','P2H'] = availability_total['Up'].at[P2H_unit,'s3_B']
    plot_availability.loc['RES min','CHP+P2H'] = 0
    plot_availability.loc['RES mid','CHP_GAS'] = availability_total['Up'].at[CHP_GAS,'s6_B']
    plot_availability.loc['RES mid','CHP_HRD'] = availability_total['Up'].at[CHP_HRD,'s6_B']
    plot_availability.loc['RES mid','CHP_OIL'] = availability_total['Up'].at[CHP_OIL,'s6_B']
    plot_availability.loc['RES mid','P2H'] = availability_total['Up'].at[P2H_unit,'s7_B']
    plot_availability.loc['RES mid','CHP+P2H'] = 0
    plot_availability.loc['RES max','CHP_GAS'] = availability_total['Up'].at[CHP_GAS,'s10_B']
    plot_availability.loc['RES max','CHP_HRD'] = availability_total['Up'].at[CHP_HRD,'s10_B']
    plot_availability.loc['RES max','CHP_OIL'] = availability_total['Up'].at[CHP_OIL,'s10_B']
    plot_availability.loc['RES max','P2H'] = availability_total['Up'].at[P2H_unit,'s11_B']
    plot_availability.loc['RES max','CHP+P2H'] = 0

    ax = plot_availability[['CHP_GAS','CHP_HRD','CHP_OIL']].plot.bar(figsize=(18, 6),stacked = True, alpha=0.2, width=0.1, position=1.5,colormap="rainbow",ylim = (0,MAX), legend = True)
    ax = plot_availability['P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=0.5,ylim = (0,MAX), legend = True)
    ax = plot_availability['CHP+P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=-0.5,ylim = (0,MAX), legend = False)
    ax.set_ylabel('availability [%]')
        
    #downwards
    MAX = 300*1.05
    plot_availability = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP_GAS','CHP_HRD','CHP_OIL','P2H','CHP+P2H'] )
    plot_availability.loc['RES min','CHP_GAS'] = availability_total['Down'].at[CHP_GAS,'s2_B']
    plot_availability.loc['RES min','CHP_HRD'] = availability_total['Down'].at[CHP_HRD,'s2_B']
    plot_availability.loc['RES min','CHP_OIL'] = availability_total['Down'].at[CHP_OIL,'s2_B']
    plot_availability.loc['RES min','P2H'] = availability_total['Down'].at[P2H_unit,'s3_B']
    plot_availability.loc['RES min','CHP+P2H'] = 0
    plot_availability.loc['RES mid','CHP_GAS'] = availability_total['Down'].at[CHP_GAS,'s6_B']
    plot_availability.loc['RES mid','CHP_HRD'] = availability_total['Down'].at[CHP_HRD,'s6_B']
    plot_availability.loc['RES mid','CHP_OIL'] = availability_total['Down'].at[CHP_OIL,'s6_B']
    plot_availability.loc['RES mid','P2H'] = availability_total['Down'].at[P2H_unit,'s7_B']
    plot_availability.loc['RES mid','CHP+P2H'] = 0
    plot_availability.loc['RES max','CHP_GAS'] = availability_total['Down'].at[CHP_GAS,'s10_B']
    plot_availability.loc['RES max','CHP_HRD'] = availability_total['Down'].at[CHP_HRD,'s10_B']
    plot_availability.loc['RES max','CHP_OIL'] = availability_total['Down'].at[CHP_OIL,'s10_B']
    plot_availability.loc['RES max','P2H'] = availability_total['Down'].at[P2H_unit,'s11_B']
    plot_availability.loc['RES max','CHP+P2H'] = 0

    ax = plot_availability[['CHP_GAS','CHP_HRD','CHP_OIL']].plot.bar(figsize=(18, 6),stacked = True, alpha=0.2, width=0.1, position=1.5,colormap="rainbow",ylim = (0,MAX), legend = True)
    ax = plot_availability['P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=0.5,ylim = (0,MAX), legend = True)
    ax = plot_availability['CHP+P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=-0.5,ylim = (0,MAX), legend = False)
    ax.set_ylabel('availability [%]')
    
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
    MAX = 300*1.05
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
    plot_availability.loc['RES max','CHP_GAS'] = availability_total['3U'].at[CHP_GAS,'s2_B']
    plot_availability.loc['RES max','CHP_HRD'] = availability_total['3U'].at[CHP_HRD,'s2_B']
    plot_availability.loc['RES max','CHP_OIL'] = availability_total['3U'].at[CHP_OIL,'s2_B']
    plot_availability.loc['RES max','P2H'] = availability_total['3U'].at[P2H_unit,'s11_B']
    plot_availability.loc['RES max','CHP+P2H'] = 0

    ax = plot_availability[['CHP_GAS','CHP_HRD','CHP_OIL']].plot.bar(figsize=(18, 6),stacked = True, alpha=0.2, width=0.1, position=1.5,colormap="rainbow",ylim = (0,MAX), legend = True)
    ax = plot_availability['P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=0.5,ylim = (0,MAX), legend = True)
    ax = plot_availability['CHP+P2H'].plot.bar(figsize=(18, 6),alpha=0.2, width=0.1, position=-0.5,ylim = (0,MAX), legend = False)
    ax.set_ylabel('availability [%]')

    
#%% shadow prices
    
    
    
    
    
    
    
    

#%% profit wait for more results 
    
    '''
    plot_availability = pd.DataFrame(0, index = ['RES min', 'RES mid', 'RES max'], columns = ['CHP', 'P2H','CHP+P2H'] )
    plot_availability.loc['RES min','CHP'] = 0
    plot_availability.loc['RES min','P2H'] = availability['s3_B']['2U'].at['[21] - IT_P2HT_OTH','sum']
    plot_availability.loc['RES min','CHP+P2H'] = 0
    plot_availability.loc['RES mid','CHP'] = 0
    plot_availability.loc['RES mid','P2H'] = availability['s7_B']['2U'].at['[21] - IT_P2HT_OTH','sum']
    plot_availability.loc['RES mid','CHP+P2H'] = 0
    plot_availability.loc['RES max','CHP'] = 0
    plot_availability.loc['RES max','P2H'] = availability['s11_B']['2U'].at['[21] - IT_P2HT_OTH','sum']
    plot_availability.loc['RES max','CHP+P2H'] = 0
    
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

