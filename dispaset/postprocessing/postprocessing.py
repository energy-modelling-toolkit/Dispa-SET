# -*- coding: utf-8 -*-
"""
Set of functions useful to analyse to DispaSET output data.

@author: Sylvain Quoilin, JRC
"""

from __future__ import division

import datetime as dt
import logging
import sys

import numpy as np
import pandas as pd

from ..common import commons
from .data_handler import ds_to_df




def get_load_data(inputs, z):
    """ 
    Get the load curve, the residual load curve, and the net residual load curve of a specific zone

    :param inputs:  DispaSET inputs (output of the get_sim_results function)
    :param z:       Zone to consider (e.g. 'BE')
    :return out:    Dataframe with the following columns:
                        Load:               Load curve of the specified zone
                        ResidualLoad:       Load minus the production of variable renewable sources
                        NetResidualLoad:    Residual netted from the interconnections with neightbouring zones
    """
    datain = inputs['param_df']
    out = pd.DataFrame(index=datain['Demand'].index)
    out['Load'] = datain['Demand']['DA', z]
    if ('Flex', z) in datain['Demand']:
        out['Load'] += datain['Demand'][('Flex', z)] 
    # Listing power plants with non-dispatchable power generation:
    VREunits = []
    VRE = np.zeros(len(out))
    for t in commons['tech_renewables']:
        for u in datain['Technology']:
            if datain['Technology'].loc[t, u]:
                VREunits.append(u)
                VRE = VRE + datain['AvailabilityFactor'][u].values * datain['PowerCapacity'].loc[u, 'PowerCapacity']
    Interconnections = np.zeros(len(out))
    for l in datain['FlowMinimum']:
        if l[:2] == z:
            Interconnections = Interconnections - datain['FlowMinimum'][l].values
        elif l[-2:] == z:
            Interconnections = Interconnections + datain['FlowMinimum'][l].values
    out['ResidualLoad'] = out['Load'] - VRE
    out['NetResidualLoad'] = out['ResidualLoad'] - Interconnections
    return out


def aggregate_by_fuel(PowerOutput, Inputs, SpecifyFuels=None):
    """
    This function sorts the power generation curves of the different units by technology

    :param PowerOutput:     Dataframe of power generationwith units as columns and time as index
    :param Inputs:          Dispaset inputs version 2.1.1
    :param SpecifyFuels:     If not all fuels should be considered, list containing the relevant ones
    :returns PowerByFuel:    Dataframe with power generation by fuel
    """
    if SpecifyFuels is None:
        if isinstance(Inputs, list):
            fuels = Inputs[0]['f']
        elif isinstance(Inputs, dict):
            fuels = Inputs['sets']['f']
        else:
            logging.error('Inputs variable no valid')
            sys.exit(1)
    else:
        fuels = SpecifyFuels
    PowerByFuel = pd.DataFrame(0, index=PowerOutput.index, columns=fuels)
    uFuel = Inputs['units']['Fuel']

    for u in PowerOutput:
        if uFuel[u] in fuels:
            PowerByFuel[uFuel[u]] = PowerByFuel[uFuel[u]] + PowerOutput[u]
        else:
            logging.warning('Fuel not found for unit ' + u + ' with fuel ' + uFuel[u])

    return PowerByFuel


def filter_by_zone(Output, inputs, z):
    """
    This function filters the dispaset Output dataframes by zone

    :param Output:          Dataframe of dispaset output with units as columns and time as index
    :param Inputs:          Dispaset inputs 
    :param z:               Selected zone (e.g. 'BE')
    :returns Power:         Dataframe with power generation by zone
    """
    loc = inputs['units']['Zone']
    ZoneOutput = Output.loc[:, [u for u in Output.columns if loc[u] == z]]
    return ZoneOutput


def filter_by_tech(PowerOutput, inputs, t):
    """
    This function filters the dispaset power output dataframe by technology

    :param PowerOutput:   Dataframe of power generation with units as columns and time as index
    :param inputs:        Dispaset inputs version 2.1.1
    :param t:             Selected tech (e.g. 'HDAM')
    :returns Power:
    """
    loc = inputs['units']['Technology']
    Power = PowerOutput.loc[:, [u for u in PowerOutput.columns if loc[u] == t]]
    return Power


def filter_by_storage(PowerOutput, Inputs, StorageSubset=None):
    """
    This function filters the power generation curves of the different storage units by storage type

    :param PowerOutput:     Dataframe of power generationwith units as columns and time as index
    :param Inputs:          Dispaset inputs version 2.1.1
    :param SpecifySubset:   If not all EES storages should be considered, list containing the relevant ones
    :returns PowerByFuel:   Dataframe with power generation by fuel
    """
    storages = Inputs['sets'][StorageSubset]
    Power = PowerOutput.loc[:, PowerOutput.columns.isin(storages)]
    return Power


def get_plot_data(inputs, results, z):
    """
    Function that reads the results dataframe of a DispaSET simulation and extract the dispatch data spedific to one zone

    :param results:         Pandas dataframe with the results (output of the GdxToDataframe function)
    :param z:               Zone to be considered (e.g. 'BE')
    :returns plotdata:      Dataframe with the dispatch data storage and outflows are negative
    """
    tmp = filter_by_zone(results['OutputPower'], inputs, z)
    plotdata = aggregate_by_fuel(tmp, inputs)

    if 'OutputStorageInput' in results:
        #onnly take the columns that correspond to storage units (StorageInput is also used for CHP plants):
        cols = [col for col in results['OutputStorageInput'] if inputs['units'].loc[col,'Technology'] in commons['tech_storage']]
        tmp = filter_by_zone(results['OutputStorageInput'][cols], inputs, z)
        bb = pd.DataFrame()
        for tech in commons['tech_storage']:
            aa = filter_by_tech(tmp, inputs, tech)
            aa = aa.sum(axis=1)
            aa = pd.DataFrame(aa,columns=[tech])
            bb = pd.concat([bb,aa],axis=1)
        bb = -bb
        plotdata = pd.concat([plotdata,bb], axis=1)
        # plotdata['Storage'] = -tmp.sum(axis=1)
    else:
        plotdata['Storage'] = 0
    plotdata.fillna(value=0, inplace=True)

    plotdata['FlowIn'] = 0
    plotdata['FlowOut'] = 0
    for col in results['OutputFlow']:
        from_node, to_node = col.split('->')
        if to_node.strip() == z:
            plotdata['FlowIn'] = plotdata['FlowIn'] + results['OutputFlow'][col]
        if from_node.strip() == z:
            plotdata['FlowOut'] = plotdata['FlowOut'] - results['OutputFlow'][col]

    # re-ordering columns:
    OrderedColumns = [col for col in commons['MeritOrder'] if col in plotdata.columns]
    plotdata = plotdata[OrderedColumns]
    
    # remove empty columns:
    for col in plotdata.columns:
        if plotdata[col].max() == 0 and plotdata[col].min()==0:
            del plotdata[col]

    return plotdata


def get_imports(flows, z):
    """ 
    Function that computes the balance of the imports/exports of a given zone

    :param flows:           Pandas dataframe with the timeseries of the exchanges
    :param z:               Zone to consider
    :returns NetImports:    Scalar with the net balance over the whole time period
    """
    NetImports = 0
    for key in flows:
        if key[:len(z)] == z:
            NetImports -= flows[key].sum()
        elif key[-len(z):] == z:
            NetImports += flows[key].sum()
    return NetImports


# %%
def get_result_analysis(inputs, results):
    """
    Reads the DispaSET results and provides useful general information to stdout

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    """

    # inputs into the dataframe format:
    dfin = inputs['param_df']

    StartDate = inputs['config']['StartDate']
    StopDate = inputs['config']['StopDate']
    index = pd.date_range(start=dt.datetime(*StartDate), end=dt.datetime(*StopDate), freq='h')

    # Aggregated values:
    demand = {}
    for z in inputs['sets']['n']:
        if 'OutputPowerConsumption' in results:
            demand_p2h = filter_by_zone(results['OutputPowerConsumption'], inputs, z)
            demand_p2h = demand_p2h.sum(axis=1)
        else:
            demand_p2h = pd.Series(0, index=results['OutputPower'].index)
        if ('Flex', z) in inputs['param_df']['Demand']:
            demand_flex = inputs['param_df']['Demand'][('Flex', z)] 
        else:
            demand_flex = pd.Series(0, index=results['OutputPower'].index)
        demand_da = inputs['param_df']['Demand'][('DA', z)]
        demand[z] = pd.DataFrame(demand_da + demand_p2h + demand_flex, columns = [('DA', z)])
    demand = pd.concat(demand, axis=1)
    demand.columns = demand.columns.droplevel(-1)

    TotalLoad = demand.sum().sum()
    PeakLoad = demand.sum(axis=1).max(axis=0)
    LoadShedding = results['OutputShedLoad'].sum().sum() / 1e6
    Curtailment = results['OutputCurtailedPower'].sum().sum()
    MaxCurtailemnt = results['OutputCurtailedPower'].sum(axis=1).max() / 1e6
    MaxLoadShedding = results['OutputShedLoad'].sum(axis=1).max()
    
    if 'OutputDemandModulation' in results:
        ShiftedLoad_net = results['OutputDemandModulation'].sum().sum() / 1E6
        ShiftedLoad_tot = results['OutputDemandModulation'].abs().sum().sum()/2 /1E6
        if ShiftedLoad_net > 0.1 * ShiftedLoad_tot:
            logging.error('The net shifted load is higher than 10% of the total shifted load, although it should be zero')
    else:
        ShiftedLoad_tot = 0
    

    # TotalLoad = dfin['Demand']['DA'].loc[index, :].sum().sum()
    # # PeakLoad = inputs['parameters']['Demand']['val'][0,:,idx].sum(axis=0).max()
    # PeakLoad = dfin['Demand']['DA'].sum(axis=1).max(axis=0)

    NetImports = -get_imports(results['OutputFlow'], 'RoW')

    Cost_kwh = results['OutputSystemCost'].sum() / (TotalLoad - NetImports)

    print ('\nAverage electricity cost : ' + str(Cost_kwh) + ' EUR/MWh')

    for key in ['LostLoad_RampUp', 'LostLoad_2D', 'LostLoad_MinPower',
                'LostLoad_RampDown', 'LostLoad_2U', 'LostLoad_3U', 'LostLoad_MaxPower', 'LostLoad_WaterSlack']:
        if key == 'LostLoad_WaterSlack':
            if isinstance(results[key], pd.Series):
                LL = results[key].sum()
            else:
                LL = results[key]
        else:
            LL = results[key].values.sum()
        if LL > 0.0001 * TotalLoad:
            logging.critical('\nThere is a significant amount of lost load for ' + key + ': ' + str(
                LL) + ' MWh. The results should be checked carefully')
        elif LL > 100:
            logging.warning('\nThere is lost load for ' + key + ': ' + str(
                LL) + ' MWh. The results should be checked')
    
    print ('\nAggregated statistics for the considered area:')
    print ('Total Consumption:' + str(TotalLoad / 1E6) + ' TWh')
    print ('Peak Load:' + str(PeakLoad) + ' MW')
    print ('Net Importations:' + str(NetImports / 1E6) + ' TWh')
    print ('Total Load Shedding:' + str(LoadShedding) + ' TWh')
    print ('Total shifted load:' + str(ShiftedLoad_tot) + ' TWh')
    print ('Maximum Load Shedding:' + str(MaxLoadShedding) + ' MW')
    print ('Total Curtailed RES:' + str(Curtailment) + ' TWh')
    print ('Maximum Curtailed RES:' + str(MaxCurtailemnt) + ' MW')

    # Zone-specific values:
    ZoneData = pd.DataFrame(index=inputs['sets']['n'])

    
    if 'Flex' in dfin['Demand']:
        ZoneData['Flexible Demand'] = inputs['param_df']['Demand']['Flex'].sum(axis=0) / 1E6
        ZoneData['Demand'] = dfin['Demand']['DA'].sum(axis=0) / 1E6 + ZoneData['Flexible Demand']
        ZoneData['PeakLoad'] = (dfin['Demand']['DA']+dfin['Demand']['Flex']).max(axis=0) 
    else:
        ZoneData['PeakLoad'] = dfin['Demand']['DA'].max(axis=0) 
        ZoneData['Demand'] = dfin['Demand']['DA'].sum(axis=0) / 1E6

    ZoneData['NetImports'] = 0
    for z in ZoneData.index:
        ZoneData.loc[z, 'NetImports'] = get_imports(results['OutputFlow'], str(z)) / 1E6

    ZoneData['LoadShedding'] = results['OutputShedLoad'].sum(axis=0) / 1E6
    ZoneData['MaxLoadShedding'] = results['OutputShedLoad'].max()
    if 'OutputDemandModulation' in results:
        ZoneData['ShiftedLoad'] = results['OutputDemandModulation'].abs().sum() / 1E6
    ZoneData['Curtailment'] = results['OutputCurtailedPower'].sum(axis=0) / 1E6
    ZoneData['MaxCurtailment'] = results['OutputCurtailedPower'].max()

    print('\nZone-Specific values (in TWh or in MW):')
    print(ZoneData)

    # Congestion:
    Congestion = {}
    if 'OutputFlow' in results:
        for l in results['OutputFlow']:
            if l[:3] != 'RoW' and l[-3:] != 'RoW':
                Congestion[l] = np.sum(
                    (results['OutputFlow'][l] == dfin['FlowMaximum'].loc[results['OutputFlow'].index, l]) & (
                    dfin['FlowMaximum'].loc[results['OutputFlow'].index, l] > 0))
    print("\nNumber of hours of congestion on each line: ")
    import pprint
    pprint.pprint(Congestion)

    # Zone-specific storage data:
    try:
        StorageData = pd.DataFrame(index=inputs['sets']['n'])
        for z in StorageData.index:
            isstorage = pd.Series(index=inputs['units'].index)
            for u in isstorage.index:
                isstorage[u] = inputs['units'].Technology[u] in commons['tech_storage']
            sto_units = inputs['units'][(inputs['units'].Zone == z) & isstorage]
            StorageData.loc[z,'Storage Capacity [MWh]'] = (sto_units.Nunits*sto_units.StorageCapacity).sum()
            StorageData.loc[z,'Storage Power [MW]'] = (sto_units.Nunits*sto_units.PowerCapacity).sum()
            StorageData.loc[z,'Peak load shifting [hours]'] = StorageData.loc[z,'Storage Capacity [MWh]']/ZoneData.loc[z,'PeakLoad']
            AverageStorageOutput = 0
            for u in results['OutputPower'].columns:
                if u in sto_units.index:
                    AverageStorageOutput += results['OutputPower'][u].mean()
            StorageData.loc[z,'Average daily cycle depth [%]'] = AverageStorageOutput*24/(1e-9+StorageData.loc[z,'Storage Capacity [MWh]'])
        print('\nZone-Specific storage data')
        print(StorageData)
    except:
        logging.error('Could compute storage data')
        StorageData = None

    co2 = results['OutputPower'].sum() * inputs['param_df']['EmissionRate'] # MWh * tCO2 / MWh = tCO2
    co2.fillna(0,inplace=True)

    UnitData = pd.DataFrame(index=inputs['sets']['u'])
    UnitData.loc[:, 'Fuel'] = inputs['units']['Fuel']
    UnitData.loc[:, 'Technology'] = inputs['units']['Technology']
    UnitData.loc[:, 'Zone'] = inputs['units']['Zone']
    UnitData.loc[:, 'CHP'] = inputs['units']['CHPType']
    UnitData.loc[:, 'Generation [TWh]'] = results['OutputPower'].sum() / 1e6
    UnitData.loc[:, 'CO2 [t]'] = co2.loc['CO2',:]
    UnitData.loc[:, 'Total Costs [EUR]'] = get_units_operation_cost(inputs, results).sum()
    print('\nUnit-Specific data')
    print(UnitData)

    FuelData = {}
    chp = {'Extraction': 'CHP', 'back-pressure': 'CHP', 'P2H': 'CHP', '': 'Non-CHP'}
    tmp = UnitData
    tmp['CHP'] = tmp['CHP'].map(chp)
    for bo in ['CHP', 'Non-CHP']:
        tmp_data = tmp.loc[tmp['CHP'] == bo]
        FuelData[bo] = {}
        for l in ['Generation [TWh]','CO2 [t]','Total Costs [EUR]']:
            FuelData[bo][l] = pd.DataFrame(index=inputs['sets']['f'], columns=inputs['sets']['t'])
            for f in inputs['sets']['f']:
                for t in inputs['sets']['t']:
                    FuelData[bo][l].loc[f,t] = tmp_data.loc[(tmp_data['Fuel'] == f) &
                                                            (tmp_data['Technology'] == t)][l].sum()

    return {'Cost_kwh': Cost_kwh, 'TotalLoad': TotalLoad, 'PeakLoad': PeakLoad, 'NetImports': NetImports,
            'Curtailment': Curtailment, 'MaxCurtailment': MaxCurtailemnt, 'ShedLoad': LoadShedding,'ShiftedLoad':ShiftedLoad_tot,
            'MaxShedLoad': MaxLoadShedding, 'ZoneData': ZoneData, 'Congestion': Congestion, 'StorageData': StorageData,
            'UnitData': UnitData, 'FuelData': FuelData}


def get_indicators_powerplant(inputs, results):
    """
    Function that analyses the dispa-set results at the power plant level
    Computes the number of startups, the capacity factor, etc

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :returns out:        Dataframe with the main power plants characteristics and the computed indicators
    """
    out = inputs['units'].loc[:, ['Nunits','PowerCapacity', 'Zone', 'Technology', 'Fuel']]

    out['startups'] = 0
    for u in out.index:
        if u in results['OutputCommitted']:
            # count the number of start-ups
            values = results['OutputCommitted'].loc[:, u].values
            diff = -(values - np.roll(values, 1))
            startups = diff > 0
            out.loc[u, 'startups'] = startups.sum()

    out['CF'] = 0
    out['Generation'] = 0
    for u in out.index:
        if u in results['OutputPower']:
            # count the number of start-ups
            out.loc[u, 'CF'] = results['OutputPower'][u].mean() / (out.loc[u, 'PowerCapacity']*out.loc[u,'Nunits'])
            out.loc[u, 'Generation'] = results['OutputPower'][u].sum()
    return out


def CostExPost(inputs,results):
    '''
    Ex post computation of the operational costs with plotting. This allows breaking down
    the cost into its different components and check that it matches with the objective 
    function from the optimization.
    
    The cost objective function is the following:
         SystemCost(i)
         =E=
         sum(u,CostFixed(u)*TimeStep*Committed(u,i))
         +sum(u,CostStartUpH(u,i) + CostShutDownH(u,i))
         +sum(u,CostRampUpH(u,i) + CostRampDownH(u,i))
         +sum(u,CostVariable(u,i) * Power(u,i)*TimeStep)
         +sum(l,PriceTransmission(l,i)*Flow(l,i)*TimeStep)
         +sum(n,CostLoadShedding(n,i)*ShedLoad(n,i)*TimeStep)
         +sum(th, CostHeatSlack(th,i) * HeatSlack(th,i)*TimeStep)
         +sum(chp, CostVariable(chp,i) * CHPPowerLossFactor(chp) * Heat(chp,i)*TimeStep)
         +Config("ValueOfLostLoad","val")*(sum(n,(LL_MaxPower(n,i)+LL_MinPower(n,i))*TimeStep))
         +0.8*Config("ValueOfLostLoad","val")*(sum(n,(LL_2U(n,i)+LL_2D(n,i)+LL_3U(n,i))*TimeStep))
         +0.7*Config("ValueOfLostLoad","val")*sum(u,(LL_RampUp(u,i)+LL_RampDown(u,i))*TimeStep)
         +Config("CostOfSpillage","val")*sum(s,spillage(s,i));

             
    :returns: tuple[0]: dataframe with all different types of costs
    :returns: tuple[1]: cumulative sum of all costs time series
    :returns: tuple[2]: costs time series disaggregated by zone
    '''
    import datetime
    
    dfin = inputs['param_df']
    timeindex = results['OutputPower'].index
    config = inputs['config']
    TimeStep = config['SimulationTimeStep']
    VOLL = inputs['config']['default']['ValueOfLostLoad']
    
    
    # Define output dataframe:
    costs = pd.DataFrame(index=timeindex)
    
    # Add possible missing columns:"
    if 'OutputCommitted' not in results:
        results['OutputCommitted'] = pd.DataFrame(1,index=results['OutputPower'].index,columns=results['OutputPower'].columns)
    
    
    #%% Fixed Costs:
    costs['FixedCosts'] = 0
    for u in results['OutputCommitted']:
        if u in dfin['CostFixed'].index:
            costs['FixedCosts'] =+ dfin['CostFixed'].loc[u,'CostFixed'] * results['OutputCommitted'][u] * TimeStep
    
            
    #%% Ramping and startup costs:
    indexinitial = timeindex[0] - datetime.timedelta(hours=1)
    powerlong = results['OutputPower'].copy()
    powerlong.loc[indexinitial,:] = 0 
    powerlong.sort_index(inplace=True)
    committedlong = results['OutputCommitted'].copy()
    for u in powerlong:
        if u in dfin['PowerInitial'].index:
            powerlong.loc[indexinitial,u] = dfin['PowerInitial'].loc[u,'PowerInitial']
            committedlong.loc[indexinitial,u] = dfin['PowerInitial'].loc[u,'PowerInitial']>0
    committedlong.sort_index(inplace=True)
    committedlong.fillna(0,inplace=True)
            
    powerlong_shifted = powerlong.copy()
    powerlong_shifted.iloc[1:,:] = powerlong.iloc[:-1,:].values
    committedlong_shifted = committedlong.copy()
    committedlong_shifted.iloc[1:,:] = committedlong.iloc[:-1,:].values
    committedlong_shifted.fillna(0,inplace=True)
    ramping = powerlong - powerlong_shifted
   
    startups = committedlong.astype(int) - committedlong_shifted.astype(int)
    ramping.drop([ramping.index[0]],inplace=True); startups.drop([startups.index[0]],inplace=True)
    
    CostStartUp = pd.DataFrame(index=startups.index,columns=startups.columns)
    for u in CostStartUp:
        if u in dfin['CostStartUp'].index:
            CostStartUp[u] = startups[startups>0][u].fillna(0) * dfin['CostStartUp'].loc[u,'CostStartUp']
        else:
            print('Unit ' + u + ' not found in input table CostStartUp!')
            
    CostShutDown = pd.DataFrame(index=startups.index,columns=startups.columns)
    for u in CostShutDown:
        if u in dfin['CostShutDown'].index:
            CostShutDown[u] = startups[startups<0][u].fillna(0) * dfin['CostShutDown'].loc[u,'CostShutDown']
        else:
            print('Unit ' + u + ' not found in input table CostShutDown!')
            
    CostRampUp = pd.DataFrame(index=ramping.index,columns=ramping.columns)
    for u in CostRampUp:
        if u in dfin['CostRampUp'].index:
            CostRampUp[u] = ramping[ramping>0][u].fillna(0) * dfin['CostRampUp'].loc[u,'CostRampUp']
        else:
            print('Unit ' + u + ' not found in input table CostRampUp!')
            
    CostRampDown = pd.DataFrame(index=ramping.index,columns=ramping.columns)
    for u in CostRampDown:
        if u in dfin['CostRampDown'].index:
            CostRampDown[u] = ramping[ramping<0][u].fillna(0) * dfin['CostRampDown'].loc[u,'CostRampDown']
        else:
            print('Unit ' + u + ' not found in input table CostRampDown!')
    
    costs['CostStartUp'] = CostStartUp.sum(axis=1).fillna(0)
    costs['CostShutDown'] = CostShutDown.sum(axis=1).fillna(0)
    costs['CostRampUp'] = CostRampUp.sum(axis=1).fillna(0)
    costs['CostRampDown'] = CostRampDown.sum(axis=1).fillna(0)
    
    #%% Variable cost:
    costs['CostVariable'] = (results['OutputPower'] * dfin['CostVariable']).fillna(0).sum(axis=1) * TimeStep
    
    #%% Transmission cost:
    costs['CostTransmission'] = (results['OutputFlow'] * dfin['PriceTransmission']).fillna(0).sum(axis=1) * TimeStep
    
    #%% Shedding cost:
    costs['CostLoadShedding'] = (results['OutputShedLoad'] * dfin['CostLoadShedding']).fillna(0).sum(axis=1) * TimeStep
    
    #%% Heating costs:
    costs['CostHeatSlack'] = (results['OutputHeatSlack'] * dfin['CostHeatSlack']).fillna(0).sum(axis=1) * TimeStep
    CostHeat = pd.DataFrame(index=results['OutputHeat'].index,columns=results['OutputHeat'].columns)
    for u in CostHeat:
        if u in dfin['CHPPowerLossFactor'].index:
            CostHeat[u] = dfin['CostVariable'][u].fillna(0) * results['OutputHeat'][u].fillna(0) * dfin['CHPPowerLossFactor'].loc[u,'CHPPowerLossFactor'] * TimeStep
        else:
            CostHeat[u] = dfin['CostVariable'][u].fillna(0) * results['OutputHeat'][u].fillna(0) * TimeStep
    costs['CostHeat'] = CostHeat.sum(axis=1).fillna(0) 
    
    #%% Lost loads:
    # NB: the value of lost load is currently hard coded. This will have to be updated
    # Locate prices for LL
    #TODO:
    costs['LostLoad'] = 0.8 * VOLL * TimeStep * (results['LostLoad_2D'].reindex(timeindex).sum(axis=1).fillna(0) + results['LostLoad_2U'].reindex(timeindex).sum(axis=1).fillna(0) + results['LostLoad_3U'].reindex(timeindex).sum(axis=1).fillna(0))  \
                       + VOLL * TimeStep * (results['LostLoad_MaxPower'].reindex(timeindex).sum(axis=1).fillna(0) + results['LostLoad_MinPower'].reindex(timeindex).sum(axis=1).fillna(0)) \
                       + 0.7 * VOLL * TimeStep * (results['LostLoad_RampDown'].reindex(timeindex).sum(axis=1).fillna(0) + results['LostLoad_RampUp'].reindex(timeindex).sum(axis=1).fillna(0))
        
    #%% Spillage:
    costs['Spillage'] = config['default']['PriceOfSpillage'] * results['OutputSpillage'].sum(axis=1).fillna(0)
    
    #%% Plotting
    # Drop na columns:
    costs.dropna(axis=1, how='all',inplace=True)
    # Delete all-zero columns:
    # costs = costs.loc[:, (costs != 0).any(axis=0)]
    
    sumcost = costs.cumsum(axis=1)
    sumcost['OutputSystemCost'] = results['OutputSystemCost']
    
    sumcost.plot(title='Cumulative sum of the cost components')
    
    #%% Warning if significant error:
    diff = (costs.sum(axis=1) - results['OutputSystemCost']).abs()
    if diff.max() > 0.01 * results['OutputSystemCost'].max():
        logging.critical('There are significant differences between the cost computed ex post and and the cost provided by the optimization results!')

    #%% Compute cost per zone
    costs_zone = pd.DataFrame(index=timeindex)
    zones = inputs['sets']['n']
    for z in zones:
        # Assign transmission costs to the importing country:
        lines = [l for l in results['OutputFlow'].columns if l.endswith(z)]
        FlowToZone = results['OutputFlow'].loc[:, lines]
        ZoneTransmissionPrice = dfin['PriceTransmission'].loc[:,lines]
        # Add possibly missing columns:
        for key in ['OutputShedLoad','LostLoad_2D','LostLoad_2U','LostLoad_3U','LostLoad_MaxPower','LostLoad_MinPower']:
            if z not in results[key]:
                results[key][z] = 0
        FixedCosts = 0
        for u in filter_by_zone(results['OutputCommitted'],inputs,z):
            if u in dfin['CostFixed'].index:
                FixedCosts =+ dfin['CostFixed'].loc[u,'CostFixed'] * results['OutputCommitted'][u] * TimeStep
    
        costs_zone[z] = FixedCosts \
                        + filter_by_zone(CostStartUp,inputs,z).sum(axis=1).fillna(0)   \
                        + filter_by_zone(CostShutDown,inputs,z).sum(axis=1).fillna(0)   \
                        + filter_by_zone(CostRampUp,inputs,z).sum(axis=1).fillna(0)  \
                        + filter_by_zone(CostRampDown,inputs,z).sum(axis=1).fillna(0) \
                        + filter_by_zone((results['OutputPower'] * dfin['CostVariable']),inputs,z).fillna(0).sum(axis=1) * TimeStep \
                        + (FlowToZone * ZoneTransmissionPrice).fillna(0).sum(axis=1) * TimeStep \
                        + (results['OutputShedLoad'][z] * dfin['CostLoadShedding'][z]).fillna(0) * TimeStep \
                        + filter_by_zone((results['OutputHeatSlack'] * dfin['CostHeatSlack']),inputs,z).fillna(0).sum(axis=1) * TimeStep \
                        + filter_by_zone(CostHeat,inputs,z).sum(axis=1).fillna(0)  \
                        + 0.8 * VOLL * TimeStep * (results['LostLoad_2D'][z].reindex(timeindex).fillna(0) + results['LostLoad_2U'][z].reindex(timeindex).fillna(0) + results['LostLoad_3U'][z].reindex(timeindex).fillna(0))  \
                        + VOLL * TimeStep * (results['LostLoad_MaxPower'][z].reindex(timeindex).fillna(0) + results['LostLoad_MinPower'][z].reindex(timeindex).fillna(0)) \
                        + 0.7 * VOLL * TimeStep * filter_by_zone((results['LostLoad_RampDown'].reindex(timeindex).fillna(0) + results['LostLoad_RampUp'].reindex(timeindex).fillna(0)),inputs,z).sum(axis=1) \
                        + config['default']['PriceOfSpillage'] * filter_by_zone(results['OutputSpillage'],inputs,z).sum(axis=1).fillna(0)

    return costs,sumcost,costs_zone


def get_units_operation_cost(inputs, results):
    """
    Function that computes the operation cost for each power unit at each instant of time from the DispaSET results
    Operation cost includes: CostFixed + CostStartUp + CostShutDown + CostRampUp + CostRampDown + CostVariable

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :returns out:       Dataframe with the the power units in columns and the operation cost at each instant in rows

    Main Author: @AbdullahAlawad
    """

    datain = ds_to_df(inputs)

    #DataFrame with startup times for each unit (1 for startup)
    StartUps = results['OutputCommitted'].copy()
    for u in StartUps:
        values = StartUps.loc[:, u].values
        diff = -(np.roll(values, 1) - values )
        diff[diff <= 0] = 0
        StartUps[u] = diff

    #DataFrame with shutdown times for each unit (1 for shutdown)
    ShutDowns = results['OutputCommitted'].copy()
    for u in ShutDowns:
        values = ShutDowns.loc[:, u].values
        diff = (np.roll(values, 1) - values )
        diff[diff <= 0] = 0
        ShutDowns[u] = diff

    #DataFrame with ramping up levels for each unit at each instant (0 for ramping-down & leveling out)
    RampUps = results['OutputPower'].copy()
    for u in RampUps:
        values = RampUps.loc[:, u].values
        diff = -(np.roll(values, 1) - values )
        diff[diff <= 0] = 0
        RampUps[u] = diff

    #DataFrame with ramping down levels for each unit at each instant (0 for ramping-up & leveling out)
    RampDowns = results['OutputPower'].copy()
    for u in RampDowns:
        values = RampDowns.loc[:, u].values
        diff = (np.roll(values, 1) - values )
        diff[diff <= 0] = 0
        RampDowns[u] = diff

    FiexedCost = results['OutputCommitted'].copy()
    StartUpCost = results['OutputCommitted'].copy()
    ShutDownCost = results['OutputCommitted'].copy()
    RampUpCost = results['OutputCommitted'].copy()
    RampDownCost = results['OutputCommitted'].copy()
    VariableCost = results['OutputCommitted'].copy()
    UnitOperationCost = results['OutputCommitted'].copy()

    OperatedUnitList = results['OutputCommitted'].columns
    for u in OperatedUnitList:
        if u not in results['OutputPower']:
            results['OutputPower'][u]=0
            RampUps[u] = 0
            RampDowns[u] = 0
            StartUps[u] = 0
            ShutDowns[u] = 0
            logging.warning('Unit ' + u + ' is in the table committed but not in the output power table')
        unit_indexNo = inputs['units'].index.get_loc(u)
        FiexedCost.loc[:,[u]] = np.array(results['OutputCommitted'].loc[:,[u]])*inputs['parameters']['CostFixed']['val'][unit_indexNo]
        StartUpCost.loc[:,[u]] = np.array(StartUps.loc[:,[u]])*inputs['parameters']['CostStartUp']['val'][unit_indexNo]
        ShutDownCost.loc[:,[u]] = np.array(ShutDowns.loc[:,[u]])*inputs['parameters']['CostShutDown']['val'][unit_indexNo]
        RampUpCost.loc[:,[u]] = np.array(RampUps.loc[:,[u]])*inputs['parameters']['CostRampUp']['val'][unit_indexNo]
        RampDownCost.loc[:,[u]] = np.array(RampDowns.loc[:,[u]])*inputs['parameters']['CostRampDown']['val'][unit_indexNo]
        VariableCost.loc[:,[u]] = np.array(datain['CostVariable'].loc[:,[u]])*np.array(results['OutputPower'][u]).reshape(-1,1)

    UnitOperationCost = FiexedCost+StartUpCost+ShutDownCost+RampUpCost+RampDownCost+VariableCost

    return UnitOperationCost
