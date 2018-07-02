# -*- coding: utf-8 -*-
"""
Set of functions useful to analyse to DispaSET output data.

@author: Sylvain Quoilin, JRC
"""

from __future__ import division

import logging
import sys

import numpy as np
import pandas as pd

from .data_handler import commons


# get color definitions:


def get_load_data(inputs, c):
    """ 
    Get the load curve, the residual load curve, and the net residual load curve of a specific country

    :param inputs:  DispaSET inputs (output of the get_sim_results function)
    :param c:       Country to consider (e.g. 'BE')
    :return out:    Dataframe with the following columns:
                        Load:               Load curve of the specified country
                        ResidualLoad:       Load minus the production of variable renewable sources
                        NetResidualLoad:    Residual netted from the interconnections with neightbouring countries
    """
    datain = inputs['param_df']
    out = pd.DataFrame(index=datain['Demand'].index)
    out['Load'] = datain['Demand']['DA', c]
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
        if l[:2] == c:
            Interconnections = Interconnections - datain['FlowMinimum'][l].values
        elif l[-2:] == c:
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
            logging.warn('Fuel not found for unit ' + u + ' with fuel ' + uFuel[u])

    return PowerByFuel


def filter_by_country(PowerOutput, inputs, c):
    """
    This function filters the dispaset Output Power dataframe by country

    :param PowerOutput:     Dataframe of power generationwith units as columns and time as index
    :param Inputs:          Dispaset inputs version 2.1.1
    :param c:               Selected country (e.g. 'BE')
    :returns Power:          Dataframe with power generation by country
    """
    loc = inputs['units']['Zone']
    Power = PowerOutput.loc[:, [u for u in PowerOutput.columns if loc[u] == c]]
    return Power


def get_plot_data(inputs, results, c):
    """
    Function that reads the results dataframe of a DispaSET simulation and extract the dispatch data spedific to one country

    :param results:         Pandas dataframe with the results (output of the GdxToDataframe function)
    :param c:               Country to be considered (e.g. 'BE')
    :returns plotdata:       Dataframe with the dispatch data storage and outflows are negative
    """
    tmp = filter_by_country(results['OutputPower'], inputs, c)
    plotdata = aggregate_by_fuel(tmp, inputs)

    if 'OutputStorageInput' in results:
        #onnly take the columns that correspond to storage units (StorageInput is also used for CHP plants):
        cols = [col for col in results['OutputStorageInput'] if inputs['units'].loc[col,'Technology'] in commons['tech_storage']]
        tmp = filter_by_country(results['OutputStorageInput'][cols], inputs, c)
        plotdata['Storage'] = -tmp.sum(axis=1)
    else:
        plotdata['Storage'] = 0
    plotdata.fillna(value=0, inplace=True)

    plotdata['FlowIn'] = 0
    plotdata['FlowOut'] = 0
    for col in results['OutputFlow']:
        if col[-2:] == c:
            plotdata['FlowIn'] = plotdata['FlowIn'] + results['OutputFlow'][col]
    for col in results['OutputFlow']:
        if col[:2] == c:
            plotdata['FlowOut'] = plotdata['FlowOut'] - results['OutputFlow'][col]

    # re-ordering columns:
    OrderedColumns = [col for col in commons['MeritOrder'] if col in plotdata.columns]
    plotdata = plotdata[OrderedColumns]
    
    # remove empty columns:
    for col in plotdata.columns:
        if plotdata[col].max() == 0 and plotdata[col].min()==0:
            del plotdata[col]

    return plotdata


def get_imports(flows, c):
    """ 
    Function that computes the balance of the imports/exports of a given zone

    :param flows:       Pandas dataframe with the timeseries of the exchanges
    :param c:           Country (zone) to consider
    :returns NetImports: Scalar with the net balance over the whole time period
    """
    NetImports = 0
    for key in flows:
        if key[:len(c)] == c:
            NetImports -= flows[key].sum()
        elif key[-len(c):] == c:
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
    index = pd.DatetimeIndex(start=pd.datetime(*StartDate), end=pd.datetime(*StopDate), freq='h')

    # Aggregated values:
    TotalLoad = dfin['Demand']['DA'].loc[index, :].sum().sum()
    # PeakLoad = inputs['parameters']['Demand']['val'][0,:,idx].sum(axis=0).max()
    PeakLoad = dfin['Demand']['DA'].sum(axis=1).max(axis=0)

    NetImports = -get_imports(results['OutputFlow'], 'RoW')

    Cost_kwh = results['OutputSystemCost'].sum() / (TotalLoad - NetImports)

    print ('\nAverage electricity cost : ' + str(Cost_kwh) + ' EUR/MWh')

    for key in ['LostLoad_RampUp', 'LostLoad_2D', 'LostLoad_MinPower',
                'LostLoad_RampDown', 'LostLoad_2U', 'LostLoad_3U', 'LostLoad_MaxPower']:
        LL = results[key].values.sum()
        if LL > 0.0001 * TotalLoad:
            logging.critical('\nThere is a significant amount of lost load for ' + key + ': ' + str(
                LL) + ' MWh. The results should be checked carefully')
        elif LL > 100:
            logging.warning('\nThere is lost load for ' + key + ': ' + str(
                LL) + ' MWh. The results should be checked')

    print ('\nAggregated statistics for the considered area:')
    print ('Total consumption:' + str(TotalLoad / 1E6) + ' TWh')
    print ('Peak load:' + str(PeakLoad) + ' MW')
    print ('Net importations:' + str(NetImports / 1E6) + ' TWh')

    # Country-specific values:
    CountryData = pd.DataFrame(index=inputs['sets']['n'])

    CountryData['Demand'] = dfin['Demand']['DA'].sum(axis=0) / 1E6
    CountryData['PeakLoad'] = dfin['Demand']['DA'].max(axis=0)

    CountryData['NetImports'] = 0
    for c in CountryData.index:
        CountryData.loc[c, 'NetImports'] = get_imports(results['OutputFlow'], str(c)) / 1E6

    CountryData['LoadShedding'] = results['OutputShedLoad'].sum(axis=0) / 1E6
    CountryData['Curtailment'] = results['OutputCurtailedPower'].sum(axis=0) / 1E6
    print('\nCountry-Specific values (in TWh or in MW):')
    print(CountryData)

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

    # Country-specific storage data:
    try:
        StorageData = pd.DataFrame(index=inputs['sets']['n'])
        List_tech_storage = ['HDAM', 'HPHS', 'BATS', 'BEVS', 'CAES', 'THMS']
        for c in StorageData.index:
            isstorage = pd.Series(index=inputs['units'].index)
            for u in isstorage.index:
                isstorage[u] = inputs['units'].Technology[u] in List_tech_storage
            sto_units = inputs['units'][(inputs['units'].Zone == c) & isstorage]
            StorageData.loc[c,'Storage Capacity [MWh]'] = (sto_units.Nunits*sto_units.StorageCapacity).sum()
            StorageData.loc[c,'Storage Power [MW]'] = (sto_units.Nunits*sto_units.PowerCapacity).sum()
            StorageData.loc[c,'Peak load shifting [hours]'] = StorageData.loc[c,'Storage Capacity [MWh]']/CountryData.loc[c,'PeakLoad']
            AverageStorageOutput = 0
            for u in results['OutputPower'].columns:
                if u in sto_units.index:
                    AverageStorageOutput += results['OutputPower'][u].mean()
            StorageData.loc[c,'Average daily cycle depth [%]'] = AverageStorageOutput*24/StorageData.loc[c,'Storage Capacity [MWh]']
        print('\nCountry-Specific storage data')
        print(StorageData)
    except:
        logging.error('Could compute storage data')
        StorageData = None

    return {'Cost_kwh': Cost_kwh, 'TotalLoad': TotalLoad, 'PeakLoad': PeakLoad, 'NetImports': NetImports,
            'CountryData': CountryData, 'Congestion': Congestion, 'StorageData': StorageData}    

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
             sum(u,CostFixed(u)*Committed(u,i))
             +sum(u,CostStartUpH(u,i) + CostShutDownH(u,i))
             +sum(u,CostRampUpH(u,i) + CostRampDownH(u,i))
             +sum(u,CostVariable(u,i) * Power(u,i))
             +sum(l,PriceTransmission(l,i)*Flow(l,i))
             +sum(n,CostLoadShedding(n,i)*ShedLoad(n,i))
             +sum(chp, CostHeatSlack(chp,i) * HeatSlack(chp,i))
             +sum(chp, CostVariable(chp,i) * CHPPowerLossFactor(chp) * Heat(chp,i))
             +100E3*(sum(n,LL_MaxPower(n,i)+LL_MinPower(n,i)))
             +80E3*(sum(n,LL_2U(n,i)+LL_2D(n,i)+LL_3U(n,i)))
             +70E3*sum(u,LL_RampUp(u,i)+LL_RampDown(u,i))
             +1*sum(s,spillage(s,i))
             
    :returns: tuple with the cost components and their cumulative sums in two dataframes.
    '''
    import datetime
    
    dfin = inputs['param_df']
    timeindex = results['OutputPower'].index
    
    costs = pd.DataFrame(index=timeindex)
    
    #%% Fixed Costs:
    costs['FixedCosts'] = 0
    for u in results['OutputCommitted']:
        if u in dfin:
            costs['FixedCosts'] =+ dfin[u] * results['OutputCommited'][u]
    
            
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
            
    powerlong_shifted = powerlong.copy()
    powerlong_shifted.iloc[1:,:] = powerlong.iloc[:-1,:].values
    committedlong_shifted = committedlong.copy()
    committedlong_shifted.iloc[1:,:] = committedlong.iloc[:-1,:].values
    
    ramping = powerlong - powerlong_shifted
    startups = committedlong - committedlong_shifted
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
    costs['CostVariable'] = (results['OutputPower'] * dfin['CostVariable']).fillna(0).sum(axis=1)
    
    #%% Transmission cost:
    costs['CostTransmission'] = (results['OutputFlow'] * dfin['PriceTransmission']).fillna(0).sum(axis=1)
    
    #%% Shedding cost:
    costs['CostLoadShedding'] = (results['OutputShedLoad'] * dfin['CostLoadShedding']).fillna(0).sum(axis=1)
    
    #%% Heating costs:
    costs['CostHeatSlack'] = (results['OutputHeatSlack'] * dfin['CostHeatSlack']).fillna(0).sum(axis=0)
    CostHeat = (results['OutputHeatSlack'] * dfin['CostHeatSlack']).fillna(0)
    CostHeat = pd.DataFrame(index=results['OutputHeat'].index,columns=results['OutputHeat'].columns)
    for u in CostHeat:
        if u in dfin['CHPPowerLossFactor'].index:
            CostHeat[u] = dfin['CostVariable'][u].fillna(0) * results['OutputHeat'][u].fillna(0) * dfin['CHPPowerLossFactor'].loc[u,'CHPPowerLossFactor']
        else:
            print('Unit ' + u + ' not found in input table CHPPowerLossFactor!')
    costs['CostHeat'] = CostHeat.sum(axis=1).fillna(0)
    
    #%% Lost loads:
    costs['LostLoad'] = 80e3* (results['LostLoad_2D'].reindex(timeindex).sum(axis=1).fillna(0) + results['LostLoad_2U'].reindex(timeindex).sum(axis=1).fillna(0) + results['LostLoad_3U'].reindex(timeindex).sum(axis=1).fillna(0))  \
                       + 100e3*(results['LostLoad_MaxPower'].reindex(timeindex).sum(axis=1).fillna(0) + results['LostLoad_MinPower'].reindex(timeindex).sum(axis=1).fillna(0)) \
                       + 70e3*(results['LostLoad_RampDown'].reindex(timeindex).sum(axis=1).fillna(0) + results['LostLoad_RampUp'].reindex(timeindex).sum(axis=1).fillna(0))
        
    #%% Spillage:
    costs['Spillage'] = 1 * results['OutputSpillage'].sum(axis=1).fillna(0)
    
    #%% Plotting
    # Drop na columns:
    costs.dropna(axis=1, how='all',inplace=True)
    # Delete all-zero columns:
    costs = costs.loc[:, (costs != 0).any(axis=0)]
    
    sumcost = costs.cumsum(axis=1)
    sumcost['OutputSystemCost'] = results['OutputSystemCost']
    
    sumcost.plot(title='Cumulative sum of the cost components')
    
    #%% Warning if significant error:
    diff = (costs.sum(axis=1) - results['OutputSystemCost']).abs()
    if diff.max() > 0.01 * results['OutputSystemCost'].max():
        logging.critical('There are significant differences between the cost computed ex post and and the cost provided by the optimization results!')
    return costs,sumcost


