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

from scipy.integrate import odeint
import matplotlib.pyplot as plt


from ..common import commons


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
    for line in datain['FlowMinimum']:
        if line[:2] == z:
            Interconnections = Interconnections - datain['FlowMinimum'][line].values
        elif line[-2:] == z:
            Interconnections = Interconnections + datain['FlowMinimum'][line].values
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

def filter_sector(InputParam, inputs):
    loc = inputs['units']['Zone']
    InputParamCopy = InputParam.copy()
    if str(InputParamCopy.columns[0]).startswith('S_'):
        for s in loc[inputs['units']['Technology'] == 'HDAMC'].index:
           InputParamCopy = InputParamCopy.rename(columns={inputs['units']['Sector1'][s]: s})
        InputParamCopy = InputParamCopy[InputParamCopy.columns.intersection(loc.index.to_list())]
    elif str(InputParamCopy.index[0]).startswith('S_'):
        for s in loc[inputs['units']['Technology'] == 'HDAMC'].index:
            InputParamCopy = InputParamCopy.rename(index={inputs['units']['Sector1'][s]: s})
        InputParamCopy = InputParamCopy.loc[InputParamCopy.index.intersection(loc.index.to_list())]
    return InputParamCopy

def filter_by_zone(PowerOutput, inputs, z, thermal = None, sector = False):
    """
    This function filters the dispaset Output dataframes by zone

    This function filters the dispaset Output Power dataframe by zone

    :param PowerOutput:     Dataframe of power generation with units as columns and time as index
    :param inputs:          Dispaset inputs version 2.1.1
    :param z:               Selected zone (e.g. 'BE')
    :returns Power:         Dataframe with power generation by zone
    """
    if thermal:
        loc = inputs['units']['Zone_th']
        PowerOutputCopy = PowerOutput.copy()
        Power = PowerOutputCopy.loc[:, [u for u in PowerOutputCopy.columns if loc[u] == z]]
    else:
        if sector != True:
            loc = inputs['units']['Zone']
            PowerOutputCopy = PowerOutput.copy()
            Power = PowerOutputCopy.loc[:, [u for u in PowerOutputCopy.columns if loc[u] == z]]
        if sector == True:
            loc = inputs['units'][['Zone', 'Sector1']]
            loc = loc[~loc['Sector1'].str.contains('nan')].dropna(how='any').drop_duplicates()
            result = loc.groupby('Sector1', as_index=False).first()[['Zone', 'Sector1']]
            result.set_index('Zone', inplace=True)
            value = result.loc[z, 'Sector1']
            # Check if the value is a pandas Series
            if isinstance(value, pd.Series):
                indices = value.tolist()
            elif isinstance(value, str):
                indices = [value]
            else:
                indices = []
            PowerOutputCopy = PowerOutput.copy()
            Power = PowerOutputCopy.loc[:, [u for u in PowerOutputCopy.columns if u in indices]]
            if PowerOutputCopy.empty:
                Power = pd.DataFrame()
    return Power


def filter_by_heating_zone(HeatOutput, inputs, z_th):
    """
    This function filters the dispaset Output Power dataframe by zone

    :param HeatOutput:     Dataframe of power generationwith units as columns and time as index
    :param inputs:          Dispaset inputs version 2.1.1
    :param z:               Selected zone (e.g. 'BE')
    :returns Heat:          Dataframe with power generation by zone
    """
    loc = inputs['units']['Zone_th']
    Heat = HeatOutput.loc[:, [u for u in HeatOutput.columns if loc[u] == z_th]]
    return Heat


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


def filter_by_tech_list(PowerOutput, inputs, tech):
    """
    This function filters the dispaset power output dataframe by technology

    :param PowerOutput:   Dataframe of power generation with units as columns and time as index
    :param inputs:        Dispaset inputs version 2.1.1
    :param t:             Selected tech (e.g. 'HDAM')
    :returns Power:
    """
    loc = inputs['units']['Technology']
    Power = pd.DataFrame()
    for t in tech:
        tmp = PowerOutput.loc[:, [u for u in PowerOutput.columns if loc[u] == t]]
        Power = pd.concat([Power, tmp], axis=1)
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
    Function that reads the results dataframe of a DispaSET simulation
    and extract the dispatch data specific to one zone

    :param results:         Pandas dataframe with the results (output of the GdxToDataframe function)
    :param z:               Zone to be considered (e.g. 'BE')
    :returns plotdata:      Dataframe with the dispatch data storage and outflows are negative
    """
    tmp = filter_by_zone(results['OutputPower'], inputs, z)
    plotdata = aggregate_by_fuel(tmp, inputs)

    if 'OutputStorageInput' in results:
        # onnly take the columns that correspond to storage units (StorageInput is also used for CHP plants):
        cols = [col for col in results['OutputStorageInput'] if
                inputs['units'].loc[col, 'Technology'] in commons['tech_storage']]
        tmp = filter_by_zone(results['OutputStorageInput'][cols], inputs, z)
        bb = pd.DataFrame()
        for tech in commons['tech_storage']:
            aa = filter_by_tech(tmp, inputs, tech)
            aa = aa.sum(axis=1)
            aa = pd.DataFrame(aa, columns=[tech])
            bb = pd.concat([bb, aa], axis=1)
        bb = -bb
        plotdata = pd.concat([plotdata, bb], axis=1)
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
        if plotdata[col].max() == 0 and plotdata[col].min() == 0 and col not in ['FlowIn', 'FlowOut']:
            del plotdata[col]

    return plotdata


def get_heat_plot_data(inputs, results, z_th):
    """
    Function that reads the results dataframe of a DispaSET simulation
    and extract the dispatch data specific to one zone

    :param results:         Pandas dataframe with the results (output of the GdxToDataframe function)
    :param z_th:            Heating zone to be considered (e.g. 'BE_th')
    :returns plotdata:      Dataframe with the dispatch data storage and outflows are negative
    """
    tmp = filter_by_heating_zone(results['OutputHeat'], inputs, z_th)
    plotdata = aggregate_by_fuel(tmp, inputs)


    if z_th not in results['OutputHeatSlack'].columns:
        plotdata.loc[:, 'HeatSlack'] = 0
    else:
        plotdata.loc[:, 'HeatSlack'] = results['OutputHeatSlack'].loc[:, z_th]

    plotdata.fillna(value=0, inplace=True)

    # re-ordering columns:
    OrderedColumns = [col for col in commons['MeritOrderHeat'] if col in plotdata.columns]
    plotdata = plotdata[OrderedColumns]

    # remove empty columns:
    for col in plotdata.columns:
        if plotdata[col].max() == 0 and plotdata[col].min() == 0:
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
def get_result_analysis(inputs, results, units='MWh'):
    """
    Reads the DispaSET results and provides useful general information to stdout

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    """

    if units == 'MWh':
        unit_small = ['MW','MWh']
        unit_medium = ['GW','GWh']
        unit_large = ['TW', 'TWh']
    elif units == 'kWh':
        unit_small = ['kW', 'MWh']
        unit_medium = ['MW', 'GWh']
        unit_large = ['GW', 'TWh']

    # inputs into the dataframe format:
    dfin = inputs['param_df']

    # Aggregated values:
    demand_zone = pd.DataFrame()
    for z in inputs['sets']['n']:
        if 'OutputPowerConsumption' in results:
            demand_p2x = filter_by_zone(results['OutputPowerConsumption'], inputs, z)
            demand_p2x = demand_p2x.sum(axis=1)
        else:
            demand_p2x = pd.Series(0, index=results['OutputPower'].index)
        if ('Flex', z) in inputs['param_df']['Demand']:
            demand_flex = inputs['param_df']['Demand'][('Flex', z)]
        else:
            demand_flex = pd.Series(0, index=results['OutputPower'].index)
        demand_da = inputs['param_df']['Demand'][('DA', z)]
        demand_temp = pd.DataFrame(demand_da + demand_p2x + demand_flex, columns=[z])
        demand_zone = pd.concat([demand_zone, demand_temp], axis=1)

    demand_total = demand_zone.sum(axis=1)

    TotalLoad = demand_total.sum().sum()
    PeakLoad = demand_total.max()
    LoadShedding = results['OutputShedLoad'].sum().sum() / 1e6
    Curtailment = results['OutputCurtailedPower'].sum().sum()
    MaxCurtailemnt = results['OutputCurtailedPower'].sum(axis=1).max() / 1e6
    MaxLoadShedding = results['OutputShedLoad'].sum(axis=1).max()

    if 'OutputDemandModulation' in results:
        ShiftedLoad_net = results['OutputDemandModulation'].sum().sum() / 1E6
        ShiftedLoad_tot = results['OutputDemandModulation'].abs().sum().sum() / 2 / 1E6
        if ShiftedLoad_net > 0.1 * ShiftedLoad_tot:
            logging.error(
                'The net shifted load is higher than 10% of the total shifted load, although it should be zero')
    else:
        ShiftedLoad_tot = 0

    if 'ShadowPrice' in results:
        if results['ShadowPrice'].isnull().values.any():
            logging.warning('Shadow Prices are not computed properly. DataFrame has nan values')
        else:
            Cost_kwh_zone = (demand_zone * results['ShadowPrice']).sum(axis=0) / demand_zone.sum(axis=0)
            Cost_kwh = (demand_zone * results['ShadowPrice']).sum(axis=1).sum() / demand_zone.sum(axis=0).sum()
            print('\nAverage electricity cost (EUR/' + unit_small[1] + '): \n'  + str(Cost_kwh_zone) + '\nEntireRegion ' + str(Cost_kwh))

    for key in ['LostLoad_RampUp', 'LostLoad_2D', 'LostLoad_MinPower',
                'LostLoad_RampDown', 'LostLoad_2U', 'LostLoad_3U', 'LostLoad_MaxPower',
                'LostLoad_StorageLevelViolation']:
        if key == 'LostLoad_StorageLevelViolation':
            if isinstance(results[key], pd.Series):
                LL = results[key].sum()
            else:
                LL = results[key]
        else:
            LL = results[key].values.sum()
        if LL > 0.0001 * TotalLoad:
            logging.critical('\nThere is a significant amount of lost load for ' + key + ': ' + str(
                LL) + unit_small[1] + '. The results should be checked carefully')
        elif LL > 100:
            logging.warning('\nThere is lost load for ' + key + ': ' + str(
                LL) + unit_small[1] + '. The results should be checked')

    NetImports = -get_imports(results['OutputFlow'], 'RoW')

    print('\nAggregated statistics for the considered area:')
    print('Total Consumption:' + str(TotalLoad / 1E6) + ' ' + unit_large[1])
    print('Peak Load:' + str(PeakLoad) + ' ' + unit_small[0])
    print('Net Importations:' + str(NetImports / 1E6) + ' ' + unit_large[1])
    print('Total Load Shedding:' + str(LoadShedding) + ' ' + unit_large[1])
    print('Total shifted load:' + str(ShiftedLoad_tot) + ' ' + unit_large[1])
    print('Maximum Load Shedding:' + str(MaxLoadShedding) + ' ' + unit_small[0])
    print('Total Curtailed RES:' + str(Curtailment) + ' ' + unit_large[1])
    print('Maximum Curtailed RES:' + str(MaxCurtailemnt) + ' ' + unit_small[0])

    # Zone-specific values:
    ZoneData = pd.DataFrame(index=inputs['sets']['n'])

    if 'Flex' in dfin['Demand']:
        ZoneData['Flexible Demand'] = inputs['param_df']['Demand']['Flex'].sum(axis=0) / 1E6
        ZoneData['Demand'] = dfin['Demand']['DA'].sum(axis=0) / 1E6 + ZoneData['Flexible Demand']
        ZoneData['PeakLoad'] = (dfin['Demand']['DA'] + dfin['Demand']['Flex']).max(axis=0)
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
        for flow in results['OutputFlow']:
            if flow[:3] != 'RoW' and flow[-3:] != 'RoW':
                Congestion[flow] = np.sum(
                    (results['OutputFlow'][flow] == dfin['FlowMaximum'].loc[results['OutputFlow'].index, flow]) & (
                            dfin['FlowMaximum'].loc[results['OutputFlow'].index, flow] > 0))
    print("\nNumber of hours of congestion on each line: ")
    import pprint
    pprint.pprint(Congestion)

    # Zone-specific storage data:
    try:
        StorageData = pd.DataFrame(index=inputs['sets']['n'])
        for z in StorageData.index:
            isstorage = pd.Series(index=inputs['units'].index, dtype='float64')
            for u in isstorage.index:
                isstorage[u] = inputs['units'].Technology[u] in commons['tech_storage']
            sto_units = inputs['units'][(inputs['units'].Zone == z) & isstorage]
            StorageData.loc[z, 'Storage Capacity [MWh]'] = (sto_units.Nunits * sto_units.StorageCapacity).sum()
            StorageData.loc[z, 'Storage Power [MW]'] = (sto_units.Nunits * sto_units.PowerCapacity).sum()
            StorageData.loc[z, 'Peak load shifting [hours]'] = StorageData.loc[z, 'Storage Capacity [MWh]'] / \
                                                               ZoneData.loc[z, 'PeakLoad']
            AverageStorageOutput = 0
            for u in results['OutputPower'].columns:
                if u in sto_units.index:
                    AverageStorageOutput += results['OutputPower'][u].mean()
            StorageData.loc[z, 'Average daily cycle depth [%]'] = AverageStorageOutput * 24 / (
                    1e-9 + StorageData.loc[z, 'Storage Capacity [MWh]'])
        print('\nZone-Specific storage data')
        print(StorageData)
    except:
        logging.error('Couldnt compute storage data')
        StorageData = None

    co2 = results['OutputPower'].sum() * inputs['param_df']['EmissionRate']  # MWh * tCO2 / MWh = tCO2
    co2.fillna(0, inplace=True)

    UnitData = pd.DataFrame(index=inputs['sets']['u'])
    UnitData.loc[:, 'Fuel'] = inputs['units']['Fuel']
    UnitData.loc[:, 'Technology'] = inputs['units']['Technology']
    UnitData.loc[:, 'Zone'] = inputs['units']['Zone']
    UnitData.loc[:, 'CHP'] = inputs['units']['CHPType']
    UnitData.loc[:, 'Generation [TWh]'] = results['OutputPower'].sum() / 1e6
    UnitData.loc[:, 'CO2 [t]'] = co2.loc['CO2', :]
    UnitData.loc[:, 'Total Costs [EUR]'] = get_units_operation_cost(inputs, results).sum(axis=1)
    UnitData.loc[:, 'WaterWithdrawal'] = results['OutputPower'].sum() * \
                                         inputs['units'].loc[:, 'WaterWithdrawal'].fillna(0)
    UnitData.loc[:, 'WaterConsumption'] = results['OutputPower'].sum() * \
                                          inputs['units'].loc[:, 'WaterConsumption'].fillna(0)
    print('\nUnit-Specific data')
    print(UnitData)

    FuelData = {}
    chp = {'Extraction': 'CHP', 'back-pressure': 'CHP', 'P2H': 'CHP', '': 'Non-CHP'}
    tmp = UnitData
    tmp['CHP'] = tmp['CHP'].map(chp)
    for bo in ['CHP', 'Non-CHP']:
        tmp_data = tmp.loc[tmp['CHP'] == bo]
        FuelData[bo] = {}
        for var in ['Generation [TWh]', 'CO2 [t]', 'Total Costs [EUR]']:
            FuelData[bo][var] = pd.DataFrame(index=inputs['sets']['f'], columns=inputs['sets']['t'])
            for f in inputs['sets']['f']:
                for t in inputs['sets']['t']:
                    FuelData[bo][var].loc[f, t] = tmp_data.loc[(tmp_data['Fuel'] == f) &
                                                               (tmp_data['Technology'] == t)][var].sum()

    WaterData = {}
    WaterData['UnitLevel'] = {}
    WaterData['ZoneLevel'] = {}
    WaterData['UnitLevel']['WaterWithdrawal'] = results['OutputPower'] * inputs['units'].loc[:, 'WaterWithdrawal']
    WaterData['UnitLevel']['WaterConsumption'] = results['OutputPower'] * inputs['units'].loc[:, 'WaterConsumption']
    WaterData['ZoneLevel']['WaterWithdrawal'] = pd.DataFrame()
    WaterData['ZoneLevel']['WaterConsumption'] = pd.DataFrame()
    for z in inputs['sets']['n']:
        tmp = pd.DataFrame(filter_by_zone(WaterData['UnitLevel']['WaterWithdrawal'], inputs, z).sum(axis=1),
                           columns=[z])
        WaterData['ZoneLevel']['WaterWithdrawal'] = pd.concat([WaterData['ZoneLevel']['WaterWithdrawal'], tmp], axis=1)
        tmp = pd.DataFrame(filter_by_zone(WaterData['UnitLevel']['WaterConsumption'], inputs, z).sum(axis=1),
                           columns=[z])
        WaterData['ZoneLevel']['WaterConsumption'] = pd.concat([WaterData['ZoneLevel']['WaterConsumption'], tmp],
                                                               axis=1)

    # if not inputs['param_df']['HeatDemand'].empty:
    #     return {'Cost_kwh': Cost_kwh, 'TotalLoad': TotalLoad, 'PeakLoad': PeakLoad, 'NetImports': NetImports,
    #             'TotalHeatLoad': TotalHeatLoad, 'PeakHeatLoad': PeakHeatLoad,
    #             'Curtailment': Curtailment, 'MaxCurtailment': MaxCurtailemnt,
    #             'HeatCurtailment': HeatCurtailment, 'MaxHeatCurtailment': MaxHeatCurtailemnt,
    #             'ShedLoad': LoadShedding, 'MaxShedLoad': MaxLoadShedding,
    #             'ShiftedLoad': ShiftedLoad_tot,
    #             'ZoneData': ZoneData, 'Congestion': Congestion, 'StorageData': StorageData,
    #             'UnitData': UnitData, 'FuelData': FuelData, 'WaterConsumptionData': WaterData}
    # else:
    return {'Cost_kwh': Cost_kwh, 'TotalLoad': TotalLoad, 'PeakLoad': PeakLoad, 'NetImports': NetImports,
            'Curtailment': Curtailment, 'MaxCurtailment': MaxCurtailemnt,
            'ShedLoad': LoadShedding, 'MaxShedLoad': MaxLoadShedding,
            'ShiftedLoad': ShiftedLoad_tot,
            'ZoneData': ZoneData, 'Congestion': Congestion, 'StorageData': StorageData,
            'UnitData': UnitData, 'FuelData': FuelData, 'WaterConsumptionData': WaterData}


def get_indicators_powerplant(inputs, results):
    """
    Function that analyses the dispa-set results at the power plant level
    Computes the number of startups, the capacity factor, etc

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :returns out:        Dataframe with the main power plants characteristics and the computed indicators
    """
    out = inputs['units'].loc[:, ['Nunits', 'PowerCapacity', 'Zone', 'Zone_th', 'Technology', 'Fuel']]

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
            out.loc[u, 'CF'] = results['OutputPower'][u].mean() / (out.loc[u, 'PowerCapacity'] * out.loc[u, 'Nunits'])
            out.loc[u, 'Generation'] = results['OutputPower'][u].sum()
        if u in results['OutputHeat']:
            out.loc[u, 'HeatGeneration'] = results['OutputHeat'][u].sum()
    # if inputs['param_df']['HeatDemand'].empty:
    #     # out.drop(['Zone_th'], axis=1, inplace=True)
    #     out.loc[:, 'Zone_th'] = ''
    return out


def CostExPost(inputs, results):
    """
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
             +Config("ValueOfLostLoad","val")*(sum(n,LL_MaxPower(n,i)+LL_MinPower(n,i)))
             +0.8*Config("ValueOfLostLoad","val")*(sum(n,LL_2U(n,i)+LL_2D(n,i)+LL_3U(n,i)))
             +0.7*Config("ValueOfLostLoad","val")*sum(u,LL_RampUp(u,i)+LL_RampDown(u,i))
             +Config("CostOfSpillage","val")*sum(s,spillage(s,i));


    :returns: tuple with the cost components and their cumulative sums in two dataframes.
    """
    import datetime

    dfin = inputs['param_df']
    timeindex = results['OutputPower'].index

    costs = pd.DataFrame(index=timeindex)

    # %% Fixed Costs:
    costs['FixedCosts'] = 0
    for u in results['OutputCommitted']:
        if u in dfin['CostFixed'].index:
            costs['FixedCosts'] = + dfin['CostFixed'].loc[u, 'CostFixed'] * results['OutputCommitted'][u]

    # %% Ramping and startup costs:
    indexinitial = timeindex[0] - datetime.timedelta(hours=1)
    powerlong = results['OutputPower'].copy()
    powerlong.loc[indexinitial, :] = 0
    powerlong.sort_index(inplace=True)
    committedlong = results['OutputCommitted'].copy()
    for u in powerlong:
        if u in dfin['PowerInitial'].index:
            powerlong.loc[indexinitial, u] = dfin['PowerInitial'].loc[u, 'PowerInitial']
            committedlong.loc[indexinitial, u] = dfin['PowerInitial'].loc[u, 'PowerInitial'] > 0
    committedlong.sort_index(inplace=True)

    powerlong_shifted = powerlong.copy()
    powerlong_shifted.iloc[1:, :] = powerlong.iloc[:-1, :].values
    committedlong_shifted = committedlong.copy()
    committedlong_shifted.iloc[1:, :] = committedlong.iloc[:-1, :].values

    ramping = powerlong - powerlong_shifted
   
    startups = committedlong.astype(int) - committedlong_shifted.astype(int)
    ramping.drop([ramping.index[0]], inplace=True)
    startups.drop([startups.index[0]], inplace=True)

    CostStartUp = pd.DataFrame(index=startups.index, columns=startups.columns)
    for u in CostStartUp:
        if u in dfin['CostStartUp'].index:
            CostStartUp[u] = startups[startups > 0][u].fillna(0) * dfin['CostStartUp'].loc[u, 'CostStartUp']
        else:
            print('Unit ' + u + ' not found in input table CostStartUp!')

    CostShutDown = pd.DataFrame(index=startups.index, columns=startups.columns)
    for u in CostShutDown:
        if u in dfin['CostShutDown'].index:
            CostShutDown[u] = startups[startups < 0][u].fillna(0) * dfin['CostShutDown'].loc[u, 'CostShutDown']
        else:
            print('Unit ' + u + ' not found in input table CostShutDown!')

    CostRampUp = pd.DataFrame(index=ramping.index, columns=ramping.columns)
    for u in CostRampUp:
        if u in dfin['CostRampUp'].index:
            CostRampUp[u] = ramping[ramping > 0][u].fillna(0) * dfin['CostRampUp'].loc[u, 'CostRampUp']
        else:
            print('Unit ' + u + ' not found in input table CostRampUp!')

    CostRampDown = pd.DataFrame(index=ramping.index, columns=ramping.columns)
    for u in CostRampDown:
        if u in dfin['CostRampDown'].index:
            CostRampDown[u] = ramping[ramping < 0][u].fillna(0) * dfin['CostRampDown'].loc[u, 'CostRampDown']
        else:
            print('Unit ' + u + ' not found in input table CostRampDown!')

    costs['CostStartUp'] = CostStartUp.sum(axis=1).fillna(0)
    costs['CostShutDown'] = CostShutDown.sum(axis=1).fillna(0)
    costs['CostRampUp'] = CostRampUp.sum(axis=1).fillna(0)
    costs['CostRampDown'] = CostRampDown.sum(axis=1).fillna(0)

    # %% Variable cost:
    costs['CostVariable'] = (results['OutputPower'] * dfin['CostVariable']).fillna(0).sum(axis=1)

    # %% Transmission cost:
    costs['CostTransmission'] = (results['OutputFlow'] * dfin['PriceTransmission']).fillna(0).sum(axis=1)

    # %% Shedding cost:
    costs['CostLoadShedding'] = (results['OutputShedLoad'] * dfin['CostLoadShedding']).fillna(0).sum(axis=1)

    # %% Heating costs:
    # costs['CostHeatSlack'] = (results['OutputHeatSlack'] * dfin['CostHeatSlack']).fillna(0).sum(axis=1)
    CostHeat = pd.DataFrame(index=results['OutputHeat'].index, columns=results['OutputHeat'].columns)
    for u in CostHeat:
        if u in dfin['CHPPowerLossFactor'].index:
            CostHeat[u] = dfin['CostVariable'][u].fillna(0) * results['OutputHeat'][u].fillna(0) * \
                          dfin['CHPPowerLossFactor'].loc[u, 'CHPPowerLossFactor']
        else:
            CostHeat[u] = dfin['CostVariable'][u].fillna(0) * results['OutputHeat'][u].fillna(0)
    costs['CostHeat'] = CostHeat.sum(axis=1).fillna(0)

    # # %% Cost H2:
    # CostH2 = pd.DataFrame(index=costs.index, columns=dfin['CostH2Slack'].columns)
    # for u in dfin['CostH2Slack'].columns:
    #     CostH2[u] = dfin['CostH2Slack'][u].fillna(0) * results['OutputStorageSlack'][u].fillna(0)
    # costs['CostH2Slack'] = CostH2.fillna(0).sum(axis=1)

    # %% Lost loads:
    # NB: the value of lost load is currently hard coded. This will have to be updated
    # Locate prices for LL
    # TODO:
    costs['LostLoad'] = 80e3 * (results['LostLoad_2D'].reindex(timeindex).sum(axis=1).fillna(0) +
                                results['LostLoad_2U'].reindex(timeindex).sum(axis=1).fillna(0) +
                                results['LostLoad_3U'].reindex(timeindex).sum(axis=1).fillna(0)) + \
                        100e3 * (results['LostLoad_MaxPower'].reindex(timeindex).sum(axis=1).fillna(0) +
                                 results['LostLoad_MinPower'].reindex(timeindex).sum(axis=1).fillna(0)) + \
                        70e3 * (results['LostLoad_RampDown'].reindex(timeindex).sum(axis=1).fillna(0) +
                                results['LostLoad_RampUp'].reindex(timeindex).sum(axis=1).fillna(0))

    # %% Spillage:
    costs['Spillage'] = 1 * results['OutputSpillage'].sum(axis=1).fillna(0)

    # %% Plotting
    # Drop na columns:
    costs.dropna(axis=1, how='all', inplace=True)
    # Delete all-zero columns:
    # costs = costs.loc[:, (costs != 0).any(axis=0)]

    sumcost = costs.cumsum(axis=1)
    sumcost['OutputSystemCost'] = results['OutputSystemCost']

    sumcost.plot(title='Cumulative sum of the cost components')

    # %% Warning if significant error:
    diff = (costs.sum(axis=1) - results['OutputSystemCost']).abs()
    if diff.max() > 0.01 * results['OutputSystemCost'].max():
        logging.critical(
            'There are significant differences between the cost computed ex post and and the cost provided by the optimization results!')
    return costs, sumcost


def get_units_operation_cost(inputs, results):
    """
    Function that computes the operation cost for each power unit at each instant of time from the DispaSET results
    Operation cost includes: CostFixed + CostStartUp + CostShutDown + CostRampUp + CostRampDown + CostVariable

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :returns out:       Dataframe with the the power units in columns and the operation cost at each instant in rows

    Main Author: @AbdullahAlawad
    """
    datain = inputs['param_df']


    # make sure that we have the same columns in the committed table and in the power table:
    for u in results['OutputCommitted']:
        if u not in results['OutputPower']:
            results['OutputPower'][u] = 0

    # DataFrame with startup times for each unit (1 for startup)
    StartUps = results['OutputCommitted'].copy()
    for u in StartUps:
        values = StartUps.loc[:, u].values
        diff = -(np.roll(values, 1) - values)
        diff[diff <= 0] = 0
        StartUps[u] = diff

    # DataFrame with shutdown times for each unit (1 for shutdown)
    ShutDowns = results['OutputCommitted'].copy()
    for u in ShutDowns:
        values = ShutDowns.loc[:, u].values
        diff = (np.roll(values, 1) - values)
        diff[diff <= 0] = 0
        ShutDowns[u] = diff

    # DataFrame with ramping up levels for each unit at each instant (0 for ramping-down & leveling out)
    RampUps = results['OutputPower'].copy()
    for u in RampUps:
        values = RampUps.loc[:, u].values
        diff = -(np.roll(values, 1) - values)
        diff[diff <= 0] = 0
        RampUps[u] = diff

    # DataFrame with ramping down levels for each unit at each instant (0 for ramping-up & leveling out)
    RampDowns = results['OutputPower'].copy()
    for u in RampDowns:
        values = RampDowns.loc[:, u].values
        diff = (np.roll(values, 1) - values)
        diff[diff <= 0] = 0
        RampDowns[u] = diff

    FiexedCost = results['OutputCommitted'].copy()
    StartUpCost = results['OutputCommitted'].copy()
    ShutDownCost = results['OutputCommitted'].copy()
    RampUpCost = results['OutputCommitted'].copy()
    RampDownCost = results['OutputCommitted'].copy()
    VariableCost = results['OutputCommitted'].copy()
    PowerUnitOperationCost = results['OutputCommitted'].copy()
    ConsumptionCost = results['OutputPowerConsumption'].copy()

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
        FiexedCost.loc[:, [u]] = np.array(results['OutputCommitted'].loc[:, [u]]) * \
                                 inputs['parameters']['CostFixed']['val'][unit_indexNo]
        StartUpCost.loc[:, [u]] = np.array(StartUps.loc[:, [u]]) * inputs['parameters']['CostStartUp']['val'][
            unit_indexNo]
        ShutDownCost.loc[:, [u]] = np.array(ShutDowns.loc[:, [u]]) * inputs['parameters']['CostShutDown']['val'][
            unit_indexNo]
        RampUpCost.loc[:, [u]] = np.array(RampUps.loc[:, [u]]) * inputs['parameters']['CostRampUp']['val'][unit_indexNo]
        RampDownCost.loc[:, [u]] = np.array(RampDowns.loc[:, [u]]) * inputs['parameters']['CostRampDown']['val'][
            unit_indexNo]
        VariableCost.loc[:, [u]] = np.array(datain['CostVariable'].loc[:, [u]]) * np.array(
            results['OutputPower'][u]).reshape(-1, 1)

    PowerConsumers = results['OutputPowerConsumption'].columns
    # Shadowprice = results['ShadowPrice'].head(inputs['config']['LookAhead']*-24+1)#reindexing
    Shadowprice = results['ShadowPrice']
    for u in PowerConsumers:#at this moment only variable costs
        z = inputs['units'].at[u,'Zone']
        ConsumptionCost.loc[:,[u]] = np.array(Shadowprice.loc[:,[z]])*np.array(results['OutputPowerConsumption'][u]).reshape(-1,1)
      
    PowerUnitOperationCost = FiexedCost + StartUpCost + ShutDownCost + RampUpCost + RampDownCost + VariableCost
    UnitOperationCost = pd.concat([PowerUnitOperationCost, ConsumptionCost], axis=1)
    
    return UnitOperationCost


def get_EFOH(inputs, results):
    """
    Function that computes the "Equivalent Full Load Operating Hours" of the Elyzers
    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :returns EFOH:      Dataframe with the EFOH of each country
    """
    EFOH = pd.DataFrame(index=inputs['sets']['p2h2'], columns=['EFOH'])
    for i in EFOH.index:
        for count, j in enumerate(inputs['sets']['au']):
            if i == j:
                Cap = inputs['parameters']['StorageChargingCapacity']['val'][count]
                StoInput = results['OutputStorageInput'].loc[:, j].sum()
                EFOH.loc[i, 'EFOH'] = StoInput / Cap


def get_power_flow_tracing(inputs, results, idx=None, type=None):
    """
    Up-stream and down-stream (not yet implemented) looking algorithms using Bialek's power flow tracing method
    (Bialek at al. DOI: 10.1049/ip-gtd:19960461)

    :param inputs:      Dispa-SET inputs
    :param results:     Dispa-SET results
    :param idx:         pandas datetime index (can be both a single and multiple time intervals)
    :return:            NxN matrix of zones (Row -> Excess generation zones, Col -> Lack of generation)
    """
    # Extract input data
    flows = results['OutputFlow']
    zones = inputs['sets']['n']
    lines = inputs['sets']['l']
    Demand = inputs['param_df']['Demand']['DA'].copy()
    ShedLoad = results['OutputShedLoad'].copy()

    # Adjust idx if not specified
    if idx is None:
        idx = pd.date_range(Demand.index[0], periods=24, freq='H')
        logging.warning('Date range not specified, Power Flow will be traced only for the first day of the simulation')

    flows = flows.reindex(columns=lines).fillna(0)
    # Predefine dataframes
    NetExports = pd.DataFrame(index=flows.index)
    Generation = pd.DataFrame(index=flows.index)
    NetExportsZ = [0]
    for z in zones:
        Generation = pd.concat([Generation, pd.DataFrame(filter_by_zone(results['OutputPower'], inputs, z).sum(axis=1)).set_axis([z], axis=1)], axis=1)
        NetExportsZ= [0]
        for key in flows:
            if key[:len(z)] == z:
                NetExportsZ += flows.loc[:, key].copy()
        if any(NetExportsZ):
            NetExports = pd.concat([NetExports, pd.DataFrame(NetExportsZ).set_axis([z], axis=1)], axis=1)
        else:
            NetExports = pd.concat([NetExports, pd.DataFrame(np.zeros((len(flows), 1))).set_axis([z], axis=1).set_axis(flows.index, axis=0)], axis=1)
    Demand = Demand.loc[:, (Demand != 0).any(axis=0)]

    if ShedLoad.empty:
        P = (Demand.loc[idx].sum() + NetExports.loc[idx].sum()).dropna()
    # TODO: Check if lost load and curtailment are messing up the P and resulting in NaN
    else:
        ShedLoad = ShedLoad.reindex(columns=zones).fillna(0)
        P = (Demand.loc[idx].sum() + NetExports.loc[idx].sum() - ShedLoad.loc[idx].sum()).dropna()

    D = pd.DataFrame(index=zones, columns=zones).fillna(0).astype(float)
    B = pd.DataFrame(index=zones, columns=zones).fillna(0).astype(float)
    np.fill_diagonal(D.values, P.values)
    for l in inputs['sets']['l']:
        if l in flows.columns:
            [from_node, to_node] = l.split(' -> ')
            if (from_node.strip() in zones) and (to_node.strip() in zones):
                D.loc[from_node, to_node] = -flows.loc[idx, l].sum()
                B.loc[from_node, to_node] = flows.loc[idx, l].sum()
        else:
            logging.info('Line ' + l + ' was probably not used. Skipping!')
    A = D.loc[Demand.columns][Demand.columns].T / P.values
    # A.fillna(0, inplace=True)
    A_T = pd.DataFrame(np.linalg.inv(A.values), A.columns, A.index)

    Share = Demand.loc[idx].sum() / P
    gen = A_T * Generation.loc[idx].sum()
    Trace = pd.DataFrame(index=Demand.columns, columns=Demand.columns)
    for z in Demand.columns:
        tmp = gen[Demand.columns].loc[z, :] * Share.loc[z]
        Trace[z] = tmp
    Trace = Trace.T
    Trace_prct = Trace.div(Trace.sum(axis=1), axis=0)
    # TODO: Implement capability for RoW flows (not sure if curent algorithm is taking it into the account properly)

    return Trace, Trace_prct


def get_from_to_flows(inputs, flows, zones, idx=None):
    """
    Helper function for braking down flows into networkx readable format

    :param inputs:      Dispa-SET inputs
    :param flows:       Flows from Dispa-SET results
    :param zones:       List of selected zones
    :param idx:         datetime index (can be a single hour or a range of hours)
    :return:
    """
    Flows = pd.DataFrame(columns=['From', 'To', 'Flow'])
    i = 0
    for l in inputs['sets']['l']:
        if l in flows.columns:
            [from_node, to_node] = l.split(' -> ')
            if (from_node.strip() in zones) and (to_node.strip() in zones):
                Flows.loc[i, 'From'] = from_node.strip()
                Flows.loc[i, 'To'] = to_node.strip()
                if isinstance(idx, pd.DatetimeIndex):
                    Flows.loc[i, 'Flow'] = flows.loc[idx, l].sum()
                else:
                    Flows.loc[i, 'Flow'] = flows.loc[0, l].sum()
                i = i + 1
    return Flows


def get_net_positions(inputs, results, zones, idx):
    """
    Helper function for calculating net positions in individual zones

    :param inputs:      Dispa-SET inputs
    :param results:     Dispa-SET results
    :param zones:       List of selected zones
    :param idx:         datetime index (can be a single hour or a range of hours)
    :return:            NetImports and Net position
    """
    NetImports = pd.DataFrame(columns=zones)
    for z in zones:
        NetImports.loc[0, z] = get_imports(results['OutputFlow'].loc[idx], z)

    Demand = inputs['param_df']['Demand']['DA'].copy()
    NetExports = -NetImports.copy()
    NetExports[NetExports < 0] = 0
    NetPosition = Demand.loc[idx].sum() + NetExports.iloc[0]
    return NetImports, NetPosition

#%%
def shadowprices(results, zone):
    """
    this function retrieves the schadowprices of DA, 2U, 2D, 3U for 1 zone
    """
    shadowprices = pd.DataFrame(0,index = results['OutputPower'].index, columns = ['DA','2U','2D','3U'])
    
    if  zone in results['ShadowPrice'].columns:
        shadowprices['DA'] = results['ShadowPrice'][zone]
    if zone in results['ShadowPrice_2U'].columns:
        shadowprices['2U'] = results['ShadowPrice_2U'][zone]
    if zone in results['ShadowPrice_2D'].columns:
        shadowprices['2D'] = results['ShadowPrice_2D'][zone]        
    if zone in results['ShadowPrice_3U'].columns:
        shadowprices['3U'] = results['ShadowPrice_3U'][zone]

    shadowprices.fillna(0,inplace=True)
    return shadowprices
#%%
def Cashflows(inputs,results,unit):
    """
    This function calculates the different cashflows (DA,2U,2D,3U,Heat,costs) for one specific unit
    returns: hourly cashflow
    """
    zone = inputs['units'].at[unit,'Zone']
    cashflows = pd.DataFrame(index = inputs['config']['idx'])
    
    TMP = shadowprices(results, zone)

    
    if TMP['DA'].max() > 10000:
        for i in TMP.index:
            if TMP.loc[i,'DA']>10000:
                TMP.loc[i,'DA']=TMP.at[i-1,'DA']
        
    if TMP['2U'].max() > 10000:
        for i in TMP.index:
            if TMP.loc[i,'2U']>10000:
                TMP.loc[i,'2U']=TMP.at[i-1,'2U']
        
    if TMP['3U'].max() > 10000:
        for i in TMP.index:
            if TMP.loc[i,'3U']>10000:
                TMP.loc[i,'3U']=TMP.at[i-1,'3U']
        
    if TMP['2D'].max() > 10000:
        for i in TMP.index:
            if TMP.loc[i,'2D']>10000:
                TMP.loc[i,'2D']=TMP.at[i-1,'2D']

    
    #positive cashflows
    if  unit in results['OutputPower'].columns and zone in results['ShadowPrice'].columns:
        cashflows['DA'] = results['OutputPower'][unit]*TMP['DA']
    else:
        cashflows['DA'] = 0
    if unit in results['OutputReserve_2U'].columns and zone in results['ShadowPrice_2U'].columns:
        cashflows['2U'] = results['OutputReserve_2U'][unit]*TMP['2U']
    else:
        cashflows['2U'] = 0
    if unit in results['OutputReserve_2D'].columns and zone in results['ShadowPrice_2D'].columns:
        cashflows['2D'] = results['OutputReserve_2D'][unit]*TMP['2D']
    else:
        cashflows['2D'] = 0
    if unit in results['OutputReserve_3U'].columns and zone in results['ShadowPrice_3U'].columns:
        cashflows['3U'] = results['OutputReserve_3U'][unit]*TMP['3U']
    else:
        cashflows['3U'] = 0
    if unit in results['OutputHeat'].columns and unit in results['HeatShadowPrice'].columns:
        cashflows['heat'] = results['OutputHeat'][unit]*results['HeatShadowPrice'][unit]
    else:
        cashflows['heat'] = 0

    #negative cashflow
    units_operation_cost = get_units_operation_cost(inputs, results)
    
    cashflows['costs'] = -units_operation_cost[unit]

    cashflows.fillna(0,inplace = True)
    return cashflows
#%%
def reserve_availability_demand(inputs,results):
    """
    this function evaluates the reserve demand and availability of all units
    returns: hourly_availability : hourly availability over the corresponding hourly reserve demand in [%]
    returns: availability        : the mean of hourly availability,
                                   total hourly availability over total corresponding hourly reserve demand in [%]
    returns: reserve_demand        mean demand over peak load
    """
    
    hourly_availability = {}
    hourly_availability['2U'] = pd.DataFrame()
    hourly_availability['3U'] = pd.DataFrame()
    hourly_availability['Down'] = pd.DataFrame()
    hourly_availability['Up'] = pd.DataFrame()
    

    total_up_reserves = pd.concat([results['OutputReserve_2U'],results['OutputReserve_3U']],axis=1)
    total_up_reserves = total_up_reserves.groupby(level=0, axis=1).sum()

    availability = {}
    availability['2U'] = pd.DataFrame(0.0,index = results['OutputReserve_2U'].columns, columns = ['mean','total'])
    availability['3U'] = pd.DataFrame(0.0,index = results['OutputReserve_3U'].columns, columns = ['mean','total'])
    availability['Down'] = pd.DataFrame(0.0,index = results['OutputReserve_2D'].columns, columns = ['mean','total'])
    availability['Up'] = pd.DataFrame(0.0,index = total_up_reserves.columns, columns = ['mean','total'])
    
    
    for unit in results['OutputReserve_2U'].columns:
        zone = inputs['units'].at[unit,'Zone']
        hourly_availability['2U'][unit] = results['OutputReserve_2U'][unit]/(inputs['param_df']['Demand']['2U',zone]/2)*100
        availability['2U'].at[unit,'mean'] = hourly_availability['2U'][unit].mean()
        availability['2U'].at[unit,'total'] = results['OutputReserve_2U'][unit].sum()/(inputs['param_df']['Demand']['2U',zone].sum()/2)*100
        
    for unit in results['OutputReserve_3U'].columns:
        zone = inputs['units'].at[unit,'Zone']
        hourly_availability['3U'][unit] = results['OutputReserve_3U'][unit]/(inputs['param_df']['Demand']['2U',zone]/2)*100
        availability['3U'].at[unit,'mean'] = hourly_availability['3U'][unit].mean()
        availability['3U'].at[unit,'total'] = results['OutputReserve_3U'][unit].sum()/(inputs['param_df']['Demand']['2U',zone].sum()/2)*100

    for unit in results['OutputReserve_2D'].columns:
        zone = inputs['units'].at[unit,'Zone']
        hourly_availability['Down'][unit] = results['OutputReserve_2D'][unit]/inputs['param_df']['Demand']['2D',zone]*100
        availability['Down'].at[unit,'mean'] = hourly_availability['Down'][unit].mean()
        availability['Down'].at[unit,'total'] = results['OutputReserve_2D'][unit].sum()/inputs['param_df']['Demand']['2D',zone].sum()*100

    for unit in total_up_reserves.columns:
        zone = inputs['units'].at[unit,'Zone']
        hourly_availability['Up'][unit] = total_up_reserves[unit]/inputs['param_df']['Demand']['2U',zone]*100
        availability['Up'].at[unit,'mean'] = hourly_availability['Up'][unit].mean()
        availability['Up'].at[unit,'total'] = total_up_reserves[unit].sum()/inputs['param_df']['Demand']['2U',zone].sum()*100
    
    reserve_demand = pd.DataFrame(0.0,index = inputs['config']['zones'], columns = ['upwards','downwards'])
    for zone in inputs['config']['zones']:
        reserve_demand.at[zone,'upwards'] = inputs['param_df']['Demand']['2U',zone].mean()/inputs['param_df']['Demand']['DA',zone].max()
        reserve_demand.at[zone,'downwards'] = inputs['param_df']['Demand']['2D',zone].mean()/inputs['param_df']['Demand']['DA',zone].max()
        
    return hourly_availability, availability, reserve_demand

#%%
def emissions(inputs,results):
    """
    function calculates emissions for each zone
    returns: emissions : dataframe with summed up emissions for each zone
    """
    emissions = pd.DataFrame(0,index=inputs['config']['zones'], columns = ['emissions'])

    unit_emissions = results['OutputPower'].sum() * inputs['param_df']['EmissionRate'] # MWh * tCO2 / MWh = tCO2
    unit_emissions.fillna(0,inplace=True)
    
    for zone in inputs['config']['zones']:
        emissions.at[zone,'emissions'] = filter_by_zone(unit_emissions, inputs, zone).sum(axis=1)
    
    return emissions
    
#%% load shedding
def load_shedding(inputs,results):
    loadshedding = pd.DataFrame(0,index = ['max','sum','amount'], columns = inputs['config']['zones'])
    for z in inputs['config']['zones']:
        if z in results['OutputShedLoad']:
            loadshedding.loc['max',z] = results['OutputShedLoad'][z].max()
            loadshedding.loc['sum',z] = results['OutputShedLoad'][z].sum()
            loadshedding.loc['amount',z] = (results['OutputShedLoad'][z]!=0).sum()

    return loadshedding

#%% curtailment
def curtailment(inputs,results):
    curtailment = pd.DataFrame(0,index = ['max','sum','amount'], columns = inputs['config']['zones'])
    for z in inputs['config']['zones']:
        if z in results['OutputCurtailedPower']:
            curtailment.loc['max',z] = results['OutputCurtailedPower'][z].max()
            curtailment.loc['sum',z] = results['OutputCurtailedPower'][z].sum()
            curtailment.loc['amount',z] = (results['OutputCurtailedPower'][z]!=0).sum()

    return curtailment

# %%
def get_frequency_security_constraints(inputs, results):
    """
    Reads the DispaSET results and provides frequency security constraints 

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    """
    
    OutputCommitted = results['OutputCommitted']
    OutputPower = results['OutputPower'] 
    units = inputs['units']
        
    # Filter "units" to keep rows with 'WAT', 'GAS', 'BIO', 'OIL' in the "Fuel" column
    convunits = units[units['Fuel'].isin(['WAT', 'GAS', 'BIO', 'OIL'])]
        
    # Get the list of generator names from the filtered "units"
    generator_list = units.Unit.tolist()
    convgenerator_list = convunits.Unit.tolist()
            
    # Filter the dataframe to keep only the columns that exist in the dataframe
    generator_list = [col for col in generator_list if col in OutputCommitted.columns]
        
    # Filter the dataframe to keep only the columns that exist in the dataframe
    convgenerator_list = [col for col in convgenerator_list if col in OutputCommitted.columns]
        
    # Keep only the columns in "OutputCommitted" that are in "convgenerator_list"
    OutputSyncCommitted = OutputCommitted[convgenerator_list]
        
    # Create an empty DataFrame for the result
    GeneratorContingency = pd.DataFrame(columns=['TIMESTAMP', 'GeneratorName', 'PowerCapacity'])
    
        
    # Iterate through the rows of OutputSyncCommitted
    for timestamp, row in OutputSyncCommitted.iterrows():
        # Find generators with commitment value of 1
        committed_generators = [col for col in row.index if row[col] == 1]
        
        # If there are committed generators, find the one with the highest PowerCapacity
        if committed_generators:
            max_committed_generator = max(committed_generators, key=lambda x: convunits.loc[x, 'PowerCapacity'])
            max_committed_power = convunits.loc[max_committed_generator, 'PowerCapacity']
            
            # Append the result to GeneratorContingency
            GeneratorContingency = GeneratorContingency.append({'TIMESTAMP': timestamp,
                                                                'GeneratorName': max_committed_generator,
                                                                'PowerCapacity': max_committed_power},
                                                                 ignore_index=True)
    
    
    # Inicializar un DataFrame vacío para almacenar los resultados
    total_power_per_time = pd.DataFrame(columns=['TIMESTAMP', 'TotalPower'])
    
    # Iterar sobre las filas de 'OutputCommitted'
    for timestamp, row in OutputCommitted.iterrows():
        # Filtrar 'relevant_units' para obtener solo las columnas de generadores encendidos (valor igual a 1)
        committed_generators = [col for col in row.index if row[col] == 1]
    
        # Crear una serie booleana indicando si cada generador está encendido
        is_committed = units['Unit'].isin(committed_generators)
        
        # Filtrar el DataFrame original para incluir solo generadores encendidos
        committed_units = units[is_committed]
        
        # Calcular la suma de PowerCapacity para generadores encendidos
        total_power = committed_units['PowerCapacity'].sum()
        total_power_per_time = total_power_per_time.append({'TIMESTAMP': timestamp, 'TotalPower': total_power}, ignore_index=True)
    
    
    # Set the 'Hour' column as the index
    GeneratorContingency.set_index('TIMESTAMP', inplace=True)
    total_power_per_time.set_index('TIMESTAMP', inplace=True)
        
    # Print the resulting DataFrame
    print(GeneratorContingency)
        
    # Merge DataFrames across columns
    dispatch = pd.concat([GeneratorContingency, total_power_per_time], axis=1)
        
    # rename columns
    dispatch = dispatch.rename(columns={'PowerCapacity': 'pcontingency','TotalPower': 'pcommitted'})
    dispatch['ploss'] = (dispatch['pcontingency']/(dispatch['pcommitted']-dispatch['pcontingency']))*-1
    # dispatch['ploss'] = (dispatch['pcontingency'])*-1
    
    dispatch = dispatch.reset_index()
    
    # Find the ploss max for each GeneratorName
    max_ploss = dispatch.groupby('GeneratorName')['ploss'].max()



# %% power swing solutions
# Create an empty dictionary to store results
    freq_security = {}
    
    # Crear un archivo Excel para almacenar los resultados
    with pd.ExcelWriter('power_swing_results.xlsx', engine='xlsxwriter') as writer:
        contingency_count = 1  # Inicializar el contador de hojas
    
        for index, row in dispatch.iterrows():
            # Perform the operations on each row
    
            def power_swing(y, t, k, H):
                deltap, f = y
    
                if t < 1:
                    ploss = 0
                    pprim = 0
                elif t < time_preparation:
                    ploss = row['ploss']
                    deltap = ploss
                    pprim = 0
                else:
                    ploss = row['ploss'] 
                    pcommitted = row['pcommitted']
                    pcontingency = row['pcontingency']
                    pprim = (-k * f) / (pcommitted-pcontingency) * (t-time_preparation)
                    deltap = ploss + pprim
    
                delf = deltap*1 / (2 * H )
                return [deltap, delf]
    
            # Initial conditions
            f0 = 0
            y0 = [0, f0]
    
            # Time vector
            t = np.arange(0, 61, 1)  # Time steps every 1 second from 0 to 60 seconds
    
            # Create an empty DataFrame to store results
            results_df = pd.DataFrame({'Time': t})
    
            # Function to calculate frequency at each time step
            def calc_frequency_at_time(k, H):
                sol = odeint(power_swing, y0, t, args=(k, H))
                f_values = sol[:, 1]  # Frequency values at each time step
                return f_values
    
            # Constraints
            def constraints_satisfied(k, H):
                f_values = calc_frequency_at_time(k, H)
                der_f_values = np.gradient(f_values, t)
                min_der_f = np.min(der_f_values)
                min_f = np.min(f_values)
                f_at_60s = f_values[-1]
    
                return min_der_f >= -0.002 and min_f >= -0.016 and f_at_60s >= -0.01
    
            # Search for the first solution by iterating through values of H and k
            lower_bound_H = 1
            upper_bound_H = 1000  # Adjust the upper bound as needed
            # lower_bound_k = 0
            # upper_bound_k = 4000  # Adjust the upper bound as needed
            step = 1
    
            found_solution = False
            # min_k = None                                  # k variable
            
            # Set the constant and time preparation values
            k = 1700  # Valor constante de k
            time_preparation = 5
            min_H = None
    
            for H in range(lower_bound_H, upper_bound_H + 1, step):
                # for k in range(lower_bound_k, upper_bound_k + 1, step):      # k variable
                    #if constraints_satisfied(k, H):                        # k variable
                    if constraints_satisfied(k, H):
                        found_solution = True
                        #min_k = k          # k variable
                        min_H = H
                        break
                # if found_solution      # k variable
                #     break              # k variable
    
            if found_solution:
                print(f"The preparation time for the Primary Reserve is set at: {time_preparation}")
                print(f"The proportional gain constant k of the speed regulators is set at: {k}")
                #print(f"Minimum k that satisfies constraints: {min/_k}")
                print(f"Minimum H that satisfies constraints: {min_H}")
                # Update the DataFrame with the new values
                dispatch.at[index, 'k'] = k  # Store in 'k'
                dispatch.at[index, 'min_H'] = min_H  # Store in 'min_H'
    
                # Solve the ODE system with the optimal k and H
                sol = odeint(power_swing, y0, t, args=(k, min_H))
    
                # Extract the results
                deltap = sol[:, 0]
                f = sol[:, 1]
                der_f = np.gradient(f, t)
                pprim = np.where(t < 5, 0, (-k * f) / row['pcommitted'] * (t-5))
                ploss = np.where(t < 5, row['ploss'], row['ploss'])
                pcommitted = row['pcommitted']
                pcontingency = row['pcontingency']
    
                # Store results in the DataFrame
                results_df['TIMESTAMP'] = row['TIMESTAMP']
                results_df['Frequency'] = f
                results_df['RoCoF'] = der_f
                results_df['DeltaP'] = deltap
                results_df['pcontingency'] = pcontingency
                results_df['PrimaryResponsePU'] = pprim
                results_df['PrimaryResponseMW'] = pprim*(pcommitted-pcontingency)
                results_df['pcommitted'] = pcommitted
                results_df['ploss'] = ploss
    
                # Add the values of k and H to the DataFrame
                results_df['K'] = k
                results_df['min_H'] = min_H
    
                # Store the results_df in the dictionary
                freq_security[f"Contingency{contingency_count}"] = results_df    
    
                # Guardar el DataFrame en una pestaña numerada
                sheet_name = f'Contingency{contingency_count}'
                results_df.to_excel(writer, sheet_name=sheet_name, index=False)
                
                contingency_count += 1  # Incrementar el contador de hojas
                
                # Plotting code remains the same
                fig, ax1 = plt.subplots(figsize=(12, 8))
            
                # Frequency Values Plot (Left Y-axis)
                ax1.plot(t, f, color='tab:blue', linestyle='-', label='$∆f_{COI}$')
                ax1.plot(t, der_f, color='tab:orange', linestyle='-', label='$d(∆f_{COI})$')
                ax1.set_xlabel('Time [s]', fontsize=16)
                ax1.set_ylabel('Frequency [p.u.]', color='tab:blue', fontsize=16)
                ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=14)
                ax1.legend(loc='upper left', fontsize=14)
            
                # Power Values Plot (Right Y-axis)
                ax2 = ax1.twinx()
                ax2.plot(t, deltap, color='tab:green', linestyle='--', label='∆Power')
                ax2.plot(t, ploss, color='tab:red', linestyle='--', label='PowerLoss')
                ax2.plot(t, pprim, color='tab:purple', linestyle='--', label='PrimaryReserve')
                ax2.set_ylabel('Power [p.u.]', color='tab:red', fontsize=16)
                ax2.tick_params(axis='y', labelcolor='tab:red', labelsize=14)
                ax2.legend(loc='upper right', fontsize=14)
            
                # Overall Plot Settings
                plt.title('Power Swing', fontsize=18)
                plt.xlim(0, 60)
                plt.tight_layout()
                plt.show()
                
            else:
                print("No solution found within the specified range of H and k that satisfies constraints.")
    
    # Crear un DataFrame vacío para almacenar los resultados
    summary = pd.DataFrame(columns=['TIMESTAMP', 'Contingency', 'TotalPowerCommitted', 'KineticEnergy', 'Inertia','PrimaryResponse'])
    
    # Para obtener el valor de Kinetic Energy
    for key, df in freq_security.items():
        # Calcula la nueva columna 'KineticEnergy'
        df['KineticEnergy'] = df['min_H'] * df['pcontingency'] - df['pcommitted']
    
    # Iterar sobre los DataFrames en el diccionario y tomar el primer valor de 'A' y 'B'
    for nombre_df, df in freq_security.items():
        primer_timestamp = df['TIMESTAMP'].iloc[0]
        primer_contingency = df['pcontingency'].iloc[0]  # Primer valor de la columna 'A'
        primer_pcommitted = df['pcommitted'].iloc[0]
        primer_KineticEnergy = df['KineticEnergy'].iloc[0]  # Primer valor de la columna 'KineticEnergy'
        primer_inertia = df['min_H'].iloc[0]  # Primer valor de la columna 'B'
        max_PrimaryResponseMW = df['PrimaryResponseMW'].max()  # Valor maximo de la columna 'PrimaryResponse'
        
        # Agregar los valores a 'Summary'
        summary = summary.append({'TIMESTAMP': primer_timestamp, 'Contingency': primer_contingency, 'TotalPowerCommitted': primer_pcommitted, 'KineticEnergy': primer_KineticEnergy, 'Inertia': primer_inertia, 'PrimaryResponse': max_PrimaryResponseMW}, ignore_index=True)
    
    dispatch.set_index('TIMESTAMP', inplace=True)
    summary.set_index('TIMESTAMP', inplace=True)
        
    return freq_security, summary, dispatch