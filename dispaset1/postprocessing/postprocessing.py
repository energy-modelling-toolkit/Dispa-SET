# -*- coding: utf-8 -*-
"""
Set of functions useful to analyse to DispaSET output data.

@author: Sylvain Quoilin, JRC
"""

from __future__ import division

import datetime as dt
import logging
import os
import pickle
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..misc.gdx_handler import gdx_to_list, gdx_to_dataframe, get_gams_path
from ..misc.str_handler import shrink_to_64, clean_strings

# 'chocolate', darkviolet, maroon, slateblue
COLORS = {'NUC': 'slateblue', 'LFG': 'aqua', 'DSL': 'grey', 'GAS': 'khaki', 'OIL': 'lightcoral',
          'MSW': 'dodgerblue', 'SUN': 'yellow', 'HFO': 'red', 'WIN': 'greenyellow',
          'FlowIn': 'green', 'FlowOut': 'green'}

def GAMSstatus(statustype,num):
    '''
    Function that returns the model status or the solve status from gams

    :param statustype: String with the type of status to retrieve ("solver" or "model")
    :param num:     Indicated termination condition (Integer)
    :returns:       String with the status
    '''
    if statustype=="model":
        msg =   {1: u'Optimal solution achieved',
                 2: u'Local optimal solution achieved',
                 3: u'Unbounded model found',
                 4: u'Infeasible model found',
                 5: u'Locally infeasible model found (in NLPs)',
                 6: u'Solver terminated early and model was infeasible',
                 7: u'Solver terminated early and model was feasible but not yet optimal',
                 8: u'Integer solution model found',
                 9: u'Solver terminated early with a non integer solution found (only in MIPs)',
                 10: u'No feasible integer solution could be found',
                 11: u'Licensing problem',
                 12: u'Error achieved, No cause known',
                 13: u'Error achieved, No solution attained',
                 14: u'No solution returned',
                 15: u'Feasible in a CNS models',
                 16: u'Locally feasible in a CNS models',
                 17: u'Singular in a CNS models',
                 18: u'Unbounded, no solution',
                 19: u'Infeasible, no solution'}
    elif statustype=="solver":
        msg =   {1: u'Normal termination',
                 2: u'Solver ran out of iterations (fix with iterlim)',
                 3: u'Solver exceeded time limit (fix with reslim)',
                 4: u'Solver quit with a problem (see LST file) found',
                 5: u'Solver quit with excessive nonlinear term evaluation errors (see LST file and fix with bounds or domlim)',
                 6: u'Solver terminated for unknown reason (see LST file)',
                 7: u'Solver terminated with preprocessor error (see LST file)',
                 8: u'User interrupt',
                 9: u'Solver terminated with some type of failure (see LST file)',
                 10: u'Solver terminated with some type of failure (see LST file)',
                 11: u'Solver terminated with some type of failure (see LST file)',
                 12: u'Solver terminated with some type of failure (see LST file)',
                 13: u'Solver terminated with some type of failure (see LST file)'}
    else:
        sys.exit('Incorrect GAMS status type')
    return str(msg[num])


def unit_location(Inputs, shrink=True):
    """ 
    Function that associates its location to each unit from the Dispaset inputs

    :param inputs:      DispaSet inputs (version 2.1.1)
    :returns loc:        Dictionary with the location of each unit
    """
    loc = {}
    if isinstance(Inputs, list):
        u = Inputs[0]['u']
        n = Inputs[0]['n']
        location = Inputs[1]['Location']['val']
    elif isinstance(Inputs, dict):
        u = Inputs['sets']['u']
        n = Inputs['sets']['n']
        location = Inputs['parameters']['Location']['val']
    else:
        logging.error('Inputs variable no valid')
        sys.exit(1)
    if shrink:
        uu = shrink_to_64(u)
    else:
        uu = u
    for i in range(len(u)):
        j = np.where(location[i, :] > 0)
        loc[uu[i]] = n[j[0][0]]
    return loc


def unit_fuel(Inputs, shrink=True):
    """ 
    Function that associates its fuel to each unit from the Dispaset inputs

    :param inputs:      DispaSet inputs (version 2.1.1)
    :param shrink:      If True, the unit name is reduced to 64, as in GAMS
    :returns fuels:        Dictionary with the fuel of each unit
    """
    fuels = {}
    if isinstance(Inputs, list):
        u = Inputs[0]['u']
        f = Inputs[0]['f']
        Fuel = Inputs[1]['Fuel']['val']
    elif isinstance(Inputs, dict):
        u = Inputs['sets']['u']
        f = Inputs['sets']['f']
        Fuel = Inputs['parameters']['Fuel']['val']
    else:
        logging.error('Inputs variable no valid')
        sys.exit(1)
    if shrink:
        uu = shrink_to_64(u)
    else:
        uu = u
    for i in range(len(u)):
        j = np.where(Fuel[i, :] > 0)
        fuels[uu[i]] = f[j[0][0]]
    return fuels


def get_load_data(datain, c):
    """ 
    Get the load curve, the residual load curve, and the net residual load curve of a specific country

    :param datain:  DispaSET inputs formatted into a list of dataframes (output of the ds_to_df function)
    :param c:       Country to consider (e.g. 'BE')
    :return out:    Dataframe with the following columns:
                        Load:               Load curve of the specified country
                        ResidualLoad:       Load minus the production of variable renewable sources
                        NetResidualLoad:    Residual netted from the interconnections with neightbouring countries
    """
    out = pd.DataFrame(index=datain['Demand'].index)
    out['Load'] = datain['Demand']['DA', c]
    # Listing power plants with non-dispatchable power generation:
    VREunits = []
    VRE = np.zeros(len(out))
    for t in ['WTON', 'WTOF', 'PHOT', 'HROR']:
        for u in datain['Technology']:
            if datain['Technology'].loc[t, u]:
                VREunits.append(u)
                VRE = VRE + datain['AvailabilityFactor'][u].values * datain['PowerCapacity'].loc[u, 'PowerCapacity']
    Interconnections = np.zeros(len(out))
    for l in datain['FlowMinimum']:
        ####### Edits #######
        from_node, to_node = l.split('->')
        #        if l[:2] == c:
        if from_node.strip() == c:
            Interconnections = Interconnections - datain['FlowMinimum'][l].values
        #        elif l[-2:] == c:
        elif to_node.strip() == c:
            Interconnections = Interconnections + datain['FlowMinimum'][l].values
        ####### Edits #######

    out['ResidualLoad'] = out['Load'] - VRE
    out['NetResidualLoad'] = out['ResidualLoad'] - Interconnections
    return out


def get_demand(Inputs, c):
    """
    Get the load curve and the residual load curve of a specific country

    :param Inputs:  DispaSET inputs
    :param c:       Country to consider (e.g. 'BE')
    """
    StartDate = Inputs['config']['StartDate']
    StopDate = Inputs['config']['StopDate']
    index = pd.DatetimeIndex(start=pd.datetime(*StartDate), end=pd.datetime(*StopDate), freq='h')

    if isinstance(Inputs, list):
        idx = Inputs[0]['n'].index(c)
        data = Inputs[1]['Demand']['val'][0, idx, range(len(index))]
    elif isinstance(Inputs, dict):
        idx = Inputs['sets']['n'].index(c)
        data = Inputs['parameters']['Demand']['val'][0, idx, range(len(index))]
    else:
        logging.error('Inputs variable no valid')
        sys.exit(1)
    return pd.Series(data, index=index, name=c)


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
    uFuel = unit_fuel(Inputs)

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
    :returns Power:          Dataframe with power generation by fuel
    """
    loc = unit_location(inputs)
    Power = PowerOutput.loc[:, [u for u in PowerOutput.columns if loc[u] == c]]
    return Power


def get_plot_data(inputs, results, c):
    """
    Function that reads the results dataframe of a DispaSET simulation and extract the dispatch data spedific to one country

    :param results:         Pandas dataframe with the results (output of the GdxToDataframe function)
    :param c:               Country to be considered (e.g. 'BE')
    :returns plotdata:       Dataframe with the dispatch data storage and outflows are negative
    """
    tmp = filter_by_country(results['UnitOutputPower'], inputs, c)
    plotdata = aggregate_by_fuel(tmp, inputs)

    if 'OutputStorageInput' in results:
        tmp = filter_by_country(results['OutputStorageInput'], inputs, c)
        plotdata['Storage'] = -tmp.sum(axis=1)
    #        tmp = AggregateByFuel(tmp,inputs,SpecifyFuels=['WAT'])
    #        for col in tmp:
    #            plotdata[col+'_pump'] = -tmp[col]
    else:
        plotdata['Storage'] = 0
    plotdata.fillna(value=0, inplace=True)

    plotdata['FlowIn'] = 0
    plotdata['FlowOut'] = 0
    for col in results['PowerFlow']:
        ####### Edits #######
        from_node, to_node = col.split('->')
        #        if col[-2:] == c:
        if to_node.strip() == c:
            plotdata['FlowIn'] = plotdata['FlowIn'] + results['PowerFlow'][col]
        #    for col in results['OutputFlow']:
        #        if col[:2] == c:
        if from_node.strip() == c:
            plotdata['FlowOut'] = plotdata['FlowOut'] - results['PowerFlow'][col]
        ####### Edits #######

    # re-ordering columns:
    plotdata = plotdata[[u'FlowOut', u'NUC', u'LFG', u'DSL', u'MSW',
                         u'GAS', u'OIL', u'SUN', u'WIN', u'FlowIn', u'HFO']]

    return plotdata


def plot_dispatch(demand, plotdata, level=None, rng=[]):
    """
    Function that plots the dispatch data and the reservoir level as a cumulative sum

    :param demand:      Pandas Series with the demand curve
    :param plotdata:    Pandas Dataframe with the data to be plotted. Negative columns should be at the beginning. Output of the function GetPlotData
    :param level:       Optional pandas series with an aggregated reservoir level for the considered zone.
    :param rng:         Indexes of the values to be plotted. If undefined, the first week is plotted
    """
    import itertools
    # loop through colors.
    colors = itertools.cycle([COLORS[x] for x in plotdata.columns])
    hatches = itertools.cycle(['x', '//', '\\', '/'])

    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    if isinstance(rng,list):
        pdrng = plotdata.index[:min(len(plotdata)-1,7*24)]
    elif rng[0] < plotdata.index[0] or rng[0] > plotdata.index[-1] or rng[-1] < plotdata.index[0] or rng[-1] > plotdata.index[-1]:
        logging.warn('Plotting range is not properly defined, considering the first simulated week')
        pdrng = plotdata.index[:min(len(plotdata)-1,7*24)]
    else:
        pdrng = rng
    alpha = '0.3'

    # find the zero line position:
    cols = plotdata.columns.tolist()
    idx_zero = 0
    tmp = plotdata.iloc[:,idx_zero].mean()
    while tmp <= 0 and idx_zero<len(cols)-1:
        idx_zero += 1
        tmp = plotdata.iloc[:,idx_zero].mean()

    tmp = plotdata[cols[:idx_zero]].sum(axis=1)
    sumplot_neg = pd.DataFrame()
    sumplot_neg['sum'] = tmp
    tmp2 = plotdata[cols[:idx_zero]]
    for col in tmp2:
        sumplot_neg[col] = - tmp2[col]
    sumplot_neg = sumplot_neg.cumsum(axis=1)

    sumplot_pos = plotdata[cols[idx_zero:]].cumsum(axis=1)
    sumplot_pos['zero'] = 0
    sumplot_pos = sumplot_pos[['zero'] + sumplot_pos.columns[:-1].tolist()]

    fig = plt.figure(figsize=(13, 7))

    # Create left axis:
    ax = fig.add_subplot(111)
    ax.plot(pdrng, demand[pdrng], color='k')
    plt.title('Power dispatch for zone ' + demand.name)

    labels = []
    patches = []
    colorlist = []

    # Plot negative values:
    for j in range(idx_zero):
        col1 = sumplot_neg.columns[j]
        col2 = sumplot_neg.columns[j + 1]
        color = next(colors)
        hatch = next(hatches)
        plt.fill_between(pdrng, sumplot_neg.loc[pdrng, col1], sumplot_neg.loc[pdrng, col2], color=color, alpha=alpha,
                         hatch=hatch)
        labels.append(col1)
        patches.append(mpatches.Patch(color=color, alpha=0.3, hatch=hatch, label=col2))
        colorlist.append(color)

    # Plot negative values:
    for j in range(len(sumplot_pos.columns) - 1):
        col1 = sumplot_pos.columns[j]
        col2 = sumplot_pos.columns[j + 1]
        color = next(colors)
        hatch = next(hatches)
        plt.fill_between(pdrng, sumplot_pos.loc[pdrng, col1], sumplot_pos.loc[pdrng, col2], color=color, alpha=alpha,
                         hatch=hatch)
        labels.append(col2)
        patches.append(mpatches.Patch(color=color, alpha=0.3, hatch=hatch, label=col2))
        colorlist.append(color)

    ax.set_ylabel('Power [MW]')
    ax.yaxis.label.set_fontsize(16)

    if level is not None:
        # Create right axis:
        ax2 = fig.add_subplot(111, sharex=ax, frameon=False, label='aa')
        ax2.plot(pdrng, level[pdrng], color='k', alpha=0.3, linestyle='--')
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.set_ylabel('Level [MWh]')
        ax2.yaxis.label.set_fontsize(16)
        line_SOC = mlines.Line2D([], [], color='black', alpha=0.3, label='Reservoir', linestyle='--')

    line_demand = mlines.Line2D([], [], color='black', label='Load')

    if level is None:
        plt.legend(handles=[line_demand] + patches[::-1], loc=4)
    else:
        plt.legend(title='Dispatch for ' + demand.name, handles=[line_demand] + [line_SOC] + patches[::-1], loc=4)


def plot_rug(df_series, on_off=False, cmap='Greys'):
    """Create multiaxis rug plot from pandas Dataframe

    Parameters:
        df_series: 2D pandas with timed index

        on_off:
            * If True all points that are above 0 will be plotted as one color.
            * If False all values will be colored based on their value.

        cmap: palette name (from colorbrewer, matplotlib etc.)

    Returns:
        plot

    Function written by K. Kavvadias
    """

    def format_axis(iax):
        # Formatting: remove all lines (not so elegant)
        iax.axes.spines['top'].set_visible(False)
        iax.spines['right'].set_visible(False)
        iax.spines['left'].set_visible(False)
        iax.spines['bottom'].set_visible(False)
        # iax.xaxis.set_ticks_position('none')
        iax.yaxis.set_ticks_position('none')
        iax.get_yaxis().set_ticks([])
        iax.yaxis.set_label_coords(-.05, -.1)

    def flag_operation(v):
        if np.isnan(v) or v == 0:
            return False
        else:
            return True

    # check if Series or dataframe
    if len(df_series.shape) == 2:
        df_series2 = df_series
        df_series2.columns = [clean_strings(x) for x in df_series.columns]
        rows = len(df_series2.columns)
    elif len(df_series.shape) == 1:
        df_series2 = df_series.to_frame()
        rows = 1
    else:
        raise AssertionError("Has to be either Series or Dataframe")

    fig, axes = plt.subplots(nrows=rows, ncols=1, sharex=True,
                             figsize=(16, 0.25 * rows), squeeze=False,
                             frameon=False, gridspec_kw={'hspace': 0.15})

    for (item, iseries), iax in zip(df_series2.iteritems(), axes.ravel()):
        format_axis(iax)
        iax.set_ylabel(item, rotation=0)
        if iseries.sum() > 0:
            if on_off:
                i_on_off = iseries.apply(flag_operation).replace(False, np.nan)
                i_on_off.plot(ax=iax, style='|', lw=.7, cmap=cmap)
            else:
                x = iseries.index
                y = np.ones(len(iseries))
                iax.scatter(x, y, marker='|', s=100,
                            c=iseries.values * 100, cmap=cmap)

    iax.spines['bottom'].set_visible(True)
    iax.set_xlim(min(x), max(x))


def plot_energy_country_fuel(datain, results, PPindicators):
    """
    Plots the generation for each country, disaggregated by fuel type

    :param datain:          Inputs values (in the dataframe format: output of the function ds_to_df)
    :param results:         Dictionnary with the outputs of the model (output of the function GetResults)
    :param PPindicators:    Por powerplant statistics (output of the function PerPowerPlantIndicators)
    """
    fuels = PPindicators.Fuel.unique()
    zones = PPindicators.Zone.unique()

    GenPerCountry = pd.DataFrame(index=zones, columns=fuels)
    # First make sure that all fuels are present. If not, initialize an empty series
    for f in ['NUC', 'LFG', 'DSL', 'MSW', 'GAS', 'OIL', 'WST', 'SUN', 'WIN', 'FlowIn', 'HFO']:
        if f not in GenPerCountry:
            GenPerCountry[f] = 0
    for c in zones:
        for f in fuels:
            tmp = PPindicators[(PPindicators.Fuel == f) & (PPindicators.Zone == c)]
            GenPerCountry.loc[c, f] = tmp.Generation.sum()
        NetImports = get_imports(results['PowerFlow'], c)
        if NetImports > 0:
            GenPerCountry.loc[c, 'FlowIn'] = NetImports

    GenPerCountry = GenPerCountry[['NUC', 'LFG', 'DSL', 'MSW', 'GAS', 'OIL',
                                   'SUN', 'WIN', 'FlowIn', 'HFO']] / 1E6
    colors = [COLORS[tech] for tech in GenPerCountry.columns]
    ax = GenPerCountry.plot(kind="bar", figsize=(12, 8), stacked=True, color=colors, alpha=0.8, legend='reverse',
                            title='Generation per country (the horizontal lines indicate the demand)')
    ax.set_ylabel('Generation [TWh]')
    demand = datain['Demand']['DA'].sum() / 1E6
    ax.barh(demand, left=ax.get_xticks() - 0.4, width=[0.8] * len(demand), height=demand.values * 0, linewidth=2,
            color='k')
    return ax


def get_sim_results(path='.', cache=False, temp_path='.pickle'):
    """
    This function reads the simulation environment folder once it has been solved and loads
    the input variables together with the results.

    :param path:                Relative path to the simulation environment folder (current path by default)
    :param cache:               If true, caches the simulation results in a pickle file for faster loading the next time
    :param temp_path:            Temporary path to store the cache file
    :returns inputs,results:    Two dictionaries with all the input and outputs
    """

    inputfile = path + '/Inputs.p'
    resultfile = path + '/Results.gdx'

    with open(inputfile, 'rb') as f:
        inputs = pickle.load(f)

    # Clean power plant names:
    inputs['sets']['u'] = clean_strings(inputs['sets']['u'])
    inputs['units'].index = clean_strings(inputs['units'].index.tolist())
    # inputs['units']['Unit'] = clean_strings(inputs['units']['Unit'].tolist())
    gams_dir = inputs['config']['GAMS_folder'].encode() # We need to pass the dir in config if we run it in clusters. PBS script fail to autolocate
    if not os.path.exists(gams_dir):
        logging.warn('The provided path for GAMS (' + gams_dir + ') does not exist. Trying to locate...')
        gams_dir = get_gams_path()
        if not os.path.exists(gams_dir):
            logging.error('GAMS path cannot be located. Simulation is stopped')
            return False
    gams_dir = gams_dir.encode()
    # Load results and store in cache file in the .pickle folder:
    if cache:
        import cPickle
        import hashlib
        m = hashlib.new('md5', resultfile)
        resultfile_hash = m.hexdigest()
        filepath_pickle = str(temp_path + os.path.sep + resultfile_hash + '.p')
        if not os.path.isdir(temp_path):
            os.mkdir(temp_path)
        if not os.path.isfile(filepath_pickle):
            time_pd = 0
        else:
            time_pd = os.path.getmtime(filepath_pickle)
        if os.path.getmtime(resultfile) > time_pd:
            ########
            results = gdx_to_dataframe(gdx_to_list(gams_dir, resultfile, varname='all', verbose=True), fixindex=True, verbose=True)
            # list_gdx = gdx_to_list(gams_dir, resultfile, varname='all', verbose=True)
            #results = gdx_to_dataframe(list_gdx, fixindex=True,
            #                           verbose=True)
            #for element in list_gdx:
            #    try:
            #        if len(list_gdx[element][0]) == 4:
            #            list_gdx[element] = [[ele[0]+'--'+ele[1],ele[2],ele[3]]for ele in list_gdx[element]]
            #    except:
            #        pass
            #results = gdx_to_dataframe(list_gdx, fixindex=True, verbose=True)
            ########
            with open(filepath_pickle, 'wb') as pfile:
                cPickle.dump(results, pfile, protocol=cPickle.HIGHEST_PROTOCOL)
        else:
            with open(filepath_pickle, 'rb') as pfile:
                results = cPickle.load(pfile)
    else:
        results = gdx_to_dataframe(gdx_to_list(gams_dir, resultfile, varname='all', verbose=True), fixindex=True, verbose=True)
            ########

        ######### If a Parameter (or Variable) has dimension of 3 in GAMS, drop the second dimension
        # If a Parameter (or Variable) has dimension of 3 in GAMS, drop the second dimension
        #results = gdx_to_dataframe(gdx_to_list(gams_dir, resultfile, varname='all', verbose=True), fixindex=True, verbose=True)
        #list_gdx = gdx_to_list(gams_dir, resultfile, varname='all', verbose=True)
        #for element in list_gdx:
        #    try:
        #        if len(list_gdx[element][0]) == 4:
        #            list_gdx[element] = [[ele[0], ele[2], ele[3]] for ele in list_gdx[element]]
        #    except:
        #        pass
        #results = gdx_to_dataframe(list_gdx, fixindex=True, verbose=True)
        ########

    # Set datetime index:
    StartDate = inputs['config']['StartDate']
    StopDate = inputs['config']['StopDate']  # last day of the simulation with look-ahead period
    StopDate_long = pd.datetime(*StopDate) + dt.timedelta(days=inputs['config']['LookAhead'])
    index = pd.DatetimeIndex(start=pd.datetime(*StartDate), end=pd.datetime(*StopDate), freq='h')
    index_long = pd.DatetimeIndex(start=pd.datetime(*StartDate), end=StopDate_long, freq='h')

    # Setting the proper index to the result dataframes:
    for key in ['UnitOutputPower', 'UnitOutputPowerInside', 'UnitOutputPowerOutside', 'UnitOutputPowerForNode', 'TotalNodeOperationCost',
                'NetNodeOperationCost', 'UnitCommitment', 'UnitOutputCurtailedPower', 'PowerFlow','NodeOutputShedLoad', 'LostLoad_Reserve2U',
                'LostLoad_MaxPower', 'LostLoad_MinPower', 'LostLoad_RampUp','LostLoad_RampDown', 'LostLoad_Reserve2D','NodeShadowPrice',
                'status', 'UnitFixedCost', 'UnitStartUpCost','UnitShutDownCost', 'UnitRampUpCost', 'UnitRampDownCost', 'UnitVariableCost',
                'UnitOperationCost', 'LocalOutputPower', 'LocalOutputPowerCost','TotalImportedPower', 'TotalImportedPowerCost',
                'ImportedPowerFromNode', 'ImportedPowerFromNodeCost','TotalExportedPower', 'TotalExportedPowerCost','ExportedPowerToNode',
                'ExportedPowerToNodeCost','PowerFlowMax', 'PowerFlowMin', 'TotalNodeOperationCost2', 'TotalSystemCost', 'NetExportedPower',
                'NodeFuelPower', 'NodeFuelCost', 'NodeFuelConsumption', 'NodeExportValue', 'NodeImportValue','LineCongestion','LineCongestion_KSA_GCC',
                'LineCongestion_KW_GCC','LineCongestion_BA_GCC','LineCongestion_QA_GCC','LineCongestion_UAE_GCC','LineCongestion_OM_GCC',
                'LineCongestion_UAE_Salwa','LineCongestion_Ghunan_Salwa','LineCongestion_Ghunan_Alfadhili']:
        if key in results:

            # Drop second level of DataFrame columns for some variables
            if key in ['UnitVariableCost', 'UnitOperationCost'] and results[key].columns.nlevels > 2:
                results[key].columns = results[key].columns.droplevel(1)
            
            if len(results[key]) == len(
                    index_long):  # Case of variables for which the look-ahead period recorded (e.g. the lost loads)
                results[key].index = index_long
            elif len(results[key]) == len(
                    index):  # Case of variables for which the look-ahead is not recorded (standard case)
                results[key].index = index
            else:  # Variables whose index is not complete (sparse formulation)
                results[key].index = index_long[results[key].index - 1]
                if key in ['UnitOutputPower', 'status', 'UnitCommitment', 'UnitOutputCurtailedPower', 'PowerFlow', 'NodeShadowPrice',
                           'NodeOutputShedLoad', 'UnitFixedCost', 'UnitStartUpCost','UnitShutDownCost', 'UnitRampUpCost', 'UnitRampDownCost',
                           'UnitVariableCost', 'UnitOperationCost', 'LocalOutputPower', 'LocalOutputPowerCost', 'UnitOutputPowerInside',
                           'UnitOutputPowerOutside', 'UnitOutputPowerForNode', 'TotalNodeOperationCost', 'NetNodeOperationCost',
                           'TotalImportedPower', 'TotalImportedPowerCost', 'ImportedPowerFromNode', 'ImportedPowerFromNodeCost',
                           'TotalExportedPower', 'TotalExportedPowerCost', 'ExportedPowerToNode', 'ExportedPowerToNodeCost',
                           'PowerFlowMax', 'PowerFlowMin', 'TotalNodeOperationCost2','TotalSystemCost', 'NetExportedPower',
                           'NodeFuelPower', 'NodeFuelCost', 'NodeFuelConsumption', 'NodeExportValue', 'NodeImportValue','LineCongestion',
                           'LineCongestion_KSA_GCC','LineCongestion_KW_GCC','LineCongestion_BA_GCC','LineCongestion_QA_GCC',
                           'LineCongestion_UAE_GCC','LineCongestion_OM_GCC','LineCongestion_UAE_Salwa','LineCongestion_Ghunan_Salwa',
                           'LineCongestion_Ghunan_Alfadhili']:
                    results[key] = results[key].reindex(index).fillna(0)
                    # results[key].fillna(0,inplace=True)
        else:
            results[key] = pd.DataFrame(index=index)

    # Clean power plant names:
    results['UnitOutputPower'].columns = clean_strings(results['UnitOutputPower'].columns.tolist())
    # Remove epsilons:
    if 'ShadowPrice' in results:
        results['NodeShadowPrice'][results['NodeShadowPrice'] == 5e300] = 0

    for key in results['UnitOutputPower']:
        if key not in inputs['units'].index:
            logging.error("Unit '" + key + "' present in the results cannot be found in the input 'units' dataframe")
        if key not in inputs['sets']['u']:
            logging.error("Unit '" + key + "' present in the results cannot be found in the set 'u' from the inputs")

    if "model" in results['status']:
        errors = results['status'][(results['status']['model'] != 1) & (results['status']['model'] != 8)]
        if len(errors) > 0:
            logging.critical('Some simulation errors were encountered. Some results could not be computed, for example at \n \
                            time ' + str(errors.index[0]) + ', with the error message: "' + GAMSstatus('model',errors['model'].iloc[0]) + '". \n \
                            The complete list is available in results["errors"] \n \
                            The optimization might be debugged by activating the Debug flag in the GAMS simulation file and running it')
            for i in errors.index:
                errors.loc[i,'Error Message'] = GAMSstatus('model',errors['model'][i])
            results['errors'] = errors
    return inputs, results


def plot_country(inputs, results, c='', rng=[]):
    """
    Generates plots from the dispa-SET results for one spedific country

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :param c:           Considered country (e.g. 'BE')
    """
    if c =='':
        Nzones = len(inputs['sets']['n'])
        c = inputs['sets']['n'][np.random.randint(Nzones)]
        print('Randomly selected zone for the detailed analysis: '+ c)

    plotdata = get_plot_data(inputs, results, c)

    if 'OutputStorageLevel' in results:
        level = filter_by_country(results['OutputStorageLevel'], inputs, c)
        level = level.sum(axis=1)
    else:
        level = pd.Series(0, index=results['UnitOutputPower'].index)

    demand = get_demand(inputs, c)
    sum_generation = plotdata.sum(axis=1)

    plot_dispatch(demand, plotdata, level, rng=rng)

    # Generation plot:
    CountryGeneration = filter_by_country(results['UnitOutputPower'], inputs, c)

    #plot_rug(CountryGeneration, on_off=False, cmap='Greys')

    return True


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
    dfin = ds_to_df(inputs)

    StartDate = inputs['config']['StartDate']
    StopDate = inputs['config']['StopDate']
    index = pd.DatetimeIndex(start=pd.datetime(*StartDate), end=pd.datetime(*StopDate), freq='h')

    # Aggregated values:
    TotalLoad = dfin['Demand']['DA'].loc[index, :].sum().sum()
    # PeakLoad = inputs['parameters']['Demand']['val'][0,:,idx].sum(axis=0).max()
    PeakLoad = dfin['Demand']['DA'].sum(axis=1).max(axis=0)

    NetImports = -get_imports(results['PowerFlow'], 'RoW')

    Cost_kwh = results['TotalSystemCost'].sum() / (TotalLoad - NetImports)

    print '\nAverage electricity cost : ' + str(Cost_kwh) + ' EUR/MWh'

    for key in ['LostLoad_RampUp', 'LostLoad_Reserve2D', 'LostLoad_MinPower',
                'LostLoad_RampDown', 'LostLoad_Reserve2U', 'LostLoad_MaxPower']:
        LL = results[key].values.sum()
        if LL > 0.0001 * TotalLoad:
            print '\nThere is a significant amount of lost load for ' + key + ': ' + str(
                LL) + ' MWh. The results should be checked carefully'

    print '\nAggregated statistics for the considered area:'
    print 'Total consumption:' + str(TotalLoad / 1E6) + ' TWh'
    print 'Peak load:' + str(PeakLoad) + ' MW'
    print 'Net importations:' + str(NetImports / 1E6) + ' TWh'

    # Country-specific values:
    CountryData = pd.DataFrame(index=inputs['sets']['n'])

    CountryData['Demand'] = dfin['Demand']['DA'].sum(axis=0) / 1E6
    CountryData['PeakLoad'] = dfin['Demand']['DA'].max(axis=0)

    CountryData['NetImports'] = 0
    for c in CountryData.index:
        CountryData.loc[c, 'NetImports'] = get_imports(results['PowerFlow'], str(c)) / 1E6

    CountryData['LoadShedding'] = results['NodeOutputShedLoad'].sum(axis=0) / 1E6
    CountryData['Curtailment'] = results['UnitOutputCurtailedPower'].sum(axis=0) / 1E6
    print '\nCountry-Specific values (in TWh or in MW):'
    print CountryData

    # Congestion:
    Congestion = {}
    if 'PowerFlow' in results:
        for l in results['PowerFlow']:
            if l[:3] != 'RoW' and l[-3:] != 'RoW':
                Congestion[l] = np.sum(
                    (results['PowerFlow'][l] == dfin['FlowMaximum'].loc[results['PowerFlow'].index, l]) & (
                            dfin['FlowMaximum'].loc[results['PowerFlow'].index, l] > 0))
    print "\nNumber of hours of congestion on each line: "
    import pprint
    pprint.pprint(Congestion)
    return {'Cost_kwh': Cost_kwh, 'TotalLoad': TotalLoad, 'PeakLoad': PeakLoad, 'NetImports': NetImports,
            'CountryData': CountryData, 'Congestion': Congestion}


def get_indicators_powerplant(inputs, results):
    """
    Function that analyses the dispa-set results at the power plant level
    Computes the number of startups, the capacity factor, etc

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :returns out:        Dataframe with the main power plants characteristics and the computed indicators
    """
    out = inputs['units'].loc[:, ['PowerCapacity', 'Zone', 'Technology', 'Fuel']]

    out['startups'] = 0
    for u in out.index:
        if u in results['UnitCommitment']:
            # count the number of start-ups
            values = results['UnitCommitment'].loc[:, u].values
            diff = -(values - np.roll(values, 1))
            startups = diff > 0
            out.loc[u, 'startups'] = startups.sum()

    out['CF'] = 0
    out['Generation'] = 0
    for u in out.index:
        if u in results['UnitOutputPower']:
            # count the number of start-ups
            out.loc[u, 'CF'] = results['UnitOutputPower'][u].mean() / out.loc[u, 'PowerCapacity']
            out.loc[u, 'Generation'] = results['UnitOutputPower'][u].sum()
    return out


def ds_to_df(inputs):
    """
    Function that converts the dispaset data format into a dictionary of dataframes

    :param inputs: input file
    :return: dictionary of dataframes
    """

    sets, parameters = inputs['sets'], inputs['parameters']

    # config = parameters['Config']['val']
    try:
        config = inputs['config']
        first_day = pd.datetime(config['StartDate'][0], config['StartDate'][1], config['StartDate'][2], 0)
        last_day = pd.datetime(config['StopDate'][0], config['StopDate'][1], config['StopDate'][2], 23)
        dates = pd.date_range(start=first_day, end=last_day, freq='1h')
    except:
        logging.warn('Could not find the start/stop date information in the inputs. Using an integer index')
        dates = range(1, len(sets['z']) + 1)
    if len(dates) > len(sets['h']):
        logging.error('The provided index has a length of ' + str(len(dates)) + ' while the data only comprises ' + str(
            len(sets['h'])) + ' time elements')
        sys.exit(1)
    elif len(dates) > len(sets['z']):
        logging.warn('The provided index has a length of ' + str(len(dates)) + ' while the simulation was designed for ' + str(
            len(sets['z'])) + ' time elements')
    elif len(dates) < len(sets['z']):
        logging.warn('The provided index has a length of ' + str(len(dates)) + ' while the simulation was designed for ' + str(
            len(sets['z'])) + ' time elements')

    idx = range(len(dates))

    out = {}
    out['sets'] = sets

    # Printing each parameter in a separate sheet and workbook:
    for p in parameters:
        var = parameters[p]
        dim = len(var['sets'])
        if var['sets'][-1] == 'h' and isinstance(dates, pd.DatetimeIndex) and dim > 1:
            # if len(dates) != var['val'].shape[-1]:
            #    sys.exit('The date range in the Config variable (' + str(len(dates)) + ' time steps) does not match the length of the time index (' + str(var['val'].shape[-1]) + ') for variable ' + p)
            var['firstrow'] = 5
        else:
            var['firstrow'] = 1
        if dim == 1:
            if var['sets'][0] == 'h':
                out[p] = pd.DataFrame(var['val'][idx], columns=[p], index=dates)
            else:
                out[p] = pd.DataFrame(var['val'], columns=[p], index=sets[var['sets'][0]])
        elif dim == 2:
            values = var['val']
            list_sets = [sets[var['sets'][0]], sets[var['sets'][1]]]
            if var['sets'][1] == 'h':
                out[p] = pd.DataFrame(values.transpose()[idx, :], index=dates, columns=list_sets[0])
            else:
                out[p] = pd.DataFrame(values.transpose(), index=list_sets[1], columns=list_sets[0])
        elif dim == 3:
            list_sets = [sets[var['sets'][0]], sets[var['sets'][1]], sets[var['sets'][2]]]
            values = var['val']
            values2 = np.zeros([len(list_sets[0]) * len(list_sets[1]), len(list_sets[2])])
            cols = np.zeros([2, len(list_sets[0]) * len(list_sets[1])])
            for i in range(len(list_sets[0])):
                values2[i * len(list_sets[1]):(i + 1) * len(list_sets[1]), :] = values[i, :, :]
                cols[0, i * len(list_sets[1]):(i + 1) * len(list_sets[1])] = i
                cols[1, i * len(list_sets[1]):(i + 1) * len(list_sets[1])] = range(len(list_sets[1]))

            columns = pd.MultiIndex([list_sets[0], list_sets[1]], cols)
            if var['sets'][2] == 'h':
                out[p] = pd.DataFrame(values2.transpose()[idx, :], index=dates, columns=columns)
            else:
                out[p] = pd.DataFrame(values2.transpose(), index=list_sets[2], columns=columns)
        else:
            logging.error('Only three dimensions currently supported. Parameter ' + p + ' has ' + str(dim) + ' dimensions.')
            sys.exit(1)
    return out

def get_units_operation_cost(inputs, results, zone = 'All'):
    """
    Function that computes the operation cost for each power unit at each instant of time from the DispaSET results
    Operation cost includes: CostFixed + CostStartUp + CostShutDown + CostRampUp + CostRampDown + CostVariable

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :param zone:        country or zone to be considered
    :returns UnitOperationCost:       Dataframe with the the power units in columns and the operatopn cost at each instant in rows
    """
    datain = ds_to_df(inputs)

    #DataFrame with startup times for each unit (1 for startup)
    StartUps = results['UnitCommitment'].copy()
    for u in StartUps:
        values = StartUps.loc[:, u].values
        diff = -(np.roll(values, 1) - values )
        diff[diff <= 0] = 0
        StartUps[u] = diff

    #DataFrame with shutdown times for each unit (1 for shutdown)
    ShutDowns = results['UnitCommitment'].copy()
    for u in ShutDowns:
        values = ShutDowns.loc[:, u].values
        diff = (np.roll(values, 1) - values )
        diff[diff <= 0] = 0
        ShutDowns[u] = diff

    #DataFrame with ramping up levels for each unit at each instant (0 for ramping-down & leveling out)
    RampUps = results['UnitOutputPower'].copy()
    for u in RampUps:
        values = RampUps.loc[:, u].values
        diff = -(np.roll(values, 1) - values )
        diff[diff <= 0] = 0
        RampUps[u] = diff

    #DataFrame with ramping down levels for each unit at each instant (0 for ramping-up & leveling out)
    RampDowns = results['UnitOutputPower'].copy()
    for u in RampDowns:
        values = RampDowns.loc[:, u].values
        diff = (np.roll(values, 1) - values )
        diff[diff <= 0] = 0
        RampDowns[u] = diff

    FiexedCost = results['UnitCommitment'].copy()
    StartUpCost = results['UnitCommitment'].copy()
    ShutDownCost = results['UnitCommitment'].copy()
    RampUpCost = results['UnitCommitment'].copy()
    RampDownCost = results['UnitCommitment'].copy()
    VariableCost = results['UnitCommitment'].copy()
    UnitOperationCost = results['UnitCommitment'].copy()

    OperatedUnitList = results['UnitCommitment'].columns
    for u in OperatedUnitList:
        unit_indexNo = inputs['units'].index.get_loc(u)
        FiexedCost.loc[:,[u]] = np.array(results['UnitCommitment'].loc[:,[u]])*inputs['parameters']['CostFixed']['val'][unit_indexNo]
        StartUpCost.loc[:,[u]] = np.array(StartUps.loc[:,[u]])*inputs['parameters']['CostStartUp']['val'][unit_indexNo]
        ShutDownCost.loc[:,[u]] = np.array(ShutDowns.loc[:,[u]])*inputs['parameters']['CostShutDown']['val'][unit_indexNo]
        RampUpCost.loc[:,[u]] = np.array(RampUps.loc[:,[u]])*inputs['parameters']['CostRampUp']['val'][unit_indexNo]
        RampDownCost.loc[:,[u]] = np.array(RampDowns.loc[:,[u]])*inputs['parameters']['CostRampDown']['val'][unit_indexNo]
        VariableCost.loc[:,[u]] = np.array(datain['CostVariable'].loc[:,[u]])*np.array(results['UnitOutputPower'][u]).reshape(-1,1)

    UnitOperationCost = FiexedCost+StartUpCost+ShutDownCost+RampUpCost+RampDownCost+VariableCost

    if zone == 'All':
        UnitOperationCost = UnitOperationCost
    else:
        UnitOperationCost = UnitOperationCost.filter(list(inputs['units'][inputs['units']['Zone'] == zone].index))

    return UnitOperationCost

def get_sorted_units(inputs, results, zone='All', sortby='CostVariable', ascending=True):
    """
    Function that sorts power units list based on specific parameter (default is variable cost)

    in case of sorting by variable cost, variable cost is assumed to be fixed by fixing efficiency
    and fuel cost (taken from the first row of each unit's variable cost time series)

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :param zone:        country or zone to be considered
    :param sortby:      parameter in which unit list is sorted based on
    :param ascending:   'True' when sorting in ascending order ('False' for descending order)
    :returns SortedInputs:       Dataframe with the the sorted power units in rows and their characteristics in columns
    """
    datain = ds_to_df(inputs)

    try:
        DataType = inputs['units'][sortby].dtype
    except:
        DataType = 'String'
        pass

    if DataType == 'float64' or DataType == 'int' or sortby is 'CostVariable':
        if zone is not 'All':
            if sortby is 'CostVariable':
                InputsNew = inputs['units'][inputs['units']['Zone']==zone].copy()
                for u in InputsNew.index:
                    InputsNew.at[u , 'CostVariable'] = datain['CostVariable'][u][0]
                SortedInputs = InputsNew.sort_values(by=['CostVariable'], ascending=ascending)
            else:
                InputsNew = inputs['units'][inputs['units']['Zone']==zone].copy()
                SortedInputs = InputsNew.sort_values(by=[sortby], ascending=ascending)

        elif zone is 'All':
            if sortby is 'CostVariable':
                InputsNew = inputs['units'].copy()
                for u in InputsNew.index:
                    InputsNew.at[u , 'CostVariable'] = datain['CostVariable'][u][0]
                SortedInputs = InputsNew.sort_values(by=['CostVariable'], ascending=ascending)
            else:
                InputsNew = inputs['units'].copy()
                SortedInputs = InputsNew.sort_values(by=[sortby], ascending=ascending)
        else:
            logging.error('zone entered does not exist')
    else:
        logging.error('column for the sorting parameter is not of dtype= "float64" nor "int" ')

    return SortedInputs

def get_committed_sorted_units(inputs, results, sortby='CostVariable', ascending=False):
    """
    Function that sorts power units that were committed to a list based on specific parameter (default is variable cost)

    in case of sorting by variable cost, variable cost is assumed to be fixed by fixing efficiency
    and fuel cost (taken from the first row of each unit's variable cost time series)

    :param inputs:                         DispaSET inputs
    :param results:                        DispaSET results
    :param sortby:                         parameter in which unit list is sorted based on
    :param ascending:                      'False' when sorting in descending order ('True' for ascending order)
    :returns committedSortedUnitList:      Dictionary of DataFrames for each zone; each DataFrame with the the sorted
                                           committed power units in rows and their characteristics in columns
    """
    SortedUnitList = {}
    committedSortedUnitList = {}
    for zone in inputs['config']['zones']:
        committedUnitMeritOrderList = []
        for unit in get_sorted_units(inputs, results, zone=zone, sortby=sortby, ascending=ascending).index:
            if unit in results['UnitCommitment'].columns:
                committedUnitMeritOrderList.append(unit)
        SortedUnitList[zone] = get_sorted_units(inputs, results, zone=zone, sortby=sortby, ascending=ascending)
        committedSortedUnitList[zone] = SortedUnitList[zone].filter(committedUnitMeritOrderList, axis=0).filter(committedUnitMeritOrderList, axis=0)

    return committedSortedUnitList

def get_nodes_generation_cost(inputs, results):
    """
    Function that computes the generation cost for each zone at each instant of time from the DispaSET results
    Generation cost includes: CostFixed + CostStartUp + CostShutDown + CostRampUp + CostRampDown + CostVariable

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :returns out:       Dataframe with the the zones in columns and the generation cost at each instant in rows
    """
    AllUnitsCosts = get_units_operation_cost(inputs, results)
    NodesGenerationCost = pd.DataFrame(index = AllUnitsCosts.index, columns = inputs['config']['zones'])
    for zone in inputs['config']['zones']:
        NodesGenerationCost[zone] = get_units_operation_cost(inputs, results, zone = zone).sum(axis=1).values

    return NodesGenerationCost

def Local_generation_minus_demand(inputs, results):
    """
    This function determines generation self sufficiency by computing (local generation - local demand)
    for each zone (country)

    :param Inputs:          Dispaset inputs
    :param results:         Pandas Dataframe with the results (output of the GdxToDataframe function)
    :returns PowerSelf:     Dataframe with local power generation minus demand in rows and zones (countries) in columns
    """
    PowerSelf = pd.DataFrame(index=results['UnitOutputPower'].index,columns=inputs['config']['zones'])
    for zone in inputs['config']['zones']:
        Local_generation_level = results['UnitOutputPower'].filter(list(inputs['units'][inputs['units']['Zone'] == zone].index)).sum(axis=1).values
        TotalDemand = get_demand(inputs, zone).values
        PowerSelfarray = Local_generation_level - TotalDemand
        PowerSelfarray[abs(PowerSelfarray) < abs(1e-8)] = 0       #Zero out very small values
        PowerSelf[zone] = PowerSelfarray

    return PowerSelf
