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
from ..common import commonvars

 # get color definitions: 
commons= commonvars()

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
                 12: u'Error achieved \u2013 No cause known',
                 13: u'Error achieved \u2013 No solution attained',
                 14: u'No solution returned',
                 15: u'Feasible in a CNS models',
                 16: u'Locally feasible in a CNS models',
                 17: u'Singular in a CNS models',
                 18: u'Unbounded \u2013 no solution',
                 19: u'Infeasible \u2013 no solution'}
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
        if col[-3:] == c:
            plotdata['FlowIn'] = plotdata['FlowIn'] + results['OutputFlow'][col]
    for col in results['OutputFlow']:
        if col[:3] == c:
            plotdata['FlowOut'] = plotdata['FlowOut'] - results['OutputFlow'][col]

    # re-ordering columns:
    OrderedColumns = [col for col in commons['MeritOrder'] if col in plotdata.columns]
    plotdata = plotdata[OrderedColumns]
    
    # remove empty columns:
    for col in plotdata.columns:
        if plotdata[col].max() == 0 and plotdata[col].min()==0:
            del plotdata[col]

    return plotdata


def plot_dispatch_safe(demand, plotdata, level=None, curtailment=None, rng=None):
    """
    Function that plots the dispatch data and the reservoir level as a cumulative sum.
    In this case, the Pandas index is not used since it can cause a bug in matplotlib
    
    :param demand:      Pandas Series with the demand curve
    :param plotdata:    Pandas Dataframe with the data to be plotted. Negative columns should be at the beginning. Output of the function GetPlotData
    :param level:       Optional pandas series with an aggregated reservoir level for the considered zone.
    :param curtailment: Optional pandas series with the value of the curtailment 
    :param rng:         Indexes of the values to be plotted. If undefined, the first week is plotted
    """
    logging.critical('Error in when drawing the dispatch plot. Trying a safer mode without x-axis')

    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    if rng is None:
        pdrng = plotdata.index[:min(len(plotdata)-1,7*24)]
    elif not type(rng) == type(demand.index):
        logging.error('The "rng" variable must be a pandas DatetimeIndex')
        sys.exit(1)
    elif rng[0] < plotdata.index[0] or rng[0] > plotdata.index[-1] or rng[-1] < plotdata.index[0] or rng[-1] > plotdata.index[-1]:
        logging.warn('Plotting range is not properly defined, considering the first simulated week')
        pdrng = plotdata.index[:min(len(plotdata)-1,7*24)]
    else:
        pdrng = rng
    alpha = '0.3'        
    idx = [d.to_pydatetime() for d in pdrng]  

    # Netting the interconnections:
    if 'FlowIn' in plotdata and 'FlowOut' in plotdata:
        plotdata['FlowIn'],plotdata['FlowOut'] = (np.maximum(0,plotdata['FlowIn']-plotdata['FlowOut']),np.maximum(0,plotdata['FlowOut']-plotdata['FlowIn']))
      
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
#    ax.set_ylim([-10000,15000])
    ax.plot(idx, demand[pdrng].values, color='k')
    plt.title('Power dispatch for country ' + demand.name[1])

    labels = []
    patches = []
    colorlist = []

#    # Plot negative values:
    for j in range(idx_zero):
        col1 = sumplot_neg.columns[j]
        col2 = sumplot_neg.columns[j + 1]
        color = commons['colors'][col2]
        hatch = commons['hatches'][col2]
        plt.fill_between(idx, sumplot_neg.loc[pdrng, col1].values, sumplot_neg.loc[pdrng, col2].values, color=color, alpha=alpha,
                         hatch=hatch)
        labels.append(col1)
        patches.append(mpatches.Patch(color=color, alpha=0.3, hatch=hatch, label=col2))
        colorlist.append(color)

    # Plot Positive values:
    for j in range(len(sumplot_pos.columns) - 1):
        col1 = sumplot_pos.columns[j]
        col2 = sumplot_pos.columns[j + 1]
        color = commons['colors'][col2]
        hatch = commons['hatches'][col2]
        plt.fill_between(idx, sumplot_pos.loc[pdrng, col1].values, sumplot_pos.loc[pdrng, col2].values, color=color, alpha=alpha,
                         hatch=hatch)
        labels.append(col2)
        patches.append(mpatches.Patch(color=color, alpha=0.3, hatch=hatch, label=col2))
        colorlist.append(color)
    
    # Plot curtailment:    
    if isinstance(curtailment,pd.Series):
        if not curtailment.index.equals(demand.index):
            logging.error('The curtailment time series must have the same index as the demand')
            sys.exit(1)
        color = 'red'
        plt.fill_between(idx, sumplot_neg.loc[pdrng, 'sum'].values-curtailment[pdrng].values, sumplot_neg.loc[pdrng, 'sum'].values, color='red', alpha=0.7)
        labels.append('Curtailment')
        patches.append(mpatches.Patch(color='red', alpha=0.7, label='Curtailment'))
   
    plt.xticks(rotation=45)
    ax.set_ylabel('Power [MW]')
    ax.yaxis.label.set_fontsize(16)

    if level is not None:
        # Create right axis:
        ax2 = fig.add_subplot(111, sharex=ax, frameon=False, label='aa')
        ax2.plot(idx, level[pdrng].values, color='k', alpha=0.3, linestyle='--')
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.set_ylabel('Level [MWh]')
        ax2.yaxis.label.set_fontsize(16)
        line_SOC = mlines.Line2D([], [], color='black', alpha=0.3, label='Reservoir', linestyle='--')

    plt.xticks(rotation=45)
    line_demand = mlines.Line2D([], [], color='black', label='Load')
    plt.legend(handles=[line_demand] + patches[::-1], loc=4)
    if level is None:
        plt.legend(handles=[line_demand] + patches[::-1], loc=4)
    else:
        plt.legend(title='Dispatch for ' + demand.name[1], handles=[line_demand] + [line_SOC] + patches[::-1], loc=4)

def plot_dispatch(demand, plotdata, level=None, curtailment=None, rng=None):
    """
    Function that plots the dispatch data and the reservoir level as a cumulative sum
    
    :param demand:      Pandas Series with the demand curve
    :param plotdata:    Pandas Dataframe with the data to be plotted. Negative columns should be at the beginning. Output of the function GetPlotData
    :param level:       Optional pandas series with an aggregated reservoir level for the considered zone.
    :param curtailment: Optional pandas series with the value of the curtailment 
    :param rng:         Indexes of the values to be plotted. If undefined, the first week is plotted
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    if rng is None:
        pdrng = plotdata.index[:min(len(plotdata)-1,7*24)]
    elif not type(rng) == type(demand.index):
        logging.error('The "rng" variable must be a pandas DatetimeIndex')
        sys.exit(1)
    elif rng[0] < plotdata.index[0] or rng[0] > plotdata.index[-1] or rng[-1] < plotdata.index[0] or rng[-1] > plotdata.index[-1]:
        logging.warn('Plotting range is not properly defined, considering the first simulated week')
        pdrng = plotdata.index[:min(len(plotdata)-1,7*24)]
    else:
        pdrng = rng
    alpha = '0.3'        
    
    # Netting the interconnections:
    if 'FlowIn' in plotdata and 'FlowOut' in plotdata:
        plotdata['FlowOut'],plotdata['FlowIn'] = (np.minimum(0,plotdata['FlowIn']+plotdata['FlowOut']),np.maximum(0,plotdata['FlowOut']+plotdata['FlowIn']))
        
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
#    ax.set_ylim([-10000,15000])
    ax.plot(pdrng, demand[pdrng], color='k')
    plt.title('Power dispatch for country ' + demand.name[1])

    labels = []
    patches = []
    colorlist = []

#    # Plot negative values:
    for j in range(idx_zero):
        col1 = sumplot_neg.columns[j]
        col2 = sumplot_neg.columns[j + 1]
        color = commons['colors'][col2]
        hatch = commons['hatches'][col2]
        plt.fill_between(pdrng, sumplot_neg.loc[pdrng, col1], sumplot_neg.loc[pdrng, col2], color=color, alpha=alpha,
                         hatch=hatch)
        labels.append(col1)
        patches.append(mpatches.Patch(color=color, alpha=0.3, hatch=hatch, label=col2))
        colorlist.append(color)

    # Plot Positive values:
    for j in range(len(sumplot_pos.columns) - 1):
        col1 = sumplot_pos.columns[j]
        col2 = sumplot_pos.columns[j + 1]
        color = commons['colors'][col2]
        hatch = commons['hatches'][col2]
        plt.fill_between(pdrng, sumplot_pos.loc[pdrng, col1], sumplot_pos.loc[pdrng, col2], color=color, alpha=alpha,
                         hatch=hatch)
        labels.append(col2)
        patches.append(mpatches.Patch(color=color, alpha=0.3, hatch=hatch, label=col2))
        colorlist.append(color)
    
    # Plot curtailment:    
    if isinstance(curtailment,pd.Series):
        if not curtailment.index.equals(demand.index):
            logging.error('The curtailment time series must have the same index as the demand')
            sys.exit(1)
        color = 'red'
        plt.fill_between(pdrng, sumplot_neg.loc[pdrng, 'sum']-curtailment[pdrng], sumplot_neg.loc[pdrng, 'sum'], color='red', alpha=0.7)
        labels.append('Curtailment')
        patches.append(mpatches.Patch(color='red', alpha=0.7, label='Curtailment'))
   

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
    plt.legend(handles=[line_demand] + patches[::-1], loc=4)
    if level is None:
        plt.legend(handles=[line_demand] + patches[::-1], loc=4)
    else:
        plt.legend(title='Dispatch for ' + demand.name[1], handles=[line_demand] + [line_SOC] + patches[::-1], loc=4)


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
                x = [d.to_pydatetime() for d in iseries.index]
                y = np.ones(len(iseries))
                iax.scatter(x, y, marker='|', s=100,
                            c=iseries.values * 100, cmap=cmap)

    iax.spines['bottom'].set_visible(True)
    plt.xticks(rotation=45)
    #iax.set_xlim(min(x), max(x))



def plot_energy_country_fuel(inputs, results, PPindicators):
    """
    Plots the generation for each country, disaggregated by fuel type

    :param results:         Dictionnary with the outputs of the model (output of the function GetResults)
    :param PPindicators:    Por powerplant statistics (output of the function get_indicators_powerplant)
    """
    fuels = PPindicators.Fuel.unique()
    countries = PPindicators.Zone.unique()

    GenPerCountry = pd.DataFrame(index=countries, columns=fuels)
    # First make sure that all fuels are present. If not, initialize an empty series
    for f in commons['Fuels'] + ['FlowIn']:
        if f not in GenPerCountry:
            GenPerCountry[f] = 0
    for c in countries:
        for f in fuels:
            tmp = PPindicators[(PPindicators.Fuel == f) & (PPindicators.Zone == c)]
            GenPerCountry.loc[c, f] = tmp.Generation.sum()
        NetImports = get_imports(results['OutputFlow'], c)
        if NetImports > 0:
            GenPerCountry.loc[c, 'FlowIn'] = NetImports

    cols = [col for col in commons['MeritOrder'] if col in GenPerCountry]
    GenPerCountry = GenPerCountry[cols] / 1E6
    colors = [commons['colors'][tech] for tech in GenPerCountry.columns]
    ax = GenPerCountry.plot(kind="bar", figsize=(12, 8), stacked=True, color=colors, alpha=0.8, legend='reverse',
                            title='Generation per country (the horizontal lines indicate the demand)')
    ax.set_ylabel('Generation [TWh]')
    demand = inputs['param_df']['Demand']['DA'].sum() / 1E6
    ax.barh(demand, left=ax.get_xticks() - 0.4, width=[0.8] * len(demand), height=ax.get_ylim()[1]*0.005, linewidth=2,
            color='k')
    return ax


def plot_country_capacities(inputs,plot=True):
    """
    Plots the installed capacity for each country, disaggregated by fuel type

    :param inputs:         Dictionnary with the inputs of the model (output of the function GetResults)
    """
    units = inputs['units']
    CountryFuels = {}
    for u in units.index:
        CountryFuels[(units.Zone[u],units.Fuel[u])] = (units.Zone[u],units.Fuel[u])
    
    PowerCapacity = pd.DataFrame(columns=inputs['sets']['f'],index=inputs['sets']['n'])
    StorageCapacity = pd.DataFrame(columns=inputs['sets']['f'],index=inputs['sets']['n'])
    for n,f in CountryFuels:
        idx = ((units.Zone == n) & (units.Fuel==f))
        PowerCapacity.loc[n,f] = (units.PowerCapacity[idx]*units.Nunits[idx]).sum()
        StorageCapacity.loc[n,f] = (units.StorageCapacity[idx]*units.Nunits[idx]).sum()

    cols = [col for col in commons['MeritOrder'] if col in PowerCapacity]
    PowerCapacity = PowerCapacity[cols]
    if plot:
        colors = [commons['colors'][tech] for tech in PowerCapacity.columns]
        ax = PowerCapacity.plot(kind="bar", figsize=(12, 8), stacked=True, color=colors, alpha=0.8, legend='reverse',
                                title='Installed capacity per country (the horizontal lines indicate the peak demand)')
        ax.set_ylabel('Capacity [MW]')
        demand = inputs['param_df']['Demand']['DA'].max() 
        ax.barh(demand, left=ax.get_xticks() - 0.4, width=[0.8] * len(demand), height=ax.get_ylim()[1]*0.005, linewidth=2,
                color='k')
    return {'PowerCapacity':PowerCapacity,'StorageCapacity':StorageCapacity}



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
    
    # Add the formated parameters in the inputs variable if not already present:
    if not 'param_df' in inputs:
        inputs['param_df'] = ds_to_df(inputs)


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
            results = gdx_to_dataframe(gdx_to_list(gams_dir, resultfile, varname='all', verbose=True), fixindex=True,
                                       verbose=True)
            with open(filepath_pickle, 'wb') as pfile:
                cPickle.dump(results, pfile, protocol=cPickle.HIGHEST_PROTOCOL)
        else:
            with open(filepath_pickle, 'rb') as pfile:
                results = cPickle.load(pfile)
    else:
        results = gdx_to_dataframe(gdx_to_list(gams_dir, resultfile, varname='all', verbose=True), fixindex=True,
                                   verbose=True)

    # Set datetime index:
    StartDate = inputs['config']['StartDate']
    StopDate = inputs['config']['StopDate']  # last day of the simulation with look-ahead period
    StopDate_long = pd.datetime(*StopDate) + dt.timedelta(days=inputs['config']['LookAhead'])
    index = pd.DatetimeIndex(start=pd.datetime(*StartDate), end=pd.datetime(*StopDate), freq='h')
    index_long = pd.DatetimeIndex(start=pd.datetime(*StartDate), end=StopDate_long, freq='h')

    # Setting the proper index to the result dataframes:
    for key in ['OutputPower', 'OutputSystemCost', 'OutputCommitted', 'OutputCurtailedPower', 'OutputFlow',
                'OutputShedLoad', 'OutputSpillage', 'OutputStorageLevel', 'OutputStorageInput', 'LostLoad_2U', 'LostLoad_3U',
                'LostLoad_MaxPower', 'LostLoad_MinPower', 'LostLoad_RampUp', 'LostLoad_RampDown', 'LostLoad_2D',
                'ShadowPrice', 'OutputHeat', 'OutputHeatSlack','status']:
        if key in results:
            if len(results[key]) == len(
                    index_long):  # Case of variables for which the look-ahead period recorded (e.g. the lost loads)
                results[key].index = index_long
            elif len(results[key]) == len(
                    index):  # Case of variables for which the look-ahead is not recorded (standard case)
                results[key].index = index
            else:  # Variables whose index is not complete (sparse formulation)
                results[key].index = index_long[results[key].index - 1]
                if key in ['OutputPower', 'OutputSystemCost', 'OutputCommitted', 'OutputCurtailedPower', 'OutputFlow',
                           'OutputShedLoad', 'OutputSpillage', 'OutputStorageLevel', 'OutputStorageInput','OutputHeat', 'OutputHeatSlack']:
                    results[key] = results[key].reindex(index).fillna(0)
                    # results[key].fillna(0,inplace=True)
        else:
            results[key] = pd.DataFrame(index=index)

    # Clean power plant names:
    results['OutputPower'].columns = clean_strings(results['OutputPower'].columns.tolist())
    # Remove epsilons:
    if 'ShadowPrice' in results:
        results['ShadowPrice'][results['ShadowPrice'] == 5e300] = 0

    for key in results['OutputPower']:
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


def plot_country(inputs, results, c='', rng=None):
    """
    Generates plots from the dispa-SET results for one specific country

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :param c:           Considered country (e.g. 'BE')
    """
    if c =='':
        Nzones = len(inputs['sets']['n'])
        c = inputs['sets']['n'][np.random.randint(Nzones)]
        print('Randomly selected zone for the detailed analysis: '+ c)
    elif c not in inputs['sets']['n']:
        logging.critical('Country ' + c + ' is not in the results')
        Nzones = len(inputs['sets']['n'])
        c = inputs['sets']['n'][np.random.randint(Nzones)]
        logging.critical('Randomly selected country: '+ c)

    plotdata = get_plot_data(inputs, results, c)

    if 'OutputStorageLevel' in results:
        level = filter_by_country(results['OutputStorageLevel'], inputs, c)
        level = level.sum(axis=1)
    else:
        level = pd.Series(0, index=results['OutputPower'].index)

    demand = inputs['param_df']['Demand'][('DA', c)]
    sum_generation = plotdata.sum(axis=1)
    diff = (sum_generation - demand).abs()
    if diff.max() > 0.01 * demand.max():
        logging.critical('There is up to ' + str(diff.max()/demand.max()*100) + '% difference in the instantaneous energy balance of country ' + c)

    try:
        if 'OutputCurtailedPower' in results and c in results['OutputCurtailedPower']:
            plot_dispatch(demand, plotdata, level, curtailment=results['OutputCurtailedPower'][c], rng=rng)
        else:
            plot_dispatch(demand, plotdata, level, rng=rng)
    except:
        if 'OutputCurtailedPower' in results and c in results['OutputCurtailedPower']:
            plot_dispatch_safe(demand, plotdata, level, curtailment=results['OutputCurtailedPower'][c], rng=rng)
        else:
            plot_dispatch_safe(demand, plotdata, level, rng=rng)

    # Generation plot: 
    CountryGeneration = filter_by_country(results['OutputPower'], inputs, c)

    plot_rug(CountryGeneration, on_off=False, cmap='Greys')

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
        timeindex=True
    except:
        logging.warn('Could not find the start/stop date information in the inputs. Using an integer index')
        dates = range(1, len(sets['z']) + 1)
        timeindex=False
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
        if var['sets'][-1] == 'h' and timeindex and dim > 1:
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


