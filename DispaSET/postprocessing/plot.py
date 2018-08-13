import logging
import sys

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from ..misc.str_handler import clean_strings
from ..common import commons


from .postprocessing import get_imports, get_plot_data, filter_by_country


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