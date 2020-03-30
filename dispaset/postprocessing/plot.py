import logging
import sys

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from ..misc.str_handler import clean_strings
from ..common import commons

from .postprocessing import get_imports, get_plot_data, filter_by_zone, filter_by_tech, filter_by_storage


def plot_dispatch(demand, plotdata, level=None, curtailment=None, shedload=None, shiftedload=None, rng=None,
                  alpha=None, figsize=(13, 6)):
    """
    Function that plots the dispatch data and the reservoir level as a cumulative sum

    :param demand:      Pandas Series with the demand curve
    :param plotdata:    Pandas Dataframe with the data to be plotted. Negative columns should be at the beginning. Output of the function GetPlotData
    :param level:       Optional pandas series/dataframe with an (dis)aggregated reservoir level for the considered zone.
    :param curtailment: Optional pandas series with the value of the curtailment
    :param shedload:    Optional pandas series with the value of the shed load
    :param rng:         Indexes of the values to be plotted. If undefined, the first week is plotted
    """
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    pd.plotting.register_matplotlib_converters()

    if rng is None:
        pdrng = plotdata.index[:min(len(plotdata) - 1, 7 * 24)]
    elif not type(rng) == type(demand.index):
        logging.error('The "rng" variable must be a pandas DatetimeIndex')
        raise ValueError()
    elif rng[0] < plotdata.index[0] or rng[0] > plotdata.index[-1] or rng[-1] < plotdata.index[0] or rng[-1] > \
            plotdata.index[-1]:
        logging.warning('Plotting range is not properly defined, considering the first simulated week')
        pdrng = plotdata.index[:min(len(plotdata) - 1, 7 * 24)]
    else:
        pdrng = rng

    # Netting the interconnections:
    if 'FlowIn' in plotdata and 'FlowOut' in plotdata:
        plotdata['FlowOut'], plotdata['FlowIn'] = (np.minimum(0, plotdata['FlowIn'] + plotdata['FlowOut']),
                                                   np.maximum(0, plotdata['FlowOut'] + plotdata['FlowIn']))

    # find the zero line position:
    cols = plotdata.columns.tolist()
    idx_zero = 0
    tmp = plotdata.iloc[:, idx_zero].mean()
    while tmp <= 0 and idx_zero < len(cols) - 1:
        idx_zero += 1
        tmp = plotdata.iloc[:, idx_zero].mean()

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

    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(figsize), frameon=True,  # 14 4*2
                             gridspec_kw={'height_ratios': [2.7, .8], 'hspace': 0.04})

    # Create left axis:
    #    ax.set_ylim([-10000,15000])
    axes[0].plot(pdrng, demand[pdrng], color='k')
    axes[0].set_xlim(pdrng[0], pdrng[-1])

    fig.suptitle('Power dispatch for zone ' + demand.name[1])

    # Define labels, patches and colors
    labels = []
    patches = []
    colorlist = []

    # Plot reservoir levels (either separated or as one value)
    if level is not None:
        if isinstance(level, pd.DataFrame):
            cols_lvl = level.columns.tolist()
            sumplot_lev = level[cols_lvl[0:]].cumsum(axis=1)
            sumplot_lev['zero'] = 0
            sumplot_lev = sumplot_lev[['zero'] + sumplot_lev.columns[:-1].tolist()]
            for j in range(len(sumplot_lev.columns) - 1):
                col3 = sumplot_lev.columns[j]
                col4 = sumplot_lev.columns[j + 1]
                rez_color = commons['colors'][col4]
                rez_hatch = commons['hatches'][col4]
                axes[1].plot(pdrng, sumplot_lev.loc[pdrng, col4], color='k', alpha=alpha, linestyle=':')
                axes[1].fill_between(pdrng, sumplot_lev.loc[pdrng, col3], sumplot_lev.loc[pdrng, col4],
                                     facecolor=rez_color, alpha=0.3)
                labels.append(col4)
                patches.append(mpatches.Patch(facecolor=rez_color, alpha=0.3, label=col4))
                colorlist.append(rez_color)
        elif isinstance(level, pd.Series):
            # Create lower axis:
            axes[1].plot(pdrng, level[pdrng], color='k', alpha=alpha, linestyle=':')
            axes[1].fill_between(pdrng, 0, level[pdrng],
                                 facecolor=commons['colors']['WAT'], alpha=.3)
        axes[1].set_ylabel('Level [GWh]')
        axes[1].yaxis.label.set_fontsize(12)
        line_SOC = mlines.Line2D([], [], color='black', alpha=alpha, label='Reservoir', linestyle=':')

    # Plot negative values:
    for j in range(idx_zero):
        col1 = sumplot_neg.columns[j]
        col2 = sumplot_neg.columns[j + 1]
        color = commons['colors'][col2]
        hatch = commons['hatches'][col2]
        axes[0].fill_between(pdrng, sumplot_neg.loc[pdrng, col1], sumplot_neg.loc[pdrng, col2], facecolor=color,
                             alpha=alpha, hatch=hatch)
        if col2 not in labels:
            labels.append(col1)
            patches.append(mpatches.Patch(facecolor=color, alpha=alpha, hatch=hatch, label=col2))
            colorlist.append(color)

    # Plot Positive values:
    for j in range(len(sumplot_pos.columns) - 1):
        col1 = sumplot_pos.columns[j]
        col2 = sumplot_pos.columns[j + 1]
        color = commons['colors'][col2]
        hatch = commons['hatches'][col2]
        axes[0].fill_between(pdrng, sumplot_pos.loc[pdrng, col1], sumplot_pos.loc[pdrng, col2], facecolor=color,
                             alpha=alpha,
                             hatch=hatch)
        labels.append(col2)
        patches.append(mpatches.Patch(facecolor=color, alpha=alpha, hatch=hatch, label=col2))
        colorlist.append(color)

    # Plot curtailment:
    if isinstance(curtailment, pd.Series):
        if not curtailment.index.equals(demand.index):
            logging.error('The curtailment time series must have the same index as the demand')
            sys.exit(1)
        axes[0].fill_between(pdrng, sumplot_neg.loc[pdrng, 'sum'] - curtailment[pdrng], sumplot_neg.loc[pdrng, 'sum'],
                             facecolor=commons['colors']['curtailment'])
        labels.append('Curtailment')
        patches.append(mpatches.Patch(facecolor=commons['colors']['curtailment'], label='Curtailment'))

    axes[0].set_ylabel('Power [GW]')
    axes[0].yaxis.label.set_fontsize(12)

    load_change = pd.Series(0, index=demand.index)
    load_changed = False
    if isinstance(shedload, pd.Series):
        if not shedload.index.equals(demand.index):
            logging.critical('The shedload time series must have the same index as the demand')
            sys.exit(1)
        load_change += -shedload
        load_changed = True
    if isinstance(shiftedload, pd.Series):
        if not shiftedload.index.equals(demand.index):
            logging.critical('The shiftedload time series must have the same index as the demand')
            sys.exit(1)
        load_change += -shiftedload
        load_changed = True
    reduced_demand = demand + load_change
    axes[0].plot(pdrng, reduced_demand[pdrng], color='k', alpha=alpha, linestyle='dashed')
    line_shedload = mlines.Line2D([], [], color='black', alpha=alpha, label='New load', linestyle='dashed')

    line_demand = mlines.Line2D([], [], color='black', label='Load')
    # plt.legend(handles=[line_demand] + patches[::-1], loc=4)

    if not load_changed and level is None:
        plt.legend(handles=[line_demand] + patches[::-1], loc=4, bbox_to_anchor=(1.2, 0.5))
    elif not load_changed:
        plt.legend(handles=[line_demand] + [line_SOC] + patches[::-1], loc=4, bbox_to_anchor=(1.2, 0.5))
    elif level is None:
        plt.legend(handles=[line_demand] + [line_shedload] + patches[::-1], loc=4, bbox_to_anchor=(1.2, 0.5))
        axes[0].fill_between(demand.index, demand, reduced_demand, facecolor="none", hatch="X", edgecolor="k",
                             linestyle='dashed')
    else:
        plt.legend(title='Dispatch for ' + demand.name[1], handles=[line_demand] + [line_shedload] + [line_SOC] +
                                                                   patches[::-1], loc=4, bbox_to_anchor=(1.2, 0.5))
        axes[0].fill_between(demand.index, demand, reduced_demand, facecolor="none", hatch="X", edgecolor="k",
                             linestyle='dashed')

    plt.subplots_adjust(right=0.8)
    plt.show()


def plot_rug(df_series, on_off=False, cmap='Greys', fig_title='', normalized=False):
    """Create multiaxis rug plot from pandas Dataframe

    Arguments:
        df_series (pd.DataFrame): 2D pandas with timed index
        on_off (bool): if True all points that are above 0 will be plotted as one color. If False all values will be colored based on their value.
        cmap (str): palette name (from colorbrewer, matplotlib etc.)
        fig_title (str): Figure title
        normalized (bool): if True, all series colormaps will be normalized based on the maximum value of the dataframe
    Returns:
        plot

    Function copied from enlopy v0.1 www.github.com/kavvkon/enlopy. Install with `pip install enlopy` for latest version.
    """

    def format_axis(iax):
        # Formatting: remove all lines (not so elegant)
        for spine in ['top', 'right', 'left', 'bottom']:
            iax.axes.spines[spine].set_visible(False)
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
    if isinstance(df_series, pd.DataFrame):
        rows = len(df_series.columns)
    elif isinstance(df_series, pd.Series):
        df_series = df_series.to_frame()
        rows = 1
    else:
        raise ValueError("Has to be either Series or Dataframe")
    if len(df_series) < 1:
        raise ValueError("Has to be non empty Series or Dataframe")

    max_color = np.nanmax(df_series.values)
    min_color = np.nanmin(df_series.values)

    __, axes = plt.subplots(nrows=rows, ncols=1, sharex=True,
                            figsize=(16, 0.25 * rows), squeeze=False,
                            frameon=True, gridspec_kw={'hspace': 0.15})

    for (item, iseries), iax in zip(df_series.iteritems(), axes.ravel()):
        format_axis(iax)
        iax.set_ylabel(str(item)[:30], rotation='horizontal',
                       rotation_mode='anchor',
                       horizontalalignment='right', x=-0.01)
        x = iseries.index

        if iseries.sum() > 0:  # if series is not empty
            if on_off:
                i_on_off = iseries.apply(flag_operation).replace(False, np.nan)
                i_on_off.plot(ax=iax, style='|', lw=.7, cmap=cmap)
            else:
                y = np.ones(len(iseries))
                # Define (truncated) colormap:
                if not normalized:  # Replace max_color (frame) with series max
                    max_color = np.nanmax(iseries.values)
                    min_color = np.nanmin(iseries.values)
                # Hack to plot max color when all series are equal
                if np.isclose(min_color, max_color):
                    min_color = min_color * 0.99

                iax.scatter(x, y,
                            marker='|', s=100,
                            c=iseries.values,
                            vmin=min_color,
                            vmax=max_color,
                            cmap=cmap)

    axes.ravel()[0].set_title(fig_title)
    axes.ravel()[-1].spines['bottom'].set_visible(True)
    axes.ravel()[-1].set_xlim(np.min(x), np.max(x))


def plot_energy_zone_fuel(inputs, results, PPindicators):
    """
    Plots the generation for each zone, disaggregated by fuel type

    :param results:         Dictionnary with the outputs of the model (output of the function GetResults)
    :param PPindicators:    Por powerplant statistics (output of the function get_indicators_powerplant)
    """
    fuels = PPindicators.Fuel.unique()
    zones = PPindicators.Zone.unique()

    GenPerZone = pd.DataFrame(index=zones, columns=fuels)
    # First make sure that all fuels are present. If not, initialize an empty series
    for f in commons['Fuels'] + ['FlowIn']:
        if f not in GenPerZone:
            GenPerZone[f] = 0
    for z in zones:
        for f in fuels:
            tmp = PPindicators[(PPindicators.Fuel == f) & (PPindicators.Zone == z)]
            GenPerZone.loc[z, f] = tmp.Generation.sum()
        NetImports = get_imports(results['OutputFlow'], z)
        if NetImports > 0:
            GenPerZone.loc[z, 'FlowIn'] = NetImports

    cols = [col for col in commons['MeritOrder'] if col in GenPerZone]
    GenPerZone = GenPerZone[cols] / 1E6
    colors = [commons['colors'][tech] for tech in GenPerZone.columns]
    ax = GenPerZone.plot(kind="bar", figsize=(12, 8), stacked=True, color=colors, alpha=0.8, legend='reverse',
                         title='Generation per zone (the horizontal lines indicate the demand)')
    ax.set_ylabel('Generation [TWh]')
    demand = inputs['param_df']['Demand']['DA'].sum() / 1E6
    ax.barh(demand, left=ax.get_xticks() - 0.4, width=[0.8] * len(demand), height=ax.get_ylim()[1] * 0.005, linewidth=2,
            color='k')
    plt.show()
    return ax


def plot_zone_capacities(inputs, plot=True):
    """
    Plots the installed capacity for each zone, disaggregated by fuel type

    :param inputs:         Dictionnary with the inputs of the model (output of the function GetResults)
    """
    units = inputs['units']
    ZoneFuels = {}
    for u in units.index:
        ZoneFuels[(units.Zone[u], units.Fuel[u])] = (units.Zone[u], units.Fuel[u])

    PowerCapacity = pd.DataFrame(columns=inputs['sets']['f'], index=inputs['sets']['n'])
    StorageCapacity = pd.DataFrame(columns=inputs['sets']['f'], index=inputs['sets']['n'])
    for n, f in ZoneFuels:
        idx = ((units.Zone == n) & (units.Fuel == f))
        PowerCapacity.loc[n, f] = (units.PowerCapacity[idx] * units.Nunits[idx]).sum()
        StorageCapacity.loc[n, f] = (units.StorageCapacity[idx] * units.Nunits[idx]).sum()

    cols = [col for col in commons['MeritOrder'] if col in PowerCapacity]
    PowerCapacity = PowerCapacity[cols]
    if plot:
        colors = [commons['colors'][tech] for tech in PowerCapacity.columns]
        ax = PowerCapacity.plot(kind="bar", figsize=(12, 8), stacked=True, color=colors, alpha=0.8, legend='reverse',
                                title='Installed capacity per zone (the horizontal lines indicate the peak demand)')
        ax.set_ylabel('Capacity [MW]')
        demand = inputs['param_df']['Demand']['DA'].max()
        ax.barh(demand, left=ax.get_xticks() - 0.4, width=[0.8] * len(demand), height=ax.get_ylim()[1] * 0.005,
                linewidth=2,
                color='k')
    return {'PowerCapacity': PowerCapacity, 'StorageCapacity': StorageCapacity}


def plot_zone(inputs, results, z='', rng=None, rug_plot=True):
    """
    Generates plots from the dispa-SET results for one specific zone

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :param z:           Considered zone (e.g. 'BE')
    """
    if z == '':
        Nzones = len(inputs['sets']['n'])
        z = inputs['sets']['n'][np.random.randint(Nzones)]
        print('Randomly selected zone for the detailed analysis: ' + z)
    elif z not in inputs['sets']['n']:
        logging.critical('Zone ' + z + ' is not in the results')
        Nzones = len(inputs['sets']['n'])
        z = inputs['sets']['n'][np.random.randint(Nzones)]
        logging.critical('Randomly selected zone: ' + z)

    plotdata = get_plot_data(inputs, results, z) / 1000  # GW

    aggregation = False
    if 'OutputStorageLevel' in results:
        lev = filter_by_zone(results['OutputStorageLevel'], inputs, z)
        lev = lev * inputs['units']['StorageCapacity'].loc[lev.columns] / 1e3  # GWh of storage
        level = filter_by_storage(lev, inputs, StorageSubset='s')
        levels = pd.DataFrame(index=results['OutputStorageLevel'].index, columns=inputs['sets']['t'])
        for t in ['HDAM', 'HPHS', 'BEVS', 'BATS']:
            temp = filter_by_tech(level, inputs, t)
            levels[t] = temp.sum(axis=1)
        levels.dropna(axis=1, inplace=True)
        for col in levels.columns:
            if levels[col].max() == 0 and levels[col].min() == 0:
                del levels[col]
        if aggregation is True:
            level = level.sum(axis=1)
        else:
            level = levels
    else:
        level = pd.Series(0, index=results['OutputPower'].index)

    if 'OutputPowerConsumption' in results:
        demand_p2h = filter_by_zone(results['OutputPowerConsumption'], inputs, z) / 1000  # GW
        demand_p2h = demand_p2h.sum(axis=1)
    else:
        demand_p2h = pd.Series(0, index=results['OutputPower'].index)
    if ('Flex', z) in inputs['param_df']['Demand']:
        demand_flex = inputs['param_df']['Demand'][('Flex', z)] / 1000
    else:
        demand_flex = pd.Series(0, index=results['OutputPower'].index)

    demand_da = inputs['param_df']['Demand'][('DA', z)] / 1000  # GW
    demand = pd.DataFrame(demand_da + demand_p2h + demand_flex, columns=[('DA', z)])
    demand = demand[('DA', z)]

    sum_generation = plotdata.sum(axis=1)
    # if 'OutputShedLoad' in results:
    if 'OutputShedLoad' in results and z in results['OutputShedLoad']:
        shed_load = results['OutputShedLoad'][z] / 1000  # GW
        shed_load = pd.Series(shed_load, index=demand.index).fillna(0)
    else:
        shed_load = pd.Series(0, index=demand.index) / 1000  # GW
    if 'OutputDemandModulation' in results and z in results['OutputDemandModulation']:
        shifted_load = -results['OutputDemandModulation'][z] / 1000  # GW
        shifted_load = pd.Series(shifted_load, index=demand.index).fillna(0)
    else:
        shifted_load = pd.Series(0, index=demand.index) / 1000  # GW
    diff = (sum_generation - demand + shifted_load + shed_load).abs()

    if diff.max() > 0.01 * demand.max():
        logging.critical('There is up to ' + str(
            diff.max() / demand.max() * 100) + '% difference in the instantaneous energy balance of zone ' + z)

    if 'OutputCurtailedPower' in results and z in results['OutputCurtailedPower']:
        curtailment = results['OutputCurtailedPower'][z] / 1000  # GW
    else:
        curtailment = None

    plot_dispatch(demand, plotdata, level, curtailment=curtailment, shedload=shed_load, shiftedload=shifted_load,
                  rng=rng, alpha=0.5)

    # Generation plot:
    if rug_plot:
        ZoneGeneration = filter_by_zone(results['OutputPower'], inputs, z)
        try:
            import enlopy as el  # try to get latest version
            el.plot_rug(ZoneGeneration, on_off=False, cmap='gist_heat_r', fig_title=z)
        except ImportError:
            plot_rug(ZoneGeneration, on_off=False, cmap='gist_heat_r', fig_title=z)

    return True


def storage_levels(inputs, results):
    """
    Reads the DispaSET results and provides the difference between the minimum storage profile and the computed storage profile

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    """
    isstorage = pd.Series(index=inputs['units'].index)
    for u in isstorage.index:
        isstorage[u] = inputs['units'].Technology[u] in commons['tech_storage']
    sto_units = inputs['units'][isstorage]
    results['OutputStorageLevel'].plot(figsize=(12, 6), title='Storage levels')
    StorageError = ((inputs['param_df']['StorageProfile'] * sto_units.StorageCapacity).subtract(
        results['OutputStorageLevel'], 'columns')).divide((sto_units.StorageCapacity), 'columns') * (-100)
    StorageError = StorageError.dropna(1)
    if not StorageError.empty:
        ax = StorageError.plot(figsize=(12, 6),
                               title='Difference between the calculated storage Levels and the (imposed) minimum level')
        ax.set_ylabel('%')

    return True
