# -*- coding: utf-8 -*-

import logging
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from .postprocessing import get_imports, get_plot_data, filter_by_zone, filter_by_tech, filter_by_storage, \
    get_power_flow_tracing, get_EFOH, get_heat_plot_data, filter_by_heating_zone, filter_by_tech_list
from ..common import commons


def plot_dispatch(demand, plotdata, y_ax='', level=None, curtailment=None, shedload=None, shiftedload=None, rng=None,
                  alpha=None, figsize=(13, 7), ntc=None):
    """
    Function that plots the dispatch data and the reservoir level as a cumulative sum

    :param demand:      Pandas Series with the demand curve
    :param plotdata:    Pandas Dataframe with the data to be plotted. Negative columns should be at the beginning.
                        Output of the function GetPlotData
    :param y_ax:        Y axis label
    :param level:       Optional pandas series/dataframe with an (dis)aggregated reservoir level for the considered zone
    :param curtailment: Optional pandas series with the value of the curtailment
    :param shedload:    Optional pandas series with the value of the shed load
    :param shiftedload: Optional pandas series with the value of the shifted load
    :param rng:         Indexes of the values to be plotted. If undefined, the first week is plotted
    :param alpha:       Alpha value for colours
    :param figsize:     Figure size in inch
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

    # NTC plot data
    if (ntc is not None) and ('FlowIn' in plotdata) and ('FlowOut' in plotdata):
        ntc['FlowIn'], ntc['FlowOut'], ntc['ZeroLine'] = plotdata['FlowIn'], plotdata['FlowOut'], 0

    # Netting the interconnections:
    if ('FlowIn' in plotdata) and ('FlowOut' in plotdata):
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

    if ntc is not None:
        fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=figsize, frameon=True,  # 14 4*2
                                 gridspec_kw={'height_ratios': [2.7, .8, .8], 'hspace': 0.04})

        axes[2].plot(pdrng, ntc.loc[pdrng, 'NTCIn'].values, color='r')
        axes[2].plot(pdrng, ntc.loc[pdrng, 'NTCOut'].values, color='g')
        axes[2].set_xlim(pdrng[0], pdrng[-1])
        axes[2].fill_between(pdrng, ntc.loc[pdrng, 'NTCIn'], ntc.loc[pdrng, 'ZeroLine'], facecolor='red',
                             alpha=alpha, hatch=commons['hatches']['FlowIn'])
        axes[2].fill_between(pdrng, ntc.loc[pdrng, 'ZeroLine'], ntc.loc[pdrng, 'NTCOut'], facecolor='green',
                             alpha=alpha, hatch=commons['hatches']['FlowOut'])
        axes[2].set_ylabel('NTC [GWh]')
    else:
        fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=figsize, frameon=True,  # 14 4*2
                                 gridspec_kw={'height_ratios': [2.7, .8], 'hspace': 0.04})

    # Create left axis:
    axes[0].plot(pdrng, demand[pdrng], color='k')
    axes[0].set_xlim(pdrng[0], pdrng[-1])

    fig.suptitle(y_ax + ' dispatch for zone ' + demand.name)

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

    axes[0].set_ylabel(y_ax + ' [GW]')
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
        plt.legend(title='Dispatch for ' + demand.name, handles=[line_demand] + [line_shedload] + [line_SOC] +
                   patches[::-1], loc=4, bbox_to_anchor=(1.2, 0.5))
        axes[0].fill_between(demand.index, demand, reduced_demand, facecolor="none", hatch="X", edgecolor="k",
                             linestyle='dashed')

    plt.subplots_adjust(right=0.8)
    plt.show()


def plot_rug(df_series, on_off=False, cmap='Greys', fig_title='', normalized=False):
    """
    Create multiaxis rug plot from pandas Dataframe

    :param df_series:    2D pandas with timed index (pd.DataFrame)
    :param on_off:       if True all points that are above 0 will be plotted as one color. If False all
                         values will be colored based on their value. (bool)
    :param cmap:         palette name (from colorbrewer, matplotlib etc.) (str)
    :param fig_title:    Figure title (str)
    :param normalized:   if True, all series colormaps will be normalized based on the maximum value
                         of the dataframe (bool)
    :return plot:

    Function copied from enlopy v0.1 www.github.com/kavvkon/enlopy. Install with `pip install enlopy` for latest version
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
    plt.show()


def plot_energy_zone_fuel(inputs, results, PPindicators):
    """
    Plots the generation for each zone, disaggregated by fuel type

    :param inputs:          Dictionary with the inputs of the model
    :param results:         Dictionary with the outputs of the model (output of the function GetResults)
    :param PPindicators:    Por powerplant statistics (output of the function get_indicators_powerplant)
    """
    fuels = PPindicators[[u in [x for x in commons['Technologies'] if x not in commons['tech_heat'] +
                                commons['tech_p2ht'] + commons['tech_thermal_storage']]
                          for u in PPindicators['Technology']]].Fuel.unique()
    zones = PPindicators.Zone.unique()
    GenPerZone = pd.DataFrame(index=zones, columns=fuels)

    # First make sure that all fuels are present. If not, initialize an empty series
    for f in commons['Fuels'] + ['FlowIn']:
        if f not in GenPerZone:
            GenPerZone[f] = 0

    # Assign generation
    power_consumption = pd.DataFrame(index=results['OutputPowerConsumption'].index)
    storage_input = pd.DataFrame(index=results['OutputPower'].index)
    for z in zones:
        for f in fuels:
            tmp = PPindicators[(PPindicators.Fuel == f) & (PPindicators.Zone == z)]
            GenPerZone.loc[z, f] = tmp.Generation.sum()
        NetImports = get_imports(results['OutputFlow'], z)
        if NetImports > 0:
            GenPerZone.loc[z, 'FlowIn'] = NetImports
        if filter_by_zone(results['OutputPowerConsumption'], inputs, z).empty:
            power_consumption.loc[:, z] = 0
        else:
            power_consumption.loc[:, z] = filter_by_zone(results['OutputPowerConsumption'], inputs, z).sum(axis=1)
        if filter_by_zone(filter_by_tech_list(results['OutputStorageInput'], inputs, commons['tech_storage']), inputs,
                          z).empty:
            storage_input.loc[:, z] = 0
        else:
            storage_input.loc[:, z] = filter_by_zone(filter_by_tech_list(results['OutputStorageInput'], inputs,
                                                                         commons['tech_storage']),
                                                     inputs, z).sum(axis=1)

    cols = [col for col in commons['MeritOrder'] if col in GenPerZone]
    GenPerZone = GenPerZone[cols] / 1E6
    GenPerZone.sort_index(inplace=True)
    colors = [commons['colors'][tech] for tech in GenPerZone.columns]
    ax = GenPerZone.plot(kind="bar", figsize=(12, 8), stacked=True, color=colors, alpha=0.8, legend='reverse',
                         title='Generation per zone (the horizontal lines indicate the demand [black], storage input '
                               '[green] & power consumption [blue])')
    ax.set_ylabel('Generation [TWh]')
    demand = inputs['param_df']['Demand']['DA'].sum() / 1E6
    demand.sort_index(inplace=True)
    power_consumption = power_consumption.sum() / 1E6
    power_consumption.sort_index(inplace=True)
    storage_input = storage_input.sum() / 1e6
    storage_input.sort_index(inplace=True)

    # Needs to be plotted in this order so that black line is always on top
    ax.barh(power_consumption + storage_input + demand, left=ax.get_xticks() - 0.4,
            width=[0.8] * len(power_consumption),
            height=ax.get_ylim()[1] * 0.005, linewidth=2, color='b')
    ax.barh(storage_input + demand, left=ax.get_xticks() - 0.4, width=[0.8] * len(storage_input),
            height=ax.get_ylim()[1] * 0.005, linewidth=2, color='g')
    ax.barh(demand, left=ax.get_xticks() - 0.4, width=[0.8] * len(demand), height=ax.get_ylim()[1] * 0.005, linewidth=2,
            color='k')

    ZonePosition = GenPerZone.copy()
    ZonePosition['ShedLoad'] = results['OutputShedLoad'].sum() / 1E6
    ShedLoad = ZonePosition.loc[:, 'ShedLoad']
    ax.bar(range(0, ShedLoad.index.size), ShedLoad.values, bottom=GenPerZone.sum(axis=1).values,
           edgecolor='black', hatch='X', color='w', width=[0.4])

    curtailed_power = -results['OutputCurtailedPower'].sum() / 1E6
    curtailed_power.sort_index(inplace=True)
    ax.bar(range(0, curtailed_power.index.size), curtailed_power.values, color='r', width=[0.5])
    ax.plot([-1, curtailed_power.index.size + 1], [0, 0], color='k')

    handles, labels = ax.get_legend_handles_labels()  # get the handles
    ax.legend(reversed(handles), reversed(labels), loc=4, bbox_to_anchor=(1.15, 0.25))
    plt.show()

    # Generation share plot
    GenPerZone_prct = GenPerZone.div(ZonePosition.sum(axis=1), axis=0)
    ShedLoad = ZonePosition.loc[:, 'ShedLoad'] / ZonePosition.sum(axis=1)
    demand_prct = demand / ZonePosition.sum(axis=1)
    power_consumption_prct = (demand + power_consumption) / ZonePosition.sum(axis=1)
    colors2 = [commons['colors'][tech] for tech in GenPerZone_prct.columns]

    ax2 = GenPerZone_prct.plot(kind="bar", figsize=(12, 8), stacked=True, color=colors2, alpha=0.8, legend='reverse',
                               title='Generation share per zone (the horizontal lines indicate the demand)')
    ax2.set_ylabel('Generation share [%]')
    ax2.bar(range(0, ShedLoad.index.size), ShedLoad.values, bottom=GenPerZone_prct.sum(axis=1).values,
            edgecolor='black', hatch='X', color='w', width=[0.4])
    ax2.barh(power_consumption_prct, left=ax2.get_xticks() - 0.4, width=[0.8] * len(power_consumption_prct),
             height=ax2.get_ylim()[1] * 0.005, linewidth=2, color='b')
    ax2.barh(demand_prct, left=ax2.get_xticks() - 0.4, width=[0.8] * len(demand_prct), height=ax2.get_ylim()[1] * 0.005,
             linewidth=2, color='k')
    handles, labels = ax2.get_legend_handles_labels()  # get the handles
    ax2.legend(reversed(handles), reversed(labels), loc=4, bbox_to_anchor=(1.14, 0.25))
    plt.show()

    # Heat generation bar chart
    fuels_heat = PPindicators[[u in [x for x in commons['Technologies']]
                               for u in PPindicators['Technology']]].Fuel.unique()
    zones_heat = PPindicators['Zone_th'].unique()
    zones_heat = np.array([i for i in zones_heat if i])

    if zones_heat.size != 0 and fuels_heat.size != 0:
        HeatPerZone_th = pd.DataFrame(index=zones_heat, columns=fuels_heat)
        for fuel in commons['Fuels']:
            if fuel not in HeatPerZone_th:
                HeatPerZone_th[fuel] = 0

        heat_storage_input = pd.DataFrame(index=results['OutputStorageInput'].index)
        for z in zones_heat:
            for f in fuels_heat:
                tmp = PPindicators[(PPindicators.Fuel == f) & (PPindicators.Zone_th == z)]
                HeatPerZone_th.loc[z, f] = tmp.HeatGeneration.sum()
            # HeatSlack = results['OutputHeatSlack'].loc[:,z].sum()
            # if HeatSlack > 0:
            #     HeatPerZone_th.loc[z, 'HeatSlack'] = HeatSlack
            z = str(z)
            if filter_by_heating_zone(filter_by_tech_list(results['OutputStorageInput'], inputs,
                                                          commons['tech_thermal_storage']), inputs, z).empty:
                heat_storage_input.loc[:, z] = 0
            else:
                heat_storage_input.loc[:, z] = filter_by_heating_zone(
                    filter_by_tech_list(results['OutputStorageInput'], inputs,
                                        commons['tech_thermal_storage']), inputs, z).sum(axis=1)

        heat_storage_input = heat_storage_input.sum() / 1e6
        heat_storage_input.sort_index(inplace=True)

        cols_heat = [col for col in commons['MeritOrderHeat'] if col in HeatPerZone_th]
        HeatPerZone_th = HeatPerZone_th[cols_heat] / 1e6
        HeatPerZone_th.sort_index(inplace=True)
        colors = [commons['colors'][tech] for tech in HeatPerZone_th.columns]
        ax = HeatPerZone_th.plot(kind="bar", figsize=(12, 8), stacked=True, color=colors, alpha=0.8, legend='reverse',
                                 title='Heat generation per heating zone (the horizontal lines indicate the demand '
                                       '[black] & storage input [green])')
        ax.set_ylabel('Generation [TWh]')

        # Check for heat slack and plot it
        Slack = results['OutputHeatSlack'].sum() / 1e6
        Slack = Slack.reindex(inputs['sets']['n_th']).fillna(0)
        Slack.sort_index(inplace=True)
        ax.bar(range(0, Slack.index.size), Slack.values, bottom=HeatPerZone_th.sum(axis=1).values,
               edgecolor='black', hatch='X', color='#943126ff', width=[0.4])

        # Identify heat demand and plot it
        demand = inputs['param_df']['HeatDemand'].sum() / 1E6
        demand.sort_index(inplace=True)
        ax.barh(demand + heat_storage_input, left=ax.get_xticks() - 0.4, width=[0.8] * len(heat_storage_input),
                height=ax.get_ylim()[1] * 0.005, linewidth=2, color='g')
        ax.barh(demand, left=ax.get_xticks() - 0.4, width=[0.8] * len(demand), height=ax.get_ylim()[1] * 0.005,
                linewidth=2, color='k')

        handles, labels = ax.get_legend_handles_labels()  # get the handles
        handles.append(Patch(facecolor='#943126ff', hatch='X'))
        labels.append('HeatSlack')
        ax.legend(reversed(handles), reversed(labels), loc=4, bbox_to_anchor=(1.14, 0.25))
        plt.show()

        return GenPerZone, HeatPerZone_th
    else:
        return GenPerZone


def plot_zone_capacities(inputs, results, plot=True):
    """
    Plots the installed capacity for each zone, disaggregated by fuel type

    :param inputs:          Dictionary with the inputs of the model
    :param results:         Dictionary with the results of the model (output of the function GetResults)
    :param plot:            Bool param to turn the plotting on/off
    """

    units = inputs['units'][[u in [x for x in commons['Technologies'] if x not in commons['tech_heat'] +
                                   commons['tech_p2ht']] for u in inputs['units']['Technology']]]
    units_heat = inputs['units'][[u in [x for x in commons['Technologies'] if x in commons['tech_heat'] +
                                        commons['tech_p2ht'] + commons['tech_thermal_storage']]
                                  for u in inputs['units']['Technology']]]
    units_chp = inputs['units'][[x.lower() in commons['types_CHP'] for x in inputs['units']['CHPType']]]

    power_consumption = pd.DataFrame(index=results['OutputPowerConsumption'].index)
    for z in inputs['sets']['n']:
        if filter_by_zone(results['OutputPowerConsumption'], inputs, z).empty:
            power_consumption.loc[:, z] = 0
        else:
            power_consumption.loc[:, z] = filter_by_zone(results['OutputPowerConsumption'], inputs, z).sum(axis=1)

    for u in units_chp.index:
        PowerCapacity = units_chp.loc[u, 'PowerCapacity']

        if units_chp.loc[u, 'CHPType'].lower() == 'p2h':
            PurePowerCapacity = PowerCapacity
        else:
            # If maximum heat is not defined, then it is defined as the intersection between two lines
            if pd.isnull(units_chp.loc[u, 'CHPMaxHeat']):
                MaxHeat = PowerCapacity / units_chp.loc[u, 'CHPPowerToHeat']
                units_chp.loc[u, 'CHPMaxHeat'] = 'inf'
            else:
                MaxHeat = units_chp.loc[u, 'CHPMaxHeat']
            PurePowerCapacity = PowerCapacity + units_chp.loc[u, 'CHPPowerLossFactor'] * MaxHeat
        units.loc[u, 'PowerCapacity'] = PurePowerCapacity
        units_chp.loc[u, 'PowerCapacity'] = MaxHeat

    # units_heat = units_heat.append(units_chp)
    units_heat = pd.concat([units_heat, units_chp])

    ZoneFuels = {}
    ZoneFuelsHT = {}
    for u in units.index:
        ZoneFuels[(units.Zone[u], units.Fuel[u])] = (units.Zone[u], units.Fuel[u])
    for u in units_heat.index:
        ZoneFuelsHT[(units_heat.Zone_th[u], units_heat.Fuel[u])] = (units_heat.Zone_th[u], units_heat.Fuel[u])

    PowerCapacity = pd.DataFrame(columns=inputs['sets']['f'], index=inputs['sets']['n'])
    StorageCapacity = pd.DataFrame(columns=inputs['sets']['f'], index=inputs['sets']['n'])
    HeatCapacity = pd.DataFrame(columns=inputs['sets']['f'], index=inputs['sets']['n_th'])
    for n, f in ZoneFuels:
        idx = ((units.Zone == n) & (units.Fuel == f))
        PowerCapacity.loc[n, f] = (units.PowerCapacity[idx] * units.Nunits[idx]).sum()
        StorageCapacity.loc[n, f] = (units.StorageCapacity[idx] * units.Nunits[idx]).sum()

    for n_th, f in ZoneFuelsHT:
        idx = ((units_heat.Zone_th == n_th) & (units_heat.Fuel == f))
        units_heat.COP.fillna(1, inplace=True)
        HeatCapacity.loc[n_th, f] = (units_heat.PowerCapacity[idx] * units_heat.Nunits[idx] * units_heat.COP[idx]).sum()

    cols = [col for col in commons['MeritOrder'] if col in PowerCapacity]
    PowerCapacity = PowerCapacity[cols]
    PowerCapacity.sort_index(inplace=True)
    if plot:
        colors = [commons['colors'][tech] for tech in PowerCapacity.columns]
        ax = PowerCapacity.plot(kind="bar", figsize=(12, 8), stacked=True, color=colors, alpha=0.8, legend='reverse',
                                title='Installed capacity per zone (the horizontal lines indicate the peak demand '
                                      '[black] & peak power consumption [blue])')
        ax.set_ylabel('Capacity [MW]')
        demand = inputs['param_df']['Demand']['DA'].max()
        demand.sort_index(inplace=True)
        power_consumption = power_consumption.max()
        power_consumption.sort_index(inplace=True)

        # Needs to be plotted in this order so that black line is always on top
        ax.barh(power_consumption + demand, left=ax.get_xticks() - 0.4, width=[0.8] * len(power_consumption),
                height=ax.get_ylim()[1] * 0.005, linewidth=2, color='b')
        ax.barh(demand, left=ax.get_xticks() - 0.4, width=[0.8] * len(demand), height=ax.get_ylim()[1] * 0.005,
                linewidth=2, color='k')

        # ax.barh(demand, left=ax.get_xticks() - 0.4, width=[0.8] * len(demand), height=ax.get_ylim()[1] * 0.005,
        #         linewidth=2, color='k')
    plt.show()

    cols_heat = [col for col in commons['MeritOrderHeat'] if col in HeatCapacity]
    HeatCapacity = HeatCapacity[cols_heat]
    HeatCapacity.sort_index(inplace=True)
    if len(HeatCapacity) > 0 and plot:
        colors = [commons['colors'][tech] for tech in HeatCapacity.columns]
        ax = HeatCapacity.plot(kind="bar", figsize=(12, 8), stacked=True, color=colors, alpha=0.8, legend='reverse',
                               title='Installed heat capacity per zone ' +
                                     '(the horizontal lines indicate the peak heating demand)')
        ax.set_ylabel('Capacity [MW]')
        demand = inputs['param_df']['HeatDemand'].max()
        demand.sort_index(inplace=True)
        ax.barh(demand, left=ax.get_xticks() - 0.4, width=[0.8] * len(demand), height=ax.get_ylim()[1] * 0.005,
                linewidth=2,
                color='k')
        plt.show()

    return {'PowerCapacity': PowerCapacity, 'StorageCapacity': StorageCapacity, 'HeatCapacity': HeatCapacity}


def plot_zone(inputs, results, z='', z_th=None, rng=None, rug_plot=True):
    """
    Generates plots from the dispa-SET results for one specific zone

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    :param z:           Considered zone (e.g. 'BE')
    :param z_th:        Considered thermal zone (e.g. 'BE_th')
    :param rng:         Date range to be considered in the plot
    :param rug_plot:    Rug plot on/off
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
        for col in lev.columns:
            if 'BEVS' in col:
                lev[col] = lev[col] * inputs['param_df']['AvailabilityFactor'][col]
        level = filter_by_storage(lev, inputs, StorageSubset='s')
        levels = pd.DataFrame(index=results['OutputStorageLevel'].index, columns=inputs['sets']['t'])
        for t in commons['tech_storage']:
            temp = filter_by_tech(level, inputs, t)
            levels[t] = temp.sum(axis=1)
        levels.dropna(axis=1, inplace=True)
        for col in levels.columns:
            if levels[col].max() == 0 and levels[col].min() == 0:
                del levels[col]

        lev_heat = filter_by_heating_zone(results['OutputStorageLevel'], inputs, z_th)
        lev_heat = lev_heat * inputs['units']['StorageCapacity'].loc[lev_heat.columns] / 1e3  # GWh of storage
        # Filter storage levels for thermal storage
        level_heat = filter_by_storage(lev_heat, inputs, StorageSubset='thms')
        levels_heat = pd.DataFrame(index=results['OutputStorageLevel'].index, columns=inputs['sets']['t'])
        for t in commons['tech_thermal_storage']:
            temp = filter_by_tech(level_heat, inputs, t)
            levels_heat[t] = temp.sum(axis=1)
        levels_heat.dropna(axis=1, inplace=True)

        for col in levels_heat.columns:
            if levels_heat[col].max() == 0 and levels_heat[col].min() == 0:
                del levels_heat[col]

        if aggregation is True:
            level = level.sum(axis=1)
            level_heat = level_heat.sum(axis=1)
        else:
            level = levels
            level_heat = levels_heat
    else:
        level = pd.Series(0, index=results['OutputPower'].index)
        level_heat = pd.Series(0, index=results['OutputPower'].index)

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

    # Assign NTC
    ntc = pd.DataFrame(0, columns=['NTCIn', 'NTCOut'], index=plotdata.index)
    for col in inputs['param_df']['FlowMaximum']:
        from_node, to_node = col.split('->')
        if to_node.strip() == z:
            ntc['NTCIn'] = ntc['NTCIn'] + inputs['param_df']['FlowMaximum'][col]
        if from_node.strip() == z:
            ntc['NTCOut'] = ntc['NTCOut'] - inputs['param_df']['FlowMaximum'][col]
    ntc = ntc / 1e3  # GW
    ntc = ntc.loc[:, (ntc != 0).any(axis=0)]

    # Plot power dispatch
    demand.rename(z, inplace=True)
    if ntc.empty:
        plot_dispatch(demand, plotdata, y_ax='Power', level=level, curtailment=curtailment, shedload=shed_load,
                      shiftedload=shifted_load, rng=rng, alpha=0.5)
    else:
        plot_dispatch(demand, plotdata, y_ax='Power', level=level, curtailment=curtailment, shedload=shed_load,
                      shiftedload=shifted_load, ntc=ntc, rng=rng, alpha=0.5)

    # Plot heat dispatch
    Nzones_th = len(inputs['sets']['n_th'])
    if Nzones_th > 0 and (z_th is None or (z_th not in inputs['sets']['n_th'])):
        z_th = inputs['sets']['n_th'][np.random.randint(Nzones_th)]
        logging.info('Randomly selected thermal zone: ' + z_th)

    if Nzones_th > 0:
        heat_demand = inputs['param_df']['HeatDemand'][z_th] / 1000
        heat_demand = pd.DataFrame(heat_demand, columns=[z_th])
        heat_demand = heat_demand[z_th]
        heat_plotdata = get_heat_plot_data(inputs, results, z_th) / 1000

        if 'OutputCurtailedHeat' in results and z_th in results['OutputCurtailedHeat']:
            heat_curtailment = results['OutputCurtailedHeat'][z_th] / 1000  # GW
        else:
            heat_curtailment = None

        plot_dispatch(heat_demand, heat_plotdata, y_ax='Heat', level=level_heat, curtailment=heat_curtailment, rng=rng,
                      alpha=0.5)

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
    Reads the DispaSET results and provides the difference between the minimum storage profile and the computed storage
    profile

    :param inputs:      DispaSET inputs
    :param results:     DispaSET results
    """
    isstorage = pd.Series(index=inputs['units'].index)
    for u in isstorage.index:
        isstorage[u] = inputs['units'].Technology[u] in commons['tech_storage']
    sto_units = inputs['units'][isstorage]
    results['OutputStorageLevel'].plot(figsize=(12, 6), title='Storage levels')
    StorageError = ((inputs['param_df']['StorageProfile'] * sto_units.StorageCapacity).subtract(
        results['OutputStorageLevel'], 'columns')).divide(sto_units.StorageCapacity, 'columns') * (-100)
    StorageError = StorageError.dropna(1)
    if not StorageError.empty:
        ax = StorageError.plot(figsize=(12, 6),
                               title='Difference between the calculated storage Levels and the (imposed) minimum level')
        ax.set_ylabel('%')

    return True


# TODO this function should be generalized
def plot_storage_levels(inputs, results, z):
    """
    This function plots the reservoir levels profiles during the year
    #TODO: Check this function and decide if necessary
    :param inputs:      Dispa-SET inputs
    :param results:     Dispa-SET results
    :param z:           zone considered
    """
    BATS_unit = [u for u in inputs['units'].index if
                 (inputs['units'].loc[u, 'Zone'] == z and inputs['units'].loc[u, 'Technology'] == 'BATS')]
    H2_unit = [u for u in inputs['units'].index if
               inputs['units'].loc[u, 'Zone'] == z and inputs['units'].loc[u, 'Technology'] == 'P2GS']
    HDAM_unit = [u for u in inputs['units'].index if
                 inputs['units'].loc[u, 'Zone'] == z and inputs['units'].loc[u, 'Technology'] == 'HDAM']

    BATS_level = results['OutputStorageLevel'][BATS_unit]
    H2_level = results['OutputStorageLevel'][H2_unit]
    HDAM_level = results['OutputStorageLevel'][HDAM_unit]

    # plot
    index_plot = np.arange(len(BATS_level.index))
    plt.rcParams.update({'font.size': 10})
    fig, (ax, ax1, ax2) = plt.subplots(3, 1, sharex=True)

    ax.fill_between(index_plot, 0, BATS_level.iloc[:, 0], color='#41A317ff')
    ax.set_ylabel('SOC - BATS [%]')
    ax.set_xlim(0, len(index_plot))

    ax1.fill_between(index_plot, 0, H2_level.iloc[:, 0], color='#A0522D')
    ax1.set_ylabel('SOC - H2[%]')

    ax2.fill_between(index_plot, 0, HDAM_level.iloc[:, 0], color='#00a0e1ff')
    ax2.set_ylabel('SOC - HDAM [%]')
    month_num = [360, 1056, 1800, 2544, 3264, 4008, 4728, 5472, 6216, 6936, 7680, 8400]
    month_name = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul',
                  'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    ax2.set_xticks(month_num)
    ax2.set_xticklabels(month_name)
    plt.show()
    return True


def plot_EFOH(inputs, results):
    """
    This function plots the equivalent full load operating hours of the electrolysers, together with the marginal price
    #TODO: Check this function and decide if necessary
    :param inputs:      Dispa-SET inputs
    :param results:     Dispa-SET results
    """
    EFOH = {}
    for i, s in enumerate(list(inputs.keys())): # TODO something with scenarios should be dropped
        EFOH[i] = get_EFOH(inputs[s], results[s])
        EFOH[i] = EFOH[i].sort_values(by=['EFOH'], ascending=False)

    ind = np.arange(len(EFOH[0].index))
    labels = list(EFOH[0].index)
    for count, i in enumerate(labels):
        labels[count] = i.split()[2].split("_")[0]

    x = np.arange(len(labels))  # the label locations
    width = 0.6 / len(list(inputs.keys()))  # the width of the bars

    fig, ax = plt.subplots(figsize=(16, 6))
    plt.rcParams.update({'font.size': 15})
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(15)
    if len(list(inputs.keys())) == 2:
        ax.bar(x - width / 2, EFOH[0].iloc[:, 0], width, label=list(inputs.keys())[0], color=(1, 0.5, 0, 1))
        ax.bar(x + width / 2, EFOH[1].iloc[:, 0], width, label=list(inputs.keys())[1], color=(0.192, 0.549, 0.905))
    elif len(list(inputs.keys())) == 3:
        ax.bar(x - 0.2, EFOH[0].iloc[:, 0], width, label=list(inputs.keys())[0], color=(1, 0.5, 0, 1))
        ax.bar(x, EFOH[1].iloc[:, 0], width, label=list(inputs.keys())[1], color=(0.192, 0.549, 0.905))
        ax.bar(x + 0.2, EFOH[2].iloc[:, 0], width, label=list(inputs.keys())[2], color=(0.192, 0.549, 0.3))

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('EFOH [h]')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

    return True


def plot_ElyserCap_vs_Utilization(inputs, results):
    # TODO: Check this function and decide if necesarry
    Data = pd.DataFrame(index=inputs['param_df']['sets']['n'], columns=['Cap', 'Max'])
    FC = pd.DataFrame(index=inputs['param_df']['sets']['n'], columns=['FC'])
    Sto = pd.DataFrame(index=inputs['param_df']['sets']['n'], columns=['Sto'])
    for u in inputs['param_df']['sets']['p2h2']:
        for z in inputs['param_df']['sets']['n']:
            if inputs['param_df']['Location'].loc[z, u]:
                Data.loc[z, 'Cap'] = inputs['param_df']['StorageChargingCapacity'].loc[
                                         u, 'StorageChargingCapacity'] / 1e3
                Data.loc[z, 'Max'] = results['OutputStorageInput'].loc[:, u].max() / 1e3
                FC.loc[z, 'FC'] = inputs['param_df']['PowerCapacity'].loc[u, 'PowerCapacity'] / 1e3
                Sto.loc[z, 'Sto'] = inputs['param_df']['StorageCapacity'].loc[u, 'StorageCapacity'] / 1e3

    labels = Data.index
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars
    ax = Data.loc[:, 'Cap'].plot(kind="bar", figsize=(12, 8), stacked=True,
                                 alpha=0.8, fontsize='medium')
    ax.set_ylabel('Capacity [GW]')
    ax.barh(Data.loc[:, 'Max'], left=ax.get_xticks() - 0.4, width=0.8, height=ax.get_ylim()[1] * 0.005,
            linewidth=2,
            color='k')
    ax2 = FC.loc[:, 'FC'].plot(kind="bar", figsize=(12, 8), stacked=True,
                               alpha=0.8, fontsize='medium')
    ax2.set_ylabel('Capacity [GW]')

    ax3 = Sto.loc[:, 'Sto'].plot(kind="bar", figsize=(12, 8), stacked=True,
                                 alpha=0.8, fontsize='medium')
    ax3.set_ylabel('Capacity [GWh]')

    plt.show()
    return True


def plot_H2_and_demand(inputs, results):
    # TODO: Check this function adn decide if necessary
    """
    This function plots the demand and the electrolyser demand as bar chart
    
    :param inputs:      Dispa-SET inputs
    :param results:     Dispa-SET results
    """

    # Get demand
    def get_demand(inputs, results):
        Demand = pd.DataFrame(index=['0'], columns=inputs['sets']['n'])
        for z in inputs['sets']['n']:
            plotdata = get_plot_data(inputs, results, z) / 1000  # GW

            aggregation = False
            if 'OutputStorageLevel' in results:
                lev = filter_by_zone(results['OutputStorageLevel'], inputs, z)
                lev = lev * inputs['units']['StorageCapacity'].loc[lev.columns] / 1e3  # GWh of storage
                for col in lev.columns:
                    if 'BEVS' in col:
                        lev[col] = lev[col] * inputs['param_df']['AvailabilityFactor'][col]
                level = filter_by_storage(lev, inputs, StorageSubset='s')
                levels = pd.DataFrame(index=results['OutputStorageLevel'].index, columns=inputs['sets']['t'])
                for t in ['HDAM', 'HPHS', 'BEVS', 'BATS', 'SCSP', 'P2GS']:
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
                demand_p2h = filter_by_zone(results['OutputPowerConsumption'], inputs, z) / 1e6  # TWh
                demand_p2h = demand_p2h.sum(axis=1)
            else:
                demand_p2h = pd.Series(0, index=results['OutputPower'].index)
            if ('Flex', z) in inputs['param_df']['Demand']:
                demand_flex = inputs['param_df']['Demand'][('Flex', z)] / 1e6
            else:
                demand_flex = pd.Series(0, index=results['OutputPower'].index)

            demand_da = inputs['param_df']['Demand'][('DA', z)] / 1e6  # TWh
            demand = pd.DataFrame(demand_da + demand_p2h + demand_flex, columns=[('DA', z)])
            demand = demand[('DA', z)]
            Demand[z] = demand.sum()
            # Get elyser consumption
        Elyser_consumption = pd.DataFrame(index=['0'], columns=inputs['sets']['n'])
        for u in results['OutputStorageInput'].columns:
            if 'P2GS' in u:
                c = u.split()[2].split('_')[0]
                Elyser_consumption[c] = results['OutputStorageInput'][u].sum() / 1e6  # TWh
                Elyser_consumption[c].fillna(0)

        return Demand, Elyser_consumption

    Demand = {}
    Elyser_consumption = {}
    scenarios = list(inputs.keys()) #TODO Remove scenarios
    for s in scenarios:
        Demand[s], Elyser_consumption[s] = get_demand(inputs[s], results[s])

    labels = inputs[scenarios[0]]['sets']['n']

    x = np.arange(len(labels))  # the label locations
    if len(scenarios) == 2:
        width = 0.3
    else:
        width = 0.2

    fig, ax = plt.subplots(figsize=(13, 5))
    plt.rcParams.update({'font.size': 15})
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(15)
    if len(scenarios) == 2:
        ax.bar(x - width / 2 - 0.03, Demand[scenarios[0]].iloc[0, :] + Elyser_consumption[scenarios[0]].iloc[0, :],
               width, color=(1, 0.5, 0, 0.7))
        ax.bar(x - width / 2 - 0.03, Demand[scenarios[0]].iloc[0, :], width, color=(1, 0.5, 0, 1),
               label=str(scenarios[0]))
        ax.bar(x + width / 2 + 0.03, Demand[scenarios[1]].iloc[0, :] + Elyser_consumption[scenarios[1]].iloc[0, :],
               width, color=(0.192, 0.549, 0.905, 0.7))
        ax.bar(x + width / 2 + 0.03, Demand[scenarios[1]].iloc[0, :], width, color=(0.192, 0.549, 0.905),
               label=str(scenarios[1]))
    elif len(scenarios) == 3:
        ax.bar(x - 0.22, Demand[scenarios[0]].iloc[0, :] + Elyser_consumption[scenarios[0]].iloc[0, :], width,
               label=str(scenarios[0]), color=(1, 0.5, 0, 0.7))
        ax.bar(x - 0.22, Demand[scenarios[0]].iloc[0, :], width, color=(1, 0.5, 0, 1))
        ax.bar(x, Demand[scenarios[1]].iloc[0, :] + Elyser_consumption[scenarios[1]].iloc[0, :], width,
               label=str(scenarios[1]), color=(0.192, 0.549, 0.905, 0.7))
        ax.bar(x, Demand[scenarios[1]].iloc[0, :], width, color=(0.192, 0.549, 0.905))
        ax.bar(x + 0.22, Demand[scenarios[2]].iloc[0, :] + Elyser_consumption[scenarios[2]].iloc[0, :], width,
               label=str(scenarios[2]), color=(0.192, 0.549, 0.3, 0.7))
        ax.bar(x + 0.22, Demand[scenarios[2]].iloc[0, :], width, color=(0.192, 0.549, 0.3))

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('[GWh]')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    plt.show()

    return True


def plot_tech_cap(inputs, plot=True, figsize=(10, 7), alpha=0.8, width=0.5):
    """
    Plot installed storage capacity
    :param inputs:      Dictionary of standard inputs from the model
    :param plot:        Bool for plotting
    :param figsize:     Size of the figure to be plotted
    :param alpha:       Alpha value for the chart
    :param width:       Width of the bars in the plot
    :return:            Storage capacities
    """

    # 2 plot storage cap
    Cap = pd.DataFrame(columns=['BEVS', 'H2', 'Hydro', 'Thermal', 'BATS'], index=inputs['sets']['n']).fillna(0)
    for i, u in enumerate(inputs['param_df']['StorageCapacity'].index):
        for z in inputs['sets']['n']:
            if inputs['param_df']['Location'].loc[z, u]:
                if inputs['param_df']['Fuel'].loc['WAT', u]:
                    Cap.loc[z, 'Hydro'] += inputs['param_df']['StorageCapacity'].fillna(0).iloc[i, 0] * \
                                           inputs['param_df']['Nunits'].iloc[i, 0]
                elif inputs['param_df']['Fuel'].loc['HYD', u]:
                    Cap.loc[z, 'H2'] += inputs['param_df']['StorageCapacity'].fillna(0).iloc[i, 0] * \
                                        inputs['param_df']['Nunits'].iloc[i, 0]
                elif inputs['param_df']['Technology'].loc['BEVS', u]:
                    Cap.loc[z, 'BEVS'] += inputs['param_df']['StorageCapacity'].fillna(0).iloc[i, 0] * \
                                          inputs['param_df']['Nunits'].iloc[i, 0]
                elif inputs['param_df']['Technology'].loc['BATS', u]:
                    Cap.loc[z, 'BATS'] += inputs['param_df']['StorageCapacity'].fillna(0).iloc[i, 0] * \
                                          inputs['param_df']['Nunits'].iloc[i, 0]
                elif u in inputs['param_df']['sets']['p2h'] or inputs['param_df']['sets']['chp'] or \
                        inputs['param_df']['Technology'].loc['SCSP', u] or inputs['param_df']['sets']['thms']:
                    Cap.loc[z, 'Thermal'] += inputs['param_df']['StorageCapacity'].fillna(0).iloc[i, 0] * \
                                             inputs['param_df']['Nunits'].iloc[i, 0]

    Cap = Cap / 1000
    x = [i for i in range(len(Cap.index))]  # the label locations
    colors = ['#57D53B', '#A0522D', '#00a0e1ff', '#C04000ff', 'black']
    fig, ax = plt.subplots(figsize=figsize)

    ax = Cap.plot(kind="bar", figsize=figsize, stacked=True, color=colors, alpha=alpha, legend='reverse', width=width)
    ax.set_ylabel('Capacity [GWh]')
    ax.set_title('Installed storage capacity')
    ax.set_xlabel('Zones')
    if plot is True:
        plt.show()

    return Cap

# TODO Check if necessary
def H2_demand_satisfaction(inputs, results):
    H2Demand = [0]
    OutputH2 = [0]
    for i in range(len(inputs.keys()) - 1):
        H2Demand.append(0)
        OutputH2.append(0)
    for i, s in enumerate(list(inputs.keys())):
        for u in inputs[s]['param_df']['sets']['p2h2']:
            H2Demand[i] += (inputs[s]['param_df']['H2Demand'].loc[:, u].sum()
                            + inputs[s]['param_df']['PtLDemandInput'].loc[:, u].sum())
            OutputH2[i] += results[s]['OutputH2Output'].loc[:, u].sum()

    labels = list(inputs.keys())
    x = [i for i in range(len(labels))]  # the label locations
    width = 0.2  # the width of the bars
    colors = ['#41afaaff']
    fig, ax = plt.subplots(figsize=(3 * len(labels), 6))
    for i, s in enumerate(list(inputs.keys())):
        if i == 0:
            ax.bar(x[i], H2Demand[i] / 1e6, width, color=colors, hatch='//', label='H2 slack')
            ax.bar(x[i], OutputH2[i] / 1e6, width, color=colors, label='H2 produced by elysers')
        else:
            ax.bar(x[i], H2Demand[i] / 1e6, width, color=colors, hatch='//')
            ax.bar(x[i], OutputH2[i] / 1e6, width, color=colors)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('H2 demand [TWh]')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.grid(axis='y', linestyle='--')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
              ncol=2)

    return True


def heatmap(data, row_labels, col_labels, ax=None, cbar_kw=None, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.
    https://matplotlib.org/3.3.1/gallery/images_contours_and_fields/image_annotated_heatmap.html

    :param data:        A 2D numpy array of shape (N, M).
    :param row_labels:  A list or array of length N with the labels for the rows.
    :param col_labels:  A list or array of length M with the labels for the columns.
    :param ax:          A `matplotlib.axes.Axes`instance to which the heatmap is plotted. If not provided, use current
                        axes or create a new one.  Optional.
    :param cbar_kw:     A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    :param cbarlabel:   The label for the colorbar.  Optional.
    """

    if cbar_kw is None:
        cbar_kw = {}
    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=0, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title('Power Flow Tracing Matrix')
    ax.set_ylabel('Supplied zones')
    ax.set_xlabel('Power Generating zones')

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}", textcolors=("black", "white"), threshold=None, **textkw):
    """
    A function to annotate a heatmap.
    https://matplotlib.org/3.3.1/gallery/images_contours_and_fields/image_annotated_heatmap.html
    :param im:          The AxesImage to be labeled.
    :param data:        Data used to annotate.  If None, the image's data is used.  Optional.
    :param valfmt:      The format of the annotations inside the heatmap. This should either use the string format
                        method, e.g. "$ {x:.2f}", or be a `matplotlib.ticker.Formatter`.  Optional.
    :param textcolors:  A pair of colors.  The first is used for values below a threshold, the second for those above.
                        Optional.
    :param threshold:   Value in data units according to which the colors from textcolors are applied.
                        If None (the default) uses the middle of the colormap as separation.  Optional.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def plot_power_flow_tracing_matrix(inputs, results, idx=None, figsize=(10, 7), **kwargs):
    """
    Plot power flow tracing matrix
    :param inputs:      Dispa-SET inputs
    :param results:     Dispa-SET results
    :param idx:         datetime index, a range of dates to analyze. By default it looks only at the first day of the
                        optimization. Optional.
    :param figsize:     Figure size. Optional.
    :param kwargs:
    :return:            Power flow tracing in MW and as a % of the total load in a particular region
    """
    data, data_prct = get_power_flow_tracing(inputs, results, idx)
    fig, ax = plt.subplots(figsize=figsize)
    im, cbar = heatmap(data_prct.values, data_prct.index, data_prct.columns,
                       cbarlabel="% of the total demand linked to one zone", **kwargs)
    texts = annotate_heatmap(im, valfmt="{x:.1f}")
    fig.tight_layout()
    plt.show()

    return data, data_prct


def plot_co2(inputs, results, idx=None, figsize=(10,7), width=0.9, alpha=0.8, points=100, facecolor='#D43F3A'):
    """

    :param inputs:      Dictionary with model inputs
    :param results:     Dictionary with model outputs
    :param idx:         Index for selecting data
    :param figsize:     Size of the plot
    :return:
    """
    data = results['OutputEmissions'].T.swaplevel().T
    data[data > 100] = 0
    data = data.reindex(sorted(data.columns), axis=1)
    labels = list(data.loc[:, 'CO2'].columns)
    df = pd.DataFrame(np.sort(data.loc[:, 'CO2'].values, axis=0), columns=labels)
    df = df.to_numpy()
    fig, ax = plt.subplots(figsize=figsize)
    parts = ax.violinplot(df, points=points, widths=width, showmeans=False, showmedians=False,
        showextrema=False, bw_method=0.5)

    for pc in parts['bodies']:
        pc.set_facecolor(facecolor)
        pc.set_edgecolor('black')
        pc.set_alpha(alpha)

    quartile1, medians, quartile3 = np.percentile(df.T, [25, 50, 75], axis=1)
    whiskers = np.array([adjacent_values(sorted_array, q1, q3) for sorted_array, q1, q3 in zip(df.T, quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

    inds = np.arange(1, len(medians) + 1)
    ax.scatter(inds, medians, marker='o', color='white', s=40, zorder=3)
    ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=6)
    ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

    set_axis_style(ax, labels, xlabel='Zone', ylabel=r'$gCO_{2eq}/kWh$')
    ax.set_title('$CO_{2}$ Emissions')
    fig.tight_layout()
    plt.show()


def set_axis_style(ax, labels, xlabel='', ylabel=''):
    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_ylim(0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value
