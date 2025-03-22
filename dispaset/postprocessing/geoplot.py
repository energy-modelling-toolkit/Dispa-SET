# -*- coding: utf-8 -*-

import logging

import cartopy
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.mpl.geoaxes
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import requests

from .postprocessing import get_from_to_flows, get_net_positions


# Insipired by Pypsa geo plot
def get_projection_from_crs(crs):
    """
    EPSG coordinate system, 4326 is the default horizontal component of 3D system. Used by the GPS satellite navigation
    system and for NATO military geodetic surveying.

    :param crs:     EPSG:4326 (WGS84 Bounds:        -180.0000, -90.0000, 180.0000, 90.0000,
                               Projected Bounds:    -180.0000, -90.0000, 180.0000, 90.0000,
                               Area:                World)
    :return:        map projection
    """
    # if data is in lat-lon system, return default map with lat-lon system
    if crs == 4326:
        return ccrs.PlateCarree()
    try:
        return ccrs.epsg(crs)
    except requests.RequestException:
        logging.warning("A connection to http://epsg.io/ is required for a projected coordinate reference system. "
                        "Falling back to lat-long.")
    except ValueError:
        logging.warning("'{crs}' does not define a projected coordinate system. "
                        "Falling back to lat-long.".format(crs=crs))
        return ccrs.PlateCarree()


def compute_bbox_with_margins(x, y, margin_type='Fixed', margin=0):
    """
    'Helper function to compute bounding box for the plot'

    :param margin:  how much percent is box beyond the ploted network
    :param x:       latitude
    :param y:       longitude
    :return:        tuples of 2 coordinates (min, max)
    """
    # set margins
    pos = np.asarray((x, y))
    minxy, maxxy = pos.min(axis=1), pos.max(axis=1)
    if margin_type == 'Fixed':
        xy1 = minxy - margin
        xy2 = maxxy + margin
    else:
        xy1 = minxy - margin * (maxxy - minxy)
        xy2 = maxxy + margin * (maxxy - minxy)

    return tuple(xy1), tuple(xy2)


def draw_map_cartopy(x, y, ax, crs=4326, boundaries=None, margin_type='Fixed', margin=0.05, geomap=True,
                     color_geomap=None, terrain=False):
    """
    This function draws a map in the background of the plot

    :param x:               geo['Longitude']
    :param y:               geo['Latitude']
    :param ax:              Plot axis
    :param crs:             Geographic projection code (i.e. 4326 is a standard EPSG coordinate system)
    :param boundaries:      if boundaries is None "compute_bbox_with_margins(x, y, margin_type, margin)" will be used,
                            otherwise boundaries need to be specified manually "x1, x2, y1, y2 = boundaries"
    :param margin_type:     'Fixed': xy1 = minxy - margin, xy2 = maxxy + margin, otherwise:
                            xy1 = minxy - margin * (maxxy - minxy), xy2 = maxxy + margin * (maxxy - minxy)
    :param margin:          Margin for extending the plotting range
    :param geomap:          Boolean for plotting the edges
    :param color_geomap:    Background colors. Should be specified in form: {'ocean': 'lightblue', 'land': 'whitesmoke'}
    :param terrain:         Boolean for plotting the background terrain, if False only country borders will be drawn
    :return:                Axis transformation
    """
    
    if boundaries is None:
        (x1, y1), (x2, y2) = compute_bbox_with_margins(x, y, margin_type, margin)
    else:
        x1, x2, y1, y2 = boundaries

    resolution = '50m' if isinstance(geomap, bool) else geomap
    assert resolution in ['10m', '50m', '110m'], ("Resolution has to be one of '10m', '50m', '110m'")
    axis_transformation = get_projection_from_crs(crs)
    ax.set_extent([x1, x2, y1, y2], crs=axis_transformation)

    if color_geomap is None:
        color_geomap = {'ocean': 'w', 'land': 'w'}
    elif color_geomap and not isinstance(color_geomap, dict):
        color_geomap = {'ocean': 'lightblue', 'land': 'whitesmoke'}

    ax.add_feature(cartopy.feature.LAND.with_scale(resolution), facecolor=color_geomap['land'])
    ax.add_feature(cartopy.feature.OCEAN.with_scale(resolution), facecolor=color_geomap['ocean'])

    ax.coastlines(linewidth=1, zorder=2, resolution=resolution)
    border = cartopy.feature.BORDERS.with_scale(resolution)
    ax.add_feature(border, linewidth=0.8)

    #FIXME: Sometimes there is an issue with connecting to the stamen server and terain can not be plotted. I dont know how to fix this!
    if terrain is True:
        # Create a Stamen terrain background instance.
        stamen_terrain = cimgt.Stamen('terrain')
        # Add the Stamen data at zoom level 8.
        ax.add_image(stamen_terrain, 8)

    return axis_transformation


def get_congestion(inputs, flows, idx):
    """
    This function computes the congestion in the interconnection lines (power system only)

    :param inputs:  Dispa-SET inputs
    :param flows:   Flows between neighbouring (interconnected) nodes
    :param idx:     Datetime index defining the range for which the congestion should be computed
    :return:        Congestion
    """

    cols = [x for x in inputs['sets']['l'] if "RoW" not in x]
    congestion = pd.DataFrame(columns=cols)
    for l in flows:
        if l[:3] != 'RoW' and l[-3:] != 'RoW':
            congestion.loc[0, l] = (
                    (flows.loc[idx, l] == inputs['param_df']['FlowMaximum'].loc[idx, l]) & (
                    inputs['param_df']['FlowMaximum'].loc[idx, l] > 0)).sum()

    congestion.fillna(0, inplace=True)
    congestion = congestion / inputs['param_df']['Demand'].loc[idx, :].index.size
    return congestion


def plot_net_flows_map(inputs, results, idx=None, crs=4326, boundaries=None, margin_type='Fixed', margin=0.20,
                       geomap=True, color_geomap=None, terrain=False, figsize=(12, 8), bublesize=5000):
    """
    Function for generating the net flows map diagram

    :param inputs:          Dispa-SET inputs
    :param results:         Dispa-SET results
    :param idx:             Datetime index range
    :param crs:             Geographic projection code (i.e. 4326 is a standard EPSG coordinate system)
    :param boundaries:      if boundaries is None "compute_bbox_with_margins(x, y, margin_type, margin)" will be used,
                            otherwise boundaries need to be specified manually "x1, x2, y1, y2 = boundaries"
    :param margin_type:     'Fixed': xy1 = minxy - margin, xy2 = maxxy + margin, otherwise:
                            xy1 = minxy - margin * (maxxy - minxy), xy2 = maxxy + margin * (maxxy - minxy)
    :param margin:          Margin for extending the plotting range
    :param geomap:          Boolean for plotting the edges
    :param color_geomap:    Background colors. Should be specified in form: {'ocean': 'lightblue', 'land': 'whitesmoke'}
    :param terrain:         Boolean for plotting the background terrain, if False only country borders will be drawn
    :param figsize:         Figure size
    :param bublesize:       Size of the bubbles. This value is impacted by the figsize.
    :return:                Net flows plot
    """

    # Preprocess input data
    flows = results['OutputFlow'].copy()
    zones = inputs['sets']['n'].copy()
    geo = inputs['geo'].copy()

    # Checking if index was selected
    if idx is None:
        idx = inputs['param_df']['Demand'].index
        logging.info('No datetime range specified, net flows map is printed for the entire optimization')
    else:
        idx = idx

    Flows = get_from_to_flows(inputs, flows, zones, idx)
    NetImports, P = get_net_positions(inputs, results, zones, idx)

    # Scale net position of the zone TODO: maybe apply some other algorithm instead of scaling based on the highest net position
    P = P / P.max()

    # Create a directed graph
    g = nx.DiGraph()
    # Add zones
    g.add_nodes_from(zones)
    # Define and add edges
    edges = Flows[['From', 'To']].values
    g.add_edges_from(edges)
    # Assign weights (between 0 - 10 seems quite reasonable), weights are sized according to the highest one
    # TODO: Not sure if this is somethign we are after, maybe there should be another method used for scaling i.e. based on max NTC
    weights = (10 * Flows['Flow'] / Flows['Flow'].max()).values

    # Define geospatial coordinates
    pos = {zone: (v['Longitude'], v['Latitude']) for zone, v in geo.to_dict('index').items()}

    # Node sizes (Based on the net position of a zone)
    sizes = [bublesize * P[i] for i in g.nodes]

    # Assign colors based on net flows (if importing/exporting/neutral)
    node_neg = NetImports.columns[(NetImports < 0).any()].tolist()
    node_pos = NetImports.columns[(NetImports > 0).any()].tolist()

    color_map = []
    for node in g:
        if node in node_neg:
            color_map.append('green')
        elif node in node_pos:
            color_map.append('red')
        else:
            color_map.append('blue')

    # Show labels only in nodes whose size is > 100
    labels = {i: i if bublesize * P[i] >= 150 else '' for i in g.nodes}

    # Define projection (FIXME: currently only 4326 possible)
    projection = get_projection_from_crs(4326)

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw=dict(projection=projection))
    title = "Power feed (red=Imports, green=Exports, blue=Neutral)"

    # Assign geo coordinates and draw them on a map
    x, y = geo['Longitude'], geo['Latitude']
    transform = draw_map_cartopy(x, y, ax, crs, boundaries, margin_type, margin, geomap, color_geomap, terrain)
    x, y, z = ax.projection.transform_points(transform, x.values, y.values).T

    x, y = pd.Series(x, geo.index), pd.Series(y, geo.index)

    # Draw networkx graph with nodes and edges
    nx.draw_networkx(g, ax=ax, font_size=16,
                     # alpha=.7,
                     width=weights,
                     node_size=sizes,
                     labels=labels,
                     pos=pos,
                     node_color=color_map,
                     cmap=plt.cm.autumn,
                     arrows=True, arrowstyle='-|>', arrowsize=20, connectionstyle='arc3, rad=0.05',
                     )

    ax.update_datalim(compute_bbox_with_margins(x, y, margin_type, margin))
    ax.autoscale_view()

    if geomap:
        for spine in ax.spines.values():
            spine.set_visible(False)
    else:
        ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title)

    plt.show()


# TODO: Generalize this function and provide descriptions
def plot_line_congestion_map(inputs, results, idx=None, crs=4326, boundaries=None, margin_type='Fixed', margin=0.20,
                             geomap=True, color_geomap=None, terrain=False, figsize=(12, 8), edge_width=10,
                             bublesize=5000):
    """
    This function plots congestion in the lines

    :param inputs:          DispaSET inputs
    :param results:         DispaSET results
    :param idx:             Datetime index range
    :param crs:             Geographic projection code (i.e. 4326 is a standard EPSG coordinate system)
    :param boundaries:      if boundaries is None "compute_bbox_with_margins(x, y, margin_type, margin)" will be used,
                            otherwise boundaries need to be specified manually "x1, x2, y1, y2 = boundaries"
    :param margin_type:     'Fixed': xy1 = minxy - margin, xy2 = maxxy + margin, otherwise:
                            xy1 = minxy - margin * (maxxy - minxy), xy2 = maxxy + margin * (maxxy - minxy)
    :param margin:          Margin for extending the plotting range
    :param geomap:          Boolean for plotting the edges
    :param color_geomap:    Background colors. Should be specified in form: {'ocean': 'lightblue', 'land': 'whitesmoke'}
    :param terrain:         Boolean for plotting the background terrain, if False only country borders will be drawn
    :param figsize:         Figure size
    :param edge_width:      Thickness of the NTC lines, 10 usually works well for the default figsize
    :param bublesize:       Size of the bubbles. This value is impacted by the figsize.
    :return:                Plot the diagram
    """

    # Preprocess input data
    zones = inputs['sets']['n'].copy()
    geo = inputs['geo'].copy()
    flows = results['OutputFlow'].copy()

    # Checking if index was selected
    if idx is None:
        idx = inputs['param_df']['Demand'].index
        logging.info('No datetime range specified, net flows map is printed for the entire optimization')
    else:
        idx = idx

    cgst = get_congestion(inputs, flows, idx)

    Congestion = get_from_to_flows(inputs, cgst, zones)

    # Create a directed graph
    g = nx.DiGraph()
    # Add zones
    g.add_nodes_from(zones)
    # Define and add edges
    edges = Congestion[['From', 'To']].values
    g.add_edges_from(edges)
    # Assign weights (between 0 - 10 seems quite reasonable), weights are sized according to the highest one
    # TODO: Not sure if this is somethign we are after, maybe there should be another method used for scaling i.e. based on max NTC
    weights = (100 * Congestion['Flow']).values

    # Define geospatial coordinates
    pos = {zone: (v['Longitude'], v['Latitude']) for zone, v in geo.to_dict('index').items()}

    # Node sizes (Based on the net position of a zone)
    sizes = [bublesize for i in g.nodes]

    # Show labels only in nodes whose size is > 100
    labels = {i: i if bublesize >= 500 else '' for i in g.nodes}

    # Define projection (FIXME: currently only 4326 possible)
    projection = get_projection_from_crs(4326)

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw=dict(projection=projection))
    title = "Line Congestion (Congestion levels: dark_red=High, green=Middle, blue=None)"

    # Assign geo coordinates and draw them on a map
    x, y = geo['Longitude'], geo['Latitude']
    transform = draw_map_cartopy(x, y, ax, crs, boundaries, margin_type, margin, geomap, color_geomap, terrain)
    x, y, z = ax.projection.transform_points(transform, x.values, y.values).T

    x, y = pd.Series(x, geo.index), pd.Series(y, geo.index)

    edge_cmap = plt.cm.jet
    edge_vmin = 0.0
    edge_vmax = 100
    cmap = plt.cm.autumn

    # Draw networkx graph with nodes and edges
    nx.draw_networkx(g, ax=ax, font_size=18,
                     # alpha=.7,
                     width=edge_width,
                     node_size=sizes,
                     node_color='black',
                     labels=labels,
                     pos=pos,
                     edge_color=weights,
                     edge_cmap=edge_cmap,
                     cmap=cmap,
                     arrows=True, arrowstyle='-|>', arrowsize=20, connectionstyle='arc3, rad=0.05',
                     edge_vmin=edge_vmin, edge_vmax=edge_vmax
                     )

    sm = plt.cm.ScalarMappable(cmap=edge_cmap, norm=plt.Normalize(vmin=edge_vmin, vmax=edge_vmax))
    sm.set_array([])
    
    # Create colorbar with proper axes
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.02)
    cbar.set_label('Congestion [%]')

    ax.update_datalim(compute_bbox_with_margins(x, y, margin_type, margin))
    ax.autoscale_view()

    if geomap:
        for spine in ax.spines.values():
            spine.set_visible(False)
    else:
        ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title)

    plt.tight_layout()
    plt.show()
