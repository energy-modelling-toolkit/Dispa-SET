"""
Functions for handling boundary sectors in Dispa-SET.

This module contains all boundary sector-related functions that were previously in utils.py and build.py.
It handles the creation, formatting, and processing of boundary sectors.
"""

import os
import logging
import numpy as np
import pandas as pd

from ..common import commons
from .data_handler import load_time_series
from .interconnections import interconnections


def process_boundary_sector_data(config, zones_bs=None):
    """
    Process boundary sector data from configuration and return relevant dataframes

    :param config:           Dispa-SET configuration dictionary
    :param zones_bs:         List of boundary sector zones (optional)
    :return:                 Dictionary containing boundary sector related dataframes
    """
    result = {}
    
    # Initialize boundary sector dataframe
    result['BoundarySector'] = pd.DataFrame()
    
    # Load boundary sector data from file
    if 'BoundarySectorData' in config and os.path.isfile(config['BoundarySectorData']):
        result['BoundarySector'] = pd.read_csv(config['BoundarySectorData'],
                                   na_values=commons['na_values'],
                                   keep_default_na=False, index_col='Sector')
        for key in ['STOCapacity', 'STOSelfDischarge', 'STOMaxPower', 'STOMinSOC', 'STOHours']:
            if key not in result['BoundarySector'].columns:
                result['BoundarySector'][key] = np.nan
                logging.warning(f'{key} is not defined in the boundary sector data table')
    else:
        for key in ['STOCapacity', 'STOSelfDischarge', 'STOMaxPower', 'STOMinSOC', 'STOHours']:
            result['BoundarySector'][key] = np.nan
    
    # Load boundary sector interconnections data
    result['BS_flows'] = pd.DataFrame(index=config['idx_long'])
    if 'BoundarySectorInterconnections' in config and os.path.isfile(config['BoundarySectorInterconnections']):
        result['BS_flows'] = load_time_series(config, config['BoundarySectorInterconnections']).fillna(0)

    # Load boundary sector NTC data
    result['BS_NTC'] = pd.DataFrame(index=config['idx_long'])
    if 'BoundarySectorNTC' in config and os.path.isfile(config['BoundarySectorNTC']):
        result['BS_NTC'] = load_time_series(config, config['BoundarySectorNTC']).fillna(0)
    
    # Load spillage capacity data
    result['BS_spillage'] = pd.DataFrame(index=config['idx_long'])
    if 'BoundarySectorMaxSpillage' in config and os.path.isfile(config['BoundarySectorMaxSpillage']):
        result['BS_spillage'] = load_time_series(config, config['BoundarySectorMaxSpillage']).fillna(0)
    else:
        logging.warning('No maximum spillage capacity provided.')
    
    # Initialize forced spillage with zeros (if spillage capacity exists)
    result['BS_forced_spillage'] = pd.DataFrame(index=config['idx_long'])
    if 'BoundarySectorMaxSpillage' in config and os.path.isfile(config['BoundarySectorMaxSpillage']):
        result['BS_forced_spillage'] = pd.DataFrame(0, index=result['BS_spillage'].index, columns=result['BS_spillage'].columns)
    
    # Load costXspillage data
    result['CostXSpillage'] = pd.DataFrame(index=config['idx_long'])
    if 'CostXSpillage' in config and config['CostXSpillage'] != '' and os.path.isfile(config['CostXSpillage']):
        result['CostXSpillage'] = load_time_series(config, config['CostXSpillage']).fillna(0)
    
    return result


def process_boundary_sector_interconnections(config, bs_data, zones_bs):
    """
    Process boundary sector interconnections data and return the processed data

    :param config:           Dispa-SET configuration dictionary
    :param bs_data:          Dictionary with boundary sector data from process_boundary_sector_data
    :param zones_bs:         List of boundary sector zones
    :return:                 Dictionary containing processed interconnection dataframes
    """
    result = {}
    
    # Process boundary sector interconnections
    [BSInterconnections_sim, BSInterconnections_RoW, BSInterconnections] = interconnections(
        zones_bs, bs_data['BS_NTC'], bs_data['BS_flows'])
    
    # Process boundary sector spillage
    [BSSpillage_sim, BSSpillage_RoW, BSSpillage] = interconnections(
        zones_bs, bs_data['BS_spillage'], bs_data['BS_forced_spillage'])
    
    # Format the results
    if len(BSInterconnections_sim.columns) > 0:
        result['BS_NTCs'] = BSInterconnections_sim.reindex(config['idx_long'])
    else:
        result['BS_NTCs'] = pd.DataFrame(index=config['idx_long'])
    
    result['BS_Inter_RoW'] = BSInterconnections_RoW.reindex(config['idx_long'])
    
    if len(BSSpillage_sim.columns) > 0:
        result['BS_Spillages'] = BSSpillage_sim.reindex(config['idx_long'])
    else:
        result['BS_Spillages'] = pd.DataFrame(index=config['idx_long'])
    
    result['BS_Spillage_RoW'] = BSSpillage_RoW.reindex(config['idx_long'])
    
    # Store the list of all interconnections
    result['BSInterconnections'] = BSInterconnections
    result['BSSpillage'] = BSSpillage
    
    return result


def get_boundary_sector_zones(boundary_sector_df):
    """
    Extract list of boundary sector zones from boundary sector dataframe

    :param boundary_sector_df: Boundary sector dataframe
    :return:                   List of boundary sector zones
    """
    if len(boundary_sector_df.index) > 0:
        zones_bs = boundary_sector_df.index.to_list()
        # Clean up the zone list
        if '' in zones_bs:
            zones_bs.remove('')
        if np.nan in zones_bs:
            zones_bs = [x for x in zones_bs if pd.isnull(x) == False]
        if 'nan' in zones_bs:
            zones_bs.remove('nan')
        return zones_bs
    else:
        return []


def zone_to_bs_mapping(plants_all_bs):
    """
    Creates a mapping function from boundary sector indices to their corresponding zones

    :param plants_all_bs:  DataFrame containing all boundary sector plants with their zone information
    :returns:              Function that maps boundary sector index to its zone
    """
    # Create a mapping dictionary from boundary sector index to zone
    bs_to_zone = {}
    for idx in plants_all_bs.index:
        # Get all columns that start with 'Sector'
        sector_cols = [col for col in plants_all_bs.columns if col.startswith('Sector')]
        for col in sector_cols:
            if not pd.isna(plants_all_bs.loc[idx, col]):
                bs_to_zone[plants_all_bs.loc[idx, col]] = plants_all_bs.loc[idx, 'Zone']

    def mapping_function(bs_idx):
        return bs_to_zone.get(bs_idx, None)

    return mapping_function


def boundary_sector_efficiency_time_series(config, plants, zones_bs):
    """
    Function that calculates boundary sector efficiency time series for each unit
    In case of generation unit, the efficiency is defined as the EfficiencySectorX
    In case of of power consumption units, the efficiency is defined as the ChargingEfficiencySectorX, (for now)

    :param config:          Dispa-SET config file
    :param plants:          Pandas dataframe with the original list of units
    :param zones_bs:        Boundary sector zones

    :returns:               Dictionary of Dataframes with a time series of the efficiency for each unit
    """
    # Get the sector columns (Sector1, Sector2, etc.)
    sector_cols = [col for col in plants if col.startswith('Sector') and not col.startswith('Efficiency')]
    Efficiencies = {}
    
    for n in zones_bs:
        Efficiencies[n] = pd.DataFrame(columns=plants.index, index=config['idx_long'])
        for s in sector_cols:
            for u in plants.index:
                if ((plants.loc[u, 'Technology'] in commons['tech_p2bs'] or plants.loc[u, 'Technology'] in commons['tech_boundary_sector']) and 
                    (plants.loc[u, s] == n)):
                    # For p2x units, use EfficiencySectorX for Power2XConversionMultiplier
                    # For x2p units, use EfficiencySectorX for X2PowerConversionMultiplier
                    eff_col = 'EfficiencySector' + s[6:]  # Convert 'Sector1' to 'EfficiencySector1'
                    Efficiencies[n][u] = plants.loc[u, eff_col]
        if n in Efficiencies:
            Efficiencies[n] = Efficiencies[n].fillna(0)
    
    return {'Efficiency': Efficiencies} 