# -*- coding: utf-8 -*-
"""
This is the main file of the DispaSET pre-processing tool. It comprises a single function that generated the DispaSET simulation environment.

@author: S. Quoilin
"""
import datetime as dt
import logging
import sys

import pandas as pd

from .build import build_single_run

from ..solve import solve_GAMS
from ..common import commons
from ..misc.gdx_handler import gdx_to_dataframe, gdx_to_list

try:
    from future.builtins import int
except ImportError:
    logging.warning("Couldn't import future package. Numeric operations may differ among different versions due to incompatible variable types")
    pass


def build_simulation(config, mts_plot=None):
    '''
    Dispa-SET function that builds different simulation environments based on the hydro scheduling option in the config file
    Hydro scheduling options:
        Off      -  Hydro scheduling turned off, normal call of BuildSimulation function
        Zonal    -  Zonal variation of hydro scheduling, if zones are not individually specified in a list (e.a. zones = ['AT','DE'])
                    hydro scheduling is imposed on all active zones from the Config file
        Regional -  Regional variation of hydro scheduling, if zones from a specific region are not individually specified in a list
                    (e.a. zones = ['AT','DE']), hydro scheduling is imposed on all active zones from the Config file simultaneously

    :config:                    Read config file
    :zones_mts:                 List of zones where new reservoir levels should be calculated eg. ['AT','BE',...'UK']
    :mts_plot:                  If ms_plot = True indicative plot with temporary computed reservoir levels is displayed
    '''
    y_start, m_start, d_start, __, __, __ = config['StartDate']
    y_stop, m_stop, d_stop, __, __, __ = config['StopDate']
    # Check existance of hydro scheduling module in the config file
    # Hydro scheduling turned off, build_simulation performed without temporary computed reservoir levels
    hydro_flag = config.get('HydroScheduling', "")  # If key does not exist it returns ""
    if (hydro_flag == "") or (hydro_flag == "Off"):
        logging.info('Simulation without mid therm scheduling')
        SimData = build_single_run(config)
    # Hydro scheduling per Zone
    else:
        # Dates for the mid term scheduling
        if config['HydroSchedulingHorizon'] == 'Annual':
            config['StartDate'] = (y_start, 1, 1, 00, 00, 00)  # updating start date to the beginning of the year
            config['StopDate'] = (y_start, 12, 31, 23, 59, 00)  # updating stopdate to the end of the year
            logging.info('Hydro scheduling is performed for the period between 01.01.' + str(y_start) + ' and 12.31.' + str(y_start))
        else:
            logging.info('Hydro scheduling is performed between Start and Stop dates!')
            # Mid term scheduling zone selection and new profile calculation
        if config['mts_zones'] is None:
            new_profiles = mid_term_scheduling(config, config['zones'])
            logging.info('Simulation with all zones selected')
        else:
            new_profiles = mid_term_scheduling(config, config['mts_zones'])
        # Plot new profiles
        if mts_plot:
            new_profiles.plot()
            logging.info('Simulation with specified zones selected')
        else:
            logging.info('No temporary profiles selected for display')
        # Build simulation data with new profiles
        config['StartDate'] = (y_start, m_start, d_start, 00, 00, 00)  # updating start date to the beginning of the year
        config['StopDate'] = (y_stop, m_stop, d_stop, 23, 59, 00)  # updating stopdate to the end of the year
        SimData = build_single_run(config, new_profiles)
    return SimData


def mid_term_scheduling(config, zones, profiles=None):
    """
    This function reads the DispaSET config file, searches for active zones,
    loads data for each zone individually and solves model using UCM_h_simple.gms

    :config:                    Read config file
    """

    # Day/hour corresponding to the first and last days of the simulation:
    # Note that the first available data corresponds to 2015.01.31 (23.00) and the
    # last day with data is 2015.12.31 (22.00)
    y_end, m_end, d_end, _, _, _ = config['StopDate']
    config['StopDate'] = (y_end, m_end, d_end, 23, 59, 00)  # updating stopdate to the end of the day

    # Indexes of the simulation:
    idx_std = pd.date_range(start=dt.datetime(*config['StartDate']),
                            end=dt.datetime(*config['StopDate']),
                            freq=commons['TimeStep'])
    idx_utc_noloc = idx_std - dt.timedelta(hours=1)
    idx = idx_utc_noloc

    # Checking which type of hydro scheduling simulation is specified in the config file:
    # Solving reservoir levels for each zone individually
    if config['HydroScheduling'] == 'Zonal':
        no_of_zones = len(zones)
        results = {}
        temp_results = {}
        i = 0
        for c in zones:
            i = i + 1
            logging.info('(Currently simulating Zone): ' + str(i) + ' out of ' + str(no_of_zones))
            temp_config = dict(config)
            temp_config['zones'] = [c]                    # Override zone that needs to be simulated
            _ = build_single_run(temp_config)       # Create temporary SimData
            r = solve_GAMS(sim_folder=temp_config['SimulationDirectory'],
                           gams_folder=temp_config['GAMS_folder'],
                           gams_file='UCM_h_simple.gms',
                           result_file='Results_simple.gdx')
            temp_results[c] = get_temp_sim_results(config)
    #        print('Zones simulated: ' + str(i) + '/' + str(no_of_zones))
        temp = pd.DataFrame()
        for c in zones:
            if 'OutputStorageLevel' not in temp_results[c]:
                logging.critical('Storage levels in zone ' + c + ' were not computed, please check that storage units '
                                'are present in the ' + c + ' power plant database! If not, unselect ' + c + ' form the '
                                'zones in the MTS module')
                sys.exit(0)
            else:
                results[c] = dict(temp_results[c]['OutputStorageLevel'])
                r = pd.DataFrame.from_dict(results[c], orient='columns')
                results_t = pd.concat([temp, r], axis = 1)
                temp = results_t
        temp = temp.set_index(idx)
        temp = temp.rename(columns={col: col.split(' - ')[1] for col in temp.columns})
    # Solving reservoir levels for all regions simultaneously
    elif config['HydroScheduling'] == 'Regional':
        if zones is None:
            temp_config = dict(config)
        else:
            temp_config = dict(config)
            temp_config['zones'] = zones               # Override zones that need to be simmulated
        _ = build_single_run(temp_config)        # Create temporary SimData
        r = solve_GAMS(sim_folder=temp_config['SimulationDirectory'],
                       gams_folder=temp_config['GAMS_folder'],
                       gams_file='UCM_h_simple.gms',
                       result_file='Results_simple.gdx')
        temp_results = get_temp_sim_results(config)
        if 'OutputStorageLevel' not in temp_results:
            logging.critical('Storage levels in the selected region were not computed, please check that storage units '
                             'are present in the power plant database! If not, unselect zones with no storage units form '
                             'the zones in the MTS module')
            sys.exit(0)
        else:
            results = dict(temp_results['OutputStorageLevel'])
        temp = pd.DataFrame()
        r = pd.DataFrame.from_dict(results, orient='columns')
        results_t = pd.concat([temp, r], axis=1)
        temp = results_t
        temp = temp.set_index(idx)
        temp = temp.rename(columns={col: col.split(' - ')[1] for col in temp.columns})
    else:
        logging.info('Mid term scheduling turned off')
    import pickle
    pickle.dump(temp, open("temp_profiles.p", "wb"))
    return temp


def get_temp_sim_results(config, gams_dir=None):
    """
    This function reads the simulation environment folder once it has been solved and loads
    the input variables together with the results.

    :param path:                Relative path to the simulation environment folder (current path by default)
    :param cache:               If true, caches the simulation results in a pickle file for faster loading the next time
    :param temp_path:           Temporary path to store the cache file
    :returns inputs,results:    Two dictionaries with all the outputs
    """

    resultfile = config['SimulationDirectory'] + '/Results_simple.gdx'
    results = gdx_to_dataframe(gdx_to_list(gams_dir, resultfile, varname='all', verbose=True), fixindex=True,
                               verbose=True)
    results['OutputStorageLevel'] = results['OutputStorageLevel'].reindex(list(range(results['OutputStorageLevel'].index.min(),
                                                         results['OutputStorageLevel'].index.max() + 1)), fill_value=0)
    return results