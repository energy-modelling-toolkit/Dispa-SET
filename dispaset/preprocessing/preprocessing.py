# -*- coding: utf-8 -*-
"""
This is the main file of the DispaSET pre-processing tool.
It comprises a single function that generates the DispaSET simulation environment.

@author: S. Quoilin
"""
import datetime as dt
import logging
import sys

import pandas as pd

from .build import build_single_run

from ..solve import solve_GAMS
from ..misc.gdx_handler import gdx_to_dataframe, gdx_to_list
from .utils import pd_timestep

try:
    from future.builtins import int
except ImportError:
    logging.warning(
        "Couldn't import future package. Numeric operations may differ among different versions due to incompatible "
        "variable types")
    pass


def build_simulation(config, mts_plot=None, MTSTimeStep=None):
    """
    A function that builds different simulation environments based on the hydro scheduling option in the config file
    Hydro scheduling options:
        Off      -  Hydro scheduling turned off, normal call of BuildSimulation function
        Zonal    -  Zonal variation of hydro scheduling, if zones are not individually specified in a list
                    (e.a. zones = ['AT','DE']) hydro scheduling is imposed on all active zones from the Config file
        Regional -  Regional variation of hydro scheduling, if zones from a specific region are not individually
                    specified in a list (e.a. zones = ['AT','DE']), hydro scheduling is imposed on all active zones
                    from the Config file simultaneously

    :param config:          Read config file
    :param mts_plot:        If ms_plot = True indicative plot with temporary computed reservoir levels is displayed
    :param MTSTimeStep:     Run the mid-term scheduling with a different (to speed things up). If unspecified,
                            the old MTS formulation is used
    :return SimData:        Simulation data for unit-commitment module
    """
    # Check existance of hydro scheduling module in the config file
    hydro_flag = config.get('HydroScheduling', "")  # If key does not exist it returns ""
    if (hydro_flag == "") or (hydro_flag == "Off"):
        logging.info('Simulation without mid therm scheduling')
        SimData = build_single_run(config)
    else:
        new_profiles = mid_term_scheduling(config, mts_plot=mts_plot, TimeStep=MTSTimeStep)
        # Build simulation data with new profiles
        SimData = build_single_run(config, new_profiles)
    return SimData


def mid_term_scheduling(config, TimeStep=None, mts_plot=None):
    """
    This function reads the DispaSET config file, searches for active zones,
    loads data for each zone individually and solves model using UCM_h_simple.gms

    :param config:          Read config file
    :param TimeStep:        Time step (1, 2, 3, 4, 6, 8, 12, 24) number of hours to be considered at once.
    :param mts_plot:        If ms_plot = True indicative plot with temporary computed reservoir levels is displayed
    :return profiles:       Newly computed profile levels
    """

    # Day/hour corresponding to the first and last days of the simulation:
    # Note that the first available data corresponds to 2015.01.31 (23.00) and the
    # last day with data is 2015.12.31 (22.00)
    import pickle
    y_start, m_start, d_start, __, __, __ = config['StartDate']
    y_end, m_end, d_end, _, _, _ = config['StopDate']

    # New configuration dictionary, specific to MTS:
    temp_config = dict(config)

    # Dates for the mid term scheduling
    if config['HydroSchedulingHorizon'] == 'Annual':
        temp_config['StartDate'] = (y_start, 1, 1, 00, 00, 00)  # updating start date to the beginning of the year
        temp_config['StopDate'] = (y_start, 12, 31, 23, 59, 00)  # updating stopdate to the end of the year
        logging.info(
            'Hydro scheduling is performed for the period between 01.01.' + str(y_start) + ' and 12.31.' + str(y_start))
    else:
        logging.info('Hydro scheduling is performed between Start and Stop dates!')

    # Indexes with the original time step:
    idx_orig = pd.date_range(start=dt.datetime(*config['StartDate']),
                             end=dt.datetime(*config['StopDate']),
                             freq=pd_timestep(config['SimulationTimeStep'])).tz_localize(None)

    if config['mts_zones'] is None:
        temp_config['mts_zones'] = config['zones']
        logging.info('Simulation with all zones selected')
    else:
        temp_config['zones'] = config['mts_zones']

    # remove look ahead:
    temp_config['LookAhead'] = 0

    # Adjust time step:
    if TimeStep is not None:
        # Indexes of the MTS simulation:
        idx = pd.date_range(start=dt.datetime(*temp_config['StartDate']),
                            end=dt.datetime(*temp_config['StopDate']),
                            freq=pd_timestep(TimeStep)).tz_localize(None)
        temp_config['SimulationTimeStep'] = TimeStep
        gams_file = 'UCM_h.gms'
        temp_config['HorizonLength'] = (idx[-1] - idx[0]).days
        resultfile = 'Results.gdx'
    else:
        idx = pd.date_range(start=dt.datetime(*temp_config['StartDate']),
                            end=dt.datetime(*temp_config['StopDate']),
                            freq=pd_timestep(temp_config['SimulationTimeStep'])).tz_localize(None)
        gams_file = 'UCM_h_simple.gms'
        resultfile = 'Results_simple.gdx'

    # Checking which type of hydro scheduling simulation is specified in the config file:
    if config['HydroScheduling'] == 'Zonal':
        if config['SimulationType'] == 'Standard':
            logging.critical("SimulationType: 'Standard' and HydroScheduling: 'Zonal' not supported! Please"
                             " choose different input options")
            sys.exit(1)
        else:
            no_of_zones = len(config['mts_zones'])
            temp_results = {}
            for i, c in enumerate(config['mts_zones']):
                logging.info('(Currently simulating Zone): ' + str(i + 1) + ' out of ' + str(no_of_zones))
                temp_config['zones'] = [c]  # Override zone that needs to be simulated
                _ = build_single_run(temp_config)  # Create temporary SimData
                r = solve_GAMS(sim_folder=temp_config['SimulationDirectory'],
                               gams_folder=temp_config['GAMS_folder'],
                               gams_file=gams_file,
                               result_file=resultfile)
                temp_results[c] = gdx_to_dataframe(
                    gdx_to_list(config['GAMS_folder'], config['SimulationDirectory'] + '/' + resultfile, varname='all',
                                verbose=True), fixindex=True, verbose=True)
            profiles = pd.DataFrame(index=idx)
            for c in config['mts_zones']:
                if 'OutputStorageLevel' not in temp_results[c]:
                    logging.critical('Storage levels in zone ' + c + ' were not computed, please check that storage units '
                                     'are present in the ' + c + ' power plant database! If not, unselect ' + c +
                                     ' form the zones in the MTS module')
                    sys.exit(0)
                elif len(temp_results[c]['OutputStorageLevel']) > len(idx):
                    logging.critical('The number of time steps in the mid-term simulation results (' + str(
                                     len(temp_results[c]['OutputStorageLevel'])) +
                                     ') does not match the length of the index (' + str(len(idx)) + ')')
                    sys.exit(0)
                elif len(temp_results[c]['OutputStorageLevel']) < len(idx):
                    for u in temp_results[c]['OutputStorageLevel']:
                        profiles[u] = temp_results[c]['OutputStorageLevel'][u].reindex(range(1, len(idx) + 1)).fillna(0).values
                else:
                    for u in temp_results[c]['OutputStorageLevel']:
                        profiles[u] = temp_results[c]['OutputStorageLevel'][u].values

    # Solving reservoir levels for all regions simultaneously
    elif config['HydroScheduling'] == 'Regional':
        _ = build_single_run(temp_config)  # Create temporary SimData
        r = solve_GAMS(sim_folder=temp_config['SimulationDirectory'],
                       gams_folder=temp_config['GAMS_folder'],
                       gams_file=gams_file,
                       result_file=resultfile)
        temp_results = gdx_to_dataframe(
            gdx_to_list(config['GAMS_folder'], config['SimulationDirectory'] + '/' + resultfile, varname='all',
                        verbose=True), fixindex=True, verbose=True)
        if 'OutputStorageLevel' not in temp_results:
            logging.critical('Storage levels in the selected region were not computed, please check that storage units '
                             'are present in the power plant database! If not, unselect zones with no storage units '
                             'form the zones in the MTS module')
            sys.exit(0)
        elif len(temp_results['OutputStorageLevel']) < len(idx):
            profiles = temp_results['OutputStorageLevel'].reindex(range(1, len(idx) + 1)).fillna(0).set_index(idx)
        else:
            profiles = temp_results['OutputStorageLevel'].set_index(idx)
    else:
        logging.error('HydroScheduling parameter should be either "Regional" or "Zonal" (case sensitive). ')
        sys.exit()

    if mts_plot:
        profiles.plot()

    # Re-index to the main simulation time step:
    if config['SimulationTimeStep'] != temp_config['SimulationTimeStep']:
        profiles = profiles.reindex(idx_orig, method='nearest')
    pickle.dump(profiles, open("temp_profiles.p", "wb"))
    return profiles


def get_temp_sim_results(config, gams_dir=None):
    """
    This function reads the simulation environment folder once it has been solved and loads
    the input variables together with the results.

    :param config:              Read config file
    :param gams_dir:            Path to GAMS directory
    :returns results:           Two dictionaries with all the outputs
    """

    resultfile = config['SimulationDirectory'] + '/Results_simple.gdx'
    results = gdx_to_dataframe(gdx_to_list(gams_dir, resultfile, varname='all', verbose=True), fixindex=True,
                               verbose=True)
    results['OutputStorageLevel'] = results['OutputStorageLevel'].reindex(
        list(range(results['OutputStorageLevel'].index.min(),
                   results['OutputStorageLevel'].index.max() + 1)), fill_value=0)
    return results
