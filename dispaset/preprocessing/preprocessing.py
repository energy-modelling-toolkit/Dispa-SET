# -*- coding: utf-8 -*-
"""
This is the main file of the DispaSET pre-processing tool.
It comprises a single function that generates the DispaSET simulation environment.

@author: S. Quoilin
"""
import datetime as dt
import logging
import sys
import os, shutil

import pandas as pd
import numpy as np

from .build import build_single_run

from ..solve import solve_GAMS
from ..misc.gdx_handler import gdx_to_dataframe, gdx_to_list
from ..postprocessing.data_handler import GAMSstatus
from .utils import pd_timestep
from ..common import commons


def build_simulation(config, mts_plot=None, MTSTimeStep=24):
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
        if config['H2FlexibleDemand'] != '':
            [new_profiles, new_PtLDemand] = mid_term_scheduling(config, mts_plot=mts_plot, TimeStep=MTSTimeStep)
            # Build simulation data with new profiles
            logging.info('\n\nBuilding final simulation\n')
            SimData = build_single_run(config, new_profiles, new_PtLDemand)
        else:
            new_profiles = mid_term_scheduling(config, mts_plot=mts_plot, TimeStep=MTSTimeStep)
            # Build simulation data with new profiles
            logging.info('\n\nBuilding final simulation\n')
            SimData = build_single_run(config, new_profiles)
        
    # Copy the log to the simulation folder:
    if os.path.isfile(commons['logfile']):
        shutil.copy(commons['logfile'], os.path.join(config['SimulationDirectory'], commons['logfile']))
    else:
        logging.error('Could not find log file in current directory')
    return SimData

def _check_results(results):
    '''
    Function that checks the gams status in the results
    '''
    if "model" in results['status']:
        errors = results['status'][(results['status']['model'] != 1) & (results['status']['model'] != 8)]
        if len(errors) > 0:
            logging.critical('Some simulation errors were encountered when running the regional MTS. Some results could not be computed, for example at  time ' + str(errors.index[0]) + ', with the error message: "' + GAMSstatus('model', errors['model'].iloc[0]) + '"')
            for i in errors.index:
                errors.loc[i,'Error Message'] = GAMSstatus('model',errors['model'][i])
            sys.exit(1)
    return True

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
        logging.info('MTS Simulation with all zones selected')
    else:
        temp_config['zones'] = config['mts_zones']

    # remove look ahead:
    temp_config['LookAhead'] = 0
    
    # Don't use any historical reservoir level:
    temp_config['ReservoirLevels'] = ''
    
    # use a LP formulation
    temp_config['SimulationType'] = 'LP clustered'

    # Adjust time step:
    if TimeStep is not None:
        # Indexes of the MTS simulation:
        idx = pd.date_range(start=dt.datetime(*temp_config['StartDate']),
                            end=dt.datetime(*temp_config['StopDate']),
                            freq=pd_timestep(TimeStep)).tz_localize(None)
        temp_config['SimulationTimeStep'] = TimeStep
        gams_file = 'UCM_h.gms'
        temp_config['HorizonLength'] = (idx[-1] - idx[0]).days+1
        resultfile = 'Results.gdx'
    else:
        idx = pd.date_range(start=dt.datetime(*temp_config['StartDate']),
                            end=dt.datetime(*temp_config['StopDate']),
                            freq=pd_timestep(temp_config['SimulationTimeStep'])).tz_localize(None)
        gams_file = 'UCM_h_simple.gms'
        resultfile = 'Results_simple.gdx'

    # Checking which type of hydro scheduling simulation is specified in the config file:
    if config['HydroScheduling'] == 'Zonal':
        no_of_zones = len(config['mts_zones'])
        temp_results = {}
        profiles = pd.DataFrame(index=idx)
        PtLDemand = pd.DataFrame(index=idx)
        for i, c in enumerate(config['mts_zones']):
            logging.info('\n\nLaunching Mid-Term Scheduling for zone '+ c + ' (Number ' + str(i + 1) + ' out of ' + str(no_of_zones) + ')\n')
            temp_config['zones'] = [c]  # Override zone that needs to be simulated
            SimData = build_single_run(temp_config, MTS=1)  # Create temporary SimData
            units = SimData['units']
            r = solve_GAMS(sim_folder=temp_config['SimulationDirectory'],
                           gams_folder=temp_config['GAMS_folder'],
                           gams_file=gams_file,
                           result_file=resultfile)
            temp_results[c] = gdx_to_dataframe(
                gdx_to_list(config['GAMS_folder'], config['SimulationDirectory'] + '/' + resultfile, varname='all',
                            verbose=True), fixindex=True, verbose=True)
            _check_results(temp_results[c])
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
                temp_results[c]['OutputStorageLevel'] = temp_results[c]['OutputStorageLevel'].reindex(range(1, len(idx) + 1)).fillna(0).values                
            for u in temp_results[c]['OutputStorageLevel']:
                if u not in units.index:
                    logging.critical('Unit "' + u + '" is reported in the reservoir levels of the result file but does not appear in the units table')
                    sys.exit(1)
                for u_old in units.loc[u,'FormerUnits']:
                    profiles[u_old] = temp_results[c]['OutputStorageLevel'][u].values
            if config['H2FlexibleDemand'] != '':
                if 'OutputPtLDemand' not in temp_results[c]:
                    logging.critical('PtL demand in zone ' + c + ' was not computed')
                    sys.exit(0)
                elif len(temp_results[c]['OutputSPtLDemand']) > len(idx):
                    logging.critical('The number of time steps in the mid-term simulation results (' + str(
                                     len(temp_results[c]['OutputPtLDemand'])) +
                                     ') does not match the length of the index (' + str(len(idx)) + ')')
                    sys.exit(0)                            
                elif len(temp_results[c]['OutputPtLDemand']) < len(idx):
                        temp_results[c]['OutputPtLDemand'] = temp_results[c]['OutputPtLDemand'].reindex(range(1, len(idx) + 1)).fillna(0).values
                for u in temp_results[c]['OutputPtLDemand']:
                    if u not in units.index:
                        logging.critical('Unit "' + u + '" is reported in the PtL demand of the result file but does not appear in the units table')
                        sys.exit(1)
                    for u_old in units.loc[u,'FormerUnits']:
                        PtLDemand[u_old] = temp_results[c]['OutputPtLDemand'][u].values  
                
    # Solving reservoir levels for all regions simultaneously
    elif config['HydroScheduling'] == 'Regional':
        logging.info('\n\nLaunching regional Mid-Term Scheduling \n')
        SimData = build_single_run(temp_config, MTS=1)  # Create temporary SimData
        units = SimData['units']
        r = solve_GAMS(sim_folder=temp_config['SimulationDirectory'],
                        gams_folder=temp_config['GAMS_folder'],
                        gams_file=gams_file,
                        result_file=resultfile)
        temp_results = gdx_to_dataframe(
            gdx_to_list(config['GAMS_folder'], config['SimulationDirectory'] + '/' + resultfile, varname='all',
                        verbose=True), fixindex=True, verbose=True)
        _check_results(temp_results)
        if 'OutputStorageLevel' not in temp_results:
            logging.critical('Storage levels in the selected region were not computed, please check that storage units '
                             'are present in the power plant database! If not, unselect zones with no storage units '
                             'form the zones in the MTS module')
            sys.exit(0)
        if len(temp_results['OutputStorageLevel']) < len(idx):
            profiles = temp_results['OutputStorageLevel'].reindex(range(1, len(idx) + 1)).fillna(0).set_index(idx)
        else:
            profiles = temp_results['OutputStorageLevel'].set_index(idx)
        
        if config['H2FlexibleDemand'] != '':
            if 'OutputPtLDemand' not in temp_results:
                logging.critical('PtL Demand in the selected region was not computed')
                sys.exit(0)   
                if len(temp_results['OutputPtLDemand']) < len(idx):
                    PtLDemand = temp_results['OutputPtLDemand'].reindex(range(1, len(idx) + 1)).fillna(0).set_index(idx)
                else:
                    PtLDemand = temp_results['OutputPtLDemand'].set_index(idx)
            
        # Updating the profiles table with the original unit names:
        for u in profiles:
            if u not in units.index:
                logging.critical('Unit "' + u + '" is reported in the reservoir levels of the result file but does not appear in the units table')
                sys.exit(1)
            for u_old in units.loc[u,'FormerUnits']:
                profiles[u_old] = profiles[u]
            profiles.drop(u,axis=1,inplace=True)
        if config['H2FlexibleDemand'] != '':
            for u in PtLDemand:
                if u not in units.index:
                    logging.critical('Unit "' + u + '" is reported in the PtL demand of the result file but does not appear in the units table')
                    sys.exit(1)
                    for u_old in units.loc[u,'FormerUnits']:
                        PtLDemand[u_old] = PtLDemand[u]
                        PtLDemand.drop(u,axis=1,inplace=True)            
    else:
        logging.error('HydroScheduling parameter should be either "Regional" or "Zonal" (case sensitive). ')
        sys.exit()
    
    # replace all 1.000000e+300 values by nan since they correspond to undefined in GAMS:
    profiles[profiles>=1E300] =  np.nan
    if config['H2FlexibleDemand'] != '':
        PtLDemand[PtLDemand>=1E300] =  np.nan

    if mts_plot:
        profiles.plot()

    # Re-index to the main simulation time step:
    if config['SimulationTimeStep'] != temp_config['SimulationTimeStep']:
        profiles = profiles.reindex(idx_orig, method='nearest')
        if config['H2FlexibleDemand'] != '':
            PtLDemand = PtLDemand.resample(pd_timestep(config['SimulationTimeStep'])).pad()
    pickle.dump(profiles, open(os.path.join(config['SimulationDirectory'],"temp_profiles.p"), "wb"))
    if config['H2FlexibleDemand'] != '':
        return profiles, PtLDemand
    else:
        return profiles

