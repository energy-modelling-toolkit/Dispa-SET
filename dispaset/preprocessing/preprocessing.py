# -*- coding: utf-8 -*-
"""
This is the main file of the DispaSET pre-processing tool.
It comprises a single function that generates the DispaSET simulation environment.

@author: S. Quoilin
"""
import datetime as dt
import logging
import os
import shutil
import sys

import numpy as np
import pandas as pd

from .build import build_single_run
from .utils import pd_timestep
from ..common import commons
from ..misc.gdx_handler import gdx_to_dataframe, gdx_to_list
from ..postprocessing.data_handler import GAMSstatus
from ..solve import solve_GAMS


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
    elif 'SectorXFlexibleDemand' and 'SectorXFlexibleSupply' in config:
        if (config['SectorXFlexibleDemand'] != '') and (config['SectorXFlexibleSupply'] == ''):
            [new_profiles, new_SectorXFlexDemand, new_profilesSectorX] = mid_term_scheduling(config, mts_plot=mts_plot,
                                                                                              TimeStep=MTSTimeStep)
            # Build simulation data with new profiles
            logging.info('\n\nBuilding final simulation\n')
            SimData = build_single_run(config, profiles=new_profiles, SectorXFlexDemand=new_SectorXFlexDemand,
                                        profilesSectorX=new_profilesSectorX)
        elif (config['SectorXFlexibleSupply'] != '') and (config['SectorXFlexibleDemand'] == ''):
            [new_profiles, new_SectorXFlexSupply, new_profilesSectorX] = mid_term_scheduling(config, mts_plot=mts_plot,
                                                                                              TimeStep=MTSTimeStep)
            # Build simulation data with new profiles
            logging.info('\n\nBuilding final simulation\n')
            SimData = build_single_run(config, profiles=new_profiles, SectorXFlexSupply=new_SectorXFlexSupply,
                                        profilesSectorX=new_profilesSectorX)
        elif (config['SectorXFlexibleSupply'] != '') and (config['SectorXFlexibleDemand'] != ''):
            [new_profiles, new_SectorXFlexDemand, new_SectorXFlexSupply, new_profilesSectorX] = mid_term_scheduling(
                config,
                mts_plot=mts_plot,
                TimeStep=MTSTimeStep)
            # Build simulation data with new profiles
            logging.info('\n\nBuilding final simulation\n')
            SimData = build_single_run(config, profiles=new_profiles, SectorXFlexDemand=new_SectorXFlexDemand,
                                        SectorXFlexSupply=new_SectorXFlexSupply, profilesSectorX=new_profilesSectorX)
        else:
            [new_profiles, new_profilesSectorX] = mid_term_scheduling(config, mts_plot=mts_plot, TimeStep=MTSTimeStep)
            # Build simulation data with new profiles
            logging.info('\n\nBuilding final simulation\n')
            SimData = build_single_run(config, profiles=new_profiles, profilesSectorX=new_profilesSectorX)
    else:
        [new_profiles, new_profilesSectorX] = mid_term_scheduling(config, mts_plot=mts_plot, TimeStep=MTSTimeStep)
        # Build simulation data with new profiles
        logging.info('\n\nBuilding final simulation\n')
        SimData = build_single_run(config, profiles=new_profiles, profilesSectorX=new_profilesSectorX)

    # Copy the log to the simulation folder:
    if os.path.isfile(commons['logfile']):
        shutil.copy(commons['logfile'], os.path.join(config['SimulationDirectory'], commons['logfile']))
    else:
        logging.error('Could not find log file in current directory')
    return SimData


def _check_results(results):
    """
    Function that checks the gams status in the results
    """
    if "model" in results['status']:
        errors = results['status'][(results['status']['model'] != 1) & (results['status']['model'] != 8)]
        if len(errors) > 0:
            logging.critical('Some simulation errors were encountered when running the regional MTS. '
                             'Some results could not be computed, for example at  time ' + str(errors.index[0]) +
                             ', with the error message: "' + GAMSstatus('model', errors['model'].iloc[0]) + '"')
            for i in errors.index:
                errors.loc[i, 'Error Message'] = GAMSstatus('model', errors['model'][i])
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
    config['StopDate'] = (y_end, m_end, d_end, 23, 59, 00)  # updating stopdate to the end of the day
    config['idx'] = pd.date_range(start=dt.datetime(*config['StartDate']),
                                  end=dt.datetime(*config['StopDate']),
                                  freq=commons['TimeStep']).tz_localize(None)
    # Indexes including the last look-ahead period
    enddate_long = config['idx'][-1] + dt.timedelta(days=config['LookAhead'])
    idx_long = pd.date_range(start=config['idx'][0], end=enddate_long, freq=commons['TimeStep'])

    # New configuration dictionary, specific to MTS:
    temp_config = dict(config)

    # Dates for the mid term scheduling
    if config['HydroSchedulingHorizon'] == 'Annual':
        temp_config['StartDate'] = (y_start, 1, 1, 00, 00, 00)  # updating start date to the beginning of the year
        temp_config['StopDate'] = (y_start, 12, 31, 23, 59, 00)  # updating stopdate to the end of the year
        logging.info('Hydro scheduling is performed for the period between 01.01.' + str(y_start) + ' and 12.31.' +
                     str(y_start))
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
    if config['InitialFinalReservoirLevel'] != 0:
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
        temp_config['HorizonLength'] = (idx[-1] - idx[0]).days + 1
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
        SectorXFlexDemand = pd.DataFrame(index=idx)
        SectorXFlexSupply = pd.DataFrame(index=idx)
        profilesSectorX = pd.DataFrame()
        for i, c in enumerate(config['mts_zones']):
            logging.info(
                '\n\nLaunching Mid-Term Scheduling for zone ' + c + ' (Number ' + str(i + 1) + ' out of ' + str(
                    no_of_zones) + ')\n')
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
                temp_results[c]['OutputStorageLevel'] = temp_results[c]['OutputStorageLevel'].reindex(
                    range(1, len(idx) + 1)).fillna(0).values
            for u in temp_results[c]['OutputStorageLevel']:
                if u not in units.index:
                    logging.critical('Unit "' + u + '" is reported in the reservoir levels of the result file but '
                                                    'does not appear in the units table')
                    sys.exit(1)
                for u_old in units.loc[u, 'FormerUnits']:
                    profiles[u_old] = temp_results[c]['OutputStorageLevel'][u].values

            # for nx in temp_results[c]['OutputSectorXStorageLevel']:
            #     if nx not in units.index:
            #         logging.critical('Unit "' + u + '" is reported in the reservoir levels of the result file but '
            #                                         'does not appear in the units table')
            #         sys.exit(1)
            #     for u_old in units.loc[u, 'FormerUnits']:
            #         profiles[u_old] = temp_results[c]['OutputStorageLevel'][u].values

            if config['SectorXFlexibleDemand'] != '':
                if 'OutputSectorXFlexDemand' not in temp_results[c]:
                    logging.critical('BS Flex demand in zone ' + c + ' was not computed')
                    sys.exit(0)
                elif len(temp_results[c]['OutputSectorXFlexDemand']) > len(idx):
                    logging.critical('The number of time steps in the mid-term simulation results (' + str(
                        len(temp_results[c]['OutputSectorXFlexDemand'])) +
                                     ') does not match the length of the index (' + str(len(idx)) + ')')
                    sys.exit(0)
                elif len(temp_results[c]['OutputSectorXFlexDemand']) < len(idx):
                    temp_results[c]['OutputSectorXFlexDemand'] = temp_results[c]['OutputSectorXFlexDemand'].reindex(
                        range(1, len(idx) + 1)).fillna(0)

            if config['SectorXFlexibleSupply'] != '':
                if 'OutputSectorXFlexSupply' not in temp_results[c]:
                    logging.critical('BS Flex demand in zone ' + c + ' was not computed')
                    sys.exit(0)
                elif len(temp_results[c]['OutputSectorXFlexSupply']) > len(idx):
                    logging.critical('The number of time steps in the mid-term simulation results (' + str(
                        len(temp_results[c]['OutputSectorXFlexSupply'])) +
                                     ') does not match the length of the index (' + str(len(idx)) + ')')
                    sys.exit(0)
                elif len(temp_results[c]['OutputSectorXFlexSupply']) < len(idx):
                    temp_results[c]['OutputSectorXFlexSupply'] = temp_results[c]['OutputSectorXFlexSupply'].reindex(
                        range(1, len(idx) + 1)).fillna(0)

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
        # Marco: create an empty df for profilesSectorX
        profilesSectorX = pd.DataFrame()
        if 'OutputStorageLevel' not in temp_results:
            logging.critical('Storage levels in the selected region were not computed, please check that storage units '
                             'are present in the power plant database! If not, unselect zones with no storage units '
                             'form the zones in the MTS module')
            sys.exit(0)
        if len(temp_results['OutputStorageLevel']) < len(idx):
            profiles = temp_results['OutputStorageLevel'].reindex(range(1, len(idx) + 1)).fillna(0).set_index(idx)
        else:
            profiles = temp_results['OutputStorageLevel'].set_index(idx)
        # Marco: New "else" condition added, check the name of the parameters 'SectorXReservoirLevels' different than 'OutputSectorXStorageLevel'
            
        if 'SectorXReservoirLevels' in config and config['SectorXReservoirLevels'] != '':
            if len(temp_results['OutputSectorXStorageLevel']) < len(idx):
                profilesSectorX = temp_results['OutputSectorXStorageLevel'].reindex(range(1, len(idx) + 1)).fillna(
                0).set_index(idx)
            else:
                profilesSectorX = temp_results['OutputSectorXStorageLevel'].set_index(idx)
        else: 
            logging.info('BS Storage Sectors were not computed')
            profilesSectorX = temp_results['OutputStorageLevel'].reindex(range(1, len(idx) + 1)).fillna(
            0).set_index(idx)

        if 'SectorXFlexibleDemand' in config and config['SectorXFlexibleDemand'] != '':
            if len(temp_results['OutputSectorXFlexDemand']) < len(idx):
                SectorXFlexDemand = temp_results['OutputSectorXFlexDemand'].reindex(range(1, len(idx) + 1)).fillna(
                    0).set_index(
                    idx)
            else:
                SectorXFlexDemand = temp_results['OutputSectorXFlexDemand'].set_index(idx)
        else:
            logging.info('BS Flexible Demand Sectors were not computed')
            SectorXFlexDemand = temp_results['OutputCommitted'].reindex(range(1, len(idx) + 1)).fillna(
            0).set_index(idx)
        

        if 'SectorXFlexibleSupply' in config and config['SectorXFlexibleSupply'] != '':
            if len(temp_results['OutputSectorXFlexSupply']) < len(idx):
                SectorXFlexSupply = temp_results['OutputSectorXFlexSupply'].reindex(range(1, len(idx) + 1)).fillna(
                    0).set_index(
                    idx)
            else:
                SectorXFlexSupply = temp_results['OutputSectorXFlexSupply'].set_index(idx)
        else:
            logging.info('BS Flexible Supply Sectors were not computed')
            SectorXFlexSupply = temp_results['OutputCommitted'].reindex(range(1, len(idx) + 1)).fillna(
            0).set_index(idx)

        # Updating the profiles table with the original unit names:
        for u in profiles:
            if u not in units.index:
                logging.critical('Unit "' + u + '" is reported in the reservoir levels of the result file but '
                                                'does not appear in the units table')
                sys.exit(1)
            for u_old in units.loc[u, 'FormerUnits']:
                profiles[u_old] = profiles[u]
            profiles.drop(u, axis=1, inplace=True)
    else:
        logging.error('HydroScheduling parameter should be either "Regional" or "Zonal" (case sensitive). ')
        sys.exit()

    # replace all 1.000000e+300 values by nan since they correspond to undefined in GAMS:
    profiles[profiles >= 1E300] = np.nan
    profilesSectorX[profilesSectorX >= 1E300]  = np.nan
    if 'SectorXFlexibleDemand' in config and config['SectorXFlexibleDemand'] != '':
        SectorXFlexDemand[SectorXFlexDemand >= 1E300] = np.nan

    if 'SectorXFlexibleSupply' in config and config['SectorXFlexibleSupply'] != '':
        SectorXFlexSupply[SectorXFlexSupply >= 1E300] = np.nan

    if mts_plot:
        profiles.plot()

    # Copy results from pre-processing
    sim_folder = config['SimulationDirectory']
    shutil.copyfile(os.path.join(sim_folder, 'Results.gdx'),
                    os.path.join(sim_folder, 'Results_MTS.gdx'))
    shutil.copyfile(os.path.join(sim_folder, 'Inputs.gdx'),
                    os.path.join(sim_folder, 'Inputs_MTS.gdx'))

    # Re-index to the main simulation time step:

    if config['SimulationTimeStep'] != temp_config['SimulationTimeStep']:
        profiles = profiles.reindex(idx_long, method='nearest')
        #Marco: Missing condition if profilesSectorX dont exist
        # profilesSectorX = profilesSectorX.reindex(idx_long, method='nearest')

        temp_config['StartDate'] = (y_start, 1, 1, 00, 00, 00)  # updating start date to the beginning of the year
        temp_config['StopDate'] = (y_start + 1, 1, 2, 00, 59, 00)
        idx_tmp = pd.date_range(start=dt.datetime(*temp_config['StartDate']),
                                end=dt.datetime(*temp_config['StopDate']),
                                freq=pd_timestep(TimeStep)).tz_localize(None)

        if 'SectorXFlexibleDemand' in config and config['SectorXFlexibleDemand'] != '':
            SectorXFlexDemand = pd.DataFrame(SectorXFlexDemand, index=idx_tmp).fillna(0)
            SectorXFlexDemand = SectorXFlexDemand.resample(pd_timestep(config['SimulationTimeStep'])).ffill()
            SectorXFlexDemand = SectorXFlexDemand.loc[idx_long, :]

        if 'SectorXFlexibleSupply' in config and config['SectorXFlexibleSupply'] != '':
            SectorXFlexSupply = pd.DataFrame(SectorXFlexSupply, index=idx_tmp).fillna(0)
            SectorXFlexSupply = SectorXFlexSupply.resample(pd_timestep(config['SimulationTimeStep'])).ffill()
            SectorXFlexSupply = SectorXFlexSupply.loc[idx_long, :]

    pickle.dump(profiles, open(os.path.join(config['SimulationDirectory'], "temp_profiles.p"), "wb"))
    # pickle.dump(profilesSectorX, open(os.path.join(config['SimulationDirectory'], "temp_profilesSectorX.p"), "wb"))
    if 'SectorXFlexibleDemand' and 'SectorXFlexibleSupply' in config:
        if (config['SectorXFlexibleSupply'] == '') and (config['SectorXFlexibleDemand'] != ''):
            return profiles, SectorXFlexDemand, profilesSectorX
        elif (config['SectorXFlexibleSupply'] != '') and (config['SectorXFlexibleDemand'] == ''):
            return profiles, SectorXFlexSupply, profilesSectorX
        elif (config['SectorXFlexibleSupply'] != '') and (config['SectorXFlexibleDemand'] != ''):
            return profiles, SectorXFlexDemand, SectorXFlexSupply, profilesSectorX
        else:
            return profiles, profilesSectorX
    else:
        return profiles, profilesSectorX
