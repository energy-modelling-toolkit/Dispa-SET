# -*- coding: utf-8 -*-
"""
This worksheet contains the main function to solve the Dispa-SET optimization problem using GAMS.

Solve with GAMS and the high level API
--------------------------------------
The high level interface is recommended for all users.

Installation:
    To install the high-level API in Python 2.x::

        cd gams24.4_linux_x64_64_sfx/apifiles/Python/api
        python gamssetup.py install

    To install the high-level API in Python 3.x::
    
        cd gams24.6_linux_x64_64_sfx/apifiles/Python/api_34
        python setup.py install    

"""

#######################################################################################################################
############################################ Dispa-SET: main model ####################################################
#######################################################################################################################


import os
import shutil
import logging
import time

from .misc.gdx_handler import get_gams_path, package_exists
from .misc.gms_handler import solve_high_level
from .common import commons


def is_sim_folder_ok(sim_folder):
    '''
    Function that checks if the provided path is a valid Dispa-SET simulation folder.
    The following files are required:

        - Inputs.gdx
        - UCM_h.gms

    :param sim_folder: path (absolute or relative) to the simulation folder
    '''
    if not os.path.exists(sim_folder):
        logging.error('The provided DispaSET simulation environment folder (' + sim_folder + ') does not exist')
        return False

    if not os.path.exists(os.path.join(sim_folder, u'Inputs.gdx')):
        logging.error(
            'There is no Inputs.gdx file within the specified DispaSET simulation environment folder (' + sim_folder + '). Check that the GDX output is activated in the option file and that no error stated during the pre-processing')
        return False

    # Check for either the dispatch GAMS file or the MTS GAMS file
    if not os.path.exists(os.path.join(sim_folder, u'UCM_h.gms')) and \
       not os.path.exists(os.path.join(sim_folder, u'UCM_MTS.gms')):
        logging.error(
            'Neither UCM_h.gms nor UCM_MTS.gms file found within the specified DispaSET simulation environment folder (' + sim_folder + ')')
        return False
    return True


def solve_GAMS(sim_folder, gams_folder=None, gams_file='UCM_h.gms', result_file='Results.gdx', output_lst=False):
    '''
    Function used to run the optimization using the GAMS engine.

    :param sim_folder: path to a valid Dispa-SET simulation folder
    :param gams_folder: optional path to the GAMS installation. If not provided, will try to detect automatically
    :param gams_file: name of the GAMS file to run (default: UCM_h.gms for dispatch, specify UCM_MTS.gms for MTS)
    :param result_file: name of the result file (default: Results.gdx)
    :param output_lst: Set to True to conserve a copy of the GAMS lst file in the simulation folder
    :return: True if simulation was successful, False otherwise
    '''
    if not package_exists('gams'):
        logging.error('GAMS API not found. Please install the GAMS Python API from your GAMS installation.')
        return False

    # Get GAMS path from provided path or try to locate it
    gams_folder = get_gams_path(gams_folder)
    if not gams_folder:
        logging.error('GAMS installation not found. Please set GAMSPATH or GAMSDIR environment variable or provide a valid path.')
        return False

    sim_folder = os.path.abspath(sim_folder)
    gams_folder = os.path.abspath(gams_folder)

    if is_sim_folder_ok(sim_folder):
        #Temporary warning for Spyder users:
        if any(['SPY_' in name for name in os.environ]): # check if spyder
            logging.info("\nIf the script seems stuck at this place \n(gams is optimizing but not output is displayed), \nit is preferable to run Dispa-SET in a \nseparate terminal (in Spyder: Preferences - Run - \nExecute in an external system terminal)")
        ret = solve_high_level(gams_folder, sim_folder, gams_file, result_file, output_lst=output_lst)
        if os.path.isfile(os.path.join(sim_folder, 'debug.gdx')):
            logging.warning('A debug file was created. There has probably been an optimization error')
        if os.path.isfile(commons['logfile']):
            shutil.copy(commons['logfile'], os.path.join(sim_folder, 'warn_solve.log'))
        return ret
    else:
        return False