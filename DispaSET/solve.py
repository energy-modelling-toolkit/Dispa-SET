# -*- coding: utf-8 -*-
"""
This worksheet contains the two main functions to solve the DispaSET optimization problem using PYOMO or GAMS.


@author: 'Sylvain Quoilin'
"""

#######################################################################################################################
############################################ Dispa-SET: main model ####################################################
#######################################################################################################################


import sys
import os
import shutil
import logging
import time

from .misc.gdx_handler import get_gams_path, import_local_lib, package_exists

def is_sim_folder_ok(sim_folder):

    if not os.path.exists(sim_folder):
        logging.error('The provided DispaSET simulation environment folder (' + sim_folder + ') does not exist')
        return False

    if not os.path.exists(os.path.join(sim_folder, 'Inputs.gdx')):
        logging.error('There is no Inputs.gdx file within the specified DispaSET simulation environment folder (' + sim_folder + '). Check that the GDX output is activated in the option file and that no error stated during the pre-processing')
        return False

    if not os.path.exists(os.path.join(sim_folder, 'UCM_h.gms')):
        logging.error('There is no UCM_h.gms file within the specified DispaSET simulation environment folder (' + sim_folder + ')')
        return False
    return True


def solve_GAMS(sim_folder, gams_folder=None, output_lst=False):
    if not package_exists('gams'):
        logging.warning('Could not import gams. Trying to automatically locate gdxcc folder')
        if not import_local_lib('gams'):
            return False

    if not os.path.exists(gams_folder):
        logging.warn('The provided path for GAMS (' + gams_folder + ') does not exist. Trying to locate...')
        gams_folder = get_gams_path()
        if not os.path.exists(gams_folder):
            logging.error('GAMS path cannot be located. Simulation is stopped')
            return False

    if is_sim_folder_ok(sim_folder):
        # create GAMS workspace:
        from gams import GamsWorkspace
        ws = GamsWorkspace(system_directory=gams_folder, debug=3)
        shutil.copy(os.path.join(sim_folder, 'UCM_h.gms'), ws.working_directory)
        shutil.copy(os.path.join(sim_folder, 'Inputs.gdx'), ws.working_directory)
        t1 = ws.add_job_from_file('UCM_h.gms')
        opt = ws.add_options()
        #Do not create .lst file
        if not output_lst:
            if sys.platform == 'win32':
                opt.output = 'nul'
            else:
                opt.output = '/dev/null'
        time0 = time.time()
        t1.run(opt)
        logging.info('Completed simulation in {0:.2f} seconds'.format(time.time() - time0))

        # copy the result file to the simulation environment folder:
        shutil.copy(os.path.join(ws.working_directory, 'Results.gdx'), sim_folder)
        if output_lst:
            shutil.copy(os.path.join(ws.working_directory, 'UCM_h.lst'), sim_folder)
        return t1
    else:
        return False


def solve_pyomo(sim_folder):
    import pickle
    from .pyomo.model import DispaSolve

    with open(os.path.join(sim_folder, 'Inputs.p'), 'rb') as pyomo_input_file:
        SimData = pickle.load(pyomo_input_file)
    if SimData['config']['SimulationType'] == 'MILP':
        LPFormulation = False
    else:
        LPFormulation = True
    time0 = time.time()
    results = DispaSolve(SimData['sets'], SimData['parameters'], LPFormulation=LPFormulation)
    logging.info('Completed simulation in {0:.2f} seconds'.format(time.time() - time0))

    return results
