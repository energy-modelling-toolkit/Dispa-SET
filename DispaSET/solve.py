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


def solve_GAMS(sim_folder, gams_folder=None, work_dir=None, output_lst=False):
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
    sim_folder = sim_folder.encode()
    gams_folder = gams_folder.encode()

    if is_sim_folder_ok(sim_folder):
        # create GAMS workspace:
        from gams import GamsWorkspace
        try:
            ws = GamsWorkspace(system_directory=gams_folder, working_directory=work_dir, debug=3)
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
        except Exception as e:
            if 'optCreateD' in str(e):
                logging.error('The GAMS solver can only be run once in the same console. Please open another console')
                sys.exit(1)
            else:
                logging.error('The following error occured when trying to solve the model in gams: ' + str(e))
                sys.exit(1)
        logging.info('Completed simulation in {0:.2f} seconds'.format(time.time() - time0))

        # copy the result file to the simulation environment folder:
        shutil.copy(os.path.join(ws.working_directory, 'Results.gdx'), sim_folder)
        for filename in ['UCM_h.lst','UCM_h.log','debug.gdx']:
            if os.path.isfile(os.path.join(ws.working_directory, 'debug.gdx')):
                shutil.copy(os.path.join(ws.working_directory, 'debug.gdx'), sim_folder)
        if os.path.isfile(os.path.join(ws.working_directory, 'debug.gdx')):
            logging.warn('A debug file was created. There has probably been an optimization error')
        if os.path.isfile('warn.log'):
            shutil.copy('warn.log', os.path.join(sim_folder, 'warn_solve.log'))
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
    
    if os.path.isfile(SimData['config']['cplex_path']):
        path_cplex = SimData['config']['cplex_path']
    else:
        path_cplex = ''
        if len(SimData['config']['cplex_path']) > 2:
            logging.warn('The specified path for cplex (' + SimData['config']['cplex_path'] + ') is not valid. It will be ignored')
        
    time0 = time.time()
    results = DispaSolve(SimData['sets'], SimData['parameters'], LPFormulation=LPFormulation, path_cplex=path_cplex)
    logging.info('Completed simulation in {0:.2f} seconds'.format(time.time() - time0))
    if os.path.isfile('warn.log'):
        shutil.copy('warn.log', os.path.join(sim_folder, 'warn_solve.log'))
    return results
