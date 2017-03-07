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

from .misc.gdx_handler import get_gams_path

def solve_GAMS(sim_folder, gams_folder=None):
    try:
        import gams
    except ImportError:
        logging.critical('The gams library is required to run the GAMS versions of DispaSET.'
                 'Please install if from the /apifiles/Python/api/ folder in the GAMS directory')
        sys.exit(1)

    if not os.path.exists(gams_folder):
        gams_folder = get_gams_path()
        if not os.path.exists(gams_folder):
            logging.critical('The provided path for GAMS (' + gams_folder + ') does not exist')
            sys.exit(1)

    if not os.path.exists(sim_folder):
        logging.critical('The provided DispaSET simulation environment folder (' + sim_folder + ') does not exist')
        sys.exit(1)

    if not os.path.exists(sim_folder + '/Inputs.gdx'):
        logging.critical('There is no Inputs.gdx file within the specified DispaSET simulation environment folder (' + sim_folder + ')')
        sys.exit(1)

    if not os.path.exists(sim_folder + '/UCM_h.gms'):
        logging.critical('There is no UCM_h.gms file within the specified DispaSET simulation environment folder (' + sim_folder + ')')
        sys.exit(1)

    # create GAMS workspace:
    ws = gams.GamsWorkspace(system_directory=gams_folder, debug=3)
    shutil.copy(os.path.join(sim_folder, 'UCM_h.gms'), ws.working_directory)
    shutil.copy(os.path.join(sim_folder, 'Inputs.gdx'), ws.working_directory)
    t1 = ws.add_job_from_file('UCM_h.gms')
    time0 = time.time()
    t1.run()
    logging.info('Completed simulation in {0:.2f} seconds'.format(time.time() - time0))

    # copy the result file to the simulation environment folder:
    shutil.copy(os.path.join(ws.working_directory, 'Results.gdx'), sim_folder)
    return t1
#    for rec in t1.out_db["x"]:
#        print "x(" + rec.key(0) + "," + rec.key(1) + "): level=" + str(rec.level) + " marginal=" + str(rec.marginal)


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
