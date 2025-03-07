import time
import logging
import sys
import os
import shutil

from .str_handler import force_str

def solve_high_level(gams_folder, sim_folder, gams_file='UCM_h.gms', result_file='Results.gdx', output_lst=False):
    """Use higher level apis to run GAMS"""
    # create GAMS workspace:
    gams_folder = force_str(gams_folder)
    from gams import GamsWorkspace
    ws = GamsWorkspace(system_directory=str(gams_folder), debug=3)
    shutil.copy(os.path.join(sim_folder, gams_file), ws.working_directory)
    shutil.copy(os.path.join(sim_folder, 'Inputs.gdx'), ws.working_directory)
    shutil.copy(os.path.join(sim_folder, 'cplex.opt'), ws.working_directory)
    t1 = ws.add_job_from_file(gams_file)
    opt = ws.add_options()
    # Do not create .lst file
    # if not output_lst:
    #     if sys.platform == 'win32':
    #         opt.output = 'nul'
    #     else:
    #         opt.output = '/dev/null'
    time0 = time.time()
    status = t1.run(opt)
    # copy the result file to the simulation environment folder:
    shutil.copy(os.path.join(ws.working_directory, result_file), sim_folder)
    gams_file_base = os.path.splitext(gams_file)[0]
    for filename in [gams_file_base + '.lst', gams_file_base + '.log', 'debug.gdx']:
        if os.path.isfile(os.path.join(ws.working_directory, filename)):
            shutil.copy(os.path.join(ws.working_directory, filename), sim_folder)
    logging.info('Completed simulation in {0:.2f} seconds'.format(time.time() - time0))
    return status