import time
import logging
import sys
import os
import shutil

from .str_handler import force_str

def solve_high_level(gams_folder, sim_folder, gams_file='UCM_h.gms', result_file='Results.gdx', output_lst=False):
    """Use higher level apis to run GAMSy"""
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


def solve_low_level(gams_folder, sim_folder, gams_file='UCM_h.gms', result_file='Results.gdx', output_lst=False, logoption=3):
    """Use lower level apis to run GAMS. Based on GAMS xpexample2.py.
     We use the same signature as in solve_high_level for consistency when we call it.
     As we define the working directory to be the simulation directory we do not need to define the output names
     and move them around the filesystem"""
    from gamsxcc import new_gamsxHandle_tp, gamsxRunExecDLL, gamsxFree, gamsxCreateD, GMS_SSSIZE
    from optcc import new_optHandle_tp, optReadDefinition, optSetStrStr, optSetIntStr, optHandleToPtr, optFree, \
        optCreateD

    sim_folder = force_str(sim_folder)
    gams_folder = force_str(gams_folder)
    model = os.path.abspath(os.path.join(sim_folder, gams_file))

    def callGams(gamsxHandle, optHandle, sysDir, model):
        deffile = force_str(os.path.join(sysDir, 'optgams.def'))

        if optReadDefinition(optHandle, deffile):
            logging.error("*** Error ReadDefinition, cannot read def file:" + deffile)
            return False

        optSetStrStr(optHandle, "SysDir", sysDir)
        optSetStrStr(optHandle, "WorkDir", sim_folder)
        optSetStrStr(optHandle, "Input", model)
        optSetIntStr(optHandle, "LogOption", logoption)
        ret = gamsxRunExecDLL(gamsxHandle, optHandleToPtr(optHandle), sysDir, 1)
        if ret[0] != 0:
            logging.error("*** Error RunExecDLL: Error in GAMS call = " + str(ret[1]))
            if 'White space' in str(ret[1]):
                logging.error("The Unix GAMS API does not accept white spaces. Move dispaset to a folder that does not contain white spaces")
            return False

        return True

    optHandle = new_optHandle_tp()
    gamsxHandle = new_gamsxHandle_tp()

    try:
        __ = gamsxCreateD(gamsxHandle, gams_folder, GMS_SSSIZE)
        __ = optCreateD(optHandle, gams_folder, GMS_SSSIZE)
        time0 = time.time()
        status = callGams(gamsxHandle, optHandle, gams_folder, model)
    except Exception as e:
        logging.error('The following error occured when trying to solve the model in gams: ' + str(e))
        return False
    finally:
        optFree(optHandle)
        gamsxFree(gamsxHandle)
    logging.info('Completed simulation in {0:.2f} seconds'.format(time.time() - time0))
    return status