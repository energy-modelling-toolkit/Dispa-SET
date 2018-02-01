import time
import logging
import sys
import os
import shutil

from .str_handler import force_str

def solve_high_level(gams_folder,sim_folder,output_lst=False):
    """Use higher level apis to run GAMSy"""
    # create GAMS workspace:
    try:
        from gams import GamsWorkspace
        ws = GamsWorkspace(system_directory=str(gams_folder), debug=3)
        shutil.copy(os.path.join(sim_folder, 'UCM_h.gms'), ws.working_directory)
        shutil.copy(os.path.join(sim_folder, 'Inputs.gdx'), ws.working_directory)
        t1 = ws.add_job_from_file('UCM_h.gms')
        opt = ws.add_options()
        # Do not create .lst file
        if not output_lst:
            if sys.platform == 'win32':
                opt.output = 'nul'
            else:
                opt.output = '/dev/null'
        time0 = time.time()
        status = t1.run(opt)
    except Exception as e:
        if 'optCreateD' in str(e):
            logging.error('The GAMS solver can only be run once in the same console. Please open another console')
            sys.exit(1)
        else:
            logging.error('The following error occured when trying to solve the model in gams: ' + str(e))
            sys.exit(1)
    # copy the result file to the simulation environment folder:
    shutil.copy(os.path.join(ws.working_directory, 'Results.gdx'), sim_folder)
    for filename in ['UCM_h.lst','UCM_h.log','debug.gdx']:
        if os.path.isfile(os.path.join(ws.working_directory, filename)):
            shutil.copy(os.path.join(ws.working_directory, filename), sim_folder)
    logging.info('Completed simulation in {0:.2f} seconds'.format(time.time() - time0))
    return status


def solve_low_level(gams_folder,sim_folder,output_lst=False,logoption=3):
    """Use lower level apis to run GAMS. BAsed on GAMS xpexample2.py"""
    from gamsxcc import new_gamsxHandle_tp, gamsxRunExecDLL, gamsxFree, gamsxCreateD, GMS_SSSIZE
    from optcc import new_optHandle_tp, optReadDefinition, optSetStrStr, optSetIntStr, optHandleToPtr, optFree, \
        optCreateD

    sim_folder = force_str(sim_folder)
    gams_folder = force_str(gams_folder)
    model = os.path.abspath(os.path.join(sim_folder, 'UCM_h.gms'))

    def callGams(gamsxHandle, optHandle, sysDir, model):
        deffile = force_str(sysDir + u'/optgams.def')

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
            return False

        return True

    optHandle = new_optHandle_tp()
    gamsxHandle = new_gamsxHandle_tp()

    try:
        rc = gamsxCreateD(gamsxHandle, gams_folder, GMS_SSSIZE)
        rc = optCreateD(optHandle, gams_folder, GMS_SSSIZE)
        time0 = time.time()
        status = callGams(gamsxHandle, optHandle, gams_folder, model)
    except Exception as e:
        logging.error('The following error occured when trying to solve the model in gams: ' + str(e))
        sys.exit(1)
    finally:
        optFree(optHandle)
        gamsxFree(gamsxHandle)
    logging.info('Completed simulation in {0:.2f} seconds'.format(time.time() - time0))
    return status
