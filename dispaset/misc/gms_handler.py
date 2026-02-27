import time
import logging
import sys
import os
import shutil

from .str_handler import force_str

def solve_high_level(gams_folder, sim_folder, gams_file='UCM_h.gms', result_file='Results.gdx', input_gdx_file='Inputs.gdx', output_lst=False):
    """Use higher level apis to run GAMS

    :param gams_folder: Path to GAMS installation
    :param sim_folder: Path to simulation folder
    :param gams_file: Name of the GAMS file to run (e.g., UCM_h.gms or UCM_MTS.gms)
    :param result_file: Name of the GAMS result GDX file (e.g., Results.gdx or Results_MTS.gdx)
    :param input_gdx_file: Name of the GAMS input GDX file (e.g., Inputs.gdx or Inputs_MTS.gdx)
    :param output_lst: Whether to keep the GAMS .lst file
    """
    # create GAMS workspace:
    gams_folder = force_str(gams_folder)
    from gams import GamsWorkspace
    ws = GamsWorkspace(working_directory=str(sim_folder), system_directory=str(gams_folder), debug=2)
    t1 = ws.add_job_from_file(gams_file)
    opt = ws.add_options()
    opt.solprint = 0  # Suppress solution printout
    opt.logoption = 0  # Minimize logging
    time0 = time.time()
    status = t1.run(opt)
    # removing _gams_py files
    for filename in os.listdir(sim_folder):
        if filename.startswith("_gams_py"):
            file_path = os.path.join(sim_folder, filename)
            if os.path.isfile(file_path):
                try:
                    os.remove(file_path)
                except PermissionError:
                    print(f"Warning: Could not delete {file_path} because it is in use.")
    logging.info('Completed simulation in {0:.2f} seconds'.format(time.time() - time0))
    return status