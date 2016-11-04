# -*- coding: utf-8 -*-
"""
This worksheet contains the two main functions to solve the Dispa-SET optimization problem using PYOMO.

Functions: 
- DispaSolve: Implements the rolling horizon
- DispOptim: Performs the optimization for each horizon

@author: 'Sylvain Quoilin'
"""

#######################################################################################################################
############################################ Dispa-SET: main model ####################################################
#######################################################################################################################


import sys
import os
import shutil

import numpy as np
import pandas as pd

try:
    import gams
except ImportError:
    sys.exit('The gams library is required to run the GAMS versions of Dispa-SET. Please install if from the /apifiles/Python/api/ folder in the GAMS directory')

def SolveMILP(sim_folder,gams_folder):
    
    if not os.path.exists(gams_folder):
        sys.exit('The provided path for GAMS (' + gams_folder + ') does not exist')    
    
    if not os.path.exists(sim_folder):
        sys.exit('The provided DispaSET simulation environment folder (' + sim_folder + ') does not exist')  
        
    if not os.path.exists(sim_folder + '/Inputs.gdx'):
        sys.exit('There is no Inputs.gdx file within the specified DispaSET simulation environment folder (' + sim_folder + ')')

    if not os.path.exists(sim_folder + '/UCM_h.gms'):
        sys.exit('There is no UCM_h.gms file within the specified DispaSET simulation environment folder (' + sim_folder + ')')
        
    # create GAMS workspace:
    ws = gams.GamsWorkspace(system_directory=gams_folder,debug=3)
    shutil.copy(os.path.join(sim_folder,'UCM_h.gms'),ws.working_directory)
    shutil.copy(os.path.join(sim_folder,'Inputs.gdx'),ws.working_directory)
    t1 = ws.add_job_from_file('UCM_h.gms')    
    t1.run()
    # copy the result file to the simulation environment folder:
    shutil.copy(os.path.join(ws.working_directory,'Results.gdx'),sim_folder)
    return t1
#    for rec in t1.out_db["x"]:
#        print "x(" + rec.key(0) + "," + rec.key(1) + "): level=" + str(rec.level) + " marginal=" + str(rec.marginal)
    
