# -*- coding: utf-8 -*-
"""

Test script for the DispaSolve functions

Loads a pickle data set and solves it.

Created on Mon Feb  8 17:15:31 2016

@author: sylvain
"""

import pickle
import logging
from DispaSolve import DispaSolve


if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.DEBUG) # filename='run.log',
    logging.info("New run started")
    # Load parameter lists and loop:
    SimData = pickle.load(open('../Simulations/simulationNL/Inputs.p','rb'))
    results = DispaSolve(SimData['sets'], SimData['parameters'], LPFormulation=False)
    logging.info("Run finished")
