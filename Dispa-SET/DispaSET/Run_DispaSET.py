# -*- coding: utf-8 -*-
"""

Test script for the DispaSolve functions

Loads a pickle data set and solves it.

Created on Mon Feb  8 17:15:31 2016

@author: sylvain
"""

import pickle
import logging
from DispaSolve import *


if __name__ == "__main__":
    # Testing the above functions with a generic input:
    logging.basicConfig(format='%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.DEBUG) # filename='run.log',
    logging.info("New run started")
    # Load parameter lists and loop:
    [sets, parameters] = pickle.load(open('../Simulation/Inputs.p','rb'))
    results = DispaSolve(sets, parameters, Mixed_Integer_LP=False)
    logging.info("Run finished")
