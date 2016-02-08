# -*- coding: utf-8 -*-
"""

Test script for the DispaSolve functions

Loads a pickle data set and solves it.

Created on Mon Feb  8 17:15:31 2016

@author: sylvain
"""

import pickle
from DispaSolve import *

# Testing the above functions with a generic input:

# Load parameter lists and loop:
[sets, parameters] = pickle.load(open('../Simulation/Inputs.p','rb'))

results = DispaSolve(sets,parameters)