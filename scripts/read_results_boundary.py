# -*- coding: utf-8 -*-
"""
Minimalist example file showing how to read the results of a Dispa-SET run and plot boundary sector dispatch

@author: Sylvain Quoilin
"""

# Add the root folder of Dispa-SET to the path so that the library can be loaded:
import sys,os
sys.path.append(os.path.abspath('..'))

# Import Dispa-SET
import dispaset as ds
import pandas as pd
import json

# Load the inputs and the results of the simulation
inputs,results = ds.get_sim_results(path='Simulations/simulation_test',cache=False)

# Get the actual date range from the results
if 'OutputPower' in results:
    rng = results['OutputPower'].index
else:
    # Fallback to a default range if PowerDemand is not available
    rng = pd.date_range(start='2016-01-01',end='2016-12-31',freq='h')


# Print information about the PowerX variables
if 'OutputPowerX' in results:
    print("\nPowerX variables:")
    print(results['OutputPowerX'].columns)
    print("Sample data:")
    print(results['OutputPowerX'].head())

# Test the new boundary sector dispatch plot
ds.plot_dispatchX(inputs, results, rng=rng, z='Z1_h2')
