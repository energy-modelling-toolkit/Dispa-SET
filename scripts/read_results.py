# -*- coding: utf-8 -*-
"""
Minimalist example file showing how to read the results of a Dispa-SET run

@author: Sylvain Quoilin
"""

# Change directory to the root folder of Dispa-SET:
import os
os.chdir('..')

# Import Dispa-SET
import DispaSET as ds

# Load the inputs and the results of the simulation
inputs,results = ds.get_sim_results(path='Simulations/simulation_test',cache=True)

#Format the inputs as a dictionary of dataframes:
datain = ds.ds_to_df(inputs)

# Generate plots
ds.plot_country(inputs,results)

# Analyse the results for each country and provide quantitative indicators:
r = ds.get_result_analysis(inputs,results)
