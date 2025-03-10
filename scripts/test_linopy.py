# -*- coding: utf-8 -*-
"""
Test file for the Linopy solver implementation.
This file runs both GAMS and Linopy solvers and compares their results.

@author: Your Name
"""

import sys,os
sys.path.append(os.path.abspath('..'))

import dispaset as ds
import numpy as np
import pandas as pd

def compare_results(gams_results, linopy_results):
    """Compare results from both solvers"""
    print("\nComparing results:")
    
    if not all(k in linopy_results['Outputs'] for k in gams_results['Outputs']):
        print("WARNING: Some outputs missing in Linopy results")
        return
        
    for key in gams_results['Outputs']:
        gams_val = gams_results['Outputs'][key]
        linopy_val = linopy_results['Outputs'][key]
        
        if isinstance(gams_val, (np.ndarray, pd.DataFrame)):
            diff = np.abs(gams_val - linopy_val).mean()
            print(f"{key}: Mean absolute difference = {diff:.6f}")
        else:
            print(f"{key}: GAMS={gams_val}, Linopy={linopy_val}")

def main():
    # Load test configuration
    config = ds.load_config('../ConfigFiles/ConfigTest.xlsx')
    
    # Build simulation environment
    SimData = ds.build_simulation(config)
    
    # Run both solvers
    print("\nSolving with GAMS...")
    config['GAMS_folder'] = '/home/sylvain/progs/GAMS/gams45.7_linux_x64_64_sfx'  # Update path
    gams_results = ds.solve_GAMS(config['SimulationDirectory'], config['GAMS_folder'])
    
    print("\nSolving with Linopy...")
    linopy_results = ds.solve_linopy(SimData)
    
    # Compare results
    compare_results(gams_results, linopy_results)

if __name__ == "__main__":
    main() 