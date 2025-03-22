# -*- coding: utf-8 -*-
"""
Test script to build the simulation environment and debug any errors
"""

import sys
import os
import dispaset as ds

# Add the root folder of Dispa-SET to the path
sys.path.append(os.path.abspath('..'))

# Load the configuration file
config = ds.load_config('ConfigFiles/ConfigTest.yml')
config['GAMS_folder'] = '/home/sylvain/progs/GAMS/gams45.7_linux_x64_64_sfx'  # Update with your GAMS path

# Build the simulation environment
print("Starting to build simulation environment...")
SimData = ds.build_simulation(config)
print("Simulation environment built successfully!") 