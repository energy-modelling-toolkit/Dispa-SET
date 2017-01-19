# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 15:37:16 2016

- add noload cost
- add 'all' header for fuel prices

@author: Sylvain Quoilin
"""

import DispaSET as ds
import argparse
import os
import sys

parser = argparse.ArgumentParser(description='Build and run the Dispaset model according to a config file')
parser.add_argument("-b", "--build", help="Path to the config file (eg ConfigFiles/Config.xlsx)")
parser.add_argument("-br", "--build_and_run", help="Path to the config file (eg ConfigFiles/Config.xlsx)")
args = parser.parse_args()

if not ds.gdxcc_ok:
    print 'WARNING: the gdxcc library could not be loaded. The gdx will not be generated'

if args.build and not args.build_and_run:
    print "Using config file " + args.build + " to build the simulation environment"
    path = args.build
    build = True
    simulate = False
elif not args.build and args.build_and_run:
    print "Using config file " + args.build_and_run + " to build the simulation environment and run the simulation"
    path = args.build_and_run
    build = True
    simulate = True
elif args.build and args.build_and_run:
    sys.exit('The "build" and "build_and_run" arguments cannot be defined simultaneously')
else:
    print "No config file specified to build the simulation environment"
    try:
        import easygui
        path = easygui.fileopenbox('Please select the config file')
        print "Using config file " + path + " to build the simulation environment"
    except:
        path = ds.InputFile('Specify the path to the config file (e.g. ../Config.xlsx): ')
        print "Using config file " + path + " to build the simulation environment"
    build = True
    simulate = False

if not os.path.isfile(path):
    sys.exit('Could not find config file ' + path + ' \n Exiting')


config = ds.LoadConfig(path)

if not os.path.isdir(config['GAMS_folder']):
    config['GAMS_folder'] = ds.get_gams_path()

if build: 
    SimData = ds.BuildSimulation(config)

if simulate:
    if ds.gams_ok:
        from DispaSET.DispaSolveGAMS import SolveMILP
        r = SolveMILP(config['SimulationDirectory'],config['GAMS_folder'])
    else:
        print('The gams library is required to run the GAMS versions of Dispa-SET. Please install if from the /apifiles/Python/api/ folder in the GAMS directory')

