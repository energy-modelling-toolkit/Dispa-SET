#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 15:37:16 2016

- add noload cost
- add 'all' header for fuel prices

@author: Sylvain Quoilin, K. Kavvadias
"""

import click
import DispaSET as ds
import logging


@click.group(chain=True)
@click.option('-c','--config', prompt='Specify the path to the config file', type=click.Path(exists=True),
              default='./ConfigFiles/ConfigTest.xlsx', help='Path to the config file (eg ConfigFiles/Config.xlsx)')
@click.option('-g','--gams', 'engine',
              flag_value='gams', default=True, help='Use GAMS version (default)') #Interface to use for solving. Gams or pyomo')
@click.option('-p','--pyomo', 'engine',
              flag_value='pyomo', help='Use pyomo version')
#@click.option('--log', help='Where to print log, screen,file?')
@click.version_option(prog_name='DispaSET', version=ds.__version__)
@click.pass_context
def cli(ctx, config, engine):
    """Build and run the Dispaset model according to a config file.
     E.g. ./dispacli.py -c ./ConfigFiles/ConfigTest.xlsx build simulate
    """
    ctx.obj['conf'] = config
    ctx.obj['engine'] = engine


@cli.command()
@click.option('-l','--loog', help='Use GAMS version (default)') #Interface to use for solving. Gams or pyomo')
@click.pass_context
def build(ctx,loog):
    """Build simulation files"""
    conf = ds.load_config_excel(ctx.obj['conf'])
    logging.info("Using config file " + ctx.obj['conf'] + " to build the simulation environment")

    if not ds.gdxcc_ok:
        click.echo('WARNING: the gdxcc library could not be loaded. The gdx will not be generated')

    logging.info("New build started")

    SimData = ds.build_simulation(conf)

    logging.info("Build finished")


@cli.command()
@click.pass_context
def simulate(ctx):
    """Run GAMS or pyomo for simulation"""
    conf = ds.load_config_excel(ctx.obj['conf'])
    engine = ctx.obj['engine']

    if engine == 'pyomo':
        if ds.pyomo_ok:
            r = ds.solve_pyomo(conf['SimulationDirectory'])
        else:
            logging.critical('Pyomo and a compatible solver is needed to run the pyomo version of DispaSET. Please install by typing "pip install pyomo"')

    if engine == 'gams':
        if ds.gams_ok:
            r = ds.solve_GAMS(conf['SimulationDirectory'], conf['GAMS_folder'])
        else:
            logging.critical('The gams library is required to run the GAMS versions of DispaSET. Please install if from the /apifiles/Python/api/ folder in the GAMS directory')

if __name__ == '__main__':
    cli(obj={})
