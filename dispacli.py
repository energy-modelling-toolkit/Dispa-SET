#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 15:37:16 2016

- add noload cost
- add 'all' header for fuel prices

@author: Sylvain Quoilin, K. Kavvadias
"""
import sys
import logging
import click
import DispaSET as ds


@click.group(chain=True)
@click.option('-c','--config', prompt='Specify the path to the config file', type=click.Path(exists=True),
              default='./ConfigFiles/ConfigTest.xlsx', help='Path to the config file (eg ConfigFiles/Config.xlsx)')
@click.option('-g','--gams', 'engine',
              flag_value='gams', default=True, help='Use GAMS version (default)') #Interface to use for solving. Gams or pyomo')
@click.option('-p','--pyomo', 'engine',
              flag_value='pyomo', help='Use pyomo version')
@click.version_option(prog_name='DispaSET', version=ds.__version__)
@click.pass_context
def cli(ctx, config, engine):
    """Build and run the Dispaset model according to a config file.
     E.g. ./dispacli.py -c ./ConfigFiles/ConfigTest.xlsx build simulate
    """
    if config.endswith(('.xlsx','.xls')):
        ctx.obj['conf'] = ds.load_config_excel(config)
    elif config.endswith(('.yml','.yaml')):
        ctx.obj['conf'] = ds.load_config_yaml(config)
    else:
        logging.error('Unrecognized file format')
        sys.exit(1)

    ctx.obj['engine'] = engine


@cli.command()
@click.pass_context
def build(ctx):
    """Build simulation files"""

    SimData = ds.build_simulation(ctx.obj['conf'] )


@cli.command()
@click.pass_context
def simulate(ctx):
    """Run GAMS or pyomo for simulation"""
    conf = ctx.obj['conf']
    engine = ctx.obj['engine']

    if engine == 'pyomo':
        r = ds.solve_pyomo(conf['SimulationDirectory'])

    if engine == 'gams':
        r = ds.solve_GAMS(conf['SimulationDirectory'], conf['GAMS_folder'])

if __name__ == '__main__':
    cli(obj={},standalone_mode=False)
