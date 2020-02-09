import sys
import logging
import click

from .preprocessing.data_handler import load_config
from .preprocessing.preprocessing import build_simulation
from .solve import solve_GAMS
from ._version import __version__

@click.group(chain=True)
@click.option('-c','--config', prompt='Specify the path to the config file', type=click.Path(exists=True),
              default='./ConfigFiles/ConfigTest.xlsx', help='Path to the config file (eg ConfigFiles/Config.xlsx)')
@click.version_option(prog_name='DispaSET', version=__version__)
@click.pass_context
def cli(ctx, config):
    """Build and run the Dispaset model according to a config file.
     E.g. dispaset -c ./ConfigFiles/ConfigTest.xlsx build simulate
    """
    ctx.obj = {}
    ctx.obj['conf'] = load_config(config)

@cli.command()
@click.pass_context
def build(ctx):
    """Build simulation files"""
    conf = ctx.obj['conf']
    __ = build_simulation(ctx.obj['conf'] )


@cli.command()
@click.pass_context
def simulate(ctx):
    """Run GAMS for simulation"""
    conf = ctx.obj['conf']

    r = solve_GAMS(conf['SimulationDirectory'], conf['GAMS_folder'])
