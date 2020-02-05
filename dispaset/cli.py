import sys
import logging
import click

from .preprocessing.data_handler import load_config_excel,load_config_yaml
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
     E.g. ./dispaset.py -c ./ConfigFiles/ConfigTest.xlsx build simulate
    """
    ctx.obj = {}
    if config.endswith(('.xlsx','.xls')):
        ctx.obj['conf'] = load_config_excel(config)
    elif config.endswith(('.yml','.yaml')):
        ctx.obj['conf'] = load_config_yaml(config)
    else:
        logging.error('Unrecognized file format')
        sys.exit(1)


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

#if __name__ == '__main__':
#    cli(obj={},standalone_mode=False)
