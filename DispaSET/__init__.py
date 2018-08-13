import logging.config
import os

# Sets the __version__ variable
from ._version import __version__

from .preprocessing.preprocessing import get_git_revision_tag

__gitversion__ = get_git_revision_tag()

# Logging: # TODO: Parametrize in dispacli or external config
_LOGCONFIG = {
    "version": 1,
    "disable_existing_loggers": False,
    'formatters': {
        'standard': {
            'format': '%(asctime)s [%(levelname)-8s] (%(funcName)s): %(message)s',
            'datefmt': '%y/%m/%d %H:%M:%S'
        },
        'notime': {
            'format': '[%(levelname)-8s] (%(funcName)s): %(message)s',
            'datefmt': '%y/%m/%d %H:%M:%S'
        },
    },
    "handlers": {
        "console": {
            "class": "DispaSET.misc.colorstreamhandler.ColorStreamHandler",
            "stream": "ext://sys.stderr",
            #             "stream": "sys.stdout",
            "level": "INFO",
            'formatter': 'notime',
        },

        "error_file": {
            "class": "logging.FileHandler",
            "level": "INFO",
            'formatter': 'standard',
            'filename': 'warn.log',
            'encoding': 'utf8'

        }
    },

    "root": {
        "level": "INFO",
        "handlers": ["console", "error_file"],
    }
}

# Importing the main Dispa-SET functions so that they can be called with "ds.function"
from .preprocessing.data_handler import load_config_excel, load_config_yaml
from .preprocessing.preprocessing import build_simulation, adjust_capacity, adjust_storage

from .solve import solve_GAMS, solve_pyomo

from .postprocessing.data_handler import get_sim_results, ds_to_df
from .postprocessing.postprocessing import get_result_analysis, get_indicators_powerplant, aggregate_by_fuel, CostExPost
from .postprocessing.plot import plot_energy_country_fuel, plot_country_capacities, plot_country

from .cli import *

# Removeing log file:
if os.path.isfile('warn.log'):
    try:
        os.remove('warn.log')
    except:
        print ('Could not erase previous log file "warn.log"')

# Setting logging configuration:
try:
    logging.config.dictConfig(_LOGCONFIG)
except:
    # if it didn't work, it might be due to ipython messing with the output
    # typical error: Unable to configure handler 'console': IOStream has no fileno
    # try without console output:
    print('WARNING: the colored console output is failing (possibly because of ipython). Switching to monochromatic output')
    _LOGCONFIG['handlers']['console']['class'] = "logging.StreamHandler"
    logging.config.dictConfig(_LOGCONFIG)