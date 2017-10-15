import logging.config
import os

__version__ = "2.2.0"

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
from .preprocessing.preprocessing import build_simulation
from .solve import solve_GAMS, solve_pyomo
from .preprocessing.data_handler import load_config_excel, load_config_yaml
from .postprocessing.postprocessing import get_sim_results
from .postprocessing.postprocessing import ds_to_df
from .postprocessing.postprocessing import plot_country
from .postprocessing.postprocessing import get_result_analysis
from .postprocessing.postprocessing import get_indicators_powerplant
from .postprocessing.postprocessing import plot_energy_country_fuel

# Removeing log file:
if os.path.isfile('warn.log'):
    try:
        os.remove('warn.log')
    except:
        print 'Could not erase previous log file "warn.log"'
        
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