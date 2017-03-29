import logging.config

__version__ = "2.2dev"

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
             "level": "WARNING",
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

# Setting logging configuration:
try: 
    logging.config.dictConfig(_LOGCONFIG)
except:
    # if it didn't work, it might be due to ipython messing with the output
    # typical error: Unable to configure handler 'console': IOStream has no fileno
    # try without console output:
    print('WARNING: the console output is failing (possibly because of ipython). Please refer to the log file')
    _LOGCONFIG = {
         "version": 1,
         "disable_existing_loggers": False,
         'formatters': {
            'standard': {
                'format': '%(asctime)s [%(levelname)-8s] (%(funcName)s): %(message)s',
                'datefmt': '%y/%m/%d %H:%M:%S'
            },
         },
         "handlers": {    
             "error_file": {
                 "class": "logging.FileHandler",
                 "level": "WARNING",
                 'formatter': 'standard',
                 'filename': 'warn.log',
                 'encoding': 'utf8'
    
             }
         },
    
         "root": {
             "level": "INFO",
             "handlers": [ "error_file"],
         }
    }    
    logging.config.dictConfig(_LOGCONFIG)