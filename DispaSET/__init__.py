import platform
import sys
import logging.config

__version__ = "2.2dev"

# Logging:
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
         "console": {
            "class": "DispaSET.misc.colorstreamhandler.ColorStreamHandler",
             "stream": "ext://sys.stderr",
             "level": "INFO",
             'formatter': 'standard',
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
         "handlers": ["console","error_file"],
     }
}


# Try to import the GAMS api and gdxcc to write gdx files:
try:
    import gams
    import gdxcc

    gams_ok = True
    gdxcc_ok = True
except:
    if sys.platform == 'linux2' and platform.architecture()[0] == '64bit':
        sys.path.append('./Externals/gams_api/linux64/')
    elif sys.platform == 'linux2' and platform.architecture()[0] == '32bit':
        sys.path.append('./Externals/gdxcc/linux32/')
    elif sys.platform == 'win32' and platform.architecture()[0] == '64bit':
        sys.path.append('./Externals/gams_api/win64/')
    elif sys.platform == 'win32' and platform.architecture()[0] == '32bit':
        sys.path.append('./Externals/gdxcc/win32/')
    elif sys.platform == 'darwin':
        sys.path.append('./Externals/gdxcc/osx64/')
    try:
        import gams

        gams_ok = True
    except ImportError:
        gams_ok = False
    try:
        import gdxcc

        gdxcc_ok = True
    except ImportError:
        gdxcc_ok = False

try:
    import pyomo

    pyomo_ok = True
except ImportError:
    pyomo_ok = False

# Importing the main Dispa-SET functions so that they can be called with "ds.function"
from .preprocessing.preprocessing import build_simulation
if gams_ok:
    from .solve import solve_GAMS
if pyomo_ok:
    from .solve import solve_pyomo
from .preprocessing.data_handler import load_config_excel
from .postprocessing.postprocessing import get_sim_results
from .postprocessing.postprocessing import ds_to_df
from .postprocessing.postprocessing import plot_country
from .postprocessing.postprocessing import get_result_analysis
from .postprocessing.postprocessing import get_indicators_powerplant
from .postprocessing.postprocessing import plot_energy_country_fuel


# Setting logging configuration:
logging.config.dictConfig(_LOGCONFIG)