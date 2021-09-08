import logging.config
import logging
import os
import importlib

from .common import commons

# Remove old log file:
for filename in (f for f in os.listdir('.') if f.endswith('.dispa.log')):
    try:
        os.remove(filename)
    except OSError:
        print('Could not erase previous log file ' + filename)

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
            "class": "dispaset.misc.colorstreamhandler.ColorStreamHandler",
            "stream": "ext://sys.stderr",
            # "stream": "sys.stdout",
            "level": "INFO",
            'formatter': 'notime',
        },

        "error_file": {
            "class": "logging.FileHandler",
            "level": "INFO",
            'formatter': 'standard',
            'filename': commons['logfile'],
            'encoding': 'utf8'

        }
    },

    "root": {
        "level": "INFO",
        "handlers": ["console", "error_file"],
    }
}

# Setting logging configuration:
try:
    logging.config.dictConfig(_LOGCONFIG)
except Exception:
    # if it didn't work, it might be due to ipython messing with the output
    # typical error: Unable to configure handler 'console': IOStream has no fileno
    # try without console output:
    logging.warning('The colored console output is failing (possibly because of ipython). '
                    'Switching to monochromatic output')
    _LOGCONFIG['handlers']['console']['class'] = "logging.StreamHandler"
    logging.config.dictConfig(_LOGCONFIG)

# Sets the __version__ variable
try:
    from setuptools_scm import get_version
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        version = get_version(version_scheme='post-release',
                              local_scheme=lambda version: version.format_choice("" if version.exact else "+{node}",
                                                                                 "+dirty"),
                              root='..', relative_to=__file__)
except (ImportError, LookupError):
    try:
        from pkg_resources import get_distribution, DistributionNotFound

        version = get_distribution(__package__).version
    except DistributionNotFound:
        logging.warning("Unable to detect version, most probably because you did not install it properly. "
                        "To avoid further errors, please install it by running 'pip install -e .'.")
        version = 'N/A'
__version__ = version

# Importing the main Dispa-SET functions so that they can be called with "ds.function"
from .preprocessing.data_handler import load_config_excel, load_config_yaml, load_config, export_yaml_config
from .preprocessing.preprocessing import build_simulation, mid_term_scheduling
from .preprocessing.utils import adjust_storage, adjust_capacity

from .solve import solve_GAMS

from .postprocessing.data_handler import get_sim_results, ds_to_df

from .postprocessing.postprocessing import get_result_analysis, get_indicators_powerplant, aggregate_by_fuel, \
    CostExPost, get_EFOH, get_units_operation_cost, get_power_flow_tracing

from .postprocessing.plot import plot_energy_zone_fuel, plot_zone_capacities, plot_zone, storage_levels, \
    plot_storage_levels, plot_EFOH, plot_H2_and_demand, plot_tech_cap, \
    H2_demand_satisfaction, plot_ElyserCap_vs_Utilization, plot_power_flow_tracing_matrix, plot_co2

if importlib.util.find_spec("cartopy"):
    from .postprocessing.geoplot import plot_net_flows_map, plot_line_congestion_map
else:
    logging.warning("Cartopy module not found. The geoplot functions will not be available")

from .cli import *
