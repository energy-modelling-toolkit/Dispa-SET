import itertools
import os
import pytest
import pandas as pd
import dispaset as ds

from dispaset.preprocessing.data_handler import NodeBasedTable, UnitBasedTable
from dispaset.preprocessing.utils import select_units

# Input data
CONFIG_FOLDER = ['./tests/']
CONFIG_FILE = ['conf.yml', 'conf.xlsx']
# Scenario creation
cartesian = [elem for elem in itertools.product(*[CONFIG_FOLDER, CONFIG_FILE])]


@pytest.fixture(name='config_file', scope='module', params=cartesian)
def config(request):
    # Set destinations for config files
    conf_file = os.path.abspath(request.param[0] + request.param[1])
    return conf_file


def test_config(config_file):
    # Test if config files can be loaded
    config = ds.load_config(config_file)
    assert isinstance(config, dict)


@pytest.fixture(name='config_wrong_title')
def config_title():
    # Set destinations for config files
    conf_file = os.path.abspath('./tests/config_wrong_title.xlsx')
    return conf_file


def test_config_title(config_wrong_title):
    # Test a dummy file and check if sys.exit is being called
    with pytest.raises(SystemExit):
        config = ds.load_config(config_wrong_title)


@pytest.fixture(name='config_wrong_extension')
def config_wrong_extension():
    # Config file with wrong extension
    conf_file = os.path.abspath('./tests/conf.ymx')
    return conf_file


def test_config_wrong_extension(config_wrong_extension):
    # Test a dummy file and check if sys.exit is being called
    with pytest.raises(SystemExit):
        config = ds.load_config(config_wrong_extension)


@pytest.fixture(name='config')
def config_ok():
    config_file = os.path.abspath('./tests/conf.yml')
    config = ds.load_config(config_file)
    config['idx'] = pd.date_range(start='2015-01-01 00:00:00', end='2015-01-01 23:00:00', freq='h')
    idx_long = pd.date_range(start='2015-01-01 00:00:00', end='2015-01-01 23:00:00', freq='h')
    Nhours_long = len(idx_long)
    config['idx_long'] = idx_long
    return config


def test_node_no_data_file(config):
    config['zones'] = ['Z1', 'Z2', 'Z3']
    with pytest.raises(SystemExit):
        NodeBasedTable('Demand',config)


def test_node_single_file(config):
    config['Demand'] = os.path.abspath('./tests/dummy_data/Load_RealTime/2015.csv')
    NodeBasedTable('Demand', config)


def test_node_no_file(config):
    config['Demand'] = os.path.abspath('./tests/dummy_data/Load_RealTime')
    with pytest.raises(SystemExit):
        NodeBasedTable('Demand',config)


@pytest.fixture(name='plants')
def plants(config):
    plants = pd.read_csv(config['PowerPlantData'])
    plants = select_units(plants, config)
    return plants


def test_unit_no_file(config, plants):
    config['Outages'] = os.path.abspath('./tests/dummy_data/buggy_data/')
    with pytest.raises(SystemExit):
        UnitBasedTable(plants, 'Outages', config, fallbacks=['Unit', 'Technology'])


def test_unit_default(config, plants):
    config['Outages'] = os.path.abspath('./tests/dummy_data/buggy_data/')
    with pytest.raises(SystemExit):
        UnitBasedTable(plants, 'Outages', config, fallbacks=['Unit', 'Technology'], default='xxx')
