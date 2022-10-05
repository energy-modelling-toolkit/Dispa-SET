import itertools
import os
import pytest
import pandas as pd
import dispaset as ds

from dispaset.preprocessing.utils import select_units
from dispaset.preprocessing.data_check import check_units, check_sto, check_AvailabilityFactors, check_heat_demand, \
                                              check_temperatures, check_clustering, isStorage, check_chp, check_p2h, \
                                              check_h2, check_df, check_MinMaxFlows, check_FlexibleDemand, \
                                              check_reserves, check_PtLDemand
from dispaset.common import commons

TECH = ['NUC', 1.5]
CONFIG_FILE = ['./tests/dummy_data/buggy_data/Units_bugNonNaNKeys1.csv',
               './tests/dummy_data/buggy_data/Units_bugNonNaNKeys2.csv',
               './tests/dummy_data/buggy_data/Units_bugKeys.csv']
cartesian_1 = [elem for elem in itertools.product(*[CONFIG_FILE, TECH])]


@pytest.fixture(name='config')
def config_ok():
    config_file = os.path.abspath('./tests/conf.yml')
    config = ds.load_config(config_file)
    config['idx'] = pd.date_range(start='2015-01-01 00:00:00', end='2015-01-01 23:00:00', freq='h')
    idx_long = pd.date_range(start='2015-01-01 00:00:00', end='2015-01-01 23:00:00', freq='h')
    Nhours_long = len(idx_long)
    config['idx_long'] = idx_long
    return config


@pytest.fixture(name='keys', params=cartesian_1)
def plants(config, request):
    config['PowerPlantData'] = os.path.abspath(request.param[0])
    plants = pd.read_csv(config['PowerPlantData'])
    plants.set_index('Unit', drop=False, inplace=True)
    plants.iloc[7, 5] = request.param[1]
    # plants = select_units(plants, config)
    plants_sto = plants[[u in commons['tech_storage'] for u in plants['Technology']]]
    plants_chp = plants[[str(x).lower() in commons['types_CHP'] for x in plants['CHPType']]]
    plants_p2h = plants[plants['Technology'] == 'P2HT']
    plants_h2 = plants[plants['Technology'] == 'P2GS']
    return config, plants, plants_sto, plants_chp, plants_p2h, plants_h2


def test_bad_keys(keys):
    with pytest.raises(SystemExit):
        check_units(keys[0], keys[1])
        check_sto(keys[0], keys[2])
        check_chp(keys[0], keys[3])
        check_p2h(keys[0], keys[4])
        check_h2(keys[0], keys[5])
