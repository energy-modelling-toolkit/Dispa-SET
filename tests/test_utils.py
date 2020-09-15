import itertools
import os
import pytest
import dispaset as ds
import numpy as np
import pandas as pd

from dispaset.preprocessing.utils import pd_timestep, _mylogspace, _find_nearest, _flatten_list, _split_list, \
    _merge_two_dicts, create_agg_dict

# TODO: test following functions: incidence_matrix, interconnections, group_plants, _get_index, clustering, adjust_storage, adjust_capacity


@pytest.mark.parametrize('TIMESTEP', [1, 0.25, 24, '24'])
def test_pd_timestep(TIMESTEP):
    if (type(TIMESTEP) == int) or (type(TIMESTEP) == float):
        pd_timestep(TIMESTEP)
    else:
        with pytest.raises(SystemExit):
            pd_timestep(TIMESTEP)


@pytest.mark.parametrize('min,max,N,expected', [(0, 1, 3, np.array([0., 0.41421356237309515, 1.]))])
def test_mylogspace(min, max, N, expected):
    results = _mylogspace(min, max, N) == expected
    assert np.sum(results) == 3


@pytest.mark.parametrize('array, value, expected', [(np.array([0, 1, 2, 3, 4]), 3.4, 3)])
def test_find_nearest(array, value, expected):
    result = _find_nearest(array, value) == expected
    assert result


@pytest.mark.parametrize('list, expected', [([1, 3, ['aa', 'bb'], 4], [1, 3, 'aa', 'bb', 4])])
def test_flatten_list(list, expected):
    result = _flatten_list(list) == expected
    assert result


@pytest.mark.parametrize('list, expected', [(['Split', 'me'], 'Split - me')])
def test_split_list(list, expected):
    result = _split_list(list) == expected
    assert result


@pytest.mark.parametrize('x,y,expected',
                         [({'a': 1, 'b': 'c'}, {'d': 2, 'e': 'f'}, {'a': 1, 'b': 'c', 'd': 2, 'e': 'f'})])
def test_merge_two_dict(x, y, expected):
    result = _merge_two_dicts(x, y) == expected
    assert result


CLUSTERING = ['LP clustered', 'MILP', 'Standard', 'Integer clustering', 'No Clustering']
PLANTS = [os.path.abspath('./tests/dummy_data/buggy_data/Clustering.csv')]
cartesian = [elem for elem in itertools.product(*[CLUSTERING, PLANTS])]


@pytest.fixture(scope='module', params=cartesian, name='clustering')
def read_plants(request):
    method = request.param[0]
    plants = pd.read_csv(request.param[1])
    plants.set_index('Unit', inplace=True, drop=False)
    return method, plants


def test_agg_dict(clustering):
    if clustering[0] == 'No Clustering':
        with pytest.raises(SystemExit):
            create_agg_dict(clustering[1], clustering[0])
    else:
        create_agg_dict(clustering[1], clustering[0])
