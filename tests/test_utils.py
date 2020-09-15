import os
import pytest
import dispaset as ds
import numpy as np

from dispaset.preprocessing.utils import pd_timestep, _mylogspace, _find_nearest, _flatten_list
# TODO: test following functions: incidence_matrix, interconnections, _get_index, group_plants, clustering, adjust_storage, adjust_capacity


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
