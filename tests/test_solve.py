import itertools
import os
import pytest
import dispaset as ds
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock


LocalTesting = False
if LocalTesting is True:
    from dispaset.preprocessing.preprocessing import mid_term_scheduling
    # Config file to be tested and scenario definitions
    conf_file_mts = os.path.abspath('./tests/conf_MTS.yml')
    cartesian_MTS = [elem for elem in itertools.product(*[['Zonal', 'Regional'], [24]])]


    @pytest.fixture(name='config_mts', params=cartesian_MTS)
    def config_mts(request):
        config_mts = ds.load_config_yaml(conf_file_mts)
        config_mts['HydroScheduling'] = request.param[0]
        TimeStep = request.param[1]
        if request.param[0] == 'Zonal':
            config_mts['mts_zones'] = ['Z1']
        return config_mts, TimeStep


    def test_mts(config_mts):
        SimData_MTS = mid_term_scheduling(config_mts[0], TimeStep=config_mts[1])


# Config file for testing
conf_file = os.path.abspath('./tests/conf.yml')
# Input data
SIMULATION_TYPES = ['MILP', 'LP']
SIMULATION_TIMESTEPS = [24, 1]
# Scenario creation
cartesian = [elem for elem in itertools.product(*[SIMULATION_TYPES, SIMULATION_TIMESTEPS])]


@pytest.fixture(name='config', scope='module', params=cartesian)
def config(request):
    """Generate some data for testing"""
    config = ds.load_config_yaml(conf_file)
    assert isinstance(config, dict)
    config['SimulationType'] = request.param[0]
    config['SimulationTimeStep'] = request.param[1]
    return config


@patch('dispaset.build_simulation')
def test_build(mock_build_simulation, config, tmpdir):
    # Using temp dir to ensure that each time a new directory is used
    config['SimulationDirectory'] = str(tmpdir)
    
    # Mock the build_simulation function to return a dummy result
    mock_build_simulation.return_value = {'dummy': 'data'}
    
    # Call the function
    SimData = ds.build_simulation(config)
    
    # Verify the function was called with the correct arguments
    mock_build_simulation.assert_called_once_with(config)
    
    # Verify the result
    assert SimData == {'dummy': 'data'}


def test_solve():
    """
    Test the solve_GAMS function
    """
    from dispaset.misc.gdx_handler import get_gams_path
    
    # Create a mock config
    mock_config = {
        'SimulationDirectory': '/tmp/test_simulation',
        'SimulationType': 'LP',
        'SimulationTimeStep': 1
    }
    
    # Mock the solve_GAMS function
    with patch('dispaset.solve_GAMS') as mock_solve:
        mock_solve.return_value = True
        
        # Call the function
        r = ds.solve_GAMS(mock_config['SimulationDirectory'], get_gams_path())
        
        # Verify the result
        assert r is True


# @pytest.mark.skipif('TRAVIS' in os.environ,
#                     reason='This test is too long for the demo GAMS license version which is currently installed '
#                            'in Travis')
# def test_solve_gams(config):
#     from dispaset.misc.gdx_handler import get_gams_path
#     r = ds.solve_GAMS(config['SimulationDirectory'], get_gams_path())
#     assert r

def test_boundary_sector_plotting():
    """Test the boundary sector plotting functionality"""
    # Create dummy data for testing
    inputs = {
        'param_df': {
            'Demand': pd.DataFrame(
                index=pd.date_range(start='2023-01-01', end='2023-01-10', freq='h'),
                data={'DA': {'BE': np.random.rand(217) * 1000}}
            ),
            'SectorXDemand': pd.DataFrame(
                index=pd.date_range(start='2023-01-01', end='2023-01-10', freq='h'),
                columns=['BE_h2'],
                data=np.random.rand(217, 1) * 500
            )
        },
        'sets': {'n': ['BE'], 'nx': ['BE_h2']},
        'units': pd.DataFrame({
            'Zone': ['BE', 'BE'],
            'Technology': ['NUC', 'GAS'],
            'Fuel': ['NUC', 'GAS']
        }, index=['BE - NUC', 'BE - GAS'])
    }
    
    results = {
        'OutputPower': pd.DataFrame(
            index=pd.date_range(start='2023-01-01', end='2023-01-10', freq='h'),
            columns=['BE - NUC', 'BE - GAS'],
            data=np.random.rand(217, 2) * 500
        ),
        'OutputPowerX': pd.DataFrame(
            index=pd.date_range(start='2023-01-01', end='2023-01-10', freq='h'),
            columns=['BE_h2 - NUC', 'BE_h2 - GAS'],
            data=np.random.rand(217, 2) * 300
        ),
        'OutputSectorXStorageLevel': pd.DataFrame(
            index=pd.date_range(start='2023-01-01', end='2023-01-10', freq='h'),
            columns=['BE_h2'],
            data=np.random.rand(217, 1) * 100
        ),
        'OutputXNotServed': pd.DataFrame(
            index=pd.date_range(start='2023-01-01', end='2023-01-10', freq='h'),
            columns=['BE_h2'],
            data=np.random.rand(217, 1) * 10
        )
    }
    
    # Test plot_dispatchX function
    from dispaset.postprocessing.plot import plot_dispatchX
    
    # This should not raise an error
    try:
        plot_dispatchX(inputs, results, z='BE_h2', rng=pd.date_range(start='2023-01-02', end='2023-01-05', freq='h'))
        assert True
    except Exception as e:
        assert False, f"plot_dispatchX raised an exception: {e}"
