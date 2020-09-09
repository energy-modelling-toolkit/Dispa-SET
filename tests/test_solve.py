import itertools
import os
import pytest
import dispaset as ds

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


def test_build(config, tmpdir):
    # Using temp dir to ensure that each time a new directory is used
    config['SimulationDirectory'] = str(tmpdir)
    SimData = ds.build_simulation(config)
    assert isinstance(SimData, dict)  # how to test if sucessful build?


@pytest.mark.skipif('TRAVIS' in os.environ,
                    reason='This test is too long for the demo GAMS license version which is currently installed '
                           'in Travis')
def test_solve_gams(config):
    from dispaset.misc.gdx_handler import get_gams_path
    r = ds.solve_GAMS(config['SimulationDirectory'], get_gams_path())
    assert r
