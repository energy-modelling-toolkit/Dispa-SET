import dispaset as ds
import os
import pytest

conf_file = os.path.abspath('./tests/conf.yml')

SIMULATION_TYPES = ['MILP', 'LP']


@pytest.fixture(scope='module',
                params=SIMULATION_TYPES,
                )
def config(request):
    """Generate some data for testing"""
    config = ds.load_config_yaml(conf_file)
    assert isinstance(config, dict)
    config['SimulationType'] = request.param
    return config

def test_build(config, tmpdir):
    # Using temp dir to ensure that each time a new directory is used
    config['SimulationDirectory'] = str(tmpdir)
    SimData = ds.build_simulation(config)
    assert isinstance(SimData, dict) #how to test if sucessful build?

@pytest.mark.skipif('TRAVIS' in os.environ,
                    reason='This test is too long for the demo GAMS license version which is currently installed in Travis')
def test_solve_gams(config):
    from dispaset.misc.gdx_handler import get_gams_path
    r = ds.solve_GAMS(config['SimulationDirectory'], get_gams_path())

    assert r

@pytest.mark.skip(reason="Skipping as pyomo version is not currently up to date")
def test_solve_pyomo(config):
    r = ds.solve_pyomo(config['SimulationDirectory'])

    assert r