import dispaset as ds
import os
import pytest

conf_file = os.path.abspath('./tests/conf.yml')


def test_load_conf():
    config = ds.load_config_yaml(conf_file)
    assert isinstance(config, dict)

def test_build():
    config = ds.load_config_yaml(conf_file)

    SimData = ds.build_simulation(config)
    assert isinstance(SimData, dict) #how to test if sucessful build?

def test_solve_gams():
    from dispaset.misc.gdx_handler import get_gams_path
    config = ds.load_config_yaml(conf_file)

    r = ds.solve_GAMS(config['SimulationDirectory'], get_gams_path())

    assert r

@pytest.mark.skip(reason="Skipping as pyomo version is not currently up to date")
def test_solve_pyomo():
    config = ds.load_config_yaml(conf_file)
    r = ds.solve_pyomo(config['SimulationDirectory'])

    assert r