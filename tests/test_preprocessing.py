import dispaset as ds
import os
import pytest

conf_file = os.path.abspath('./tests/conf.yml')
@pytest.fixture(scope='module',
                params=SIMULATION_TYPES,
                )
def config(request):
    """Generate some data for testing"""
    config = ds.load_config_yaml(conf_file)
    assert isinstance(config, dict)
    config['SimulationType'] = request.param
    return config

