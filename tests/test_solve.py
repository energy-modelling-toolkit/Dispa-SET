import dispaset as ds


def test_load_conf():
    test_conf = './tests/conf.yml'
    #data = ds.(data_filename)
    data = ds.load_config_yaml('')
    assert True