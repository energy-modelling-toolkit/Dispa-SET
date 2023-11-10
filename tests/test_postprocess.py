import dispaset as ds
import os
import sys
import pytest
import xarray as xr

SIM_DIR = os.path.abspath('./tests/dummy_results')

@pytest.mark.skipif(sys.version_info < (3, 5), reason="requires python3.5 or higher due to incompatible pickle file in tests.")
def test_read_results_dicts():
    # Using temp dir to ensure that each time a new directory is used
    inputs, results = ds.get_sim_results(path=SIM_DIR, return_xarray=False, return_status=False)
    assert isinstance(inputs, dict)
    assert isinstance(results, dict)

# @pytest.mark.skipif(sys.version_info < (3, 5), reason="requires python3.5 or higher due to incompatible pickle file in tests.")
# def test_read_results_xarray():
#     # Using temp dir to ensure that each time a new directory is used
#     inputs, results = ds.get_sim_results(path=SIM_DIR, return_xarray=True, return_status=False)
#     assert isinstance(inputs, xr.Dataset)
#     assert isinstance(results, xr.Dataset)

