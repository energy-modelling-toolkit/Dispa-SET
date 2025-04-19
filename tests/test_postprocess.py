import dispaset as ds
import os
import sys
import pytest
import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from dispaset.postprocessing.postprocessing import get_plot_data

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

@pytest.mark.skipif(sys.version_info < (3, 5), reason="requires python3.5 or higher due to incompatible pickle file in tests.")
def test_time_range_validation():
    """Test the time range validation improvements in plotting functions"""
    # Create dummy data for testing
    start_date = datetime(2023, 1, 1)
    end_date = datetime(2023, 1, 10)
    dates = pd.date_range(start=start_date, end=end_date, freq='h')
    
    # Create dummy inputs and results
    inputs = {
        'param_df': {
            'Demand': pd.DataFrame(
                index=dates,
                data={'DA': {'BE': np.random.rand(len(dates)) * 1000}}
            )
        },
        'sets': {'n': ['BE'], 'f': ['NUC', 'GAS']},
        'units': pd.DataFrame({
            'Zone': ['BE', 'BE'],
            'Technology': ['NUC', 'GAS'],
            'Fuel': ['NUC', 'GAS']
        }, index=['BE - NUC', 'BE - GAS'])
    }
    
    results = {
        'OutputPower': pd.DataFrame(
            index=dates,
            columns=['BE - NUC', 'BE - GAS'],
            data=np.random.rand(len(dates), 2) * 500
        )
    }
    
    # Test with valid time range
    valid_start = datetime(2023, 1, 2)
    valid_end = datetime(2023, 1, 5)
    plot_data = get_plot_data(inputs, results, 'BE')
    assert isinstance(plot_data, pd.DataFrame)
    assert not plot_data.empty
    assert 'NUC' in plot_data.columns
    assert 'Demand' not in plot_data
    
    # Test with missing time range (should use full range)
    plot_data = get_plot_data(inputs, results, 'BE')
    assert isinstance(plot_data, pd.DataFrame)
    assert not plot_data.empty
    assert plot_data.index[0] == start_date
    assert plot_data.index[-1] == end_date

