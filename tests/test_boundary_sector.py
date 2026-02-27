import os
import pytest
import pandas as pd
import numpy as np
import dispaset as ds
from dispaset.postprocessing.plot import plot_dispatchX, plot_zone
from dispaset.postprocessing.postprocessing import get_plot_data, filter_sector, filter_by_zone


def test_boundary_sector_data_handling():
    """Test the handling of boundary sector data in postprocessing functions"""
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
            ),
            'SectorXStorageCapacity': pd.DataFrame(
                index=['BE_h2 - BATS'],
                columns=['StorageCapacity'],
                data=[1000]
            )
        },
        'sets': {'n': ['BE'], 'nx': ['BE_h2'], 'f': ['NUC', 'GAS']},
        'units': pd.DataFrame({
            'Zone': ['BE', 'BE'],
            'Technology': ['NUC', 'GAS'],
            'Fuel': ['NUC', 'GAS'],
            'Unit': ['BE - NUC', 'BE - GAS'],
            'Sector1': ['BE', 'BE'],
            'PowerCapacity': [1000, 800],
            'PartLoadMin': [0.3, 0.4],
            'RampUp': [100, 200],
            'RampDown': [100, 200],
            'StartUpCost': [1000, 500],
            'NoLoadCost': [50, 30],
            'MinUpTime': [4, 2],
            'MinDownTime': [4, 2]
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
        'OutputSectorXStorageInput': pd.DataFrame(
            index=pd.date_range(start='2023-01-01', end='2023-01-10', freq='h'),
            columns=['BE_h2 - BATS'],
            data=np.random.rand(217, 1) * 50 - 25  # Mix of positive and negative values
        ),
        'OutputSectorXFlow': pd.DataFrame(
            index=pd.date_range(start='2023-01-01', end='2023-01-10', freq='h'),
            columns=['BE->BE_h2', 'BE_h2->BE'],
            data=np.random.rand(217, 2) * 100
        ),
        'OutputXNotServed': pd.DataFrame(
            index=pd.date_range(start='2023-01-01', end='2023-01-10', freq='h'),
            columns=['BE_h2'],
            data=np.random.rand(217, 1) * 10
        )
    }
    
    # Test filter_sector function
    sector_data = filter_sector(inputs['param_df']['SectorXDemand'], inputs)
    assert isinstance(sector_data, pd.DataFrame)
    assert not sector_data.empty
    
    # Test filter_by_zone with sector=True
    zone_data = filter_by_zone(results['OutputPowerX'], inputs, 'BE_h2', sector=True)
    assert isinstance(zone_data, pd.DataFrame)
    assert not zone_data.empty
    
    # Test get_plot_data with boundary sector
    plot_data = get_plot_data(inputs, results, 'BE_h2')
    assert isinstance(plot_data, pd.DataFrame)
    assert not plot_data.empty


def test_boundary_sector_plotting_functions():
    """Test the plotting functions for boundary sector data"""
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
            ),
            'SectorXStorageCapacity': pd.DataFrame(
                index=['BE_h2 - BATS'],
                columns=['StorageCapacity'],
                data=[1000]
            )
        },
        'sets': {'n': ['BE'], 'nx': ['BE_h2'], 'f': ['NUC', 'GAS']},
        'units': pd.DataFrame({
            'Zone': ['BE', 'BE'],
            'Technology': ['NUC', 'GAS'],
            'Fuel': ['NUC', 'GAS'],
            'Unit': ['BE - NUC', 'BE - GAS'],
            'Sector1': ['BE', 'BE'],
            'PowerCapacity': [1000, 800],
            'PartLoadMin': [0.3, 0.4],
            'RampUp': [100, 200],
            'RampDown': [100, 200],
            'StartUpCost': [1000, 500],
            'NoLoadCost': [50, 30],
            'MinUpTime': [4, 2],
            'MinDownTime': [4, 2]
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
        'OutputSectorXStorageInput': pd.DataFrame(
            index=pd.date_range(start='2023-01-01', end='2023-01-10', freq='h'),
            columns=['BE_h2 - BATS'],
            data=np.random.rand(217, 1) * 50 - 25
        ),
        'OutputSectorXFlow': pd.DataFrame(
            index=pd.date_range(start='2023-01-01', end='2023-01-10', freq='h'),
            columns=['BE->BE_h2', 'BE_h2->BE'],
            data=np.random.rand(217, 2) * 100
        ),
        'OutputXNotServed': pd.DataFrame(
            index=pd.date_range(start='2023-01-01', end='2023-01-10', freq='h'),
            columns=['BE_h2'],
            data=np.random.rand(217, 1) * 10
        )
    }

    # Test plot_dispatchX with valid time range
    valid_start = pd.Timestamp('2023-01-01')
    valid_end = pd.Timestamp('2023-01-05')
    valid_rng = pd.date_range(start=valid_start, end=valid_end, freq='h') # Create DatetimeIndex
    try:
        fig = plot_dispatchX(inputs, results, 'BE_h2', rng=valid_rng)
        assert fig is not None
    except Exception as e:
        pytest.fail(f"plot_dispatchX failed with valid time range: {str(e)}")

    # Test plot_dispatchX with time range outside available data
    invalid_start = pd.Timestamp('2022-12-01')
    invalid_end = pd.Timestamp('2022-12-31')
    invalid_rng = pd.date_range(start=invalid_start, end=invalid_end, freq='h') # Create DatetimeIndex
    try:
        fig = plot_dispatchX(inputs, results, 'BE_h2', rng=invalid_rng)
        assert fig is not None  # Should still work but use available data range
    except Exception as e:
        pytest.fail(f"plot_dispatchX failed with invalid time range: {str(e)}")

    # Test plot_dispatchX with custom colors
    try:
        fig = plot_dispatchX(inputs, results, 'BE_h2', colors={'NUC': 'red', 'GAS': 'blue'})
        assert fig is not None
    except Exception as e:
        pytest.fail(f"plot_dispatchX failed with custom colors: {str(e)}")

    # Test plot_dispatchX with dispatch and storage limits
    try:
        fig = plot_dispatchX(inputs, results, 'BE_h2')
        assert fig is not None
    except Exception as e:
        pytest.fail(f"plot_dispatchX failed with dispatch and storage limits: {str(e)}") 