"""
Linopy-based solver for the Dispa-SET optimization problem.
This is an alternative to the GAMS-based solver.

@author: Your Name
"""

import logging
import pandas as pd
import numpy as np
import xarray as xr
from linopy import Model
import datetime as dt

def create_parameter_array(param_data, sets_data):
    """
    Create an xarray DataArray from parameter data using set information
    
    :param param_data: Dictionary containing 'sets' and 'val' keys
    :param sets_data: Dictionary containing the actual set values
    :return: xarray DataArray with proper dimensions and coordinates
    """
    # Create coords dictionary from sets
    coords = {set_name: sets_data[set_name] for set_name in param_data['sets']}
    
    # Verify dimensions match
    expected_shape = tuple(len(sets_data[set_name]) for set_name in param_data['sets'])
    if param_data['val'].shape != expected_shape:
        logging.error(f"Parameter shape mismatch: expected {expected_shape}, got {param_data['val'].shape}")
        raise ValueError(f"Parameter dimensions do not match set sizes")
    
    # Create and return DataArray
    return xr.DataArray(
        data=param_data['val'],
        dims=param_data['sets'],
        coords=coords
    )

def solve_linopy(SimData):
    """
    Function that runs the Dispa-SET optimization using linopy
    
    :param SimData:      Dictionary with all simulation data from the preprocessing
    :return:            Dictionary with the optimization results
    """
    logging.info('Building optimization problem using Linopy')
    
    # Extract key data from SimData
    sets = SimData['sets']
    parameters = SimData['parameters']
    config = SimData['config']
    
    # Log configuration details
    logging.info(f"Start date: {config['StartDate']}")
    logging.info(f"End date: {config['StopDate']}")
    logging.info(f"Horizon length: {config['HorizonLength']} days")
    logging.info(f"Look ahead: {config['LookAhead']} days")
    
    # Create time indices
    start_date = dt.datetime(config['StartDate'][0], config['StartDate'][1], 
                           config['StartDate'][2], config['StartDate'][3])
    end_date = dt.datetime(config['StopDate'][0], config['StopDate'][1], 
                          config['StopDate'][2], config['StopDate'][3])
    
    horizon_hours = config['HorizonLength'] * 24
    overlap_hours = config['LookAhead'] * 24
    
    # Create full time index
    full_timeindex = pd.date_range(start=start_date, end=end_date, freq='h')
    logging.info(f"Created time index with {len(full_timeindex)} hours")
    
    # Create base dimensions dictionary
    dims = {
        'h': full_timeindex,          # hours
        'u': list(sets['u']),         # units
        'n': list(sets['n']),         # nodes
        'f': list(sets['f']),         # fuel types
        't': list(sets['t']),         # technologies
        'l': list(sets['l']),         # lines
        'mk': list(sets['mk']),       # markets
        'p': list(sets['p']),         # pollutants
        's': list(sets['s']),         # storage units
    }
    
    # Log dimensions
    for dim_name, dim_values in dims.items():
        logging.info(f"Dimension {dim_name}: {len(dim_values)} elements")
    
    # Create parameter arrays
    param_arrays = {}
    for param_name, param_data in parameters.items():
        logging.info(f"Processing parameter: {param_name}")
        param_arrays[param_name] = create_parameter_array(param_data, sets)
    
    # Initialize results arrays with proper dimensions
    results = {
        'OutputPower': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputCommitted': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputStorageLevel': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['s']))),
            dims=['h', 's'],
            coords={'h': dims['h'], 's': dims['s']}
        ),
        'OutputStorageInput': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['s']))),
            dims=['h', 's'],
            coords={'h': dims['h'], 's': dims['s']}
        ),
        'OutputFlow': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['l']))),
            dims=['h', 'l'],
            coords={'h': dims['h'], 'l': dims['l']}
        ),
        'OutputSystemCost': xr.DataArray(
            np.zeros(len(full_timeindex)),
            dims=['h'],
            coords={'h': dims['h']}
        ),
        'ShadowPrice': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputCurtailedPower': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputShedLoad': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputReserve_2U': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputReserve_2D': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputReserve_3U': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputEmissions': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']), len(sets['p']))),
            dims=['h', 'n', 'p'],
            coords={'h': dims['h'], 'n': dims['n'], 'p': dims['p']}
        ),
        'LostLoad_MaxPower': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'LostLoad_MinPower': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'LostLoad_2U': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'LostLoad_2D': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'LostLoad_3U': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputDemand_2U': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputDemand_3U': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputDemand_2D': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputMaxOutageUp': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputMaxOutageDown': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputPowerX': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']), len(sets['u']))),
            dims=['h', 'n', 'u'],
            coords={'h': dims['h'], 'n': dims['n'], 'u': dims['u']}
        ),
        'OutputPowerConsumption': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputResidualLoad': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputHeat': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputHeatSlack': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputH2Slack': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputXNotServed': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'OutputH2Output': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputPowerMustRun': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputStartUp': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputShutDown': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputRampRate': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'CapacityMargin': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['n']))),
            dims=['h', 'n'],
            coords={'h': dims['h'], 'n': dims['n']}
        ),
        'UnitHourlyPowerRevenue': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'UnitHourlyFixedCost': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'UnitHourlyVariableCost': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'UnitHourlyStartUpCost': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'UnitHourlyShutDownCost': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'UnitHourlyRampingCost': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'UnitHourlyProductionCost': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'UnitHourlyProfit': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputSysInertia': xr.DataArray(
            np.zeros(len(full_timeindex)),
            dims=['h'],
            coords={'h': dims['h']}
        ),
        'OutputSystemGain': xr.DataArray(
            np.zeros(len(full_timeindex)),
            dims=['h'],
            coords={'h': dims['h']}
        ),
        'OutputPrimaryReserve_Available': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputFFRGain': xr.DataArray(
            np.zeros(len(full_timeindex)),
            dims=['h'],
            coords={'h': dims['h']}
        ),
        'OutputFFR_Available': xr.DataArray(
            np.zeros((len(full_timeindex), len(sets['u']))),
            dims=['h', 'u'],
            coords={'h': dims['h'], 'u': dims['u']}
        ),
        'OutputPowerLoss': xr.DataArray(
            np.zeros(len(full_timeindex)),
            dims=['h'],
            coords={'h': dims['h']}
        ),
        'OutputTotalDemand_2U': xr.DataArray(
            np.zeros(len(full_timeindex)),
            dims=['h'],
            coords={'h': dims['h']}
        )
    }
    
    # Rolling horizon loop
    for start_idx in range(0, len(full_timeindex), horizon_hours - overlap_hours):
        # Get time slice for this optimization
        end_idx = min(start_idx + horizon_hours, len(full_timeindex))
        time_slice = slice(start_idx, end_idx)
        current_times = full_timeindex[time_slice]
        
        logging.info(f'Optimizing period from {current_times[0]} to {current_times[-1]}')
        
        # Create Linopy model for this horizon
        model = Model()
        
        # TODO: Add variables, constraints and objective function
        # This will be implemented in subsequent iterations
        
        # For now, just fill with dummy values (ones)
        results['OutputPower'].loc[dict(h=current_times)] = 1.0
        results['OutputCommitted'].loc[dict(h=current_times)] = 1.0
        results['OutputStorageLevel'].loc[dict(h=current_times)] = 1.0
        results['OutputSystemCost'].loc[dict(h=current_times)] = 1.0
        results['ShadowPrice'].loc[dict(h=current_times)] = 1.0
        
        # Only keep non-overlap period results
        keep_until = min(start_idx + horizon_hours - overlap_hours, len(full_timeindex))
    
    logging.info('Optimization completed successfully')
    

    
    return results
