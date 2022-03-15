import datetime as dt
import logging
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd

from ..misc.gdx_handler import get_gams_path, gdx_to_dataframe, gdx_to_list
from ..misc.str_handler import clean_strings

col_keys = {'OutputCommitted': ('u', 'h'),
            'OutputFlow': ('l', 'h'),
            'OutputPower': ('u', 'h'),
            'OutputPowerConsumption': ('u', 'h'),
            'OutputHeat': ('u', 'h'),
            'OutputHeatSlack': ('n_th', 'h'),
            'OutputStorageInput': ('u', 'h'),
            'OutputStorageLevel': ('u', 'h'),
            'OutputSystemCost': ('h'),
            'OutputSpillage': ('u', 'h'),
            'OutputShedLoad': ('n', 'h'),
            'OutputCurtailedPower': ('n', 'h'),
            'OutputCurtailedHeat': ('n_th','h'),
            'OutputDemandModulation': ('n', 'h'),
            'LostLoad_MaxPower': ('n', 'h'),
            'LostLoad_MinPower': ('n', 'h'),
            'LostLoad_2D': ('n', 'h'),
            'LostLoad_2U': ('n', 'h'),
            'LostLoad_3U': ('n', 'h'),
            'LostLoad_RampUp': ('n', 'h'),
            'LostLoad_RampDown': ('n', 'h'),
            'ShadowPrice': ('n', 'h'),
            'StorageShadowPrice': ('u', 'h'),
            'LostLoad_WaterSlack': ('u'),
            'OutputH2Output': ('u', 'h'),
            'OutputStorageSlack': ('u', 'h'),
            'OutputPtLDemand': ('u', 'h'),
            'HeatShadowPrice': ('n_th', 'h'),
            'H2ShadowPrice': ('u', 'h'),
            'ShadowPrice_2U': ('u', 'h'),
            'ShadowPrice_2D': ('u', 'h'),
            'ShadowPrice_3U': ('u', 'h'),
            'OutputReserve_2U': ('u', 'h'),
            'OutputReserve_2D': ('u', 'h'),
            'OutputReserve_3U': ('u', 'h'),
            'ShadowPrice_RampUp_TC': ('u', 'h'),
            'ShadowPrice_RampDown_TC': ('u', 'h'),
            'OutputRampRate': ('u', 'h'),
            'OutputStartUp': ('u', 'h'),
            'OutputShutDown': ('u', 'h'),
            'OutputPowerMustRun': ('u', 'h'),
            'OutputCostStartUpH': ('u', 'h'),
            'OutputCostRampUpH': ('u', 'h'),
            'OutputCostRampDownH': ('u', 'h'),
            'status': tuple(),
            '*': tuple()
            }


def GAMSstatus(statustype, num):
    """
    Function that returns the model status or the solve status from gams

    :param statustype: String with the type of status to retrieve ("solver" or "model")
    :param num:     Indicated termination condition (Integer)
    :returns:       String with the status
    """
    if statustype == "model":
        msg = {1: u'Optimal solution achieved',
               2: u'Local optimal solution achieved',
               3: u'Unbounded model found',
               4: u'Infeasible model found',
               5: u'Locally infeasible model found (in NLPs)',
               6: u'Solver terminated early and model was infeasible',
               7: u'Solver terminated early and model was feasible but not yet optimal',
               8: u'Integer solution model found',
               9: u'Solver terminated early with a non integer solution found (only in MIPs)',
               10: u'No feasible integer solution could be found',
               11: u'Licensing problem',
               12: u'Error achieved \u2013 No cause known',
               13: u'Error achieved \u2013 No solution attained',
               14: u'No solution returned',
               15: u'Feasible in a CNS models',
               16: u'Locally feasible in a CNS models',
               17: u'Singular in a CNS models',
               18: u'Unbounded \u2013 no solution',
               19: u'Infeasible \u2013 no solution'}
    elif statustype == "solver":
        msg = {1: u'Normal termination',
               2: u'Solver ran out of iterations (fix with iterlim)',
               3: u'Solver exceeded time limit (fix with reslim)',
               4: u'Solver quit with a problem (see LST file) found',
               5: u'Solver quit with excessive nonlinear term evaluation errors (see LST file and fix with bounds or domlim)',
               6: u'Solver terminated for unknown reason (see LST file)',
               7: u'Solver terminated with preprocessor error (see LST file)',
               8: u'User interrupt',
               9: u'Solver terminated with some type of failure (see LST file)',
               10: u'Solver terminated with some type of failure (see LST file)',
               11: u'Solver terminated with some type of failure (see LST file)',
               12: u'Solver terminated with some type of failure (see LST file)',
               13: u'Solver terminated with some type of failure (see LST file)'}
    else:
        sys.exit('Incorrect GAMS status type')
    return str(msg[num])


def get_sim_results(path='.', cache=None, temp_path=None, return_xarray=False, return_status=False):
    """
    This function reads the simulation environment folder once it has been solved and loads
    the input variables together with the results.

    :param path:                Relative path to the simulation environment folder (current path by default)
    :param return_xarray:       If true the results are returned as a multidimensional xarray.
                                Otherwise a dict of dataframes will be returned.
    :param return_status:       IF true the output of this function is a tuple containng the following:
                                (inputs, results, status). The latter is a dictionary with diagnostic messages.
    :returns inputs,results:    Two dictionaries with all the input and outputs
    """

    inputfile = path + '/Inputs.p'
    resultfile = path + '/Results.gdx'
    if cache is not None or temp_path is not None:
        logging.warning(
            'Caching option has been removed. Try to save manually the results, e.g. results.to_netcdf("res.nc")')

    inputs = pd.read_pickle(inputfile)

    # Clean power plant names:
    inputs['sets']['u'] = clean_strings(inputs['sets']['u'])
    inputs['units'].index = clean_strings(inputs['units'].index.tolist())

    # Add the formated parameters in the inputs variable if not already present:
    if not 'param_df' in inputs:
        inputs['param_df'] = ds_to_df(inputs)

    # We need to pass the dir in config if we run it in clusters. PBS script fail to autolocate
    gams_dir = get_gams_path(gams_dir=inputs['config']['GAMS_folder'].encode())
    if not gams_dir:  # couldn't locate
        logging.error('GAMS path cannot be located. Cannot parse gdx files')
        return False

    results = gdx_to_dataframe(gdx_to_list(gams_dir, resultfile, varname='all', verbose=True),
                               fixindex=True, verbose=True)

    # Set datetime index:
    StartDate = inputs['config']['StartDate']
    StopDate = inputs['config']['StopDate']  # last day of the simulation with look-ahead period
    StopDate_long = dt.datetime(*StopDate) + dt.timedelta(days=inputs['config']['LookAhead'])
    index = pd.date_range(start=dt.datetime(*StartDate), end=dt.datetime(*StopDate), freq='h')
    index_long = pd.date_range(start=dt.datetime(*StartDate), end=StopDate_long, freq='h')

    keys = ['LostLoad_2U', 'LostLoad_3U', 'LostLoad_MaxPower', 'LostLoad_MinPower', 'LostLoad_RampUp',
            'LostLoad_RampDown', 'LostLoad_2D', 'ShadowPrice', 'StorageShadowPrice',
            'OutputCostStartUpH', 'OutputCostRampUpH', 'ShadowPrice_2U', 'ShadowPrice_2D', 'ShadowPrice_3U', 'status']  # 'status'

    keys_sparse = ['OutputPower', 'OutputPowerConsumption', 'OutputSystemCost', 'OutputCommitted',
                   'OutputCurtailedPower', 'OutputFlow', 'OutputShedLoad', 'OutputSpillage', 'OutputStorageLevel',
                   'OutputStorageInput', 'OutputHeat', 'OutputHeatSlack', 'OutputDemandModulation',
                   'OutputStorageSlack', 'OutputPtLDemand', 'OutputH2Output', 'OutputPowerMustRun',
                   'OutputReserve_2U', 'OutputReserve_2D', 'OutputReserve_3U', 'ShadowPrice_RampUp_TC',
                   'ShadowPrice_RampDown_TC', 'OutputRampRate', 'OutputStartUp', 'OutputShutDown', 'HeatShadowPrice',
                   'H2ShadowPrice', 'OutputCurtailedHeat', 'OutputEmissions', 'CapacityMargin']

    # Setting the proper index to the result dataframes:
    from itertools import chain
    for key in chain(keys, keys_sparse):
        if key in results:
            # Case of variables for which the look-ahead period is recorded (e.g. the lost loads)
            if len(results[key]) == len(index_long):
                results[key].index = index_long
            # Case of variables for which the look-ahead is not recorded (standard case)
            elif len(results[key]) == len(index):
                results[key].index = index
            # Variables whose index is not complete (sparse formulation)
            else:
                results[key].index = index_long[results[key].index - 1]
                if key in keys_sparse:
                    results[key] = results[key].reindex(index).fillna(0)
        else:
            results[key] = pd.DataFrame(index=index)
        if isinstance(results[key], type(pd.DataFrame())):
            results[key] = results[key].reindex(sorted(results[key].columns), axis=1)

    # Include water slack in the results (only one number)
    if 'LostLoad_WaterSlack' in results:
        results['LostLoad_WaterSlack'] = results['LostLoad_WaterSlack']
    else:
        results['LostLoad_WaterSlack'] = 0

    # Clean power plant names:
    results['OutputPower'].columns = clean_strings(results['OutputPower'].columns.tolist())
    # Remove epsilons:
    if 'ShadowPrice' in results:
        results['ShadowPrice'][results['ShadowPrice'] >= 1e300] = 0
        results['StorageShadowPrice'][results['StorageShadowPrice'] >= 1e300] = 0
        results['ShadowPrice_2D'][results['ShadowPrice_2D'] >= 1e300] = 0
        results['ShadowPrice_2U'][results['ShadowPrice_2U'] >= 1e300] = 0
        results['ShadowPrice_3U'][results['ShadowPrice_3U'] >= 1e300] = 0
        results['HeatShadowPrice'][results['HeatShadowPrice'] >= 1e300] = 0
        results['H2ShadowPrice'][results['H2ShadowPrice'] >= 1e300] = 0

    status = {}
    if "model" in results['status']:
        errors = results['status'][(results['status']['model'] != 1) & (results['status']['model'] != 8)]
        if len(errors) > 0:
            logging.critical('Some simulation errors were encountered. Some results could not be computed, for example at \n \
                            time ' + str(errors.index[0]) + ', with the error message: "' + GAMSstatus('model', errors[
                'model'].iloc[0]) + '". \n \
                            The complete list is available in results["errors"] \n \
                            The optimization might be debugged by activating the Debug flag in the GAMS simulation file and running it')
            for i in errors.index:
                errors.loc[i, 'Error Message'] = GAMSstatus('model', errors['model'][i])
            status['errors'] = errors

    if return_xarray:
        if return_status:
            status['*'] = results.pop('*')
            status['status'] = results.pop('status')
        results = results_to_xarray(results)
        inputs = inputs_to_xarray(inputs)
        out = (inputs, results)
        if return_status:
            return out + (status)
        else:
            return out
    else:
        return (inputs, results)


def results_to_xarray(results):
    """
    Convert results to xarray format

    :param results:     Dispa-SET results
    :return:            results in xarray format
    """

    try:
        import xarray as xr
        all_ds = []

        for k, v in results.items():
            df = v
            var_name = k
            ind = col_keys[k]
            if len(ind) == 2:
                ds = xr.DataArray(df.values,
                                  coords={ind[1]: df.index.values,
                                          ind[0]: df.columns.values},
                                  dims=[ind[1], ind[0]],
                                  name=var_name)
                all_ds.append(ds)
            elif len(ind) == 1:
                if k == 'LostLoad_WaterSlack':
                    ds = xr.DataArray(df,
                                      name=var_name)
                else:
                    ds = xr.DataArray(df.values,
                                      coords={ind[0]: df.index.values, },
                                      dims=[ind[0]],
                                      name=var_name)
                all_ds.append(ds)
            else:
                pass  # print('Ignoring ', var_name)
        results = xr.merge(all_ds)
        return results
    except ImportError:
        logging.warn('Cannot find xarray package. Falling back to dict of dataframes')


def inputs_to_xarray(inputs):
    """Read crude input dispaset dictionary and return a structured xarray dataset

    :param inputs:  Dispa-SET inputs
    :return:        Inputs in xarray format
    """

    try:
        import xarray as xr
        # build set dictionary {parameter:set}
        in_keys = {}
        for k, v in inputs['parameters'].items():
            in_keys[k] = tuple(v['sets'])
        # Remove config from dict and add it later as dataset attribute (metadata)
        __ = in_keys.pop('Config')
        config = inputs['config']
        version = inputs.get('version', '')
        all_ds = []
        # Iterate all and build values nad coordinates
        for k, v in in_keys.items():
            val = inputs['parameters'][k]['val']
            var_name = k
            ind = v
            coords = OrderedDict()
            # Build dictionary for the specific parameter {sets:set elements}
            for sett in ind:
                coords[sett] = inputs['sets'][sett]
            dims = coords.keys()

            ds = xr.DataArray(val,
                              coords=coords,
                              dims=dims,
                              name=var_name)

            all_ds.append(ds)
        inputs = xr.merge(all_ds)
        inputs.attrs['version'] = version
        for key in config:
            if isinstance(config[key], (float, int, str)):
                inputs.attrs[key] = config[key]
        # Replace h with DateTimeIndex
        StartDate = config['StartDate']
        StopDate = config['StopDate']  # last day of the simulation with look-ahead period
        StopDate_long = dt.datetime(*StopDate) + dt.timedelta(days=config['LookAhead'])
        index_long = pd.date_range(start=dt.datetime(*StartDate), end=StopDate_long, freq='h')
        inputs.coords['h'] = index_long

    except ImportError:
        logging.warn('Cannot find xarray package. Falling back to dict of dataframes')
    return inputs


def ds_to_df(inputs):
    """
    Function that converts the dispaset data format into a dictionary of dataframes

    :param inputs:      input file
    :return:            dictionary of dataframes
    """

    sets, parameters = inputs['sets'], inputs['parameters']

    # config = parameters['Config']['val']
    try:
        config = inputs['config']
        first_day = dt.datetime(config['StartDate'][0], config['StartDate'][1], config['StartDate'][2], 0)
        last_day = dt.datetime(config['StopDate'][0], config['StopDate'][1], config['StopDate'][2], 23)
        dates = pd.date_range(start=first_day, end=last_day,
                              freq=str(int(inputs['config']['SimulationTimeStep'])) + 'h')
        timeindex = True
    except KeyError:
        dates = range(1, len(sets['z']) + 1)
        timeindex = False
    if len(dates) > len(sets['h']):
        logging.error('The provided index has a length of ' + str(len(dates)) + ' while the data only comprises ' + str(
            len(sets['h'])) + ' time elements')
        sys.exit(1)
    elif len(dates) > len(sets['z']):
        logging.warn(
            'The provided index has a length of ' + str(len(dates)) + ' while the simulation was designed for ' + str(
                len(sets['z'])) + ' time elements')
    elif len(dates) < len(sets['z']):
        logging.warn(
            'The provided index has a length of ' + str(len(dates)) + ' while the simulation was designed for ' + str(
                len(sets['z'])) + ' time elements')

    idx = range(len(dates))

    out = {}
    out['sets'] = sets

    # Printing each parameter in a separate sheet and workbook:
    for p in parameters:
        var = parameters[p]
        dim = len(var['sets'])
        if var['sets'][-1] == 'h' and timeindex and dim > 1:
            # if len(dates) != var['val'].shape[-1]:
            #    sys.exit('The date range in the Config variable (' + str(len(dates)) + ' time steps) does not match the length of the time index (' + str(var['val'].shape[-1]) + ') for variable ' + p)
            var['firstrow'] = 5
        else:
            var['firstrow'] = 1
        if dim == 1:
            if var['sets'][0] == 'h':
                out[p] = pd.DataFrame(var['val'][idx], columns=[p], index=dates)
            else:
                out[p] = pd.DataFrame(var['val'], columns=[p], index=sets[var['sets'][0]])
        elif dim == 2:
            values = var['val']
            list_sets = [sets[var['sets'][0]], sets[var['sets'][1]]]
            if var['sets'][1] == 'h':
                out[p] = pd.DataFrame(values.transpose()[idx, :], index=dates, columns=list_sets[0])
            else:
                out[p] = pd.DataFrame(values.transpose(), index=list_sets[1], columns=list_sets[0])
        elif dim == 3:
            list_sets = [sets[var['sets'][0]], sets[var['sets'][1]], sets[var['sets'][2]]]
            values = var['val']
            values2 = np.zeros([len(list_sets[0]) * len(list_sets[1]), len(list_sets[2])])
            cols = np.zeros([2, len(list_sets[0]) * len(list_sets[1])])
            for i in range(len(list_sets[0])):
                values2[i * len(list_sets[1]):(i + 1) * len(list_sets[1]), :] = values[i, :, :]
                cols[0, i * len(list_sets[1]):(i + 1) * len(list_sets[1])] = i
                cols[1, i * len(list_sets[1]):(i + 1) * len(list_sets[1])] = range(len(list_sets[1]))

            columns = pd.MultiIndex([list_sets[0], list_sets[1]], cols)
            if var['sets'][2] == 'h':
                out[p] = pd.DataFrame(values2.transpose()[idx, :], index=dates, columns=columns)
            else:
                out[p] = pd.DataFrame(values2.transpose(), index=list_sets[2], columns=columns)
        else:
            logging.error(
                'Only three dimensions currently supported. Parameter ' + p + ' has ' + str(dim) + ' dimensions.')
            sys.exit(1)
    return out
