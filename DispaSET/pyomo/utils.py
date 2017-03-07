import sys

import numpy as np
import pandas as pd
import logging

def get_sets(instance, varname):
    """Get sets that belong to a pyomo Variable or Param

    :param instance: Pyomo Instance
    :param varname: Name of the Pyomo Variable (string)
    :return: A list with the sets that belong to this Param
    """
    var = getattr(instance, varname)

    if var.dim() > 1:
        sets = [pset.cname() for pset in var._index.set_tuple]
    else:
        sets = [var._index.name]
    return sets


def get_set_members(instance, sets):
    """Get set members relative to a list of sets

    :param instance: Pyomo Instance
    :param sets: List of strings with the set names
    :return: A list with the set members
    """
    sm = []
    for s in sets:
        sm.append([v for v in getattr(instance, s).value])
    return sm


def pyomo_to_pandas(instance, varname):
    """
    Function converting a pyomo variable or parameter into a pandas dataframe.
    The variable must have one or two dimensions and the sets must be provided as a list of lists

    :param instance: Pyomo Instance
    :param varname: Name of the Pyomo Variable (string)
    """
    setnames = get_sets(instance, varname)
    sets = get_set_members(instance, setnames)
    var = getattr(instance, varname)  # Previous script used model.var instead of var
    ####
    if len(sets) != var.dim():
        logging.error('The number of provided set lists (' + str(
            len(sets)) + ') does not match the dimensions of the variable (' + str(var.dim()) + ')')
        sys.exit(1)
    if var.dim() == 1:
        [SecondSet] = sets
        out = pd.DataFrame(columns=[var.name], index=SecondSet)
        data = var.get_values()
        for idx in data:
            out[var.name][idx] = data[idx]
        return out

    elif var.dim() == 2:
        [FirstSet, SecondSet] = sets
        out = pd.DataFrame(columns=FirstSet, index=SecondSet)
        data = var.get_values()
        for idx in data:
            out[idx[0]][idx[1]] = data[idx]
        return out
    else:
        logging.warn('the pyomo_to_pandas function currently only accepts one or two-dimensional variables')
        return []


def pyomo_format(sets, param):
    """
    Function that flattens the multidimensional dispaset input data into the pyomo format: a dictionnary with a tuple and the parameter value.
    The tuple contains the strings of the corresponding set values
    """

    param_new = {}
    ndims = len(param['sets'])
    for i in range(ndims):
        if param['val'].shape[i] != len(sets[param['sets'][i]]):
            logging.error('The size of the matrix and the number of set values do no match for set ' + param['sets'][i])
            sys.exit(1)
    if ndims == 1:
        for k in range(len(param['val'])):
            param_new[(sets[param['sets'][0]][k])] = param['val'][k]
    elif ndims == 2:
        X, Y = np.meshgrid(np.arange(param['val'].shape[1]), np.arange(param['val'].shape[0]))
        X_flat = X.flatten()
        Y_flat = Y.flatten()
        array_flat = param['val'].flatten()
        for k in range(len(array_flat)):
            i = Y_flat[k]
            j = X_flat[k]
            param_new[(sets[param['sets'][0]][i], sets[param['sets'][1]][j])] = array_flat[k]
    elif ndims == 3:
        X, Y, Z = np.meshgrid(np.arange(param['val'].shape[1]), np.arange(param['val'].shape[0]),
                              np.arange(param['val'].shape[2]))
        X_flat = X.flatten()
        Y_flat = Y.flatten()
        Z_flat = Z.flatten()
        array_flat = param['val'].flatten()
        for k in range(len(array_flat)):
            i = Y_flat[k]
            j = X_flat[k]
            z = Z_flat[k]
            param_new[(sets[param['sets'][0]][i], sets[param['sets'][1]][j], sets[param['sets'][2]][z])] = array_flat[k]
    else:
        logging.error('Maximum 3 dimensions')
        sys.exit(1)
    return param_new