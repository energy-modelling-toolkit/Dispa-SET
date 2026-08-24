# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 12:02:10 2020

@author: Chiara Magni
"""

import pandas as pd

from ..common import commons


def compute_fast_reserve_eligible(res, res_U, au, ba, time_startup, full_activation_time=None):
    """
    Compute the FastReserveEligible(res, au) table: whether a unit can deliver a given
    upward reserve product while NOT committed (i.e. from an offline/cold state).

    A unit is eligible for an upward reserve category ``r`` if its start-up time is either
    0 (no data / instant-start) or does not exceed ``r``'s full-activation-time requirement.
    Battery units (``ba``) are always eligible for every reserve category (including
    downward ones), independent of start-up time.

    :param res:                  Iterable of all reserve category names (e.g. sets['res']).
    :param res_U:                Iterable of upward reserve category names (sets['res_U']).
    :param au:                   Iterable of all unit names (sets['au']).
    :param ba:                   Iterable of battery unit names (sets['ba']).
    :param time_startup:         Mapping (dict or pandas Series) from unit name to
                                 TimeStartUp [h]. A unit absent from the mapping is treated
                                 as 0, matching the "no data" convention used upstream.
    :param full_activation_time: Mapping {res: hours}; defaults to
                                 commons['FullActivationTime']. Pass
                                 config.get('FullActivationTime', commons['FullActivationTime'])
                                 to respect a case-level override.
    :return:                     dict {(r, u): 0 or 1} for every (r, u) in res x au.
    """
    if full_activation_time is None:
        full_activation_time = commons['FullActivationTime']
    ba_set = set(ba)
    eligible = {}
    for u in au:
        time_startup_u = time_startup.get(u, 0)
        for r in res:
            eligible[(r, u)] = 0
        for r in res_U:
            fat = full_activation_time.get(r, 0)
            if time_startup_u == 0 or (0 < time_startup_u <= fat):
                eligible[(r, u)] = 1
        if u in ba_set:
            for r in res:
                eligible[(r, u)] = 1
    return eligible


def percentage_reserve(config, allunits, load, AvailabilityFactors, zone):
    """
    Reserve Requirements : (3+5) rule

    dynamic method for the evaluation of reserve requirements based on forecasted
    load and forecasted available wind and solar power.

    :param config:              Dictionary with all the configuration fields loaded from the excel file. Output of the
                                'LoadConfig' function.
    :param allunits:            Dataframe of all power plants
    :param load:                Dataframe of load timeseries
    :param AvailabilityFactors: Dataframe of availability factors
    :param zone:                Zone for reserve calculation
    """

    c = zone
    func = lambda col: col * units.loc[col.name, 'PowerCapacity'] * units.loc[col.name, 'Nunits']

    # 3% of the load forecast
    load = load[c]

    units = allunits[allunits['Zone'] == c]

    # 5% of the windon forecast
    if units['Technology'].str.contains('WTON').any():
        tmp = units.index[units['Technology'] == 'WTON']
        wton = AvailabilityFactors[tmp].agg(func).sum(axis=1)
    else:
        wton = pd.Series(0, index=load.index)

    # 5% of the windoff forecast

    if units['Technology'].str.contains('WTOF').any():
        tmp = units.index[units['Technology'] == 'WTOF']
        wtof = AvailabilityFactors[tmp].agg(func).sum(axis=1)
    else:
        wtof = pd.Series(0, index=load.index)

    # 5% of solar forecast

    if units['Technology'].str.contains('PHOT').any():
        tmp = units.index[units['Technology'] == 'PHOT']
        phot = AvailabilityFactors[tmp].agg(func).sum(axis=1)
    else:
        phot = pd.Series(0, index=load.index)

    # total reserve demand
    rr = 0.03 * load + 0.05 * wton + 0.05 * wtof + 0.05 * phot

    return rr, rr


def probabilistic_reserve(config, allunits, load, AvailabilityFactors, zone,
                          std={'WTON': 0.2, 'WTOF': 0.2, 'PHOT': 0.045, 'LOAD': 0.02}):
    """
    dynamic method for the evaluation of reserve requirements based on load
    (secondary reserves requiements) and error forecasting for load, wind and
    solar power (tertiary reserves requirements).

    :param config:              Dictionary with all the configuration fields loaded from the excel file. Output of the
                                'LoadConfig' function.
    :param allunits:            Dataframe of all power plants
    :param load:                Dataframe of load timeseries
    :param AvailabilityFactors: Dataframe of availability factors
    :param zone:                Zone for reserve calculation
    :param std:                 Dictionary with multiplication factors for the calculation of the standard deviation
    """

    c = zone
    func = lambda col: col * units.loc[col.name, 'PowerCapacity'] * units.loc[col.name, 'Nunits']

    load = load[c]

    units = allunits[allunits['Zone'] == c]

    # 5% of the windon forecast
    if units['Technology'].str.contains('WTON').any():
        tmp = units.index[units['Technology'] == 'WTON']
        wton = AvailabilityFactors[tmp].agg(func).sum(axis=1)
    else:
        wton = pd.Series(0, index=load.index)

    # 5% of the windoff forecast
    if units['Technology'].str.contains('WTOF').any():
        tmp = units.index[units['Technology'] == 'WTOF']
        wtof = AvailabilityFactors[tmp].agg(func).sum(axis=1)
    else:
        wtof = pd.Series(0, index=load.index)

    # 5% of solar forecast
    if units['Technology'].str.contains('PHOT').any():
        tmp = units.index[units['Technology'] == 'PHOT']
        phot = AvailabilityFactors[tmp].agg(func).sum(axis=1)
    else:
        phot = pd.Series(0, index=load.index)

    # Standard deviations
    load_std = load * std.get('LOAD', 0)
    wton_std = wton * std.get('WTON', 0)
    wtof_std = wtof * std.get('WTOF', 0)
    phot_std = phot * std.get('PHOT', 0)

    srr = (10 * load + 150 ** 2) ** 0.5 - 150
    trr = 2.74 * (load_std ** 2 + (wton_std + wtof_std) ** 2 + phot_std ** 2) ** 0.5
    rr = srr + trr

    return rr, rr


def generic_reserve(load):
    """
    generic reserve sizing based on the load curve.
    Not adapted to high shares of renewables

    :param load:                Dataframe of load timeseries
    """
    # up = (10 * load.max() + 150 ** 2) ** 0.5 - 150
    # down = 0.5 * up
    up = load.apply(lambda col: ((10 * col + 150 ** 2) ** 0.5) - 150)
    down = 0.5 * up   
    return pd.Series(up, index=load.index), pd.Series(down, index=load.index)
