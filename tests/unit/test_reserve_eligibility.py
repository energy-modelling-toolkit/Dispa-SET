# -*- coding: utf-8 -*-
"""
Unit tests for the frequency-stability reserve-eligibility fixes described in
UCED_ISSUES.md (ISSUE-UCED-001) and UCED_DECISIONS.md (DECISION-UCED-001/002):

* ``compute_quick_start_eligible`` (dispaset/preprocessing/reserves.py) -- the Python
  replacement for the previously GAMS-hardcoded ``QuickStartEligible(res,au)`` assignment.
* The technology-based default ``StartUpTime`` fallback in
  ``dispaset/preprocessing/data_check.py::check_units`` -- replaces the previous blind
  ``fillna(0)`` (in ``utils.py::clustering``) that made every unit with missing data
  unconditionally fast-reserve eligible. Moved into ``check_units`` (mirroring its existing
  "Nunits missing -> default to 1" precedent) so the default applies regardless of clustering
  method, and runs before clustering() instead of inside it.

How to run
----------
1. ``pytest tests/unit/test_reserve_eligibility.py``
2. ``python tests/unit/test_reserve_eligibility.py``
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from dispaset.common import commons  # noqa: E402
from dispaset.preprocessing.reserves import compute_quick_start_eligible  # noqa: E402
from dispaset.preprocessing.data_check import check_units  # noqa: E402

RES = ['VIRU', 'FFRU', 'FCRU', 'aFRRU', 'mFRRU', 'FFRD', 'FCRD', 'aFRRD']
RES_U = ['VIRU', 'FFRU', 'FCRU', 'aFRRU', 'mFRRU']


# --------------------------------------------------------------------------- #
# compute_quick_start_eligible
# --------------------------------------------------------------------------- #
def test_battery_eligible_for_every_reserve_regardless_of_startup_time():
    eligible = compute_quick_start_eligible(
        res=RES, res_U=RES_U, au=['BAT1'], ba=['BAT1'],
        time_startup={'BAT1': 999})  # even an absurd start-up time must not matter
    assert all(eligible[(r, 'BAT1')] == 1 for r in RES)


def test_zero_startup_time_is_eligible_for_all_upward_reserves():
    # A real-world data case: StartUpTime == 0 for
    # (previously) every non-battery unit, which the formula treats as "no start-up
    # penalty" -- eligible for everything upward, but NOT for downward reserves.
    eligible = compute_quick_start_eligible(
        res=RES, res_U=RES_U, au=['STUR1'], ba=[], time_startup={'STUR1': 0})
    assert all(eligible[(r, 'STUR1')] == 1 for r in RES_U)
    assert all(eligible[(r, 'STUR1')] == 0 for r in RES if r not in RES_U)


def test_slow_thermal_unit_not_eligible_for_fast_products_when_startup_time_is_realistic():
    # A steam turbine with a realistic 4h start-up time must not be credited with the
    # ability to deliver FCR (15s) or aFRR (5min) while offline -- this is exactly the
    # physical-realism problem found in real data
    # (where StartUpTime was 0 for every thermal unit).
    eligible = compute_quick_start_eligible(
        res=RES, res_U=RES_U, au=['STUR1'], ba=[],
        time_startup={'STUR1': commons['default_StartUpTime']['STUR']})
    assert eligible[('FFRU', 'STUR1')] == 0
    assert eligible[('FCRU', 'STUR1')] == 0
    assert eligible[('aFRRU', 'STUR1')] == 0
    assert eligible[('mFRRU', 'STUR1')] == 0


def test_fast_thermal_unit_eligible_only_up_to_its_own_startup_time():
    # A gas turbine (default 0.17h = ~10min) clears FCR (15s), aFRR (5min)?? aFRRU full
    # activation time is 5min == 0.0833h < 0.17h -- so a 0.17h start-up unit should NOT
    # clear aFRRU or mFRRU, only... check against the actual FullActivationTime table.
    fat = commons['FullActivationTime']
    tsu = 0.17
    eligible = compute_quick_start_eligible(
        res=RES, res_U=RES_U, au=['GTUR1'], ba=[], time_startup={'GTUR1': tsu})
    for r in RES_U:
        if r == 'VIRU':
            continue  # VIRU has no FullActivationTime entry (defaults to 0, like in GAMS)
        expected = 1 if tsu <= fat[r] else 0
        assert eligible[(r, 'GTUR1')] == expected, r


def test_unit_missing_from_time_startup_mapping_defaults_to_zero_like_gams():
    # A unit absent from the mapping (e.g. a data-alignment bug) must behave like
    # TimeStartUp=0, matching the previous GAMS default and the preprocessing fallback --
    # not silently raise or silently become ineligible.
    eligible = compute_quick_start_eligible(
        res=RES, res_U=RES_U, au=['GHOST1'], ba=[], time_startup={})
    assert all(eligible[(r, 'GHOST1')] == 1 for r in RES_U)


def test_reserve_duration_and_full_activation_time_defaults_are_registered():
    # Sanity check that the single-source-of-truth tables moved into commons.py still
    # cover every reserve category referenced by UCM.gms's res_U/res_D sets.
    for r in RES_U + ['FFRD', 'FCRD', 'aFRRD']:
        if r == 'VIRU':
            continue
        assert r in commons['ReserveDuration']
        assert r in commons['FullActivationTime']


# --------------------------------------------------------------------------- #
# Technology-based default StartUpTime (data_check.py::check_units)
#
# Deliberately independent of clustering(): check_units runs before clustering() on the raw
# plants table (build.py:250 vs. build.py:580) and must apply the default regardless of
# SimulationType/clustering method, mirroring the existing "Nunits missing -> default to 1"
# precedent already in check_units (data_check.py:438-441).
# --------------------------------------------------------------------------- #
def _minimal_clean_units(startup_time=0.0, technology="GTUR"):
    """Minimal valid power-plants DataFrame accepted by check_units (mirrors
    tests/failure/test_invalid_units_data.py::_minimal_clean_units)."""
    return pd.DataFrame([{
        "Unit": "U1", "PowerCapacity": 100.0, "Nunits": 1,
        "Zone": "Z1", "Technology": technology, "Fuel": "GAS",
        "Efficiency": 0.4, "MinUpTime": 1, "MinDownTime": 1,
        "RampUpRate": 1.0, "RampDownRate": 1.0,
        "StartUpCost": 0.0, "NoLoadCost": 0.0,
        "RampingCost": 0.0, "PartLoadMin": 0.0,
        "MinEfficiency": 0.0, "StartUpTime": startup_time,
        "CO2Intensity": 0.2, "CHPType": "", "CHPPowerToHeat": 0,
        "CHPPowerLossFactor": 0, "CHPMaxHeat": 0,
        "STOCapacity": 0, "STOSelfDischarge": 0,
        "STOMaxChargingPower": 0, "STOChargingEfficiency": 0,
        "InertiaConstant": 0, "WaterWithdrawal": 0, "WaterConsumption": 0,
    }])


def _config():
    return {"HorizonLength": 1, "DataTimeStep": 1.0,
            "SimulationTimeStep": 1.0, "SimulationType": "LP clustered"}


def test_missing_startup_time_uses_technology_default_not_zero(caplog):
    plants = _minimal_clean_units(startup_time=np.nan, technology="STUR")
    check_units(_config(), plants)
    assert plants.loc[0, "StartUpTime"] == commons['default_StartUpTime']['STUR']
    assert any("U1" in r.message and "default" in r.message.lower() for r in caplog.records)


def test_explicit_zero_startup_time_is_never_overridden():
    # An explicitly-entered 0 (real modeler data) must be left
    # untouched -- only genuinely missing (NaN) values get a technology default.
    plants = _minimal_clean_units(startup_time=0.0, technology="STUR")
    check_units(_config(), plants)
    assert plants.loc[0, "StartUpTime"] == 0.0


def test_missing_startup_time_for_unlisted_technology_uses_fallback_and_warns(caplog):
    plants = _minimal_clean_units(startup_time=np.nan, technology="NOT_A_REAL_TECH")
    check_units(_config(), plants)
    assert plants.loc[0, "StartUpTime"] == commons['default_StartUpTime_fallback']
    assert any("no registered default" in r.message for r in caplog.records)


def test_battery_and_renewable_defaults_are_at_or_near_zero():
    for tech in ['BATS', 'BEVS', 'PHOT', 'WTON']:
        assert commons['default_StartUpTime'][tech] == 0


def test_thermal_defaults_are_strictly_positive_and_ordered_by_flexibility():
    d = commons['default_StartUpTime']
    assert 0 < d['ICEN'] < d['GTUR'] < d['COMC'] < d['STUR']


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
