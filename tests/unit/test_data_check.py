# -*- coding: utf-8 -*-
"""
Unit test: dispaset.preprocessing.data_check
============================================

What this test does
-------------------
Dispa-SET ships a battery of `check_*` functions that validate the input
power plant table, availability factors, flexibility data, and so on.
They are critical for surfacing user errors early and they all rely on
``sys.exit(1)`` to abort.

We test the most important paths:

* `check_units` raises ``SystemExit`` when:
    - the unit names contain duplicate keys.
    - a NonNaN key contains a string instead of a number.
    - a NonNaN key is missing.
* `check_AvailabilityFactors` does not raise on a clean dataset.
* `check_FlexibleDemand` raises if the values are outside [0, 1].
* `check_clustering` accepts identical capacity totals.

How to run
----------

1. ``pytest tests/unit/test_data_check.py``
2. ``python tests/unit/test_data_check.py``

This test does NOT require GAMS.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from dispaset.preprocessing.data_check import (  # noqa: E402
    check_units, check_AvailabilityFactors, check_FlexibleDemand,
    check_clustering,
)
from dispaset.common import commons  # noqa: E402


BUGGY_DIR = Path(__file__).resolve().parents[1] / "data" / "buggy"
DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def _config_with_horizon(h: int = 7):
    return {
        "HorizonLength": h,
        "SimulationTimeStep": 1,
        "SimulationType": "LP clustered",
        "HydroScheduling": "Off",
    }


# --------------------------------------------------------------------------- #
# check_units
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("buggy_csv", [
    "Units_bugNonNaNKeys1.csv",
    "Units_bugNonNaNKeys2.csv",
    "Units_bugKeys.csv",
])
def test_check_units_buggy_inputs_exit(buggy_csv):
    """A buggy power plant CSV must trigger SystemExit in check_units."""
    plants = pd.read_csv(BUGGY_DIR / buggy_csv)
    plants.set_index("Unit", drop=False, inplace=True)
    config = _config_with_horizon(7)
    with pytest.raises(SystemExit):
        check_units(config, plants)


def test_check_units_clean_inputs_pass():
    """The reference power plant CSV should pass `check_units` cleanly.

    The reference CSV uses ``NoLoadCost_pu`` rather than ``NoLoadCost`` and
    omits ``InertiaConstant``/``RampingCost``. These columns are normally
    derived in :func:`build_single_run`. Here we re-create the same minimal
    transformation so the table can be validated in isolation.
    """
    plants = pd.read_csv(DATA_DIR / "Units_testcase.csv")
    plants = plants[~plants["Unit"].str.startswith("Bug ")].copy()
    plants.set_index("Unit", drop=False, inplace=True)

    if "NoLoadCost" not in plants.columns and "NoLoadCost_pu" in plants.columns:
        plants["NoLoadCost"] = plants["NoLoadCost_pu"].fillna(0) * plants["PowerCapacity"]
    for col in ("InertiaConstant", "RampingCost"):
        if col not in plants.columns:
            plants[col] = 0

    config = _config_with_horizon(96)
    assert check_units(config, plants) is True


# --------------------------------------------------------------------------- #
# check_AvailabilityFactors
# --------------------------------------------------------------------------- #
def test_check_availability_factors_handles_renewables():
    """`check_AvailabilityFactors` must not raise for a valid renewable AF."""
    units = pd.DataFrame({
        "Unit":       ["WT1", "PV1"],
        "Technology": ["WTON", "PHOT"],
    })
    af = pd.DataFrame({
        "WT1": [0.1, 0.2, 0.3],
        "PV1": [0.0, 0.5, 0.9],
    })
    # Should run silently:
    check_AvailabilityFactors(units, af)


def test_check_availability_factors_raises_on_missing_renewable():
    """A renewable unit with no AF row should raise a ValueError."""
    units = pd.DataFrame({
        "Unit":       ["WT1"],
        "Technology": ["WTON"],
    })
    af = pd.DataFrame({"OTHER": [0.1, 0.2]})
    with pytest.raises(ValueError):
        check_AvailabilityFactors(units, af)


# --------------------------------------------------------------------------- #
# check_FlexibleDemand
# --------------------------------------------------------------------------- #
def test_check_flexible_demand_within_bounds():
    """Values inside [0, 1] should pass."""
    flex = pd.DataFrame({"Z1": [0.1, 0.2, 0.5]})
    check_FlexibleDemand(flex)


@pytest.mark.parametrize("bad_value", [-0.1, 1.1])
def test_check_flexible_demand_outside_bounds_exits(bad_value):
    """A value outside [0, 1] should trigger SystemExit."""
    flex = pd.DataFrame({"Z1": [0.0, bad_value, 0.5]})
    with pytest.raises(SystemExit):
        check_FlexibleDemand(flex)


# --------------------------------------------------------------------------- #
# check_clustering
# --------------------------------------------------------------------------- #
def test_check_clustering_accepts_equal_capacity():
    """When the merged capacity equals the original, the check passes."""
    plants = pd.DataFrame({
        "Unit": ["A", "B"], "Technology": ["GTUR", "GTUR"], "Fuel": ["GAS", "GAS"],
        "PowerCapacity": [100, 200], "Nunits": [1, 1], "STOCapacity": [0, 0],
    })
    plants_merged = pd.DataFrame({
        "Unit": ["AB"], "Technology": ["GTUR"], "Fuel": ["GAS"],
        "PowerCapacity": [300], "Nunits": [1], "STOCapacity": [0],
    })
    assert check_clustering(plants, plants_merged) is True


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
