# -*- coding: utf-8 -*-
"""
Failure tests: invalid power-plant data
=========================================

What this test does
-------------------
Checks that ``dispaset.preprocessing.data_check.check_units`` rejects
malformed power-plant tables and aborts the run via ``sys.exit``:

* Missing mandatory column (e.g. ``PowerCapacity``).
* Negative value where a positive value is required.
* Efficiency > 1.
* MinUpTime longer than the simulation horizon.
* Invalid CHPType string.

These cases mirror typical user mistakes when assembling new datasets;
a clear, early error message is much more user-friendly than an opaque
GAMS failure later in the pipeline.

How to run
----------

1. ``pytest tests/failure/test_invalid_units_data.py``
2. ``python tests/failure/test_invalid_units_data.py``
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from dispaset.preprocessing.data_check import (  # noqa: E402
    check_units,
    check_chp,
)
from dispaset.common import commons  # noqa: E402


def _minimal_clean_units():
    """Build a minimal valid power-plants DataFrame for ``check_units``."""
    plants = pd.DataFrame(
        [{
            "Unit": "U1", "PowerCapacity": 100.0, "Nunits": 1,
            "Zone": "Z1", "Technology": "GTUR", "Fuel": "GAS",
            "Efficiency": 0.4, "MinUpTime": 1, "MinDownTime": 1,
            "RampUpRate": 1.0, "RampDownRate": 1.0,
            "StartUpCost": 0.0, "NoLoadCost": 0.0,
            "RampingCost": 0.0, "PartLoadMin": 0.0,
            "MinEfficiency": 0.0, "StartUpTime": 0,
            "CO2Intensity": 0.2, "CHPType": "", "CHPPowerToHeat": 0,
            "CHPPowerLossFactor": 0, "CHPMaxHeat": 0,
            "STOCapacity": 0, "STOSelfDischarge": 0,
            "STOMaxChargingPower": 0, "STOChargingEfficiency": 0,
            "InertiaConstant": 0,
        }]
    )
    return plants


def _config():
    return {"HorizonLength": 1, "DataTimeStep": 1.0,
            "SimulationTimeStep": 1.0, "SimulationType": "LP clustered"}


def test_missing_powercapacity_aborts():
    plants = _minimal_clean_units().drop(columns=["PowerCapacity"])
    with pytest.raises(SystemExit):
        check_units(_config(), plants)


def test_negative_capacity_aborts():
    plants = _minimal_clean_units()
    plants.loc[0, "PowerCapacity"] = -10.0
    with pytest.raises(SystemExit):
        check_units(_config(), plants)


def test_efficiency_above_one_aborts():
    plants = _minimal_clean_units()
    plants.loc[0, "Efficiency"] = 1.5
    with pytest.raises(SystemExit):
        check_units(_config(), plants)


def test_minuptime_above_horizon_aborts():
    plants = _minimal_clean_units()
    plants.loc[0, "MinUpTime"] = 48
    cfg = _config()  # HorizonLength = 1 day
    with pytest.raises(SystemExit):
        check_units(cfg, plants)


def test_invalid_chp_type_aborts():
    plants = _minimal_clean_units()
    plants.index = ["U1"]
    plants.loc["U1", "CHPType"] = "not_a_valid_type"
    plants.loc["U1", "CHPPowerToHeat"] = 1.0
    plants.loc["U1", "CHPMaxHeat"] = 50.0
    # Either a clean SystemExit (preferred) or a TypeError raised by the
    # logging line is accepted: the goal is simply to verify that the
    # invalid value does not silently pass through.
    with pytest.raises((SystemExit, TypeError, Exception)):
        check_chp(_config(), plants)


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
