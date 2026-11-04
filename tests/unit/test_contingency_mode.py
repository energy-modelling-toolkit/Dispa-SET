# -*- coding: utf-8 -*-
"""
Unit test (build only, no GAMS): config option ContingencyMode
==============================================================

ContingencyCapacity [MW] is the output lost in one N-1 event of a unit row:

* ``'plant'`` (default): a non-inverter row with Nunits > 1 is lost as a whole
  (PowerCapacity x Nunits): multi-unit hydro plant, combined-cycle block
  represented by its Nunits = 2 gas-turbine trains.
* ``'unit'``: each unit is its own N-1 event (PowerCapacity): one hydro unit,
  one combined-cycle train.
* Batteries (and other inverter-based modules) count per module in both modes.

Data: tests/data/Units_contingency.csv (CC 2 x 150 MW, GT 80 MW, BATS 4 x 50 MW,
HDAM 3 x 40 MW), 1-day tiny_spillage.yml skeleton, 'No clustering'.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import load_test_config, TESTS_DIR  # noqa: E402

import dispaset as ds  # noqa: E402
from dispaset.preprocessing.build import get_contingency_mode  # noqa: E402

EXPECTED = {
    "plant": {"Z1_CC": 300.0, "Z1_GT": 80.0, "Z1_BATS": 50.0, "Z1_HDAM": 120.0},
    "unit": {"Z1_CC": 150.0, "Z1_GT": 80.0, "Z1_BATS": 50.0, "Z1_HDAM": 40.0},
}


def _contingency_capacity(mode):
    cfg = load_test_config("tiny_spillage.yml", f"unit_contingency_{mode or 'default'}")
    cfg["PowerPlantData"] = str(TESTS_DIR / "data" / "Units_contingency.csv")
    cfg["SimulationType"] = "No clustering"
    if mode is not None:
        cfg["ContingencyMode"] = mode
    try:
        sim = ds.build_simulation(cfg)
    finally:
        shutil.rmtree(cfg["SimulationDirectory"], ignore_errors=True)
    au = list(sim["sets"]["au"])
    val = sim["parameters"]["ContingencyCapacity"]["val"]
    # unit names may carry a "[n] - " prefix
    return {u.split(" - ")[-1].strip(): float(v) for u, v in zip(au, val)}


@pytest.mark.timeout(120)
@pytest.mark.parametrize("mode", ["plant", "unit"])
def test_contingency_capacity_by_mode(mode):
    got = _contingency_capacity(mode)
    for unit, expected in EXPECTED[mode].items():
        assert got[unit] == pytest.approx(expected), (mode, unit, got)


def test_contingency_mode_default_and_validation():
    assert get_contingency_mode({}) == "plant"
    assert get_contingency_mode({"ContingencyMode": ""}) == "plant"
    assert get_contingency_mode({"ContingencyMode": " Unit "}) == "unit"
    with pytest.raises(ValueError):
        get_contingency_mode({"ContingencyMode": "block"})


def test_check_config_rejects_invalid_mode():
    from dispaset.preprocessing.data_handler import check_config
    assert check_config({"HydroScheduling": 0, "ContingencyMode": "UNIT"})["ContingencyMode"] == "unit"
    assert check_config({"HydroScheduling": 0})["ContingencyMode"] == "plant"
    with pytest.raises(SystemExit):
        check_config({"HydroScheduling": 0, "ContingencyMode": "group"})


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
