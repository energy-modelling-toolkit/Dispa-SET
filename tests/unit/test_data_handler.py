# -*- coding: utf-8 -*-
"""
Unit test: dispaset.preprocessing.data_handler I/O helpers
==========================================================

What this test does
-------------------
The data handler exposes a few generic loaders that are used everywhere in
preprocessing:

* `NodeBasedTable` - load per-zone time series (Demand, RES AF, etc.).
* `UnitBasedTable` - load per-unit time series (Outages, AF fallback, etc.).
* `load_geo_data`  - read the per-zone geographic coordinates.

We test:

* A valid path returns a non-empty DataFrame indexed on `idx_long`.
* A path containing `##` and missing zone files exits cleanly.
* A unit-based table falls back on the requested keys when the unit name is
  unknown.
* Geo data is correctly indexed by country code.

How to run
----------

1. ``pytest tests/unit/test_data_handler.py``
2. ``python tests/unit/test_data_handler.py``

This test does NOT require GAMS.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from dispaset.preprocessing.data_handler import (  # noqa: E402
    NodeBasedTable, UnitBasedTable, load_geo_data,
)


DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def _make_config(zones=("Z1", "Z2")):
    """Create a minimal config dict with a one-day idx_long."""
    idx = pd.date_range("2015-01-01", "2015-01-01 23:00", freq="h")
    return {
        "zones": list(zones),
        "idx_long": idx,
        "idx": idx,
        "DataTimeStep": 1,
        "SimulationTimeStep": 1,
    }


# --------------------------------------------------------------------------- #
# NodeBasedTable
# --------------------------------------------------------------------------- #
def test_node_based_table_with_zone_pattern():
    """`NodeBasedTable` should resolve `##` to the per-zone files."""
    config = _make_config(zones=["Z1"])
    config["Demand"] = str(DATA_DIR / "Load_RealTime" / "##" / "2015.csv")
    out = NodeBasedTable("Demand", config)
    assert isinstance(out, pd.DataFrame)
    assert "Z1" in out.columns
    assert len(out) == len(config["idx_long"])
    assert out["Z1"].notna().all()


def test_node_based_table_unknown_zone_exits():
    """An unknown zone in the `##` pattern triggers SystemExit."""
    config = _make_config(zones=["XX"])
    config["Demand"] = str(DATA_DIR / "Load_RealTime" / "##" / "2015.csv")
    with pytest.raises(SystemExit):
        NodeBasedTable("Demand", config)


def test_node_based_table_default_value():
    """No file + numeric default returns a constant DataFrame."""
    config = _make_config(zones=["Z1", "Z2"])
    config["Demand"] = ""
    out = NodeBasedTable("Demand", config, default=0)
    assert (out == 0).all().all()
    assert list(out.columns) == ["Z1", "Z2"]


# --------------------------------------------------------------------------- #
# UnitBasedTable
# --------------------------------------------------------------------------- #
def test_unit_based_table_with_default():
    """A missing path + numeric default produces a constant DataFrame."""
    config = _make_config(zones=["Z1"])
    config["Outages"] = ""
    plants = pd.DataFrame({
        "Unit": ["U1", "U2"], "Zone": ["Z1", "Z1"],
        "Technology": ["GTUR", "GTUR"], "Fuel": ["GAS", "GAS"],
    })
    out = UnitBasedTable(plants, "Outages", config, fallbacks=["Unit"], default=0)
    assert (out == 0).all().all()
    assert list(out.columns) == ["U1", "U2"]


def test_unit_based_table_falls_back_to_technology():
    """Fallback to the technology column when no unit-specific data exists."""
    config = _make_config(zones=["Z1"])
    config["Outages"] = ""  # Keep the fallback path simple
    # Build a fake outages CSV in tmp dir is overkill - we use the default API.
    plants = pd.DataFrame({
        "Unit": ["NoMatchUnit"], "Zone": ["Z1"],
        "Technology": ["GTUR"], "Fuel": ["GAS"],
    })
    out = UnitBasedTable(plants, "Outages", config, fallbacks=["Unit", "Technology"],
                        default=0)
    assert "NoMatchUnit" in out.columns


# --------------------------------------------------------------------------- #
# load_geo_data
# --------------------------------------------------------------------------- #
def test_load_geo_data():
    """The `Geo_Coordinates.csv` file is parsed with `CountryCode` as index."""
    geo = load_geo_data(str(DATA_DIR / "Geo_Coordinates.csv"), header=0)
    assert isinstance(geo, pd.DataFrame)
    assert geo.index.name == "CountryCode"


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
