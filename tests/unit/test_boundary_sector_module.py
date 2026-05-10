# -*- coding: utf-8 -*-
"""
Unit test: dispaset.preprocessing.boundary_sector module
========================================================

What this test does
-------------------
The boundary-sector module groups the helpers that translate the legacy
"Heat" / "H2" formulations into the unified "Boundary Sector" abstraction.
The module is *new* and not exhaustively documented yet, so the tests below
serve a double purpose: validating behaviour AND documenting expectations.

Covered:

* `process_boundary_sector_data`: returns a dict with the expected dataframe
  keys and reads the input CSVs without crashing.
* `get_boundary_sector_zones`: extracts the boundary sector zones from a
  populated `BoundarySector` DataFrame and filters out NaN-like entries.
* `zone_to_bs_mapping`: returns a callable that maps a boundary sector zone
  back to its parent power zone.
* `boundary_sector_efficiency_time_series`: produces a per-unit efficiency
  time series indexed on `config['idx_long']`.

How to run
----------

1. ``pytest tests/unit/test_boundary_sector_module.py``
2. ``python tests/unit/test_boundary_sector_module.py``

This test does NOT require GAMS.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from dispaset.preprocessing.boundary_sector import (  # noqa: E402
    process_boundary_sector_data, get_boundary_sector_zones,
    zone_to_bs_mapping, boundary_sector_efficiency_time_series,
)


DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def _make_config():
    idx = pd.date_range("2015-01-01", "2015-01-01 23:00", freq="h").tz_localize(None)
    return {
        "idx_long": idx,
        "idx": idx,
        "DataTimeStep": 1,
        "SimulationTimeStep": 1,
        "BoundarySectorData": str(DATA_DIR / "boundary_sector" / "BS_Inputs.csv"),
        "BoundarySectorNTC":  str(DATA_DIR / "boundary_sector" / "BS_NTCs.csv"),
        "BoundarySectorInterconnections": str(DATA_DIR / "boundary_sector" / "BS_NTCs.csv"),
        "BoundarySectorMaxSpillage": str(DATA_DIR / "boundary_sector" / "BS_Spillage_Capacity.csv"),
        "CostXSpillage": str(DATA_DIR / "boundary_sector" / "BS_Spillage_cost.csv"),
    }


def test_process_boundary_sector_data_keys():
    """`process_boundary_sector_data` returns the expected set of keys."""
    cfg = _make_config()
    out = process_boundary_sector_data(cfg)
    expected = {"BoundarySector", "BS_flows", "BS_NTC",
                "BS_spillage", "BS_forced_spillage", "CostXSpillage"}
    assert expected.issubset(out.keys()), \
        f"Missing keys: {expected - out.keys()}"
    assert isinstance(out["BoundarySector"], pd.DataFrame)
    # The shipped BS_Inputs.csv defines at least one sector.
    assert not out["BoundarySector"].empty


def test_get_boundary_sector_zones_clean():
    """NaN/empty values must be discarded from the boundary sector zone list."""
    bs_df = pd.DataFrame(
        index=pd.Index(["Z1_h2", "Z2_th_2", "", np.nan, "nan"], name="Sector"),
        data={"STOCapacity": [1000, 0, 0, 0, 0]},
    )
    out = get_boundary_sector_zones(bs_df)
    assert sorted(out) == sorted(["Z1_h2", "Z2_th_2"])


def test_zone_to_bs_mapping():
    """Mapping should resolve `BS zone -> parent power zone`."""
    plants = pd.DataFrame({
        "Zone":    ["Z1", "Z2"],
        "Sector1": ["Z1_h2", "Z2_th_2"],
        "Sector2": [np.nan, np.nan],
    }, index=["unit_a", "unit_b"])
    mapping = zone_to_bs_mapping(plants)
    assert mapping("Z1_h2") == "Z1"
    assert mapping("Z2_th_2") == "Z2"
    assert mapping("unknown_sector") is None


def test_boundary_sector_efficiency_time_series_shape():
    """Output dict is indexed on `idx_long` with one column per unit."""
    cfg = _make_config()
    plants = pd.DataFrame({
        "Technology":        ["P2GS"],
        "Sector1":           ["Z1_h2"],
        "EfficiencySector1": [0.65],
    }, index=["electrolyzer_Z2"])
    eff = boundary_sector_efficiency_time_series(cfg, plants, ["Z1_h2"])
    assert "Efficiency" in eff
    assert "Z1_h2" in eff["Efficiency"]
    series = eff["Efficiency"]["Z1_h2"]
    assert len(series.index) == len(cfg["idx_long"])
    assert "electrolyzer_Z2" in series.columns
    assert (series["electrolyzer_Z2"] == 0.65).all()


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
