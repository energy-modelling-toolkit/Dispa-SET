# -*- coding: utf-8 -*-
"""
Unit test: dispaset.postprocessing filtering helpers
====================================================

What this test does
-------------------
The postprocessing module exposes a series of small filter helpers that
take a result DataFrame (one column per unit) and return a sub-frame
filtered by zone, technology, fuel, etc. They are pure-Python and can be
tested in isolation with synthetic dictionaries.

Covered:

* `filter_by_zone`     - filter by zone name (string match on `units['Zone']`).
* `filter_by_tech`     - filter by single technology code.
* `filter_by_tech_list`- filter by a list of technology codes.
* `aggregate_by_fuel`  - aggregate per-unit power by fuel type.
* `get_plot_data`      - prepare per-zone dispatch plot data.

How to run
----------

1. ``pytest tests/unit/test_postprocess_filters.py``
2. ``python tests/unit/test_postprocess_filters.py``

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

from dispaset.postprocessing.postprocessing import (  # noqa: E402
    filter_by_zone, filter_by_tech, filter_by_tech_list,
    aggregate_by_fuel, get_plot_data,
)


def _make_inputs_results(periods=24):
    """Create a synthetic mini Dispa-SET input/results dict."""
    idx = pd.date_range("2015-01-01", periods=periods, freq="h")
    units = pd.DataFrame({
        "Zone":       ["Z1", "Z1", "Z2"],
        "Technology": ["GTUR", "WTON", "GTUR"],
        "Fuel":       ["GAS",  "WIN",  "GAS"],
    }, index=["Z1 - GAS", "Z1 - WIN", "Z2 - GAS"])

    inputs = {
        "sets": {
            "n": ["Z1", "Z2"],
            "f": ["GAS", "WIN"],
            "t": ["GTUR", "WTON"],
        },
        "units": units,
        "param_df": {
            "Demand": pd.DataFrame(
                {("DA", "Z1"): np.full(periods, 100.0),
                 ("DA", "Z2"): np.full(periods, 50.0)},
                index=idx,
            ),
            "Technology": pd.DataFrame(  # tech-by-unit indicator
                np.eye(3, dtype=int),
                index=["GTUR", "WTON", "GTUR_2"],
                columns=units.index,
            ),
        },
    }
    results = {
        "OutputPower": pd.DataFrame(
            {"Z1 - GAS": np.full(periods, 70.0),
             "Z1 - WIN": np.full(periods, 30.0),
             "Z2 - GAS": np.full(periods, 50.0)},
            index=idx,
        ),
    }
    return inputs, results


def test_filter_by_zone_returns_only_unit_in_zone():
    inputs, results = _make_inputs_results()
    out = filter_by_zone(results["OutputPower"], inputs, "Z1")
    assert sorted(out.columns) == ["Z1 - GAS", "Z1 - WIN"]


def test_filter_by_tech():
    inputs, results = _make_inputs_results()
    out = filter_by_tech(results["OutputPower"], inputs, "GTUR")
    assert sorted(out.columns) == ["Z1 - GAS", "Z2 - GAS"]


def test_filter_by_tech_list():
    inputs, results = _make_inputs_results()
    out = filter_by_tech_list(results["OutputPower"], inputs, ["GTUR", "WTON"])
    assert sorted(out.columns) == sorted(["Z1 - GAS", "Z1 - WIN", "Z2 - GAS"])


def test_aggregate_by_fuel():
    inputs, results = _make_inputs_results()
    by_fuel = aggregate_by_fuel(results["OutputPower"], inputs)
    assert isinstance(by_fuel, pd.DataFrame)
    # GAS is generated in both zones (70 + 50 = 120)
    assert by_fuel["GAS"].iloc[0] == pytest.approx(70 + 50)
    assert by_fuel["WIN"].iloc[0] == pytest.approx(30)


def test_get_plot_data_returns_non_empty():
    inputs, results = _make_inputs_results()
    plot_data = get_plot_data(inputs, results, "Z1")
    assert isinstance(plot_data, pd.DataFrame)
    assert not plot_data.empty
    assert "Demand" not in plot_data.columns
    # WIN/GAS columns should be present (one per fuel/tech depending on impl)
    assert any(col in plot_data.columns for col in ["GAS", "WIN", "GTUR", "WTON"])


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
