# -*- coding: utf-8 -*-
"""
Unit test: dispaset.preprocessing.utils
========================================

What this test does
-------------------
The `utils` module collects many small helpers that drive the heavy lifting
of preprocessing. They are pure-Python and easy to test in isolation.

Covered helpers:

* `pd_timestep`           - converts an int/float into a pandas frequency.
* `_mylogspace`           - log-spaced numbers (used to slice piecewise
                            partial-load curves).
* `_find_nearest`         - find index of the nearest value.
* `_flatten_list`         - flatten a list with nested sub-lists.
* `_split_list`           - join a list of strings with " - ".
* `_merge_two_dicts`      - merge two dicts.
* `create_agg_dict`       - clustering helper that builds the aggregator dict.
* `select_units`          - filter the units table per zone list.

How to run
----------

1. ``pytest tests/unit/test_utils.py``
2. ``python tests/unit/test_utils.py``

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

from dispaset.preprocessing.utils import (  # noqa: E402
    pd_timestep, _mylogspace, _find_nearest, _flatten_list, _split_list,
    _merge_two_dicts, create_agg_dict, select_units,
)


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("step", [1, 0.25, 24])
def test_pd_timestep_valid(step):
    """`pd_timestep` should return a valid pandas-compatible frequency."""
    out = pd_timestep(step)
    assert out is not None


@pytest.mark.parametrize("step", ["24", None])
def test_pd_timestep_invalid_exits(step):
    """A non-numeric input should trigger a SystemExit (logged.critical)."""
    with pytest.raises(SystemExit):
        pd_timestep(step)


def test_mylogspace_endpoints():
    """`_mylogspace` must include the bounds."""
    out = _mylogspace(0, 1, 3)
    assert out[0] == pytest.approx(0)
    assert out[-1] == pytest.approx(1)
    assert len(out) == 3


def test_find_nearest():
    arr = np.array([0, 1, 2, 3, 4])
    assert _find_nearest(arr, 3.4) == 3


def test_flatten_list():
    assert _flatten_list([1, 3, ["aa", "bb"], 4]) == [1, 3, "aa", "bb", 4]


def test_split_list():
    assert _split_list(["Split", "me"]) == "Split - me"


def test_merge_two_dicts():
    assert _merge_two_dicts({"a": 1}, {"b": 2}) == {"a": 1, "b": 2}


# --------------------------------------------------------------------------- #
# Clustering
# --------------------------------------------------------------------------- #
def _read_clustering_plants() -> pd.DataFrame:
    path = Path(__file__).resolve().parents[1] / "data" / "buggy" / "Clustering.csv"
    plants = pd.read_csv(path)
    plants.set_index("Unit", inplace=True, drop=False)
    return plants


@pytest.mark.parametrize("method", ["LP clustered", "MILP", "Standard",
                                    "Integer clustering"])
def test_create_agg_dict_supported_modes(method):
    """All supported clustering modes should produce an aggregator dict."""
    plants = _read_clustering_plants()
    agg = create_agg_dict(plants, method)
    assert isinstance(agg, dict)


def test_create_agg_dict_no_clustering_exits():
    """The `'No Clustering'` mode is not supported by `create_agg_dict`."""
    plants = _read_clustering_plants()
    with pytest.raises(SystemExit):
        create_agg_dict(plants, "No Clustering")


# --------------------------------------------------------------------------- #
# select_units
# --------------------------------------------------------------------------- #
def test_select_units_filters_zones(tmp_path):
    """`select_units` should drop units whose zone is not in `config['zones']`."""
    csv = Path(__file__).resolve().parents[1] / "data" / "Units_testcase.csv"
    plants = pd.read_csv(csv)
    config = {"zones": ["Z1"]}
    out = select_units(plants, config)
    assert set(out["Zone"].unique()).issubset({"Z1"}), \
        f"Got non-Z1 zones in output: {out['Zone'].unique()}"
    assert not out.empty


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
