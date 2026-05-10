# -*- coding: utf-8 -*-
"""
Property-based tests: dispaset.preprocessing.utils
====================================================

What this test does
-------------------
Uses Hypothesis to check algebraic / invariant properties of the small
helper functions in ``dispaset.preprocessing.utils``:

* ``_mylogspace``       – log-spaced range: length, monotonicity, endpoints.
* ``_split_list``       – round-trip: every element appears in the result,
                          elements are separated by `` - ``.
* ``_merge_two_dicts``  – commutativity of the union (last-wins), identity,
                          all keys present.
* ``select_units``      – units outside *config['zones']* are always removed;
                          units inside *config['zones']* with non-zero capacity
                          are always kept.

How to run
----------

1. ``pytest tests/unit/test_utils_property.py``
2. ``python tests/unit/test_utils_property.py``

This test does NOT require GAMS.
Hypothesis is used for property-based generation.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import given, assume, settings
from hypothesis import strategies as st

from dispaset.preprocessing.utils import (  # noqa: E402
    _mylogspace,
    _split_list,
    _merge_two_dicts,
    select_units,
)
from dispaset.common import commons  # noqa: E402


# ---------------------------------------------------------------------------
# _mylogspace
# ---------------------------------------------------------------------------

@given(
    low=st.floats(min_value=0.0, max_value=1e4, allow_nan=False, allow_infinity=False),
    high=st.floats(min_value=0.0, max_value=1e4, allow_nan=False, allow_infinity=False),
    N=st.integers(min_value=2, max_value=200),
)
@settings(max_examples=200)
def test_mylogspace_length(low, high, N):
    """The output always has exactly N elements."""
    assume(high >= low)
    result = _mylogspace(low, high, N)
    assert len(result) == N


@given(
    low=st.floats(min_value=0.0, max_value=1e4, allow_nan=False, allow_infinity=False),
    high=st.floats(min_value=0.0, max_value=1e4, allow_nan=False, allow_infinity=False),
    N=st.integers(min_value=2, max_value=200),
)
@settings(max_examples=200)
def test_mylogspace_monotone(low, high, N):
    """The output must be non-decreasing."""
    assume(high >= low)
    result = _mylogspace(low, high, N)
    diffs = np.diff(result)
    assert (diffs >= -1e-9).all(), f"Not monotone: {result}"


@given(
    low=st.floats(min_value=0.0, max_value=1e4, allow_nan=False, allow_infinity=False),
    high=st.floats(min_value=0.0, max_value=1e4, allow_nan=False, allow_infinity=False),
    N=st.integers(min_value=2, max_value=200),
)
@settings(max_examples=200)
def test_mylogspace_endpoints(low, high, N):
    """Last element ≈ high; first element ≈ 0 when low=0 (actual formula: space[0] = -low)."""
    assume(high >= low)
    result = _mylogspace(low, high, N)
    # The formula is np.logspace(0, log10(high+low+1), N) - (low+1),
    # so result[-1] == high but result[0] == -low (not low).
    assert abs(result[-1] - high) < 1e-6, f"Last element {result[-1]} != high {high}"
    assert abs(result[0] - (-low)) < 1e-6, f"First element {result[0]} != -low {-low}"


# ---------------------------------------------------------------------------
# _split_list
# ---------------------------------------------------------------------------

_text = st.text(
    alphabet=st.characters(whitelist_categories=("Lu", "Ll", "Nd"), min_codepoint=32, max_codepoint=127),
    min_size=1, max_size=20,
)


@given(elements=st.lists(_text, min_size=1, max_size=10))
def test_split_list_all_elements_present(elements):
    """Every non-empty element of the list must appear in the result string."""
    result = _split_list(elements)
    for el in elements:
        if str(el) not in ("nan", ""):
            assert str(el) in result


@given(elements=st.lists(_text, min_size=2, max_size=10))
def test_split_list_separator(elements):
    """Elements are joined by \" - \" so the result contains at least one \" - \"."""
    # Filter out empty/nan entries the way _split_list does
    filtered = [e for e in elements if str(e) not in ("nan", "") and str(e).strip() != ""]
    assume(len(filtered) >= 2)
    # _split_list detects the last element by value equality (l != newlist[-1]),
    # so duplicate-valued elements can be treated as the last prematurely.
    # Require all filtered elements are distinct to avoid this edge case.
    assume(len(set(filtered)) == len(filtered))
    result = _split_list(elements)
    assert " - " in result


@given(single=_text)
def test_split_list_single_element_no_separator(single):
    """A one-element list must not contain the \" - \" separator."""
    assume(str(single) not in ("nan", "") and str(single).strip() != "")
    result = _split_list([single])
    assert " - " not in result
    assert result == str(single)


# ---------------------------------------------------------------------------
# _merge_two_dicts
# ---------------------------------------------------------------------------

_key_strategy = st.text(min_size=1, max_size=5, alphabet="abcde")
_val_strategy = st.integers(min_value=0, max_value=100)
_dict_strategy = st.dictionaries(_key_strategy, _val_strategy, max_size=6)


@given(x=_dict_strategy, y=_dict_strategy)
def test_merge_two_dicts_all_keys(x, y):
    """All keys from both dicts must appear in the merged result."""
    merged = _merge_two_dicts(x, y)
    for k in x:
        assert k in merged
    for k in y:
        assert k in merged


@given(x=_dict_strategy, y=_dict_strategy)
def test_merge_two_dicts_y_wins(x, y):
    """For keys in both dicts, y's value wins (last-write wins)."""
    merged = _merge_two_dicts(x, y)
    for k in y:
        assert merged[k] == y[k]


@given(x=_dict_strategy)
def test_merge_two_dicts_identity(x):
    """Merging with an empty dict returns an equal dict."""
    assert _merge_two_dicts(x, {}) == x
    assert _merge_two_dicts({}, x) == x


@given(x=_dict_strategy, y=_dict_strategy)
def test_merge_two_dicts_does_not_mutate(x, y):
    """The original dicts must not be modified."""
    x_copy = dict(x)
    y_copy = dict(y)
    _merge_two_dicts(x, y)
    assert x == x_copy
    assert y == y_copy


# ---------------------------------------------------------------------------
# select_units
# ---------------------------------------------------------------------------

# Known valid zone name for tests (first in commons zones list from config)
_VALID_ZONE = "Z1"
_TECH = "GTUR"
_FUEL = "GAS"


def _make_units(zones: list[str], capacities: list[float]) -> pd.DataFrame:
    """Build a minimal units DataFrame for select_units."""
    rows = []
    for i, (z, cap) in enumerate(zip(zones, capacities)):
        rows.append({
            "Unit": f"U{i}",
            "Zone": z,
            "Technology": _TECH,
            "Fuel": _FUEL,
            "PowerCapacity": cap,
            "STOMaxChargingPower": 0.0,
            "Nunits": 1,
        })
    df = pd.DataFrame(rows)
    df.index = range(len(df))
    return df


@given(
    extra_zones=st.lists(
        st.text(min_size=2, max_size=4, alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ"),
        min_size=1, max_size=4,
    ),
    n_bad=st.integers(min_value=1, max_value=4),
)
def test_select_units_removes_unknown_zones(extra_zones, n_bad):
    """Units whose zone is not in config['zones'] must be removed."""
    assume(all(z != _VALID_ZONE for z in extra_zones))
    bad_zones = extra_zones[:n_bad]
    units = _make_units(
        [_VALID_ZONE] + bad_zones,
        [100.0] * (1 + len(bad_zones)),
    )
    config = {"zones": [_VALID_ZONE]}
    result = select_units(units.copy(), config)
    assert (result["Zone"] == _VALID_ZONE).all()
    assert len(result) == 1


@given(
    n_good=st.integers(min_value=1, max_value=5),
    capacity=st.floats(min_value=1.0, max_value=1000.0, allow_nan=False, allow_infinity=False),
)
def test_select_units_keeps_nonzero_capacity_in_zone(n_good, capacity):
    """Units with non-zero capacity in a known zone are always kept."""
    units = _make_units([_VALID_ZONE] * n_good, [capacity] * n_good)
    config = {"zones": [_VALID_ZONE]}
    result = select_units(units.copy(), config)
    assert len(result) == n_good


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
