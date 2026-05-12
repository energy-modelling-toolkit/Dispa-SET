# -*- coding: utf-8 -*-
"""
Unit tests: solver KPI helpers (get_solver_kpis, assert_feasible)
=================================================================

What this test does
-------------------
Tests the two KPI helper functions added to ``tests/_helpers.py``:

* ``get_solver_kpis(results)`` — returns ``total_cost``,
  ``total_lostload_mwh`` and ``storage_violation_mwh`` from a results dict.
* ``assert_feasible(results, max_lostload_mwh)`` — raises ``AssertionError``
  when lost load exceeds the threshold or when the cost is non-finite.

All tests use synthetic (mock) results dicts and require NO GAMS installation.

How to run
----------

1. ``pytest tests/unit/test_kpi_helper.py``
2. ``python tests/unit/test_kpi_helper.py``
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd
import pytest

# Make both the repo root and the tests/ directory importable.
if __package__ is None or __package__ == "":
    _tests_dir = Path(__file__).resolve().parents[1]
    _repo_root = _tests_dir.parent
    sys.path.insert(0, str(_repo_root))
    sys.path.insert(0, str(_tests_dir))

from _helpers import get_solver_kpis, assert_feasible  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_scalar_results(**kwargs) -> dict:
    """Build a minimal results dict using scalar DataFrames."""
    idx = pd.date_range("2015-01-01", periods=3, freq="h")
    results = {}
    for key, value in kwargs.items():
        results[key] = pd.DataFrame({"Z1": [value] * 3}, index=idx)
    return results


# ---------------------------------------------------------------------------
# get_solver_kpis
# ---------------------------------------------------------------------------

class TestGetSolverKpis:

    def test_returns_total_cost(self):
        r = _make_scalar_results(OutputSystemCost=100.0)
        kpis = get_solver_kpis(r)
        assert "total_cost" in kpis
        assert abs(kpis["total_cost"] - 300.0) < 1e-6  # 100 * 3 hours

    def test_cost_absent_key_omitted(self):
        r = {}
        kpis = get_solver_kpis(r)
        assert "total_cost" not in kpis

    def test_zero_lostload_when_no_keys(self):
        r = _make_scalar_results(OutputSystemCost=50.0)
        kpis = get_solver_kpis(r)
        assert kpis["total_lostload_mwh"] == 0.0

    def test_lostload_maxpower_counted(self):
        r = _make_scalar_results(LostLoad_MaxPower=2.0)
        kpis = get_solver_kpis(r)
        assert kpis["total_lostload_mwh"] == pytest.approx(6.0)  # 2 * 3 hours

    def test_multiple_lostload_keys_summed(self):
        r = _make_scalar_results(
            LostLoad_MaxPower=1.0,
            LostLoad_aFRRU=1.0,
            LostLoad_RampUp=0.5,
        )
        kpis = get_solver_kpis(r)
        # Each key contributes 3 h × value
        assert kpis["total_lostload_mwh"] == pytest.approx(3.0 + 3.0 + 1.5)

    def test_storage_violation_separate(self):
        r = _make_scalar_results(
            LostLoad_MaxPower=0.0,
            LostLoad_StorageLevelViolation=5.0,
        )
        kpis = get_solver_kpis(r)
        assert kpis["total_lostload_mwh"] == 0.0
        assert kpis["storage_violation_mwh"] == pytest.approx(15.0)  # 5 * 3 hours

    def test_negative_lostload_clamped_to_zero(self):
        """Numerical noise can produce tiny negative values; they are clamped."""
        r = _make_scalar_results(LostLoad_MaxPower=-1e-9)
        kpis = get_solver_kpis(r)
        assert kpis["total_lostload_mwh"] == 0.0


# ---------------------------------------------------------------------------
# assert_feasible
# ---------------------------------------------------------------------------

class TestAssertFeasible:

    def test_no_lostload_passes(self):
        r = _make_scalar_results(OutputSystemCost=50.0)
        kpis = assert_feasible(r)
        assert isinstance(kpis, dict)

    def test_small_lostload_within_tolerance_passes(self):
        r = _make_scalar_results(LostLoad_MaxPower=0.1)  # 0.3 MWh total
        assert_feasible(r, max_lostload_mwh=1.0)

    def test_large_lostload_raises(self):
        r = _make_scalar_results(LostLoad_MaxPower=5.0)  # 15 MWh total
        with pytest.raises(AssertionError, match="Lost load"):
            assert_feasible(r, max_lostload_mwh=1.0)

    def test_zero_threshold_strict(self):
        r = _make_scalar_results(LostLoad_MinPower=0.01)  # 0.03 MWh
        with pytest.raises(AssertionError):
            assert_feasible(r, max_lostload_mwh=0.0)

    def test_negative_cost_raises(self):
        r = _make_scalar_results(OutputSystemCost=-10.0)
        with pytest.raises(AssertionError, match="Negative"):
            assert_feasible(r, max_lostload_mwh=0.0)

    def test_non_finite_cost_raises(self):
        idx = pd.date_range("2015-01-01", periods=2, freq="h")
        r = {"OutputSystemCost": pd.DataFrame({"Z1": [math.inf, 1.0]}, index=idx)}
        with pytest.raises(AssertionError, match="Non-finite"):
            assert_feasible(r, max_lostload_mwh=0.0)

    def test_returns_kpis_dict(self):
        r = _make_scalar_results(OutputSystemCost=200.0)
        kpis = assert_feasible(r)
        assert "total_cost" in kpis
        assert "total_lostload_mwh" in kpis
        assert "storage_violation_mwh" in kpis


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
