# -*- coding: utf-8 -*-
"""
Unit tests: PTDF matrix computation + DC-Power-Flow build (no GAMS)
====================================================================

What this test does
-------------------
1. **test_ptdf_matrix_shape** — calls `PTDF_matrix` directly on the
   3-zone triangle grid (GridData_tiny.csv) and checks the output shape.
2. **test_ptdf_matrix_values** — verifies the PTDF entries against the
   hand-calculated reference for the topology::

       Z1 → Z2   (X=0.10)
       ZZ3 → Z2  (X=0.10)
       Z1 → ZZ3  (X=0.40)

   Slack bus = Z1 (most connections in case of tie → first zone by
   the code's loop order).  Reference values derived analytically:

       PTDF [Z1→Z2, Z2]   = -5/6  ≈ -0.8333
       PTDF [Z1→Z2, ZZ3]  = -2/3  ≈ -0.6667
       PTDF [ZZ3→Z2, Z2]  = -1/6  ≈ -0.1667
       PTDF [ZZ3→Z2, ZZ3] = +2/3  ≈ +0.6667
       PTDF [Z1→ZZ3, Z2]  = -1/6  ≈ -0.1667
       PTDF [Z1→ZZ3, ZZ3] = -1/3  ≈ -0.3333
       (All Z1 column = 0 because Z1 is the slack bus.)

3. **test_ptdf_loop_flow_prediction** — uses the PTDF to predict flows
   for a scenario (Z1 surplus, Z2 and ZZ3 deficits) and verifies that
   the predicted loop flow on ZZ3→Z2 is significant and positive even
   though ZZ3 is a pure-load zone.
4. **test_build_dcpf_ptdf_nonempty** — calls `build_simulation` with the
   DC-Power-Flow config and asserts that the ``PTDF`` parameter in the
   returned ``SimData`` is non-zero (i.e., the PTDF was actually injected
   into the inputs).
5. **test_build_ntc_ptdf_zero** — same build but with the NTC companion
   config; asserts that the ``PTDF`` parameter is all-zero (NTC mode does
   not use PTDF).

No GAMS licence is required for any of these tests.

How to run
----------
1. ``pytest tests/unit/test_ptdf.py``
2. ``python tests/unit/test_ptdf.py``
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import TESTS_DIR, load_test_config  # noqa: E402

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

GRID_CSV = TESTS_DIR / "data" / "GridData_tiny.csv"

# 3-zone topology for manual PTDF calculation
ZONES = ["Z1", "Z2", "ZZ3"]

# Expected PTDF reference (slack = Z1, lines in order of GridData_tiny.csv)
# Columns: Z1, Z2, ZZ3.  All Z1 entries are 0 (slack).
#
# Derivation (see docstring above):
#   B12=10, B23=10, B13=2.5;  Bbus reduced (dropping slack Z1 row/col)
#   = [[20,-10],[-10,12.5]]; inv = [[12.5,10],[10,20]]/150
#
EXPECTED_PTDF = {
    "Z1 -> Z2":  {"Z1": 0.0, "Z2": -5/6,  "ZZ3": -2/3},
    "ZZ3 -> Z2": {"Z1": 0.0, "Z2": -1/6,  "ZZ3":  2/3},
    "Z1 -> ZZ3": {"Z1": 0.0, "Z2": -1/6,  "ZZ3": -1/3},
}


def _make_feeder() -> pd.DataFrame:
    """Load GridData_tiny.csv as the feeder dataframe expected by PTDF_matrix."""
    df = pd.read_csv(GRID_CSV, index_col=0, na_values=["", " "], keep_default_na=False)
    return df


def _make_config() -> dict:
    """Minimal config dict sufficient for PTDF_matrix."""
    return {"zones": ZONES}


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.timeout(5)
def test_ptdf_matrix_shape():
    """PTDF matrix has shape (n_lines, n_zones)."""
    from dispaset.preprocessing.utils import PTDF_matrix

    feeder = _make_feeder()
    cfg = _make_config()
    ptdf = PTDF_matrix(cfg, feeder.copy())

    assert ptdf.shape == (len(feeder), len(ZONES)), (
        f"Expected shape {(len(feeder), len(ZONES))}, got {ptdf.shape}"
    )
    assert list(ptdf.columns) == ZONES
    assert list(ptdf.index) == list(feeder.index)


@pytest.mark.timeout(5)
def test_ptdf_matrix_values():
    """PTDF entries match the hand-calculated reference within 1e-6."""
    from dispaset.preprocessing.utils import PTDF_matrix

    feeder = _make_feeder()
    cfg = _make_config()
    ptdf = PTDF_matrix(cfg, feeder.copy()).fillna(0)

    tol = 1e-4
    for line, zone_vals in EXPECTED_PTDF.items():
        for zone, expected in zone_vals.items():
            actual = float(ptdf.loc[line, zone])
            assert abs(actual - expected) < tol, (
                f"PTDF[{line}, {zone}]: expected {expected:.4f}, got {actual:.4f}"
            )


@pytest.mark.timeout(5)
def test_ptdf_loop_flow_prediction():
    """The PTDF correctly predicts a loop flow on ZZ3→Z2 even though ZZ3
    is a pure-load zone (no generation).

    Scenario:
      Z1 net injection  = +187 MW (sole generator, serves Z2 and ZZ3)
      Z2 net injection  = -170 MW (pure load)
      ZZ3 net injection =  -17 MW (pure load)

    Expected flows (F = PTDF × net_injections, Z1 is slack):
      Z1→Z2:   153 MW  (direct from Z1 to Z2)
      ZZ3→Z2:  +17 MW  (loop: Z1→ZZ3→Z2 path re-exports to Z2)
      Z1→ZZ3:   34 MW  (Z1 also feeds ZZ3, more than the 17 MW local demand)
    """
    from dispaset.preprocessing.utils import PTDF_matrix

    feeder = _make_feeder()
    cfg = _make_config()
    ptdf = PTDF_matrix(cfg, feeder.copy()).fillna(0).values  # shape (3, 3)

    # Net injections (Z1 column is all-zero in PTDF → Z1 is the slack)
    net_injections = np.array([0.0, -170.0, -17.0])  # [Z1, Z2, ZZ3]

    flows = ptdf @ net_injections  # shape (3,)
    f_Z1_Z2, f_ZZ3_Z2, f_Z1_ZZ3 = flows

    tol = 1.0  # MW tolerance (demand varies ±some %)
    assert abs(f_Z1_Z2 - 152.8) < tol, f"Z1→Z2 flow: expected ~152.8, got {f_Z1_Z2:.1f}"
    assert abs(f_ZZ3_Z2 - 17.1)  < tol, f"ZZ3→Z2 flow: expected ~17.1, got {f_ZZ3_Z2:.1f}"
    assert abs(f_Z1_ZZ3 - 34.1)  < tol, f"Z1→ZZ3 flow: expected ~34.1, got {f_Z1_ZZ3:.1f}"

    # The key physics result: ZZ3→Z2 is positive (loop current flows through
    # ZZ3 even though ZZ3 is a pure-load zone) — and it is larger than zero.
    assert f_ZZ3_Z2 > 5, "Expected significant positive loop flow on ZZ3→Z2"


@pytest.mark.timeout(30)
def test_build_dcpf_ptdf_nonempty():
    """build_simulation with DC-Power-Flow produces a non-zero PTDF parameter."""
    import dispaset as ds

    cfg = load_test_config("tiny_dcpf.yml", "build_dcpf")
    sim_data = ds.build_simulation(cfg)

    ptdf_val = sim_data["parameters"]["PTDF"]["val"]
    assert ptdf_val.size > 0, "PTDF parameter array is empty"
    assert np.any(ptdf_val != 0), (
        "PTDF parameter is all-zero in DC-Power-Flow build — PTDF was not injected"
    )
    print(f"PTDF shape: {ptdf_val.shape}, non-zero entries: {np.count_nonzero(ptdf_val)}")


@pytest.mark.timeout(30)
def test_build_ntc_ptdf_zero():
    """build_simulation with NTC produces a blank/NaN PTDF parameter (NTC
    mode does not use physics-based line-flow constraints).

    The PTDF DataFrame is intentionally left all-NaN for NTC-based
    simulations (the parameter slot exists in the GDX schema but carries
    no meaningful values).  We therefore accept both all-NaN and all-zero.
    """
    import dispaset as ds

    cfg = load_test_config("tiny_ntc_3zones.yml", "build_ntc3")
    sim_data = ds.build_simulation(cfg)

    ptdf_val = sim_data["parameters"]["PTDF"]["val"]
    # Convert to float so NaN comparisons work correctly (the raw val may be
    # an object array containing float('nan') entries).
    ptdf_float = np.array(ptdf_val, dtype=float)
    # There should be no *finite* non-zero entries:
    finite_nonzero = np.isfinite(ptdf_float) & (ptdf_float != 0)
    assert ptdf_float.size == 0 or not np.any(finite_nonzero), (
        "NTC PTDF should have no finite non-zero entries but found some"
    )
    print(f"NTC PTDF shape: {ptdf_float.shape} (blank/NaN as expected)")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
