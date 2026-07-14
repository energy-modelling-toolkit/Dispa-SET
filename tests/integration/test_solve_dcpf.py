# -*- coding: utf-8 -*-
"""
Integration test: DC Power Flow vs NTC transmission
=====================================================

What this test does
-------------------
Tests both DC power flow and NTC transmission on the same 3-zone
triangle topology (Z1 → Z2, ZZ3 → Z2, Z1 → ZZ3).

Zone roles:
  * Z1  – sole generator zone (500 MW gas unit)
  * Z2  – pure-load zone (~170 MW at modifier 0.02)
  * ZZ3 – pure-load zone (~ 17 MW at modifier 0.02)

**Grid topology** (``tests/data/GridData_tiny.csv``)::

    Z1 ──(X=0.10)──► Z2
    │                 ▲
    (X=0.40)    (X=0.10)
    ▼                 │
    ZZ3 ──────────────┘

Asymmetric reactances create loop flows under DC-PF:
  With NTC, the optimizer sends power directly (Z1→Z2 ≈ 170 MW,
  Z1→ZZ3 ≈ 17 MW, ZZ3→Z2 ≈ 0 MW).
  With DC-PF, physics forces loop flows (ZZ3→Z2 ≈ +17 MW even
  though ZZ3 is a pure load zone).

Tests
-----
1. **test_solve_dcpf_3zones** – builds and solves with
   ``DC-Power Flow``; verifies solve succeeds and OutputFlow exists.
2. **test_ntc_vs_dcpf_loop_flows** – runs both NTC and DC-PF on the
   identical system and asserts that the ZZ3→Z2 line carries a
   significantly larger mean flow under DC-PF than NTC (demonstrating
   the loop-flow effect).

Performance budget
------------------
Both GAMS solves are 1-day horizons with one unit → typically < 10 s each.
Pytest timeout: 120 s.

How to run
----------
1. ``pytest tests/integration/test_solve_dcpf.py``
2. ``python tests/integration/test_solve_dcpf.py``

Skipped automatically when no GAMS installation is detected.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import build_solve, load_test_config, skip_if_no_gams  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mean_flow(results, line_name: str) -> float:
    """Return the mean hourly flow on *line_name* (0 if line absent)."""
    flows = results.get("OutputFlow", None)
    if flows is None or flows.empty:
        return 0.0
    if line_name in flows.columns:
        return float(flows[line_name].mean())
    return 0.0


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.timeout(120)
def test_solve_dcpf_3zones():
    """DC-Power-Flow 3-zone build+solve succeeds and the PTDF parameter is
    non-zero in the GDX inputs."""
    skip_if_no_gams()

    cfg = load_test_config("tiny_dcpf.yml", "dcpf_solve")
    out = build_solve(cfg)

    inputs  = out["inputs"]
    results = out["results"]

    # --- Basic solve checks ---
    assert "OutputPower" in results, "OutputPower missing from results"
    assert "OutputFlow"  in results, "OutputFlow missing from results (should exist for multi-zone)"
    assert not results["OutputFlow"].empty, "OutputFlow DataFrame is unexpectedly empty"

    # --- PTDF is non-zero (DC-PF was actually applied) ---
    ptdf_val = inputs["parameters"]["PTDF"]["val"]
    assert np.any(ptdf_val != 0), (
        "PTDF parameter is all-zero in DC-Power-Flow result — "
        "the physics-based line constraints were not applied"
    )

    # --- All three zones appear in the result ---
    zones_in_result = inputs["sets"]["n"]
    for z in ("Z1", "Z2", "ZZ3"):
        assert z in zones_in_result, f"Zone {z} missing from simulation sets"

    # --- Flows are present on expected lines ---
    flow_cols = list(results["OutputFlow"].columns)
    assert any("Z1" in c for c in flow_cols), "No Z1 line found in OutputFlow"

    print(f"DC-PF build {out['build_time']:.2f}s  solve {out['solve_time']:.2f}s")
    print(f"PTDF non-zero entries: {np.count_nonzero(ptdf_val)}")
    print(f"OutputFlow columns: {flow_cols}")


@pytest.mark.timeout(240)
def test_ntc_vs_dcpf_loop_flows():
    """DC-Power-Flow shows loop flows absent from the NTC solution.

    The key observable: under DC-PF the ZZ3→Z2 line carries a substantial
    positive mean flow (~17 MW), while under NTC the optimizer does not use
    that line (mean flow ≈ 0 MW).

    Both simulations use the identical load data, units, and demand
    modifier so the total generation is the same; only the transmission
    model differs.
    """
    skip_if_no_gams()

    # Run NTC
    cfg_ntc  = load_test_config("tiny_ntc_3zones.yml", "ntc3_compare")
    out_ntc  = build_solve(cfg_ntc)

    # Run DC-PF
    cfg_dcpf = load_test_config("tiny_dcpf.yml", "dcpf_compare")
    out_dcpf = build_solve(cfg_dcpf)

    results_ntc  = out_ntc["results"]
    results_dcpf = out_dcpf["results"]

    # --- Both simulations solved ---
    assert "OutputPower" in results_ntc,  "NTC solve failed (no OutputPower)"
    assert "OutputPower" in results_dcpf, "DC-PF solve failed (no OutputPower)"

    # --- PTDF is non-zero for DC-PF and blank/NaN for NTC ---
    ptdf_dcpf = np.array(out_dcpf["inputs"]["parameters"]["PTDF"]["val"], dtype=float)
    ptdf_ntc  = np.array(out_ntc["inputs"]["parameters"]["PTDF"]["val"],  dtype=float)
    assert np.any(ptdf_dcpf != 0), "DC-PF PTDF should be non-zero"
    # NTC PTDF is all-NaN (blank) – no finite non-zero entries:
    assert not np.any(np.isfinite(ptdf_ntc) & (ptdf_ntc != 0)), (
        "NTC PTDF should have no finite non-zero entries"
    )

    # --- Loop flow on ZZ3→Z2 ---
    # Under DC-PF: physics forces ~+17 MW loop flow on ZZ3→Z2.
    # Under NTC:   optimizer routes directly; ZZ3→Z2 is unused (~0 MW).
    flow_dcpf_zzz_z2 = _mean_flow(results_dcpf, "ZZ3 -> Z2")
    flow_ntc_zzz_z2  = _mean_flow(results_ntc,  "ZZ3 -> Z2")

    print(f"NTC   ZZ3→Z2 mean flow: {flow_ntc_zzz_z2:+.1f} MW")
    print(f"DC-PF ZZ3→Z2 mean flow: {flow_dcpf_zzz_z2:+.1f} MW")

    # DC-PF must produce a significantly larger loop flow on ZZ3→Z2.
    # Expected: DC-PF ≈ +17 MW, NTC ≈ 0 MW.
    assert flow_dcpf_zzz_z2 > flow_ntc_zzz_z2 + 5, (
        f"DC-PF loop flow on ZZ3→Z2 ({flow_dcpf_zzz_z2:.1f} MW) should be "
        f"at least 5 MW larger than NTC ({flow_ntc_zzz_z2:.1f} MW)"
    )

    # --- Z1→Z2 flows differ between NTC and DC-PF ---
    # Under DC-PF: ~153 MW (physics splits power via ZZ3 path)
    # Under NTC:   ~170 MW (direct routing from Z1 to Z2)
    flow_dcpf_z1_z2 = _mean_flow(results_dcpf, "Z1 -> Z2")
    flow_ntc_z1_z2  = _mean_flow(results_ntc,  "Z1 -> Z2")

    print(f"NTC   Z1→Z2 mean flow: {flow_ntc_z1_z2:+.1f} MW")
    print(f"DC-PF Z1→Z2 mean flow: {flow_dcpf_z1_z2:+.1f} MW")

    # Under NTC the direct Z1→Z2 flow should exceed the DC-PF value
    # (because DC-PF routes some power through ZZ3 instead).
    assert flow_ntc_z1_z2 > flow_dcpf_z1_z2 + 5, (
        f"NTC Z1→Z2 ({flow_ntc_z1_z2:.1f} MW) should exceed DC-PF "
        f"({flow_dcpf_z1_z2:.1f} MW) by at least 5 MW"
    )

    print(f"NTC  build {out_ntc['build_time']:.2f}s  solve {out_ntc['solve_time']:.2f}s")
    print(f"DCPF build {out_dcpf['build_time']:.2f}s  solve {out_dcpf['solve_time']:.2f}s")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
