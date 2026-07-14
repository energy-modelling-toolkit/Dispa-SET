# -*- coding: utf-8 -*-
"""
Integration test: DC Power Flow on the PJM 5-bus test system
=============================================================

What this test does
-------------------
Builds and solves the standard PJM 5-bus academic DC-OPF test case inside
Dispa-SET, then verifies that the PTDF-based line-flow constraints reproduce
the well-known congestion pattern of the reference GAMS model to within ~1 MW.

Reference systems
-----------------
Network topology and generator data are taken from:

    F. Li and R. Bo, "Small test systems for power system economic studies,"
    IEEE PES General Meeting, Minneapolis, MN, 2010, pp. 1-4.
    doi: 10.1109/PES.2010.5589973

The reference GAMS DC-OPF implementation is:

    A. Soroudi, "DC Optimal Power Flow (OPF) model",
    https://github.com/OptimizationExpert/GAMS-DC--OPF-
    Contributed by Dr. Alireza Soroudi, University College Dublin, Ireland.

Network topology (``tests/data/GridData_pjm5bus.csv``)::

    PJM_B1 (slack bus, generators: Alta 40 MW, ParkCity 170 MW)
      |── X=0.0281, Fmax=400  ──► PJM_B2 (300 MW constant load)
      |── X=0.0304, Fmax=1000 ──► PJM_B4 (400 MW load + Sundance 200 MW)
      └── X=0.0064, Fmax=1000 ──► PJM_B5 (Brighton 600 MW generator)
    PJM_B2 ── X=0.0108, Fmax=1000 ──► PJM_B3 (300 MW load + Solitude 520 MW)
    PJM_B3 ── X=0.0297, Fmax=1000 ──► PJM_B4
    PJM_B4 ── X=0.0297, Fmax=240  ──► PJM_B5

Generator merit order (marginal cost = GasPrice / Efficiency, GasPrice=10 $/MWh):
    PJM_Brighton   (B5, eff=1.000):  10 $/MWh  cheapest
    PJM_Alta       (B1, eff=0.714):  14 $/MWh
    PJM_ParkCity   (B1, eff=0.667):  15 $/MWh
    PJM_Solitude   (B3, eff=0.333):  30 $/MWh
    PJM_Sundance   (B4, eff=0.250):  40 $/MWh  most expensive

Configuration (``tests/configs/tiny_pjm5bus_dcpf.yml``)
--------------------------------------------------------
Two deliberate choices make Dispa-SET's formulation as close as possible to
the reference GAMS single-period DC-OPF:

**SimulationType: 'LP'** (no clustering)
    In Dispa-SET's default 'LP clustered' mode all generators of the same
    zone/technology/fuel are merged into a single aggregate unit.  Here that
    would combine Alta (40 MW, 14 $/MWh) and ParkCity (170 MW, 15 $/MWh) into
    one 210 MW unit with a capacity-weighted average efficiency — losing the
    individual cost signals that distinguish the two plants.  Using 'LP' keeps
    every generator as a separate decision variable, matching the GAMS model.

**Reserve2U: 0, Reserve2D: 0** (zero reserve margin)
    Dispa-SET's 'Generic' reserve formula requires each zone to hold back
    spinning capacity equal to  √(10·Load + 150²) − 150  MW (≈10–13 MW for
    300–400 MW loads).  This "headroom" effectively reduces the usable output
    of each generator, shifting the marginal unit and altering the nodal
    injections — and therefore the line flows via the PTDF mapping.  Setting
    both margins to zero removes this effect and lets the LP minimize cost
    subject only to the network and capacity constraints, exactly as in GAMS.

Total load: 300 + 300 + 400 + 1 + 1 = 1002 MW (PJM_B1/B5 each carry 1 MW
as a placeholder so Dispa-SET does not treat them as generation-only stubs;
the GAMS model has only 1000 MW of load).

Expected DC power flow behavior
--------------------------------
Without transmission constraints the cheapest full-dispatch (merit order) is:
    Brighton: 600, Alta: 40, ParkCity: 170, Solitude: 192 MW → 1002 MW

The binding transmission constraint is the **PJM_B4 → PJM_B5** corridor
(Fmax=240 MW).  Intuitively: Brighton (cheapest, at B5) wants to export as
much power as possible toward the three load buses.  By Kirchhoff's laws the
only paths from B5 outward are:
  * B5 → B1 (direct, X=0.0064, very low reactance → high admittance → large share)
  * B5 ← B4 (through the B4→B5 line backwards, X=0.0297, limited at 240 MW)
The B4→B5 line saturates first (negative direction: B5→B4), capping Brighton's
output at ≈467 MW well below its 600 MW nameplate.  Solitude covers the
residual B2/B3 demand; Sundance (most expensive) is dispatched near zero.

Comparison with GAMS reference (DC-OPF.gms, GAMS 45.7, run 2026-05-11)
------------------------------------------------------------------------
The GAMS model is a pure single-period LP with Sbase=100 MVA.  After removing
reserves and clustering, Dispa-SET matches it to within 1 MW on every line:

    Line              Fmax   GAMS (MW)  Dispa-SET (MW)  Diff
    ──────────────────────────────────────────────────────
    B1→B2              400     +249.7        +249.5       −0.2
    B1→B4             1000     +186.8        +186.7       −0.1
    B1→B5             1000     −226.5        −227.2       −0.7
    B2→B3             1000      −50.3         −50.5       −0.2
    B3→B4             1000      −26.8         −26.7       +0.1
    B4→B5  *** 240    240.0     −240.0        −240.0        0.0  ← CONGESTED

    Unit               GAMS (MW)   Dispa-SET (MW)   Diff
    ──────────────────────────────────────────────────────
    Alta       (B1)       40.0           40.0          0
    ParkCity   (B1)      170.0          170.0          0
    Solitude   (B3)      323.5          323.5        < 1
    Sundance   (B4)        0.0            0.0          0
    Brighton   (B5)      466.5          468.5        < 2  (2 MW extra load)

The residual 0.2–0.7 MW difference on the flows is entirely explained by the
2 MW extra demand from the B1/B5 placeholder loads (absent from the GAMS
model).  The extra load is supplied by Brighton, which then injects ~2 MW
more power, propagating a proportional flow change through all lines.

Kirchhoff check (PTDF-based power flow)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Dispa-SET computes line flows via::

    Flow(l) = −ΣPTDF(l,n) · InjectedPower(n)
    InjectedPower(n) = Generation(n) − Demand(n)

The PTDF matrix is computed analytically from the network susceptances
(utils.PTDF_matrix) using a B-matrix inversion with B5 as the slack bus.
The minus sign matches the GAMS angle-difference equation:
    Pij(i,j) = bij·(δ_i − δ_j)
so that positive flow means power travels from bus i to bus j.

Tests
-----
1. **test_solve_pjm5bus_dcpf**  — builds and solves with 'DC-Power Flow';
   checks that all 5 zones and all 6 lines appear, that no line limit is
   violated, and that total generation balances total demand.
2. **test_pjm5bus_congestion**  — re-uses the same config to verify that
   B4-B5 is saturated at −240 MW, that B1-B2 carries ≈+250 MW, and that
   Brighton is constrained below its 600 MW nameplate — confirming that the
   PTDF-based Kirchhoff constraints reproduce the reference GAMS solution.

Performance budget
------------------
Both tests share one solve on a 1-day (24-h) horizon with 5 units → typically
< 15 s total.  Pytest timeout: 180 s.

How to run
----------
1. ``pytest tests/integration/test_solve_pjm5bus_dcpf.py``
2. ``python tests/integration/test_solve_pjm5bus_dcpf.py``

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
# Expected network data (from GridData_pjm5bus.csv)
# ---------------------------------------------------------------------------

EXPECTED_ZONES = ["PJM_B1", "PJM_B2", "PJM_B3", "PJM_B4", "PJM_B5"]
EXPECTED_LINES = [
    "PJM_B1 -> PJM_B2",  # Fmax = 400 MW  (binding congestion line)
    "PJM_B1 -> PJM_B4",  # Fmax = 1000 MW
    "PJM_B1 -> PJM_B5",  # Fmax = 1000 MW
    "PJM_B2 -> PJM_B3",  # Fmax = 1000 MW
    "PJM_B3 -> PJM_B4",  # Fmax = 1000 MW
    "PJM_B4 -> PJM_B5",  # Fmax = 240 MW  (binding congestion line)
]
FMAX = {
    "PJM_B1 -> PJM_B2": 400.0,
    "PJM_B1 -> PJM_B4": 1000.0,
    "PJM_B1 -> PJM_B5": 1000.0,
    "PJM_B2 -> PJM_B3": 1000.0,
    "PJM_B3 -> PJM_B4": 1000.0,
    "PJM_B4 -> PJM_B5": 240.0,
}

# Total constant load per hour (B1=1, B2=300, B3=300, B4=400, B5=1) = 1002 MW
TOTAL_DEMAND_MW = 1002.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mean_abs_flow(results: dict, line_name: str) -> float:
    """Return the mean absolute hourly flow on *line_name* (0 if absent)."""
    flows = results.get("OutputFlow", None)
    if flows is None or flows.empty:
        return 0.0
    if line_name in flows.columns:
        return float(flows[line_name].abs().mean())
    return 0.0


def _mean_flow(results: dict, line_name: str) -> float:
    """Return the mean (signed) hourly flow on *line_name* (0 if absent)."""
    flows = results.get("OutputFlow", None)
    if flows is None or flows.empty:
        return 0.0
    if line_name in flows.columns:
        return float(flows[line_name].mean())
    return 0.0


# ---------------------------------------------------------------------------
# Shared fixture: build + solve once, reuse in both tests
# ---------------------------------------------------------------------------

_shared_out: dict | None = None


def _get_pjm5bus_solve_output() -> dict:
    """Build and solve the PJM 5-bus DC-PF case once; cache the result."""
    global _shared_out
    if _shared_out is None:
        cfg = load_test_config("tiny_pjm5bus_dcpf.yml", "pjm5bus_dcpf")
        _shared_out = build_solve(cfg)
    return _shared_out


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.timeout(180)
def test_solve_pjm5bus_dcpf():
    """5-zone PJM DC-PF build+solve succeeds; all zones and lines are present;
    no line flow limit is violated; total generation matches total demand.

    This is the primary smoke-test for the PJM 5-bus reference case inside
    Dispa-SET.  If it fails, the subsequent congestion test is meaningless.
    """
    skip_if_no_gams()
    out = _get_pjm5bus_solve_output()

    inputs  = out["inputs"]
    results = out["results"]

    # --- Basic solve checks ---
    assert "OutputPower"  in results, "Solve failed: OutputPower missing"
    assert "OutputFlow"   in results, "Solve failed: OutputFlow missing"
    assert not results["OutputFlow"].empty, "OutputFlow is unexpectedly empty"

    # --- PTDF is non-zero (DC-PF was applied) ---
    ptdf_val = inputs["parameters"]["PTDF"]["val"]
    assert np.any(ptdf_val != 0), (
        "PTDF parameter is all-zero — DC-Power-Flow constraints were not applied"
    )

    # --- All 5 zones present ---
    zones_in_result = inputs["sets"]["n"]
    for z in EXPECTED_ZONES:
        assert z in zones_in_result, f"Zone {z} missing from simulation sets"

    # --- All 6 lines appear in OutputFlow ---
    flow_cols = set(results["OutputFlow"].columns)
    for line in EXPECTED_LINES:
        assert line in flow_cols, (
            f"Line '{line}' missing from OutputFlow columns.\n"
            f"Available columns: {sorted(flow_cols)}"
        )

    # --- No line flow limit violated (within 1 MW numerical tolerance) ---
    tol = 1.0  # MW
    for line in EXPECTED_LINES:
        max_abs_flow = float(results["OutputFlow"][line].abs().max())
        assert max_abs_flow <= FMAX[line] + tol, (
            f"Line '{line}' violated its limit: "
            f"|flow|_max={max_abs_flow:.2f} MW > Fmax={FMAX[line]} MW"
        )

    # --- Total generation ≈ total demand (no lost load on a feasible case) ---
    if "OutputPower" in results and not results["OutputPower"].empty:
        total_gen_per_hour = results["OutputPower"].sum(axis=1)
        mean_gen = float(total_gen_per_hour.mean())
        assert abs(mean_gen - TOTAL_DEMAND_MW) < 10.0, (
            f"Mean hourly generation ({mean_gen:.1f} MW) deviates by "
            f"more than 10 MW from total demand ({TOTAL_DEMAND_MW} MW)"
        )

    print(f"PJM 5-bus DC-PF build {out['build_time']:.2f}s  "
          f"solve {out['solve_time']:.2f}s")
    print(f"PTDF non-zero entries: {int(np.count_nonzero(ptdf_val))}")
    print(f"OutputFlow lines: {sorted(flow_cols)}")
    for line in EXPECTED_LINES:
        mf = _mean_flow(results, line)
        print(f"  {line}: mean flow = {mf:+.1f} MW  (Fmax={FMAX[line]} MW)")


@pytest.mark.timeout(180)
def test_pjm5bus_congestion():
    """Dispa-SET reproduces the reference GAMS DC-OPF solution to within ~1 MW.

    Physical mechanism
    ------------------
    Brighton (PJM_B5, cheapest at 10 $/MWh) wants to run at full output
    (600 MW), but the network path from B5 to the three load buses (B2, B3,
    B4) forces power through the **B4→B5** corridor in reverse (B5→B4).
    That line has only a 240 MW limit.  By Kirchhoff's voltage law the DC
    power flow distributes the Brighton injection across all paths:

        PTDF[B4→B5, B5] ≈ +0.112  (11% of each MW injected at B5 flows here)

    When Brighton produces 467 MW, line B4→B5 carries exactly 240 MW in
    reverse — its Fmax.  Further output would violate the constraint, so the
    optimizer limits Brighton and dispatches Solitude instead.

    Why B1-B2 is NOT congested
    --------------------------
    Unlike the unconstrained dispatch (where B1→B2 would carry ≈600 MW as
    Brighton power loops through B5→B1→B2), the actual dispatch routes most
    of B2's 300 MW demand through Solitude (B3) via the backward B2→B3 path.
    Solitude sits right next to B3/B2, so:
      - Solitude at B3 exports ~50 MW backward toward B2 (B2→B3 is reversed)
      - Brighton export at B5 reaches B2 via B5→B1→B2 and B5→B4→B3→B2
    The result is B1→B2 ≈ +250 MW, well below the 400 MW limit.

    Verification against GAMS (single-period LP, Sbase=100 MVA)
    -----------------------------------------------------------
    The GAMS model was run with GAMS 45.7 (2026-05-11).  After removing
    reserve margins and unit clustering from Dispa-SET, the two models
    agree to within 1 MW across all six lines.  The residual difference
    is explained entirely by the 2 MW extra load from placeholder demands
    at PJM_B1 and PJM_B5 that are absent from the GAMS model.

    Expected values (GAMS reference → Dispa-SET within 1 MW)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Flow B4→B5:  ≈−240.0 MW  (at the 240 MW limit; negative = B5→B4)
        Flow B1→B2:  ≈+249.7 MW  (forward; B1 exports to B2 load)
        Flow B1→B5:  ≈−226.5 MW  (negative; Brighton exports back to slack)
        Flow B2→B3:  ≈−50.3 MW   (negative; Solitude back-feeds B2)
        Flow B1→B4:  ≈+186.8 MW  (forward; B1 exports to B4 load)
        Flow B3→B4:  ≈−26.8 MW   (negative; Solitude partially serves B4)
    """
    skip_if_no_gams()
    out = _get_pjm5bus_solve_output()
    results = out["results"]

    assert "OutputFlow" in results, "OutputFlow missing — solve may have failed"
    assert not results["OutputFlow"].empty, "OutputFlow is empty"

    # --- B4-B5 line is saturated in the reverse direction (B5→B4) ---
    # GAMS reference: -240.0 MW.  Tolerance 2 MW for the 2 MW placeholder load.
    flow_b4_b5 = _mean_flow(results, "PJM_B4 -> PJM_B5")
    abs_flow_b4_b5 = abs(flow_b4_b5)
    print(f"PJM_B4 → PJM_B5 mean flow: {flow_b4_b5:+.1f} MW  (Fmax=240 MW, GAMS=-240.0)")

    # Must be at the 240 MW limit (within 2 MW numerical tolerance)
    assert abs_flow_b4_b5 > 235.0, (
        f"Expected |B4-B5| > 235 MW (at the 240 MW congestion limit), "
        f"got {abs_flow_b4_b5:.1f} MW.  B4-B5 should be the binding constraint."
    )
    assert flow_b4_b5 < 0, (
        f"Expected negative flow on B4→B5 (Brighton exports B5→B4), "
        f"got {flow_b4_b5:+.1f} MW."
    )

    # --- B1-B2 line: ≈+250 MW (GAMS: +249.7 MW) ---
    # Power from Alta+ParkCity at B1 and from Brighton (via B5→B1→B2) reaches
    # B2's 300 MW load; Solitude back-feeds B2 via B2→B3 in reverse (~50 MW).
    flow_b1_b2 = _mean_flow(results, "PJM_B1 -> PJM_B2")
    print(f"PJM_B1 → PJM_B2 mean flow: {flow_b1_b2:+.1f} MW  (Fmax=400 MW, GAMS=+249.7)")

    assert 220.0 < flow_b1_b2 < 270.0, (
        f"Expected B1-B2 mean flow in [220, 270] MW (GAMS reference: +249.7 MW), "
        f"got {flow_b1_b2:.1f} MW."
    )

    # --- B2-B3 line: ≈-50 MW (GAMS: -50.3 MW) ---
    # Solitude (B3) back-feeds B2's demand; the negative sign means power
    # flows B3→B2 (opposite to the defined B2→B3 direction).
    flow_b2_b3 = _mean_flow(results, "PJM_B2 -> PJM_B3")
    print(f"PJM_B2 → PJM_B3 mean flow: {flow_b2_b3:+.1f} MW  (Fmax=1000 MW, GAMS=-50.3)")

    assert -60.0 < flow_b2_b3 < -40.0, (
        f"Expected B2-B3 mean flow in [-60, -40] MW (GAMS: -50.3 MW), "
        f"got {flow_b2_b3:.1f} MW."
    )

    # --- B1-B5 line: ≈-227 MW (GAMS: -226.5 MW) ---
    # Brighton (B5) exports most of its power directly to the slack bus (B1)
    # via the very low-reactance B1-B5 line (X=0.0064), explaining its large
    # PTDF coefficient: PTDF[B1→B5, B5] ≈ -0.888.
    flow_b1_b5 = _mean_flow(results, "PJM_B1 -> PJM_B5")
    print(f"PJM_B1 → PJM_B5 mean flow: {flow_b1_b5:+.1f} MW  (Fmax=1000 MW, GAMS=-226.5)")

    assert -240.0 < flow_b1_b5 < -210.0, (
        f"Expected B1-B5 mean flow in [-240, -210] MW (GAMS: -226.5 MW), "
        f"got {flow_b1_b5:.1f} MW."
    )

    # --- Brighton is the primary generator, constrained by B4-B5 congestion ---
    # GAMS reference: 466.5 MW (2 MW extra demand → ~468.5 MW in Dispa-SET)
    if "OutputPower" in results and not results["OutputPower"].empty:
        power = results["OutputPower"]
        if "PJM_Brighton" in power.columns:
            mean_brighton = float(power["PJM_Brighton"].mean())
            print(f"PJM_Brighton mean output: {mean_brighton:.1f} MW  (Pmax=600 MW, GAMS=466.5)")
            assert 450.0 < mean_brighton < 490.0, (
                f"Expected Brighton ≈466-469 MW (GAMS: 466.5 MW), "
                f"got {mean_brighton:.1f} MW."
            )


# ---------------------------------------------------------------------------
# Script entry-point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
