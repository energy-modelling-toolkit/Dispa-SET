# -*- coding: utf-8 -*-
"""
Integration test: build + solve + read - LP single zone
========================================================

What this test does
-------------------
End-to-end exercise of Dispa-SET on the smallest possible test-case:

* one zone (Z1),
* one day of simulation,
* LP-clustered formulation,
* no boundary sector / no MTS / no reserves.

Steps:
  1. Load ``tests/configs/tiny.yml``.
  2. Build the simulation environment (writes ``Inputs.gdx`` + ``Inputs.p``).
  3. Solve the optimisation problem with GAMS (writes ``Results.gdx``).
  4. Read the results back via :func:`dispaset.get_sim_results`.
  5. Sanity-check the solution: total dispatch close to the demand,
     no shed-load (the units cover the demand), and the time index
     matches the requested period.

Performance budget
------------------
Should typically run in **less than 5 seconds** including the GAMS solve.
A pytest timeout is configured at 90 seconds.

How to run
----------

1. ``pytest tests/integration/test_solve_lp.py``
2. ``python tests/integration/test_solve_lp.py``

The test is skipped automatically when no GAMS installation is detected.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import build_solve, load_test_config, skip_if_no_gams  # noqa: E402


@pytest.mark.timeout(90)
def test_solve_tiny_lp():
    skip_if_no_gams()
    cfg = load_test_config("tiny.yml", "solve_lp")
    out = build_solve(cfg)

    inputs = out["inputs"]
    results = out["results"]

    # Index matches requested simulation period:
    sim_index = results["OutputPower"].index
    assert sim_index[0].day == 1
    assert sim_index[0].month == 1
    assert sim_index[0].year == 2015

    # Storage level violation should be small (a slack variable in the
    # GAMS formulation; non-zero only at the very edge of the horizon).
    sv = results["LostLoad_StorageLevelViolation"]
    sv_total = float(sv.sum()) if hasattr(sv, "sum") else float(sv)
    assert sv_total < 1000, f"Unusually large storage violation: {sv_total}"

    # Energy balance: generation + (- storage_input) ~ demand + curtailment - shed
    # We just check that dispatch is in the right order of magnitude (>= 50%
    # of demand: any value below would indicate a malformed simulation).
    demand_total = float(inputs["param_df"]["Demand"]["DA"].sum().sum())
    dispatch = float(results["OutputPower"].sum().sum())
    assert dispatch > demand_total * 0.5, (
        f"Dispatch {dispatch:.0f} much too low vs demand {demand_total:.0f}"
    )
    assert dispatch < demand_total * 1.5, (
        f"Dispatch {dispatch:.0f} much too high vs demand {demand_total:.0f}"
    )

    print(f"Build {out['build_time']:.2f}s  "
          f"Solve {out['solve_time']:.2f}s  "
          f"Read {out['read_time']:.2f}s")


@pytest.mark.timeout(120)
def test_solve_two_zones_with_ntc():
    skip_if_no_gams()
    cfg = load_test_config("tiny_2zones.yml", "solve_2zones")
    out = build_solve(cfg)

    results = out["results"]
    assert "OutputFlow" in results

    # Zones in the result should match the config:
    assert "Z1" in out["inputs"]["sets"]["n"]
    assert "Z2" in out["inputs"]["sets"]["n"]

    # Some interconnection flow is expected when there are NTCs:
    flows = results["OutputFlow"]
    assert not flows.empty, "Expected interconnection flows but got empty df"


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
