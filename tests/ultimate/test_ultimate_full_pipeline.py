# -*- coding: utf-8 -*-
"""
Ultimate end-to-end test
==========================

What this test does
-------------------
This is the "everything together" test: it runs the full Dispa-SET
pipeline (load_config -> build -> solve_GAMS -> get_sim_results ->
postprocess -> plots) on the realistic ``ConfigTest`` fictitious test
case but on a **3-day window** instead of a full year.

It exercises:

* Multi-zone configuration with NTC interconnections.
* CHP units with thermal boundary sector.
* Boundary sector (H2) with electrolyzer / fuel cell units.
* Hydro storage and reservoirs (mid-term scheduling typically off here).
* All the headline plotting and result-analysis APIs.

The test is intentionally tolerant: it checks that each step *runs*
without raising and that the most important outputs are non-empty.
Detailed value assertions are left to the smaller, focused integration
tests in ``tests/integration``.

How to run
----------

1. ``pytest tests/ultimate/test_ultimate_full_pipeline.py``
2. ``python tests/ultimate/test_ultimate_full_pipeline.py``

Skipped automatically if no GAMS installation is detected.

Expected runtime: under ~30 seconds on a recent laptop.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import (  # noqa: E402
    skip_if_no_gams,
    load_test_config,
    fresh_simdir,
    solve_succeeded,
)

import dispaset as ds  # noqa: E402

OUTPUT_DIR = REPO_ROOT / "tests" / "output" / "ultimate_plots"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def _save_current_fig(name: str) -> None:
    """Save the currently active matplotlib figure (off-screen) and close
    every open figure to avoid accumulating memory across the test."""
    fig = plt.gcf()
    fig.savefig(OUTPUT_DIR / f"{name}.png", dpi=80, bbox_inches="tight")
    plt.close("all")


@pytest.mark.timeout(180)
def test_ultimate_pipeline():
    skip_if_no_gams()
    cfg = load_test_config("ultimate.yml", "ultimate_pipeline")

    SimData = ds.build_simulation(cfg)
    assert isinstance(SimData, dict)
    assert "sets" in SimData and "parameters" in SimData

    ds.solve_GAMS(cfg["SimulationDirectory"], cfg["GAMS_folder"])
    assert solve_succeeded(cfg["SimulationDirectory"]), (
        f"GAMS did not produce Results.gdx in {cfg['SimulationDirectory']}"
    )

    inputs, results = ds.get_sim_results(
        path=cfg["SimulationDirectory"], cache=False
    )
    assert isinstance(inputs, dict)
    assert isinstance(results, dict)
    assert len(results) > 0

    expected_outputs = ["OutputPower", "OutputSystemCost"]
    for key in expected_outputs:
        assert key in results, f"Missing expected output {key}"
        assert not results[key].empty

    # ---- Energy balance (all zones) ------------------------------------
    balance = ds.check_energy_balance(inputs, results)
    for zone, rel_err in balance.items():
        assert rel_err <= 0.01, (
            f"Energy balance exceeds 1 % for zone {zone}: {rel_err * 100:.3f}%"
        )

    # ---- Plotting (all off-screen) ---------------------------------------
    ds.plot_zone(inputs, results)
    _save_current_fig("plot_zone")

    if any(z.startswith("Z1_h2") or "_h2" in z for z in inputs.get("sets", {}).get("nx", [])):
        try:
            ds.plot_dispatchX(inputs, results, z="Z1_h2")
            _save_current_fig("plot_dispatchX")
        except Exception as e:  # noqa: BLE001
            # We don't want to fail the ultimate test on optional plot
            # of a possibly missing boundary sector.
            print(f"plot_dispatchX skipped: {e}")

    ds.plot_zone_capacities(inputs, results)
    _save_current_fig("plot_zone_capacities")

    indicators = ds.get_indicators_powerplant(inputs, results)
    assert indicators is not None
    ds.plot_energy_zone_fuel(inputs, results, indicators)
    _save_current_fig("plot_energy_zone_fuel")

    analysis = ds.get_result_analysis(inputs, results)
    assert isinstance(analysis, dict)


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
