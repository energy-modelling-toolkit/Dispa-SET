# -*- coding: utf-8 -*-
"""
Integration test: postprocess + plot a freshly-solved simulation
================================================================

What this test does
-------------------
Most postprocessing helpers are exercised through ``test_postprocess_*``
unit tests using mocked dataframes. This test goes one step further and
runs them on the **actual** output of a small build+solve so that any
incompatibility between the GDX content and the postprocessing code
surfaces immediately.

The test:

1. Builds and solves the tiny 2-zone case.
2. Runs ``ds.get_indicators_powerplant``, ``ds.get_result_analysis``,
   ``ds.aggregate_by_fuel``, ``ds.CostExPost`` and
   ``ds.get_units_operation_cost``.
3. Renders ``ds.plot_zone``, ``ds.plot_zone_capacities`` and
   ``ds.plot_energy_zone_fuel`` off-screen and saves the PNGs to
   ``tests/output/plots_solve/``.

How to run
----------

1. ``pytest tests/integration/test_postprocess_after_solve.py``
2. ``python tests/integration/test_postprocess_after_solve.py``

Skipped automatically if no GAMS installation is detected.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import build_solve, load_test_config, skip_if_no_gams  # noqa: E402


PLOT_DIR = Path(__file__).resolve().parents[1] / "output" / "plots_solve"
PLOT_DIR.mkdir(parents=True, exist_ok=True)


def _save_active_figs(name: str):
    out = PLOT_DIR / f"{name}.png"
    for fignum in plt.get_fignums():
        plt.figure(fignum).savefig(out, bbox_inches="tight")
    plt.close("all")
    return out


@pytest.mark.timeout(180)
def test_postprocess_after_solve():
    skip_if_no_gams()
    cfg = load_test_config("tiny_2zones.yml", "postprocess_after_solve")
    out = build_solve(cfg)
    inputs, results = out["inputs"], out["results"]

    import dispaset as ds

    # --- numerical postprocessing helpers ---
    indi = ds.get_indicators_powerplant(inputs, results)
    assert "Generation" in indi.columns

    analysis = ds.get_result_analysis(inputs, results)
    assert isinstance(analysis, dict)

    by_fuel = ds.aggregate_by_fuel(results["OutputPower"], inputs)
    assert by_fuel.shape[0] == len(results["OutputPower"].index)

    cep = ds.CostExPost(inputs, results)
    assert cep is not None

    uoc = ds.get_units_operation_cost(inputs, results)
    assert uoc is not None

    # --- plotting helpers (off-screen) ---
    z = inputs["sets"]["n"][0]
    ds.plot_zone(inputs, results, z=z)
    out1 = _save_active_figs(f"plot_zone_{z}")
    assert out1.exists() and out1.stat().st_size > 0

    ds.plot_zone_capacities(inputs, results, plot=True)
    out2 = _save_active_figs("plot_zone_capacities")
    assert out2.exists() and out2.stat().st_size > 0

    ds.plot_energy_zone_fuel(inputs, results, indi)
    out3 = _save_active_figs("plot_energy_zone_fuel")
    assert out3.exists() and out3.stat().st_size > 0


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
