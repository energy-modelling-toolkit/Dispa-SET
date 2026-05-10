# -*- coding: utf-8 -*-
"""
Unit test: read & analyse pre-built reference results
======================================================

What this test does
-------------------
The repository ships a pre-built simulation result under
``tests/data/reference/`` (Inputs.gdx + Inputs.p + Results.gdx). This lets
us exercise the full postprocessing pipeline (gdx -> dataframe ->
indicators -> plots) WITHOUT having to re-run a GAMS optimisation.

The test:

1. loads the reference simulation through ``ds.get_sim_results``;
2. checks that the standard input and result keys are present;
3. computes ``get_indicators_powerplant``, ``get_result_analysis`` and
   ``aggregate_by_fuel`` and checks the output types;
4. invokes ``plot_zone``, ``plot_zone_capacities`` and
   ``plot_energy_zone_fuel`` off-screen and stores the PNG outputs in
   ``tests/output/plots_reference/``.

NOTE
----
This test still requires a working GAMS API to read the GDX file
(``gdxcc``). On a machine without GAMS the test is skipped.

How to run
----------

1. ``pytest tests/unit/test_postprocess_reference.py``
2. ``python tests/unit/test_postprocess_reference.py``
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import dispaset as ds  # noqa: E402
from dispaset.misc.gdx_handler import package_exists, get_gams_path  # noqa: E402


REFERENCE_DIR = Path(__file__).resolve().parents[1] / "data" / "reference"
PLOT_DIR = Path(__file__).resolve().parents[1] / "output" / "plots_reference"
PLOT_DIR.mkdir(parents=True, exist_ok=True)


def _skip_if_no_gams():
    if not package_exists("gams"):
        pytest.skip("GAMS python API not installed")
    if get_gams_path() is None:
        pytest.skip("GAMS installation not detected")


@pytest.fixture(scope="module")
def sim_results():
    _skip_if_no_gams()
    if not (REFERENCE_DIR / "Inputs.p").is_file():
        pytest.skip(f"Reference simulation not found at {REFERENCE_DIR}")
    return ds.get_sim_results(path=str(REFERENCE_DIR))


def test_get_sim_results_returns_dicts(sim_results):
    inputs, results = sim_results
    assert isinstance(inputs, dict)
    assert isinstance(results, dict)
    for key in ("sets", "units", "param_df"):
        assert key in inputs, f"Missing inputs key {key!r}"
    assert "OutputPower" in results, "Missing OutputPower in results"


def test_get_indicators_powerplant(sim_results):
    inputs, results = sim_results
    indi = ds.get_indicators_powerplant(inputs, results)
    assert hasattr(indi, "Generation")  # DataFrame column
    assert "Zone" in indi.columns
    assert "Fuel" in indi.columns


def test_get_result_analysis(sim_results):
    """`get_result_analysis` returns a dict and prints summary stats."""
    inputs, results = sim_results
    out = ds.get_result_analysis(inputs, results)
    assert isinstance(out, dict)


def test_aggregate_by_fuel(sim_results):
    inputs, results = sim_results
    by_fuel = ds.aggregate_by_fuel(results["OutputPower"], inputs)
    assert by_fuel.shape[0] == len(results["OutputPower"].index)


# --- Plot tests ---------------------------------------------------------- #
def test_plot_zone(sim_results):
    inputs, results = sim_results
    z = inputs["sets"]["n"][0]
    ds.plot_zone(inputs, results, z=z)
    out = PLOT_DIR / f"plot_zone_{z}.png"
    for fignum in plt.get_fignums():
        plt.figure(fignum).savefig(out, bbox_inches="tight")
    plt.close("all")
    assert out.exists() and out.stat().st_size > 0


def test_plot_zone_capacities(sim_results):
    inputs, results = sim_results
    ds.plot_zone_capacities(inputs, results, plot=True)
    out = PLOT_DIR / "plot_zone_capacities.png"
    for fignum in plt.get_fignums():
        plt.figure(fignum).savefig(out, bbox_inches="tight")
    plt.close("all")
    assert out.exists() and out.stat().st_size > 0


def test_plot_energy_zone_fuel(sim_results):
    inputs, results = sim_results
    indi = ds.get_indicators_powerplant(inputs, results)
    ds.plot_energy_zone_fuel(inputs, results, indi)
    out = PLOT_DIR / "plot_energy_zone_fuel.png"
    for fignum in plt.get_fignums():
        plt.figure(fignum).savefig(out, bbox_inches="tight")
    plt.close("all")
    assert out.exists() and out.stat().st_size > 0


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
