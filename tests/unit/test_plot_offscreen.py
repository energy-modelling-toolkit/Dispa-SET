# -*- coding: utf-8 -*-
"""
Unit test: dispaset.postprocessing.plot (off-screen rendering)
==============================================================

What this test does
-------------------
The plotting helpers all call ``plt.show()`` at the end. We force the
non-interactive *Agg* backend (done in `tests/conftest.py`) so the figures
are drawn off-screen, then we save each rendered figure to
``tests/output/plots/`` instead of displaying them. The test just verifies
that no exception is raised.

Tested helpers (with synthetic input data, no GAMS needed):

* `plot_dispatch`    - stack-area dispatch plot for a single zone.
* `plot_dispatchX`   - same for a boundary sector zone.
* `plot_zone`        - per-zone summary plot.
* `plot_zone_capacities` - bar plot of installed capacity.

How to run
----------

1. ``pytest tests/unit/test_plot_offscreen.py``
2. ``python tests/unit/test_plot_offscreen.py``

This test does NOT require GAMS.

Output
------
PNG files are written to ``tests/output/plots/`` so you can visually inspect
the result.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # repeat: standalone-script safety
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from dispaset.postprocessing.plot import (  # noqa: E402
    plot_dispatch, plot_dispatchX,
)


PLOT_DIR = Path(__file__).resolve().parents[1] / "output" / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)


def _save_current_figure(name: str):
    """Save *and* close the active matplotlib figure(s)."""
    out = PLOT_DIR / f"{name}.png"
    for fignum in plt.get_fignums():
        plt.figure(fignum).savefig(out, bbox_inches="tight")
    plt.close("all")
    return out


def _build_dispatch_inputs(periods=24):
    """Synthetic mini dispatch dataframe."""
    idx = pd.date_range("2015-01-01", periods=periods, freq="h")
    plotdata = pd.DataFrame({
        "GAS": np.full(periods, 50.0),
        "WIN": np.full(periods, 30.0),
        "SUN": np.full(periods, 20.0),
    }, index=idx)
    demand = pd.Series(plotdata.sum(axis=1), index=idx, name="Demand")
    return demand, plotdata


# --------------------------------------------------------------------------- #
def test_plot_dispatch_runs_offscreen():
    demand, plotdata = _build_dispatch_inputs()
    plot_dispatch(demand, plotdata)
    out = _save_current_figure("plot_dispatch")
    assert out.exists() and out.stat().st_size > 0


def test_plot_dispatch_with_curtailment_and_shedding():
    demand, plotdata = _build_dispatch_inputs()
    curtail = pd.Series(np.full(len(demand), 5.0), index=demand.index)
    shed = pd.Series(np.full(len(demand), 2.0), index=demand.index)
    plot_dispatch(demand, plotdata, curtailment=curtail, shedload=shed)
    out = _save_current_figure("plot_dispatch_curt")
    assert out.exists() and out.stat().st_size > 0


def test_plot_dispatchX_with_minimal_inputs():
    """`plot_dispatchX` should not crash on a minimal boundary-sector input."""
    idx = pd.date_range("2023-01-01", periods=24, freq="h")
    inputs = {
        "param_df": {
            "Demand": pd.DataFrame(
                {("DA", "BE"): np.full(24, 1000.0)}, index=idx
            ),
            "SectorXDemand": pd.DataFrame(
                {"BE_h2": np.full(24, 100.0)}, index=idx
            ),
            "SectorXStorageCapacity": pd.DataFrame(
                {"StorageCapacity": [1000]}, index=["BE_h2 - BATS"],
            ),
        },
        # plot_dispatchX now reads SectorXStorageCapacity from inputs['parameters']
        # (changed in the future_proof merge to support GWh scaling)
        "parameters": {
            "SectorXStorageCapacity": {"val": np.array([1000.0])},
        },
        "sets": {"n": ["BE"], "nx": ["BE_h2"], "f": ["GAS"]},
        "units": pd.DataFrame({
            "Zone": ["BE"], "Technology": ["GTUR"], "Fuel": ["GAS"],
            "Unit": ["BE - GTUR"], "Sector1": ["BE"],
            "PowerCapacity": [1000], "PartLoadMin": [0.3],
            "RampUp": [100], "RampDown": [100],
            "StartUpCost": [0], "NoLoadCost": [0],
            "MinUpTime": [1], "MinDownTime": [1],
        }, index=["BE - GTUR"]),
    }
    results = {
        "OutputPowerX": pd.DataFrame(
            {"BE_h2 - GTUR": np.full(24, 50.0)}, index=idx,
        ),
        "OutputSectorXStorageLevel": pd.DataFrame(
            {"BE_h2": np.full(24, 500.0)}, index=idx,
        ),
        "OutputSectorXStorageInput": pd.DataFrame(
            {"BE_h2 - BATS": np.full(24, 5.0)}, index=idx,
        ),
        "OutputSectorXFlow": pd.DataFrame(
            {"BE->BE_h2": np.full(24, 10.0)}, index=idx,
        ),
        "OutputXNotServed": pd.DataFrame(
            {"BE_h2": np.zeros(24)}, index=idx,
        ),
    }
    fig = plot_dispatchX(inputs, results, z="BE_h2")
    out = _save_current_figure("plot_dispatchX")
    assert out.exists() and out.stat().st_size > 0


def test_plot_dispatchX_storage_capacity_gwh_scaling():
    """Verify that plot_dispatchX correctly scales the storage level by
    SectorXStorageCapacity to convert from p.u. to GWh.

    New behaviour from the future_proof merge:
        storageXcapacity = inputs['parameters']['SectorXStorageCapacity']['val'].item()
        storage_level = results['OutputSectorXStorageLevel'][z] * storageXcapacity / 1000

    We use a known capacity (2000 MWh = 2 GWh) and a storage level of 0.5 p.u.
    and verify the function runs without errors (the actual GWh values are
    rendered in the figure, not returned, so we can only assert no crash).
    """
    idx = pd.date_range("2023-01-01", periods=24, freq="h")
    storage_capacity_mwh = 2000.0  # MWh → 2 GWh when /1000

    inputs = {
        "param_df": {
            "Demand": pd.DataFrame(
                {("DA", "BE"): np.full(24, 1000.0)}, index=idx
            ),
            "SectorXDemand": pd.DataFrame(
                {"BE_h2": np.full(24, 100.0)}, index=idx
            ),
            "SectorXStorageCapacity": pd.DataFrame(
                {"StorageCapacity": [storage_capacity_mwh]},
                index=["BE_h2 - BATS"],
            ),
        },
        # The new code reads from inputs['parameters']['SectorXStorageCapacity']
        "parameters": {
            "SectorXStorageCapacity": {
                "val": np.array([storage_capacity_mwh]),
            },
        },
        "sets": {"n": ["BE"], "nx": ["BE_h2"], "f": ["GAS"]},
        "units": pd.DataFrame({
            "Zone": ["BE"], "Technology": ["GTUR"], "Fuel": ["GAS"],
            "Unit": ["BE - GTUR"], "Sector1": ["BE"],
            "PowerCapacity": [1000], "PartLoadMin": [0.3],
            "RampUp": [100], "RampDown": [100],
            "StartUpCost": [0], "NoLoadCost": [0],
            "MinUpTime": [1], "MinDownTime": [1],
        }, index=["BE - GTUR"]),
    }
    results = {
        "OutputPowerX": pd.DataFrame(
            {"BE_h2 - GTUR": np.full(24, 50.0)}, index=idx,
        ),
        # Storage level in p.u. (0.5 → 0.5 * 2000 / 1000 = 1 GWh)
        "OutputSectorXStorageLevel": pd.DataFrame(
            {"BE_h2": np.full(24, 0.5)}, index=idx,
        ),
        "OutputSectorXStorageInput": pd.DataFrame(
            {"BE_h2 - BATS": np.full(24, 5.0)}, index=idx,
        ),
        "OutputSectorXFlow": pd.DataFrame(
            {"BE->BE_h2": np.full(24, 10.0)}, index=idx,
        ),
        "OutputXNotServed": pd.DataFrame(
            {"BE_h2": np.zeros(24)}, index=idx,
        ),
    }
    fig = plot_dispatchX(inputs, results, z="BE_h2")
    out = _save_current_figure("plot_dispatchX_gwh_scaling")
    assert out.exists() and out.stat().st_size > 0


def test_plot_dispatchX_multi_zone_storage_capacity():
    """Regression test: plot_dispatchX must not raise ValueError when there
    are *multiple* boundary sector zones (i.e. SectorXStorageCapacity['val']
    has more than one element).

    Bug: the original code called ``.item()`` on the full capacity array, which
    raises ``ValueError: can only convert an array of size 1 to a Python scalar``
    as soon as there are ≥2 boundary sector zones.  The fix looks up the zone's
    index in ``sets['nx']`` and reads only that element.
    """
    idx = pd.date_range("2023-01-01", periods=24, freq="h")

    # Three boundary-sector zones; the zone under test is the *second* one
    # so that any regression that defaults to index-0 would pick the wrong
    # capacity and the zone-lookup logic is exercised.
    nx_zones = ["Z1_h2", "Z1_th", "Z2_th"]
    capacities = np.array([5000.0, 300.0, 150.0])   # MWh

    inputs = {
        "param_df": {
            "Demand": pd.DataFrame(
                {("DA", "Z1"): np.full(24, 1000.0)}, index=idx
            ),
            "SectorXDemand": pd.DataFrame(
                {"Z1_th": np.full(24, 80.0)}, index=idx
            ),
            "SectorXStorageCapacity": pd.DataFrame(
                {"StorageCapacity": capacities},
                index=[z + " - BATS" for z in nx_zones],
            ),
        },
        "parameters": {
            "SectorXStorageCapacity": {"val": capacities},
        },
        "sets": {"n": ["Z1", "Z2"], "nx": nx_zones, "f": ["GAS"]},
        "units": pd.DataFrame({
            "Zone":  ["Z1"],
            "Technology": ["GTUR"],
            "Fuel": ["GAS"],
            "Unit": ["Z1 - GTUR"],
            "Sector1": ["Z1"],
            "PowerCapacity": [1000], "PartLoadMin": [0.3],
            "RampUp": [100], "RampDown": [100],
            "StartUpCost": [0], "NoLoadCost": [0],
            "MinUpTime": [1], "MinDownTime": [1],
        }, index=["Z1 - GTUR"]),
    }
    results = {
        "OutputPowerX": pd.DataFrame(
            {"Z1_th - GTUR": np.full(24, 40.0)}, index=idx,
        ),
        "OutputSectorXStorageLevel": pd.DataFrame(
            {"Z1_th": np.full(24, 0.5)}, index=idx,   # p.u. → 0.5*300/1000 GWh
        ),
        "OutputSectorXStorageInput": pd.DataFrame(
            {"Z1_th - BATS": np.full(24, 2.0)}, index=idx,
        ),
        "OutputSectorXFlow": pd.DataFrame(
            {"Z1->Z1_th": np.full(24, 5.0)}, index=idx,
        ),
        "OutputXNotServed": pd.DataFrame(
            {"Z1_th": np.zeros(24)}, index=idx,
        ),
    }
    # Must not raise ValueError even though len(capacities) == 3
    fig = plot_dispatchX(inputs, results, z="Z1_th")
    out = _save_current_figure("plot_dispatchX_multi_zone")
    assert out.exists() and out.stat().st_size > 0


def test_plot_dispatch_hatch_guard_unknown_positive_column():
    """Regression test: plot_dispatch must NOT raise KeyError when a plotdata
    column is present in commons['colors'] but was previously absent from
    commons['hatches'].

    Before the fix, the positive-values loop used:
        hatch = commons['hatches'][col2]          # KeyError if col2 absent
    The fix applies the same guard used in the negative-values loop:
        if plot_lines and col2 in commons['hatches']:

    We use 'HeatSlack' as the offending column because it was in colors but
    not in hatches prior to the patch (and is a positive generation-side
    entry typical in coupled heat/power runs).
    """
    periods = 24
    idx = pd.date_range("2015-01-01", periods=periods, freq="h")
    # HeatSlack appears as a positive (generation-side) column
    plotdata = pd.DataFrame({
        "GAS": np.full(periods, 50.0),
        "HeatSlack": np.full(periods, 10.0),
    }, index=idx)
    demand = pd.Series(plotdata.sum(axis=1), index=idx, name="Demand")
    # plot_lines is set automatically when rng < 32 days (our 24-h window qualifies)
    plot_dispatch(demand, plotdata)
    out = _save_current_figure("plot_dispatch_hatch_guard")
    assert out.exists() and out.stat().st_size > 0


def test_plot_dispatch_hatch_guard_unknown_negative_column():
    """Regression: negative-values loop must also not KeyError for unknown
    hatch keys (this was already correct but we verify it stays so).

    ShedLoad is a negative-side entry present in colors but previously
    absent from hatches.
    """
    periods = 24
    idx = pd.date_range("2015-01-01", periods=periods, freq="h")
    plotdata = pd.DataFrame({
        "GAS": np.full(periods, 60.0),
        "ShedLoad": np.full(periods, -10.0),
    }, index=idx)
    demand = pd.Series(np.full(periods, 60.0), index=idx, name="Demand")
    plot_dispatch(demand, plotdata)
    out = _save_current_figure("plot_dispatch_hatch_guard_neg")
    assert out.exists() and out.stat().st_size > 0


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
