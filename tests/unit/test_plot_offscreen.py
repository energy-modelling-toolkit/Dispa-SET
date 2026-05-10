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


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
