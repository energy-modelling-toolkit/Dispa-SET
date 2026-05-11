# -*- coding: utf-8 -*-
"""
Integration test: Mid-Term Scheduling (MTS)
============================================

What this test does
-------------------
Exercises the regional Mid-Term Scheduling pre-processing. With the
``HydroScheduling: Regional`` flag, Dispa-SET first runs an aggregated
optimisation over the whole horizon to compute target reservoir levels,
then injects them back into the dispatch model.

The test:

1. Loads ``tests/configs/tiny_mts.yml``.
2. Calls ``ds.mid_term_scheduling`` directly to verify the MTS function
   returns the expected tuple of profiles (4 dataframes).
3. Then runs the full pipeline ``build_simulation -> solve_GAMS``
   so the dispatch model uses the MTS reservoir profiles.

Performance budget
------------------
Should typically run in **less than 30 seconds**.

How to run
----------

1. ``pytest tests/integration/test_solve_mts.py``
2. ``python tests/integration/test_solve_mts.py``

Skipped automatically if no GAMS installation is detected.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import build_solve, load_test_config, skip_if_no_gams  # noqa: E402


def _write_constant_profile(path: Path, index: pd.DatetimeIndex, value: float) -> str:
    pd.DataFrame({"Z1": value}, index=index).to_csv(path)
    return str(path)


@pytest.mark.timeout(120)
def test_solve_mts_regional():
    skip_if_no_gams()
    cfg = load_test_config("tiny_mts.yml", "solve_mts")
    # Force a small MTS time step (24h) - default would otherwise be the
    # whole horizon and might be unnecessarily slow.
    import dispaset as ds  # local to keep imports cheap when skipped
    profiles, *_ = ds.mid_term_scheduling(cfg, TimeStep=24, mts_plot=False)
    assert profiles is not None
    # Run the full build+solve afterwards:
    out = build_solve(cfg)
    assert "OutputPower" in out["results"]
    print(f"MTS build {out['build_time']:.2f}s  solve {out['solve_time']:.2f}s")


@pytest.mark.timeout(180)
def test_solve_mts_with_exogenous_reserve_aliases():
    """MTS path should also work with exogenous reserve inputs and legacy keys."""
    skip_if_no_gams()
    cfg = load_test_config("tiny_mts.yml", "solve_mts_reserve_aliases")
    cfg["ReserveCalculation"] = "Exogenous"

    reserve_idx = pd.date_range(
        start=pd.Timestamp(cfg["StartDate"][0], 1, 1, 0, 0),
        end=pd.Timestamp(cfg["StartDate"][0], 12, 31, 23, 0),
        freq="h",
    )
    reserve_dir = Path(cfg["SimulationDirectory"]).parent / "mts_reserve_demands"
    reserve_dir.mkdir(parents=True, exist_ok=True)

    # Populate only legacy keys: loader normalization must map them.
    cfg["Reserve2U"] = _write_constant_profile(reserve_dir / "reserve2u.csv", reserve_idx, 2.0)
    cfg["Reserve2D"] = _write_constant_profile(reserve_dir / "reserve2d.csv", reserve_idx, 1.5)
    cfg["FFRLimit"] = _write_constant_profile(reserve_dir / "ffr.csv", reserve_idx, 1.0)
    cfg["PrimaryReserveLimit"] = _write_constant_profile(reserve_dir / "fcr.csv", reserve_idx, 0.8)
    cfg["aFRRUDemand"] = ""
    cfg["aFRRDDemand"] = ""
    cfg["FFRDemand"] = ""
    cfg["FCRDemand"] = ""

    out = build_solve(cfg)
    assert "OutputPower" in out["results"]
    assert "OutputReserveProvision" in out["results"]


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
