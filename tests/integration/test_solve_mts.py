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

import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import build_solve, load_test_config, skip_if_no_gams  # noqa: E402


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


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
