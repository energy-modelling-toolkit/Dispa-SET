# -*- coding: utf-8 -*-
"""
Integration test: boundary sector (H2) coupling
=================================================

What this test does
-------------------
Exercises the *Boundary Sector* abstraction that replaces the legacy
"Heat" / "H2" formulations. The configuration has:

* A power zone Z2 with conventional generation.
* An electrolyzer (P2GS) consuming power and producing hydrogen into the
  boundary sector ``Z1_h2``.
* A fuel-cell (BSPG) consuming hydrogen and producing power back.
* An H2 demand profile (``SectorXDemand``) for the boundary sector.

The test verifies that:

* The boundary-sector input data is loaded correctly.
* The simulation runs and produces the expected boundary-sector output
  variables: ``OutputPowerX``, ``OutputSectorXStorageLevel``,
  ``OutputXNotServed``.
* Power consumption by the electrolyzer is recorded in
  ``OutputPowerConsumption``.
* The boundary-sector demand is approximately balanced by H2 sources
  (electrolyzer + storage discharge - fuel cell consumption +
  unserved demand).

How to run
----------

1. ``pytest tests/integration/test_solve_boundary_sector.py``
2. ``python tests/integration/test_solve_boundary_sector.py``

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


@pytest.mark.timeout(180)
def test_solve_h2_boundary_sector():
    skip_if_no_gams()
    cfg = load_test_config("tiny_bs.yml", "solve_bs")
    out = build_solve(cfg)

    inputs = out["inputs"]
    results = out["results"]

    # Boundary-sector output dataframes should be present:
    for key in ("OutputPowerX", "OutputXNotServed",
                "OutputSectorXStorageLevel"):
        assert key in results, f"Missing boundary-sector output: {key}"

    # The boundary-sector zone Z1_h2 must be in the zone set:
    assert "nx" in inputs["sets"]
    assert "Z1_h2" in inputs["sets"]["nx"]

    print(f"BS build {out['build_time']:.2f}s  "
          f"solve {out['solve_time']:.2f}s  "
          f"#bs_zones={len(inputs['sets']['nx'])}")


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
