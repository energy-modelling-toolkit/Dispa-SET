# -*- coding: utf-8 -*-
"""
Integration test: CHP coupling through the boundary sector
==========================================================

What this test does
-------------------
The boundary sector formulation is also used to model heat zones for CHP
plants. The configuration `tiny_bs.yml` includes:

* the **back-pressure** CHP "Diemen" (Z2) coupled to the heat zone Z2_th_1,
* the **extraction** CHP "Rijnmond II" (Z2) coupled to Z2_th_2,
* the **P2H** unit "Ground source HP" (Z1) coupled to Z2_th_2.

The test verifies that:

* The thermal sectors are part of the boundary-sector zone set.
* The CHP units appear in the inputs.
* The simulation produces a non-empty ``OutputPowerX`` (power-equivalent
  flow into the boundary sector) for the heat zones.

How to run
----------

1. ``pytest tests/integration/test_solve_chp.py``
2. ``python tests/integration/test_solve_chp.py``

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
def test_solve_chp_via_boundary_sector():
    skip_if_no_gams()
    cfg = load_test_config("tiny_bs.yml", "solve_chp")
    out = build_solve(cfg)

    inputs = out["inputs"]
    results = out["results"]

    # The CHP units should be in the units table:
    units = inputs["units"]
    assert any("Diemen" in str(u) for u in units.index)

    # The thermal boundary-sector zones should be present:
    assert "nx" in inputs["sets"]
    bs_zones = set(inputs["sets"]["nx"])
    expected_th = {"Z2_th_1", "Z2_th_2"}
    assert expected_th.intersection(bs_zones), (
        f"Expected at least one of {expected_th} in boundary-sector zones, "
        f"got {bs_zones}"
    )

    # OutputPowerX should be defined when boundary sectors are used:
    assert "OutputPowerX" in results


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
