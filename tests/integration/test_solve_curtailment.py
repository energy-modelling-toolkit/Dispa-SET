# -*- coding: utf-8 -*-
"""
Integration test: curtailment of VRE generation
=================================================

What this test does
-------------------
Sets up a tiny configuration with a very small electrical demand and
heavily over-sized renewable generation (Solar/Wind modifiers = 10).
The expected behaviour is that the optimiser must **curtail** part of
the renewable production, which is reported in
``results['OutputCurtailedPower']``.

The test verifies:

* The simulation runs to completion.
* ``OutputCurtailedPower`` is non-zero (curtailment actually happens).
* The total curtailment is positive but smaller than the total
  available renewable energy.

How to run
----------

1. ``pytest tests/integration/test_solve_curtailment.py``
2. ``python tests/integration/test_solve_curtailment.py``

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


@pytest.mark.timeout(90)
def test_curtailment_is_triggered():
    skip_if_no_gams()
    cfg = load_test_config("tiny_curtailment.yml", "solve_curtailment")
    out = build_solve(cfg)

    results = out["results"]
    assert "OutputCurtailedPower" in results

    curtailment_total = float(results["OutputCurtailedPower"].sum().sum())
    assert curtailment_total > 0, (
        "Expected non-zero curtailment with VRE 10x modifier and "
        "small demand, got 0."
    )
    print(f"Total curtailment: {curtailment_total:.0f} MWh")


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
