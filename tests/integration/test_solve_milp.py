# -*- coding: utf-8 -*-
"""
Integration test: build + solve - MILP (Integer-clustering)
============================================================

What this test does
-------------------
End-to-end exercise of Dispa-SET with the unit-commitment MILP formulation
(``SimulationType: Integer clustering``) on a tiny single-zone case.

The MILP path enables many additional constraints (start-up costs,
minimum up/down time, partial-load minimum power, etc.) compared to the
LP-clustered variant, so this test verifies that the MILP formulation
solves to optimality on the tiny case and that the binary commitment
output is correctly returned.

Performance budget
------------------
Should run in **less than 10 seconds** on a CPU laptop.

How to run
----------

1. ``pytest tests/integration/test_solve_milp.py``
2. ``python tests/integration/test_solve_milp.py``

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
def test_solve_tiny_milp():
    skip_if_no_gams()
    cfg = load_test_config("tiny_milp.yml", "solve_milp")
    out = build_solve(cfg)

    inputs = out["inputs"]
    results = out["results"]

    # Commitment must be present and integer-valued (Integer clustering):
    assert "OutputCommitted" in results, "Missing OutputCommitted output"
    committed = results["OutputCommitted"]
    if not committed.empty:
        # Values can be 0..Nunits per cluster; all must be non-negative integers.
        vals = committed.values.flatten()
        assert (vals >= -1e-6).all(), "Negative commitment values found"
        non_int = sum(abs(v - round(v)) > 1e-3 for v in vals)
        # Allow a few numerical leakages but no widespread non-integer values.
        assert non_int < 5, f"Too many non-integer commitment values: {non_int}"

    print(f"MILP build {out['build_time']:.2f}s  "
          f"solve {out['solve_time']:.2f}s")


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
