# -*- coding: utf-8 -*-
"""
Integration test: build the simulation environment (no GAMS solve)
==================================================================

What this test does
-------------------
Runs ``ds.build_simulation`` on the tiny configuration and verifies that
the standard set of output artefacts is created in the simulation folder:

* ``Inputs.gdx`` - the GDX file passed to GAMS for the dispatch run.
* ``Inputs.p``   - the pickle representation of the inputs (used by
                   ``get_sim_results``).
* ``UCM_h.gms``  - copied GAMS source used by the dispatch model.

This test exercises:

* the YAML config loader,
* the per-zone time-series I/O (Demand / RES AF / etc.),
* the GDX writer (``write_variables``),
* the file copy of the GAMS model.

It does NOT call GAMS, so it works on systems where the GAMS binary is
not installed but the GAMS python API and gdxcc library are available.

How to run
----------

1. ``pytest tests/integration/test_build_only.py``
2. ``python tests/integration/test_build_only.py``
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import load_test_config  # noqa: E402

import dispaset as ds  # noqa: E402


def _check_simdir_artefacts(simdir: Path):
    expected = ["Inputs.gdx", "Inputs.p", "UCM_h.gms"]
    for name in expected:
        f = simdir / name
        assert f.is_file(), f"Missing expected build artefact: {f}"
        assert f.stat().st_size > 0


@pytest.mark.timeout(60)
def test_build_tiny():
    """Build the smallest possible simulation environment (1 zone, 1 day)."""
    cfg = load_test_config("tiny.yml", "build_only_tiny")
    simdir = Path(cfg["SimulationDirectory"])
    ds.build_simulation(cfg)
    _check_simdir_artefacts(simdir)


@pytest.mark.timeout(60)
def test_build_two_zones_with_ntc():
    """Build a two-zone simulation environment with an NTC matrix."""
    cfg = load_test_config("tiny_2zones.yml", "build_only_2zones")
    simdir = Path(cfg["SimulationDirectory"])
    ds.build_simulation(cfg)
    _check_simdir_artefacts(simdir)


@pytest.mark.timeout(120)
def test_build_with_boundary_sector():
    """Build a simulation environment with a hydrogen boundary sector."""
    cfg = load_test_config("tiny_bs.yml", "build_only_bs")
    simdir = Path(cfg["SimulationDirectory"])
    ds.build_simulation(cfg)
    _check_simdir_artefacts(simdir)


@pytest.mark.timeout(60)
def test_build_milp():
    """Build a MILP (Integer-clustering) environment."""
    cfg = load_test_config("tiny_milp.yml", "build_only_milp")
    simdir = Path(cfg["SimulationDirectory"])
    ds.build_simulation(cfg)
    _check_simdir_artefacts(simdir)


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
