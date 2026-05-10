# -*- coding: utf-8 -*-
"""
Helpers shared by every Dispa-SET test.

The functions in this module avoid duplicating boilerplate in the individual
test files. Every test file imports its helpers from here.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import matplotlib

# Force a non-interactive backend so ``plt.show`` becomes a no-op everywhere.
matplotlib.use("Agg")

import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
TESTS_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = TESTS_DIR / "output"

# Make sure dispaset is importable when the test file is run as a script.
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def gams_available() -> bool:
    """True if a GAMS installation can be located."""
    from dispaset.misc.gdx_handler import package_exists, get_gams_path
    if not package_exists("gams"):
        return False
    try:
        return get_gams_path() is not None
    except Exception:
        return False


def skip_if_no_gams():
    """Pytest helper - skip the running test if no GAMS is detected."""
    if not gams_available():
        pytest.skip("GAMS python API or installation not detected")


def fresh_simdir(suffix: str) -> Path:
    """Return an empty ``tests/output/<suffix>`` folder (creating it)."""
    p = OUTPUT_DIR / suffix
    if p.exists():
        shutil.rmtree(p)
    p.mkdir(parents=True, exist_ok=True)
    return p


def load_test_config(config_name: str, simdir_suffix: str | None = None):
    """
    Load a YAML configuration from ``tests/configs/`` and assign a fresh
    SimulationDirectory under ``tests/output/`` so concurrent tests do not
    collide.
    """
    import dispaset as ds
    config_path = TESTS_DIR / "configs" / config_name
    cfg = ds.load_config(str(config_path))
    if simdir_suffix is None:
        simdir_suffix = Path(cfg["SimulationDirectory"]).name
    simdir = fresh_simdir(simdir_suffix)
    cfg["SimulationDirectory"] = str(simdir)
    return cfg


def solve_succeeded(simdir: str | Path) -> bool:
    """Return True if a Results.gdx was created (regardless of return code)."""
    return (Path(simdir) / "Results.gdx").is_file()


def build_solve(config: dict, *, mts: bool = False) -> dict:
    """
    Build the simulation environment, run GAMS, return the inputs and results.

    Returns a dict with keys ``simdir``, ``inputs``, ``results``,
    ``build_time``, ``solve_time``, ``read_time``.
    """
    import time
    import dispaset as ds

    out: dict = {}
    out["simdir"] = config["SimulationDirectory"]

    t0 = time.time()
    ds.build_simulation(config)
    out["build_time"] = time.time() - t0

    t0 = time.time()
    if mts:
        ds.solve_GAMS(config["SimulationDirectory"], config["GAMS_folder"],
                      gams_file="UCM_MTS.gms",
                      result_file="Results_MTS.gdx")
    else:
        ds.solve_GAMS(config["SimulationDirectory"], config["GAMS_folder"])
    out["solve_time"] = time.time() - t0

    if not solve_succeeded(config["SimulationDirectory"]):
        raise RuntimeError(
            f"GAMS did not produce Results.gdx in {config['SimulationDirectory']}"
        )

    t0 = time.time()
    inputs, results = ds.get_sim_results(path=config["SimulationDirectory"])
    out["read_time"] = time.time() - t0
    out["inputs"] = inputs
    out["results"] = results

    # Systematic energy-balance check for every simulated zone.
    # check_energy_balance logs CRITICAL for any zone with > 1 % imbalance,
    # which will be caught by the fail_on_critical_log fixture in
    # tests/integration/ and tests/ultimate/.
    out["energy_balance"] = ds.check_energy_balance(inputs, results)

    return out
