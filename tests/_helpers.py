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


# Lost-load variables produced by the GAMS model.  These capture real
# violations of physical constraints (load shedding, reserve shortfall,
# ramp violations).  StorageLevelViolation is kept separate because it is
# a soft numerical penalty, not a dispatch infeasibility.
_LOSTLOAD_KEYS = [
    'LostLoad_MaxPower',
    'LostLoad_MinPower',
    'LostLoad_RampUp',
    'LostLoad_RampDown',
    'LostLoad_aFRRU',
    'LostLoad_aFRRD',
    'LostLoad_mFRRU',
    'LostLoad_Inertia',
    'LostLoad_FFRU',
    'LostLoad_FFRD',
    'LostLoad_FCRU',
    'LostLoad_FCRD',
]


def _sum_result(val) -> float:
    """Return the scalar total of a results value (DataFrame, Series or scalar)."""
    if val is None:
        return 0.0
    try:
        return float(val.sum().sum())
    except AttributeError:
        try:
            return float(val.sum())
        except (AttributeError, TypeError):
            return float(val)


def get_solver_kpis(results: dict) -> dict:
    """Return a compact dict of solution-quality KPIs.

    Keys
    ----
    total_cost : float
        Sum of OutputSystemCost (all hours × zones).
    total_lostload_mwh : float
        Sum of all LostLoad_* variables (excluding storage violation).
    storage_violation_mwh : float
        Sum of LostLoad_StorageLevelViolation.
    """
    kpis: dict = {}

    sc = results.get("OutputSystemCost")
    if sc is not None:
        kpis["total_cost"] = _sum_result(sc)

    kpis["total_lostload_mwh"] = sum(
        max(0.0, _sum_result(results.get(k))) for k in _LOSTLOAD_KEYS
    )

    sv = results.get("LostLoad_StorageLevelViolation")
    kpis["storage_violation_mwh"] = max(0.0, _sum_result(sv))

    return kpis


def assert_feasible(results: dict, max_lostload_mwh: float = 1.0) -> dict:
    """Assert that a solved simulation has no significant load shedding.

    Parameters
    ----------
    results : dict
        As returned by ``ds.get_sim_results``.
    max_lostload_mwh : float
        Tolerance in MWh.  Use 0.0 for a strict no-lostload assertion.

    Returns
    -------
    dict
        The KPI dict from :func:`get_solver_kpis` for optional further checks.
    """
    import math

    kpis = get_solver_kpis(results)

    total_ll = kpis.get("total_lostload_mwh", 0.0)
    assert total_ll <= max_lostload_mwh, (
        f"Lost load {total_ll:.3f} MWh exceeds threshold {max_lostload_mwh:.3f} MWh"
    )

    cost = kpis.get("total_cost")
    if cost is not None:
        assert cost >= -1e-3, f"Negative total system cost: {cost:.2f}"
        assert math.isfinite(cost), f"Non-finite total system cost: {cost}"

    return kpis


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

    # Solver KPIs: objective value and lost-load totals.
    out["kpis"] = get_solver_kpis(results)

    return out
