# -*- coding: utf-8 -*-
"""
Automates the ODE-SFR -> UCED hand-off for the NC/DC/SC study design: turns
``get_frequency_stability_reserves``'s ``summary_reserves`` output into the
exogenous reserve/inertia CSV files a second (DC or SC) Dispa-SET run reads
as input, replaces a manual copy/rename step with ``clone_config_for_mode``,
and orchestrates the NC -> ODE-SFR -> DC (optionally re-sized DC2) -> compare
pipeline via ``run_uced_odesfr_link``, plus the separate SC benchmark via
``run_uced_sc``.

Pure, file-writing/config-cloning helper functions (``write_odesfr_exogenous_
csvs``, ``compute_static_worst_case``, ``clone_config_for_mode``, ``extract_
contingency``, ``compare_contingency``) are unit-testable without GAMS (see
``tests/unit/test_reserve_link.py``). The two orchestrators (``run_uced_
odesfr_link``, ``run_uced_sc``) take an injected ``build_solve_fn`` callable so
they, too, are unit-testable with a fake solver (no GAMS) -- the real,
GAMS-dependent integration path is exercised in
``tests/integration/test_uced_odesfr_link.py``.

Naming and units
------------------
Dispa-SET's exogenous reserve/inertia inputs are five single-column, hourly
CSVs per run: ``TIMESTAMP,<KEY>_val`` where ``KEY`` is one of ``H`` (system
inertia, MWs, read via ``SystemInertiaDemand``), ``FCR``, ``FFR``, ``aFRR``
(read via both ``aFRRUDemand`` and ``aFRRDDemand`` -- ODE-SFR sizes one
symmetric aFRR value, not separate up/down), and ``VIR``. ``get_frequency_
stability_reserves`` itself already writes CSVs as a side effect (``path +
f"{col}.csv"`` for each of ``Hs_val, Hv_val, SIR_val, VIR_val, FFR_val,
FCR_val, aFRR_val, mFRR_val``, frequency_stability.py:597-600) -- but to the
wrong location (concatenated onto the run's own ``SimulationDirectory``
string, not the ``Reserves``/``InertiaLimit`` subfolders a *different* run's
config points to) and the wrong column name for inertia (``Hs_val``, not the
``H_val`` the ``SystemInertiaDemand`` file format uses). ``write_odesfr_
exogenous_csvs`` below reads the returned ``summary_reserves`` DataFrame
directly instead of relying on that side effect, and writes correctly
named/located files from it.

``H_val = Hs_val`` ONLY, not ``Hs_val + Hv_val`` (DECISION-ODESFR-012 corrects
an earlier, wrong version of this docstring): ``EQ_Inertia_Balance`` in
``UCM.gms`` checks ``InertiaDemand`` against ``sum(au, InertiaProvision(au,i))``,
but ``InertiaProvision`` itself is defined as ``(...)$(cu(au))``
(UCM.gms:1142-1146) -- gated to ``cu``, "Conventional units only"
(UCM.gms:75; ``sets['cu']`` = ``commons['tech_conventional']``, which
excludes ``BATS``). So the sum is syntactically over every unit but is zero
for each non-conventional one: virtual/synthetic inertia from batteries
never enters this balance at all. It is instead procured as the ``VIRU``
reserve product (``commons['ReserveEligibleTechnologies']['VIRU'] ==
commons['tech_batteries']``, common.py:124, with the explicit comment "Real
rotating-mass inertia from synchronous machines is handled separately via
EQ_Inertia_Delivery_Limit/cu(au); VIRU specifically represents *synthetic*
inertia from non-synchronous resources") -- through the ordinary
reserve-balance machinery (``ReserveDemand``/``EQ_UpwardReserves_balance``),
the same path FFR/FCR/aFRR/mFRR use, not through ``SystemInertiaDemand``.
This is exactly what the ``VIR`` entry below (``VIR_val``, ODE-SFR's own
delivered-MW output for this service, mirroring how ``FFR_val``/``FCR_val``
are also delivered-MW quantities -- DECISION-ODESFR-009) is for: ``Hv_val``
(the *Hv_cap* inertia-equivalent search variable ODE-SFR's internal ODE physics
combines with Hs_cap, DECISION-ODESFR-006 -- a statement about the swing-
equation physics, not about which UCED config key carries the result) has no
corresponding exogenous config key and is not written to any file here.

mFRR is intentionally NOT written to its own file: ``build.py`` (preprocessing/
build.py:792) unconditionally sets ``finalTS['mFRRU'] = aFRRUDemand_tot``,
i.e. it overwrites whatever mFRR series would be supplied with the aFRR one
regardless -- no config key reads a distinct mFRR file, so writing one would
be silently ignored (``ISSUE-UCED-014``). ``DEMAND_FILE_SPEC`` below covers
only the 5 keys that are actually wired into a config.
"""
from __future__ import annotations

import copy
import datetime as dt
import json
import logging
import os
from typing import Callable, Dict, List, Optional, Sequence, Union

import pandas as pd

from .frequency_stability import patch_degenerate_reserves

PathLike = Union[str, "os.PathLike[str]"]

# Reserve keys written as their own exogenous CSV -- the summary_reserves
# column each is read from, and the Dispa-SET config key(s) it is wired into.
# ``H`` = Hs_val only (synchronous inertia; see module docstring for why
# Hv_val is deliberately excluded -- it has no corresponding config key).
DEMAND_FILE_SPEC = {
    "H": {"source_columns": ("Hs_val",), "config_keys": ("SystemInertiaDemand",)},
    "VIR": {"source_columns": ("VIR_val",), "config_keys": ("VIRDemand",)},
    "FFR": {"source_columns": ("FFR_val",), "config_keys": ("FFRDemand",)},
    "FCR": {"source_columns": ("FCR_val",), "config_keys": ("FCRDemand",)},
    "aFRR": {"source_columns": ("aFRR_val",), "config_keys": ("aFRRUDemand", "aFRRDDemand")},
}

# Every config key that carries an exogenous reserve/inertia file (used by
# clone_config_for_mode to blank them out for NC).
ALL_DEMAND_CONFIG_KEYS = sorted({k for spec in DEMAND_FILE_SPEC.values() for k in spec["config_keys"]})

# Machine-readable status file stage 'nc' writes into simulation_dir_nc, and
# stage 'dc'/'sc' (run standalone, e.g. as a separate cluster array task) read
# before doing any GAMS work -- see run_stage_nc/read_nc_status/
# _require_secure_nc_status and run_all_simulations.py --stage dc/sc.
NC_STATUS_FILENAME = "ode_sfr_status.json"


# ---------------------------------------------------------------------------
# Pure helpers (no GAMS, no filesystem side effects beyond writing CSVs)
# ---------------------------------------------------------------------------

def compute_static_worst_case(summary_reserves: pd.DataFrame) -> pd.DataFrame:
    """
    The SC ("Static Constraints") series: each column's own annual maximum,
    held constant across every hour -- "static, worst-case-constant sizing"
    Returned frame has the same
    columns and index as ``summary_reserves`` (a broadcast constant per
    column), so it can be passed to ``write_odesfr_exogenous_csvs`` exactly
    like the raw, dynamic ``summary_reserves`` frame.
    """
    # Insecure hours (status False) are excluded: SC secures the most demanding *secure* hour.
    # (run_uced_odesfr_link only builds DC/SC when the screening has no insecure hour.)
    secure = summary_reserves
    if "status" in summary_reserves.columns:
        secure = summary_reserves[summary_reserves["status"].astype(bool)]
    maxima = secure.fillna(0).max()
    return pd.DataFrame({col: maxima[col] for col in summary_reserves.columns}, index=summary_reserves.index)


def config_time_index(config: dict) -> pd.DatetimeIndex:
    """Hourly simulation index of a Dispa-SET config (StartDate..StopDate), as built by build.py (config['idx'])."""
    return pd.date_range(start=dt.datetime(*config["StartDate"]), end=dt.datetime(*config["StopDate"]), freq="h")


def _dc_index(config: dict) -> Optional[pd.DatetimeIndex]:
    """``config_time_index(config)``, or ``None`` when the config has no StartDate/StopDate (test doubles)."""
    if "StartDate" in config and "StopDate" in config:
        return config_time_index(config)
    return None


def _on_index(series: pd.Series, index: pd.DatetimeIndex) -> pd.Series:
    """``series`` re-labelled onto ``index``: by timestamp (nearest, gaps 0) when its index is a datetime index,
    by position when it has exactly ``len(index)`` entries and no usable timestamps."""
    try:
        src = pd.DatetimeIndex(pd.to_datetime(series.index))
        if src.tz is not None:
            src = src.tz_localize(None)
    except (TypeError, ValueError):
        src = None
    if src is None or len(src.intersection(index)) == 0:
        if len(series) != len(index):
            raise ValueError(f"cannot align {len(series)} reserve values on a {len(index)}-step simulation index")
        return pd.Series(series.values, index=index)
    s = pd.Series(series.values, index=src)
    s = s[~s.index.duplicated(keep="first")].sort_index()
    return s.reindex(index, method="nearest").fillna(0)


def write_odesfr_exogenous_csvs(summary_reserves: pd.DataFrame, out_dir: PathLike,
                               case_prefix: str, index: Optional[pd.DatetimeIndex] = None) -> Dict[str, str]:
    """
    Writes the 5 exogenous files a DC or SC config actually reads, from a
    ``summary_reserves``-shaped frame (either ODE-SFR's own dynamic output, or
    the constant frame from ``compute_static_worst_case``).

    Creates ``<out_dir>/Reserves/`` and ``<out_dir>/InertiaLimit/`` (so
    passing a run's own ``SimulationDirectory`` as ``out_dir`` keeps every
    exogenous input specific to that run self-contained inside its own
    folder), and writes ``<case_prefix><KEY>.csv`` into the appropriate one,
    e.g. ``out_dir="Simulations/BAU_2030_DC"``, ``case_prefix="BAU_2030_DC"``
    -> ``Simulations/BAU_2030_DC/Reserves/BAU_2030_DCFCR.csv``.

    :param summary_reserves: DataFrame indexed by TIMESTAMP with at least
                              columns ``Hs_val``, ``FCR_val``, ``FFR_val``,
                              ``aFRR_val``, ``VIR_val`` (``Hv_val`` is not
                              read -- see module docstring).
    :param out_dir:           The run's own simulation directory.
    :param case_prefix:       Filename prefix, e.g. ``"BAU_2030_DC"``.
    :param index:             The DC/SC config's own simulation index (``config_time_index(config)``). Every CSV
                              is written on exactly this index (values aligned by timestamp), so
                              ``load_time_series`` reads it for any horizon, including partial-year runs.
                              ``None``: the index of ``summary_reserves`` is kept.
    :return:                  {config_key: path} for every config key in
                              ``ALL_DEMAND_CONFIG_KEYS`` (``aFRRUDemand`` and
                              ``aFRRDDemand`` both map to the same file).
    """
    reserves_dir = os.path.join(str(out_dir), "Reserves")
    inertia_dir = os.path.join(str(out_dir), "InertiaLimit")
    os.makedirs(reserves_dir, exist_ok=True)
    os.makedirs(inertia_dir, exist_ok=True)

    paths: Dict[str, str] = {}
    for key, spec in DEMAND_FILE_SPEC.items():
        series = sum(summary_reserves[col].fillna(0) for col in spec["source_columns"])
        if index is not None:
            series = _on_index(series, index)
        frame = series.to_frame(name=f"{key}_val")
        frame.index.name = "TIMESTAMP"

        target_dir = inertia_dir if key == "H" else reserves_dir
        path = os.path.join(target_dir, f"{case_prefix}{key}.csv")
        frame.to_csv(path)
        for config_key in spec["config_keys"]:
            paths[config_key] = path

    logging.info("write_odesfr_exogenous_csvs(): wrote %d exogenous reserve/inertia CSVs for '%s' to %s",
                 len(DEMAND_FILE_SPEC), case_prefix, out_dir)
    return paths


def clone_config_for_mode(base_config: dict, mode: str, csv_paths: Optional[Dict[str, str]],
                           simulation_dir: PathLike) -> dict:
    """
    In-memory config clone (no new YAML file written to disk): copies
    ``base_config``, points ``SimulationDirectory`` at ``simulation_dir``,
    and sets the 5 exogenous demand config keys.

    :param mode:            One of ``"NC"``, ``"DC"``, ``"SC"`` -- used only
                             for validation/logging, not to branch behavior
                             beyond the ``csv_paths is None`` check below.
    :param csv_paths:       The dict returned by ``write_odesfr_exogenous_csvs``,
                             or ``None`` for NC. When ``None``, every demand
                             config key is set to ``''`` -- Dispa-SET's own
                             ``build.py`` treats a missing/empty path as an
                             all-zero requirement (``os.path.isfile('') ==
                             False``), which is how "No Constraints" is
                             actually implemented: a real all-zero demand,
                             not a special code path.
    :param simulation_dir:  This run's own ``SimulationDirectory``.
    """
    if mode not in ("NC", "DC", "SC"):
        raise ValueError(f"mode must be one of 'NC', 'DC', 'SC', got {mode!r}")
    if mode == "NC" and csv_paths is not None:
        raise ValueError("csv_paths must be None for mode='NC' (NC uses an all-zero, "
                          "empty-path exogenous requirement, not a written CSV)")
    if mode != "NC" and csv_paths is None:
        raise ValueError(f"csv_paths is required for mode={mode!r}")

    config = copy.deepcopy(base_config)
    config["SimulationDirectory"] = str(simulation_dir)

    if csv_paths is None:
        for key in ALL_DEMAND_CONFIG_KEYS:
            config[key] = ""
    else:
        for key in ALL_DEMAND_CONFIG_KEYS:
            config[key] = csv_paths[key]

    return config


# ---------------------------------------------------------------------------
# Stage-split support (Array 1 / Array 2 cluster design): a machine-readable
# status file stage 'nc' writes, and a pure re-derivation of the exogenous CSV
# paths stage 'nc' already wrote -- so stage 'dc'/'sc' can run as a separate
# process/cluster task without re-running NC or ODE-SFR. See run_stage_nc/
# run_stage_dc/run_stage_sc below and run_all_simulations.py --stage.
# ---------------------------------------------------------------------------

def nc_status_path(simulation_dir_nc: PathLike) -> str:
    """Path of the status JSON ``run_stage_nc`` writes into ``simulation_dir_nc``."""
    return os.path.join(str(simulation_dir_nc), NC_STATUS_FILENAME)


def write_nc_status(simulation_dir_nc: PathLike, *, secure: bool, n_insecure: int,
                     n_hours: int, base_case_name: str) -> str:
    """
    Writes the status JSON a standalone stage 'dc'/'sc' task reads to decide
    whether to run at all: ``{"secure": bool, "n_insecure": int, "n_hours": int,
    "base_case_name": str, "written_at": iso-timestamp}``. ``secure`` is the
    single source of truth stage 'dc'/'sc' act on (``n_insecure == 0``).
    """
    os.makedirs(str(simulation_dir_nc), exist_ok=True)
    payload = {
        "base_case_name": base_case_name,
        "secure": bool(secure),
        "n_insecure": int(n_insecure),
        "n_hours": int(n_hours),
        "written_at": dt.datetime.now().isoformat(timespec="seconds"),
    }
    path = nc_status_path(simulation_dir_nc)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    logging.info("write_nc_status(): wrote %s (secure=%s, %d/%d insecure hour(s))",
                 path, payload["secure"], n_insecure, n_hours)
    return path


def read_nc_status(simulation_dir_nc: PathLike) -> Optional[dict]:
    """Reads back ``write_nc_status``'s JSON, or ``None`` if stage 'nc' has not
    run (or not finished) for this directory yet -- never raises for a missing
    file, since "NC not run yet" is an ordinary, expected state for a cluster
    stage-2 task whose stage-1 dependency has not completed."""
    path = nc_status_path(simulation_dir_nc)
    if not os.path.isfile(path):
        return None
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _require_secure_nc_status(simulation_dir_nc: PathLike) -> dict:
    """Used by ``run_stage_dc``/``run_stage_sc``: raises ``RuntimeError`` with a
    clear, specific message (missing vs. insecure) if stage 'dc'/'sc' should
    not run -- the caller (``run_all_simulations.py --stage dc/sc``) is
    expected to check ``read_nc_status`` itself first and exit 0 cleanly
    instead of letting this raise; this is the defence-in-depth guard for any
    other caller (e.g. interactive/Spyder use)."""
    status = read_nc_status(simulation_dir_nc)
    if status is None:
        raise RuntimeError(
            f"no NC status file at {nc_status_path(simulation_dir_nc)!r} -- run stage 'nc' first")
    if not status.get("secure"):
        raise RuntimeError(
            f"NC stage reported {status.get('n_insecure')} insecure hour(s) of "
            f"{status.get('n_hours')} for {status.get('base_case_name')!r} -- "
            "stage 'dc'/'sc' do not run for an insecure fleet (battery storage sizing required first)")
    return status


def _require_nc_status_exists(simulation_dir_nc: PathLike) -> dict:
    """Lenient counterpart to ``_require_secure_nc_status``, used only when a
    caller passes ``require_secure=False`` (TESTING ONLY -- see
    ``run_stage_dc``/``run_stage_sc``'s ``require_secure`` parameter): still
    requires stage 'nc' to have run at all (a status file to exist), but does
    not require ``secure``. Logs a warning either way so an insecure-but-
    forced run is never silent."""
    status = read_nc_status(simulation_dir_nc)
    if status is None:
        raise RuntimeError(
            f"no NC status file at {nc_status_path(simulation_dir_nc)!r} -- run stage 'nc' first")
    if not status.get("secure"):
        logging.warning("require_secure=False (TESTING ONLY): proceeding despite NC reporting %s insecure "
                         "hour(s) of %s for %r", status.get("n_insecure"), status.get("n_hours"),
                         status.get("base_case_name"))
    return status


def _demand_csv_paths(out_dir: PathLike, case_prefix: str, require_exists: bool = True) -> Dict[str, str]:
    """
    Re-derives the ``{config_key: path}`` mapping ``write_odesfr_exogenous_csvs``
    would return for the same ``(out_dir, case_prefix)`` -- pure string
    construction from ``DEMAND_FILE_SPEC``, no filesystem write. Used by
    ``run_stage_dc``/``run_stage_sc`` (run standalone, e.g. as a separate
    cluster array task) to locate the CSVs ``run_stage_nc`` already wrote,
    without re-running ODE-SFR. The two path-construction rules (here and in
    ``write_odesfr_exogenous_csvs``) must stay in sync -- a unit test
    (``test_demand_csv_paths_matches_write_odesfr_exogenous_csvs``) guards this.

    :param require_exists: Raises ``FileNotFoundError`` (fail-fast, not a
                            silently-wrong all-zero config) if a file is
                            missing. ``False`` is used only for the
                            cross-check unit test above.
    """
    reserves_dir = os.path.join(str(out_dir), "Reserves")
    inertia_dir = os.path.join(str(out_dir), "InertiaLimit")
    paths: Dict[str, str] = {}
    for key, spec in DEMAND_FILE_SPEC.items():
        target_dir = inertia_dir if key == "H" else reserves_dir
        path = os.path.join(target_dir, f"{case_prefix}{key}.csv")
        if require_exists and not os.path.isfile(path):
            raise FileNotFoundError(
                f"expected exogenous CSV not found: {path} (run stage 'nc' first, "
                "or check base_case_name/out_dir)")
        for config_key in spec["config_keys"]:
            paths[config_key] = path
    return paths


def extract_contingency(results: dict) -> pd.Series:
    """
    Thin, documented wrapper around ``results['OutputContingency']`` (the
    system-level N-1 contingency series a run actually solved with). ODE-SFR
    itself never consumes ``OutputContingencyPerZone`` (see
    ``frequency_stability.py``), so comparing per-zone contingency would
    test something ODE-SFR does not actually use -- an intentional non-goal,
    not an oversight.
    """
    return results["OutputContingency"]


def compare_contingency(series_a: pd.Series, series_b: pd.Series, atol: float = 1.0) -> dict:
    """
    Compares two hourly contingency series (e.g. NC vs. DC1) on their common
    timestamps.

    :param atol: Absolute tolerance in MW below which two hours' contingency
                 are considered "the same" (default 1 MW).
    :return: ``{"n_hours", "n_discrepant", "pct_discrepant", "max_abs_diff",
              "discrepant_hours"}`` -- ``discrepant_hours`` is a list of the
              timestamps where ``|series_a - series_b| > atol``.
    """
    a, b = series_a.align(series_b, join="inner")
    diff = (a - b).abs()
    discrepant = diff[diff > atol]
    n_hours = len(diff)
    n_discrepant = len(discrepant)
    return {
        "n_hours": n_hours,
        "n_discrepant": n_discrepant,
        "pct_discrepant": (100.0 * n_discrepant / n_hours) if n_hours else float("nan"),
        "max_abs_diff": float(diff.max()) if n_hours else float("nan"),
        "discrepant_hours": list(discrepant.index),
    }


# ---------------------------------------------------------------------------
# Stage functions (GAMS-dependent via injected build_solve_fn) -- the
# building blocks of the NC -> ODE-SFR -> DC/SC pipeline, usable either
# chained in-process (run_uced_odesfr_link, "--stage all") or standalone, one
# per cluster array task ("--stage nc/dc/sc" in run_all_simulations.py: Array
# 1 = run_stage_nc per scenario-year; Array 2 = run_stage_dc and run_stage_sc,
# one SLURM array each, started via --dependency=aftercorr: once the matching
# Array-1 task finishes).
# ---------------------------------------------------------------------------

def run_stage_nc(base_config: dict,
                  build_solve_fn: Callable[[dict], "tuple[dict, dict]"],
                  get_frequency_stability_reserves_fn: Callable[..., "tuple"],
                  simulation_dir_nc: PathLike,
                  simulation_dir_dc: PathLike,
                  simulation_dir_sc: Optional[PathLike],
                  base_case_name: str,
                  odesfr_kwargs: Optional[dict] = None,
                  force_write_csvs: bool = False) -> dict:
    """
    Stage 'nc': solve NC, size ODE-SFR on it (round 1 only -- see note below),
    write DC1's and, if ``simulation_dir_sc`` is given, SC's exogenous CSVs
    from that same round-1 result, and write the machine-readable status JSON
    (``write_nc_status``) stage 'dc'/'sc' read before doing any GAMS work.

    Only round 1 is computed here: the split-stage cluster design lets DC and
    SC start as soon as this stage finishes (``--dependency=aftercorr``), so
    by construction neither can depend on a DC1 result first. The optional
    DC1-based re-sizing round (``rounds=2``) stays available only through the
    sequential, non-split ``run_uced_odesfr_link``/``--stage all`` path, which
    calls this function for round 1 and then runs round 2 inline -- see its
    docstring.

    :param simulation_dir_dc: DC1's own simulation directory (where its
                               exogenous CSVs are written -- NOT solved here).
    :param simulation_dir_sc: SC's own simulation directory, or ``None`` to
                               skip writing SC's CSVs here (e.g. the ``--stage
                               all`` path, which writes them later itself via
                               ``run_uced_sc`` using the *final* round's
                               summary_reserves -- round 2's if it ran).
    :param force_write_csvs:  TESTING ONLY (``run_all_simulations.py
                               --force-continue``, never set by the cluster
                               job files) -- writes DC1's/SC's exogenous CSVs
                               even when ``secure`` is ``False``, so the
                               mechanical DC/SC pipeline can be exercised on a
                               short/insecure test horizon without a real,
                               fully-secure BESS fleet. ``secure``/
                               ``n_insecure`` and the status JSON still report
                               the true screening result either way -- this
                               never hides an insecure outcome, it only lets a
                               caller that already knows that and wants to
                               test anyway bypass the normal "don't build
                               DC/SC for an insecure fleet" gate.
    :return: dict with keys ``nc_inputs``, ``nc_results``, ``nc_contingency``,
             ``summary_reserves``, ``patch_log``, ``n_insecure``, ``secure``,
             ``csv_paths_dc`` and ``csv_paths_sc`` (``None`` if not written).
    """
    odesfr_kwargs = odesfr_kwargs or {}

    nc_config = clone_config_for_mode(base_config, "NC", None, simulation_dir_nc)
    nc_inputs, nc_results = build_solve_fn(nc_config)
    nc_contingency = extract_contingency(nc_results)

    _, summary_reserves, _ = get_frequency_stability_reserves_fn(
        str(simulation_dir_nc), nc_inputs, nc_results, **odesfr_kwargs)
    # See run_uced_odesfr_link's matching comment: repairs search-artifact
    # degenerate rows a full-bound, jointly-delivered re-check finds feasible.
    summary_reserves, patch_log = patch_degenerate_reserves(summary_reserves)

    n_insecure = 0
    if "status" in summary_reserves.columns:
        n_insecure = int((~summary_reserves["status"].astype(bool)).sum())
    secure = n_insecure == 0

    if not secure and force_write_csvs:
        logging.warning("run_stage_nc(): %s has %d insecure hour(s) but force_write_csvs=True (TESTING ONLY) -- "
                         "writing DC1/SC exogenous CSVs anyway; status JSON still records secure=False",
                         base_case_name, n_insecure)

    csv_paths_dc = None
    csv_paths_sc = None
    if secure or force_write_csvs:
        csv_paths_dc = write_odesfr_exogenous_csvs(
            summary_reserves, simulation_dir_dc, f"{base_case_name}_DC1", index=_dc_index(base_config))
        if simulation_dir_sc is not None:
            static_frame = compute_static_worst_case(summary_reserves)
            csv_paths_sc = write_odesfr_exogenous_csvs(
                static_frame, simulation_dir_sc, f"{base_case_name}_SC", index=_dc_index(base_config))

    write_nc_status(simulation_dir_nc, secure=secure, n_insecure=n_insecure,
                     n_hours=len(summary_reserves), base_case_name=base_case_name)

    return {
        "nc_inputs": nc_inputs,
        "nc_results": nc_results,
        "nc_contingency": nc_contingency,
        "summary_reserves": summary_reserves,
        "patch_log": patch_log,
        "n_insecure": n_insecure,
        "secure": secure,
        "csv_paths_dc": csv_paths_dc,
        "csv_paths_sc": csv_paths_sc,
    }


def run_stage_dc(base_config: dict,
                  build_solve_fn: Callable[[dict], "tuple[dict, dict]"],
                  simulation_dir_nc: PathLike,
                  simulation_dir_dc: PathLike,
                  base_case_name: str,
                  atol: float = 1.0,
                  nc_results: Optional[dict] = None,
                  load_nc_results_fn: Optional[Callable[[str], "tuple[dict, dict]"]] = None,
                  require_secure: bool = True) -> dict:
    """
    Stage 'dc': requires stage 'nc' to have already run and reported secure
    (checked via ``_require_secure_nc_status`` -- raises otherwise; the
    standalone cluster caller, ``run_all_simulations.py --stage dc``, is
    expected to check ``read_nc_status`` itself first and exit 0 cleanly
    instead of hitting this exception). Re-derives DC1's exogenous CSV paths
    (``_demand_csv_paths`` -- the files stage 'nc' already wrote), solves DC1,
    and produces the NC-vs-DC1 contingency comparison report.

    :param nc_results:        NC's own (already in memory) results, to avoid
                               re-reading its GDX when this runs right after
                               stage 'nc' in the same process (the ``--stage
                               all`` path). ``None`` (the standalone/cluster
                               ``--stage dc`` path) loads them from
                               ``simulation_dir_nc`` via ``load_nc_results_fn``.
    :param load_nc_results_fn: ``path -> (inputs, results)``, defaults to
                               ``dispaset.postprocessing.data_handler.
                               get_sim_results``. Injected so this is
                               unit-testable without GAMS/a real GDX file.
    :param require_secure:    TESTING ONLY (``run_all_simulations.py
                               --force-continue``, never set by the cluster
                               job files) -- ``False`` still requires stage
                               'nc' to have run (a status file to exist) but
                               does not require it to have reported secure,
                               so the mechanical DC pipeline can be exercised
                               against CSVs stage 'nc' wrote with its own
                               ``force_write_csvs=True``.
    :return: dict with keys ``config``, ``inputs``, ``results``,
             ``contingency``, ``comparison`` (``compare_contingency``'s dict
             vs. NC), ``converged`` (``comparison["n_discrepant"] == 0``).
    """
    if require_secure:
        _require_secure_nc_status(simulation_dir_nc)
    else:
        _require_nc_status_exists(simulation_dir_nc)

    try:
        csv_paths = _demand_csv_paths(simulation_dir_dc, f"{base_case_name}_DC1")
    except FileNotFoundError:
        # DC1's CSVs are missing when the NC run came from --stage all and DC never started: rebuild
        # them from NC's round-1 (patched) ODE-SFR summary, exactly as run_stage_nc writes them.
        summary_path = os.path.join(str(simulation_dir_nc), "ODESFR_summary_reserves.csv")
        logging.warning("%s: DC1 exogenous CSVs not found -- rebuilding them from %s", base_case_name, summary_path)
        summary = pd.read_csv(summary_path, index_col=0, parse_dates=True)
        csv_paths = write_odesfr_exogenous_csvs(summary, simulation_dir_dc, f"{base_case_name}_DC1",
                                                index=_dc_index(base_config))
    dc_config = clone_config_for_mode(base_config, "DC", csv_paths, simulation_dir_dc)
    dc_inputs, dc_results = build_solve_fn(dc_config)
    dc_contingency = extract_contingency(dc_results)

    if nc_results is None:
        load_nc_results_fn = load_nc_results_fn or _default_load_nc_results_fn()
        _, nc_results = load_nc_results_fn(str(simulation_dir_nc))
    nc_contingency = extract_contingency(nc_results)

    comparison = compare_contingency(nc_contingency, dc_contingency, atol=atol)

    return {
        "config": dc_config,
        "inputs": dc_inputs,
        "results": dc_results,
        "contingency": dc_contingency,
        "comparison": comparison,
        "converged": comparison["n_discrepant"] == 0,
    }


def run_stage_sc(base_config: dict,
                  build_solve_fn: Callable[[dict], "tuple[dict, dict]"],
                  simulation_dir_nc: PathLike,
                  simulation_dir_sc: PathLike,
                  base_case_name: str,
                  require_secure: bool = True) -> dict:
    """
    Stage 'sc': requires stage 'nc' to have already run and reported secure
    (same guard as ``run_stage_dc``). Re-derives SC's exogenous CSV paths
    (the static, worst-case-constant files ``run_stage_nc`` wrote when called
    with a ``simulation_dir_sc``) and solves SC. Independent of stage 'dc' --
    does not read ``simulation_dir_dc`` or anything DC1 produced, only
    ``simulation_dir_nc``'s status and ``simulation_dir_sc``'s own CSVs.

    :param require_secure: TESTING ONLY -- see ``run_stage_dc``'s parameter
                            of the same name.
    :return: ``{"config", "inputs", "results"}`` for the solved SC run (same
             shape as ``run_uced_sc``).
    """
    if require_secure:
        _require_secure_nc_status(simulation_dir_nc)
    else:
        _require_nc_status_exists(simulation_dir_nc)

    try:
        csv_paths = _demand_csv_paths(simulation_dir_sc, f"{base_case_name}_SC")
    except FileNotFoundError:
        # SC's CSVs are missing when the NC run came from --stage all and SC never started (e.g. the
        # task was stopped during DC): rebuild them from NC's round-1 ODE-SFR summary, as run_uced_sc does.
        summary_path = os.path.join(str(simulation_dir_nc), "ODESFR_summary_reserves.csv")
        logging.warning("%s: SC exogenous CSVs not found -- rebuilding them from %s", base_case_name, summary_path)
        summary = pd.read_csv(summary_path, index_col=0, parse_dates=True)
        csv_paths = write_odesfr_exogenous_csvs(compute_static_worst_case(summary), simulation_dir_sc,
                                                f"{base_case_name}_SC", index=_dc_index(base_config))
    sc_config = clone_config_for_mode(base_config, "SC", csv_paths, simulation_dir_sc)
    inputs, results = build_solve_fn(sc_config)
    return {"config": sc_config, "inputs": inputs, "results": results}


def _default_load_nc_results_fn() -> Callable[[str], "tuple[dict, dict]"]:
    """Lazy import of ``get_sim_results`` (avoids a hard, module-level
    dependency on ``data_handler`` for the pure/no-GAMS parts of this module,
    and keeps ``run_stage_dc`` trivially fake-able in unit tests)."""
    from .data_handler import get_sim_results

    def _load(path: str):
        return get_sim_results(path=path, cache=False)

    return _load


# ---------------------------------------------------------------------------
# Orchestrators (GAMS-dependent via injected build_solve_fn)
# ---------------------------------------------------------------------------

def run_uced_odesfr_link(base_config: dict,
                        build_solve_fn: Callable[[dict], "tuple[dict, dict]"],
                        get_frequency_stability_reserves_fn: Callable[..., "tuple"],
                        simulation_dir_nc: PathLike,
                        simulation_dirs_dc: Sequence[PathLike],
                        base_case_name: str,
                        atol: float = 1.0,
                        odesfr_kwargs: Optional[dict] = None) -> dict:
    """
    Orchestrates NC -> ODE-SFR -> DC1 [-> re-size ODE-SFR -> DC2] -> compare.

    1. Solve NC (``build_solve_fn(nc_config)``; NC's exogenous demand is
       all-zero via ``clone_config_for_mode(..., 'NC', None, ...)``).
    2. Size ODE-SFR on NC's own (inputs, results) -> ``summary_reserves_round1``.
    3. Write DC1's exogenous CSVs, solve DC1.
    4. Compare NC vs. DC1's system-level contingency (``round_1_comparison``).
    5. If round 1 is discrepant AND a second directory is given in
       ``simulation_dirs_dc``: re-size ODE-SFR on DC1's own (inputs, results) ->
       ``summary_reserves_round2``, write DC2, solve DC2, compare DC1 vs.
       DC2 (``round_2_comparison``). Otherwise stop after round 1 --
       ``round_2_comparison`` stays ``None`` and ``rounds_run`` stays 1; this
       is never silently assumed to be convergence (``converged`` reports
       exactly what the last comparison found).

    :param build_solve_fn:   ``config -> (inputs, results)``, e.g.
                              ``build_simulation`` + ``solve_GAMS`` +
                              ``get_sim_results`` composed by the caller --
                              injected so this function is unit-testable
                              with a fake solver (no GAMS required).
    :param get_frequency_stability_reserves_fn: Normally
                              ``dispaset.get_frequency_stability_reserves``,
                              injected for the same reason.
    :param simulation_dirs_dc: One directory per round to attempt (length 1
                              -> single pass, no re-sizing; length 2 -> a
                              second round is attempted if round 1 disagrees).
    :param odesfr_kwargs:       Extra kwargs forwarded to
                              ``get_frequency_stability_reserves_fn``
                              (e.g. custom ``limit_freq``).
    :return: dict with keys ``rounds_run``, ``converged``,
             ``round_1_comparison``, ``round_2_comparison`` (``None`` if not
             run), ``summary_reserves_round1``, ``summary_reserves_round2``
             (``None`` if not run), plus the resolved directories/names used.
    """
    odesfr_kwargs = odesfr_kwargs or {}
    simulation_dirs_dc = list(simulation_dirs_dc)
    if not simulation_dirs_dc:
        raise ValueError("simulation_dirs_dc must have at least one entry")
    dir_dc1 = simulation_dirs_dc[0]

    report: dict = {
        "simulation_dir_nc": str(simulation_dir_nc),
        "simulation_dirs_dc": [str(d) for d in simulation_dirs_dc],
        "base_case_name": base_case_name,
        "rounds_run": 0,
        "converged": False,
        "round_1_comparison": None,
        "round_2_comparison": None,
        "summary_reserves_round1": None,
        "summary_reserves_round2": None,
        "patch_log_round1": None,
        "patch_log_round2": None,
        "n_insecure_round1": None,
        "stopped_insecure": False,
    }

    # --- Stage 'nc': NC + round-1 ODE-SFR sizing (shared with --stage nc) ---
    # simulation_dir_sc=None: SC's CSVs are written later, from the *final*
    # round's summary_reserves (round 2's if it runs below, else round 1's --
    # see run_uced_sc, called separately by this function's own callers, e.g.
    # run_all_simulations.py's --stage all path). This keeps this function's
    # numeric results byte-identical to before this was factored into stages.
    nc_stage = run_stage_nc(base_config, build_solve_fn, get_frequency_stability_reserves_fn,
                             simulation_dir_nc, dir_dc1, None, base_case_name,
                             odesfr_kwargs=odesfr_kwargs)
    report["patch_log_round1"] = nc_stage["patch_log"]
    report["summary_reserves_round1"] = nc_stage["summary_reserves"]
    report["n_insecure_round1"] = nc_stage["n_insecure"]

    # Frequency-adequacy screening: DC (and SC) are only built for a fleet that secures every hour.
    # With insecure hours left, the fleet needs battery storage sizing first; NC + ODE-SFR are rerun
    # with the battery fleet before any frequency-constrained run.
    if not nc_stage["secure"]:
        report["stopped_insecure"] = True
        return report

    # --- Stage 'dc': DC1 solved against round-1's requirements + NC-vs-DC1 compare (shared with --stage dc) ---
    dc_stage = run_stage_dc(base_config, build_solve_fn, simulation_dir_nc, dir_dc1, base_case_name,
                             atol=atol, nc_results=nc_stage["nc_results"])
    report["rounds_run"] = 1
    report["round_1_comparison"] = dc_stage["comparison"]

    if dc_stage["comparison"]["n_discrepant"] == 0 or len(simulation_dirs_dc) < 2:
        report["converged"] = dc_stage["comparison"]["n_discrepant"] == 0
        return report

    # --- Round 2: re-size ODE-SFR on DC1, solved as DC2 -- NOT part of the
    # split-stage cluster design (Array 1/2 only ever run round 1; see
    # run_stage_nc's docstring), so kept inline here rather than as its own
    # stage function. ---
    _, summary_reserves_2, _ = get_frequency_stability_reserves_fn(
        str(dir_dc1), dc_stage["inputs"], dc_stage["results"], **odesfr_kwargs)
    summary_reserves_2, patch_log_2 = patch_degenerate_reserves(summary_reserves_2)
    report["patch_log_round2"] = patch_log_2
    report["summary_reserves_round2"] = summary_reserves_2

    dir_dc2 = simulation_dirs_dc[1]
    csv_paths_2 = write_odesfr_exogenous_csvs(summary_reserves_2, dir_dc2, f"{base_case_name}_DC2",
                                              index=_dc_index(base_config))
    dc2_config = clone_config_for_mode(base_config, "DC", csv_paths_2, dir_dc2)
    dc2_inputs, dc2_results = build_solve_fn(dc2_config)
    report["rounds_run"] = 2

    dc2_contingency = extract_contingency(dc2_results)
    round_2_comparison = compare_contingency(dc_stage["contingency"], dc2_contingency, atol=atol)
    report["round_2_comparison"] = round_2_comparison
    report["converged"] = round_2_comparison["n_discrepant"] == 0

    return report


def run_uced_sc(base_config: dict,
                 summary_reserves_for_worst_case: pd.DataFrame,
                 build_solve_fn: Callable[[dict], "tuple[dict, dict]"],
                 simulation_dir_sc: PathLike,
                 base_case_name: str) -> dict:
    """
    Builds and solves the SC ("Static Constraints") run: the static,
    worst-case-constant requirement derived from ``summary_reserves_for_
    worst_case`` (by convention, the final round's ``summary_reserves`` from
    ``run_uced_odesfr_link`` -- round 2's if it ran, otherwise round 1's).
    Separate from the NC/DC convergence loop, since SC is a fixed benchmark
    by construction, not something that gets iteratively refined.

    :return: ``{"config", "inputs", "results"}`` for the solved SC run.
    """
    static_frame = compute_static_worst_case(summary_reserves_for_worst_case)
    csv_paths = write_odesfr_exogenous_csvs(static_frame, simulation_dir_sc, f"{base_case_name}_SC",
                                            index=_dc_index(base_config))
    sc_config = clone_config_for_mode(base_config, "SC", csv_paths, simulation_dir_sc)
    inputs, results = build_solve_fn(sc_config)
    return {"config": sc_config, "inputs": inputs, "results": results}
