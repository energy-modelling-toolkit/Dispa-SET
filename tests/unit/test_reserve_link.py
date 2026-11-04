# -*- coding: utf-8 -*-
"""
Unit tests for dispaset.postprocessing.reserve_link -- the ODE-SFR -> UCED
exogenous-CSV hand-off and the NC/DC/SC orchestration. No GAMS dependency:
all fixtures are synthetic summary_reserves-shaped DataFrames, and
run_uced_odesfr_link/run_uced_sc are exercised with a fake build_solve_fn /
get_frequency_stability_reserves_fn, matching the pattern used in
test_odesfr.py.
"""
import pandas as pd
import pytest

from dispaset.postprocessing.reserve_link import (
    DEMAND_FILE_SPEC,
    ALL_DEMAND_CONFIG_KEYS,
    compute_static_worst_case,
    write_odesfr_exogenous_csvs,
    clone_config_for_mode,
    extract_contingency,
    compare_contingency,
    run_uced_odesfr_link,
    run_uced_sc,
    nc_status_path,
    write_nc_status,
    read_nc_status,
    _demand_csv_paths,
    _require_secure_nc_status,
    run_stage_nc,
    run_stage_dc,
    run_stage_sc,
)


_BASE_VALUES = {
    "Hs_val": [1.0, 2.0, 3.0, 0.5],
    "Hv_val": [0.0, 1.0, 0.0, 2.0],
    "FCR_val": [10.0, 20.0, 5.0, 15.0],
    "FFR_val": [1.0, 2.0, 3.0, 4.0],
    "aFRR_val": [50.0, 60.0, 40.0, 55.0],
    "VIR_val": [0.0, 0.0, 5.0, 0.0],
    "mFRR_val": [100.0, 90.0, 110.0, 95.0],
}


def _make_summary_reserves(index=None):
    idx = index if index is not None else pd.date_range("2030-01-01", periods=4, freq="h", tz="UTC")
    n = len(idx)
    data = {col: [vals[i % len(vals)] for i in range(n)] for col, vals in _BASE_VALUES.items()}
    return pd.DataFrame(data, index=idx)


def _make_base_config(tmp_path):
    return {
        "SimulationDirectory": str(tmp_path / "base"),
        "FCRDemand": "", "FFRDemand": "", "aFRRUDemand": "", "aFRRDDemand": "", "VIRDemand": "",
        "SystemInertiaDemand": "",
        "some_other_key": "unchanged",
    }


# ---------------------------------------------------------------------------
# compute_static_worst_case
# ---------------------------------------------------------------------------

def test_static_worst_case_is_constant_at_the_annual_max():
    sr = _make_summary_reserves()
    frame = compute_static_worst_case(sr)
    assert (frame.nunique() == 1).all()
    assert frame["Hs_val"].iloc[0] == 3.0
    assert frame["Hv_val"].iloc[0] == 2.0
    assert frame["FCR_val"].iloc[0] == 20.0
    assert frame["aFRR_val"].iloc[0] == 60.0
    assert frame["VIR_val"].iloc[0] == 5.0
    assert list(frame.columns) == list(sr.columns)
    assert len(frame) == len(sr)


def test_static_worst_case_at_least_as_large_as_dynamic_everywhere():
    sr = _make_summary_reserves()
    static = compute_static_worst_case(sr)
    for col in sr.columns:
        assert (static[col].fillna(0) >= sr[col].fillna(0)).all()


# ---------------------------------------------------------------------------
# write_odesfr_exogenous_csvs
# ---------------------------------------------------------------------------

def test_write_odesfr_exogenous_csvs_creates_expected_files_and_headers(tmp_path):
    sr = _make_summary_reserves()
    out_dir = tmp_path / "Simulations" / "BAU_2030_DC"

    paths = write_odesfr_exogenous_csvs(sr, str(out_dir), case_prefix="BAU_2030_DC")

    assert set(paths.keys()) == set(ALL_DEMAND_CONFIG_KEYS)
    assert paths["aFRRUDemand"] == paths["aFRRDDemand"]  # same symmetric file for both

    inertia_file = out_dir / "InertiaLimit" / "BAU_2030_DCH.csv"
    assert inertia_file.exists()
    assert inertia_file.read_text().splitlines()[0] == "TIMESTAMP,H_val"
    h_read = pd.read_csv(inertia_file)
    assert list(h_read["H_val"]) == [1.0, 2.0, 3.0, 0.5]  # Hs_val only -- see module docstring

    fcr_file = out_dir / "Reserves" / "BAU_2030_DCFCR.csv"
    assert fcr_file.exists()
    assert list(pd.read_csv(fcr_file)["FCR_val"]) == list(sr["FCR_val"])

    # mFRR is NOT written -- no config key consumes it (build.py overwrites
    # finalTS['mFRRU'] with aFRRU regardless, ISSUE-UCED-014)
    assert not (out_dir / "Reserves" / "BAU_2030_DCmFRR.csv").exists()


def test_write_odesfr_exogenous_csvs_H_ignores_Hv_val_regression_guard(tmp_path):
    """Hv_val must NOT leak into the H (SystemInertiaDemand) file: virtual
    inertia is procured as the VIRU reserve (batteries only), not through
    EQ_Inertia_Balance, which is gated to cu(au) (conventional units) --
    DECISION-ODESFR-012. A large Hv_val with a small/zero Hs_val must produce
    an H_val equal to Hs_val alone."""
    idx = pd.date_range("2030-01-01", periods=2, freq="h", tz="UTC")
    sr = pd.DataFrame(
        {
            "Hs_val": [1.0, 0.0],
            "Hv_val": [1000.0, 1000.0],  # deliberately huge, must be ignored
            "FCR_val": [0.0, 0.0], "FFR_val": [0.0, 0.0],
            "aFRR_val": [0.0, 0.0], "VIR_val": [0.0, 0.0], "mFRR_val": [0.0, 0.0],
        },
        index=idx,
    )
    paths = write_odesfr_exogenous_csvs(sr, str(tmp_path), case_prefix="BAU_2030_DC")
    h_read = pd.read_csv(paths["SystemInertiaDemand"])
    assert list(h_read["H_val"]) == [1.0, 0.0]


def test_write_odesfr_exogenous_csvs_accepts_static_worst_case_frame(tmp_path):
    """compute_static_worst_case's output must be a drop-in input to
    write_odesfr_exogenous_csvs, since run_uced_sc chains the two."""
    sr = _make_summary_reserves()
    static = compute_static_worst_case(sr)
    paths = write_odesfr_exogenous_csvs(static, str(tmp_path), case_prefix="BAU_2030_SC")
    fcr_read = pd.read_csv(paths["FCRDemand"])
    assert (fcr_read["FCR_val"] == 20.0).all()


# ---------------------------------------------------------------------------
# clone_config_for_mode
# ---------------------------------------------------------------------------

def test_clone_config_for_mode_nc_blanks_all_demand_keys(tmp_path):
    base = _make_base_config(tmp_path)
    sim_dir = tmp_path / "Simulations" / "BAU_2030_NC"
    config = clone_config_for_mode(base, "NC", None, sim_dir)

    assert config["SimulationDirectory"] == str(sim_dir)
    for key in ALL_DEMAND_CONFIG_KEYS:
        assert config[key] == ""
    assert config["some_other_key"] == "unchanged"
    assert base["SimulationDirectory"] != str(sim_dir)  # base_config untouched


def test_clone_config_for_mode_dc_uses_given_csv_paths(tmp_path):
    base = _make_base_config(tmp_path)
    sr = _make_summary_reserves()
    sim_dir = tmp_path / "Simulations" / "BAU_2030_DC"
    csv_paths = write_odesfr_exogenous_csvs(sr, str(sim_dir), "BAU_2030_DC")

    config = clone_config_for_mode(base, "DC", csv_paths, sim_dir)
    assert config["FCRDemand"] == csv_paths["FCRDemand"]
    assert config["SystemInertiaDemand"] == csv_paths["SystemInertiaDemand"]


def test_clone_config_for_mode_rejects_invalid_mode(tmp_path):
    base = _make_base_config(tmp_path)
    with pytest.raises(ValueError):
        clone_config_for_mode(base, "BOGUS", None, tmp_path)


def test_clone_config_for_mode_nc_rejects_csv_paths(tmp_path):
    base = _make_base_config(tmp_path)
    sr = _make_summary_reserves()
    csv_paths = write_odesfr_exogenous_csvs(sr, str(tmp_path), "x")
    with pytest.raises(ValueError):
        clone_config_for_mode(base, "NC", csv_paths, tmp_path)


def test_clone_config_for_mode_dc_requires_csv_paths(tmp_path):
    base = _make_base_config(tmp_path)
    with pytest.raises(ValueError):
        clone_config_for_mode(base, "DC", None, tmp_path)


# ---------------------------------------------------------------------------
# extract_contingency / compare_contingency
# ---------------------------------------------------------------------------

def test_extract_contingency_reads_output_contingency():
    series = pd.Series([100.0, 200.0])
    results = {"OutputContingency": series, "OutputContingencyPerZone": "should not be read"}
    assert extract_contingency(results) is series


def test_compare_contingency_identical_series():
    idx = pd.date_range("2030-01-01", periods=3, freq="h")
    a = pd.Series([100.0, 200.0, 300.0], index=idx)
    result = compare_contingency(a, a.copy())
    assert result["n_hours"] == 3
    assert result["n_discrepant"] == 0
    assert result["max_abs_diff"] == 0.0
    assert result["discrepant_hours"] == []


def test_compare_contingency_detects_discrepancy_beyond_tolerance():
    idx = pd.date_range("2030-01-01", periods=3, freq="h")
    a = pd.Series([100.0, 200.0, 300.0], index=idx)
    b = pd.Series([100.5, 205.0, 300.0], index=idx)  # hour 0 within 1MW atol, hour 1 not
    result = compare_contingency(a, b, atol=1.0)
    assert result["n_hours"] == 3
    assert result["n_discrepant"] == 1
    assert result["discrepant_hours"] == [idx[1]]
    assert result["max_abs_diff"] == 5.0
    assert 0 < result["pct_discrepant"] < 100


# ---------------------------------------------------------------------------
# run_uced_odesfr_link (fake build_solve_fn / ODE-SFR, no GAMS)
# ---------------------------------------------------------------------------

def _fake_odesfr_factory(summary_reserves_by_call):
    """Returns a fake get_frequency_stability_reserves_fn that pops the next
    pre-built summary_reserves frame from summary_reserves_by_call each call."""
    calls = list(summary_reserves_by_call)

    def _fake(path, inputs, results, **kwargs):
        sr = calls.pop(0)
        return None, sr, None

    return _fake


def test_run_uced_odesfr_link_single_round_when_only_one_dc_dir_given(tmp_path):
    idx = pd.date_range("2030-01-01", periods=3, freq="h", tz="UTC")
    sr1 = _make_summary_reserves(idx)

    def _build_solve(config):
        # Contingency deliberately differs between NC and DC1 to exercise
        # a discrepant round 1 that nonetheless stops (only 1 dir given).
        if config["SimulationDirectory"].endswith("_NC"):
            contingency = pd.Series([500.0, 500.0, 500.0], index=idx)
        else:
            contingency = pd.Series([600.0, 600.0, 600.0], index=idx)
        return {"config": config}, {"OutputContingency": contingency}

    report = run_uced_odesfr_link(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve,
        get_frequency_stability_reserves_fn=_fake_odesfr_factory([sr1]),
        simulation_dir_nc=str(tmp_path / "Simulations" / "BAU_2030_NC"),
        simulation_dirs_dc=[str(tmp_path / "Simulations" / "BAU_2030_DC1")],
        base_case_name="BAU_2030",
        atol=1.0,
    )

    assert report["rounds_run"] == 1
    assert report["round_1_comparison"]["n_discrepant"] == 3  # 100 MW gap, all hours
    assert report["round_2_comparison"] is None
    assert report["summary_reserves_round2"] is None
    assert report["converged"] is False


def test_run_uced_odesfr_link_runs_round_2_when_discrepant_and_second_dir_given(tmp_path):
    idx = pd.date_range("2030-01-01", periods=3, freq="h", tz="UTC")
    sr1 = _make_summary_reserves(idx)
    sr2 = _make_summary_reserves(idx)

    contingencies = {"_NC": [500.0, 500.0, 500.0], "_DC1": [600.0, 600.0, 600.0],
                      "_DC2": [600.0, 600.0, 600.0]}  # DC1 == DC2 -> round 2 converges

    def _build_solve(config):
        for suffix, vals in contingencies.items():
            if config["SimulationDirectory"].endswith(suffix):
                return {"config": config}, {"OutputContingency": pd.Series(vals, index=idx)}
        raise AssertionError(f"unexpected SimulationDirectory {config['SimulationDirectory']}")

    report = run_uced_odesfr_link(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve,
        get_frequency_stability_reserves_fn=_fake_odesfr_factory([sr1, sr2]),
        simulation_dir_nc=str(tmp_path / "Simulations" / "BAU_2030_NC"),
        simulation_dirs_dc=[str(tmp_path / "Simulations" / "BAU_2030_DC1"),
                             str(tmp_path / "Simulations" / "BAU_2030_DC2")],
        base_case_name="BAU_2030",
        atol=1.0,
    )

    assert report["rounds_run"] == 2
    assert report["round_1_comparison"]["n_discrepant"] == 3
    assert report["round_2_comparison"]["n_discrepant"] == 0
    assert report["converged"] is True
    assert report["summary_reserves_round2"] is not None


def test_run_uced_odesfr_link_skips_round_2_when_round_1_converges(tmp_path):
    idx = pd.date_range("2030-01-01", periods=2, freq="h", tz="UTC")
    sr1 = _make_summary_reserves(idx)

    def _build_solve(config):
        contingency = pd.Series([500.0, 500.0], index=idx)  # identical everywhere
        return {"config": config}, {"OutputContingency": contingency}

    report = run_uced_odesfr_link(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve,
        get_frequency_stability_reserves_fn=_fake_odesfr_factory([sr1]),
        simulation_dir_nc=str(tmp_path / "Simulations" / "BAU_2030_NC"),
        simulation_dirs_dc=[str(tmp_path / "Simulations" / "BAU_2030_DC1"),
                             str(tmp_path / "Simulations" / "BAU_2030_DC2")],
        base_case_name="BAU_2030",
    )

    assert report["rounds_run"] == 1
    assert report["converged"] is True
    assert report["round_2_comparison"] is None


def test_run_uced_odesfr_link_requires_at_least_one_dc_dir(tmp_path):
    with pytest.raises(ValueError):
        run_uced_odesfr_link(
            base_config=_make_base_config(tmp_path),
            build_solve_fn=lambda c: (None, None),
            get_frequency_stability_reserves_fn=lambda *a, **k: (None, None, None),
            simulation_dir_nc=str(tmp_path / "NC"),
            simulation_dirs_dc=[],
            base_case_name="BAU_2030",
        )


# ---------------------------------------------------------------------------
# run_uced_sc
# ---------------------------------------------------------------------------

def test_run_uced_sc_builds_static_config_and_solves(tmp_path):
    idx = pd.date_range("2030-01-01", periods=3, freq="h", tz="UTC")
    sr = _make_summary_reserves(idx)
    seen_configs = []

    def _build_solve(config):
        seen_configs.append(config)
        return {"marker": "sc_inputs"}, {"marker": "sc_results"}

    sim_dir = tmp_path / "Simulations" / "BAU_2030_SC"
    result = run_uced_sc(
        base_config=_make_base_config(tmp_path),
        summary_reserves_for_worst_case=sr,
        build_solve_fn=_build_solve,
        simulation_dir_sc=str(sim_dir),
        base_case_name="BAU_2030",
    )

    assert result["inputs"] == {"marker": "sc_inputs"}
    assert result["results"] == {"marker": "sc_results"}
    assert result["config"]["SimulationDirectory"] == str(sim_dir)
    fcr_read = pd.read_csv(result["config"]["FCRDemand"])
    assert (fcr_read["FCR_val"] == sr["FCR_val"].max()).all()
    assert len(seen_configs) == 1


# ---------------------------------------------------------------------------
# NC status JSON (nc_status_path / write_nc_status / read_nc_status /
# _require_secure_nc_status) -- the stage 'nc' -> stage 'dc'/'sc' hand-off
# for the split-stage cluster design (run_all_simulations.py --stage).
# ---------------------------------------------------------------------------

def test_read_nc_status_returns_none_when_not_written(tmp_path):
    assert read_nc_status(tmp_path / "BAU_2030_NC") is None


def test_write_then_read_nc_status_roundtrip(tmp_path):
    sim_dir = tmp_path / "BAU_2030_NC"
    path = write_nc_status(sim_dir, secure=True, n_insecure=0, n_hours=8760, base_case_name="BAU_2030")
    assert path == nc_status_path(sim_dir)
    status = read_nc_status(sim_dir)
    assert status["secure"] is True
    assert status["n_insecure"] == 0
    assert status["n_hours"] == 8760
    assert status["base_case_name"] == "BAU_2030"
    assert "written_at" in status


def test_require_secure_nc_status_raises_when_missing(tmp_path):
    with pytest.raises(RuntimeError, match="run stage 'nc' first"):
        _require_secure_nc_status(tmp_path / "BAU_2030_NC")


def test_require_secure_nc_status_raises_when_insecure(tmp_path):
    sim_dir = tmp_path / "BAU_2030_NC"
    write_nc_status(sim_dir, secure=False, n_insecure=5, n_hours=8760, base_case_name="BAU_2030")
    with pytest.raises(RuntimeError, match="insecure"):
        _require_secure_nc_status(sim_dir)


def test_require_secure_nc_status_passes_when_secure(tmp_path):
    sim_dir = tmp_path / "BAU_2030_NC"
    write_nc_status(sim_dir, secure=True, n_insecure=0, n_hours=8760, base_case_name="BAU_2030")
    status = _require_secure_nc_status(sim_dir)  # must not raise
    assert status["secure"] is True


# ---------------------------------------------------------------------------
# _demand_csv_paths -- must stay in sync with write_odesfr_exogenous_csvs
# ---------------------------------------------------------------------------

def test_demand_csv_paths_matches_write_odesfr_exogenous_csvs(tmp_path):
    sr = _make_summary_reserves()
    out_dir = tmp_path / "Simulations" / "BAU_2030_DC"
    written = write_odesfr_exogenous_csvs(sr, str(out_dir), case_prefix="BAU_2030_DC1")
    rederived = _demand_csv_paths(str(out_dir), "BAU_2030_DC1")
    assert written == rederived


def test_demand_csv_paths_raises_if_missing_by_default(tmp_path):
    with pytest.raises(FileNotFoundError):
        _demand_csv_paths(str(tmp_path / "nothing_written_here"), "BAU_2030_DC1")


def test_demand_csv_paths_allows_missing_when_not_required(tmp_path):
    paths = _demand_csv_paths(str(tmp_path / "nothing_written_here"), "BAU_2030_DC1", require_exists=False)
    assert set(paths.keys()) == set(ALL_DEMAND_CONFIG_KEYS)


# ---------------------------------------------------------------------------
# run_stage_nc / run_stage_dc / run_stage_sc -- the split-stage building
# blocks (fake build_solve_fn / ODE-SFR, no GAMS), matching the Array 1
# (run_stage_nc) / Array 2 (run_stage_dc, run_stage_sc) cluster design.
# ---------------------------------------------------------------------------

def _secure_summary_reserves(idx):
    sr = _make_summary_reserves(idx)
    sr["status"] = True
    return sr


def _insecure_summary_reserves(idx):
    """One genuinely-infeasible-even-at-full-bounds row (status=False), with
    the extra fleet-level columns patch_degenerate_reserves's docstring says
    must be retained ('SystemInertia', 'VIR SystemInertia', 'FFR Gain',
    'FCR Gain', 'Contingency', 'Damping') -- all bounds 0 and a huge
    Contingency so the repair's own full-bound re-check still finds it
    infeasible (status stays False; this is a real adequacy gap, not a
    search artifact, which is exactly the case run_stage_nc must treat as
    insecure)."""
    sr = _make_summary_reserves(idx)
    sr["status"] = [True, False, True, True][: len(idx)]
    sr["SystemInertia"] = 0.0
    sr["VIR SystemInertia"] = 0.0
    sr["FFR Gain"] = 0.0
    sr["FCR Gain"] = 0.0
    sr["Contingency"] = 1.0e6
    sr["Damping"] = 1.0
    return sr


def _build_solve_by_suffix(contingencies):
    """config -> (inputs, results) with OutputContingency picked by
    SimulationDirectory suffix, e.g. {'_NC': [...], '_DC1': [...]}."""
    def _build_solve(config):
        for suffix, vals in contingencies.items():
            if config["SimulationDirectory"].endswith(suffix):
                idx = pd.date_range("2030-01-01", periods=len(vals), freq="h", tz="UTC")
                return {"config": config}, {"OutputContingency": pd.Series(vals, index=idx)}
        raise AssertionError(f"unexpected SimulationDirectory {config['SimulationDirectory']}")
    return _build_solve


def test_run_stage_nc_secure_writes_dc_and_sc_csvs_and_status(tmp_path):
    idx = pd.date_range("2030-01-01", periods=4, freq="h", tz="UTC")
    sr = _secure_summary_reserves(idx)
    dir_nc = tmp_path / "BAU_2030_NC"
    dir_dc = tmp_path / "BAU_2030_DC"
    dir_sc = tmp_path / "BAU_2030_SC"

    result = run_stage_nc(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve_by_suffix({"_NC": [500.0] * 4}),
        get_frequency_stability_reserves_fn=_fake_odesfr_factory([sr]),
        simulation_dir_nc=str(dir_nc), simulation_dir_dc=str(dir_dc), simulation_dir_sc=str(dir_sc),
        base_case_name="BAU_2030",
    )

    assert result["secure"] is True
    assert result["n_insecure"] == 0
    assert result["csv_paths_dc"] is not None
    assert result["csv_paths_sc"] is not None
    assert (dir_dc / "Reserves" / "BAU_2030_DC1FCR.csv").exists()
    assert (dir_sc / "Reserves" / "BAU_2030_SCFCR.csv").exists()

    status = read_nc_status(dir_nc)
    assert status["secure"] is True
    assert status["n_insecure"] == 0


def test_run_stage_nc_insecure_skips_dc_sc_csvs(tmp_path):
    idx = pd.date_range("2030-01-01", periods=4, freq="h", tz="UTC")
    sr = _insecure_summary_reserves(idx)
    dir_nc = tmp_path / "BAU_2030_NC"
    dir_dc = tmp_path / "BAU_2030_DC"
    dir_sc = tmp_path / "BAU_2030_SC"

    result = run_stage_nc(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve_by_suffix({"_NC": [500.0] * 4}),
        get_frequency_stability_reserves_fn=_fake_odesfr_factory([sr]),
        simulation_dir_nc=str(dir_nc), simulation_dir_dc=str(dir_dc), simulation_dir_sc=str(dir_sc),
        base_case_name="BAU_2030",
    )

    assert result["secure"] is False
    assert result["n_insecure"] == 1
    assert result["csv_paths_dc"] is None
    assert result["csv_paths_sc"] is None
    assert not dir_dc.exists()
    assert not dir_sc.exists()

    status = read_nc_status(dir_nc)
    assert status["secure"] is False
    assert status["n_insecure"] == 1


def test_run_stage_nc_without_sc_dir_skips_sc_csvs_only(tmp_path):
    """simulation_dir_sc=None (the run_uced_odesfr_link/--stage all path): DC1's
    CSVs are still written, SC's are left for a later, separate call."""
    idx = pd.date_range("2030-01-01", periods=3, freq="h", tz="UTC")
    sr = _secure_summary_reserves(idx)
    dir_nc, dir_dc = tmp_path / "BAU_2030_NC", tmp_path / "BAU_2030_DC"

    result = run_stage_nc(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve_by_suffix({"_NC": [500.0] * 3}),
        get_frequency_stability_reserves_fn=_fake_odesfr_factory([sr]),
        simulation_dir_nc=str(dir_nc), simulation_dir_dc=str(dir_dc), simulation_dir_sc=None,
        base_case_name="BAU_2030",
    )
    assert result["csv_paths_dc"] is not None
    assert result["csv_paths_sc"] is None


def test_run_stage_dc_raises_when_nc_stage_not_run(tmp_path):
    with pytest.raises(RuntimeError, match="run stage 'nc' first"):
        run_stage_dc(
            base_config=_make_base_config(tmp_path),
            build_solve_fn=lambda c: (None, None),
            simulation_dir_nc=str(tmp_path / "BAU_2030_NC"),
            simulation_dir_dc=str(tmp_path / "BAU_2030_DC"),
            base_case_name="BAU_2030",
        )


def test_run_stage_dc_raises_when_nc_insecure(tmp_path):
    dir_nc = tmp_path / "BAU_2030_NC"
    write_nc_status(dir_nc, secure=False, n_insecure=2, n_hours=10, base_case_name="BAU_2030")
    with pytest.raises(RuntimeError, match="insecure"):
        run_stage_dc(
            base_config=_make_base_config(tmp_path),
            build_solve_fn=lambda c: (None, None),
            simulation_dir_nc=str(dir_nc),
            simulation_dir_dc=str(tmp_path / "BAU_2030_DC"),
            base_case_name="BAU_2030",
        )


def test_run_stage_dc_solves_and_compares_using_injected_nc_results(tmp_path):
    idx = pd.date_range("2030-01-01", periods=3, freq="h", tz="UTC")
    sr = _secure_summary_reserves(idx)
    dir_nc, dir_dc = tmp_path / "BAU_2030_NC", tmp_path / "BAU_2030_DC"

    nc_stage = run_stage_nc(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve_by_suffix({"_NC": [500.0, 500.0, 500.0]}),
        get_frequency_stability_reserves_fn=_fake_odesfr_factory([sr]),
        simulation_dir_nc=str(dir_nc), simulation_dir_dc=str(dir_dc), simulation_dir_sc=None,
        base_case_name="BAU_2030",
    )

    dc_stage = run_stage_dc(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve_by_suffix({"_DC": [500.0, 500.0, 600.0]}),
        simulation_dir_nc=str(dir_nc), simulation_dir_dc=str(dir_dc), base_case_name="BAU_2030",
        nc_results=nc_stage["nc_results"],
    )
    assert dc_stage["comparison"]["n_discrepant"] == 1
    assert dc_stage["converged"] is False


def test_run_stage_dc_loads_nc_results_from_disk_when_not_given(tmp_path):
    """Standalone/cluster --stage dc path: nc_results=None, loaded via the
    injected load_nc_results_fn instead of a real GDX read."""
    idx = pd.date_range("2030-01-01", periods=2, freq="h", tz="UTC")
    sr = _secure_summary_reserves(idx)
    dir_nc, dir_dc = tmp_path / "BAU_2030_NC", tmp_path / "BAU_2030_DC"

    run_stage_nc(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve_by_suffix({"_NC": [500.0, 500.0]}),
        get_frequency_stability_reserves_fn=_fake_odesfr_factory([sr]),
        simulation_dir_nc=str(dir_nc), simulation_dir_dc=str(dir_dc), simulation_dir_sc=None,
        base_case_name="BAU_2030",
    )

    loaded = {"path": None}

    def _fake_loader(path):
        loaded["path"] = path
        return {"marker": "nc_inputs"}, {"OutputContingency": pd.Series([500.0, 500.0], index=idx)}

    dc_stage = run_stage_dc(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve_by_suffix({"_DC": [500.0, 500.0]}),
        simulation_dir_nc=str(dir_nc), simulation_dir_dc=str(dir_dc), base_case_name="BAU_2030",
        nc_results=None, load_nc_results_fn=_fake_loader,
    )
    assert loaded["path"] == str(dir_nc)
    assert dc_stage["converged"] is True


def test_run_stage_sc_raises_when_nc_stage_not_run(tmp_path):
    with pytest.raises(RuntimeError, match="run stage 'nc' first"):
        run_stage_sc(
            base_config=_make_base_config(tmp_path),
            build_solve_fn=lambda c: (None, None),
            simulation_dir_nc=str(tmp_path / "BAU_2030_NC"),
            simulation_dir_sc=str(tmp_path / "BAU_2030_SC"),
            base_case_name="BAU_2030",
        )


def test_run_stage_sc_independent_of_dc_solves_from_nc_written_csvs(tmp_path):
    idx = pd.date_range("2030-01-01", periods=3, freq="h", tz="UTC")
    sr = _secure_summary_reserves(idx)
    dir_nc, dir_dc, dir_sc = tmp_path / "BAU_2030_NC", tmp_path / "BAU_2030_DC", tmp_path / "BAU_2030_SC"

    run_stage_nc(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve_by_suffix({"_NC": [500.0] * 3}),
        get_frequency_stability_reserves_fn=_fake_odesfr_factory([sr]),
        simulation_dir_nc=str(dir_nc), simulation_dir_dc=str(dir_dc), simulation_dir_sc=str(dir_sc),
        base_case_name="BAU_2030",
    )

    seen = []

    def _build_solve(config):
        seen.append(config["SimulationDirectory"])
        return {"marker": "sc_inputs"}, {"marker": "sc_results"}

    sc_stage = run_stage_sc(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve,
        simulation_dir_nc=str(dir_nc), simulation_dir_sc=str(dir_sc), base_case_name="BAU_2030",
    )
    assert sc_stage["results"] == {"marker": "sc_results"}
    assert seen == [str(dir_sc)]  # never touches dir_dc
    fcr_read = pd.read_csv(sc_stage["config"]["FCRDemand"])
    assert (fcr_read["FCR_val"] == sr["FCR_val"].max()).all()


# ---------------------------------------------------------------------------
# force_write_csvs / require_secure -- TESTING-ONLY escape hatches for
# exercising the DC/SC mechanics on a short, possibly-insecure test horizon
# (run_all_simulations.py --force-continue; never set by the cluster job files)
# ---------------------------------------------------------------------------

def test_run_stage_nc_force_write_csvs_writes_despite_insecure(tmp_path):
    idx = pd.date_range("2030-01-01", periods=4, freq="h", tz="UTC")
    sr = _insecure_summary_reserves(idx)
    dir_nc, dir_dc, dir_sc = tmp_path / "BAU_2030_NC", tmp_path / "BAU_2030_DC", tmp_path / "BAU_2030_SC"

    result = run_stage_nc(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve_by_suffix({"_NC": [500.0] * 4}),
        get_frequency_stability_reserves_fn=_fake_odesfr_factory([sr]),
        simulation_dir_nc=str(dir_nc), simulation_dir_dc=str(dir_dc), simulation_dir_sc=str(dir_sc),
        base_case_name="BAU_2030", force_write_csvs=True,
    )
    assert result["secure"] is False  # the true screening result is still reported...
    assert result["csv_paths_dc"] is not None  # ...but the CSVs are written anyway
    assert result["csv_paths_sc"] is not None
    status = read_nc_status(dir_nc)
    assert status["secure"] is False  # status JSON is never lied to


def test_run_stage_dc_require_secure_false_proceeds_when_insecure(tmp_path):
    idx = pd.date_range("2030-01-01", periods=3, freq="h", tz="UTC")
    sr = _insecure_summary_reserves(idx)
    sr["status"] = False
    dir_nc, dir_dc = tmp_path / "BAU_2030_NC", tmp_path / "BAU_2030_DC"

    run_stage_nc(
        base_config=_make_base_config(tmp_path),
        build_solve_fn=_build_solve_by_suffix({"_NC": [500.0] * 3}),
        get_frequency_stability_reserves_fn=_fake_odesfr_factory([sr]),
        simulation_dir_nc=str(dir_nc), simulation_dir_dc=str(dir_dc), simulation_dir_sc=None,
        base_case_name="BAU_2030", force_write_csvs=True,
    )

    # Without require_secure=False this must raise (the normal, safe default):
    with pytest.raises(RuntimeError, match="insecure"):
        run_stage_dc(
            base_config=_make_base_config(tmp_path), build_solve_fn=_build_solve_by_suffix({"_DC": [500.0] * 3}),
            simulation_dir_nc=str(dir_nc), simulation_dir_dc=str(dir_dc), base_case_name="BAU_2030",
        )

    # With require_secure=False (TESTING ONLY) it proceeds:
    dc_stage = run_stage_dc(
        base_config=_make_base_config(tmp_path), build_solve_fn=_build_solve_by_suffix({"_DC": [500.0] * 3}),
        simulation_dir_nc=str(dir_nc), simulation_dir_dc=str(dir_dc), base_case_name="BAU_2030",
        nc_results={"OutputContingency": pd.Series([500.0] * 3, index=idx)}, require_secure=False,
    )
    assert dc_stage["converged"] is True


def test_run_stage_dc_require_secure_false_still_requires_status_file(tmp_path):
    with pytest.raises(RuntimeError, match="run stage 'nc' first"):
        run_stage_dc(
            base_config=_make_base_config(tmp_path),
            build_solve_fn=lambda c: (None, None),
            simulation_dir_nc=str(tmp_path / "BAU_2030_NC"),
            simulation_dir_dc=str(tmp_path / "BAU_2030_DC"),
            base_case_name="BAU_2030", require_secure=False,
        )
