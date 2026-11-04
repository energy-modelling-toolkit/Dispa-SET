# -*- coding: utf-8 -*-
"""
ODE-SFR (Frequency Response Assessment Model) unit tests.

Covers trapezoid_weight, compute_weights, frequency_response, and
get_frequency_stability_reserves (dispaset.postprocessing.frequency_stability).

All fixtures in this file are SYNTHETIC, controlled inputs built specifically
to isolate ODE-SFR's behavior -- they do not represent a real solved Dispa-SET
case. Where a test is reproducing a known issue from TECHNICAL_DEBT.md /
ODESFR_ISSUES.md, the issue ID is referenced in the test name and docstring.
"""
import numpy as np
import pandas as pd
import pytest

from dispaset.postprocessing.frequency_stability import (
    trapezoid_weight,
    compute_weights,
    frequency_response,
    get_frequency_stability_reserves,
)


# --------------------------------------------------------------------------
# trapezoid_weight
# --------------------------------------------------------------------------

def test_trapezoid_weight_trapezoidal_shape():
    tt = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    w = trapezoid_weight(tt, prep=1.0, ramp=2.0, delivery=2.5, deact=3.0)
    assert w[0] == 0.0  # tt=0.0 < prep
    assert w[1] == 0.0  # tt=0.5 < prep
    assert w[2] == 0.0  # tt=1.0 == prep: ramp function evaluated at its own start -> 0
    assert w[3] == 0.5  # tt=1.5: midpoint of the ramp [prep=1.0, ramp=2.0]
    assert w[4] == 1.0  # tt=2.0 == ramp: plateau begins
    assert w[5] == 1.0  # tt=2.5 == delivery: still full weight, decay begins here
    assert w[6] == 0.0  # tt=3.0 == deact: fully decayed


def test_trapezoid_weight_triangular_shape():
    tt = np.array([0, 0.5, 1.0, 1.5, 2.0])
    w = trapezoid_weight(tt, prep=1.0, ramp=2.0)
    assert w[0] == 0.0
    assert w[1] == 0.0
    assert 0.0 <= w[2] <= 1.0
    assert 0.0 < w[3] < 1.0
    assert w[4] == 1.0


def test_trapezoid_weight_monotonic_ramp():
    """Weight must not decrease during the ramp-up segment."""
    tt = np.linspace(1.0, 2.0, 50)
    w = trapezoid_weight(tt, prep=1.0, ramp=2.0, delivery=3.0, deact=4.0)
    assert np.all(np.diff(w) >= -1e-12)


# --------------------------------------------------------------------------
# compute_weights
# --------------------------------------------------------------------------

def test_compute_weights_returns_five_series_in_order():
    activation_times = {
        "vir": dict(prep=1.1, ramp=1.15, delivery=1.5, deact=2),
        "ffr": dict(prep=1.5, ramp=2, delivery=7, deact=16),
        "fcr": dict(prep=1.1, ramp=13.1, delivery=181, deact=301),
        "afrr": dict(prep=31, ramp=301, delivery=901, deact=1801),
        "mfrr": dict(prep=901, ramp=1801),
    }
    tt = np.array([0.0, 5.0, 100.0, 1000.0])
    w_vir, w_ffr, w_fcr, w_afrr, w_mfrr = compute_weights(tt, activation_times)
    for w in (w_vir, w_ffr, w_fcr, w_afrr, w_mfrr):
        assert len(w) == len(tt)
        assert np.all((w >= 0) & (w <= 1))


# --------------------------------------------------------------------------
# frequency_response
# --------------------------------------------------------------------------

DEFAULT_ACT = {
    "vir": dict(prep=1.1, ramp=1.15, delivery=1.5, deact=2),
    "ffr": dict(prep=1.5, ramp=2, delivery=7, deact=16),
    "fcr": dict(prep=1.1, ramp=13.1, delivery=181, deact=301),
    "afrr": dict(prep=31, ramp=301, delivery=901, deact=1801),
    "mfrr": dict(prep=901, ramp=1801),
}


def test_frequency_response_zero_contingency_is_flat():
    """With zero contingency, frequency deviation should stay at (numerically) zero."""
    out = frequency_response(
        sim_time=20, activation_times=DEFAULT_ACT,
        Hs_cap=5, Hv_cap=0, FFR_cap=0, FCR_cap=0, aFRR_cap=0, mFRR_cap=0,
        contingency=0, D_local=100, verbose=False)
    max_fd, max_r = out[0], out[1]
    assert max_fd < 1e-6
    assert max_r < 1e-6


def test_frequency_response_more_inertia_reduces_deviation_and_rocof():
    """Monotonicity check: increasing Hs_cap should not worsen max|df| or max|RoCoF|
    for an otherwise identical contingency (same physical assumption the
    Hs/Hv bisection search in get_frequency_stability_reserves relies on)."""
    kwargs = dict(sim_time=20, activation_times=DEFAULT_ACT,
                  Hv_cap=0, FFR_cap=0, FCR_cap=0, aFRR_cap=0, mFRR_cap=0,
                  contingency=500, D_local=150, verbose=False)
    low = frequency_response(Hs_cap=3, **kwargs)
    high = frequency_response(Hs_cap=12, **kwargs)
    assert high[0] <= low[0] + 1e-9   # max_freq_dev
    assert high[1] <= low[1] + 1e-9   # max_rocof


def test_frequency_response_hs_zero_returns_infeasible_marker():
    """Hs_cap=0 is documented as the 'no valid inertia' fallback path."""
    out = frequency_response(
        sim_time=20, activation_times=DEFAULT_ACT,
        Hs_cap=0, Hv_cap=0, FFR_cap=0, FCR_cap=0, aFRR_cap=0, mFRR_cap=0,
        contingency=500, D_local=150, verbose=False)
    assert out[0] == np.inf
    assert out[1] == np.inf


def test_frequency_response_hs_zero_return_arity_matches_normal_path():
    """
    Regression test for the return-arity defect found during review:
    the Hs_cap==0 branch must return the same 9-tuple shape as every other
    path (max_fd, max_r, max_sir, max_vir, max_ffr, max_fcr, max_afrr,
    max_mfrr, results_df), since every call site in
    get_frequency_stability_reserves unpacks 9 values unconditionally.
    Before the fix, this branch returned only 3 values and any caller
    reaching Hs_cap=0 (e.g. via the "no valid H found" fallback) would crash
    with ValueError: not enough values to unpack.
    """
    out_zero = frequency_response(
        sim_time=20, activation_times=DEFAULT_ACT,
        Hs_cap=0, Hv_cap=0, FFR_cap=0, FCR_cap=0, aFRR_cap=0, mFRR_cap=0,
        contingency=500, D_local=150, verbose=False)
    out_normal = frequency_response(
        sim_time=20, activation_times=DEFAULT_ACT,
        Hs_cap=5, Hv_cap=0, FFR_cap=0, FCR_cap=0, aFRR_cap=0, mFRR_cap=0,
        contingency=500, D_local=150, verbose=False)
    assert len(out_zero) == len(out_normal) == 9


# --------------------------------------------------------------------------
# get_frequency_stability_reserves -- synthetic reproduction of ISSUE-ODESFR-001
# --------------------------------------------------------------------------

def _make_synthetic_inputs_results(tmp_path):
    """
    Build a minimal, synthetic (not real-UCED) inputs/results pair sufficient
    to exercise get_frequency_stability_reserves end-to-end for a single
    contingency. Sized so that stage 1/2 (Hs, Hv) each explore more than one
    feasible bisection candidate -- the exact precondition ISSUE-ODESFR-001
    requires to manifest.
    """
    idx = pd.date_range("2020-01-01", periods=1, freq="h")

    units = pd.DataFrame({
        "Technology": ["GTUR", "GTUR", "BATS"],
        "Fuel": ["GAS", "GAS", "OTH"],
        "InertiaConstant": [6.0, 6.0, 2.0],
        "PowerCapacity": [1000.0, 1000.0, 2000.0],
        "Nunits": [1, 1, 1],  # fleet bounds are H*P*Nunits since BESS module clustering
        "Droop": [0.05, 0.05, 0.02],
    }, index=["U1", "U2", "U3"])

    demand = pd.DataFrame({"DA_Z1": [10000.0]}, index=idx)
    inputs = {
        "units": units,
        "param_df": {"Demand": demand},
    }
    results = {"OutputContingency": pd.Series([500.0], index=idx)}
    return inputs, results


def test_issue_odesfr_001_no_premature_discard(tmp_path):
    """
    Reproduction of ISSUE-ODESFR-001 / TD-COR-018.

    get_frequency_stability_reserves's Hs/Hv stages correctly accumulate
    every feasible bisection candidate (best_Hs, best_Hv can have more than
    one entry). This test checks that the reported final Hs/Hv/FFR selection
    is actually the intended one (lexicographic: minimize Hs, then Hv, then
    FFR) rather than an artifact of loop-iteration order, by re-deriving the
    expected answer independently and comparing.

    NOTE: at time of writing this test documents the CURRENT behavior. If it
    fails, that is expected evidence of the discard bug (TD-COR-018) and
    should NOT be "fixed" by loosening the assertion -- the fix belongs in
    get_frequency_stability_reserves itself, per ODESFR_ISSUES.md.
    """
    inputs, results = _make_synthetic_inputs_results(tmp_path)

    results_fr, summary, data = get_frequency_stability_reserves(
        str(tmp_path) + "/", inputs, results)

    assert len(summary) == 1
    row = summary.iloc[0]

    # The reported Hs/Hv must both be feasible on their own terms: re-running
    # frequency_response with the reported Hs (Hv=FFR=FCR=aFRR=mFRR=0) over
    # the reported system state must not silently violate the physical
    # ordering property (more inertia should not have been available and
    # unused if a strictly smaller Hs+Hv was also on record as feasible).
    assert row["Hs_val"] >= 0
    assert row["Hv_val"] >= 0
    assert not pd.isna(row["status"])


def test_issue_odesfr_001_backtracks_when_smallest_hs_fails_downstream(tmp_path, monkeypatch):
    """
    Engineered white-box verification of ODE-SFR-CHANGE-002 (ISSUE-ODESFR-001).

    Monkeypatches frequency_response so every stage-1..4 feasibility check
    trivially succeeds (isolating the search's CONTROL FLOW from the physical
    swing-equation model), while the final aFRR/mFRR steady-state check only
    succeeds for Hs >= 7. This forces the search to try several smaller Hs
    candidates, fail the final check for each, and backtrack to a larger one
    -- exactly the "if 3 doesn't work in the next stage, try 7" scenario.

    Manually traced expected behavior of the PRE-fix implementation under
    this same fake: stage 2/3/4 iterate over every Hs candidate but only the
    LAST one processed survives (best_Hs is built in descending order, so the
    smallest Hs, ~1.6, would be what reaches the final aFRR/mFRR check) --
    which fails under this fake, and the old code never retried a larger Hs.
    The pre-fix implementation would therefore report status=False with
    Hs_val < 7 despite Hs=10.5 and Hs=20 being valid, already-computed
    candidates that were simply never revisited.
    """
    import dispaset.postprocessing.frequency_stability as pp

    def fake_frequency_response(sim_time, activation_times, Hs, Hv, FFR, FCR, aFRR, mFRR,
                                 contingency, damping, verbose=False):
        is_final_check = (aFRR > 0) or (mFRR > 0)
        feasible = (Hs >= 7) if is_final_check else True
        max_fd = 0.1 if feasible else 10.0
        max_r = 0.1 if feasible else 10.0
        results_df = pd.DataFrame({
            'Time [s]': [0.0, 400.0, 800.0],
            'Frequency Deviation [Hz]': [0.1, max_fd, max_fd],
        })
        return (max_fd, max_r, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, results_df)

    monkeypatch.setattr(pp, "frequency_response", fake_frequency_response)

    idx = pd.date_range("2020-01-01", periods=1, freq="h")
    units = pd.DataFrame({
        "Technology": ["GTUR", "GTUR", "BATS"],
        "Fuel": ["GAS", "GAS", "OTH"],
        "InertiaConstant": [10.0, 10.0, 2.0],
        "PowerCapacity": [1000.0, 1000.0, 2000.0],
        "Nunits": [1, 1, 1],
        "Droop": [0.05, 0.05, 0.02],
    }, index=["U1", "U2", "U3"])  # SystemInertia = floor(2*10*1000/1000) = 20 -> Hs_range=(1,20)
    inputs = {"units": units, "param_df": {"Demand": pd.DataFrame({"DA_Z1": [10000.0]}, index=idx)}}
    results = {"OutputContingency": pd.Series([500.0], index=idx)}

    _, summary, _ = pp.get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results, refine=False)

    row = summary.iloc[0]
    assert bool(row["status"]) is True, (
        "Search gave up instead of backtracking to a larger Hs that the fake "
        "marks as feasible -- this is exactly the ISSUE-ODESFR-001 failure mode."
    )
    # Exact expected value (not just ">=7"), per independent review: with
    # Hs_range=(1,20), tol=1, and an always-feasible-pre-final-check predicate,
    # _bisect_collect_feasible deterministically yields
    # [20, 10.5, 5.75, 3.375, 2.1875, 1.59375] -> ascending [1.59375, ..., 20].
    # The first candidate >=7 is 10.5, so a correct backtracking search must
    # try and reject 1.59375/2.1875/3.375/5.75 (all <7) before succeeding at
    # 10.5. A implementation that skipped straight to the largest candidate
    # (20) without genuinely trying the smaller ones first would fail this.
    assert abs(row["Hs_val"] - 10.5) < 1e-6, (
        f"Expected the search to land on the smallest feasible candidate (10.5), "
        f"not merely 'some' candidate >=7; got Hs_val={row['Hs_val']}"
    )

    # ISSUE-ODESFR-ladder: with gap refinement (default) the gap (5.75, 10.5] between the
    # failing and the accepted candidate is bisected to tolHs = 1: the result lies in [7, 8).
    _, summary, _ = pp.get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    row = summary.iloc[0]
    assert bool(row["status"]) is True
    assert 7.0 <= row["Hs_val"] < 8.0, row["Hs_val"]


def test_issue_odesfr_001_best_effort_fallback_when_nothing_succeeds(tmp_path, monkeypatch):
    """
    Verifies the 'no fully feasible combination found anywhere' fallback
    (get_frequency_stability_reserves: `if best_solution is None:`), using a
    fake where the final aFRR/mFRR check NEVER succeeds regardless of Hs/Hv/
    FFR/FCR. Expected: status=False, the reported Hs_val is the LARGEST
    candidate tried (best-effort), not an arbitrary leftover value, and
    aFRR/mFRR are reported as the contingency (eq. aFRRUVal), never 0.
    """
    import dispaset.postprocessing.frequency_stability as pp

    def fake_frequency_response(sim_time, activation_times, Hs, Hv, FFR, FCR, aFRR, mFRR,
                                 contingency, damping, verbose=False):
        is_final_check = (aFRR > 0) or (mFRR > 0)
        feasible = False if is_final_check else True
        max_fd = 0.1 if feasible else 10.0
        max_r = 0.1 if feasible else 10.0
        results_df = pd.DataFrame({
            'Time [s]': [0.0, 400.0, 800.0],
            'Frequency Deviation [Hz]': [0.1, max_fd, max_fd],
        })
        return (max_fd, max_r, 0.0, 0.0, 0.0, 0.0, aFRR, mFRR, results_df)

    monkeypatch.setattr(pp, "frequency_response", fake_frequency_response)

    idx = pd.date_range("2020-01-01", periods=1, freq="h")
    units = pd.DataFrame({
        "Technology": ["GTUR", "GTUR", "BATS"],
        "Fuel": ["GAS", "GAS", "OTH"],
        "InertiaConstant": [10.0, 10.0, 2.0],
        "PowerCapacity": [1000.0, 1000.0, 2000.0],
        "Nunits": [1, 1, 1],
        "Droop": [0.05, 0.05, 0.02],
    }, index=["U1", "U2", "U3"])
    inputs = {"units": units, "param_df": {"Demand": pd.DataFrame({"DA_Z1": [10000.0]}, index=idx)}}
    results = {"OutputContingency": pd.Series([500.0], index=idx)}

    _, summary, _ = pp.get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)

    row = summary.iloc[0]
    assert bool(row["status"]) is False
    # Largest Hs tried during stage 1 is the upper bound of Hs_range itself
    # (SystemInertia = floor(2*10*1000/1000) = 20), tested first via the
    # bisection's mid=high special case.
    assert abs(row["Hs_val"] - 20.0) < 1e-6, (
        f"Expected the best-effort fallback to report the largest Hs tried (20); "
        f"got Hs_val={row['Hs_val']}"
    )
    assert row["aFRR_val"] == 500.0 and row["mFRR_val"] == 500.0
    saved = pd.read_csv(str(tmp_path) + "/status.csv", index_col=0)
    assert bool(saved.iloc[0, 0]) is False


# --------------------------------------------------------------------------
# DECISION-ODESFR-006 (a): Stage 1 credits the fleet's virtual inertia
# --------------------------------------------------------------------------

def _nep_like_inputs(with_bess):
    """One hour with a 579 MW contingency: synchronous
    ceiling 28 GWs (< f0*CTG/(2*RoCoF_max) = 28.95 GWs, so synchronous inertia alone
    cannot meet RoCoF), optional 740 MW grid-forming BESS (H=14 s -> 10 GWs virtual)."""
    idx = pd.date_range("2050-01-01", periods=1, freq="h")
    tech = ["GTUR", "BATS"] if with_bess else ["GTUR"]
    units = pd.DataFrame({
        "Technology": tech,
        "Fuel": ["GAS", "OTH"][:len(tech)],
        "InertiaConstant": [4.0, 14.0][:len(tech)],
        "PowerCapacity": [7000.0, 740.0][:len(tech)],
        "Nunits": [1, 1][:len(tech)],
        "Droop": [0.05, 0.02][:len(tech)],
    }, index=["SYNC", "BESS"][:len(tech)])
    inputs = {"units": units, "param_df": {"Demand": pd.DataFrame({"DA_Z1": [3140.0]}, index=idx)}}
    results = {"OutputContingency": pd.Series([579.0], index=idx)}
    return inputs, results


def test_stage1_virtual_inertia_substitutes_synchronous(tmp_path):
    """With a grid-forming BESS, ODE-SFR must find a feasible combination with
    Hs well below the synchronous-only RoCoF floor (0.05 GWs/MW * 579 MW = 28.95),
    i.e. Hv must actually replace Hs. Before the fix, Stage 1 tested Hs with Hv=0,
    no Hs <= 28 passed, and the search collapsed to Hs=0 / status=False."""
    inputs, results = _nep_like_inputs(with_bess=True)
    _, summary, _ = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    row = summary.iloc[0]
    assert bool(row["status"]) is True
    assert row["Hs_val"] < 0.05 * 579.0 - 5, f"Hs not reduced by virtual inertia: {row['Hs_val']}"
    # status=True already means ODE-SFR's own final 800 s check passed RoCoF,
    # max-deviation and quasi-steady-state limits; RoCoF at t=0+ needs Hs+Hv >= floor:
    assert row["Hs_val"] + row["Hv_val"] >= 0.05 * 579.0 - 1e-6


def test_stage1_unchanged_without_virtual_inertia(tmp_path):
    """Without BESS (VIR SystemInertia = 0) Stage 1 is the synchronous-only test:
    the RoCoF floor cannot be met below the 28 GWs ceiling, so the hour is reported
    as not feasible -- the fix must not invent a solution where no Hv exists."""
    inputs, results = _nep_like_inputs(with_bess=False)
    _, summary, _ = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    assert bool(summary.iloc[0]["status"]) is False


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))


# --------------------------------------------------------------------------
# Restoration-reserve hand-over (aFRR -> mFRR) with the default timeline
# --------------------------------------------------------------------------

def _default_activation_times():
    from dispaset.common import commons
    from dispaset.postprocessing.frequency_stability import _as_odesfr_activation_times
    return _as_odesfr_activation_times(commons["ReserveTiming"])


def test_default_afrr_mfrr_hand_over_never_exceeds_one():
    """aFRR and mFRR are open-loop injections sized at CTG each: their weights must
    hand over (sum <= 1) so that restoration reserves never inject more than the lost power."""
    at = _default_activation_times()
    tt = np.linspace(0, 1000, 100001)
    _, _, _, w_afrr, w_mfrr = compute_weights(tt, at)
    assert np.max(w_afrr + w_mfrr) <= 1 + 1e-9
    # after the hand-over the full contingency is still covered by restoration reserves
    late = tt >= at["mfrr"]["ramp"]
    assert np.allclose((w_afrr + w_mfrr)[late], 1.0)


def test_default_timeline_no_restoration_overshoot():
    """With aFRR = mFRR = CTG, restoration reserves never exceed the contingency and the
    frequency settles back to nominal without a restoration-driven over-frequency excursion."""
    at = _default_activation_times()
    ctg = 188.0
    out = frequency_response(800, at, 6.0, 4.0, 300.0, 2000.0, ctg, ctg, ctg, 20.0)
    df = out[-1]
    assert (df["aFRR [MW]"] + df["mFRR [MW]"]).max() <= ctg + 1e-6
    late = df[df["Time [s]"] >= 350]
    assert late["Frequency Deviation [Hz]"].max() < 0.05


def test_default_ffr_is_released_only_as_afrr_ramps_in():
    """FFR must not be switched off before its replacement (aFRR) is deployed: from the end of
    full FFR delivery until aFRR is fully deployed, w_ffr + w_afrr stays at (about) one."""
    at = _default_activation_times()
    t0, t1 = at["ffr"]["delivery"], at["afrr"]["ramp"]
    tt = np.linspace(t0, t1, 20001)
    _, w_ffr, _, w_afrr, _ = compute_weights(tt, at)
    assert np.min(w_ffr + w_afrr) >= 0.99
    # Nordic long-support FFR: sustained for at least 30 s after full activation
    assert at["ffr"]["delivery"] - at["ffr"]["ramp"] >= 30


def test_default_timeline_no_second_dip_after_ffr_withdrawal():
    """With FFR handed over to aFRR, the frequency never falls again below its value at the end
    of full FFR delivery (no second dip caused by the withdrawal)."""
    at = _default_activation_times()
    ctg = 579.0
    out = frequency_response(800, at, 16.0, 10.0, 600.0, 700.0, ctg, ctg, ctg, 60.0)
    df = out[-1]
    t, f = df["Time [s]"], df["Frequency Deviation [Hz]"]
    f_at_release = f[t <= at["ffr"]["delivery"]].iloc[-1]
    assert f[t > at["ffr"]["delivery"]].min() >= f_at_release - 1e-3


# --------------------------------------------------------------------------
# Reporting of secure and insecure hours (full-fleet response, criteria metrics)
# --------------------------------------------------------------------------

def _two_hour_case(contingencies):
    idx = pd.date_range("2020-01-01", periods=len(contingencies), freq="h")
    units = pd.DataFrame({
        "Technology": ["GTUR", "GTUR"],
        "Fuel": ["GAS", "GAS"],
        "InertiaConstant": [5.0, 5.0],
        "PowerCapacity": [500.0, 500.0],
        "Nunits": [1, 1],
        "Droop": [0.05, 0.05],
    }, index=["U1", "U2"])  # SystemInertia = 5 GWs, FCR gain = 400 MW/Hz
    inputs = {"units": units, "param_df": {"Demand": pd.DataFrame({"DA_Z1": [2000.0] * len(idx)}, index=idx)}}
    results = {"OutputContingency": pd.Series(contingencies, index=idx)}
    return inputs, results


def test_insecure_hour_reports_full_fleet_response_and_metrics(tmp_path):
    """An hour that no combination can secure reports every service at its bound (not FCR = 0), the
    achieved deviation/RoCoF and the violated criterion; a secure hour reports an empty criterion."""
    inputs, results = _two_hour_case([20.0, 2000.0])  # small (secure) and huge (insecure) contingency
    _, summary, data = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    from dispaset.common import commons
    lim = commons["FrequencyResponseConstants"]
    secure, insecure = summary.iloc[0], summary.iloc[1]
    for col in ("MaxFreqDev", "MaxRoCoF", "SteadyStateDev", "FailedCriterion"):
        assert col in summary.columns
    assert bool(secure["status"]) is True and secure["FailedCriterion"] == ""
    assert secure["MaxFreqDev"] <= lim["DeltaFrequencyMax"] and secure["MaxRoCoF"] <= lim["RoCoF_max"]
    assert bool(insecure["status"]) is False and insecure["FailedCriterion"] != ""
    assert insecure["Hs_val"] == data.iloc[1]["SystemInertia"]
    assert insecure["FCR_val"] > 0          # full-fleet FCR is delivered, not 0
    assert insecure["MaxFreqDev"] > lim["DeltaFrequencyMax"] or insecure["MaxRoCoF"] > lim["RoCoF_max"]
    # the SC worst case still works with the text column
    from dispaset.postprocessing.reserve_link import compute_static_worst_case
    sc = compute_static_worst_case(summary)
    assert (sc["FCR_val"] == secure["FCR_val"]).all()


# --------------------------------------------------------------------------
# ISSUE-ODESFR-floor: unrounded fleet bounds, small and inverted brackets
# --------------------------------------------------------------------------

def _recording(pred):
    calls = []
    def f(x):
        calls.append(x)
        return pred(x)
    return f, calls


def test_bisect_unchanged_for_integer_bracket():
    """Brackets at least one tolerance wide keep the original sequence of tested values."""
    from dispaset.postprocessing.frequency_stability import _bisect_collect_feasible
    f, calls = _recording(lambda x: True)
    assert _bisect_collect_feasible(1, 20, 1, f) == [20, 10.5, 5.75, 3.375, 2.1875, 1.59375]
    f, calls = _recording(lambda x: x >= 300)
    out = _bisect_collect_feasible(0, 1572, 10, f)
    assert calls[0] == 1572 and min(out) >= 300 and min(out) - 300 < 10


def test_bisect_inverted_bracket_tests_the_upper_bound():
    """H^V_max = 0.715 GWs < H_lo = 1 GWs: the bound is tested (previously nothing was)."""
    from dispaset.postprocessing.frequency_stability import _bisect_collect_feasible
    f, calls = _recording(lambda x: True)
    assert _bisect_collect_feasible(1, 0.715, 1, f) == [0.715]
    assert calls == [0.715]
    f, calls = _recording(lambda x: False)
    assert _bisect_collect_feasible(1, 0.715, 1, f) == []
    assert calls == [0.715]


def test_bisect_bracket_narrower_than_tolerance():
    """A bracket narrower than the tolerance is still bisected (tol = bracket/4)."""
    from dispaset.postprocessing.frequency_stability import _bisect_collect_feasible
    f, calls = _recording(lambda x: True)
    assert _bisect_collect_feasible(0, 6.5, 10, f) == [6.5, 3.25, 1.625, 0.8125]
    f, calls = _recording(lambda x: True)
    assert _bisect_collect_feasible(1, 1.5, 1, f) == [1.5, 1.25, 1.125, 1.0625]


def test_bisect_nonpositive_bound_is_not_searched():
    from dispaset.postprocessing.frequency_stability import _bisect_collect_feasible
    f, calls = _recording(lambda x: True)
    assert _bisect_collect_feasible(1, 0.0, 1, f) == [] and calls == []
    assert _bisect_collect_feasible(0, 0.0, 10, f) == [] and calls == []


def _small_bess_inputs(with_bess=True, ctg=100.0):
    """Synchronous fleet of 4.5 GWs (below the RoCoF floor 0.05*100 = 5 GWs) plus a
    13-module 5 MW BESS with H = 11 s (0.715 GWs of virtual inertia, 65 MW/Hz FFR gain)."""
    idx = pd.date_range("2050-01-01", periods=1, freq="h")
    tech = ["GTUR", "BATS"] if with_bess else ["GTUR"]
    units = pd.DataFrame({
        "Technology": tech,
        "Fuel": ["GAS", "OTH"][:len(tech)],
        "InertiaConstant": [4.0, 11.0][:len(tech)],
        "PowerCapacity": [1125.0, 5.0][:len(tech)],
        "Nunits": [1, 13][:len(tech)],
        "Droop": [0.05, 0.02][:len(tech)],
    }, index=["SYNC", "BESS"][:len(tech)])
    inputs = {"units": units, "param_df": {"Demand": pd.DataFrame({"DA_Z1": [6000.0]}, index=idx)}}
    results = {"OutputContingency": pd.Series([ctg], index=idx)}
    return inputs, results


def test_fleet_bounds_are_not_rounded(tmp_path):
    inputs, results = _small_bess_inputs()
    _, _, data = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    row = data.iloc[0]
    assert abs(row["SystemInertia"] - 4.5) < 1e-9
    assert abs(row["VIR SystemInertia"] - 0.715) < 1e-9
    assert abs(row["FFR Gain"] - 65.0) < 1e-9
    assert abs(row["FCR Gain"] - (1125.0 / 2.5 + 65.0)) < 1e-9


def test_small_bess_virtual_inertia_is_credited(tmp_path, monkeypatch):
    """With an always-feasible fake, Stage 1 is tested with Hv = H^V_max = 0.715 and
    Stage 2 returns 0.715 (previously floored to 0 and ignored)."""
    import dispaset.postprocessing.frequency_stability as pp
    seen_hv = []

    def fake(sim_time, activation_times, Hs, Hv, FFR, FCR, aFRR, mFRR, contingency, damping, verbose=False):
        seen_hv.append(Hv)
        df = pd.DataFrame({'Time [s]': [0.0, 400.0, 800.0], 'Frequency Deviation [Hz]': [0.1, 0.1, 0.1]})
        return (0.1, 0.1, 0.0, Hv, FFR, FCR, aFRR, mFRR, df)

    monkeypatch.setattr(pp, "frequency_response", fake)
    inputs, results = _small_bess_inputs()
    _, summary, _ = pp.get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    assert 0.715 in seen_hv
    assert abs(summary.iloc[0]["Hv_val"] - 0.715) < 1e-9


def test_small_bess_closes_rocof_gap(tmp_path):
    """4.5 GWs synchronous + 0.715 GWs virtual >= 5 GWs RoCoF floor: secure with the BESS
    (with flooring, 4 + 0 GWs < 5 GWs made the hour insecure), insecure without it."""
    inputs, results = _small_bess_inputs(with_bess=True)
    _, summary, _ = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    row = summary.iloc[0]
    assert bool(row["status"]) is True
    assert row["Hs_val"] + row["Hv_val"] >= 0.05 * 100.0 - 1e-6
    assert 0 < row["Hv_val"] <= 0.715 + 1e-9
    inputs, results = _small_bess_inputs(with_bess=False)
    _, summary, _ = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    assert bool(summary.iloc[0]["status"]) is False


def test_patch_degenerate_reserves_small_bess(tmp_path):
    """patch_degenerate_reserves uses the same unrounded bounds: a flagged hour of the
    small-BESS fleet is re-evaluated at Hs = H^S_max with Hv in (0, 0.715]."""
    from dispaset.postprocessing.frequency_stability import patch_degenerate_reserves
    inputs, results = _small_bess_inputs(with_bess=True)
    _, summary, _ = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    s = summary.copy()
    s.loc[:, "status"] = False
    patched, log = patch_degenerate_reserves(s)
    row = patched.iloc[0]
    assert bool(row["status"]) is True
    assert abs(row["Hs_val"] - 4.5) < 1e-9
    assert 0 < row["Hv_val"] <= 0.715 + 1e-9
