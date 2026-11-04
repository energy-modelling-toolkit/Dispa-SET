# -*- coding: utf-8 -*-
"""
ISSUE-ODESFR-avail / ISSUE-ODESFR-ladder: ODE-SFR search bounds from the hourly AVAILABLE
fleet (nameplate x AF x (1 - OF), the UCED's LoadMaximum) and gap refinement of the
ordered-backtracking search. Synthetic fleets only.
"""
import numpy as np
import pandas as pd
import pytest

import dispaset.postprocessing.frequency_stability as fs
from dispaset.postprocessing.frequency_stability import (
    _ordered_backtracking,
    frequency_response,
    get_frequency_stability_reserves,
    patch_degenerate_reserves,
    resolve_activation_times,
    unit_hourly_availability,
)
from dispaset.common import commons


def _fleet(n_hours=3, af=None, of=None, ctg=300.0, demand=8000.0):
    """Two synchronous units (G1: 2 x 1000 MW, H = 5 s; G2: 1000 MW, H = 4 s -> 14 GWs) and a
    1000 MW PV plant. AF/OF are optional hourly DataFrames (hours x units)."""
    idx = pd.date_range("2050-01-01", periods=n_hours, freq="h")
    units = pd.DataFrame({
        "Technology": ["STUR", "GTUR", "PHOT"],
        "Fuel": ["GAS", "GAS", "SUN"],
        "InertiaConstant": [5.0, 4.0, np.nan],
        "PowerCapacity": [1000.0, 1000.0, 1000.0],
        "Nunits": [2, 1, 1],
        "Droop": [0.05, 0.05, 0.06],
    }, index=["G1", "G2", "PV"])
    pdf = {"Demand": pd.DataFrame({"DA_Z1": [demand] * n_hours}, index=idx)}
    if af is not None:
        pdf["AvailabilityFactor"] = af.set_axis(idx)
    if of is not None:
        pdf["OutageFactor"] = of.set_axis(idx)
    inputs = {"units": units, "param_df": pdf}
    results = {"OutputContingency": pd.Series([ctg] * n_hours, index=idx)}
    return inputs, results


def _af_of():
    af = pd.DataFrame({"G1": [1.0, 1.0, 1.0], "G2": [1.0, 1.0, 1.0], "PV": [0.0, 0.8, 0.0]})
    of = pd.DataFrame({"G1": [0.0, 0.5, 0.0], "G2": [0.0, 0.0, 1.0], "PV": [0.0, 0.0, 0.0]})
    return af, of


# --------------------------------------------------------------------------- availability

def test_unit_hourly_availability_is_af_times_one_minus_of():
    af, of = _af_of()
    inputs, _ = _fleet(af=af, of=of)
    idx = inputs["param_df"]["Demand"].index
    lm = unit_hourly_availability(inputs, ["G1", "G2", "PV"], idx)
    np.testing.assert_allclose(lm.values, (af * (1 - of)).values)


def test_unit_hourly_availability_defaults_to_installed_fleet():
    inputs, _ = _fleet()                      # no AF/OF in the inputs
    idx = inputs["param_df"]["Demand"].index
    lm = unit_hourly_availability(inputs, ["G1", "G2", "PV"], idx)
    assert (lm.values == 1.0).all()
    af, of = _af_of()                         # unit missing from the data -> available
    inputs, _ = _fleet(af=af[["G1", "PV"]], of=of[["G1", "PV"]])
    lm = unit_hourly_availability(inputs, ["G1", "G2", "PV"], idx)
    assert (lm["G2"] == 1.0).all() and lm["G1"].iloc[1] == 0.5


def test_unit_hourly_availability_from_gams_parameter_arrays():
    """Raw Inputs.p layout (parameters[name]['val'] with sets ['au', 'h']), aligned by position."""
    af, of = _af_of()
    idx = pd.date_range("2050-01-01", periods=3, freq="h")
    inputs = {"sets": {"au": ["G1", "G2", "PV"]},
              "parameters": {"AvailabilityFactor": {"sets": ["au", "h"], "val": af.values.T},
                             "OutageFactor": {"sets": ["au", "h"], "val": of.values.T}}}
    lm = unit_hourly_availability(inputs, ["G1", "G2", "PV"], idx)
    np.testing.assert_allclose(lm.values, (af * (1 - of)).values)


# --------------------------------------------------------------------------- hourly bounds

def test_hourly_bounds_follow_available_fleet(tmp_path):
    af, of = _af_of()
    inputs, results = _fleet(af=af, of=of)
    _, _, data = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    # H: G1 2*1000*5 = 10 GWs, G2 4 GWs; hour 2: G1 half out; hour 3: G2 out
    np.testing.assert_allclose(data["SystemInertia"].values, [14.0, 9.0, 10.0])
    kg = 1000.0 / (0.05 * 50)                 # 400 MW/Hz per 1000 MW at 5 % droop
    kpv = 1000.0 / (0.06 * 50)
    np.testing.assert_allclose(data["FCR Gain"].values,
                               [3 * kg, 2 * kg + 0.8 * kpv, 2 * kg])   # PV only in daylight
    # previous behaviour: constant nameplate bounds
    _, _, data_i = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results,
                                                    fleet_bounds="installed")
    np.testing.assert_allclose(data_i["SystemInertia"].values, [14.0] * 3)
    np.testing.assert_allclose(data_i["FCR Gain"].values, [3 * kg + kpv] * 3)


def test_invalid_fleet_bounds_option(tmp_path):
    inputs, results = _fleet()
    with pytest.raises(ValueError):
        get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results, fleet_bounds="foo")


def test_requirement_never_exceeds_available_bound(tmp_path):
    """Every reported requirement is within the hour's available bound; hours that the installed
    fleet secures only with units on outage (or PV at night) are flagged insecure."""
    af, of = _af_of()
    inputs, results = _fleet(af=af, of=of, ctg=100.0, demand=3000.0)
    _, s_av, _ = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    _, s_in, _ = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results,
                                                  fleet_bounds="installed")
    assert (s_av["Hs_val"].astype(float) <= s_av["SystemInertia"] + 1e-9).all()
    assert s_in["status"].astype(bool).all()
    # hours 2-3: G1 half out / G2 out -> not deliverable by the available fleet
    assert s_av["status"].astype(bool).tolist() == [True, False, False]
    # hour 1: PV unavailable at night -> less FCR gain, more synchronous inertia needed
    assert float(s_av["Hs_val"].iloc[0]) >= float(s_in["Hs_val"].iloc[0]) - 1e-9


def test_available_equals_installed_without_availability_data(tmp_path):
    inputs, results = _fleet(ctg=100.0, demand=3000.0)
    _, s_av, _ = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    _, s_in, _ = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results,
                                                  fleet_bounds="installed")
    np.testing.assert_allclose(s_av["Hs_val"].astype(float), s_in["Hs_val"].astype(float))
    np.testing.assert_allclose(s_av["FCR_val"].astype(float), s_in["FCR_val"].astype(float))


def test_patch_uses_hourly_bound(tmp_path):
    """patch_degenerate_reserves re-evaluates a flagged hour at that hour's available H^S_max."""
    af, of = _af_of()
    inputs, results = _fleet(af=af, of=of)
    _, s, _ = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    assert set(fs._FCR_SAT_BY_INDEX) == set(s.index)
    s2 = s.copy()
    s2.loc[:, "status"] = False
    patched, _ = patch_degenerate_reserves(s2)
    ok = patched["status"].astype(bool)
    np.testing.assert_allclose(patched.loc[ok, "Hs_val"].astype(float), s.loc[ok, "SystemInertia"])


# --------------------------------------------------------------------------- gap refinement

def test_ordered_backtracking_refines_the_gap():
    ladder = [1.59375, 2.1875, 3.375, 5.75, 10.5, 20.0]
    calls = []

    def cands(k, q):
        return ladder if k == 0 else [0]

    def accept(q):
        calls.append(q[0])
        return "ok" if q[0] >= 7.0 else None

    v, p = _ordered_backtracking(2, cands, accept, [None, None])
    assert v[0] == 10.5 and p == "ok"                       # previous behaviour
    v, p = _ordered_backtracking(2, cands, accept, [1.0, None])
    assert 7.0 <= v[0] < 8.0 and p == "ok"
    assert all(c >= 7.0 for c in [v[0]])                     # returned value passed accept
    assert _ordered_backtracking(2, cands, lambda q: None, [1.0, None]) is None


def test_ordered_backtracking_first_candidate_is_not_refined():
    """When the smallest candidate passes, there is no gap to refine (no extra evaluations)."""
    n = []

    def accept(q):
        n.append(1)
        return "ok"

    v, _ = _ordered_backtracking(1, lambda k, q: [3.0, 6.0], accept, [0.5])
    assert v == [3.0] and len(n) == 1


def test_refined_requirement_is_secure_and_not_larger(tmp_path):
    """Real swing equation: the refined requirement passes every criterion when re-simulated and is
    never larger than the unrefined one."""
    inputs, results = _fleet(n_hours=1, ctg=1000.0, demand=3000.0)
    _, s0, _ = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results, refine=False)
    _, s1, _ = get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results, refine=True)
    r0, r1 = s0.iloc[0], s1.iloc[0]
    assert bool(r1["status"]) is True or bool(r0["status"]) is False
    assert float(r1["Hs_val"]) <= float(r0["Hs_val"]) + 1e-9
    if bool(r1["status"]):
        c = commons["FrequencyResponseConstants"]
        act, _ = resolve_activation_times(None)
        assert float(r1["MaxFreqDev"]) <= c["DeltaFrequencyMax"] + 1e-9
        assert float(r1["MaxRoCoF"]) <= c["RoCoF_max"] + 1e-9
        assert r1["FailedCriterion"] == ""
