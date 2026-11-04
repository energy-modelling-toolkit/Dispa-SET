# -*- coding: utf-8 -*-
"""ISSUE-ODESFR-fallback: a stage (Hv, FFR, FCR) with no feasible value on its own window falls back
to the candidates [0, bound] (previously [0]), so that services that only secure the hour jointly
(e.g. FFR + FCR) are still tested together at their bounds by Stage 5, while 0 is still tried
first. Synthetic fleets only."""
import numpy as np
import pandas as pd

import dispaset.postprocessing.frequency_stability as pp


def _fleet():
    idx = pd.date_range("2050-03-07 02:00", periods=1, freq="h")
    units = pd.DataFrame({"Technology": ["GTUR", "BATS"], "Fuel": ["GAS", "OTH"],
                          "InertiaConstant": [4.0, 0.0], "PowerCapacity": [1000.0, 5.0],
                          "Nunits": [1, 13], "Droop": [0.05, 0.02]}, index=["SYNC", "BESS"])
    inputs = {"units": units, "param_df": {"Demand": pd.DataFrame({"DA_Z1": [6000.0]}, index=idx)}}
    return inputs, {"OutputContingency": pd.Series([100.0], index=idx)}


def _fake_joint(seen):
    """Stage windows: Hs alone passes, any FFR/FCR test on its own window fails; Stage 5 passes only
    with FFR and FCR both at their bounds (65 MW/Hz and 1000/2.5 + 65 MW/Hz)."""
    ffr_max, fcr_max = 65.0, 1000.0 / 2.5 + 65.0

    def fake(sim_time, act, Hs, Hv, FFR, FCR, aFRR, mFRR, ctg, damping, verbose=False):
        if aFRR > 0 or mFRR > 0:
            seen.append((FFR, FCR))
            ok = FFR >= ffr_max - 1e-9 and FCR >= fcr_max - 1e-9
        else:
            ok = FFR == 0 and FCR == 0
        fd = 0.1 if ok else 10.0
        df = pd.DataFrame({'Time [s]': [0.0, 400.0, 800.0], 'Frequency Deviation [Hz]': [0.1, fd, fd]})
        return (fd, fd, 0.0, 0.0, FFR, FCR, aFRR, mFRR, df)
    return fake


def test_ffr_and_fcr_are_tested_jointly_at_their_bounds(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(pp, "frequency_response", _fake_joint(seen))
    inputs, results = _fleet()
    _, s, _ = pp.get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    row = s.iloc[0]
    assert bool(row["status"]) is True
    assert abs(row["FFR_val"] - 65.0) < 1e-9 and abs(row["FCR_val"] - 465.0) < 1e-9
    assert (65.0, 465.0) in seen
    assert seen[0] == (0, 0)           # 0 still tried first (the only test of the previous code)


def test_patch_degenerate_reserves_falls_back_to_bounds(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(pp, "frequency_response", _fake_joint(seen))
    inputs, results = _fleet()
    _, s, _ = pp.get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    s2 = s.copy(); s2.loc[:, "status"] = False
    patched, log = pp.patch_degenerate_reserves(s2)
    assert bool(patched.iloc[0]["status"]) is True
    assert abs(patched.iloc[0]["FFR_val"] - 65.0) < 1e-9


def test_unused_service_still_zero(tmp_path, monkeypatch):
    """Without a battery fleet (FFR/VIR bounds 0) the FFR and Hv stages stay at 0."""
    seen = []

    def fake(sim_time, act, Hs, Hv, FFR, FCR, aFRR, mFRR, ctg, damping, verbose=False):
        seen.append((Hv, FFR))
        df = pd.DataFrame({'Time [s]': [0.0, 400.0, 800.0], 'Frequency Deviation [Hz]': [0.1, 0.1, 0.1]})
        return (0.1, 0.1, 0.0, Hv, FFR, FCR, aFRR, mFRR, df)
    monkeypatch.setattr(pp, "frequency_response", fake)
    inputs, results = _fleet()
    inputs["units"] = inputs["units"].loc[["SYNC"]]
    _, s, _ = pp.get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    assert all(hv == 0 and ffr == 0 for hv, ffr in seen)
    assert s.iloc[0]["Hv_val"] == 0


def test_zero_still_preferred_when_it_suffices(tmp_path, monkeypatch):
    """Lexicographic order kept: if FFR = 0 with FCR at its bound passes Stage 5, FFR stays 0."""
    seen = []

    def fake(sim_time, act, Hs, Hv, FFR, FCR, aFRR, mFRR, ctg, damping, verbose=False):
        if aFRR > 0 or mFRR > 0:
            seen.append((FFR, FCR)); ok = FCR >= 465.0 - 1e-9
        else:
            ok = FFR == 0 and FCR == 0
        fd = 0.1 if ok else 10.0
        df = pd.DataFrame({'Time [s]': [0.0, 400.0, 800.0], 'Frequency Deviation [Hz]': [0.1, fd, fd]})
        return (fd, fd, 0.0, 0.0, FFR, FCR, aFRR, mFRR, df)
    monkeypatch.setattr(pp, "frequency_response", fake)
    inputs, results = _fleet()
    _, s, _ = pp.get_frequency_stability_reserves(str(tmp_path) + "/", inputs, results)
    assert bool(s.iloc[0]["status"]) is True
    assert s.iloc[0]["FFR_val"] == 0 and abs(s.iloc[0]["FCR_val"] - 465.0) < 1e-9
