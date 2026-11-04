"""Mapping of the mid-term scheduling (MTS) storage levels to the simulation time steps.

The MTS level of step d is the level at the END of step d (StorageLevel in UCM.gms) but is labelled with
the start of the step. The hourly target at the last hour of each MTS step must equal that level, for any
horizon / look-ahead, and the first row must be the level before the first simulated hour."""
import numpy as np
import pandas as pd
import pytest

from dispaset.preprocessing.preprocessing import mts_levels_to_simulation_steps


def _mts(mts_hours, n, seed=0):
    idx = pd.date_range("2030-01-01", periods=n, freq=pd.Timedelta(hours=mts_hours))
    rng = np.random.default_rng(seed)
    return pd.DataFrame({"A": np.clip(0.5 + np.cumsum(rng.normal(0, 0.02, n)), 0, 1),
                         "B": np.linspace(0.2, 0.9, n)}, index=idx)


@pytest.mark.parametrize("mts_hours,sim_hours", [(24, 1), (24, 0.25), (6, 1)])
@pytest.mark.parametrize("start", ["2030-01-01", "2030-09-15"])
def test_end_of_step_levels_match_mts(mts_hours, sim_hours, start):
    mts = _mts(mts_hours, int(365 * 24 / mts_hours))
    idx_long = pd.date_range(start, periods=int(30 * 24 / sim_hours), freq=pd.Timedelta(hours=sim_hours))
    out = mts_levels_to_simulation_steps(mts, idx_long, mts_hours, sim_hours, 0.5)
    step, dt = pd.Timedelta(hours=mts_hours), pd.Timedelta(hours=sim_hours)
    for t0, row in mts.iterrows():
        t_last = t0 + step - dt          # last simulation step of the MTS step
        if t_last in out.index:
            assert np.allclose(out.loc[t_last].values, row.values)


def test_initial_level_and_first_row():
    mts = _mts(24, 365)
    idx_long = pd.date_range("2030-01-01", periods=72, freq="h")
    out = mts_levels_to_simulation_steps(mts, idx_long, 24, 1, 0.5)
    assert out.index[0] == pd.Timestamp("2029-12-31 23:00")
    assert np.allclose(out.iloc[0].values, 0.5)                          # level before the first hour
    assert np.allclose(out.loc["2030-01-01 23:00"].values, mts.iloc[0].values)
    # mid-year start: level before the first hour = MTS level at the end of the previous day
    idx_long = pd.date_range("2030-09-15", periods=72, freq="h")
    out = mts_levels_to_simulation_steps(mts, idx_long, 24, 1, 0.5)
    assert np.allclose(out.iloc[0].values, mts.loc["2030-09-14"].values)


def test_cyclic_initial_level_and_lookahead_beyond_mts():
    mts = _mts(24, 365)
    idx_long = pd.date_range("2030-12-29", periods=24 * 5, freq="h")   # look-ahead past the MTS horizon
    out = mts_levels_to_simulation_steps(mts, idx_long, 24, 1, mts.iloc[-1])
    assert np.allclose(out.loc["2031-01-02 23:00"].values, mts.iloc[-1].values)
    assert not out.isna().any().any()
    out = mts_levels_to_simulation_steps(mts, pd.date_range("2030-01-01", periods=24, freq="h"), 24, 1,
                                         mts.iloc[-1])
    assert np.allclose(out.iloc[0].values, mts.iloc[-1].values)
