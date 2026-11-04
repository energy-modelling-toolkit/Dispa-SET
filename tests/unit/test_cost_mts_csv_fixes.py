# -*- coding: utf-8 -*-
"""
Unit tests (no GAMS, seconds) for three small fixes:

* postprocessing.CostExPost: FixedCosts is the sum over ALL units of CostFixed x Committed
  (it used to keep only the last unit: ``= +`` instead of ``+=``).
* preprocessing.mid_term_scheduling: when no storage unit is left in the MTS problem (short-duration
  storage excluded) the MTS warns and continues (no GAMS solve, empty level targets) instead of sys.exit.
* reserve_link.write_odesfr_exogenous_csvs: the exogenous CSVs are written on exactly the DC config's
  index, and load_time_series reads them for a partial (3-day) horizon.
"""
from __future__ import annotations

import datetime as dt
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import _helpers  # noqa: E402,F401  (Agg backend)


# --------------------------------------------------------------------------- #
# CostExPost FixedCosts
# --------------------------------------------------------------------------- #
def test_costexpost_fixed_costs_sum_over_units():
    from dispaset.postprocessing.postprocessing import CostExPost
    t = pd.date_range("2030-01-01", periods=3, freq="h")
    units = ["A", "B"]

    def unit_param(name, values):
        return pd.DataFrame({name: values}, index=units)

    empty = pd.DataFrame(index=t)
    results = {
        "OutputPower": pd.DataFrame({"A": [10.0] * 3, "B": [5.0] * 3}, index=t),
        "OutputCommitted": pd.DataFrame({"A": [1, 1, 1], "B": [2, 2, 0]}, index=t),
        "OutputFlow": empty, "OutputShedLoad": empty, "OutputHeat": empty, "OutputSpillage": empty,
        "OutputSystemCost": pd.Series(0.0, index=t),
    }
    for k in ("aFRRD", "aFRRU", "mFRRU", "Inertia", "MaxPower", "MinPower", "RampDown", "RampUp"):
        results["LostLoad_" + k] = empty
    inputs = {"param_df": {
        "CostFixed": unit_param("CostFixed", [10.0, 20.0]),
        "PowerInitial": unit_param("PowerInitial", [10.0, 5.0]),
        "CostStartUp": unit_param("CostStartUp", [0.0, 0.0]),
        "CostShutDown": unit_param("CostShutDown", [0.0, 0.0]),
        "CostRampUp": unit_param("CostRampUp", [0.0, 0.0]),
        "CostRampDown": unit_param("CostRampDown", [0.0, 0.0]),
        "CostVariable": pd.DataFrame(0.0, index=t, columns=units),
        "PriceTransmission": pd.DataFrame(), "CostLoadShedding": pd.DataFrame(),
        "CHPPowerLossFactor": pd.DataFrame(columns=["CHPPowerLossFactor"]),
    }}
    costs, _ = CostExPost(inputs, results)
    # A: 10 x 1, B: 20 x committed (2, 2, 0)
    assert list(costs["FixedCosts"]) == [50.0, 50.0, 10.0]


# --------------------------------------------------------------------------- #
# MTS without storage: warn and continue
# --------------------------------------------------------------------------- #
def test_mts_without_storage_warns_and_continues(tmp_path, monkeypatch, caplog):
    from dispaset.preprocessing import preprocessing as pp

    def fake_build(config, MTS=0, **kw):
        return {"units": pd.DataFrame(), "sets": {"s": [], "nx": [], "au": []}, "parameters": {}}

    def no_solve(*a, **kw):
        raise AssertionError("GAMS must not be called when no storage is left in the MTS")

    monkeypatch.setattr(pp, "build_single_run", fake_build)
    monkeypatch.setattr(pp, "solve_GAMS", no_solve)
    config = {"StartDate": (2030, 1, 1, 0, 0, 0), "StopDate": (2030, 1, 3, 23, 0, 0), "LookAhead": 0,
              "SimulationTimeStep": 1, "zones": ["Z1"], "mts_zones": ["Z1"], "InitialFinalReservoirLevel": 1,
              "ReservoirLevels": "", "SimulationType": "LP", "HydroScheduling": 1, "GAMS_folder": "",
              "SimulationDirectory": str(tmp_path), "SectorXFlexibleDemand": "", "SectorXFlexibleSupply": "",
              "default": {"ReservoirLevelInitial": 0.5}}
    with caplog.at_level(logging.WARNING):
        profiles, *_ = pp.mid_term_scheduling(config, TimeStep=24)
    assert profiles.shape[1] == 0
    assert "no storage unit" in caplog.text


# --------------------------------------------------------------------------- #
# Exogenous CSVs on the DC config index, partial horizon
# --------------------------------------------------------------------------- #
def test_exogenous_csvs_on_dc_index_partial_horizon(tmp_path):
    from dispaset.postprocessing.reserve_link import write_odesfr_exogenous_csvs, config_time_index
    from dispaset.preprocessing.data_handler import load_time_series
    config = {"StartDate": (2030, 3, 1, 0, 0, 0), "StopDate": (2030, 3, 3, 23, 0, 0)}
    idx = config_time_index(config)
    assert len(idx) == 72
    # ODE-SFR summary on a tz-aware index with one missing hour (as read back from NC results)
    src = idx.drop(idx[10]).tz_localize("UTC")
    summary = pd.DataFrame({c: np.arange(len(src), dtype=float) for c in
                            ("Hs_val", "VIR_val", "FFR_val", "FCR_val", "aFRR_val")}, index=src)
    paths = write_odesfr_exogenous_csvs(summary, tmp_path, "T_DC", index=idx)
    written = pd.read_csv(paths["FCRDemand"], index_col=0, parse_dates=True)
    assert list(written.index) == list(idx)
    cfg = {"idx": idx, "idx_long": pd.date_range(idx[0], idx[-1] + dt.timedelta(days=1), freq="h")}
    for key in ("SystemInertiaDemand", "VIRDemand", "FFRDemand", "FCRDemand", "aFRRUDemand"):
        data = load_time_series(cfg, paths[key])
        assert len(data) == len(cfg["idx_long"]) and not data.isna().any().any()
    assert load_time_series(cfg, paths["FCRDemand"]).iloc[:72, 0].tolist() == written.iloc[:, 0].tolist()


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
