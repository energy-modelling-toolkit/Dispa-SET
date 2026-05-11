# -*- coding: utf-8 -*-
"""
Integration test: recursive reserve-demand activation on tiny system
====================================================================

What this test does
-------------------
Runs a sequence of very small exogenous-reserve scenarios on ``tiny.yml``.
Each step cumulatively adds one reserve-demand input file:

1. base (no exogenous reserve demand),
2. + aFRRU,
3. + aFRRD,
4. + FFR,
5. + FCR.

The test verifies that reserve activation has an observable impact on
dispatch and that battery units participate in reserve provision.

How to run
----------

1. ``pytest tests/integration/test_reserves_generic.py``
2. ``python tests/integration/test_reserves_generic.py``
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import build_solve, load_test_config, skip_if_no_gams  # noqa: E402


def _write_reserve_profile(path: Path, index: pd.DatetimeIndex, value_mw: float) -> str:
    df = pd.DataFrame({"Z1": value_mw}, index=index)
    df.to_csv(path)
    return str(path)


def _battery_power_series(results: dict) -> pd.Series:
    cols = [c for c in results["OutputPower"].columns if "BATS" in str(c)]
    if not cols:
        return pd.Series(dtype=float)
    return results["OutputPower"][cols].sum(axis=1)


def _battery_reserve_total(results: dict) -> float:
    if "OutputReserveProvision" not in results:
        return 0.0
    reserve_df = results["OutputReserveProvision"]
    if reserve_df.empty:
        return 0.0
    cols = [c for c in reserve_df.columns if "BATS" in str(c)]
    if not cols:
        return 0.0
    return float(reserve_df[cols].sum().sum())


@pytest.mark.timeout(240)
def test_recursive_reserve_demands_affect_dispatch_with_battery():
    skip_if_no_gams()

    cfg0 = load_test_config("tiny.yml", "reserves_recursive_base")
    base_index = pd.date_range(
        start=pd.Timestamp(*cfg0["StartDate"]),
        end=pd.Timestamp(*cfg0["StopDate"]) + pd.Timedelta(hours=23) + pd.Timedelta(days=cfg0["LookAhead"]),
        freq="h",
    )

    demand_dir = Path(cfg0["SimulationDirectory"]).parent / "reserve_demands"
    demand_dir.mkdir(parents=True, exist_ok=True)

    files = {
        "aFRRUDemand": _write_reserve_profile(demand_dir / "aFRRUDemand.csv", base_index, 3.0),
        "aFRRDDemand": _write_reserve_profile(demand_dir / "aFRRDDemand.csv", base_index, 2.0),
        "FFRDemand": _write_reserve_profile(demand_dir / "FFRDemand.csv", base_index, 1.0),
        "FCRDemand": _write_reserve_profile(demand_dir / "FCRDemand.csv", base_index, 1.5),
    }

    cumulative_steps = [
        [],
        ["aFRRUDemand"],
        ["aFRRUDemand", "aFRRDDemand"],
        ["aFRRUDemand", "aFRRDDemand", "FFRDemand"],
        ["aFRRUDemand", "aFRRDDemand", "FFRDemand", "FCRDemand"],
    ]

    run_data = []
    for i, enabled in enumerate(cumulative_steps):
        cfg = load_test_config("tiny.yml", f"reserves_recursive_{i}")
        cfg["ReserveCalculation"] = "Exogenous"
        for k in ("aFRRUDemand", "aFRRDDemand", "FFRDemand", "FCRDemand"):
            cfg[k] = files[k] if k in enabled else ""

        out = build_solve(cfg)
        results = out["results"]

        run_data.append(
            {
                "enabled": tuple(enabled),
                "system_cost": float(results["OutputSystemCost"].sum().sum()),
                "battery_power": _battery_power_series(results),
                "battery_reserve_total": _battery_reserve_total(results),
            }
        )

    base = run_data[0]
    final = run_data[-1]

    # Reserve activation should not reduce total objective cost.
    assert final["system_cost"] >= base["system_cost"] - 1e-6

    # Reserve provision by battery should be effectively activated in cumulative scenario.
    assert final["battery_reserve_total"] > 0.0, "Battery reserve provision stayed zero"
    assert final["battery_reserve_total"] > base["battery_reserve_total"] + 1e-6, (
        "Battery reserve provision did not increase versus base scenario"
    )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
