# -*- coding: utf-8 -*-
"""
Failure tests: insufficient generation capacity
================================================

What these tests do
-------------------
These tests deliberately build and solve configurations where installed
generation capacity is far too low to serve all demand.  When that happens:

* GAMS activates ``LostLoad_MaxPower`` (involuntary load shedding).
* ``check_energy_balance`` detects the resulting energy imbalance and emits
  a ``logging.CRITICAL`` message.

The tests **pass** when both of those things happen correctly, confirming
that the framework's diagnostic layer reliably catches infeasible or
severely under-resourced dispatch cases.

Contrast with the integration tests in ``tests/integration/``, which are
expected to produce *clean* solves (zero CRITICAL messages).

How to run
----------

1. ``pytest tests/failure/test_insufficient_generation.py``
2. ``python tests/failure/test_insufficient_generation.py``

Skipped automatically if no GAMS installation is detected.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import build_solve, load_test_config, skip_if_no_gams  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
# Demand multiplier that guarantees infeasibility.  The tiny test cases have
# ≈ 700 MW installed in Z1 and ≈ 450 MW in Z2.  A 20× multiplier drives
# peak load to many GW, making LostLoad_MaxPower unavoidable.
_OVERLOAD_FACTOR = 20.0


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------
class _CriticalLogCatcher(logging.Handler):
    """Collects every CRITICAL-level log record emitted while attached."""

    def __init__(self):
        super().__init__(level=logging.CRITICAL)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def _build_and_capture_critical(cfg: dict) -> tuple[dict, list[logging.LogRecord]]:
    """Run build_solve and return (out, critical_records)."""
    handler = _CriticalLogCatcher()
    logging.getLogger().addHandler(handler)
    try:
        out = build_solve(cfg)
    finally:
        logging.getLogger().removeHandler(handler)
    return out, handler.records


def _lostload_sum(results: dict, zone: str) -> float:
    """Total LostLoad_MaxPower [MWh] for *zone* in these results."""
    ll_df = results.get("LostLoad_MaxPower")
    if ll_df is None or zone not in ll_df.columns:
        return 0.0
    return float(ll_df[zone].sum())


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
@pytest.mark.timeout(180)
def test_single_zone_lostload_triggers_critical():
    """
    Single-zone LP: demand 20× normal → LostLoad_MaxPower > 0 in Z1, and
    ``check_energy_balance`` must emit at least one CRITICAL message.
    """
    skip_if_no_gams()

    cfg = load_test_config("tiny.yml", "insufficient_gen_lp_z1")
    cfg["modifiers"]["Demand"] = _OVERLOAD_FACTOR

    out, crits = _build_and_capture_critical(cfg)
    results = out["results"]
    balance = out["energy_balance"]

    # 1. LostLoad must be positive — generation genuinely fell short
    ll = _lostload_sum(results, "Z1")
    assert ll > 0, (
        f"Expected positive LostLoad_MaxPower in Z1 with demand×{_OVERLOAD_FACTOR}, "
        f"got {ll} MWh"
    )

    # 2. check_energy_balance must report > 1 % imbalance for Z1
    assert balance["Z1"] > 0.01, (
        f"Expected energy-balance error > 1 % for Z1, got {balance['Z1']:.4%}"
    )

    # 3. At least one CRITICAL log message must mention energy balance
    assert crits, "Expected at least one CRITICAL log message but none were captured"
    messages = " ".join(r.getMessage() for r in crits)
    assert "energy balance" in messages.lower(), (
        f"Expected 'energy balance' in a CRITICAL message; got:\n{messages}"
    )


@pytest.mark.timeout(180)
def test_two_zone_lostload_triggers_critical_per_zone():
    """
    Two-zone LP: demand 20× normal → LostLoad in both Z1 and Z2, and
    ``check_energy_balance`` must emit a CRITICAL message for each zone.
    """
    skip_if_no_gams()

    cfg = load_test_config("tiny_2zones.yml", "insufficient_gen_lp_2z")
    cfg["modifiers"]["Demand"] = _OVERLOAD_FACTOR

    out, crits = _build_and_capture_critical(cfg)
    results = out["results"]
    balance = out["energy_balance"]

    for zone in cfg["zones"]:
        ll = _lostload_sum(results, zone)
        assert ll > 0, (
            f"Expected positive LostLoad_MaxPower in {zone} with demand×{_OVERLOAD_FACTOR}, "
            f"got {ll} MWh"
        )
        assert balance[zone] > 0.01, (
            f"Expected energy-balance error > 1 % for {zone}, "
            f"got {balance[zone]:.4%}"
        )

    # At least one CRITICAL per zone
    assert len(crits) >= len(cfg["zones"]), (
        f"Expected >= {len(cfg['zones'])} CRITICAL messages (one per zone), "
        f"got {len(crits)}"
    )


@pytest.mark.timeout(180)
def test_milp_lostload_triggers_critical():
    """
    MILP formulation: demand 20× normal → same diagnostic behaviour as LP.
    Validates that the MILP commitment constraints do not prevent LostLoad
    from being reported.
    """
    skip_if_no_gams()

    cfg = load_test_config("tiny_milp.yml", "insufficient_gen_milp")
    cfg["modifiers"]["Demand"] = _OVERLOAD_FACTOR

    out, crits = _build_and_capture_critical(cfg)
    results = out["results"]
    balance = out["energy_balance"]

    ll = _lostload_sum(results, "Z1")
    assert ll > 0, (
        f"Expected positive LostLoad_MaxPower in Z1 (MILP), got {ll} MWh"
    )
    assert balance["Z1"] > 0.01, (
        f"Expected energy-balance error > 1 % for Z1 (MILP), "
        f"got {balance['Z1']:.4%}"
    )
    assert crits, "Expected at least one CRITICAL log message (MILP)"


@pytest.mark.timeout(180)
def test_balance_error_bounded_by_one():
    """
    Sanity check: even in a completely infeasible case the relative
    energy-balance error returned by ``check_energy_balance`` must stay
    strictly below 1.0 (100 %).

    A value ≥ 1.0 would indicate a double-counting bug (e.g. adding
    shifted load instead of subtracting it, or including LostLoad twice).
    """
    skip_if_no_gams()

    cfg = load_test_config("tiny.yml", "insufficient_gen_bounded")
    cfg["modifiers"]["Demand"] = _OVERLOAD_FACTOR

    out, crits = _build_and_capture_critical(cfg)
    balance = out["energy_balance"]

    # Must be non-zero (infeasible) but also strictly bounded
    assert balance["Z1"] > 0.01, (
        f"Expected imbalance > 1 % in Z1, got {balance['Z1']:.4%}"
    )
    assert balance["Z1"] < 1.0, (
        f"Imbalance for Z1 is {balance['Z1']:.2%} — this is suspiciously high "
        "and may indicate a double-counting bug in check_energy_balance"
    )
    assert crits, "Expected CRITICAL log message about energy balance"


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
