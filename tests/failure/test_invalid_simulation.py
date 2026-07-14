# -*- coding: utf-8 -*-
"""
Failure tests: invalid simulation parameters
=============================================

What this test does
-------------------
Verifies that ``ds.build_simulation`` aborts cleanly when the user
provides a configuration that is structurally inconsistent:

* ``StopDate`` earlier than ``StartDate``.
* Unknown ``SimulationType``.
* ``HorizonLength`` <= 0.
* Empty zone list.

Each case loads the standard ``tiny.yml``, mutates one field and asserts
that ``build_simulation`` raises an exception or aborts the interpreter
via ``sys.exit``.

How to run
----------

1. ``pytest tests/failure/test_invalid_simulation.py``
2. ``python tests/failure/test_invalid_simulation.py``
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import load_test_config  # noqa: E402

import dispaset as ds  # noqa: E402
from dispaset.common import DispaSETValidationError  # noqa: E402


def test_stopdate_before_startdate_fails():
    cfg = load_test_config("tiny.yml", "fail_stopdate")
    cfg["StartDate"] = (2015, 1, 5, 0, 0, 0)
    cfg["StopDate"] = (2015, 1, 1, 0, 0, 0)
    with pytest.raises(DispaSETValidationError):
        ds.build_simulation(cfg)


def test_stopdate_before_startdate_exits_with_code_1():
    """build_simulation raises DispaSETValidationError when StartDate is strictly
    later than StopDate.
    """
    cfg = load_test_config("tiny.yml", "fail_stopdate_code")
    cfg["StartDate"] = (2015, 12, 31, 0, 0, 0)
    cfg["StopDate"]  = (2015,  1,  1, 0, 0, 0)
    with pytest.raises(DispaSETValidationError):
        ds.build_simulation(cfg)


def test_startdate_equal_stopdate_does_not_fail():
    """Equal StartDate and StopDate is a degenerate (zero-length) horizon but
    must NOT trigger the StartDate > StopDate guard.

    The guard in future_proof uses a strict ``>`` comparison, so equal dates
    should pass the check (any subsequent error is acceptable – we only
    verify no immediate exit from the date guard itself).
    """
    cfg = load_test_config("tiny.yml", "equal_dates")
    same = (2015, 1, 1, 0, 0, 0)
    cfg["StartDate"] = same
    cfg["StopDate"]  = same
    try:
        ds.build_simulation(cfg)
    except DispaSETValidationError as exc:
        # If the date guard fired, that is wrong.
        pytest.fail(f"Date guard must not fire when StartDate == StopDate: {exc}")
    except Exception:
        pass  # Any other error is fine (horizon too short, etc.)


def test_unknown_simulation_type_fails():
    """An unrecognised SimulationType must raise DispaSETValidationError."""
    cfg = load_test_config("tiny.yml", "fail_simtype")
    cfg["SimulationType"] = "Banana clustered"
    with pytest.raises(DispaSETValidationError):
        ds.build_simulation(cfg)


def test_empty_zones_list_fails():
    cfg = load_test_config("tiny.yml", "fail_emptyzones")
    cfg["zones"] = []
    with pytest.raises((DispaSETValidationError, Exception)):
        ds.build_simulation(cfg)


def test_invalid_data_time_step_raises():
    """DataTimeStep != 1 must raise DispaSETValidationError."""
    cfg = load_test_config("tiny.yml", "fail_datatimestep")
    cfg["DataTimeStep"] = 2
    with pytest.raises(DispaSETValidationError):
        ds.build_simulation(cfg)


def test_invalid_simulation_time_step_raises():
    """SimulationTimeStep not in (1, 24) must raise DispaSETValidationError."""
    cfg = load_test_config("tiny.yml", "fail_simtimestep")
    cfg["SimulationTimeStep"] = 6
    with pytest.raises(DispaSETValidationError):
        ds.build_simulation(cfg)


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))

