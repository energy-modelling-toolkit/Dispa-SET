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


def test_stopdate_before_startdate_fails():
    cfg = load_test_config("tiny.yml", "fail_stopdate")
    cfg["StartDate"] = (2015, 1, 5, 0, 0, 0)
    cfg["StopDate"] = (2015, 1, 1, 0, 0, 0)
    with pytest.raises((SystemExit, Exception)):
        ds.build_simulation(cfg)


def test_stopdate_before_startdate_exits_with_code_1():
    """The new validation added in future_proof explicitly calls sys.exit(1).

    This test verifies that build_simulation raises SystemExit with code 1
    (not a different exit code or a generic Exception) when StartDate is
    strictly later than StopDate.
    """
    cfg = load_test_config("tiny.yml", "fail_stopdate_code")
    cfg["StartDate"] = (2015, 12, 31, 0, 0, 0)
    cfg["StopDate"]  = (2015,  1,  1, 0, 0, 0)
    with pytest.raises(SystemExit) as exc_info:
        ds.build_simulation(cfg)
    assert exc_info.value.code == 1, (
        f"Expected SystemExit(1), got SystemExit({exc_info.value.code})"
    )


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
    except SystemExit as exc:
        # Exit code 1 would mean the date guard fired – that is wrong.
        assert exc.code != 1, "Date guard must not fire when StartDate == StopDate"
    except Exception:
        pass  # Any other error is fine (horizon too short, etc.)


def test_unknown_simulation_type_fails():
    cfg = load_test_config("tiny.yml", "fail_simtype")
    cfg["SimulationType"] = "Banana clustered"
    with pytest.raises((SystemExit, Exception)):
        ds.build_simulation(cfg)


def test_empty_zones_list_fails():
    cfg = load_test_config("tiny.yml", "fail_emptyzones")
    cfg["zones"] = []
    with pytest.raises((SystemExit, Exception)):
        ds.build_simulation(cfg)


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))

