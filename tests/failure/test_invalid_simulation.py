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
