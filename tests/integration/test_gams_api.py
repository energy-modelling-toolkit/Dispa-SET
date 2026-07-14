# -*- coding: utf-8 -*-
"""
Integration test: GAMS / GDX low-level API
============================================

What this test does
-------------------
Verifies the low-level interaction between Dispa-SET and the GAMS Python
API:

* ``package_exists('gams')`` returns True (the package is installed).
* ``get_gams_path()`` returns a directory containing a ``gams`` executable.
* A trivial GDX file can be written using ``write_variables`` and read
  back using ``gdx_to_list`` + ``gdx_to_dataframe``.

This test isolates problems with the GAMS link from problems with the
Dispa-SET model itself.

How to run
----------

1. ``pytest tests/integration/test_gams_api.py``
2. ``python tests/integration/test_gams_api.py``

Skipped automatically if no GAMS installation is detected.
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import skip_if_no_gams  # noqa: E402

from dispaset.misc.gdx_handler import (  # noqa: E402
    package_exists, get_gams_path, write_variables,
    gdx_to_list, gdx_to_dataframe,
)


def test_gams_package_available():
    skip_if_no_gams()
    assert package_exists("gams")


def test_gams_path_resolves():
    skip_if_no_gams()
    p = get_gams_path()
    assert p is not None
    assert os.path.isdir(p)
    # Linux/macOS: the executable must exist
    if os.name != "nt":
        assert os.path.isfile(os.path.join(p, "gams"))


@pytest.mark.timeout(30)
def test_gdx_round_trip():
    """Write a trivial GDX file and read it back."""
    skip_if_no_gams()
    sets = {"i": ["a", "b", "c"]}
    parameters = {
        "p": {
            "sets": ["i"],
            "val": np.array([1.0, 2.0, 3.0]),
        }
    }
    config = {"GAMS_folder": get_gams_path()}
    with tempfile.TemporaryDirectory() as tmp:
        gdx_path = os.path.join(tmp, "trip.gdx")
        write_variables(config, gdx_path, [sets, parameters])
        assert os.path.isfile(gdx_path)

        raw = gdx_to_list(get_gams_path(), gdx_path, varname="all")
        assert "p" in raw
        assert "i" in raw

        dfs = gdx_to_dataframe(raw, fixindex=False)
        assert "p" in dfs


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
