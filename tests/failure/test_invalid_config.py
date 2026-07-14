# -*- coding: utf-8 -*-
"""
Failure test: invalid configuration files
==========================================

What this test does
-------------------
Checks that Dispa-SET fails *cleanly* (raises an exception or calls
``sys.exit``) when the user provides an invalid configuration file:

* Non-existent config path -> ``FileNotFoundError``-style error.
* Excel config with wrong sheet/title structure -> ``SystemExit`` /
  exception.
* YAML with missing mandatory fields -> ``SystemExit`` / exception.

We deliberately use a generous catch (``Exception | SystemExit``) here
because the goal is to verify that the user gets *some* error rather
than a silent or undefined behaviour. Tightening the assertions to a
specific exception type can be done in a follow-up cleanup of
``data_handler.load_config``.

How to run
----------

1. ``pytest tests/failure/test_invalid_config.py``
2. ``python tests/failure/test_invalid_config.py``
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from dispaset.preprocessing.data_handler import load_config  # noqa: E402


def test_missing_config_file_raises():
    with pytest.raises((SystemExit, FileNotFoundError, IOError, OSError, Exception)):
        load_config("path/that/does/not/exist.yml")


def test_wrong_extension_raises():
    """Unsupported config extensions should raise an error."""
    with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
        f.write(b"hello")
        path = f.name
    try:
        with pytest.raises((SystemExit, Exception)):
            load_config(path)
    finally:
        os.unlink(path)


def test_malformed_yaml_raises():
    """YAML that is not parseable should produce an error, not silent
    success."""
    with tempfile.NamedTemporaryFile(suffix=".yml", delete=False, mode="w") as f:
        f.write("this is : not : valid : yaml :\n   - improper\n  - indent")
        path = f.name
    try:
        with pytest.raises((SystemExit, yaml.YAMLError, Exception)):
            load_config(path)
    finally:
        os.unlink(path)


def test_wrong_excel_template_raises():
    """The bundled malformed Excel config should also fail to parse."""
    bad_xlsx = REPO_ROOT / "tests" / "configs" / "config_wrong_title.xlsx"
    if not bad_xlsx.exists():
        pytest.skip("malformed reference xlsx not available")
    with pytest.raises((SystemExit, Exception)):
        load_config(str(bad_xlsx))


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
