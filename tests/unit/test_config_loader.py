# -*- coding: utf-8 -*-
"""
Unit test: Dispa-SET configuration loader
=========================================

What this test does
-------------------
The loader (`dispaset.preprocessing.data_handler.load_config`) accepts both
`.xlsx` and `.yml` configuration files. Several invariants must hold:

* All four "tiny" YAML configurations under `tests/configs/` load to a
  ``dict`` with at least the keys ``zones``, ``StartDate``, ``StopDate``,
  ``SimulationType`` and ``default``.
* All ``Path``-shaped configuration values are converted to absolute paths
  (or to the empty string when not provided).
* A YAML file that uses Excel-style configurations (the legacy `.xlsx`)
  also loads correctly.
* A wrong file extension (``.ymx``) makes the loader exit cleanly.
* A garbled Excel file raises ``SystemExit`` instead of crashing.

How to run
----------

1. ``pytest tests/unit/test_config_loader.py``
2. ``python tests/unit/test_config_loader.py``

This test does NOT require GAMS.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import dispaset as ds  # noqa: E402

CONFIG_DIR = Path(__file__).resolve().parents[1] / "configs"


# --------------------------------------------------------------------------- #
# Successful loads
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("yml", [
    "tiny.yml", "tiny_2zones.yml", "tiny_milp.yml", "tiny_bs.yml",
    "tiny_mts.yml", "ultimate.yml",
])
def test_yaml_configs_load(yml):
    """Every tiny YAML config must load to a ``dict`` with required keys."""
    cfg = ds.load_config(str(CONFIG_DIR / yml))
    assert isinstance(cfg, dict)
    for key in ("zones", "StartDate", "StopDate", "SimulationType", "default"):
        assert key in cfg, f"Missing key {key!r} after loading {yml}"


def test_paths_are_absolute_or_empty():
    """Every path-typed value in tiny.yml is either '' or absolute."""
    cfg = ds.load_config(str(CONFIG_DIR / "tiny.yml"))
    path_keys = ["Demand", "RenewablesAF", "PowerPlantData", "NTC",
                 "Interconnections", "ReservoirLevels",
                 "ReservoirScaledInflows", "GeoData"]
    for key in path_keys:
        v = cfg.get(key, "")
        assert v == "" or os.path.isabs(v), \
            f"Config[{key!r}] is not absolute and not empty: {v!r}"


def test_legacy_xlsx_config_loads():
    """The legacy .xlsx config should still load (backwards compat)."""
    xlsx = CONFIG_DIR / "conf.xlsx"
    if not xlsx.exists():
        pytest.skip(f"{xlsx} not present")
    cfg = ds.load_config(str(xlsx))
    assert isinstance(cfg, dict)
    assert "zones" in cfg


# --------------------------------------------------------------------------- #
# Failure modes
# --------------------------------------------------------------------------- #
def test_wrong_extension_exits():
    """A `.ymx` extension must trigger a SystemExit."""
    bad = CONFIG_DIR / "conf.ymx"
    if not bad.exists():
        pytest.skip(f"{bad} not present")
    with pytest.raises(SystemExit):
        ds.load_config(str(bad))


def test_wrong_xlsx_title_exits():
    """An xlsx with the wrong root title must trigger a SystemExit."""
    bad = CONFIG_DIR / "config_wrong_title.xlsx"
    if not bad.exists():
        pytest.skip(f"{bad} not present")
    with pytest.raises(SystemExit):
        ds.load_config(str(bad))


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
