# -*- coding: utf-8 -*-
"""
Shared pytest fixtures and helpers for the Dispa-SET test suite.

The fixtures here keep tests short, isolated and reproducible:

* `dispaset_root`            : absolute path to the Dispa-SET repository root.
* `tests_dir`                : absolute path to the `tests/` directory.
* `output_dir`               : absolute path to a per-test output folder
                               (created on demand, removed at the end of
                               the session).
* `gams_path`                : auto-detected GAMS installation path. If GAMS
                               is not available, integration tests should be
                               skipped.

In addition, we set the `Agg` matplotlib backend at import time so that any
plot rendered during the tests is non-interactive (no blocking `plt.show`).

This file is automatically loaded by pytest. To run individual tests as
plain Python scripts (the user-facing requirement), each test file has its
own `if __name__ == "__main__"` block that performs the same imports.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import matplotlib

# Use a non-interactive backend everywhere: every plt.show() becomes a no-op,
# which is essential for fully-automated tests.
matplotlib.use("Agg")

import pytest


# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent
OUTPUT_DIR = TESTS_DIR / "output"

# Make sure the repo root is on sys.path so `import dispaset` always works
# regardless of where pytest is invoked from.
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


@pytest.fixture(scope="session")
def dispaset_root() -> Path:
    """Absolute path to the Dispa-SET repository root."""
    return REPO_ROOT


@pytest.fixture(scope="session")
def tests_dir() -> Path:
    """Absolute path to the tests directory."""
    return TESTS_DIR


@pytest.fixture(scope="session")
def output_dir() -> Path:
    """Absolute path to the output/ folder used for test artefacts."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUTPUT_DIR


@pytest.fixture(scope="session")
def gams_path() -> str | None:
    """
    Try to detect a working GAMS installation without prompting the user.
    Returns None if no GAMS is found - integration tests should call
    `pytest.skip(...)` in that case.
    """
    # 1. Honour user-set environment variables first
    for env_var in ("GAMSPATH", "GAMSDIR"):
        p = os.environ.get(env_var)
        if p and os.path.exists(os.path.join(p, "gams")):
            return p
        if p and os.path.exists(os.path.join(p, "gams.exe")):
            return p

    # 2. Try `which gams`
    from shutil import which
    g = which("gams")
    if g:
        return os.path.dirname(g)

    # 3. Common Linux/macOS install paths
    for base in ["/opt/gams", "/usr/local/gams",
                 os.path.expanduser("~/gams"),
                 os.path.expanduser("~/progs/GAMS")]:
        if os.path.isdir(base):
            for sub in sorted(os.listdir(base), reverse=True):
                cand = os.path.join(base, sub)
                if os.path.isfile(os.path.join(cand, "gams")) or \
                   os.path.isfile(os.path.join(cand, "gams.exe")):
                    return cand
    return None


# --------------------------------------------------------------------------- #
# Helpers used both by pytest fixtures and by the standalone __main__ blocks
# --------------------------------------------------------------------------- #
def shorten_horizon(config: dict, days: int = 1) -> dict:
    """Truncate the simulation period of a Dispa-SET config to `days` days."""
    import datetime as dt
    start = dt.datetime(*config["StartDate"])
    end = start + dt.timedelta(days=days)
    config["StopDate"] = (end.year, end.month, end.day, 0, 0, 0)
    config["LookAhead"] = 0
    config["HorizonLength"] = max(1, days)
    return config


def get_clean_simulation_dir(config: dict, suffix: str = "") -> str:
    """Return an absolute, empty SimulationDirectory under tests/output/."""
    base = OUTPUT_DIR / (Path(config["SimulationDirectory"]).name + suffix)
    if base.exists():
        shutil.rmtree(base)
    base.mkdir(parents=True, exist_ok=True)
    config["SimulationDirectory"] = str(base)
    return str(base)


def skip_if_no_gams():
    """Convenience: pytest.skip if no GAMS is detected."""
    g = None
    for env_var in ("GAMSPATH", "GAMSDIR"):
        if os.environ.get(env_var):
            g = os.environ[env_var]
            break
    if not g:
        from shutil import which
        if which("gams"):
            g = "found"
    if not g:
        pytest.skip("GAMS installation not detected, skipping integration test")
