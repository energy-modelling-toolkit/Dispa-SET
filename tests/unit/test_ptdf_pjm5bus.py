# -*- coding: utf-8 -*-
"""
Unit tests: PTDF matrix for the PJM 5-bus test system
=======================================================

What this test does
-------------------
Verifies that ``PTDF_matrix`` correctly computes the Power Transfer Distribution
Factor matrix for the standard PJM 5-bus academic test case from:

    F. Li and R. Bo, "Small test systems for power system economic studies,"
    IEEE PES General Meeting, Minneapolis, MN, 2010, pp. 1-4.
    doi: 10.1109/PES.2010.5589973

The GAMS reference implementation is:

    A. Soroudi, "DC Optimal Power Flow (OPF) model",
    https://github.com/OptimizationExpert/GAMS-DC--OPF-
    Contributed by Dr. Alireza Soroudi, University College Dublin, Ireland.

Network topology (``tests/data/GridData_pjm5bus.csv``)::

    PJM_B1 (slack, most connections: 3)
      |── X=0.0281, Fmax=400  ──► PJM_B2 (300 MW load)
      |── X=0.0304, Fmax=1000 ──► PJM_B4 (400 MW load + Sundance)
      └── X=0.0064, Fmax=1000 ──► PJM_B5 (Brighton generator)
    PJM_B2 ── X=0.0108, Fmax=1000 ──► PJM_B3 (300 MW load + Solitude)
    PJM_B3 ── X=0.0297, Fmax=1000 ──► PJM_B4
    PJM_B4 ── X=0.0297, Fmax=240  ──► PJM_B5

Dispa-SET automatically selects PJM_B1 as the slack bus (highest connectivity).
The PTDF matrix has shape (6 lines, 5 buses), with the PJM_B1 column all-zero.

PTDF reference values (computed analytically from the network reactances):

    Line              | PJM_B1 | PJM_B2  | PJM_B3  | PJM_B4  | PJM_B5
    B1 -> B2 (400 MW) |  0.000 | -0.6698 | -0.5429 | -0.1939 | -0.0344
    B1 -> B4 (1000MW) |  0.000 | -0.1792 | -0.2481 | -0.4376 | -0.0776
    B1 -> B5 (1000MW) |  0.000 | -0.1509 | -0.2090 | -0.3685 | -0.8880
    B2 -> B3 (1000MW) |  0.000 | +0.3302 | -0.5429 | -0.1939 | -0.0344
    B3 -> B4 (1000MW) |  0.000 | +0.3302 | +0.4571 | -0.1939 | -0.0344
    B4 -> B5 (240 MW) |  0.000 | +0.1509 | +0.2090 | +0.3685 | -0.1120

Key DC power flow physics:
  - Brighton (B5, cheapest at 10 $/MWh) has the largest PTDF on B1-B5 (-0.888),
    meaning its power is mostly exported via the direct B1-B5 line.
  - Loading B5 forces 88.8% of the flow back through the high-susceptance
    B1-B5 line, with 11.2% looping via B4-B5.
  - The B1-B2 line (X=0.0281, highest susceptance) carries the most per-MW
    injected anywhere, showing it is the main power highway.

Tests
-----
1. ``test_ptdf_pjm5bus_shape``        — matrix shape is (6, 5)
2. ``test_ptdf_pjm5bus_slack_zero``   — PJM_B1 column is all-zero
3. ``test_ptdf_pjm5bus_key_entries``  — 8 specific entries within 1e-4
4. ``test_ptdf_pjm5bus_b5_injection`` — B5 injection flows mainly on B1-B5
5. ``test_build_pjm5bus_nonempty``    — build_simulation PTDF non-zero

No GAMS licence required for any of these tests.

How to run
----------
1. ``pytest tests/unit/test_ptdf_pjm5bus.py``
2. ``python tests/unit/test_ptdf_pjm5bus.py``
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _helpers import TESTS_DIR, load_test_config  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GRID_CSV = TESTS_DIR / "data" / "GridData_pjm5bus.csv"

ZONES = ["PJM_B1", "PJM_B2", "PJM_B3", "PJM_B4", "PJM_B5"]
LINES = [
    "PJM_B1 -> PJM_B2",
    "PJM_B1 -> PJM_B4",
    "PJM_B1 -> PJM_B5",
    "PJM_B2 -> PJM_B3",
    "PJM_B3 -> PJM_B4",
    "PJM_B4 -> PJM_B5",
]

# Reference PTDF values (tolerance 1e-4) derived analytically from the
# PJM 5-bus network reactances (Sbase=100 MVA, Li & Bo 2010 Table I).
# Column order: PJM_B1(slack), PJM_B2, PJM_B3, PJM_B4, PJM_B5
EXPECTED_PTDF = {
    # format: {line_name: {zone_name: expected_value}}
    "PJM_B1 -> PJM_B2": {
        "PJM_B1": 0.0,
        "PJM_B2": -0.669811,
        "PJM_B3": -0.542906,
        "PJM_B4": -0.193917,
        "PJM_B5": -0.034379,
    },
    "PJM_B1 -> PJM_B5": {
        "PJM_B1": 0.0,
        "PJM_B2": -0.150943,
        "PJM_B3": -0.208957,
        "PJM_B4": -0.368495,
        "PJM_B5": -0.888043,
    },
    "PJM_B4 -> PJM_B5": {
        "PJM_B1": 0.0,
        "PJM_B2": 0.150943,
        "PJM_B3": 0.208957,
        "PJM_B4": 0.368495,
        "PJM_B5": -0.111957,
    },
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_grid() -> pd.DataFrame:
    return pd.read_csv(GRID_CSV, index_col=0, na_values=["", " "], keep_default_na=False)


def _make_config() -> dict:
    return {"zones": ZONES}


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.timeout(5)
def test_ptdf_pjm5bus_shape():
    """PTDF matrix has shape (n_lines=6, n_zones=5) for the PJM 5-bus grid."""
    from dispaset.preprocessing.utils import PTDF_matrix

    ptdf = PTDF_matrix(_make_config(), _load_grid())

    assert ptdf.shape == (6, 5), f"Expected (6, 5), got {ptdf.shape}"
    assert list(ptdf.columns) == ZONES, "Column order does not match zone list"
    assert list(ptdf.index) == LINES, "Row index does not match GridData lines"


@pytest.mark.timeout(5)
def test_ptdf_pjm5bus_slack_zero():
    """The PJM_B1 column (slack bus) is identically zero in the PTDF matrix.

    Dispa-SET's PTDF_matrix selects the slack bus automatically as the node
    with the most connections. PJM_B1 has three lines (B1-B2, B1-B4, B1-B5)
    while every other bus has at most two, so B1 is always chosen as slack.
    """
    from dispaset.preprocessing.utils import PTDF_matrix

    ptdf = PTDF_matrix(_make_config(), _load_grid()).fillna(0)

    assert np.allclose(ptdf["PJM_B1"], 0.0), (
        "Slack-bus column PJM_B1 is not zero:\n" + str(ptdf["PJM_B1"])
    )


@pytest.mark.timeout(5)
def test_ptdf_pjm5bus_key_entries():
    """Specific PTDF entries match the analytically derived reference (tol=1e-4).

    Reference values are derived from the PJM 5-bus network data of Li & Bo
    (2010) using standard DC power flow analysis:
      B = A * Bp * A^T  (bus susceptance matrix)
      PTDF = Bp * A_red * inv(B_red)
    where the reduced matrices drop the slack-bus row and column.
    """
    from dispaset.preprocessing.utils import PTDF_matrix

    ptdf = PTDF_matrix(_make_config(), _load_grid()).fillna(0)

    tol = 1e-4
    for line, zone_vals in EXPECTED_PTDF.items():
        for zone, expected in zone_vals.items():
            actual = float(ptdf.loc[line, zone])
            assert abs(actual - expected) < tol, (
                f"PTDF[{line}, {zone}]: expected {expected:.6f}, got {actual:.6f}"
            )


@pytest.mark.timeout(5)
def test_ptdf_pjm5bus_b5_injection():
    """B5 injection (Brighton) flows predominantly on the direct B1-B5 line.

    Physical interpretation: Brighton at bus 5 exports most power directly
    to bus 1 (via the low-reactance X=0.0064 line), with only 11% looping
    through the B4-B5 line.  This PTDF ratio is a key indicator that the
    B1-B5 path dominates compared to the B4-B5 path.

    The scenario injects 1 MW at PJM_B5, absorbed by the slack PJM_B1:
      F[B1->B5]  = -0.888  (88.8% flows directly B1<-B5, negative = counter to
                            defined direction means power from B5 to B1)
      F[B4->B5]  = -0.112  (11.2% loops via B4-B5 backward)
    The B1-B5 path carries |0.888 / 0.112| ≈ 7.9× more power than B4-B5.
    """
    from dispaset.preprocessing.utils import PTDF_matrix

    ptdf = PTDF_matrix(_make_config(), _load_grid()).fillna(0)

    ptdf_b1b5_b5 = float(ptdf.loc["PJM_B1 -> PJM_B5", "PJM_B5"])   # ≈ -0.888
    ptdf_b4b5_b5 = float(ptdf.loc["PJM_B4 -> PJM_B5", "PJM_B5"])   # ≈ -0.112

    assert abs(ptdf_b1b5_b5) > abs(ptdf_b4b5_b5), (
        "B1-B5 direct path should carry more power per MW from B5 than B4-B5, "
        f"but PTDF[B1->B5,B5]={ptdf_b1b5_b5:.4f} is not larger in magnitude "
        f"than PTDF[B4->B5,B5]={ptdf_b4b5_b5:.4f}"
    )
    # Specifically, about 7-8x more
    ratio = abs(ptdf_b1b5_b5) / (abs(ptdf_b4b5_b5) + 1e-9)
    assert ratio > 5, (
        f"Expected B1-B5 / B4-B5 ratio > 5, got {ratio:.2f}. "
        "The low-reactance B1-B5 line (X=0.0064) should dominate."
    )


@pytest.mark.timeout(30)
def test_build_pjm5bus_nonempty():
    """build_simulation with DC-Power-Flow config produces a non-zero PTDF.

    Builds the 5-bus PJM simulation environment (no GAMS solve) and checks
    that the ``PTDF`` parameter in the returned SimData is populated with
    non-zero values.  This verifies that the PTDF was correctly injected into
    the GDX inputs and that the full build pipeline handles 5 zones correctly.
    """
    import dispaset as ds

    cfg = load_test_config("tiny_pjm5bus_dcpf.yml", "build_pjm5bus")
    sim_data = ds.build_simulation(cfg)

    ptdf_val = sim_data["parameters"]["PTDF"]["val"]
    assert ptdf_val.size > 0, "PTDF parameter array is empty after build"
    assert np.any(ptdf_val != 0), (
        "PTDF parameter is all-zero in DC-Power-Flow build — "
        "the PTDF was not injected into the simulation inputs"
    )
    # Expect 6 lines × 5 buses = 30 entries, and at least 24 non-zero
    # (the 6 PJM_B1=slack column entries are zero, leaving ≥24 non-zero)
    n_nonzero = int(np.count_nonzero(ptdf_val))
    assert n_nonzero >= 24, (
        f"Expected at least 24 non-zero PTDF entries (6 lines × 4 non-slack "
        f"buses), got {n_nonzero}"
    )
    print(f"PTDF shape: {ptdf_val.shape}, non-zero entries: {n_nonzero}")


# ---------------------------------------------------------------------------
# Script entry-point (allows: python tests/unit/test_ptdf_pjm5bus.py)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import pytest as _pt

    _pt.main([__file__, "-v"])
