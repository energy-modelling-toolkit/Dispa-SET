# -*- coding: utf-8 -*-
"""
Unit test: closed-loop storage units (batteries, closed-loop pumped storage)
cannot spill
=============================================================================

Context
-------
Dispa-SET only allows water spillage for storage units that have a natural
inflow (see ``dispaset/GAMS/UCM.gms``):

    sp(au) = s(au)$(smax(h, StorageInflow(au,h)) > 0);
    ...
    EQ_Storage_balance(au,i)..
        ...
        =E=
        StorageLevel(au,i)$(s(au))
        +StorageOutflow(au,i)$(s(au))*Nunits(au)$(s(au))*TimeStep
        +spillage(au,i)$sp(au)
        + ...

The ``spillage`` term only enters the storage balance for units in the
``sp`` subset. Batteries (BATS) and closed-loop pumped-storage hydro (HPHS)
units have no natural inflow, so ``StorageInflow(au,h)`` is identically zero
for them and they are therefore never members of ``sp``: structurally, they
cannot spill (their only loss paths are the round-trip (dis)charging
efficiency and ``StorageSelfDischarge``).

``StorageInflow`` is built in ``dispaset/preprocessing/build.py`` from the
``ReservoirScaledInflows`` ("ScaledInflows") table via ``UnitBasedTable``,
which falls back, in order, on a column named after the unit, then the
*Technology*, then the *Zone* (``fallbacks=['Unit', 'Technology', 'Zone']``).
If no column matches, the inflow defaults to 0. This test checks both:

1. (always runs, no GAMS needed) That ``ds.build_simulation`` assigns a
   ``StorageInflow`` of exactly 0 to BATS and HPHS units when the inflow
   table has no matching column for them (the realistic case, confirmed in
   a real database), while a reservoir with a matching
   ``Technology`` column receives the (nonzero) scaled inflow. This is the
   Python-side precondition that makes ``sp(au)`` false for these units in
   GAMS, i.e. it mirrors/validates the structural GAMS rule above.

2. (requires GAMS; skipped otherwise) A tiny 1-day LP solve with one
   battery, one closed-loop HPHS unit and one small HDAM reservoir forced
   to overflow (inflow = 5x its capacity): asserts ``OutputSpillage`` is
   exactly 0 for the battery and the HPHS unit, and strictly positive for
   the reservoir.

How to run
----------
1. ``pytest tests/unit/test_storage_spillage.py``
2. ``python tests/unit/test_storage_spillage.py``

Both tests are fast (a few seconds each); the GAMS-based one is skipped
automatically when no GAMS installation is detected.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from conftest import skip_if_no_gams  # noqa: E402
from _helpers import load_test_config  # noqa: E402

import dispaset as ds  # noqa: E402


# --------------------------------------------------------------------------- #
# 1. Build-level test (no GAMS): StorageInflow must be 0 for BATS/HPHS
# --------------------------------------------------------------------------- #
@pytest.mark.timeout(60)
def test_no_inflow_for_batteries_and_closed_loop_phs():
    """
    Build the tiny spillage-test config (Z1_BATS, Z1_HPHS, Z1_HDAM) and check
    the StorageInflow parameter that feeds the GAMS ``sp(au)`` rule:

    * Z1_BATS and Z1_HPHS have no matching column in ScaledInflows_spillage.csv
      (no Unit, Technology or Zone match) -> StorageInflow must be all-zero.
    * Z1_HDAM has a matching "HDAM" Technology column with a large constant
      value -> StorageInflow must be strictly positive at every timestep.

    Since GAMS sets ``sp(au) = s(au)$(smax(h, StorageInflow(au,h)) > 0)``,
    an all-zero StorageInflow for Z1_BATS/Z1_HPHS means they can never enter
    ``sp`` and the ``spillage(au,i)$sp(au)`` term can never apply to them.
    """
    cfg = load_test_config("tiny_spillage.yml", "unit_spillage_build")
    sim_data = ds.build_simulation(cfg)

    sets = sim_data["sets"]
    storage_units = list(sets["asu"])

    def find_unit(name: str) -> str:
        # LP clustering prefixes clustered unit names with "[n] - ", e.g.
        # "[1] - Z1_BATS", so match by suffix rather than exact equality.
        matches = [u for u in storage_units if u == name or u.endswith(" - " + name)]
        assert len(matches) == 1, f"Expected exactly one storage unit matching {name!r}, found {matches}"
        return matches[0]

    bats_name = find_unit("Z1_BATS")
    hphs_name = find_unit("Z1_HPHS")
    hdam_name = find_unit("Z1_HDAM")

    inflow = sim_data["parameters"]["StorageInflow"]["val"]  # shape (n_storage, n_h)

    def inflow_for(unit: str):
        i = storage_units.index(unit)
        return inflow[i, :]

    bats_inflow = inflow_for(bats_name)
    hphs_inflow = inflow_for(hphs_name)
    hdam_inflow = inflow_for(hdam_name)

    assert (bats_inflow == 0).all(), (
        "Z1_BATS (battery) must have zero StorageInflow at every timestep "
        "(no inflow data path -> cannot be in the GAMS sp(au) spillage set)."
    )
    assert (hphs_inflow == 0).all(), (
        "Z1_HPHS (closed-loop pumped storage) must have zero StorageInflow "
        "at every timestep (no inflow data path -> cannot be in sp(au))."
    )
    assert (hdam_inflow > 0).all(), (
        "Z1_HDAM (open reservoir) was given a nonzero 'HDAM' Technology "
        "column in ScaledInflows_spillage.csv and must receive it, so that "
        "it correctly becomes a member of sp(au) and can spill."
    )


def test_battery_inflow_from_fallback_is_removed(tmp_path):
    """
    A column named after the technology ("BATS") in the inflow table would reach
    the battery through the Unit -> Technology -> Zone fallback. Batteries have no
    natural inflow, so build.py must drop it (commons['tech_storage_no_inflow']),
    keeping the battery out of the spillage set sp(au).
    """
    import pandas as pd
    cfg = load_test_config("tiny_spillage.yml", "unit_spillage_fallback")
    src = Path(cfg["ReservoirScaledInflows"])
    if not src.is_absolute():
        src = Path(__file__).resolve().parents[1] / src
    df = pd.read_csv(src, index_col=0)
    df["BATS"] = 0.5
    path = tmp_path / "ScaledInflows_with_bats.csv"
    df.to_csv(path)
    cfg["ReservoirScaledInflows"] = str(path)
    sim_data = ds.build_simulation(cfg)
    storage_units = list(sim_data["sets"]["asu"])
    bats = [u for u in storage_units if u == "Z1_BATS" or u.endswith(" - Z1_BATS")][0]
    inflow = sim_data["parameters"]["StorageInflow"]["val"][storage_units.index(bats), :]
    assert (inflow == 0).all(), "A battery must never receive inflow, even through the Technology fallback."


# --------------------------------------------------------------------------- #
# 2. Mirror of the one-line GAMS rule (pure logic, documents the contract)
# --------------------------------------------------------------------------- #
def test_sp_set_rule_mirror():
    """
    Pure-Python mirror of the GAMS line:

        sp(au) = s(au)$(smax(h, StorageInflow(au,h)) > 0);

    This does not replace the GAMS behaviour (which is exercised for real by
    the build-level and solve-level tests above/below); it documents the
    exact boolean rule so a future change to this line is caught even by
    someone only reading the Python test suite.
    """
    import numpy as np

    def sp_rule(storage_inflow_by_unit: dict) -> set:
        return {u for u, series in storage_inflow_by_unit.items()
                if np.max(series) > 0}

    data = {
        "BATS1": np.zeros(24),
        "HPHS1": np.zeros(24),
        "HDAM1": np.full(24, 5.0),
    }
    sp = sp_rule(data)
    assert sp == {"HDAM1"}
    assert "BATS1" not in sp
    assert "HPHS1" not in sp


# --------------------------------------------------------------------------- #
# 3. Solve-level test (requires GAMS): OutputSpillage must be 0 for
#    BATS/HPHS and > 0 for the forced-overflow reservoir
# --------------------------------------------------------------------------- #
@pytest.mark.timeout(120)
def test_spillage_solve_bats_hphs_never_spill():
    """
    Tiny 1-day LP solve. Z1_HDAM has STOCapacity=5 MWh and a constant inflow
    of 5x its PowerCapacity with no pumping (STOMaxChargingPower=0), so it is
    forced to spill well above its generation + storage capability. Z1_BATS
    and Z1_HPHS have no inflow and must show exactly zero spillage.
    """
    skip_if_no_gams()

    cfg = load_test_config("tiny_spillage.yml", "unit_spillage_solve")
    ds.build_simulation(cfg)
    ds.solve_GAMS(cfg["SimulationDirectory"], cfg["GAMS_folder"])

    from pathlib import Path as _Path
    results_gdx = _Path(cfg["SimulationDirectory"]) / "Results.gdx"
    if not results_gdx.is_file():
        pytest.fail(f"GAMS did not produce Results.gdx in {cfg['SimulationDirectory']}")

    inputs, results = ds.get_sim_results(path=cfg["SimulationDirectory"])
    spillage = results["OutputSpillage"]

    def find_col(name: str):
        """
        Return the matching OutputSpillage column, or None if absent.

        Because ``spillage(au,i)`` only ever appears in the GAMS model
        inside the ``$sp(au)``-guarded term of EQ_Storage_balance, GAMS's
        sparse model generation never instantiates the variable for a unit
        that is not in ``sp`` (no natural inflow): there is no record for
        it at all in Results.gdx, i.e. the column is entirely absent from
        the OutputSpillage dataframe (stronger evidence than a value of 0,
        since the decision variable itself was never created for that
        unit).
        """
        matches = [c for c in spillage.columns if c == name or str(c).endswith(" - " + name)]
        assert len(matches) <= 1, f"Expected at most one spillage column matching {name!r}, found {matches}"
        return matches[0] if matches else None

    bats_col = find_col("Z1_BATS")
    hphs_col = find_col("Z1_HPHS")
    hdam_col = find_col("Z1_HDAM")

    if bats_col is None:
        pass  # no spillage variable was ever generated for the battery: strongest possible "cannot spill"
    else:
        assert (spillage[bats_col] == 0).all(), "Battery spillage must be exactly 0"

    if hphs_col is None:
        pass  # idem for the closed-loop pumped-storage unit
    else:
        assert (spillage[hphs_col] == 0).all(), "Closed-loop PHS spillage must be exactly 0"

    assert hdam_col is not None, (
        "Z1_HDAM (the forced-overflow reservoir) should be a member of sp(au) "
        "and have an OutputSpillage column"
    )
    assert (spillage[hdam_col] > 0).any(), (
        "The tiny, excess-inflow reservoir was expected to spill at least "
        "once over the horizon"
    )


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
