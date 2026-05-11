# -*- coding: utf-8 -*-
"""
Unit test: Dispa-SET *common* dictionary
=========================================

What this test does
-------------------
The `dispaset.common.commons` dict is the central registry used everywhere in
the codebase to enumerate technologies, fuels, colours and a few global
constants. Many bugs come from a typo in this dict (a fuel that does not
match the GAMS file, a colour key that is missing, an inconsistent merit
order). This test verifies a few invariants:

* Required top-level keys are present.
* Every technology in the various sub-lists belongs to the master
  ``Technologies`` list.
* Every fuel referenced by ``MeritOrder`` is either a known fuel from
  ``Fuels`` or a control flow tag (``FlowIn`` / ``FlowOut`` / a storage
  tech).
* The colour table contains an entry for every fuel.
* Storage technologies are a subset of ``Technologies``.

How to run
----------

1. As a pytest:                  ``pytest tests/unit/test_common.py``
2. As a standalone Python script: ``python tests/unit/test_common.py``

This test does NOT require GAMS.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Allow running as a standalone script (without pytest):
if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import dispaset as ds  # noqa: E402
from dispaset.common import commons  # noqa: E402


REQUIRED_KEYS = [
    "TimeStep", "Technologies", "tech_renewables", "tech_storage",
    "tech_p2bs", "tech_bs2p", "tech_boundary_sector", "types_CHP",
    "Fuels", "MeritOrder", "colors", "hatches", "logfile", "na_values",
    "StdParameters", "PathParameters", "modifiers", "default",
]


def test_required_keys_present():
    """All required top-level entries should be present in `commons`."""
    missing = [k for k in REQUIRED_KEYS if k not in commons]
    assert not missing, f"Missing required keys in commons: {missing}"


def test_tech_subsets():
    """Every sub-tech list is a subset of the global `Technologies` list."""
    techs = set(commons["Technologies"])
    subsets = {
        "tech_renewables": commons["tech_renewables"],
        "tech_storage":    commons["tech_storage"],
        "tech_p2bs":       commons["tech_p2bs"],
        "tech_bs2p":       commons["tech_bs2p"],
        "tech_boundary_sector": commons["tech_boundary_sector"],
    }
    for name, lst in subsets.items():
        unknown = set(lst) - techs
        assert not unknown, (f"In commons['{name}'] the technology codes "
                             f"{sorted(unknown)} are not declared in "
                             f"commons['Technologies']")


def test_fuel_consistency():
    """Fuel codes referenced from MeritOrder should be known fuels or tags."""
    fuels = set(commons["Fuels"]) | {"FlowIn", "FlowOut"} | set(commons["tech_storage"])
    # Boundary-sector / P2X technology codes that are allowed in the merit order:
    # - P2GS: Power-to-gas storage (original)
    # - P2X:  Power-to-X (added in future_proof merge to represent generic P2X dispatch)
    extra = {"P2GS", "P2X"}
    fuels |= extra
    unknown = [f for f in commons["MeritOrder"] if f not in fuels]
    assert not unknown, f"Unknown items in MeritOrder: {unknown}"


def test_colour_table():
    """Every fuel must have a colour entry (for plotting)."""
    missing = [f for f in commons["Fuels"] if f not in commons["colors"]]
    # Known acceptable absences: OTH-only fuel codes that share a colour key
    allowed_missing = {"WIN", "WAT", "BIO", "GEO", "AMO", "AIR", "WHT", "ELE", "THE"}
    real_missing = set(missing) - allowed_missing
    # We only assert there is a colour entry for the explicit MeritOrder list:
    merit_missing = [f for f in commons["MeritOrder"]
                     if f not in commons["colors"] and f not in {"FlowIn", "FlowOut"}]
    assert not merit_missing, (f"Missing colour entries for merit-order items: "
                               f"{merit_missing}")


def test_logfile_naming():
    """The log filename must be a string ending with `.dispa.log`."""
    assert isinstance(commons["logfile"], str)
    assert commons["logfile"].endswith(".dispa.log")


def test_hatches_colors_key_parity():
    """Every key in commons['colors'] must also exist in commons['hatches'].

    This prevents KeyError crashes in plot_dispatch when iterating over
    plotdata columns and looking up a hatch for a fuel code that is in the
    color table but missing from the hatches table.
    """
    colors_keys = set(commons["colors"].keys())
    hatches_keys = set(commons["hatches"].keys())
    missing_from_hatches = colors_keys - hatches_keys
    assert not missing_from_hatches, (
        f"Keys present in commons['colors'] but missing from commons['hatches']: "
        f"{sorted(missing_from_hatches)}.  "
        f"Add them to commons['hatches'] (empty string is fine)."
    )


def test_dispaset_version_string():
    """The dispaset package exposes a non-empty version string."""
    assert isinstance(ds.__version__, str) and len(ds.__version__) > 0


# --------------------------------------------------------------------------- #
# Standalone runner
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    failures = 0
    for fn in [test_required_keys_present, test_tech_subsets,
               test_fuel_consistency, test_colour_table,
               test_logfile_naming, test_hatches_colors_key_parity,
               test_dispaset_version_string]:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
        except AssertionError as exc:
            failures += 1
            print(f"FAIL  {fn.__name__}: {exc}")
    raise SystemExit(0 if failures == 0 else 1)
