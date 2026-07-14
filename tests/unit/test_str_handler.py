# -*- coding: utf-8 -*-
"""
Unit test: dispaset.misc.str_handler
====================================

What this test does
-------------------
The `str_handler` module is responsible for two things that frequently break
the GAMS handshake:

1. `shrink_to_64` truncates strings to fit the GDX 64-character limit while
   keeping uniqueness. We test that:
   - Short strings are returned unchanged.
   - Long strings are shortened to <= 63 characters.
   - Different inputs produce different outputs (uniqueness).
2. `force_str` ensures bytes / numpy scalar objects are converted to plain
   ``str`` so the GAMS API does not crash.
3. `clean_strings` strips spaces and disallowed characters that would mess
   with column-header / unit-name matching.

How to run
----------

1. ``pytest tests/unit/test_str_handler.py``
2. ``python tests/unit/test_str_handler.py``

This test does NOT require GAMS.
"""
from __future__ import annotations

import sys
from pathlib import Path

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from dispaset.misc.str_handler import shrink_to_64, force_str, clean_strings  # noqa: E402


def test_force_str_from_bytes():
    """Bytes objects must be returned as decoded ``str``."""
    assert force_str(b"hello") == "hello"


def test_force_str_idempotent():
    """A ``str`` should be returned unchanged."""
    assert force_str("Z1 - GAS") == "Z1 - GAS"


def test_shrink_to_64_short_strings():
    """A list of strings <= 63 chars must be returned untouched."""
    inp = ["abc", "Z1 - GAS", "Wind_Z1"]
    out = shrink_to_64(inp)
    assert out == inp


def test_shrink_to_64_long_strings():
    """Long strings are shortened to fit the GDX 64-char rule."""
    long1 = "VeryVeryLongUnitName_" + "x" * 60
    long2 = "VeryVeryLongUnitName_" + "y" * 60
    out = shrink_to_64([long1, long2])
    assert all(len(s) <= 63 for s in out), f"Output too long: {out}"
    # Uniqueness must be preserved
    assert out[0] != out[1], "Different inputs must remain different"


def test_clean_strings_removes_invalid_chars():
    """`clean_strings` should strip leading/trailing whitespace."""
    raw = [" Z1 - GAS ", " Z2 - WIN"]
    cleaned = clean_strings(raw)
    assert all(not s.startswith(" ") and not s.endswith(" ") for s in cleaned)


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    failures = 0
    for fn in [test_force_str_from_bytes, test_force_str_idempotent,
               test_shrink_to_64_short_strings, test_shrink_to_64_long_strings,
               test_clean_strings_removes_invalid_chars]:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
        except AssertionError as exc:
            failures += 1
            print(f"FAIL  {fn.__name__}: {exc}")
    raise SystemExit(0 if failures == 0 else 1)
