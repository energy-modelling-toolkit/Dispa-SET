# -*- coding: utf-8 -*-
"""
Unit tests: Reserve configuration normalization and validation (Phase 1).

Tests cover:
- Backward-compat flat list → nested dict conversion
- New nested dict format: round-trip and value access
- Invalid reserve type names produce a warning (not an error)
- Out-of-range factors produce a warning
- Legacy demand aliases still work
- mFRRUDemand present in loaded configs
- Empty ReserveParticipation is handled gracefully

How to run
----------
    pytest tests/unit/test_reserve_config.py
    python tests/unit/test_reserve_config.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from dispaset.preprocessing.data_handler import (
    normalize_reserve_config,
    VALID_RESERVE_TYPES,
    _flatlist_to_reserve_dict,
    _RESERVE_FLATLIST_DEFAULT_TYPES,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _base_config():
    """Minimal config dict sufficient for normalize_reserve_config()."""
    return {
        'ReserveParticipation': [],
        'ReserveParticipation_CHP': [],
    }


# ---------------------------------------------------------------------------
# Flat-list conversion
# ---------------------------------------------------------------------------

class TestFlatListConversion:
    def test_flat_list_converts_to_nested_dict(self):
        cfg = _base_config()
        cfg['ReserveParticipation'] = ['BATS', 'COMC', 'GTUR']
        normalize_reserve_config(cfg)
        rp = cfg['ReserveParticipation']
        assert isinstance(rp, dict), "ReserveParticipation should become a dict"

    def test_flat_list_default_types_populated(self):
        cfg = _base_config()
        cfg['ReserveParticipation'] = ['BATS', 'COMC']
        normalize_reserve_config(cfg)
        rp = cfg['ReserveParticipation']
        for res_type in _RESERVE_FLATLIST_DEFAULT_TYPES:
            assert res_type in rp, f"Default type {res_type!r} should be present after flat-list conversion"

    def test_flat_list_factor_is_one(self):
        cfg = _base_config()
        cfg['ReserveParticipation'] = ['BATS']
        normalize_reserve_config(cfg)
        rp = cfg['ReserveParticipation']
        for res_type in _RESERVE_FLATLIST_DEFAULT_TYPES:
            assert rp[res_type]['BATS'] == 1.0

    def test_empty_flat_list_gives_empty_dict(self):
        cfg = _base_config()
        cfg['ReserveParticipation'] = []
        normalize_reserve_config(cfg)
        assert cfg['ReserveParticipation'] == {}

    def test_chp_flat_list_also_converted(self):
        cfg = _base_config()
        cfg['ReserveParticipation_CHP'] = ['COMC']
        normalize_reserve_config(cfg)
        rp_chp = cfg['ReserveParticipation_CHP']
        assert isinstance(rp_chp, dict)
        for res_type in _RESERVE_FLATLIST_DEFAULT_TYPES:
            assert res_type in rp_chp

    def test_missing_key_defaults_to_empty_dict(self):
        cfg = {}
        normalize_reserve_config(cfg)
        assert cfg.get('ReserveParticipation', {}) == {}
        assert cfg.get('ReserveParticipation_CHP', {}) == {}


# ---------------------------------------------------------------------------
# Nested dict round-trip
# ---------------------------------------------------------------------------

class TestNestedDictFormat:
    def test_nested_dict_preserved_unchanged(self):
        cfg = _base_config()
        cfg['ReserveParticipation'] = {
            'FFRU': {'BATS': 1.0},
            'aFRRU': {'BATS': 0.8, 'COMC': 1.0},
            'aFRRD': {'COMC': 0.5},
            'mFRRU': {'GTUR': 1.0, 'BATS': 1.0},
        }
        normalize_reserve_config(cfg)
        rp = cfg['ReserveParticipation']
        assert rp['FFRU']['BATS'] == 1.0
        assert rp['aFRRU']['BATS'] == 0.8
        assert rp['aFRRD']['COMC'] == 0.5
        assert rp['mFRRU']['GTUR'] == 1.0

    def test_all_valid_reserve_types_accepted(self):
        cfg = _base_config()
        cfg['ReserveParticipation'] = {res: {'BATS': 1.0} for res in VALID_RESERVE_TYPES}
        # Should not raise
        normalize_reserve_config(cfg)

    def test_unknown_reserve_type_logged_but_no_exception(self, caplog):
        cfg = _base_config()
        cfg['ReserveParticipation'] = {'UNKNOWN_TYPE': {'BATS': 1.0}}
        import logging
        with caplog.at_level(logging.WARNING):
            normalize_reserve_config(cfg)
        assert any('unknown reserve types' in m for m in caplog.messages)

    def test_out_of_range_factor_logged_but_no_exception(self, caplog):
        cfg = _base_config()
        cfg['ReserveParticipation'] = {'aFRRU': {'COMC': 2.5}}
        import logging
        with caplog.at_level(logging.WARNING):
            normalize_reserve_config(cfg)
        assert any('outside [0, 1]' in m for m in caplog.messages)

    def test_invalid_inner_type_raises(self):
        cfg = _base_config()
        cfg['ReserveParticipation'] = {'aFRRU': ['BATS', 'COMC']}  # list instead of dict
        with pytest.raises(ValueError, match='must be a dict'):
            normalize_reserve_config(cfg)


# ---------------------------------------------------------------------------
# Legacy demand aliases
# ---------------------------------------------------------------------------

class TestLegacyDemandAliases:
    def test_reserve2u_maps_to_afrrudemand(self):
        cfg = _base_config()
        cfg['Reserve2U'] = 'path/to/reserve2u.csv'
        cfg['aFRRUDemand'] = ''
        normalize_reserve_config(cfg)
        assert cfg['aFRRUDemand'] == 'path/to/reserve2u.csv'

    def test_reserve2d_maps_to_afrrdemand(self):
        cfg = _base_config()
        cfg['Reserve2D'] = 'path/to/reserve2d.csv'
        cfg['aFRRDDemand'] = ''
        normalize_reserve_config(cfg)
        assert cfg['aFRRDDemand'] == 'path/to/reserve2d.csv'

    def test_primaryreservelimit_maps_to_fcrdemand(self):
        cfg = _base_config()
        cfg['PrimaryReserveLimit'] = 'path/to/fcr.csv'
        cfg['FCRDemand'] = ''
        normalize_reserve_config(cfg)
        assert cfg['FCRDemand'] == 'path/to/fcr.csv'

    def test_existing_new_key_takes_precedence(self, caplog):
        cfg = _base_config()
        cfg['aFRRUDemand'] = 'new_path.csv'
        cfg['Reserve2U'] = 'old_path.csv'
        import logging
        with caplog.at_level(logging.WARNING):
            normalize_reserve_config(cfg)
        assert cfg['aFRRUDemand'] == 'new_path.csv'
        assert any('"aFRRUDemand"' in m and '"Reserve2U"' in m for m in caplog.messages)


# ---------------------------------------------------------------------------
# Full load_config integration
# ---------------------------------------------------------------------------

class TestLoadConfigIntegration:
    CONFIG_DIR = Path(__file__).resolve().parents[1] / "configs"

    def test_mfRRUDemand_in_loaded_config(self):
        """mFRRUDemand must be present (as empty string) in any loaded config."""
        import dispaset as ds
        cfg = ds.load_config(str(self.CONFIG_DIR / "tiny.yml"))
        assert 'mFRRUDemand' in cfg, "mFRRUDemand should be present after loading config"

    def test_tiny_flat_list_converted_on_load(self):
        """tiny.yml uses a flat list; after loading it must be a nested dict."""
        import dispaset as ds
        cfg = ds.load_config(str(self.CONFIG_DIR / "tiny.yml"))
        rp = cfg['ReserveParticipation']
        assert isinstance(rp, dict), "ReserveParticipation should be a nested dict after loading"

    def test_tiny_participation_contains_default_types(self):
        """After flat-list conversion, the default types must be present."""
        import dispaset as ds
        cfg = ds.load_config(str(self.CONFIG_DIR / "tiny.yml"))
        for res_type in _RESERVE_FLATLIST_DEFAULT_TYPES:
            assert res_type in cfg['ReserveParticipation'], \
                f"Default reserve type {res_type!r} missing after flat-list conversion"

    def test_valid_reserve_types_constant(self):
        assert set(VALID_RESERVE_TYPES) == {'FFRU', 'FCRU', 'aFRRU', 'mFRRU', 'FFRD', 'FCRD', 'aFRRD'}


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
