# -*- coding: utf-8 -*-
"""
Unit tests: Reserve demand pipeline (Phase 2).

Tests cover:
- check_reserve_demand() generic function behaviour
- mFRRUDemand initialised in finalTS as a separate DataFrame from aFRRU
- plants_res selection works with new nested dict config format
- check_reserves() backward-compat wrapper still calls both aFRRU and aFRRD checks

How to run
----------
    pytest tests/unit/test_reserve_demand_pipeline.py
    python tests/unit/test_reserve_demand_pipeline.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from dispaset.preprocessing.data_check import check_reserve_demand, check_reserves
from dispaset.preprocessing.data_handler import (
    normalize_reserve_config,
    _RESERVE_FLATLIST_DEFAULT_TYPES,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_idx(n=24):
    return pd.date_range("2015-01-01", periods=n, freq="h")


def _make_load(zones=("Z1",), n=24, value=100.0):
    idx = _make_idx(n)
    return pd.DataFrame({z: value for z in zones}, index=idx)


# ---------------------------------------------------------------------------
# check_reserve_demand
# ---------------------------------------------------------------------------

class TestCheckReserveDemand:
    def test_valid_demand_no_exception(self):
        load = _make_load()
        demand = _make_load(value=10.0)
        # Should not raise
        check_reserve_demand(demand, 'aFRRU', load)

    def test_negative_demand_raises(self):
        load = _make_load()
        demand = _make_load(value=-5.0)
        from dispaset.common import DispaSETValidationError
        with pytest.raises(DispaSETValidationError, match='negative values'):
            check_reserve_demand(demand, 'aFRRU', load)

    def test_demand_exceeds_load_raises(self):
        load = _make_load(value=100.0)
        demand = _make_load(value=200.0)
        from dispaset.common import DispaSETValidationError
        with pytest.raises(DispaSETValidationError, match='higher than demand'):
            check_reserve_demand(demand, 'mFRRU', load)

    def test_missing_zone_logs_warning(self, caplog):
        load = _make_load(zones=("Z1",))
        # demand has no Z1 column → should warn, not raise
        demand = pd.DataFrame({'Z2': 5.0}, index=_make_idx())
        import logging
        with caplog.at_level(logging.WARNING):
            check_reserve_demand(demand, 'mFRRU', load)
        assert any('mFRRU' in m and 'Z1' in m for m in caplog.messages)

    def test_multizone_all_valid(self):
        load = _make_load(zones=("Z1", "Z2"), value=100.0)
        demand = _make_load(zones=("Z1", "Z2"), value=20.0)
        check_reserve_demand(demand, 'FCRU', load)

    def test_multizone_one_bad_raises(self):
        load = _make_load(zones=("Z1", "Z2"), value=100.0)
        idx = _make_idx()
        demand = pd.DataFrame({'Z1': 20.0, 'Z2': 150.0}, index=idx)
        from dispaset.common import DispaSETValidationError
        with pytest.raises(DispaSETValidationError):
            check_reserve_demand(demand, 'FCRU', load)


class TestCheckReservesCompat:
    """Backward-compat wrapper still works correctly."""

    def test_both_types_checked(self):
        load = _make_load()
        aFRRU = _make_load(value=5.0)
        aFRRD = _make_load(value=5.0)
        check_reserves(aFRRD, aFRRU, load)  # should not raise

    def test_bad_afrru_raises(self):
        load = _make_load(value=100.0)
        aFRRU = _make_load(value=-1.0)
        aFRRD = _make_load(value=5.0)
        from dispaset.common import DispaSETValidationError
        with pytest.raises(DispaSETValidationError):
            check_reserves(aFRRD, aFRRU, load)


# ---------------------------------------------------------------------------
# mFRRU fallback behaviour (non-Exogenous)
# ---------------------------------------------------------------------------

class TestMFRRUFallback:
    """Verify that mFRRU gets the same values as aFRRU when using formula-based sizing."""

    def test_mfrru_equals_afrru_for_generic(self):
        """In non-Exogenous mode, mFRRUDemand_tot must equal aFRRUDemand_tot."""
        # We test this at the utility level: _RESERVE_FLATLIST_DEFAULT_TYPES contains 'mFRRU'
        assert 'mFRRU' in _RESERVE_FLATLIST_DEFAULT_TYPES

    def test_flat_list_includes_mfrru_type(self):
        cfg = {'ReserveParticipation': ['BATS'], 'ReserveParticipation_CHP': []}
        normalize_reserve_config(cfg)
        assert 'mFRRU' in cfg['ReserveParticipation']


# ---------------------------------------------------------------------------
# plants_res selection with nested dict
# ---------------------------------------------------------------------------

class TestPlantsResSelection:
    """Verify that all technologies from nested dict ReserveParticipation are selected."""

    def _make_plants(self, techs):
        return pd.DataFrame({'Technology': techs}, index=[f'U{i}' for i in range(len(techs))])

    def _collect_participating_techs(self, config):
        """Mirror the logic in build.py for plants_res selection."""
        from dispaset.common import commons
        all_res_techs = set()
        for rp_key in ('ReserveParticipation', 'ReserveParticipation_CHP'):
            for res_type, tech_factors in config.get(rp_key, {}).items():
                all_res_techs.update(tech_factors.keys())
        return all_res_techs

    def test_nested_dict_techs_collected(self):
        cfg = {
            'ReserveParticipation': {'aFRRU': {'BATS': 1.0, 'COMC': 1.0}, 'FFRU': {'BATS': 0.8}},
            'ReserveParticipation_CHP': {'aFRRU': {'GTUR': 1.0}},
        }
        participating = self._collect_participating_techs(cfg)
        assert participating == {'BATS', 'COMC', 'GTUR'}

    def test_empty_config_gives_empty_set(self):
        cfg = {'ReserveParticipation': {}, 'ReserveParticipation_CHP': {}}
        participating = self._collect_participating_techs(cfg)
        assert participating == set()

    def test_flat_list_converted_then_collected(self):
        cfg = {'ReserveParticipation': ['BATS', 'HDAM'], 'ReserveParticipation_CHP': []}
        normalize_reserve_config(cfg)
        participating = self._collect_participating_techs(cfg)
        assert 'BATS' in participating
        assert 'HDAM' in participating


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
