# -*- coding: utf-8 -*-
"""
Unit tests: Reserve participation matrix population (Phase 3).

Tests that the 3-D ReserveParticipation array [res × unit × hour] is correctly
populated from the nested-dict config format, and that renewables receive
default participation for non-FFR reserve types.

How to run
----------
    pytest tests/unit/test_reserve_participation_matrix.py
    python tests/unit/test_reserve_participation_matrix.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from dispaset.common import commons


# ---------------------------------------------------------------------------
# Mirror of the participation matrix logic from build.py
# (without importing the full build pipeline)
# ---------------------------------------------------------------------------

def _build_participation_matrix(config, sets, plants_res):
    """
    Mirror of the ReserveParticipation build logic in build.py.
    Returns a 3D array: [len(res), len(units), len(hours)]
    """
    n_res = len(sets['res'])
    n_units = len(sets['au'])
    n_hours = len(sets['h'])

    values = np.zeros((n_res, n_units, n_hours))
    res_index = {r: i for i, r in enumerate(sets['res'])}
    au_index = {u: i for i, u in enumerate(sets['au'])}

    _participation_lookup = {}
    for rp_key in ('ReserveParticipation', 'ReserveParticipation_CHP'):
        for res_type, tech_factors in config.get(rp_key, {}).items():
            if res_type not in _participation_lookup:
                _participation_lookup[res_type] = {}
            _participation_lookup[res_type].update(tech_factors)

    _vre_non_ffr = [r for r in sets['res'] if r not in ('FFRU', 'FFRD')]

    for u in sets['au']:
        if u not in plants_res:
            continue
        i = au_index[u]
        tech = plants_res[u]

        for r in sets['res']:
            j = res_index[r]
            factor = 0.0
            if tech in commons['tech_renewables'] and r in _vre_non_ffr:
                factor = 1.0
            if r in _participation_lookup and tech in _participation_lookup[r]:
                factor = float(_participation_lookup[r][tech])
            values[j, i, :] = factor

    return values


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

RESERVE_TYPES = ['FFRU', 'FCRU', 'aFRRU', 'mFRRU', 'FFRD', 'FCRD', 'aFRRD']

def _default_sets(units=None, hours=3):
    units = units or ['BAT1', 'GAS1', 'SOL1']
    return {
        'res': RESERVE_TYPES,
        'res_U': ['FFRU', 'FCRU', 'aFRRU', 'mFRRU'],
        'res_D': ['FFRD', 'FCRD', 'aFRRD'],
        'au': units,
        'h': list(range(hours)),
    }


def _default_plants_res(mapping):
    """mapping: {unit_name: technology}"""
    return mapping


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestParticipationMatrixFromNestedDict:
    def test_battery_ffr_factor(self):
        cfg = {
            'ReserveParticipation': {'FFRU': {'BATS': 0.9}, 'FFRD': {'BATS': 0.9}},
            'ReserveParticipation_CHP': {},
        }
        sets = _default_sets(['BAT1'])
        plants_res = {'BAT1': 'BATS'}
        val = _build_participation_matrix(cfg, sets, plants_res)
        ffru_idx = RESERVE_TYPES.index('FFRU')
        ffrd_idx = RESERVE_TYPES.index('FFRD')
        assert val[ffru_idx, 0, 0] == pytest.approx(0.9)
        assert val[ffrd_idx, 0, 0] == pytest.approx(0.9)

    def test_unit_not_in_plants_res_stays_zero(self):
        cfg = {
            'ReserveParticipation': {'aFRRU': {'COMC': 1.0}},
            'ReserveParticipation_CHP': {},
        }
        sets = _default_sets(['GAS1', 'BAT1'])
        # BAT1 not in plants_res
        plants_res = {'GAS1': 'GTUR'}
        val = _build_participation_matrix(cfg, sets, plants_res)
        bat_idx = sets['au'].index('BAT1')
        assert np.all(val[:, bat_idx, :] == 0.0)

    def test_config_factor_overrides_renewable_default(self):
        cfg = {
            'ReserveParticipation': {'aFRRU': {'PHOT': 0.5}},
            'ReserveParticipation_CHP': {},
        }
        sets = _default_sets(['SOL1'])
        plants_res = {'SOL1': 'PHOT'}
        val = _build_participation_matrix(cfg, sets, plants_res)
        afrru_idx = RESERVE_TYPES.index('aFRRU')
        # Explicit config factor 0.5 overrides the renewable default of 1.0
        assert val[afrru_idx, 0, 0] == pytest.approx(0.5)

    def test_renewable_gets_default_for_non_ffr(self):
        cfg = {'ReserveParticipation': {}, 'ReserveParticipation_CHP': {}}
        sets = _default_sets(['SOL1'])
        plants_res = {'SOL1': 'PHOT'}
        val = _build_participation_matrix(cfg, sets, plants_res)
        ffru_idx = RESERVE_TYPES.index('FFRU')
        fcru_idx = RESERVE_TYPES.index('FCRU')
        afrru_idx = RESERVE_TYPES.index('aFRRU')
        # FFR types: renewables should NOT participate by default
        assert val[ffru_idx, 0, 0] == pytest.approx(0.0)
        # Non-FFR types: renewables DO participate at 1.0 by default
        assert val[fcru_idx, 0, 0] == pytest.approx(1.0)
        assert val[afrru_idx, 0, 0] == pytest.approx(1.0)

    def test_multiple_units_multiple_types(self):
        cfg = {
            'ReserveParticipation': {
                'FFRU': {'BATS': 1.0},
                'FCRU': {'BATS': 1.0, 'COMC': 1.0},
                'aFRRU': {'BATS': 0.8, 'COMC': 1.0, 'GTUR': 0.5},
                'mFRRU': {'GTUR': 1.0},
            },
            'ReserveParticipation_CHP': {'aFRRU': {'COMC': 0.9}},  # override
        }
        sets = _default_sets(['BAT1', 'GAS1', 'GT1'])
        plants_res = {'BAT1': 'BATS', 'GAS1': 'COMC', 'GT1': 'GTUR'}
        val = _build_participation_matrix(cfg, sets, plants_res)

        ffru_idx = RESERVE_TYPES.index('FFRU')
        afrru_idx = RESERVE_TYPES.index('aFRRU')
        mfrru_idx = RESERVE_TYPES.index('mFRRU')

        bat_idx, gas_idx, gt_idx = 0, 1, 2

        assert val[ffru_idx, bat_idx, 0] == pytest.approx(1.0)
        assert val[ffru_idx, gas_idx, 0] == pytest.approx(0.0)
        # CHP override: COMC aFRRU factor = 0.9 (from ReserveParticipation_CHP)
        assert val[afrru_idx, gas_idx, 0] == pytest.approx(0.9)
        assert val[afrru_idx, bat_idx, 0] == pytest.approx(0.8)
        assert val[mfrru_idx, gt_idx, 0] == pytest.approx(1.0)
        assert val[mfrru_idx, gas_idx, 0] == pytest.approx(0.0)

    def test_all_zeros_for_empty_config(self):
        cfg = {'ReserveParticipation': {}, 'ReserveParticipation_CHP': {}}
        sets = _default_sets(['GAS1'])
        plants_res = {'GAS1': 'COMC'}
        val = _build_participation_matrix(cfg, sets, plants_res)
        assert np.all(val == 0.0)

    def test_hours_dimension_uniform(self):
        cfg = {'ReserveParticipation': {'aFRRU': {'BATS': 0.7}}, 'ReserveParticipation_CHP': {}}
        sets = _default_sets(['BAT1'], hours=10)
        plants_res = {'BAT1': 'BATS'}
        val = _build_participation_matrix(cfg, sets, plants_res)
        afrru_idx = RESERVE_TYPES.index('aFRRU')
        assert val.shape[2] == 10
        assert np.all(val[afrru_idx, 0, :] == pytest.approx(0.7))

    def test_factor_zero_disables_participation(self):
        cfg = {'ReserveParticipation': {'aFRRU': {'BATS': 0.0}}, 'ReserveParticipation_CHP': {}}
        sets = _default_sets(['BAT1'])
        plants_res = {'BAT1': 'BATS'}
        val = _build_participation_matrix(cfg, sets, plants_res)
        afrru_idx = RESERVE_TYPES.index('aFRRU')
        assert val[afrru_idx, 0, 0] == pytest.approx(0.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
