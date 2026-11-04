# -*- coding: utf-8 -*-
"""
Unit tests for the frequency-response plotting functions in
dispaset.postprocessing.frequency_stability: _resolve_contingency_key,
plot_frequency_response_curve, plot_frequency_response_for_hours.

No GAMS, no ODE solve -- results_frequency_response is a synthetic dict of
hand-built results_df frames, matching frequency_response()'s real column
layout, exactly what a real get_frequency_stability_reserves() call returns.
matplotlib is forced to the non-interactive Agg backend via tests/_helpers.py
(imported for its side effect, matching test_plots_ode_sfr.py's own pattern).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import _helpers  # noqa: F401,E402  -- forces matplotlib Agg backend

from dispaset.postprocessing.frequency_stability import (
    _resolve_contingency_key,
    _draw_reserve_service_traces,
    plot_frequency_response_curve,
    plot_frequency_response_for_hours,
)


def _make_results_df(peak=1.0):
    t = np.arange(-1.0, 10.0, 0.1)
    n = len(t)
    return pd.DataFrame({
        'Time [s]': t,
        'Frequency Deviation [Hz]': np.linspace(0, peak, n),
        'RoCoF [Hz/s]': np.linspace(0, -peak, n),
        'SIR [MW]': np.linspace(0, 50, n),
        'VIR [MW]': np.linspace(0, 20, n),
        'FFR [MW]': np.linspace(0, 10, n),
        'FCR [MW]': np.linspace(0, 30, n),
        'aFRR [MW]': np.linspace(0, 40, n),
        'mFRR [MW]': np.linspace(0, 5, n),
        'Contingency [MW]': np.full(n, 100.0),
        'deltap [MW]': np.linspace(0, -5, n),
    })


@pytest.fixture
def frr_and_index():
    idx = pd.date_range("2030-01-01", periods=3, freq="h", tz="UTC")
    results_frequency_response = {
        "Contingency1": _make_results_df(peak=0.3),
        "Contingency2": _make_results_df(peak=0.9),
        "Contingency3": _make_results_df(peak=0.5),
    }
    return results_frequency_response, idx


# ---------------------------------------------------------------------------
# _resolve_contingency_key
# ---------------------------------------------------------------------------

def test_resolve_by_integer_position(frr_and_index):
    frr, idx = frr_and_index
    label, key = _resolve_contingency_key(1, idx, frr)
    assert label == idx[1]
    assert key == "Contingency2"


def test_resolve_by_label(frr_and_index):
    frr, idx = frr_and_index
    label, key = _resolve_contingency_key(idx[2], idx, frr)
    assert label == idx[2]
    assert key == "Contingency3"


def test_resolve_by_literal_key(frr_and_index):
    frr, idx = frr_and_index
    label, key = _resolve_contingency_key("Contingency1", idx, frr)
    assert label == "Contingency1"
    assert key == "Contingency1"


def test_resolve_unknown_label_raises(frr_and_index):
    frr, idx = frr_and_index
    bogus = pd.Timestamp("2099-01-01", tz="UTC")
    with pytest.raises(KeyError):
        _resolve_contingency_key(bogus, idx, frr)


def test_resolve_position_out_of_range_raises(frr_and_index):
    frr, idx = frr_and_index
    with pytest.raises(IndexError):
        _resolve_contingency_key(10, idx, frr)  # index[10] doesn't exist (len(idx) == 3)


# ---------------------------------------------------------------------------
# plot_frequency_response_curve
# ---------------------------------------------------------------------------

def test_plot_frequency_response_curve_returns_fig_and_axes():
    df = _make_results_df()
    fig, (ax1, ax2) = plot_frequency_response_curve(df)
    assert fig is not None
    assert ax1.get_ylim() == (-1.3, 1.3)
    assert ax2.get_ylim() == (-300, 300)
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_plot_frequency_response_curve_respects_xlim():
    df = _make_results_df()
    fig, (ax1, ax2) = plot_frequency_response_curve(df, xlim=(0, 60))
    assert ax1.get_xlim() == (0, 60)
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_plot_frequency_response_curve_custom_ylims():
    df = _make_results_df()
    fig, (ax1, ax2) = plot_frequency_response_curve(df, freq_ylim=(-2, 2), power_ylim=(-500, 500))
    assert ax1.get_ylim() == (-2, 2)
    assert ax2.get_ylim() == (-500, 500)
    import matplotlib.pyplot as plt
    plt.close(fig)


# ---------------------------------------------------------------------------
# plot_frequency_response_for_hours
# ---------------------------------------------------------------------------

def test_plot_for_single_hour_produces_full_and_zoom(frr_and_index):
    frr, idx = frr_and_index
    out = plot_frequency_response_for_hours(frr, idx, hours=idx[0], show=False)
    assert list(out.keys()) == [idx[0]]
    assert out[idx[0]]['full'] is not None
    assert out[idx[0]]['zoom'] is not None
    import matplotlib.pyplot as plt
    plt.close('all')


def test_plot_for_multiple_hours_produces_one_entry_each(frr_and_index):
    frr, idx = frr_and_index
    hours = [idx[0], idx[2]]
    out = plot_frequency_response_for_hours(frr, idx, hours=hours, show=False)
    assert set(out.keys()) == set(hours)
    for h in hours:
        assert out[h]['full'][0] is not None
        assert out[h]['zoom'][0] is not None
    import matplotlib.pyplot as plt
    plt.close('all')


def test_plot_without_zoom_skips_zoom_figure(frr_and_index):
    frr, idx = frr_and_index
    out = plot_frequency_response_for_hours(frr, idx, hours=idx[0], include_zoom=False, show=False)
    assert out[idx[0]]['zoom'] is None
    import matplotlib.pyplot as plt
    plt.close('all')


def test_plot_accepts_mixed_hour_specifiers(frr_and_index):
    """hours can mix integer positions, labels, and literal keys in one call."""
    frr, idx = frr_and_index
    out = plot_frequency_response_for_hours(frr, idx, hours=[0, idx[1], "Contingency3"], show=False)
    assert idx[0] in out
    assert idx[1] in out
    assert "Contingency3" in out
    import matplotlib.pyplot as plt
    plt.close('all')
