# -*- coding: utf-8 -*-
"""
Frequency Response Assessment Model (ODE-SFR): sequential, ordered-backtracking
sizing of synchronous inertia, virtual inertia, and frequency reserves (FFR,
FCR, aFRR, mFRR) against RoCoF, nadir, and steady-state frequency-deviation
limits, per contingency.

Extracted out of postprocessing.py into its own module, mirroring how
reserves.py, boundary_sector.py, geoplot.py, and streamlit_functions.py are
each split out of the larger preprocessing/postprocessing files. Pure move,
no behavior change.

See ODESFR_ISSUES.md and ODESFR_DECISIONS.md at the repository root for the full
change history, known issues, and open decisions for this subsystem.
"""
import logging
import time

import numpy as np
import pandas as pd
from scipy.integrate import odeint

from ..common import commons


# %% trapezoid weights
def trapezoid_weight(tt, prep, ramp, delivery=None, deact=None):
    """
    this function builds a trapezoidal (0→1→0) or triangular (0→1) profile
    """
    if delivery is None or deact is None:
        # Triangular case: ramp only (e.g. mFRR)
        return np.piecewise(tt,
            [tt < prep,
             (prep <= tt) & (tt < ramp),
             tt >= ramp],
            [0,
             lambda t: (t - prep) / (ramp - prep),
             1]
        )
    else:
        # Trapezoidal case: ramp → delivery → deact
        return np.piecewise(tt,
            [tt < prep,
             (prep <= tt) & (tt < ramp),
             (ramp <= tt) & (tt < delivery),
             (delivery <= tt) & (tt < deact),
             tt >= deact],
            [0,
             lambda t: (t - prep) / (ramp - prep),
             1,
             lambda t: 1 - (t - delivery) / (deact - delivery),
             0]
        )
# %% compute weights
def compute_weights(tt, activation_times):
    """
    This function calculates and returns w_vir, w_ffr, w_fcr, w_afrr, w_mfrr for vector tt.
    """
    w_vir = trapezoid_weight(tt, **activation_times["vir"])
    w_ffr = trapezoid_weight(tt, **activation_times["ffr"])
    w_fcr = trapezoid_weight(tt, **activation_times["fcr"])
    w_afrr = trapezoid_weight(tt, **activation_times["afrr"])
    w_mfrr = trapezoid_weight(tt, **activation_times["mfrr"])
    return w_vir, w_ffr, w_fcr, w_afrr, w_mfrr


# %% frequency_response
# Saturation of the fleet FCR response (set by get_frequency_stability_reserves from the units):
# (K_g [MW/Hz], d_g [Hz], K_max [MW/Hz]) per group of units with the same saturation deviation d = S/K.
_FCR_SAT = None


def _fcr_shape(f):
    """Effective deviation of the procured FCR: FCR(df) = FCR_cap * _fcr_shape(df). Each unit gives
    K_u*df up to its capability S_u; procuring a share FCR_cap/K_max of the fleet pro rata gives
    FCR_cap/K_max * sum_u K_u*min(df, d_u). Identity when no saturation is defined."""
    if _FCR_SAT is None:
        return f
    Kg, dg, Kmax = _FCR_SAT
    f = np.asarray(f, dtype=float)
    return (np.sum(Kg[:, None] * np.clip(f.reshape(1, -1), -dg[:, None], dg[:, None]), axis=0) / Kmax).reshape(f.shape)


def _steady_state_ok(results_df, limit_freq_steady_state):
    """Quasi-steady deviation after full primary activation and recovery to the normal band (NOp 00)."""
    c = commons['FrequencyResponseConstants']
    t = results_df['Time [s]']; df = results_df['Frequency Deviation [Hz]'].abs()
    ss = df[t >= c.get('SteadyStateStart', 350.0)]
    rec = df[t >= c.get('RecoveryTime', np.inf)]
    ok_ss = (not ss.empty) and ss.max() <= limit_freq_steady_state
    ok_rec = rec.empty or rec.max() <= c.get('RecoveryDeviationMax', np.inf)
    return ok_ss and ok_rec


def _response_metrics(results_df, limit_freq, limit_rocof, limit_freq_steady_state):
    """Maximum deviation, maximum RoCoF and quasi-steady deviation of a frequency response, and the
    criteria it violates ('' if none), so that each hour reports how close it is to the limits."""
    c = commons['FrequencyResponseConstants']
    t = results_df['Time [s]']
    dev = results_df['Frequency Deviation [Hz]'].abs()
    max_fd = float(dev.max())
    max_r = float(results_df['RoCoF [Hz/s]'].abs().max()) if 'RoCoF [Hz/s]' in results_df else np.nan
    ss = dev[t >= c.get('SteadyStateStart', 350.0)]
    rec = dev[t >= c.get('RecoveryTime', np.inf)]
    ss_dev = float(ss.max()) if not ss.empty else np.nan
    failed = []
    if max_fd > limit_freq:
        failed.append('nadir')
    if np.isfinite(max_r) and max_r > limit_rocof:
        failed.append('RoCoF')
    if ss.empty or ss_dev > limit_freq_steady_state:
        failed.append('steady-state')
    if not rec.empty and rec.max() > c.get('RecoveryDeviationMax', np.inf):
        failed.append('recovery')
    return max_fd, max_r, ss_dev, '+'.join(failed)


def frequency_response(sim_time, activation_times, Hs_cap, Hv_cap, FFR_cap, FCR_cap, aFRR_cap, mFRR_cap,
                   contingency, D_local, verbose=False):
    """
    This function solves the power swing differential equation, for each combination
    of inertia and frequency reserves desired.
    param sim_time:             Time to evaluate the differential equation in seconds [s]
    param activation_times:     Activation times for each reserve, measured from the contingency onset (t=0)
    param Hs_cap:                System Inertia value
    param FFR_cap:              Fast Frequency Reserve capacity
    param FCR_cap:              Frequency Containment Reserve capacity
    param aFRR_cap:             Automatic Frequency Restoration Reserve capacity
    param mFRR_cap:             Manual Frequency Restoration Reserve capacity
    param contingency:          Contingency size assumed by the n-1 hypothesis
    param D_local:              Damping value assumed to be 1.5% of the load

    returns dataframe results_df:
                                - Time [s]
                                - Frequency Deviation [Hz]
                                - RoCoF [Hz/s]
                                - Inertia [GWs]
                                - FFR [MW]
                                - FCR [MW]
                                - aFRR [MW]
                                - mFRR [MW]
                                - Damping [MW/Hz]'
                                - Contingency [MW]'
                                - deltap [MW]
    """
    # If the value of Inertia is zero H=0
    if Hs_cap == 0:
        # ODE-SFR-CHANGE-001 (ISSUE-ODESFR-003): this branch must return the same
        # 9-tuple shape as every other path below, since every call site
        # unpacks 9 values unconditionally (including the "no valid H found"
        # fallback in get_frequency_stability_reserves, which sets Hs_cap=0
        # and would otherwise crash here with ValueError on unpack).
        if verbose:
            print("Hs_cap=0: unstable system, returning infinite deviations.")
        results_df = pd.DataFrame({c: [] for c in [
            'Time [s]', 'Frequency Deviation [Hz]', 'RoCoF [Hz/s]',
            'Synchronous Inertia Constant [s]', 'Virtual Inertia Constant [s]',
            'SIR [MW]', 'VIR [MW]', 'FFR [MW]', 'FCR [MW]', 'aFRR [MW]', 'mFRR [MW]',
            'Damping [MW/Hz]', 'Contingency [MW]', 'deltap [MW]']})
        return (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, results_df)

    # Simulation horizon and time step; t=0 is the contingency onset itself.
    t = np.arange(0, sim_time, 0.1)

    # Precompute weights vector for plotting and for fixed deployment logic.
    # The "vir" weight (first return value) is intentionally discarded: it is
    # never used below either -- H_total combines Hs_cap+Hv_cap unweighted
    # (see ISSUE-ODESFR-004 / DECISION-ODESFR-006 for the separate, unresolved
    # question of whether Hv should be gated by its own activation profile;
    # this fix does not touch that).
    _, W_ffr_vec, W_fcr_vec, W_afrr_vec, W_mfrr_vec = compute_weights(t, activation_times)

    f_0 = 50.0  # Nominal frequency [Hz]
    S_base = 1000  # System base [MW]

    def state(y, tt):
        f = y[0]

        # ODE-SFR-CHANGE-003 (ISSUE-ODESFR-002): interpolate the already-precomputed
        # weight vectors instead of recomputing compute_weights() from scratch
        # at every internal ODE step (measured ~1,470 recomputations per call
        # before this fix; interpolation gives numerically equivalent output
        # -- max diff ~1e-8 to 1e-12 Hz across benchmarked cases -- at a
        # measured 13-16x speedup, see ODESFR_ISSUES.md ISSUE-ODESFR-002).
        w_ffr = np.interp(tt, t, W_ffr_vec)
        w_fcr = np.interp(tt, t, W_fcr_vec)
        w_afrr = np.interp(tt, t, W_afrr_vec)
        w_mfrr = np.interp(tt, t, W_mfrr_vec)

        ffr = FFR_cap * w_ffr * f
        fcr = FCR_cap * w_fcr * _fcr_shape(f)
        afrr = aFRR_cap * w_afrr
        mfrr = mFRR_cap * w_mfrr

        H_total = Hs_cap + Hv_cap

        deltap = contingency - (ffr + fcr + afrr + mfrr) - D_local * f
        dfdt = (deltap / (2 * H_total * S_base)) * f_0

        return [dfdt]

    sol = odeint(state, [0.0], t)
    f = sol[:, 0]

    contingency_vec = np.full(t.shape, float(contingency))

    max_freq_dev = np.max(np.abs(f))
    rocof = np.gradient(f, t)
    max_rocof = np.max(np.abs(rocof))

    sir = 2 * Hs_cap * S_base * rocof/ f_0
    vir = 2 * Hv_cap * S_base * rocof/ f_0
    ffr = FFR_cap * f * W_ffr_vec
    fcr = FCR_cap * _fcr_shape(f) * W_fcr_vec
    afrr = aFRR_cap * W_afrr_vec
    mfrr= mFRR_cap * W_mfrr_vec


    max_sir = np.max(np.abs(2 * Hs_cap * S_base * rocof/ f_0))
    max_vir = np.max(np.abs(2 * Hv_cap * S_base * rocof/ f_0))
    max_ffr = np.max(ffr)
    max_fcr = np.max(fcr)
    max_afrr = np.max(afrr)
    max_mfrr = np.max(mfrr)

    deltap = contingency_vec - (ffr + fcr + afrr + mfrr) - D_local * f

    results_df = pd.DataFrame({
        'Time [s]': t,
        'Frequency Deviation [Hz]': -f,
        'RoCoF [Hz/s]': -rocof,
        'Synchronous Inertia Constant [s]': Hs_cap,
        'Virtual Inertia Constant [s]': Hv_cap,
        'SIR [MW]': sir,
        'VIR [MW]': vir,
        'FFR [MW]': ffr,
        'FCR [MW]': fcr,
        'aFRR [MW]': afrr,
        'mFRR [MW]': mfrr,
        'Damping [MW/Hz]': D_local,
        'Contingency [MW]': -contingency_vec,
        'deltap [MW]': -deltap
    })

    # Cosmetic 1s pre-roll of nominal (pre-contingency) state, for plotting only.
    preroll_t = np.arange(-1.0, 0, 0.1)
    preroll_df = pd.DataFrame({
        'Time [s]': preroll_t,
        'Frequency Deviation [Hz]': 0.0,
        'RoCoF [Hz/s]': 0.0,
        'Synchronous Inertia Constant [s]': Hs_cap,
        'Virtual Inertia Constant [s]': Hv_cap,
        'SIR [MW]': 0.0,
        'VIR [MW]': 0.0,
        'FFR [MW]': 0.0,
        'FCR [MW]': 0.0,
        'aFRR [MW]': 0.0,
        'mFRR [MW]': 0.0,
        'Damping [MW/Hz]': D_local,
        'Contingency [MW]': 0.0,
        'deltap [MW]': 0.0,
    })
    results_df = pd.concat([preroll_df, results_df], ignore_index=True)

    if verbose:
        print(f"Sim finished: H={Hs_cap:.2f},H={Hv_cap:.2f},SIR={max_sir:.2f},VIR={max_vir:.2f},FFR={max_ffr:.2f},FCR={max_fcr:.2f},aFRR={max_afrr:.2f},mFRR={max_mfrr:.2f}")
        print(f" -> max_freq_dev={max_freq_dev:.4f} Hz, max_rocof={max_rocof:.4f} Hz/s")

    return max_freq_dev, max_rocof, max_sir, max_vir, max_ffr, max_fcr, max_afrr, max_mfrr, results_df


# %% ordered-backtracking bisection helper
def _bisect_collect_feasible(low, high, tol, is_feasible):
    """
    Binary-search [low, high]; return every tested value that satisfied
    is_feasible, in the order they were found (which narrows toward the
    smallest feasible value, matching the "test the upper bound first" style
    already used by ODE-SFR's stage bisections).

    Used by get_frequency_stability_reserves (ODE-SFR-CHANGE-002, ISSUE-ODESFR-001 /
    DECISION-ODESFR-002) to generate the ordered candidate list for each of the
    four bisected stages (Hs, Hv, FFR, FCR), so the calling loop can try
    candidates smallest-first and backtrack to a larger one only if every
    downstream possibility for the smaller one is exhausted.

    Bracket handling for small or inverted ranges (ISSUE-ODESFR-floor): the
    lower bound is clipped to the upper bound (``low = min(low, high)``), so a
    fleet bound below the nominal lower bound (e.g. a small BESS fleet with
    H^V_max = 0.715 GWs < H_lo = 1 GWs) still has its upper bound tested
    instead of producing an inverted, never-searched bracket. When the bracket
    is narrower than ``tol`` the tolerance is reduced to a quarter of the
    bracket, so the bound and three midpoints are still tested; a zero-width
    bracket tests the bound once. For brackets at least ``tol`` wide (every
    case with integer bounds >= low + tol) the behaviour is unchanged.
    A non-positive upper bound returns no candidate (nothing to search).
    """
    feasible = []
    if not high > 0:
        return feasible
    low = min(low, high)
    if high - low < tol:
        tol = (high - low) / 4.0
    if tol <= 0:
        return [high] if is_feasible(high) else feasible
    first_iter = True
    while low + tol <= high:
        mid = high if first_iter else (low + high) / 2
        if is_feasible(mid):
            feasible.append(mid)
            high = mid
        elif not first_iter:
            low = mid
        first_iter = False
    return feasible


_RESERVE_TIMING_TO_ODESFR_KEY = {'VIRU': 'vir', 'FFRU': 'ffr', 'FCRU': 'fcr', 'aFRRU': 'afrr', 'mFRRU': 'mfrr'}


def _as_odesfr_activation_times(reserve_timing):
    """Translate commons['ReserveTiming'] keys (VIRU/FFRU/...) into ODE-SFR's lowercase-key shape
    (vir/ffr/...), backfilling any missing category from commons['ReserveTiming']."""
    defaults = commons['ReserveTiming']
    return {
        odesfr_key: reserve_timing.get(res_key, defaults[res_key])
        for res_key, odesfr_key in _RESERVE_TIMING_TO_ODESFR_KEY.items()
    }


def resolve_activation_times(activation_times, config=None):
    """Resolve activation_times to use and its source: 'explicit' (caller-provided, warns if
    inconsistent with the case's own timing), 'config' (config['ReserveTiming']), or 'default'."""
    case_reserve_timing = (config or {}).get('ReserveTiming', commons['ReserveTiming'])

    if activation_times is not None:
        if activation_times != _as_odesfr_activation_times(case_reserve_timing):
            logging.warning(
                "get_frequency_stability_reserves(): activation_times differs from this "
                "case's own reserve-timing assumptions. Treat these results as a what-if "
                "sensitivity, not a validation of the solved case's reserve procurement. To "
                "make UCED consistent, set config['ReserveTiming'] and rebuild + resolve."
            )
        return activation_times, 'explicit'

    if 'ReserveTiming' in (config or {}):
        return _as_odesfr_activation_times(case_reserve_timing), 'config'

    return _as_odesfr_activation_times(commons['ReserveTiming']), 'default'


# %% frequency stability reserves
# %% hourly availability of the fleet (search bounds) and gap-refined ordered backtracking
# Saturation table of the FCR response per hour (set by get_frequency_stability_reserves).
_FCR_SAT_BY_INDEX = {}


def _hourly_param(inputs, name, units_index, index, default):
    """Hourly unit parameter (hours x units) from inputs['param_df'][name] (DataFrame indexed by
    date, as built by ds_to_df) or, failing that, inputs['parameters'][name] (GAMS-style array,
    sets ['au', 'h']). Missing units/hours take ``default``."""
    df = None
    pdf = inputs.get('param_df') if isinstance(inputs, dict) else None
    if isinstance(pdf, dict) and isinstance(pdf.get(name), pd.DataFrame):
        df = pdf[name]
    elif isinstance(inputs, dict) and name in inputs.get('parameters', {}):
        par = inputs['parameters'][name]
        cols = list(inputs['sets'][par['sets'][0]])
        df = pd.DataFrame(np.asarray(par['val'], dtype=float).T, columns=cols)
    if df is None:
        return None
    df = df.loc[:, ~df.columns.duplicated()]
    if len(index.intersection(df.index)) == 0:
        # no common time index (e.g. integer-indexed array): align by position from the first hour
        if len(df) < len(index):
            logging.warning("unit_hourly_availability(): %s has %d rows < %d hours; using default %s",
                            name, len(df), len(index), default)
            return None
        df = df.iloc[:len(index)].set_axis(index, axis=0)
    out = df.reindex(index=index, columns=units_index)
    n_missing = int(out.isna().to_numpy().sum())
    if n_missing:
        logging.info("unit_hourly_availability(): %d missing %s values set to %s", n_missing, name, default)
    return out.astype(float).fillna(default)


def unit_hourly_availability(inputs, units_index, index):
    """Hourly available share of each unit, LoadMaximum(u,h) = AF(u,h) * (1 - OF(u,h)), the same
    factor the UCED applies to PowerCapacity (UCM.gms). Clipped to [0, 1]; 1 (installed fleet)
    where the inputs carry no availability data."""
    units_index = pd.Index(units_index)
    af = _hourly_param(inputs, 'AvailabilityFactor', units_index, index, 1.0)
    of = _hourly_param(inputs, 'OutageFactor', units_index, index, 0.0)
    if af is None and of is None:
        logging.warning("unit_hourly_availability(): no AvailabilityFactor/OutageFactor in inputs -- "
                        "search bounds fall back to the installed fleet")
    af = af if af is not None else pd.DataFrame(1.0, index=index, columns=units_index)
    of = of if of is not None else pd.DataFrame(0.0, index=index, columns=units_index)
    return (af * (1.0 - of)).clip(lower=0.0, upper=1.0)


def _ordered_backtracking(n_stages, candidates, accept, refine_tols=None):
    """Ordered backtracking over ``n_stages`` bisected stages (ODE-SFR-CHANGE-002) with gap
    refinement (ISSUE-ODESFR-ladder).

    ``candidates(k, q)`` returns the ascending candidate list of stage k given the values ``q``
    fixed at stages 0..k-1; ``accept(q)`` is the final (Stage-5) check of a complete combination
    and returns a payload, or None if the combination fails. Candidates are tried smallest-first
    and a stage advances only after the whole downstream chain failed for the current candidate.

    Gap refinement: the candidates of a stage are the values visited by its bisection, so
    consecutive candidates can differ by a factor of two. When candidate c_j passes downstream
    after c_{j-1} failed, the gap (c_{j-1}, c_j] is bisected (feasibility = "some downstream chain
    passes") down to ``refine_tols[k]``; every returned combination has passed ``accept``.

    :return: (values, payload) or None when no combination passes.
    """
    refine_tols = list(refine_tols) if refine_tols is not None else [None] * n_stages

    def solve(k, q):
        if k == n_stages:
            payload = accept(q)
            return None if payload is None else (list(q), payload)
        prev = None
        for c in candidates(k, q):
            sol = solve(k + 1, q + [c])
            if sol is None:
                prev = c
                continue
            tol = refine_tols[k]
            if tol and prev is not None:
                lo, hi = prev, c
                while hi - lo > tol:
                    mid = 0.5 * (lo + hi)
                    s = solve(k + 1, q + [mid])
                    if s is None:
                        lo = mid
                    else:
                        hi, sol = mid, s
            return sol
        return None

    return solve(0, [])


def get_frequency_stability_reserves(path, inputs, results, activation_times=None,
                                     limit_freq=None, limit_rocof=None, limit_freq_steady_state=None,
                                     use_vir=None, use_ffr=None, use_fcr=None, use_afrr=None, use_mfrr=None,
                                     fleet_bounds='available', refine=True,
                                     feasibility_only=False, stop_on_insecure=False):
    """
    Performs a sequential binary search over different combinations of.

    This function iteratively explores combinations of inertia and reserves (H, FFR, FCR, aFRR, and mFRR)
    throughtout a binary search and considering the specific activation time windows of each service.
    For each candidate solution, the function calls `frequency_response` to simulate the system's dynamic
    frequency behavior after a contingency, and verifies whether the minimum stability
    requirements (limit_freq, limit_rocof, limit_freq_stead_state) are satisfied.


    param inputs:              DispaSET inputs, needed to compute the damping value and the maximum system inertia available
    param results:             DispaSET results, needed to compute the size of contingency for the preliminar simulation
    param activation_times:    Reserve activation times
    limit_freq:                Maximum frequency deviation allowed by each TSO (default value 0.8 Hz)
    limit_rocof:               Maximum rate of change of frequency allowed by each TSO (default value 0.5 Hz/s)
    limit_freq_steady_state:   Maximum frequency deviation allowed by the TSO in the steady state
                               (default value 0.2 Hz, assumed to be 200 seconds after the contingency)
    fleet_bounds:              'available' (default): hourly search bounds from the available fleet,
                               nameplate x AF x (1 - OF) per unit and hour, so an hour is secure only
                               if its requirement can be delivered by available units;
                               'installed': constant nameplate bounds (previous behaviour)
    refine:                    True (default): bisect the gap between a candidate that failed
                               downstream and the next larger one that passed (see
                               _ordered_backtracking); False: previous behaviour
    feasibility_only:          True: skip the search and only check whether the hour is secure with
                               every service at its hourly bound (the criteria are monotone in the
                               resources, so this is the same secure/insecure classification);
                               the reported values are the bounds. Used for storage sizing.
    stop_on_insecure:          True: stop at the first insecure hour (rows after it stay empty).

    Returns
    -------
    dict results_frequency_response:    Results of frequency response simulation solved
                                         with the optimal combination of frequency services for each contingency.

    dataframe summary_reserves:         The optimal combination of (H, FFR, FCR, aFRR, mFRR) that ensures frequency security for each contingency.
    dataframe data:                     Contingency, and system data related to the reserve sizing.
    """
    start_time = time.time()  # begin execution timer

    # Parámetros del sistema (ajustar según tus datos)
    f_0 = 50.0  # Frecuencia nominal [Hz]
    S_base = 1000  # Base del sistema [MW]

    # Definition of default settings for the function
    activation_times, activation_times_source = resolve_activation_times(
        activation_times, inputs.get('config'))
    logging.info("get_frequency_stability_reserves(): activation_times source = '%s'" % activation_times_source)

    # Definition of Safe operational limits. limit_freq/limit_rocof default to the same
    # TSO reference constants used in preprocessing to size droop-based reserve
    # participation (build.py -> commons['FrequencyResponseConstants']), so the
    # post-contingency compliance check and the droop-sizing reference are a single
    # source of truth rather than two independently-maintained literals.
    freq_constants = commons['FrequencyResponseConstants']
    if limit_freq is None:
        limit_freq = freq_constants['DeltaFrequencyMax']
    if limit_rocof is None:
        limit_rocof = freq_constants['RoCoF_max']
    if limit_freq_steady_state is None:
        limit_freq_steady_state = freq_constants.get('SteadyStateDeviationMax', 0.2)

    # Definition of activated reserves

    if use_vir is None:
        use_vir = True
    if use_ffr is None:
        use_ffr = True
    if use_fcr is None:
        use_fcr = True
    if use_afrr is None:
        use_afrr = True
    if use_mfrr is None:
        use_mfrr = True

    # if use_vir == True:
    #     limit_rocof = 10000

    # Function settings
    print(limit_freq, limit_rocof, limit_freq_steady_state, use_vir, use_ffr, use_fcr, use_afrr, use_mfrr)

    Damping = inputs["param_df"]["Demand"].filter(like="DA").sum(axis=1).to_frame("Damping")*0.015
    Contingency = results['OutputContingency'].to_frame("Contingency")
    data = pd.concat([Contingency, Damping], axis=1)

    # # Build reduced DataFrame
    # # data_grouped = group_contingencies_data(data, "Contingency", "Damping", tol_max=10.0, tol_min=50.0)
    # data_grouped = group_contingencies_data(data, "Contingency", "Damping", method="hierarchical", tol_max=1.0)
    # group_cols = ["Contingency_group", "Damping_group"]
    # data_reduced = (
    #     data_grouped.groupby(group_cols)
    #     .size()
    #     .reset_index(name="count")
    # )

    units = inputs['units']
    # --- Search bounds (ISSUE-ODESFR-avail) ---------------------------------------------------
    # Hourly bounds from the AVAILABLE fleet: each unit's nameplate contribution is scaled by
    # the same hourly availability the UCED applies to it, LoadMaximum(u,h) = AF(u,h)*(1-OF(u,h))
    # (UCM.gms). A unit on outage, or a PV plant at night, cannot provide inertia or a droop
    # response in that hour, so a requirement sized against it is not deliverable. With
    # fleet_bounds='installed' (previous behaviour) every factor is 1 and the bounds are the
    # nameplate values of Eqs. sir_bound..fcr_bound, constant over the year.
    # Bounds are NOT rounded (ISSUE-ODESFR-floor). Small/inverted brackets are handled in
    # _bisect_collect_feasible.
    if fleet_bounds not in ('available', 'installed'):
        raise ValueError("fleet_bounds must be 'available' or 'installed', got %r" % (fleet_bounds,))
    if fleet_bounds == 'available':
        avail = unit_hourly_availability(inputs, units.index, data.index)
    else:
        avail = pd.DataFrame(1.0, index=data.index, columns=units.index)
    # Synchronous inertia fleet: same Technology classification GAMS uses for
    # cu(au) in EQ_Inertia_Delivery_Limit (build.py -> commons['tech_conventional']).
    sync = units['Technology'].isin(commons['tech_conventional'])
    # Virtual inertia / FFR fleet: batteries only (build.py VIRU/FFRU/FFRD eligibility).
    ibr = units['Technology'].isin(commons['tech_batteries'])
    # FCR fleet: batteries + conventional + renewables (build.py FCR/aFRR/mFRR eligibility).
    fcr_mask = units['Technology'].isin(commons['tech_batteries'] + commons['tech_conventional'] + commons['tech_renewables'])
    # PowerCapacity is a per-unit rating: x Nunits gives the fleet value of a clustered row.
    P = units["PowerCapacity"].astype(float) * units["Nunits"].astype(float)
    HP = units["InertiaConstant"].astype(float).fillna(0.0) * P / S_base     # GWs per unit row
    K = P / (units["Droop"].astype(float) * f_0)                                # MW/Hz per unit row
    K = K.where(np.isfinite(K), 0.0)
    sync_idx, ibr_idx, fcr_idx = sync[sync].index, ibr[ibr].index, fcr_mask[fcr_mask].index
    data["SystemInertia"] = avail.loc[:, sync_idx].mul(HP[sync_idx], axis=1).sum(axis=1).astype(float)
    data["VIR SystemInertia"] = avail.loc[:, ibr_idx].mul(HP[ibr_idx], axis=1).sum(axis=1).astype(float)
    data["FFR Gain"] = avail.loc[:, ibr_idx].mul(K[ibr_idx], axis=1).sum(axis=1).astype(float)
    data["FCR Gain"] = avail.loc[:, fcr_idx].mul(K[fcr_idx], axis=1).sum(axis=1).astype(float)
    # Per-unit saturation of the FCR response (commons['DroopSaturation']): unit u gives K_u*df up
    # to S_u = share_u*P_u, i.e. it saturates at d_u = S_u/K_u = share_u*Droop_u*f_0, which does not
    # depend on availability. Units are grouped by d_u; the group gains K_g(h) follow the hourly
    # availability, so the saturation shape is set hour by hour (kept per index for
    # patch_degenerate_reserves).
    global _FCR_SAT
    cap_tab = commons.get('DroopSaturation')
    sat_by_index = {}
    if cap_tab and len(fcr_idx):
        isbat = units.loc[fcr_idx, 'Technology'].isin(commons['tech_batteries'])
        share = pd.Series([cap_tab['battery'] if b else cap_tab.get(t, cap_tab['default'])
                           for t, b in zip(units.loc[fcr_idx, 'Technology'], isbat)], index=fcr_idx)
        d_u = np.round(share * units.loc[fcr_idx, "Droop"].astype(float) * f_0, 6)
        groups = np.sort(d_u.unique())
        G = np.column_stack([(d_u.values == g).astype(float) for g in groups])     # units x groups
        Kg_h = avail.loc[:, fcr_idx].mul(K[fcr_idx], axis=1).values @ G              # hours x groups
        for pos, ix in enumerate(data.index):
            kmax = float(Kg_h[pos].sum())
            sat_by_index[ix] = (Kg_h[pos].copy(), groups.astype(float), kmax) if kmax > 0 else None
    _FCR_SAT_BY_INDEX.clear()
    _FCR_SAT_BY_INDEX.update(sat_by_index)
    _FCR_SAT = None

    total_contingencies = len(data)  # Count the total number of contingencies"
    contingency_counter = 1  # Initialize the contingnecy counter counter
    tolHs = 1                  # Set a tolerance for the Hs_cap binary search
    tolHv = 1                  # Set a tolerance for the Hv_cap binary search
    tolReserves = 10          # Set a tolerance for the Reserves binary search
    # (tolerances shrink to a quarter of the bracket when it is narrower -- see _bisect_collect_feasible)
    # Gap refinement (ISSUE-ODESFR-ladder): when the smallest candidate of a stage fails downstream,
    # the gap between it and the next larger candidate that passes is bisected to these tolerances
    # (None = no refinement, previous behaviour: accept the next visited value, up to 2x larger).
    refine_tols = [tolHs, tolHv, tolReserves, tolReserves] if refine else [None] * 4
    print(f"Total Contingencies: {total_contingencies}")

    # Create an empty dictionary to store results
    results_frequency_response = {}
    summary_reserves = data.copy()
    columns=['Hs_val','Hv_val','SIR_val','VIR_val','FFR_val','FCR_val','aFRR_val','mFRR_val', 'status',
             'MaxFreqDev', 'MaxRoCoF', 'SteadyStateDev', 'FailedCriterion']
    for col in columns:
        summary_reserves[col] = pd.NA
    print("Initializing sequential binary search over each reserve timeframe...")

    # Perform the binary search for each row of the dataframe data
    for index, row in data.iterrows():
        use_vir = bool(row['VIR SystemInertia'] > 0)
        use_ffr = bool(row['FFR Gain'] > 0)
        use_fcr = bool(row['FCR Gain'] > 0)
        _FCR_SAT = _FCR_SAT_BY_INDEX.get(index)

        # Perform the operations on each row
        print(f"Calculating frequency reserves for Contingency {contingency_counter}")
        # Binary search range and step
        Hs_range = (1, (row['SystemInertia']))
        Hv_range = (1, (row['VIR SystemInertia']))
        FFR_range = (0, row['FFR Gain'])
        FCR_range = (0, row['FCR Gain'])

        def _feasible(sim_time, Hs, Hv, FFR, FCR, aFRR, mFRR):
            max_fd, max_r = frequency_response(
                sim_time, activation_times, Hs, Hv, FFR, FCR, aFRR, mFRR,
                row['Contingency'], row['Damping'])[:2]
            return (max_fd <= limit_freq) and (max_r <= limit_rocof)

        # adjust frequency response simulation settings per stage window
        if use_vir:
            t1 = activation_times["vir"]["prep"]
        elif use_ffr:
            t1 = activation_times["ffr"]["prep"]
        else:
            t1 = activation_times["fcr"]["prep"]
        t2 = activation_times["ffr"]["prep"] if use_ffr else activation_times["fcr"]["prep"]
        # FFR stage window ends when FFR starts to be withdrawn (end of full delivery); a gradual
        # withdrawal after that is covered by the FCR stage and the final Stage-5 check.
        t3 = activation_times["ffr"]["delivery"]
        t4 = activation_times["fcr"]["delivery"]
        t5 = max(activation_times["mfrr"]["ramp"], commons['FrequencyResponseConstants'].get('RecoveryTime', 0)) + 50

        if feasibility_only:
            caps = (row['SystemInertia'], row['VIR SystemInertia'] if use_vir else 0,
                    row['FFR Gain'] if use_ffr else 0, row['FCR Gain'] if use_fcr else 0)
            out = frequency_response(t5, activation_times, *caps, row['Contingency'], row['Contingency'],
                                     row['Contingency'], row['Damping'], True)
            status = bool(use_afrr and use_mfrr and out[0] <= limit_freq and out[1] <= limit_rocof
                          and _steady_state_ok(out[8], limit_freq_steady_state))
            metrics = _response_metrics(out[8], limit_freq, limit_rocof, limit_freq_steady_state)
            summary_reserves.loc[index] = [row['Contingency'], row['Damping'], row['SystemInertia'],
                                           row['VIR SystemInertia'], row['FFR Gain'], row['FCR Gain'],
                                           caps[0], caps[1], out[2], out[3], out[4], out[5], out[6], out[7],
                                           status, *metrics]
            contingency_counter += 1
            if stop_on_insecure and not status:
                break
            continue

        # --- Ordered-backtracking search across Hs -> Hv -> FFR -> FCR -> Stage 5 ---
        # (ODE-SFR-CHANGE-002; DECISION-ODESFR-006 option (a) for Stage 1: each Hs candidate is
        # tested with the full virtual inertia available, Hv = VIR SystemInertia.)
        hv_stage1 = row['VIR SystemInertia'] if use_vir else 0

        def _candidates(k, q):
            if k == 0:      # Stage 1: Hs over [0, t1]
                c = sorted(_bisect_collect_feasible(
                    *Hs_range, tolHs, lambda Hs: _feasible(t1, Hs, hv_stage1, 0, 0, 0, 0)))
                if not c:
                    # No synchronous inertia within the bound secures [0, t1]: continue with the
                    # bound H^S_max (Hs = 0 would return infinite deviations).
                    print("No valid H was found -- using the fleet bound H^S_max")
                    c = [Hs_range[1]]
                return c
            if k == 1:      # Stage 2: Hv over [0, t2]
                c = sorted(_bisect_collect_feasible(
                    *Hv_range, tolHv, lambda Hv: _feasible(t2, q[0], Hv, 0, 0, 0, 0))) if use_vir else []
                # ISSUE-ODESFR-fallback: a stage with no feasible value on its own window falls back to the
                # candidates [0, bound], smallest first (previously [0]): 0 is tried first, as before, then
                # the bound, so services that only secure the hour jointly (e.g. FFR + FCR) are still tested
                # together at their bounds by the later stages and Stage 5.
                return c or ([0, Hv_range[1]] if use_vir else [0])
            if k == 2:      # Stage 3: FFR over [0, t3]
                c = sorted(_bisect_collect_feasible(
                    *FFR_range, tolReserves, lambda FFR: _feasible(t3, q[0], q[1], FFR, 0, 0, 0))) if use_ffr else []
                return c or ([0, FFR_range[1]] if use_ffr else [0])
            c = sorted(_bisect_collect_feasible(  # Stage 4: FCR over [0, t4]
                *FCR_range, tolReserves, lambda FCR: _feasible(t4, q[0], q[1], q[2], FCR, 0, 0))) if use_fcr else []
            return c or ([0, FCR_range[1]] if use_fcr else [0])

        def _accept(q):     # Stage 5: aFRR + mFRR (deterministic, full-contingency test)
            if not (use_afrr and use_mfrr):
                return (0, 0, False, None)
            out = frequency_response(t5, activation_times, q[0], q[1], q[2], q[3],
                                     row['Contingency'], row['Contingency'], row['Contingency'], row['Damping'])
            ok = ((out[0] <= limit_freq) and (out[1] <= limit_rocof)
                  and _steady_state_ok(out[8], limit_freq_steady_state))
            return (row['Contingency'], row['Contingency'], True, out) if ok else None

        found = _ordered_backtracking(4, _candidates, _accept, refine_tols)
        best_df = None
        best_solution = None
        if found is not None:
            (Hs_cap, Hv_cap, FFR_val, FCR_val), (aFRR_val, mFRR_val, status, out) = found
            best_solution = (Hs_cap, Hv_cap, FFR_val, FCR_val, aFRR_val, mFRR_val, status)
            if out is not None:
                max_fd, max_r, max_sir, max_vir, max_ffr, max_fcr, max_afrr, max_mfrr, best_df = out

        if not (use_afrr and use_mfrr):
            print("Reserves aFRR and mFRR are not activated")
        if best_solution is None:
            # Insecure hour: report the best response the fleet can deliver, every service at its bound
            print(f"No fully feasible combination found for Contingency {contingency_counter} "
                  "— reporting the full-fleet response (every service at its bound)")
            best_solution = (row['SystemInertia'], row['VIR SystemInertia'] if use_vir else 0,
                             row['FFR Gain'] if use_ffr else 0, row['FCR Gain'] if use_fcr else 0,
                             row['Contingency'], row['Contingency'], False)

        Hs_cap, Hv_cap, FFR_val, FCR_val, aFRR_val, mFRR_val, status = best_solution
        if best_df is None:
            # Either the aFRR/mFRR check was skipped (not activated) or no
            # combination cleared it -- run the chosen combination once,
            # verbosely, to obtain the diagnostic values and results_df.
            max_fd, max_r, max_sir, max_vir, max_ffr, max_fcr, max_afrr, max_mfrr, best_df = frequency_response(
                t5, activation_times, Hs_cap, Hv_cap, FFR_val, FCR_val, aFRR_val, mFRR_val,
                row['Contingency'], row['Damping'], True)

        # Store the results_df in the dictionary
        results_frequency_response[f"Contingency{contingency_counter}"] = best_df

        print(f"The best reserves combination for Contingency {contingency_counter} is: Hs={Hs_cap:.2f}, Hv={Hv_cap:.2f}, SIR={max_sir:.2f}, VIR={max_vir:.2f}, FFR={max_ffr:.2f}, FCR={max_fcr:.2f}, aFRR={max_afrr:.2f}, mFRR={max_mfrr:.2f}")
        metrics = _response_metrics(best_df, limit_freq, limit_rocof, limit_freq_steady_state)
        summary_reserves.loc[index] = [row['Contingency'], row['Damping'], row['SystemInertia'], row['VIR SystemInertia'], row['FFR Gain'], row['FCR Gain'], Hs_cap, Hv_cap, max_sir, max_vir, max_ffr, max_fcr, max_afrr, max_mfrr, status, *metrics]

        contingency_counter += 1

    end_time = time.time()  # fin de la función
    elapsed_time = end_time - start_time
    print(f"Tiempo total de optimize_search: {elapsed_time:.2f} segundos")
    # # Map back to full time series
    # summary_reserves = summary_reserves.merge(
    #     summary_reserves_reduced,
    #     left_on=group_cols,
    #     right_on=group_cols,
    #     how="left"
    # )
    # We convert the 'index' column into the DataFrame's index.

    # summary_reserves = summary_reserves.set_index('index').sort_index()
    summary_reserves = summary_reserves.reset_index(drop=True)
    summary_reserves.index = data.index
    # data_grouped = data_grouped.set_index('index').sort_index()

    # # Save the results of the swing equation solutions for each contingency
    # with pd.ExcelWriter(path +'results_frecuency_response.xlsx', engine='xlsxwriter') as writer:
    #     for contingency, df in results_frequency_response.items():
    #         df.to_excel(writer, sheet_name=contingency, index=False)

    # Save the results (reserve size) of the frequency security constraints analisys
    columns_to_save  = ['Hs_val','Hv_val','SIR_val','VIR_val', 'FFR_val', 'FCR_val', 'aFRR_val', 'mFRR_val', 'status']
    for col in columns_to_save:
        filename = path +f"{col}.csv"
        summary_reserves[[col]].to_csv(filename, index=True)

    # Save the data containing the Contingency, Damping, and max System Inertia
    # data_grouped.to_csv(path +'Contingency.csv', index=True)
    data.to_csv(path +'Contingency.csv', index=True)

    # # Save the results in pickle file
    # with open(path +"freq_stab_results.pkl", "wb") as f:
    #     # pickle.dump((results_frequency_response, summary_reserves, data_grouped), f)
    #     pickle.dump((results_frequency_response, summary_reserves, data), f)


    # return results_frequency_response, summary_reserves, data_grouped
    return results_frequency_response, summary_reserves, data


def patch_degenerate_reserves(summary_reserves, activation_times=None, config=None,
                               limit_freq=None, limit_rocof=None, limit_freq_steady_state=None):
    """
    Repairs rows where ``get_frequency_stability_reserves``'s own sequential
    search collapsed to a degenerate result (``status=False``, ``Hs_val=0``)
    before ever crediting the fleet's real synchronous inertia, leaving one
    or more of ``VIR_val``/``FFR_val``/``FCR_val`` as ``inf`` -- non-finite
    values that ``write_odesfr_exogenous_csvs`` would otherwise pass straight
    into the DC/SC exogenous CSVs, which then fail Dispa-SET's own
    ``check_VIRDemand``/``check_FFRDemand`` validation (``"... contains
    values higher than demand"``) as soon as a build is attempted.

    Root cause: Stage 1
    of the search tests synchronous inertia Hs ALONE (Hv/FFR/FCR/aFRR/mFRR
    all held at 0) at the very first activation checkpoint (VIR's own
    ``prep`` time, ~0.5s post-contingency by default) -- before virtual
    inertia has been decided at all. For a contingency large enough relative
    to the fleet's own inertia ceiling, the instantaneous RoCoF at that
    checkpoint can exceed the limit even at ``Hs = H^S_max``, by a narrow
    margin, purely because Hv is still zero at this stage (RoCoF depends on
    ``Hs+Hv`` combined). When that happens, ``Hs_candidates`` collapses to
    ``[0]`` and every later stage inherits ``Hs_cap=0``, discarding the
    fleet's real inertia for the rest of the search -- not because the true
    combined system (Hs at its ceiling + Hv + FFR + FCR + aFRR/mFRR, all
    delivered together as a real operator would) is actually infeasible, but
    because the search's own staged, non-joint architecture never gets to
    test that combination. This is a different (and more specific) failure
    mode than the already-documented FFR/FCR cross-stage limitation, though
    both stem from the same staged, non-joint search design (see the module
    docstring's reference to the deliberately-deferred search-loop refactor
    -- this function is a targeted downstream repair, not that refactor).

    For each row with ``status == False`` (or any non-finite delivered
    value), this directly re-evaluates ``frequency_response`` ONCE more with
    every reserve family simultaneously at its own fleet-derived maximum
    bound -- ``Hs=SystemInertia``, ``Hv='VIR SystemInertia'``,
    ``FFR='FFR Gain'``, ``FCR='FCR Gain'`` (unrounded fleet bounds),
    ``aFRR=mFRR=Contingency`` (matching the model's own Stage-5 deterministic
    convention) -- and, if that combination is genuinely feasible, replaces
    the row's degenerate ``Hs_val``/``Hv_val``/``VIR_val``/``FFR_val``/
    ``FCR_val``/``aFRR_val``/``mFRR_val``/``status`` with the real values
    from that check. A row that is still infeasible even with every bound
    used simultaneously is left untouched (a genuine adequacy gap, not a
    search artifact) -- ``write_odesfr_exogenous_csvs`` will still receive
    ``inf`` for it and the downstream build will still (correctly) fail,
    which is the right behavior for a row this function cannot honestly fix.

    :param summary_reserves: DataFrame as returned by
        ``get_frequency_stability_reserves`` (must retain its
        ``'SystemInertia'``, ``'VIR SystemInertia'``, ``'FFR Gain'``,
        ``'FCR Gain'``, ``'Contingency'``, ``'Damping'``, ``'status'``
        columns -- call this on the full frame before any column
        subsetting).
    :param activation_times, config, limit_freq, limit_rocof,
        limit_freq_steady_state: forwarded to ``resolve_activation_times``/
        ``frequency_response`` with the same defaults
        ``get_frequency_stability_reserves`` itself uses, so the re-check is
        evaluated against identical criteria.
    :return: (patched_summary_reserves, patch_log) -- ``patch_log`` is a
        DataFrame indexed by every row this function actually re-evaluated,
        with columns ``['max_fd', 'max_rocof', 'fixed']`` (``fixed=True``
        rows had their values replaced; ``fixed=False`` rows are genuinely
        infeasible even at full bounds and were left untouched), for the
        caller to log what changed rather than have it happen silently.
    """
    activation_times, _source = resolve_activation_times(activation_times, config)
    freq_constants = commons['FrequencyResponseConstants']
    if limit_freq is None:
        limit_freq = freq_constants['DeltaFrequencyMax']
    if limit_rocof is None:
        limit_rocof = freq_constants['RoCoF_max']
    if limit_freq_steady_state is None:
        limit_freq_steady_state = freq_constants.get('SteadyStateDeviationMax', 0.2)
    tolHs, tolHv, tolReserves = 1, 1, 10  # matches get_frequency_stability_reserves's own tolerances
    t2_ffr = activation_times["ffr"]["prep"]
    t2_fcr = activation_times["fcr"]["prep"]
    t3 = activation_times["ffr"]["delivery"]  # end of full FFR delivery (see main search)
    t4 = activation_times["fcr"]["delivery"]
    t5 = max(activation_times["mfrr"]["ramp"], commons['FrequencyResponseConstants'].get('RecoveryTime', 0)) + 50

    patched = summary_reserves.copy()
    degenerate_mask = (~patched['status'].astype(bool)) if 'status' in patched.columns else pd.Series(False, index=patched.index)
    for col in ('Hs_val', 'Hv_val', 'VIR_val', 'FFR_val', 'FCR_val', 'aFRR_val'):
        if col in patched.columns:
            degenerate_mask = degenerate_mask | ~np.isfinite(pd.to_numeric(patched[col], errors='coerce'))

    log_rows = {}
    n_degenerate = int(degenerate_mask.sum())
    if n_degenerate:
        logging.warning("patch_degenerate_reserves(): %d/%d rows are degenerate (status=False or non-finite "
                         "delivered value) -- re-searching each with Hs fixed at the fleet's own ceiling "
                         "(bypassing the isolated, Hv-blind Stage-1 test that discarded it) and Hv/FFR/FCR "
                         "still bisected for their own minimum sufficient value, not jumped to their ceiling.",
                         n_degenerate, len(patched))

    global _FCR_SAT
    _sat_default = _FCR_SAT
    for idx in patched.index[degenerate_mask]:
        row = patched.loc[idx]
        # FCR saturation of this hour's available fleet (set by the search), if known
        _FCR_SAT = _FCR_SAT_BY_INDEX.get(idx, _sat_default)
        contingency = row['Contingency']
        damping = row['Damping']
        hs_cap = row['SystemInertia']  # fixed at the fleet's own ceiling -- see docstring for why
        hv_max = row['VIR SystemInertia'] if pd.notna(row.get('VIR SystemInertia')) else 0.0
        ffr_max = row['FFR Gain'] if pd.notna(row.get('FFR Gain')) else 0.0
        fcr_max = row['FCR Gain'] if pd.notna(row.get('FCR Gain')) else 0.0

        def _feasible(sim_time, Hv, FFR, FCR, aFRR, mFRR):
            max_fd, max_r = frequency_response(
                sim_time, activation_times, hs_cap, Hv, FFR, FCR, aFRR, mFRR, contingency, damping)[:2]
            return (max_fd <= limit_freq) and (max_r <= limit_rocof)

        # --- Stage 2: Hv (virtual inertia), Hs already fixed at ceiling ---
        t2 = t2_ffr if ffr_max > 0 else t2_fcr
        Hv_candidates = sorted(_bisect_collect_feasible(
            1, hv_max, tolHv, lambda Hv: _feasible(t2, Hv, 0, 0, 0, 0))) if hv_max > 0 else []
        if not Hv_candidates:
            # ISSUE-ODESFR-fallback: a stage with no feasible value on its own window falls back to the
            # candidates [0, bound], smallest first (previously [0]): 0 is tried first, as before, then
            # the bound, so services that only secure the hour jointly (e.g. FFR + FCR) are still tested
            # together at their bounds by the later stages and Stage 5.
            Hv_candidates = [0, hv_max] if hv_max > 0 else [0]

        best = None
        for Hv_cap in Hv_candidates:
            # --- Stage 3: FFR, for this Hv_cap ---
            FFR_candidates = sorted(_bisect_collect_feasible(
                0, ffr_max, tolReserves, lambda FFR: _feasible(t3, Hv_cap, FFR, 0, 0, 0))) if ffr_max > 0 else []
            if not FFR_candidates:
                FFR_candidates = [0, ffr_max] if ffr_max > 0 else [0]  # ISSUE-ODESFR-fallback
            for FFR_val in FFR_candidates:
                # --- Stage 4: FCR, for this (Hv_cap, FFR_val) ---
                FCR_candidates = sorted(_bisect_collect_feasible(
                    0, fcr_max, tolReserves,
                    lambda FCR: _feasible(t4, Hv_cap, FFR_val, FCR, 0, 0))) if fcr_max > 0 else []
                if not FCR_candidates:
                    FCR_candidates = [0, fcr_max] if fcr_max > 0 else [0]  # ISSUE-ODESFR-fallback
                for FCR_val in FCR_candidates:
                    # --- Stage 5: deterministic aFRR=mFRR=contingency check ---
                    max_fd, max_r, max_sir, max_vir, max_ffr, max_fcr, max_afrr, max_mfrr, results_df = frequency_response(
                        t5, activation_times, hs_cap, Hv_cap, FFR_val, FCR_val, contingency, contingency, contingency,
                        damping, verbose=False)
                    if (max_fd <= limit_freq) and (max_r <= limit_rocof) and _steady_state_ok(results_df, limit_freq_steady_state):
                        best = (Hv_cap, max_sir, max_vir, max_ffr, max_fcr, max_afrr, max_mfrr, max_fd, max_r)
                        break
                if best is not None:
                    break
            if best is not None:
                break

        if best is None:
            # Smallest-first bisection at every stage still didn't clear Stage
            # 5 (should be rare, given max-everything is already known
            # feasible from the fallback below) -- fall back to every bound
            # simultaneously, the same safe-but-conservative combination this
            # function used before this refinement, rather than leave the row
            # degenerate.
            max_fd, max_r, max_sir, max_vir, max_ffr, max_fcr, max_afrr, max_mfrr, results_df = frequency_response(
                t5, activation_times, hs_cap, hv_max, ffr_max, fcr_max, contingency, contingency, contingency,
                damping, verbose=False)
            fixed = bool((max_fd <= limit_freq) and (max_r <= limit_rocof)
                         and _steady_state_ok(results_df, limit_freq_steady_state))
            if fixed:
                best = (hv_max, max_sir, max_vir, max_ffr, max_fcr, max_afrr, max_mfrr, max_fd, max_r)
                logging.warning("patch_degenerate_reserves(): row %s needed the full-bound fallback (no smaller "
                                 "Hv/FFR/FCR combination cleared Stage 5) -- this one row's requirement is "
                                 "genuinely conservative, not just this function being lazy.", idx)

        if best is None:
            log_rows[idx] = {'max_fd': max_fd, 'max_rocof': max_r, 'fixed': False}
            continue

        hv_cap, max_sir, max_vir, max_ffr, max_fcr, max_afrr, max_mfrr, max_fd, max_r = best
        log_rows[idx] = {'max_fd': max_fd, 'max_rocof': max_r, 'fixed': True}
        patched.loc[idx, 'Hs_val'] = hs_cap
        patched.loc[idx, 'Hv_val'] = hv_cap
        patched.loc[idx, 'VIR_val'] = max_vir
        patched.loc[idx, 'FFR_val'] = max_ffr
        patched.loc[idx, 'FCR_val'] = max_fcr
        patched.loc[idx, 'aFRR_val'] = max_afrr
        patched.loc[idx, 'mFRR_val'] = max_mfrr
        if 'status' in patched.columns:
            patched.loc[idx, 'status'] = True
        if 'MaxFreqDev' in patched.columns:
            # results_df is the response of the combination kept above (the loops stop at it)
            m = _response_metrics(results_df, limit_freq, limit_rocof, limit_freq_steady_state)
            for col, val in zip(('MaxFreqDev', 'MaxRoCoF', 'SteadyStateDev', 'FailedCriterion'), m):
                patched.loc[idx, col] = val

    _FCR_SAT = _sat_default
    patch_log = pd.DataFrame.from_dict(log_rows, orient='index')
    if n_degenerate:
        n_fixed = int(patch_log['fixed'].sum()) if not patch_log.empty else 0
        logging.warning("patch_degenerate_reserves(): %d/%d degenerate rows resolved to a genuinely feasible "
                         "(minimum-sufficient-Hv/FFR/FCR, not blanket-maximum) combination and were patched; "
                         "%d remain infeasible even at every bound simultaneously and were left untouched "
                         "(a real adequacy gap, not a search artifact).",
                         n_fixed, n_degenerate, n_degenerate - n_fixed)
    return patched, patch_log
