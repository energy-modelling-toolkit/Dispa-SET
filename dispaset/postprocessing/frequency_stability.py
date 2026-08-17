# -*- coding: utf-8 -*-
"""
Frequency Response Assessment Model (FRAM): sequential, ordered-backtracking
sizing of synchronous inertia, virtual inertia, and frequency reserves (FFR,
FCR, aFRR, mFRR) against RoCoF, nadir, and steady-state frequency-deviation
limits, per contingency.

Extracted out of postprocessing.py into its own module, mirroring how
reserves.py, boundary_sector.py, geoplot.py, and streamlit_functions.py are
each split out of the larger preprocessing/postprocessing files. Pure move,
no behavior change.

See FRAM_ISSUES.md and FRAM_DECISIONS.md at the repository root for the full
change history, known issues, and open decisions for this subsystem.
"""
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
def frequency_response(sim_time, activation_times, Hs_cap, Hv_cap, FFR_cap, FCR_cap, aFRR_cap, mFRR_cap,
                   contingency, D_local, verbose=False):
    """
    This function solves the power swing differential equation, for each combination
    of inertia and frequency reserves desired.
    param sim_time:             Time to evaluate the differential equation in seconds [s]
    param activation_times:     Activation times for each reserve in absolut values from sim_time = 0
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
        # FRAM-CHANGE-001 (ISSUE-FRAM-003): this branch must return the same
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

    # definition of simulation horizon and time step
    t = np.arange(0, sim_time, 0.1)

    # Precompute weights vector for plotting and for fixed deployment logic.
    # The "vir" weight (first return value) is intentionally discarded: it is
    # never used below either -- H_total combines Hs_cap+Hv_cap unweighted
    # (see ISSUE-FRAM-004 / DECISION-FRAM-006 for the separate, unresolved
    # question of whether Hv should be gated by its own activation profile;
    # this fix does not touch that).
    _, W_ffr_vec, W_fcr_vec, W_afrr_vec, W_mfrr_vec = compute_weights(t, activation_times)

    f_0 = 50.0  # Nominal frequency [Hz]
    S_base = 1000  # System base [MW]

    def state(y, tt):
        f = y[0]
        contingency_local = contingency if tt >= 1.0 else 0.0

        # FRAM-CHANGE-003 (ISSUE-FRAM-002): interpolate the already-precomputed
        # weight vectors instead of recomputing compute_weights() from scratch
        # at every internal ODE step (measured ~1,470 recomputations per call
        # before this fix; interpolation gives numerically equivalent output
        # -- max diff ~1e-8 to 1e-12 Hz across benchmarked cases -- at a
        # measured 13-16x speedup, see FRAM_ISSUES.md ISSUE-FRAM-002).
        w_ffr = np.interp(tt, t, W_ffr_vec)
        w_fcr = np.interp(tt, t, W_fcr_vec)
        w_afrr = np.interp(tt, t, W_afrr_vec)
        w_mfrr = np.interp(tt, t, W_mfrr_vec)

        ffr = FFR_cap * w_ffr * f
        fcr = FCR_cap * w_fcr * f
        afrr = aFRR_cap * w_afrr
        mfrr = mFRR_cap * w_mfrr

        H_total = Hs_cap + Hv_cap

        deltap = contingency_local - (ffr + fcr + afrr + mfrr) - D_local * f
        dfdt = (deltap / (2 * H_total * S_base)) * f_0

        return [dfdt]

    sol = odeint(state, [0.0], t)
    f = sol[:, 0]

    contingency_vec = np.where(t >= 1.0, contingency, 0.0)

    max_freq_dev = np.max(np.abs(f))
    rocof = np.gradient(f, t)
    max_rocof = np.max(np.abs(rocof))

    sir = 2 * Hs_cap * S_base * rocof/ f_0
    vir = 2 * Hv_cap * S_base * rocof/ f_0
    ffr = FFR_cap * f * W_ffr_vec
    fcr = FCR_cap * f * W_fcr_vec
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
    already used by FRAM's stage bisections).

    Used by get_frequency_stability_reserves (FRAM-CHANGE-002, ISSUE-FRAM-001 /
    DECISION-FRAM-002) to generate the ordered candidate list for each of the
    four bisected stages (Hs, Hv, FFR, FCR), so the calling loop can try
    candidates smallest-first and backtrack to a larger one only if every
    downstream possibility for the smaller one is exhausted.
    """
    feasible = []
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


# %% frequency stability reserves
def get_frequency_stability_reserves(path, inputs, results, activation_times=None,
                                     limit_freq=None, limit_rocof=None, limit_freq_steady_state=None,
                                     use_vir=None, use_ffr=None, use_fcr=None, use_afrr=None, use_mfrr=None):
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
    if activation_times is None:
        activation_times = {
                    "vir": dict(prep=1.1, ramp=1.15, delivery=1.5, deact=2),
                    "ffr": dict(prep=1.5, ramp=2, delivery=7, deact=16),
                    "fcr": dict(prep=1.1, ramp=13.1, delivery=181, deact=301),
                    "afrr": dict(prep=31, ramp=301, delivery=901, deact=1801),
                    "mfrr": dict(prep=901, ramp=1801)  # mfrr leght is considered for the whole timestep
                }

    # Definition of Safe operational limits
    if limit_freq is None:
        limit_freq = 0.8
    if limit_rocof is None:
        limit_rocof = 0.5
    if limit_freq_steady_state is None:
        limit_freq_steady_state = 0.2

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
    # find max system inertia possible from synchronous generators
    # filter synchronous generators by technology eligibility
    # Synchronous inertia fleet: same Technology classification GAMS uses for
    # cu(au) in EQ_Inertia_Delivery_Limit (build.py:679, 889 -> commons['tech_conventional']).
    sync_mask = units['Technology'].isin(commons['tech_conventional'])
    syncunits = units[sync_mask]
    system_inertia = np.floor((syncunits["InertiaConstant"] * syncunits["PowerCapacity"]).sum() / S_base)
    data["SystemInertia"] = system_inertia

    # find max system inertia possible from IBR
    # filter synchronous generators by technology eligibility
    # Virtual inertia / FFR fleet: batteries only, matching build.py's VIRU/FFRU/FFRD
    # eligibility (build.py:682, 890, 1489, 1494 -> commons['tech_batteries']).
    ibr_mask = units['Technology'].isin(commons['tech_batteries'])
    ibrunits = units[ibr_mask]
    vir_system_inertia = np.floor((ibrunits["InertiaConstant"] * ibrunits["PowerCapacity"]).sum() / S_base)
    data["VIR SystemInertia"] = vir_system_inertia


    # find max FFR gain possible from IBR
    # filter synchronous generators by technology eligibility
    ffr_gain = np.floor((ibrunits["PowerCapacity"] / (ibrunits["Droop"] * f_0)).sum())
    data["FFR Gain"] = ffr_gain

    # find max FCR gain possible from IBR+CONV
    # filter synchronous generators by technology eligibility
    # FCR fleet: batteries + conventional + renewables, matching build.py's
    # FCR/aFRR/mFRR eligibility (build.py:1497-1502).
    fcr_mask = units['Technology'].isin(commons['tech_batteries'] + commons['tech_conventional'] + commons['tech_renewables'])
    fcrunits = units[fcr_mask]
    fcr_gain = np.floor((fcrunits["PowerCapacity"] / (fcrunits["Droop"] * f_0)).sum())
    data["FCR Gain"] = fcr_gain

    # total_contingencies = len(data_reduced)  # Count the total number of contingencies"
    total_contingencies = len(data)  # Count the total number of contingencies"
    contingency_counter = 1  # Initialize the contingnecy counter counter
    tolHs = 1                  # Set a tolerance for the Hs_cap binary search
    tolHv = 1                  # Set a tolerance for the Hv_cap binary search
    tolReserves = 10          # Set a tolerance for the Reserves binary search
    print(f"Total Contingencies: {total_contingencies}")

    # Create an empty dictionary to store results
    results_frequency_response = {}
    # Create an empty dataframe to store the reserves found
    # summary_reserves = data_grouped.copy()
    # summary_reserves_reduced  = data_reduced.copy()
    summary_reserves = data.copy()
    # summary_reserves_reduced  = data.copy()
    columns=['Hs_val','Hv_val','SIR_val','VIR_val','FFR_val','FCR_val','aFRR_val','mFRR_val', 'status']
    for col in columns:
        summary_reserves[col] = pd.NA
    print("Initializing sequential binary search over each reserve timeframe...")

    # Perform the binary search for each row of the dataframe data
    # for index, row in data_reduced.iterrows():
    for index, row in data.iterrows():
        if (row['VIR SystemInertia']) > 0:
            use_vir = True
        else:
            use_vir = False
        if (row['FFR Gain']) > 0:
            use_ffr = True
        else:
            use_ffr = False
        if (row['FCR Gain']) > 0:
            use_fcr = True
        else:
            use_fcr = False

        # Perform the operations on each row
        print(f"Calculating frequency reserves for Contingency {contingency_counter}")
        # Binary search range and step
        Hs_range = (1, (row['SystemInertia']))
        Hv_range = (1, (row['VIR SystemInertia']))
        # reserve_range = (0, row['Contingency_group']*1.3)
        # reserve_range = (0, row['Contingency']*1.3)
        FFR_range = (0, row['FFR Gain'])
        FCR_range = (0, row['FCR Gain'])

        # --- Ordered-backtracking search across Hs -> Hv -> FFR -> FCR ---
        # FRAM-CHANGE-002 (ISSUE-FRAM-001, DECISION-FRAM-002): candidates at
        # each stage are generated by bisection (as before) and tried
        # smallest-first. The search only advances to the next candidate at a
        # given stage once the ENTIRE downstream chain -- every later stage,
        # AND the final aFRR/mFRR steady-state check -- has been proven
        # infeasible for every candidate tried so far at that stage. This
        # replaces the previous behavior, where every downstream stage was
        # evaluated for every upstream candidate and the result was kept via
        # plain reassignment inside a loop, so only whichever candidate
        # happened to be processed last survived (not necessarily the
        # smallest, and not deliberately selected by any criterion).
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
        t3 = activation_times["ffr"]["deact"]
        t4 = activation_times["fcr"]["delivery"]
        t5 = activation_times["mfrr"]["ramp"] + 50

        # --- Stage 1: Hs (synchronous inertia) ---
        Hs_candidates = sorted(_bisect_collect_feasible(
            *Hs_range, tolHs, lambda Hs: _feasible(t1, Hs, 0, 0, 0, 0, 0)))
        if not Hs_candidates:
            print("No valid H was found")
            Hs_candidates = [0]  # assume 0 as the default value

        best_solution = None
        best_df = None
        Hv_candidates = FFR_candidates = FCR_candidates = [0]

        for Hs_cap in Hs_candidates:
            # --- Stage 2: Hv (virtual inertia), for this Hs_cap ---
            if use_vir:
                Hv_candidates = sorted(_bisect_collect_feasible(
                    *Hv_range, tolHv, lambda Hv: _feasible(t2, Hs_cap, Hv, 0, 0, 0, 0)))
                if not Hv_candidates:
                    Hv_candidates = [0]
            else:
                Hv_candidates = [0]

            for Hv_cap in Hv_candidates:
                # --- Stage 3: FFR, for this (Hs_cap, Hv_cap) ---
                if use_ffr:
                    FFR_candidates = sorted(_bisect_collect_feasible(
                        *FFR_range, tolReserves,
                        lambda FFR: _feasible(t3, Hs_cap, Hv_cap, FFR, 0, 0, 0)))
                    if not FFR_candidates:
                        FFR_candidates = [0]
                else:
                    FFR_candidates = [0]

                for FFR_val in FFR_candidates:
                    # --- Stage 4: FCR, for this (Hs_cap, Hv_cap, FFR_val) ---
                    if use_fcr:
                        FCR_candidates = sorted(_bisect_collect_feasible(
                            *FCR_range, tolReserves,
                            lambda FCR: _feasible(t4, Hs_cap, Hv_cap, FFR_val, FCR, 0, 0)))
                        if not FCR_candidates:
                            FCR_candidates = [0]
                    else:
                        FCR_candidates = [0]

                    for FCR_val in FCR_candidates:
                        # --- Stage 5: aFRR + mFRR (deterministic, full-contingency test) ---
                        if use_afrr and use_mfrr:
                            max_fd, max_r, max_sir, max_vir, max_ffr, max_fcr, max_afrr, max_mfrr, results_df = frequency_response(
                                t5, activation_times, Hs_cap, Hv_cap, FFR_val, FCR_val,
                                row['Contingency'], row['Contingency'], row['Contingency'], row['Damping'])
                            freq_window = results_df.loc[results_df['Time [s]'] > 350]
                            freq_steady_state = freq_window['Frequency Deviation [Hz]'].abs().max()
                            status = ((max_fd <= limit_freq) and (max_r <= limit_rocof)
                                      and (freq_steady_state <= limit_freq_steady_state))
                            if status:
                                best_solution = (Hs_cap, Hv_cap, FFR_val, FCR_val,
                                                  row['Contingency'], row['Contingency'], True)
                                best_df = results_df
                                break
                        else:
                            best_solution = (Hs_cap, Hv_cap, FFR_val, FCR_val, 0, 0, False)
                            break
                    if best_solution is not None:
                        break
                if best_solution is not None:
                    break
            if best_solution is not None:
                break

        if not (use_afrr and use_mfrr):
            print("Reserves aFRR and mFRR are not activated")
        if best_solution is None:
            print(f"No fully feasible combination found for Contingency {contingency_counter} "
                  "— reporting max-effort attempt (largest candidate tried at every stage)")
            best_solution = (Hs_candidates[-1], Hv_candidates[-1], FFR_candidates[-1],
                              FCR_candidates[-1], 0, 0, False)

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
        summary_reserves.loc[index] = [row['Contingency'], row['Damping'], row['SystemInertia'], row['VIR SystemInertia'], row['FFR Gain'], row['FCR Gain'], Hs_cap, Hv_cap, max_sir, max_vir, max_ffr, max_fcr, max_afrr, max_mfrr, status]

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
    columns_to_save  = ['Hs_val','Hv_val','SIR_val','VIR_val', 'FFR_val', 'FCR_val', 'aFRR_val', 'mFRR_val']
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
