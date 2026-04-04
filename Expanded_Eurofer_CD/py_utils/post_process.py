"""
post_process.py — Derived quantities for Expanded_Eurofer_CD.

Computes macroscopic observables from the ODE solution, including:
  - Total SIA / vacancy / He contents
  - Mean cluster sizes
  - Void swelling via the swelling identity (Eq. 161)
  - Frenkel pair conservation diagnostic δ_FP (Eq. 164)
  - He conservation diagnostic δ_He (Eq. 165)

Physics reference
-----------------
Ghoniem, N.M. (2026), Section 7 (Rate_Equations.pdf):
  Eqs. 161-165.
"""

import numpy as np

_kB = 8.617333262e-5


def calculate_derived_quantities(t, y, input_data, rate_eq_obj,
                                 xmax_history=None):
    """
    Compute macroscopic quantities from ODE solution.

    Parameters
    ----------
    t           : ndarray [n_t]
    y           : ndarray [N_eq, n_t]
    input_data  : InputData
    rate_eq_obj : RateEquations or BinMomentRateEquations
    xmax_history: ignored (kept for API compatibility)

    Returns
    -------
    results : dict
    """
    d       = input_data.derived
    Omega   = d['Omega']
    inv_Omega = 1.0 / Omega   # [m^-3] — converts at.frac → m^-3
    r0      = (3.0 * Omega / (4.0 * np.pi)) ** (1.0 / 3.0)

    N   = input_data.N
    M   = input_data.M
    G   = d['G']
    G_He = d['G_He']

    n_t = len(t)
    ns  = np.arange(1, N + 1, dtype=float)
    ms  = np.arange(1, M + 1, dtype=float)

    # Identify state vector layout from rate_eq_obj
    he_mode = getattr(rate_eq_obj, 'he_mode', 'case2')
    is_bin  = hasattr(rate_eq_obj, 'bins')  # BinMomentRateEquations
    qss_He  = getattr(rate_eq_obj, 'qss_He', False)

    i_SIA  = getattr(rate_eq_obj, 'i_SIA', 0)
    i_VAC  = getattr(rate_eq_obj, 'i_VAC', N)
    i_He   = getattr(rate_eq_obj, 'i_He',  None)   # None when qss_He=True
    if i_He is None and not qss_He:
        i_He = y.shape[0] - 1

    # For bin-moment: reconstruct c_n from moments
    K_bins = getattr(rate_eq_obj, 'K', 0)

    # Allocate output arrays
    C_SIA_tot = np.zeros(n_t)   # total SIA content Σ n·c_n
    C_VAC_tot = np.zeros(n_t)   # total vacancy content Σ m·c_m
    C_He_tot  = np.zeros(n_t)   # total He (free + trapped)
    mean_n_i  = np.zeros(n_t)   # mean SIA cluster size
    mean_n_v  = np.zeros(n_t)   # mean vacancy cluster size
    swelling  = np.zeros(n_t)   # void swelling S(t) [fraction]
    N_loops   = np.zeros(n_t)   # number density of SIA loops [a.frac]
    N_voids   = np.zeros(n_t)   # number density of voids [a.frac]
    C_i1      = np.zeros(n_t)   # free SIA monomer
    C_v1      = np.zeros(n_t)   # free vacancy monomer
    C_He_free = np.zeros(n_t)   # free He

    delta_FP  = np.zeros(n_t)   # Frenkel pair conservation (Eq. 164)
    delta_He  = np.zeros(n_t)   # He conservation (Eq. 165)

    # Cascade survival fraction η — used only in δ_FP accounting
    eta = input_data.production_fission.get('eta', 0.30) \
          if 'fiss' in d['spectrum'] \
          else input_data.production_fusion.get('eta', 0.28)

    for j in range(n_t):
        yj = np.maximum(y[:, j], 0.0)

        if is_bin:
            # Reconstruct c_n from moments
            mu0_j = yj[0::2][:K_bins]
            mu1_j = yj[1::2][:K_bins]
            from .bin_moment_rates import distribution_from_moments_pc
            c_n = distribution_from_moments_pc(mu0_j, mu1_j, rate_eq_obj.bins, N)
            c_v = yj[i_VAC:i_VAC + M]
        else:
            c_n = yj[i_SIA:i_VAC]     # [N]
            c_v = yj[i_VAC:i_VAC + M] # [M]

        # Q (He in voids) — extract based on he_mode
        if he_mode == 'case1':
            i_Q = getattr(rate_eq_obj, 'i_Q', i_VAC + M)
            # When qss_He, i_He is None so state ends at i_Q + M
            Q_end = (i_He if i_He is not None else i_Q + M)
            Q_m   = yj[i_Q:Q_end]    # [M] per-class He
            Q_tot = np.sum(Q_m)
        else:
            i_Qtot = getattr(rate_eq_obj, 'i_Qtot', i_VAC + M)
            Q_tot  = yj[i_Qtot]
            Q_m    = None

        # Free He: read from state or reconstruct from QSS
        if qss_He:
            c_v_pp = yj[i_VAC:i_VAC + M]
            c_h = rate_eq_obj.compute_c_h_qss(c_v_pp, Q_tot=Q_tot, Q_m=Q_m)
        else:
            c_h = yj[i_He]

        # SIA content
        C_SIA_tot[j] = np.dot(ns, c_n)
        C_VAC_tot[j] = np.dot(ms, c_v)
        C_He_tot[j]  = c_h + Q_tot

        # Mean cluster sizes (weighted average)
        sum_cni = np.sum(c_n[1:])   # clusters n ≥ 2
        mean_n_i[j] = np.dot(ns[1:], c_n[1:]) / max(sum_cni, 1e-300)

        sum_cvi = np.sum(c_v[1:])   # clusters m ≥ 2
        mean_n_v[j] = np.dot(ms[1:], c_v[1:]) / max(sum_cvi, 1e-300)

        # Number densities (clusters n,m ≥ 2)
        N_loops[j] = sum_cni
        N_voids[j] = sum_cvi

        # Void swelling (Eq. 161): S = Σ_m m·c_m
        # V_m = m·Ω per cluster; c_m is cluster number density in at.frac;
        # ΔV/V = Σ_m (m·c_m·Ω)/(N_atoms·Ω) = Σ_m m·c_m = C_VAC_tot
        swelling[j] = C_VAC_tot[j]

        C_i1[j]      = c_n[0]
        C_v1[j]      = c_v[0]
        C_He_free[j] = c_h

        # Frenkel pair balance diagnostic (Eq. 164)
        # δ_FP = |Σn·n·c_n − Σm·m·c_m| / (η·G·t)
        # Tests imbalance between SIA and vacancy cluster content.
        # Should be near 0 when equal numbers of SIA and vacancies are in clusters
        # (unbiased sinks); rises toward 1 only if one species is severely depleted
        # relative to the other.  This is NOT a sink-fraction — see post_process note.
        surviving = eta * G * max(t[j], 1e-20)
        delta_FP[j] = abs(C_SIA_tot[j] - C_VAC_tot[j]) / max(surviving, 1e-300)

        # He conservation diagnostic (Eq. 165)
        # δ_He = |c_h + Q_tot − G_He·t| / (G_He·t)
        # Approaches 0 when He is fully retained in voids/free state;
        # approaches 1 when He is dominated by fixed-sink absorption.
        denom_he = G_He * max(t[j], 1e-20)
        delta_He[j] = abs(C_He_tot[j] - G_He * t[j]) / max(denom_he, 1e-300)

    # Dose axis [dpa]
    dose = G * t

    # Convert concentrations from atomic fraction to m^-3
    # c [at.frac] / Omega [m^3/atom] = N [m^-3]
    # swelling, mean sizes, and conservation diagnostics remain dimensionless.
    C_SIA_tot  *= inv_Omega
    C_VAC_tot  *= inv_Omega
    C_He_tot   *= inv_Omega
    C_He_free  *= inv_Omega
    N_loops    *= inv_Omega
    N_voids    *= inv_Omega
    C_i1       *= inv_Omega
    C_v1       *= inv_Omega

    # Concentration names for quality checks
    conc_names = ([f'Ci_{n}' for n in range(1, N + 1)] +
                  [f'Cv_{m}' for m in range(1, M + 1)] + ['C_He'])

    results = {
        't':          t,
        'dose':       dose,
        'y':          y,           # raw ODE state, still in at.frac
        'Omega':      Omega,       # [m^3] — for downstream unit conversions
        'C_SIA_tot':  C_SIA_tot,   # [m^-3]
        'C_VAC_tot':  C_VAC_tot,   # [m^-3]
        'C_He_tot':   C_He_tot,    # [m^-3]
        'C_He_free':  C_He_free,   # [m^-3]
        'mean_n_i':   mean_n_i,    # [atoms]
        'mean_n_v':   mean_n_v,    # [atoms]
        'N_loops':    N_loops,     # [m^-3]
        'N_voids':    N_voids,     # [m^-3]
        'swelling':   swelling,    # [dimensionless fraction]
        'C_i1':       C_i1,        # [m^-3]
        'C_v1':       C_v1,        # [m^-3]
        'delta_FP':   delta_FP,    # [dimensionless]
        'delta_He':   delta_He,    # [dimensionless]
        'conc_names': conc_names,
    }
    return results


def summary_csv_row(results, input_data, solver_label=''):
    """
    Build a summary dict for a single run (for output/summary.csv).
    """
    d = input_data.derived
    t_end = results['t'][-1]
    row = {
        'solver':        solver_label,
        'physics_option':input_data.physics_option,
        'T_K':           d['T'],
        'G_dpa_s':       d['G'],
        't_end_s':       t_end,
        'dose_dpa':      results['dose'][-1],
        'C_SIA_tot_m3':  results['C_SIA_tot'][-1],
        'C_VAC_tot_m3':  results['C_VAC_tot'][-1],
        'C_He_tot_m3':   results['C_He_tot'][-1],
        'mean_n_i':      results['mean_n_i'][-1],
        'mean_n_v':      results['mean_n_v'][-1],
        'swelling':      results['swelling'][-1],
        'delta_FP':      results['delta_FP'][-1],
        'delta_He':      results['delta_He'][-1],
    }
    return row
