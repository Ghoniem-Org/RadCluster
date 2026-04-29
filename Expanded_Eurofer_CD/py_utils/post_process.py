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

    I   = input_data.I
    V   = input_data.V
    G   = d['G']
    G_He = d['G_He']

    n_t = len(t)
    ns  = np.arange(1, I + 1, dtype=float)
    ms  = np.arange(1, V + 1, dtype=float)

    # Identify state vector layout from rate_eq_obj
    he_mode = getattr(rate_eq_obj, 'he_mode', 'case2')
    is_bin  = hasattr(rate_eq_obj, 'bins')  # BinMomentRateEquations
    qss_He  = getattr(rate_eq_obj, 'qss_He', False)

    i_SIA  = getattr(rate_eq_obj, 'i_SIA', 0)
    i_VAC  = getattr(rate_eq_obj, 'i_VAC', I)
    i_He   = getattr(rate_eq_obj, 'i_He',  None)   # None when qss_He=True
    if i_He is None and not qss_He:
        i_He = y.shape[0] - 1

    # For bin-moment: reconstruct c_n from SIA moments
    I_bin    = getattr(rate_eq_obj, 'I_bin', getattr(rate_eq_obj, 'K', 0))
    # Vacancy bins (V_bin > 0 means vacancies are also binned)
    V_bin    = getattr(rate_eq_obj, 'V_bin', getattr(rate_eq_obj, 'K_v', 0))

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

    delta_FP     = np.zeros(n_t)   # Frenkel pair conservation (Eq. 164) — max(sia, vac)
    delta_FP_sia = np.zeros(n_t)   # SIA arm of FP conservation
    delta_FP_vac = np.zeros(n_t)   # VAC arm of FP conservation
    delta_He  = np.zeros(n_t)   # He conservation (Eq. 165)

    # (Conservation fluxes are now tracked by CVODE as extra state variables)

    # Rate constants for sink terms
    rr = getattr(rate_eq_obj, 'rr', None)
    v_mobile = d.get('v_mobile', d.get('m_max_v', 1))

    # Cascade survival fraction η
    eta = input_data.production_fission.get('eta', 0.30) \
          if 'fiss' in d['spectrum'] \
          else input_data.production_fusion.get('eta', 0.28)

    # Vacancy bin midpoints (for algebraic mu1 when vacancies are binned)
    vac_mid = getattr(rate_eq_obj, 'vac_mid', None)

    for j in range(n_t):
        yj = np.maximum(y[:, j], 0.0)

        if is_bin:
            from .bin_moment_rates import reconstruct_distribution
            i_d = getattr(rate_eq_obj, 'i_discrete', 0)
            v_d = getattr(rate_eq_obj, 'v_discrete', 0)
            P   = getattr(rate_eq_obj, 'n_mom', 2)
            sf  = getattr(rate_eq_obj, 'shape_function', 'linear')

            # ── SIA: discrete + binned ───────────────────────────────
            # Discrete sizes contribute directly to content
            ns_disc = np.arange(1, i_d + 1, dtype=float)
            SIA_content_from_mu1 = np.dot(ns_disc, yj[:i_d])
            # Binned: use tracked first moments μ₁ directly (or approximate)
            if I_bin > 0:
                mom = yj[i_d:i_d + P * I_bin]
                mu0_j = mom[0::P][:I_bin]
                if P >= 2:
                    mu1_j = mom[1::P][:I_bin]
                    SIA_content_from_mu1 += np.sum(mu1_j)
                else:
                    mu1_j = None
                    # Approximate μ₁ as μ₀ × midpoint
                    for kb, (nlo, nhi) in enumerate(rate_eq_obj.bins):
                        SIA_content_from_mu1 += mu0_j[kb] * (nlo + nhi - 1) / 2.0
                mu2_j = mom[2::P][:I_bin] if P >= 3 else None
                # Reconstruct c_n for mean-size and number-density calcs
                c_n = reconstruct_distribution(sf, mu0_j, mu1_j, mu2_j,
                                               rate_eq_obj.bins, I)
                c_n[:i_d] = yj[:i_d]  # overwrite discrete region
            else:
                c_n = np.zeros(I)
                c_n[:i_d] = yj[:i_d]

            # ── Vacancy: discrete + binned ───────────────────────────
            ms_disc = np.arange(1, v_d + 1, dtype=float)
            c_v = np.zeros(V)
            c_v[:v_d] = yj[i_VAC:i_VAC + v_d]
            VAC_content_from_mu1 = np.dot(ms_disc, c_v[:v_d])
            if V_bin > 0:
                vac_start = i_VAC + v_d
                vmom = yj[vac_start:vac_start + P * V_bin]
                vmu0 = vmom[0::P][:V_bin]
                if P >= 2:
                    vmu1 = vmom[1::P][:V_bin]
                    VAC_content_from_mu1 += np.sum(vmu1)
                else:
                    vmu1 = None
                    for kb, (mlo, mhi) in enumerate(rate_eq_obj.vac_bins):
                        VAC_content_from_mu1 += vmu0[kb] * (mlo + mhi - 1) / 2.0
                vmu2 = vmom[2::P][:V_bin] if P >= 3 else None
                c_v_binned = reconstruct_distribution(sf, vmu0, vmu1, vmu2,
                                                      rate_eq_obj.vac_bins, V)
                c_v[v_d:] = c_v_binned[v_d:]
        else:
            c_n = yj[i_SIA:i_VAC]     # [N]
            c_v = yj[i_VAC:i_VAC + V] # [V]
            SIA_content_from_mu1 = None
            VAC_content_from_mu1 = None

        # Q (He in voids) — extract based on he_mode
        if he_mode == 'case1':
            i_Q = getattr(rate_eq_obj, 'i_Q', i_VAC + V)
            n_Q = V_bin if V_bin > 0 else V
            Q_m   = yj[i_Q:i_Q + n_Q]
            Q_tot = np.sum(Q_m)
        else:
            i_Qtot = getattr(rate_eq_obj, 'i_Qtot', i_VAC + V)
            Q_tot  = yj[i_Qtot]
            Q_m    = None

        # Free He: read from state or reconstruct from QSS
        if qss_He:
            c_h = rate_eq_obj.compute_c_h_qss(c_v, Q_tot=Q_tot, Q_m=Q_m)
        else:
            c_h = yj[i_He]

        # SIA and vacancy total content
        # Prefer tracked moments (exact) over PC reconstruction (approximate)
        C_SIA_tot[j] = SIA_content_from_mu1 if SIA_content_from_mu1 is not None \
                        else np.dot(ns, c_n)
        # Tracked mu1 gives exact VAC content for both discrete and binned
        C_VAC_tot[j] = VAC_content_from_mu1 if VAC_content_from_mu1 is not None \
                        else np.dot(ms, c_v)
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


    # Dose axis [dpa]
    dose = G * t

    # ── Conservation diagnostics — read cumulative sinks from state vector ──
    # The last 4 entries of y are CVODE-integrated cumulative sink fluxes:
    #   y[N_eq-5] = J_SIA_fixed:  cumulative SIA content to fixed sinks
    #   y[N_eq-4] = J_SIA_mutual: cumulative SIA to recombination + cavity
    #   y[N_eq-3] = J_VAC_fixed:  cumulative VAC content to fixed sinks
    #   y[N_eq-2] = J_VAC_mutual: cumulative VAC to recombination + cavity
    #   y[N_eq-1] = J_He_sink:    cumulative He to sinks
    N_eq = y.shape[0]
    J_SIA_fixed  = np.maximum(y[N_eq - 5, :], 0.0)
    J_SIA_mutual = np.maximum(y[N_eq - 4, :], 0.0)
    J_VAC_fixed  = np.maximum(y[N_eq - 3, :], 0.0)
    J_VAC_mutual = np.maximum(y[N_eq - 2, :], 0.0)
    J_He_sink    = np.maximum(y[N_eq - 1, :], 0.0)

    # FP conservation:
    #   η·G·t = C_SIA + J_SIA_fixed + J_SIA_mutual  (SIA balance)
    #   η·G·t = C_VAC + J_VAC_fixed + J_VAC_mutual  (VAC balance)
    # δ_FP = max of relative errors
    for j in range(n_t):
        prod = eta * G * t[j]
        if prod > 1e-300:
            err_sia = abs(prod - C_SIA_tot[j] - J_SIA_fixed[j] - J_SIA_mutual[j]) / prod
            err_vac = abs(prod - C_VAC_tot[j] - J_VAC_fixed[j] - J_VAC_mutual[j]) / prod
            delta_FP[j]     = max(err_sia, err_vac)
            delta_FP_sia[j] = err_sia
            delta_FP_vac[j] = err_vac
        else:
            delta_FP[j] = delta_FP_sia[j] = delta_FP_vac[j] = 0.0

    # He conservation: G_He·t = C_He_tot + J_He_sink
    # NOTE: This identity is exact only in dynamic (non-QSS) He mode.
    # In quasi_steady_state mode, the algebraic c_h eliminates the transient
    # storage term (dc_h/dt), so the identity does not hold exactly.
    # When qss_He=True, delta_He is set to NaN to indicate it is not applicable.
    for j in range(n_t):
        he_prod = G_He * t[j]
        if qss_He:
            delta_He[j] = float('nan')
        elif he_prod > 1e-300:
            delta_He[j] = abs(he_prod - C_He_tot[j] - J_He_sink[j]) / he_prod
        else:
            delta_He[j] = 0.0

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
    conc_names = ([f'Ci_{n}' for n in range(1, I + 1)] +
                  [f'Cv_{m}' for m in range(1, V + 1)] + ['C_He'])

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
        'delta_FP':     delta_FP,     # FP conservation: max(sia, vac) relative error
        'delta_FP_sia': delta_FP_sia, # SIA arm of FP conservation
        'delta_FP_vac': delta_FP_vac, # VAC arm of FP conservation
        'delta_He':   delta_He,    # He conservation (NaN when qss_He=True)
        'J_SIA_fixed':  J_SIA_fixed,   # cumulative SIA to fixed sinks [at.frac]
        'J_SIA_mutual': J_SIA_mutual,  # cumulative SIA to recomb + cavity [at.frac]
        'J_VAC_fixed':  J_VAC_fixed,   # cumulative VAC to fixed sinks [at.frac]
        'J_VAC_mutual': J_VAC_mutual,  # cumulative VAC to recomb + cavity [at.frac]
        'J_He_sink':    J_He_sink,     # cumulative He to sinks [at.frac]
        'eta_G':        eta * G,       # FP production rate η·G [at.frac/s]
        'conc_names': conc_names,
        # Layout indices for cross-segment y merging during domain doubling
        '_y_i_VAC':  i_VAC,
        '_y_i_He':   i_VAC + (getattr(rate_eq_obj, 'v_discrete', 0)
                              + getattr(rate_eq_obj, 'n_mom', 2)
                              * getattr(rate_eq_obj, 'V_bin',
                                        getattr(rate_eq_obj, 'K_v', 0))
                              if is_bin else V),
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
