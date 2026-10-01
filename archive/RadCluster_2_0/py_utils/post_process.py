"""
post_process.py — Derived quantities for RadCluster_2_0.

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
                                 xmax_history=None, y_sia100=None):
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
    C_floor = float(input_data.reactions.get('C_floor', 1e-15))

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

    # Cascade survival fraction η — select on the canonical cascade axis
    # ('fission' | 'fusion'), not a substring match on the spectrum label.
    cascade = getattr(input_data, 'cascade', None)
    if cascade is None:
        from .input_data import split_physics_option
        cascade = split_physics_option(input_data.physics_option)[1]
    eta = input_data.production_fission.get('eta', 0.30) \
          if cascade == 'fission' \
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

        # Mean cluster sizes and number densities — computed from tracked
        # moments (μ₀, μ₁) to stay consistent with C_SIA_tot / C_VAC_tot.
        # Floor subtraction is applied at the moment level (exact, no
        # reconstruction): for a bin spanning [nlo, nhi) of width Δ, the
        # uniform C_floor IC contributes Δ·C_floor to μ₀ and
        # C_floor·Δ·(nlo+nhi-1)/2 to μ₁.
        if is_bin:
            # SIA: discrete (n=2..i_d) + binned
            disc_i_eff   = np.maximum(yj[1:i_d] - C_floor, 0.0)
            count_i      = float(np.sum(disc_i_eff))
            content_i    = float(np.dot(ns_disc[1:], disc_i_eff))
            if I_bin > 0:
                for kb, (nlo, nhi) in enumerate(rate_eq_obj.bins):
                    width   = float(nhi - nlo)
                    sum_n   = width * (nlo + nhi - 1) / 2.0
                    mu0_eff = max(float(mu0_j[kb]) - width * C_floor, 0.0)
                    if P >= 2:
                        mu1_eff = max(float(mu1_j[kb]) - C_floor * sum_n, 0.0)
                    else:
                        mu1_eff = mu0_eff * (nlo + nhi - 1) / 2.0
                    count_i   += mu0_eff
                    content_i += mu1_eff
            N_loops[j]  = count_i
            mean_n_i[j] = content_i / count_i if count_i > 0.0 else 0.0

            # VAC: discrete (m=2..v_d) + binned
            disc_v_eff   = np.maximum(yj[i_VAC + 1:i_VAC + v_d] - C_floor, 0.0)
            count_v      = float(np.sum(disc_v_eff))
            content_v    = float(np.dot(ms_disc[1:], disc_v_eff))
            if V_bin > 0:
                for kb, (mlo, mhi) in enumerate(rate_eq_obj.vac_bins):
                    width   = float(mhi - mlo)
                    sum_m   = width * (mlo + mhi - 1) / 2.0
                    mu0_eff = max(float(vmu0[kb]) - width * C_floor, 0.0)
                    if P >= 2:
                        mu1_eff = max(float(vmu1[kb]) - C_floor * sum_m, 0.0)
                    else:
                        mu1_eff = mu0_eff * (mlo + mhi - 1) / 2.0
                    count_v   += mu0_eff
                    content_v += mu1_eff
            N_voids[j]  = count_v
            mean_n_v[j] = content_v / count_v if count_v > 0.0 else 0.0
        else:
            c_n_eff = np.maximum(c_n[1:] - C_floor, 0.0)
            sum_cni_eff = np.sum(c_n_eff)
            mean_n_i[j] = (np.dot(ns[1:], c_n_eff) / sum_cni_eff) \
                          if sum_cni_eff > 0.0 else 0.0
            c_v_eff = np.maximum(c_v[1:] - C_floor, 0.0)
            sum_cvi_eff = np.sum(c_v_eff)
            mean_n_v[j] = (np.dot(ms[1:], c_v_eff) / sum_cvi_eff) \
                          if sum_cvi_eff > 0.0 else 0.0
            N_loops[j] = np.sum(c_n[1:])
            N_voids[j] = np.sum(c_v[1:])

        # Void swelling (Eq. 161): S = Σ_m m·c_m
        # V_m = m·Ω per cluster; c_m is cluster number density in at.frac;
        # ΔV/V = Σ_m (m·c_m·Ω)/(N_atoms·Ω) = Σ_m m·c_m = C_VAC_tot
        swelling[j] = C_VAC_tot[j]

        C_i1[j]      = c_n[0]
        C_v1[j]      = c_v[0]
        C_He_free[j] = c_h


    # ── Loop-conversion: split SIA loops into ½⟨111⟩ and ⟨100⟩ populations ────
    # y_sia100 [I, n_t] is the appended sessile ⟨100⟩ block (discrete); the
    # ½⟨111⟩ block is y[0:I] (discrete when conversion is on).  We expose, per
    # population: number density (N_loops_111/100), mean size (mean_n_111/100),
    # and the content-weighted ½⟨111⟩ fraction f₁₁₁(t) — directly comparable to
    # the experimental loop fraction.  The combined SIA totals (C_SIA_tot,
    # N_loops, mean_n_i) are also updated so the conservation diagnostics and
    # legacy plots see the full SIA content.  When conversion is off the ⟨100⟩
    # arrays are zero and f₁₁₁ ≡ 1 (all SIA loops are ½⟨111⟩).
    N_loops_111 = N_loops.copy()        # ½⟨111⟩ density (loops, n ≥ 2)
    N_loops_100 = np.zeros(n_t)
    mean_n_111  = mean_n_i.copy()        # ½⟨111⟩ mean size
    mean_n_100  = np.zeros(n_t)
    f_111_loop  = np.ones(n_t)           # ½⟨111⟩ loop fraction (content-weighted)
    if y_sia100 is not None and np.asarray(y_sia100).size:
        # y_sia100 is per-size [I, n_t] in BOTH discrete and bin_moment modes
        # (the bridge reconstructs the bin-moment ⟨100⟩ block before this point).
        # The ½⟨111⟩ content is C_SIA_tot as accumulated above (BEFORE the ⟨100⟩
        # contribution is folded in) — exact from the tracked μ₁ moments in
        # bin-moment mode, and Σ n·c_n in discrete mode — and the ½⟨111⟩ loop
        # density is N_loops_111.  Using these makes the split mode-agnostic.
        cont111 = C_SIA_tot.copy()                           # ½⟨111⟩ content
        cnt111  = N_loops_111                                # ½⟨111⟩ density
        c100    = np.maximum(np.asarray(y_sia100, dtype=float), 0.0)   # [I, n_t]
        cont100 = ns @ c100                                  # ⟨100⟩ content [at.frac]
        N_loops_100 = np.sum(c100[1:, :], axis=0)            # ⟨100⟩ density
        C_SIA_tot   = C_SIA_tot + cont100                    # full SIA inventory
        cnt100  = N_loops_100
        mean_n_100 = np.divide(cont100, cnt100,
                               out=np.zeros(n_t), where=cnt100 > 0)
        tot = cont111 + cont100
        f_111_loop = np.divide(cont111, tot,
                               out=np.ones(n_t), where=tot > 0)
        cnt_tot  = cnt111 + cnt100
        mean_n_i = np.divide(cont111 + cont100, cnt_tot,
                             out=mean_n_i.copy(), where=cnt_tot > 0)
        N_loops = N_loops_111 + N_loops_100                  # combined density

    # Dose axis [dpa]
    dose = G * t

    # ── Conservation diagnostics — read cumulative sinks from state vector ──
    # The cumulative sink fluxes are CVODE-integrated extra state variables.
    # Their row indices are exposed explicitly by the rate-equation object;
    # fall back to the historical "last 5" offsets only if an attribute is
    # missing.
    #   i_J_SIA_fixed  -> J_SIA_fixed:  cumulative SIA content to fixed sinks
    #   i_J_SIA_mutual -> J_SIA_mutual: cumulative SIA to recombination + cavity
    #   i_J_VAC_fixed  -> J_VAC_fixed:  cumulative VAC content to fixed sinks
    #   i_J_VAC_mutual -> J_VAC_mutual: cumulative VAC to recombination + cavity
    #   i_J_He_sink    -> J_He_sink:    cumulative He to sinks
    N_eq = y.shape[0]
    i_J_SIA_fixed  = getattr(rate_eq_obj, 'i_J_SIA_fixed',  N_eq - 5)
    i_J_SIA_mutual = getattr(rate_eq_obj, 'i_J_SIA_mutual', N_eq - 4)
    i_J_VAC_fixed  = getattr(rate_eq_obj, 'i_J_VAC_fixed',  N_eq - 3)
    i_J_VAC_mutual = getattr(rate_eq_obj, 'i_J_VAC_mutual', N_eq - 2)
    i_J_He_sink    = getattr(rate_eq_obj, 'i_J_He_sink',    N_eq - 1)
    J_SIA_fixed  = np.maximum(y[i_J_SIA_fixed,  :], 0.0)
    J_SIA_mutual = np.maximum(y[i_J_SIA_mutual, :], 0.0)
    J_VAC_fixed  = np.maximum(y[i_J_VAC_fixed,  :], 0.0)
    J_VAC_mutual = np.maximum(y[i_J_VAC_mutual, :], 0.0)
    J_He_sink    = np.maximum(y[i_J_He_sink,    :], 0.0)

    # FP conservation — swelling-identity diagnostic (Eq. 96):
    #   δ_FP = |S - S_I - J_SIA_fixed + J_VAC_fixed|
    #          / (S + S_I + J_SIA_fixed + J_VAC_fixed)
    # with S = C_VAC_tot (vacancy inventory) and S_I = C_SIA_tot (SIA
    # inventory). J_SIA_fixed enters with a MINUS sign, J_VAC_fixed with a
    # PLUS sign. J_*_mutual does NOT appear: mutual annihilation destroys
    # equal SIA and vacancy atoms, so it cancels exactly in d(S - S_I)/dt.
    # delta_FP_sia / delta_FP_vac retain the per-species residuals for
    # diagnostic inspection.
    for j in range(n_t):
        num = abs(C_VAC_tot[j] - C_SIA_tot[j]
                  - J_SIA_fixed[j] + J_VAC_fixed[j])
        den = (C_VAC_tot[j] + C_SIA_tot[j]
               + J_SIA_fixed[j] + J_VAC_fixed[j])
        delta_FP[j] = num / den if den > 1e-300 else 0.0

        # Production referenced to the first sample t[0] (see δ_He note below):
        # defects and the cumulative-flux integrals both start from 0 at t[0],
        # so the FP production balance must use ∫_{t0}^t = η·G·(t−t0).
        prod = eta * G * (t[j] - t[0])
        if prod > 1e-300:
            sia0 = C_SIA_tot[0] + J_SIA_fixed[0] + J_SIA_mutual[0]
            vac0 = C_VAC_tot[0] + J_VAC_fixed[0] + J_VAC_mutual[0]
            delta_FP_sia[j] = abs(
                prod - (C_SIA_tot[j] + J_SIA_fixed[j] + J_SIA_mutual[j] - sia0)) / prod
            delta_FP_vac[j] = abs(
                prod - (C_VAC_tot[j] + J_VAC_fixed[j] + J_VAC_mutual[j] - vac0)) / prod
        else:
            delta_FP_sia[j] = delta_FP_vac[j] = 0.0

    # He conservation diagnostic (Eq. 97):
    #   δ_He = |ΔS_He + ΔJ_He_sink - ∫G_He dt'| / (S_He + J_He_sink + ∫G_He dt')
    # The fixed-sink loss J_He_sink must be tracked explicitly; S_He is NOT
    # conserved against ∫G_He dt alone.
    #
    # The integrated balance depends on the He kinetics mode:
    #   dynamic:  c_h is a state variable, so S_He = c_h + Q (free + trapped)
    #             satisfies d(S_He)/dt = G_He - (sink losses) exactly.
    #   QSS:      c_h is algebraic (dc_h/dt = 0 substituted), so only the
    #             trapped inventory obeys dQ/dt = G_He - k2_He*c_h - sink,
    #             and the conserved combination is S_He = Q alone.  Including
    #             the algebraic c_h would book the entire instantaneous
    #             free-He concentration as an apparent conservation error.
    #
    # Production and the accounting integrals are referenced to the FIRST
    # sample t[0], not absolute t=0: CVODE applies the IC at t[0]=t_span[0]
    # (e.g. 1e-6 s) and the cumulative-flux state variables start from 0
    # there.  Using ∫G_He = G_He·t[j] (from 0) while the inventory only
    # accumulates from t[0] leaves a constant G_He·t[0] offset that reads as
    # a spurious ~1e-2 early-time imbalance decaying as t[0]/t.  Differencing
    # against the t[0] baseline removes it and is exact for a constant G_He.
    He_bal = (C_He_tot - C_He_free) if qss_He else C_He_tot
    S_He0 = He_bal[0] + J_He_sink[0]   # conserved combination at t[0]
    for j in range(n_t):
        he_prod = G_He * (t[j] - t[0])   # ∫_{t0}^t G_He dt' (constant G_He)
        num = abs(He_bal[j] + J_He_sink[j] - he_prod - S_He0)
        den = He_bal[j] + J_He_sink[j] + he_prod
        delta_He[j] = num / den if den > 1e-300 else 0.0

    # Convert concentrations from atomic fraction to m^-3
    # c [at.frac] / Omega [m^3/atom] = N [m^-3]
    # swelling, mean sizes, and conservation diagnostics remain dimensionless.
    C_SIA_tot  *= inv_Omega
    C_VAC_tot  *= inv_Omega
    C_He_tot   *= inv_Omega
    C_He_free  *= inv_Omega
    N_loops    *= inv_Omega
    N_loops_111 *= inv_Omega
    N_loops_100 *= inv_Omega
    N_voids    *= inv_Omega
    C_i1       *= inv_Omega
    C_v1       *= inv_Omega

    # n_active_sia / n_active_vac / n_active are populated by cpp_bridge from
    # the C++ window-bounds sidecar (<bin>.window.csv).  We can't infer them
    # from y because out-of-window indices remain at C_floor (the initial
    # value), so they look identical to evolved-but-still-low in-window slots.

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
        'N_loops':    N_loops,     # [m^-3]  (combined SIA loops)
        'N_loops_111': N_loops_111, # [m^-3]  ½⟨111⟩ loops
        'N_loops_100': N_loops_100, # [m^-3]  ⟨100⟩ loops (loop conversion)
        'mean_n_111':  mean_n_111,  # [atoms] ½⟨111⟩ mean loop size
        'mean_n_100':  mean_n_100,  # [atoms] ⟨100⟩ mean loop size
        'f_111_loop':  f_111_loop,  # ½⟨111⟩ loop fraction (content-weighted)
        'N_voids':    N_voids,     # [m^-3]
        'swelling':   swelling,    # [dimensionless fraction]
        'C_i1':       C_i1,        # [m^-3]
        'C_v1':       C_v1,        # [m^-3]
        'delta_FP':     delta_FP,     # FP conservation: max(sia, vac) relative error
        'delta_FP_sia': delta_FP_sia, # SIA arm of FP conservation
        'delta_FP_vac': delta_FP_vac, # VAC arm of FP conservation
        'delta_He':   delta_He,    # He conservation (Eq. 97; dynamic + QSS)
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


def _last_finite(arr):
    """Return the last finite (non-NaN, non-inf) value of an array.

    Falls back to an empty string when the array holds no finite value,
    so the CSV never carries a literal 'nan'/'inf'.
    """
    a = np.asarray(arr, dtype=float)
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return ''
    return finite[-1]


def _safe_scalar(value):
    """Guard a single scalar against NaN/inf for CSV output."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return value
    return v if np.isfinite(v) else ''


def summary_csv_row(results, input_data, solver_label=''):
    """
    Build a summary dict for a single run (for output/summary.csv).

    Every numeric value is guarded against NaN/inf: array metrics emit the
    last finite element (or empty string if none), scalar metrics emit an
    empty string when non-finite.  This prevents literal 'nan' in
    summary.csv.
    """
    d = input_data.derived
    t_end = results['t'][-1]
    row = {
        'solver':        solver_label,
        'physics_option':input_data.physics_option,
        'T_K':           _safe_scalar(d['T']),
        'G_dpa_s':       _safe_scalar(d['G']),
        't_end_s':       _safe_scalar(t_end),
        'dose_dpa':      _last_finite(results['dose']),
        'C_SIA_tot_m3':  _last_finite(results['C_SIA_tot']),
        'C_VAC_tot_m3':  _last_finite(results['C_VAC_tot']),
        'C_He_tot_m3':   _last_finite(results['C_He_tot']),
        'mean_n_i':      _last_finite(results['mean_n_i']),
        'mean_n_v':      _last_finite(results['mean_n_v']),
        'swelling':      _last_finite(results['swelling']),
        'delta_FP':      _last_finite(results['delta_FP']),
        'delta_He':      _last_finite(results['delta_He']),
    }
    return row
