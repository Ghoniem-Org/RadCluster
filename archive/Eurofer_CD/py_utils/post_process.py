"""
post_process.py — Derived quantities and solution quality checks for Eurofer_CD.

Called after the ODE integration completes.

Quantities computed
-------------------
- Total SIA content (all cluster sizes weighted by n)
- Total vacancy content (all cluster sizes weighted by m)
- Total He content (free + trapped estimate)
- Mean SIA and vacancy cluster sizes
- Void swelling estimate: S = (4π/3) · r_m^3 · N_v (volumetric)
- Active interstitial band (cluster sizes above threshold)

Differences from Full_CD/post_process.py
-----------------------------------------
- State vector has He species at index N-1: y[Ni+Nv]
- 'Nv' dimension refers to vacancy clusters, 'Ni' to SIA clusters
- Swelling calculation uses r_m = r_0 · m^(1/3)
"""

import numpy as np

_kB = 8.617333262e-5   # eV K^-1


# ── Solution quality ──────────────────────────────────────────────────────────

def check_solution_quality(t, y, concentration_names):
    """
    Inspect the ODE solution for numerical issues.

    Parameters
    ----------
    t                   : ndarray [n_time]
    y                   : ndarray [N, n_time]
    concentration_names : list[str]

    Returns
    -------
    warnings_found : list[str]
    """
    warnings_found = []

    min_vals = np.min(y, axis=1)
    neg_idx  = np.where(min_vals < -1e-12)[0]
    if len(neg_idx):
        neg_names = [concentration_names[i] for i in neg_idx[:5]]
        warnings_found.append(f"Significant negative concentrations: {neg_names}")

    if y.shape[1] > 10:
        early = np.mean(y[0, :5])
        late  = np.mean(y[0, -5:])
        if early > 1e-30 and abs(late - early) / early > 20.0:
            warnings_found.append(
                f"Ci1 changed by > {100*abs(late-early)/early:.0f}% over run"
            )

    if warnings_found:
        print("Solution quality warnings:")
        for w in warnings_found:
            print(f"  - {w}")
    else:
        print("Solution quality: Good")

    return warnings_found


# ── Derived quantities ─────────────────────────────────────────────────────────

def calculate_derived_quantities(time, concentrations, input_data, rate_equations,
                                 xmax_history=None):
    """
    Convert raw ODE output into the standard results dict.

    Parameters
    ----------
    time             : ndarray [n_time]
    concentrations   : ndarray [N, n_time]
    input_data       : InputData
    rate_equations   : RateEquations
    xmax_history     : list, optional

    Returns
    -------
    dict with keys: time, concentrations, totals, mean_sizes, swelling,
                    active_band, He_content, xmax_history, Nv, Ni
    """
    print("Calculating derived quantities…")

    re   = rate_equations
    Ni   = re.Ni
    Nv   = re.Nv
    i0   = re.i_SIA   # = 0
    iv   = re.i_VAC   # = Ni
    iHe  = re.i_He    # = Ni + Nv
    names = re.concentration_names

    conc = np.maximum(concentrations, 0.0)

    # ── SIA cluster totals ───────────────────────────────────────────────────
    ns = np.arange(1, Ni + 1, dtype=float)[:, None]    # (Ni, 1)
    Ci_arr  = conc[i0:iv, :]                            # (Ni, n_time)
    total_i = np.sum(ns * Ci_arr, axis=0)               # weighted sum

    # ── Vacancy cluster totals ───────────────────────────────────────────────
    ms = np.arange(1, Nv + 1, dtype=float)[:, None]    # (Nv, 1)
    Cv_arr  = conc[iv:iHe, :]                           # (Nv, n_time)
    total_v = np.sum(ms * Cv_arr, axis=0)

    # ── He content ───────────────────────────────────────────────────────────
    C_He_arr = conc[iHe, :]                             # (n_time,)

    # ── Mean cluster sizes (clusters only: n≥2, m≥2) ─────────────────────────
    # Exclude monomers: Ci1 >> Ci_n≥2, so including n=1 pins the mean at ≈1.
    sum_nCi = np.sum(ns[1:] * Ci_arr[1:], axis=0)
    sum_Ci  = np.sum(Ci_arr[1:], axis=0)
    mean_n  = np.where(sum_Ci > 1e-30, sum_nCi / sum_Ci, 1.0)

    sum_mCv = np.sum(ms[1:] * Cv_arr[1:], axis=0)
    sum_Cv  = np.sum(Cv_arr[1:], axis=0)
    mean_m  = np.where(sum_Cv > 1e-30, sum_mCv / sum_Cv, 1.0)

    # ── Void swelling ─────────────────────────────────────────────────────────
    # S = (4π/3) · Σ_m r_m^3 · C_v(m)  (volume fraction, dimensionless)
    # r_m = r_0 · m^(1/3); r_0 = (3Ω/4π)^(1/3)
    Omega = float(input_data.derived['Omega'])
    r0    = float(input_data.derived['r0'])
    r_m   = r0 * np.arange(1, Nv + 1, dtype=float)**(1.0/3.0)   # (Nv,)
    vol_m = (4.0 * np.pi / 3.0) * r_m**3                         # m^3 per cluster

    # Convert from atom fractions to number density [m^-3]:
    # N_m [m^-3] = C_v(m) / Omega
    # V_fraction = Σ_m N_m · vol_m = Σ_m C_v(m) * vol_m / Omega
    swelling = np.sum(
        (vol_m[:, None] / Omega) * Cv_arr, axis=0
    )  # (n_time,)

    # ── Active interstitial band (sizes above threshold) ─────────────────────
    THRESH = 1e-18
    x_max_plot = []
    x_min_plot = []
    for j in range(conc.shape[1]):
        appeared = np.where(Ci_arr[:, j] > THRESH)[0]
        x_max_plot.append(int(appeared[-1] + 1) if len(appeared) else 0)
        x_min_plot.append(int(appeared[0]  + 1) if len(appeared) else 1)

    print("Derived quantities calculated.")

    return {
        'time':           time,
        'concentrations': {
            name: concentrations[idx, :]
            for idx, name in enumerate(names)
        },
        'totals': {
            'total_i':  total_i,
            'total_v':  total_v,
            'C_He':     C_He_arr,
        },
        'mean_sizes': {
            'mean_n':   mean_n,    # mean SIA cluster size
            'mean_m':   mean_m,    # mean vacancy cluster size
        },
        'swelling':      swelling,
        'active_band': {
            'x_max':    np.array(x_max_plot),
            'x_min':    np.array(x_min_plot),
        },
        'xmax_history':  xmax_history if xmax_history is not None else [],
        'Ni':            Ni,
        'Nv':            Nv,
    }
