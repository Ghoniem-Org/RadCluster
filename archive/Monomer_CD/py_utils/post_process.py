"""
post_process.py – Derived quantities and solution quality checks.

Called after the ODE integration completes.  All physics transformations
of the raw concentrations into summary statistics live here.

This module is intentionally decoupled from the solver so that the C++
back-end reuses the same post-processing pipeline.
"""

import numpy as np


# ── Solution quality ──────────────────────────────────────────────────────────

def check_solution_quality(t, y, concentration_names):
    """
    Inspect the ODE solution for numerical issues.

    Parameters
    ----------
    t                   : numpy.ndarray  [n_time]
    y                   : numpy.ndarray  [N, n_time]
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

    # Large overall change in Cv1 + Ci1
    if y.shape[1] > 10:
        Nv = sum(1 for n in concentration_names if n.startswith('Cv'))
        early = np.mean(y[0, :5]  + y[Nv, :5])
        late  = np.mean(y[0, -5:] + y[Nv, -5:])
        if early > 1e-30 and abs(late - early) / early > 10.0:
            warnings_found.append(
                f"Cv/Ci changed by >{100*(abs(late-early)/early):.0f}% over run"
            )

    if warnings_found:
        print("⚠️  Solution quality warnings:")
        for w in warnings_found:
            print(f"    - {w}")
    else:
        print("✓ Solution quality: Good")

    return warnings_found


# ── Derived quantities ─────────────────────────────────────────────────────────

def calculate_derived_quantities(time, concentrations, input_data, rate_equations,
                                 xmax_history=None):
    """
    Convert raw ODE output into the standard results dict.

    Parameters
    ----------
    time             : numpy.ndarray  [n_time]
    concentrations   : numpy.ndarray  [N, n_time]
    input_data       : InputData
    rate_equations   : RateEquations
    xmax_history     : list, optional – x_max per segment (Python solver only)

    Returns
    -------
    dict with keys: time, concentrations, totals, mean_sizes, xmax_history
    """
    print("Calculating derived quantities (vectorized)…")

    Nv = rate_equations.Nv
    Ni = rate_equations.Ni
    N  = rate_equations.N
    names = rate_equations.concentration_names

    conc = np.maximum(concentrations, 0.0)

    # ── Total vacancy / interstitial atom content (weighted by cluster size) ─
    xs_v = np.arange(1, Nv + 1, dtype=float)[:, None]   # (Nv, 1)
    xs_i = np.arange(1, Ni + 1, dtype=float)[:, None]   # (Ni, 1)

    total_v = np.sum(xs_v * conc[:Nv, :], axis=0)        # (n_time,)
    total_i = np.sum(xs_i * conc[Nv:, :], axis=0)        # (n_time,)

    # ── Mean cluster sizes ────────────────────────────────────────────────────
    sum_nCv = np.sum(xs_v * conc[:Nv, :], axis=0)
    sum_Cv  = np.sum(conc[:Nv, :], axis=0)
    mean_v  = np.where(sum_Cv > 1e-30, sum_nCv / sum_Cv, 1.0)

    sum_nCi = np.sum(xs_i * conc[Nv:, :], axis=0)
    sum_Ci  = np.sum(conc[Nv:, :], axis=0)
    mean_i  = np.where(sum_Ci > 1e-30, sum_nCi / sum_Ci, 1.0)

    # ── Active interstitial band (sizes above threshold) ─────────────────────
    THRESH = 1e-18
    x_max_plot = []
    x_min_plot = []
    for j in range(conc.shape[1]):
        appeared = np.where(conc[Nv:, j] > THRESH)[0]
        x_max_plot.append(int(appeared[-1] + 1) if len(appeared) else 0)
        x_min_plot.append(int(appeared[0]  + 1) if len(appeared) else 1)

    print("✓ Derived quantities calculated.")

    return {
        'time': time,
        'concentrations': {
            name: concentrations[idx, :]
            for idx, name in enumerate(names)
        },
        'totals': {
            'total_v': total_v,
            'total_i': total_i,
        },
        'mean_sizes': {
            'mean_v': mean_v,
            'mean_i': mean_i,
        },
        'active_band': {
            'x_max': np.array(x_max_plot),
            'x_min': np.array(x_min_plot),
        },
        'xmax_history': xmax_history if xmax_history is not None else [],
        'Nv': Nv,
        'Ni': Ni,
    }
