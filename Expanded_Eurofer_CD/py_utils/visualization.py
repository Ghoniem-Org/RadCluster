"""
visualization.py — Plotting routines for Expanded_Eurofer_CD.

Generates standardized figures from the ODE post-processing results.
All concentration quantities are in m^-3 (converted in post_process.py).
"""

import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

_CONC_LABEL = r'Concentration (m$^{-3}$)'


def _check_mpl():
    if not _HAS_MPL:
        raise ImportError("matplotlib is required for visualization.")


def plot_point_defects(results, out_path=None, title=''):
    """Free SIA (Ci1) and vacancy (Cv1) monomers vs. dose."""
    _check_mpl()
    fig, ax = plt.subplots(figsize=(7, 4))
    dose = results['dose']
    ax.loglog(dose, results['C_i1'], label='C_i1 (free SIA)', color='royalblue')
    ax.loglog(dose, results['C_v1'], label='C_v1 (free vac)', color='firebrick')
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel(_CONC_LABEL)
    ax.set_title(f'Point Defects {title}')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(left=1e-4)
    ax.set_ylim(bottom=1e8)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_totals(results, out_path=None, title=''):
    """Total SIA, vacancy, and He contents vs. dose."""
    _check_mpl()
    fig, ax = plt.subplots(figsize=(7, 4))
    dose = results['dose']
    ax.loglog(dose, results['C_SIA_tot'], label=r'SIA total ($\Sigma n\cdot c_n$)', color='steelblue')
    ax.loglog(dose, results['C_VAC_tot'], label=r'Vac total ($\Sigma m\cdot c_m$)', color='tomato')
    ax.loglog(dose, results['C_He_tot'],  label='He total',                          color='green')
    ax.set_xlim(left=1e-6)
    ax.set_ylim(bottom=1e8)
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel(_CONC_LABEL)
    ax.set_title(f'Defect Contents {title}')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_swelling(results, out_path=None, title=''):
    """Void swelling S(t) vs. dose."""
    _check_mpl()
    fig, ax = plt.subplots(figsize=(7, 4))
    dose = results['dose']
    ax.semilogx(dose, results['swelling'] * 100, color='darkorange', lw=2)
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel('Swelling S(t) (%)')
    ax.set_title(f'Void Swelling (Eq. 161) {title}')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_mean_sizes(results, out_path=None, title=''):
    """Mean SIA and vacancy cluster sizes vs. dose."""
    _check_mpl()
    fig, ax = plt.subplots(figsize=(7, 4))
    dose = results['dose']
    ax.semilogx(dose, results['mean_n_i'], label=r'$\langle n \rangle$ SIA loops', color='steelblue')
    ax.semilogx(dose, results['mean_n_v'], label=r'$\langle m \rangle$ voids',     color='tomato')
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel('Mean cluster size (defects)')
    ax.set_title(f'Mean Cluster Sizes {title}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=1e-6)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_he_content(results, out_path=None, title=''):
    """Free and total He vs. dose."""
    _check_mpl()
    fig, ax = plt.subplots(figsize=(7, 4))
    dose = results['dose']
    ax.loglog(dose, results['C_He_tot'],  label='He total',       color='darkgreen')
    ax.loglog(dose, results['C_He_free'], label=r'He free ($c_h$)', color='limegreen', ls='--')
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel(_CONC_LABEL)
    ax.set_title(f'Helium Content {title}')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(left=1e-6)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_conservation(results, out_path=None, title=''):
    """
    FP balance and He retention diagnostics (Eqs. 164-165).

    δ_FP = |Σn·n·c_n − Σm·m·c_m| / (η·G·t)  — SIA/VAC imbalance; near 0 = balanced.
    δ_He = |c_h + Q_tot − G_He·t| / (G_He·t) — He retention; near 0 = all He in clusters/free.
    """
    _check_mpl()
    fig, ax = plt.subplots(figsize=(7, 4))
    dose = results['dose']
    ax.semilogy(dose, np.maximum(results['delta_FP'], 1e-20),
                label=r'$\delta_{FP}$ (SIA–VAC imbalance, Eq. 164)', color='navy')
    ax.semilogy(dose, np.maximum(results['delta_He'], 1e-20),
                label=r'$\delta_{He}$ (He retention, Eq. 165)', color='darkred')
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel('Diagnostic (near 0 = good)')
    ax.set_title(f'Balance Diagnostics {title}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=1e-6)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_size_distribution(results, input_data, t_idx=-1, out_path=None, title=''):
    """SIA and vacancy cluster size distributions at a given time index."""
    _check_mpl()
    N = input_data.N
    M = input_data.M
    Omega = results.get('Omega', input_data.derived['Omega'])
    inv_Omega = 1.0 / Omega

    y = results['y'][:, t_idx]
    if y.shape[0] < N + M:
        return None

    # Raw ODE state is in at.frac; convert to m^-3 for display
    c_n = np.maximum(y[:N], 0.0) * inv_Omega
    c_v = np.maximum(y[N:N + M], 0.0) * inv_Omega

    ns = np.arange(1, N + 1)
    ms = np.arange(1, M + 1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    ax1.semilogy(ns, c_n + 1e-10, color='steelblue')
    ax1.set_xlabel('SIA cluster size n')
    ax1.set_ylabel(_CONC_LABEL)
    ax1.set_title(f'SIA Distribution  t={results["t"][t_idx]:.2e} s {title}')
    ax1.grid(True, alpha=0.3)

    ax2.semilogy(ms, c_v + 1e-10, color='tomato')
    ax2.set_xlabel('Vacancy cluster size m')
    ax2.set_ylabel(_CONC_LABEL)
    ax2.set_title(f'Void Distribution  t={results["t"][t_idx]:.2e} s {title}')
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def _extract_cluster_arrays(results, input_data):
    """
    Return (c_n, c_v, dose) arrays [size x n_t] in m^-3.

    Works only for full_CD modes (y has per-cluster rows).
    Returns None for bin-moment runs where y rows are moments.
    """
    N     = input_data.N
    M     = input_data.M
    y     = results['y']          # [N_eq, n_t], at.frac
    Omega = results.get('Omega', input_data.derived['Omega'])
    inv_O = 1.0 / Omega

    # Sanity check: full_CD has at least N+M rows
    if y.shape[0] < N + M:
        return None, None, None

    c_n = np.maximum(y[:N, :],     0.0) * inv_O   # [N, n_t] m^-3
    c_v = np.maximum(y[N:N + M, :], 0.0) * inv_O  # [M, n_t] m^-3
    return c_n, c_v, results['dose']


def plot_sia_clusters(results, input_data, out_path=None, title=''):
    """SIA cluster concentrations vs. dose, split into small / mid / large."""
    _check_mpl()
    c_n, c_v, dose = _extract_cluster_arrays(results, input_data)
    if c_n is None:
        print("plot_sia_clusters: skipped (bin-moment mode).")
        return None

    N = input_data.N
    groups = [
        ('SIA clusters: n = 1–5',      range(1, min(6,   N + 1)),  'sia_small'),
        ('SIA clusters: n = 6–20',     range(6, min(21,  N + 1)),  'sia_mid'),
        ('SIA clusters: n = 20–100',   range(20, min(101, N + 1)), 'sia_large'),
    ]
    figs = []
    for gtitle, rng, fname in groups:
        idx = [n - 1 for n in rng if n - 1 < N]
        if not idx:
            continue
        fig, ax = plt.subplots(figsize=(7, 4))
        cmap = plt.cm.viridis(np.linspace(0, 1, len(idx)))
        for k, color in zip(idx, cmap):
            ax.loglog(dose, np.maximum(c_n[k, :], 1e-10),
                      color=color, label=f'n={k + 1}')
        ax.set_xlabel('Dose (dpa)')
        ax.set_ylabel(_CONC_LABEL)
        ax.set_title(f'{gtitle} {title}')
        ax.grid(True, which='both', alpha=0.3)
        ax.set_xlim(left=1e-6)
        if len(idx) <= 10:
            ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        if out_path:
            # Replace file stem with group-specific name
            from pathlib import Path as _P
            p = _P(out_path)
            fig.savefig(str(p.parent / f'{fname}.png'), dpi=150)
        figs.append(fig)
    return figs


def plot_vac_clusters(results, input_data, out_path=None, title=''):
    """Vacancy cluster concentrations vs. dose, split into small / mid / large."""
    _check_mpl()
    c_n, c_v, dose = _extract_cluster_arrays(results, input_data)
    if c_v is None:
        print("plot_vac_clusters: skipped (bin-moment mode).")
        return None

    M = input_data.M
    groups = [
        ('Vacancy clusters: m = 1–5',      range(1, min(6,   M + 1)),  'vac_small'),
        ('Vacancy clusters: m = 6–20',     range(6, min(21,  M + 1)),  'vac_mid'),
        ('Vacancy clusters: m = 20–100',   range(20, min(101, M + 1)), 'vac_large'),
    ]
    figs = []
    for gtitle, rng, fname in groups:
        idx = [m - 1 for m in rng if m - 1 < M]
        if not idx:
            continue
        fig, ax = plt.subplots(figsize=(7, 4))
        cmap = plt.cm.plasma(np.linspace(0, 1, len(idx)))
        for k, color in zip(idx, cmap):
            ax.loglog(dose, np.maximum(c_v[k, :], 1e-10),
                      color=color, label=f'm={k + 1}')
        ax.set_xlabel('Dose (dpa)')
        ax.set_ylabel(_CONC_LABEL)
        ax.set_title(f'{gtitle} {title}')
        ax.grid(True, which='both', alpha=0.3)
        ax.set_xlim(left=1e-6)
        if len(idx) <= 10:
            ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        if out_path:
            from pathlib import Path as _P
            p = _P(out_path)
            fig.savefig(str(p.parent / f'{fname}.png'), dpi=150)
        figs.append(fig)
    return figs


def _log_snapshot_indices(dose, n_times=10):
    """Return n_times indices log-spaced across the dose array."""
    d = np.asarray(dose)
    d_pos = d[d > 0]
    if len(d_pos) < 2:
        return np.linspace(0, len(d) - 1, min(n_times, len(d)), dtype=int)
    log_min = np.log10(d_pos[0])
    log_max = np.log10(d[-1])
    targets = np.logspace(log_min, log_max, n_times)
    return [int(np.argmin(np.abs(d - t))) for t in targets]


def plot_sia_distribution_evolution(results, input_data, n_times=10,
                                    out_path=None, title=''):
    """
    SIA cluster size distribution c_n(n) at n_times log-spaced dose snapshots.

    Each curve is coloured from light to dark (viridis) and labelled by dose.
    """
    _check_mpl()
    N     = input_data.N
    Omega = results.get('Omega', input_data.derived['Omega'])
    y     = results['y']           # [N_eq, n_t] at.frac
    dose  = results['dose']

    if y.shape[0] < N:
        print("plot_sia_distribution_evolution: skipped (bin-moment mode).")
        return None

    inv_O   = 1.0 / Omega
    ns      = np.arange(1, N + 1)
    indices = _log_snapshot_indices(dose, n_times)
    cmap    = plt.cm.viridis(np.linspace(0.15, 0.95, len(indices)))

    fig, ax = plt.subplots(figsize=(8, 5))
    for idx, color in zip(indices, cmap):
        c_n = np.maximum(y[:N, idx], 0.0) * inv_O
        ax.semilogy(ns, c_n + 1e-10, color=color,
                    label=f'{dose[idx]:.2e} dpa')

    ax.set_xlabel('SIA cluster size n')
    ax.set_ylabel(_CONC_LABEL)
    ax.set_title(f'SIA Size Distribution Evolution {title}')
    ax.legend(fontsize=7, ncol=2, title='Dose')
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_void_distribution_evolution(results, input_data, n_times=10,
                                     out_path=None, title=''):
    """
    Void/vacancy cluster size distribution c_m(m) at n_times log-spaced dose snapshots.

    Each curve is coloured from light to dark (plasma) and labelled by dose.
    """
    _check_mpl()
    N     = input_data.N
    M     = input_data.M
    Omega = results.get('Omega', input_data.derived['Omega'])
    y     = results['y']           # [N_eq, n_t] at.frac
    dose  = results['dose']

    if y.shape[0] < N + M:
        print("plot_void_distribution_evolution: skipped (bin-moment mode).")
        return None

    inv_O   = 1.0 / Omega
    ms      = np.arange(1, M + 1)
    indices = _log_snapshot_indices(dose, n_times)
    cmap    = plt.cm.plasma(np.linspace(0.15, 0.85, len(indices)))

    fig, ax = plt.subplots(figsize=(8, 5))
    for idx, color in zip(indices, cmap):
        c_v = np.maximum(y[N:N + M, idx], 0.0) * inv_O
        ax.semilogy(ms, c_v + 1e-10, color=color,
                    label=f'{dose[idx]:.2e} dpa')

    ax.set_xlabel('Vacancy cluster size m')
    ax.set_ylabel(_CONC_LABEL)
    ax.set_title(f'Void Size Distribution Evolution {title}')
    ax.legend(fontsize=7, ncol=2, title='Dose')
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_number_densities(results, out_path=None, title=''):
    """Loop and void number densities (clusters n,m ≥ 2) vs. dose."""
    _check_mpl()
    fig, ax = plt.subplots(figsize=(7, 4))
    dose = results['dose']
    ax.loglog(dose, np.maximum(results['N_loops'], 1e-10),
              label='Loop density (SIA, n≥2)', color='steelblue')
    ax.loglog(dose, np.maximum(results['N_voids'], 1e-10),
              label='Void density (vac, m≥2)', color='tomato', ls='--')
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel(r'Number density (m$^{-3}$)')
    ax.set_title(f'Loop and Void Number Densities {title}')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(left=1e-6)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def save_all_plots(results, input_data, out_dir, label=''):
    """Save all standard plots to out_dir/."""
    from pathlib import Path
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    opts = dict(title=label)
    plot_point_defects(results,      out_path=f"{out_dir}/point_defects.png",    **opts)
    plot_totals(results,             out_path=f"{out_dir}/totals.png",            **opts)
    plot_swelling(results,           out_path=f"{out_dir}/swelling.png",          **opts)
    plot_mean_sizes(results,         out_path=f"{out_dir}/mean_sizes.png",        **opts)
    plot_he_content(results,         out_path=f"{out_dir}/he_content.png",        **opts)
    plot_conservation(results,       out_path=f"{out_dir}/conservation.png",      **opts)
    plot_number_densities(results,   out_path=f"{out_dir}/number_densities.png",  **opts)
    plot_size_distribution(results, input_data, out_path=f"{out_dir}/size_distributions.png", **opts)
    plot_sia_distribution_evolution(results, input_data, out_path=f"{out_dir}/sia_dist_evolution.png", **opts)
    plot_void_distribution_evolution(results, input_data, out_path=f"{out_dir}/void_dist_evolution.png", **opts)
    # Per-cluster time-series (full_CD modes only; skipped for bin_moment)
    plot_sia_clusters(results, input_data, out_path=f"{out_dir}/sia_small.png",  **opts)
    plot_vac_clusters(results, input_data, out_path=f"{out_dir}/vac_small.png",  **opts)
    print(f"Saved plots to {out_dir}/")
