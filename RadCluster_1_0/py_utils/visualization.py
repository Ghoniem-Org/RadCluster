"""
visualization.py — Plotting routines for RadCluster_1_0.

Generates standardized figures from the ODE post-processing results.
All concentration quantities are in m^-3 (converted in post_process.py).
"""

import numpy as np

try:
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

_CONC_LABEL = r'Concentration (m$^{-3}$)'

_PLOT_FONTSIZE = 16
if _HAS_MPL:
    plt.rcParams.update({
        'axes.titlesize':  _PLOT_FONTSIZE,
        'axes.labelsize':  _PLOT_FONTSIZE,
        'xtick.labelsize': _PLOT_FONTSIZE,
        'ytick.labelsize': _PLOT_FONTSIZE,
    })


# ── User-configurable axis controls ──────────────────────────────────────────
# Plots are split into three groups so a small set of knobs covers the whole
# suite.  Each group entry accepts:
#   xlim / ylim : (min, max) tuple — None for auto on either bound.
#   xscale / yscale : 'log' | 'linear' | None (None = keep the plot's default).
# A global `legend_fontsize` controls the size of the figure-level legends
# placed underneath cluster-concentration plots.
#
# Override from the notebook with:
#     viz.set_plot_config({'concentration': {'ylim': (1e16, 1e22)},
#                          'legend_fontsize': 5})
_PLOT_CONFIG = {
    # Concentration vs dose: point defects, totals, He content, number
    # densities (incl. TEM), SIA/cavity cluster concentrations.
    'concentration': {'xlim': (None, None), 'ylim': (None, None),
                      'xscale': None,       'yscale': None},
    # Scalar metrics vs dose: swelling, mean sizes, fraction balances,
    # conservation diagnostics.
    'scalar':        {'xlim': (None, None), 'ylim': (None, None),
                      'xscale': None,       'yscale': None},
    # Size-distribution snapshots: concentration vs cluster size n/m or
    # vs loop/void diameter (nm).
    'size_dist':     {'xlim': (None, None), 'ylim': (None, None),
                      'xscale': None,       'yscale': None},
    'legend_fontsize': 5,
}


def set_plot_config(cfg):
    """Update the module-level plot config; partial overrides are merged."""
    for k, v in cfg.items():
        cur = _PLOT_CONFIG.get(k)
        if isinstance(cur, dict) and isinstance(v, dict):
            cur.update(v)
        else:
            _PLOT_CONFIG[k] = v


def _apply_axis_config(ax, group):
    """Apply the named group's scale and limits to a single axis."""
    cfg = _PLOT_CONFIG.get(group)
    if not cfg:
        return
    xs, ys = cfg.get('xscale'), cfg.get('yscale')
    if xs:
        ax.set_xscale(xs)
    if ys:
        ax.set_yscale(ys)
    xlim = cfg.get('xlim') or (None, None)
    ylim = cfg.get('ylim') or (None, None)
    if xlim[0] is not None or xlim[1] is not None:
        ax.set_xlim(left=xlim[0], right=xlim[1])
    if ylim[0] is not None or ylim[1] is not None:
        ax.set_ylim(bottom=ylim[0], top=ylim[1])


def _apply_axis_config_to_fig(fig, group):
    for ax in fig.axes:
        _apply_axis_config(ax, group)


def _legend_fs():
    return _PLOT_CONFIG.get('legend_fontsize', 5)


def _legend_below(fig, axes, ncol=8, fontsize=None, bottom=0.28):
    """
    Place a single shared legend underneath the plot area.

    Removes any existing per-axes legends, then adds a figure-level legend
    using the handles from the first axis (cluster-concentration plots
    duplicate the same labels across panels).
    """
    if not _HAS_MPL:
        return
    if fontsize is None:
        fontsize = _legend_fs()
    if not isinstance(axes, (list, tuple, np.ndarray)):
        axes = [axes]
    for ax in axes:
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()
    handles, labels = axes[0].get_legend_handles_labels()
    if not handles:
        return
    fig.legend(handles, labels, loc='lower center',
               bbox_to_anchor=(0.5, 0.0),
               ncol=ncol, fontsize=fontsize, frameon=True)
    fig.subplots_adjust(bottom=bottom)


def _check_mpl():
    if not _HAS_MPL:
        raise ImportError("matplotlib is required for visualization.")


def _dose_xlim(dose):
    """
    Return the left x-limit for a dose axis.

    Uses the minimum positive dose value in the data so the plot is never
    reversed when the run ends at a dose below any hard-coded threshold.
    """
    pos = dose[dose > 0]
    return float(pos.min()) if len(pos) > 0 else 1e-20


def _align_dose_to_y(results):
    """
    Return (y, dose) with matching time dimensions.

    After adaptive domain-doubling, ``results['y']`` may contain only the
    last segment (its row count changed), while ``results['dose']`` spans
    the full run.  When the two lengths differ, take the *tail* of dose
    to match ``y``.
    """
    y = results['y']
    dose = results['dose']
    n_y = y.shape[1]
    if len(dose) != n_y:
        dose = dose[-n_y:]
    return y, dose


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
    ax.set_xlim(left=_dose_xlim(dose))
    ax.set_ylim(bottom=1e14)
    _apply_axis_config(ax, 'concentration')
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
    ax.set_xlim(left=_dose_xlim(dose))
    ax.set_ylim(bottom=1e11)
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel(_CONC_LABEL)
    ax.set_title(f'Defect Contents {title}')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    _apply_axis_config(ax, 'concentration')
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
    _apply_axis_config(ax, 'scalar')
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
    ax.set_xlim(left=_dose_xlim(dose))
    ax.set_ylim(bottom=0)
    _apply_axis_config(ax, 'scalar')
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
    ax.set_xlim(left=_dose_xlim(dose))
    _apply_axis_config(ax, 'concentration')
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
    ax.set_xlim(left=_dose_xlim(dose))
    ax.set_ylim(bottom=1e-5, top=10)
    _apply_axis_config(ax, 'scalar')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_size_distribution(results, input_data, t_idx=-1, out_path=None, title=''):
    """SIA and vacancy cluster size distributions at a given time index."""
    _check_mpl()
    N = input_data.I
    M = input_data.V
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

    # Diameter conversions
    r0 = results.get('r0', input_data.derived['r0'])         # Wigner-Seitz radius [m]
    b_111 = results.get('b_111', input_data.derived['b_111'])  # Burgers vector [m]
    d_loop = 2.0 * np.sqrt(ns * Omega / (np.pi * b_111)) * 1e9   # SIA loop diameter [nm]
    d_void = 2.0 * r0 * ms ** (1.0 / 3.0) * 1e9                  # void diameter [nm]

    time_label = f't={results["t"][t_idx]:.2e} s'

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # --- Row 1: log-scale vs cluster size ---
    ax1 = axes[0, 0]
    ax1.semilogy(ns, c_n + 1e-10, color='steelblue')
    ax1.set_xlabel('SIA cluster size n')
    ax1.set_ylabel(_CONC_LABEL)
    ax1.set_title(f'SIA Distribution  {time_label} {title}')
    ax1.grid(True, alpha=0.3)
    c_n_pos = c_n[c_n > 0]
    if c_n_pos.size > 0:
        ax1.set_ylim(bottom=c_n_pos.min() * 0.1)

    ax2 = axes[0, 1]
    ax2.semilogy(ms, c_v + 1e-10, color='tomato')
    ax2.set_xlabel('Vacancy cluster size m')
    ax2.set_ylabel(_CONC_LABEL)
    ax2.set_title(f'Void Distribution  {time_label} {title}')
    ax2.grid(True, alpha=0.3)
    c_v_pos = c_v[c_v > 0]
    if c_v_pos.size > 0:
        ax2.set_ylim(bottom=c_v_pos.min() * 0.1)

    # --- Row 2: linear-scale vs diameter (nm) ---
    ax3 = axes[1, 0]
    ax3.plot(d_loop, c_n, color='steelblue')
    ax3.set_xlabel('SIA loop diameter (nm)')
    ax3.set_ylabel(_CONC_LABEL)
    ax3.set_title(f'SIA Loop Size  {time_label} {title}')
    ax3.grid(True, alpha=0.3)

    ax4 = axes[1, 1]
    ax4.plot(d_void, c_v, color='tomato')
    ax4.set_xlabel('Void diameter (nm)')
    ax4.set_ylabel(_CONC_LABEL)
    ax4.set_title(f'Void Size  {time_label} {title}')
    ax4.grid(True, alpha=0.3)

    _apply_axis_config_to_fig(fig, 'size_dist')
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
    N     = input_data.I
    M     = input_data.V
    y, dose = _align_dose_to_y(results)
    Omega = results.get('Omega', input_data.derived['Omega'])
    inv_O = 1.0 / Omega

    # Sanity check: full_CD has at least N+M rows
    if y.shape[0] < N + M:
        return None, None, None

    c_n = np.maximum(y[:N, :],     0.0) * inv_O   # [N, n_t] m^-3
    c_v = np.maximum(y[N:N + M, :], 0.0) * inv_O  # [M, n_t] m^-3
    return c_n, c_v, dose


def plot_sia_clusters(results, input_data, out_path=None, title=''):
    """SIA cluster concentrations vs. dose, split into small / mid / large."""
    _check_mpl()
    c_n, c_v, dose = _extract_cluster_arrays(results, input_data)
    if c_n is None:
        print("plot_sia_clusters: skipped (bin-moment mode).")
        return None

    N = input_data.I
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
        fig, ax = plt.subplots(figsize=(7, 5))
        cmap = plt.cm.viridis(np.linspace(0, 1, len(idx)))
        for k, color in zip(idx, cmap):
            ax.loglog(dose, np.maximum(c_n[k, :], 1e-10),
                      color=color, label=f'n={k + 1}')
        ax.set_xlabel('Dose (dpa)')
        ax.set_ylabel(_CONC_LABEL)
        full_dose = results['dose']
        ax.set_title(f'{gtitle} {title}')
        ax.grid(True, which='both', alpha=0.3)
        ax.set_xlim(left=_dose_xlim(full_dose), right=full_dose[-1])
        _apply_axis_config(ax, 'concentration')
        fig.tight_layout()
        ncol = min(8, max(2, len(idx)))
        _legend_below(fig, ax, ncol=ncol, bottom=0.25)
        if out_path:
            # Replace file stem with group-specific name
            from pathlib import Path as _P
            p = _P(out_path)
            fig.savefig(str(p.parent / f'{fname}.png'), dpi=150,
                        bbox_inches='tight')
        figs.append(fig)
    return figs


def plot_vac_clusters(results, input_data, out_path=None, title=''):
    """Vacancy cluster concentrations vs. dose, split into small / mid / large."""
    _check_mpl()
    c_n, c_v, dose = _extract_cluster_arrays(results, input_data)
    if c_v is None:
        print("plot_vac_clusters: skipped (bin-moment mode).")
        return None

    M = input_data.V
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
        fig, ax = plt.subplots(figsize=(7, 5))
        cmap = plt.cm.plasma(np.linspace(0, 1, len(idx)))
        for k, color in zip(idx, cmap):
            ax.loglog(dose, np.maximum(c_v[k, :], 1e-10),
                      color=color, label=f'm={k + 1}')
        full_dose = results['dose']
        ax.set_xlabel('Dose (dpa)')
        ax.set_ylabel(_CONC_LABEL)
        ax.set_title(f'{gtitle} {title}')
        ax.grid(True, which='both', alpha=0.3)
        ax.set_xlim(left=_dose_xlim(full_dose), right=full_dose[-1])
        _apply_axis_config(ax, 'concentration')
        fig.tight_layout()
        ncol = min(8, max(2, len(idx)))
        _legend_below(fig, ax, ncol=ncol, bottom=0.25)
        if out_path:
            from pathlib import Path as _P
            p = _P(out_path)
            fig.savefig(str(p.parent / f'{fname}.png'), dpi=150,
                        bbox_inches='tight')
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

    Top panel: log-scale concentration vs cluster size n.
    Bottom panel: linear-scale concentration vs SIA loop diameter (nm).
    """
    _check_mpl()
    N     = input_data.I
    Omega = results.get('Omega', input_data.derived['Omega'])
    b_111 = results.get('b_111', input_data.derived['b_111'])
    y, dose = _align_dose_to_y(results)

    if y.shape[0] < N:
        print("plot_sia_distribution_evolution: skipped (bin-moment mode).")
        return None

    inv_O   = 1.0 / Omega
    ns      = np.arange(1, N + 1)
    d_loop  = 2.0 * np.sqrt(ns * Omega / (np.pi * b_111)) * 1e9  # nm
    indices = _log_snapshot_indices(dose, n_times)
    cmap    = plt.cm.viridis(np.linspace(0.15, 0.95, len(indices)))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 9))

    for idx, color in zip(indices, cmap):
        c_n = np.maximum(y[:N, idx], 0.0) * inv_O
        ax1.semilogy(ns, c_n + 1e-10, color=color,
                     label=f'{dose[idx]:.2e} dpa')
        ax2.plot(d_loop, c_n, color=color,
                 label=f'{dose[idx]:.2e} dpa')

    ax1.set_xlabel('SIA cluster size n')
    ax1.set_ylabel(_CONC_LABEL)
    ax1.set_title(f'SIA Size Distribution {title}')
    ax1.legend(fontsize=7, ncol=2, title='Dose')
    ax1.grid(True, which='both', alpha=0.3)

    ax2.set_xlabel('SIA loop diameter (nm)')
    ax2.set_ylabel(_CONC_LABEL)
    ax2.set_title(f'SIA Loop Size Distribution {title}')
    ax2.legend(fontsize=7, ncol=2, title='Dose')
    ax2.grid(True, alpha=0.3)

    _apply_axis_config_to_fig(fig, 'size_dist')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_void_distribution_evolution(results, input_data, n_times=10,
                                     out_path=None, title=''):
    """
    Void/vacancy cluster size distribution c_m(m) at n_times log-spaced dose snapshots.

    Top panel: log-scale concentration vs vacancy cluster size m.
    Bottom panel: linear-scale concentration vs void diameter (nm).
    """
    _check_mpl()
    N     = input_data.I
    M     = input_data.V
    Omega = results.get('Omega', input_data.derived['Omega'])
    r0    = results.get('r0', input_data.derived['r0'])
    y, dose = _align_dose_to_y(results)

    if y.shape[0] < N + M:
        print("plot_void_distribution_evolution: skipped (bin-moment mode).")
        return None

    inv_O   = 1.0 / Omega
    ms      = np.arange(1, M + 1)
    d_void  = 2.0 * r0 * ms ** (1.0 / 3.0) * 1e9  # nm
    indices = _log_snapshot_indices(dose, n_times)
    cmap    = plt.cm.plasma(np.linspace(0.15, 0.85, len(indices)))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 9))

    for idx, color in zip(indices, cmap):
        c_v = np.maximum(y[N:N + M, idx], 0.0) * inv_O
        ax1.semilogy(ms, c_v + 1e-10, color=color,
                     label=f'{dose[idx]:.2e} dpa')
        ax2.plot(d_void, c_v, color=color,
                 label=f'{dose[idx]:.2e} dpa')

    ax1.set_xlabel('Vacancy cluster size m')
    ax1.set_ylabel(_CONC_LABEL)
    ax1.set_title(f'Void Size Distribution {title}')
    ax1.legend(fontsize=7, ncol=2, title='Dose')
    ax1.grid(True, which='both', alpha=0.3)

    ax2.set_xlabel('Void diameter (nm)')
    ax2.set_ylabel(_CONC_LABEL)
    ax2.set_title(f'Void Diameter Distribution {title}')
    ax2.legend(fontsize=7, ncol=2, title='Dose')
    ax2.grid(True, alpha=0.3)

    _apply_axis_config_to_fig(fig, 'size_dist')
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
    ax.set_xlim(left=_dose_xlim(dose))
    _apply_axis_config(ax, 'concentration')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# TEM-visible cluster plots (n,m ≥ N_MIN_TEM)
# ══════════════════════════════════════════════════════════════════════════════

_N_MIN_TEM = 10       # atoms — minimum TEM-visible cluster size (~1 nm)
_OMEGA     = 1.18e-29 # m^3 — atomic volume (bcc Fe)
_B_111     = 2.482e-10 # m — Burgers vector 1/2<111>


def _tem_filter(c_arr, n_min=_N_MIN_TEM):
    """Return mask for TEM-visible clusters (size index >= n_min - 1)."""
    N = len(c_arr)
    return np.arange(1, N + 1) >= n_min


def plot_number_densities_tem(results, input_data, rate_eq_obj,
                              out_path=None, title=''):
    """
    Loop and void number densities vs dose — TEM-visible only (n,m ≥ 10).
    Plotted alongside the total (all sizes) for comparison.
    """
    _check_mpl()

    c_n_all, c_v_all, dose = _reconstruct_distributions(
        results, input_data, rate_eq_obj)
    if c_n_all is None:
        return None

    N, M = c_n_all.shape[0], c_v_all.shape[0]
    n_t = c_n_all.shape[1]
    mask_i = np.arange(1, N + 1) >= _N_MIN_TEM
    mask_v = np.arange(1, M + 1) >= _N_MIN_TEM

    N_loops_all = np.sum(c_n_all[1:, :], axis=0)   # n ≥ 2
    N_loops_tem = np.sum(c_n_all[mask_i, :], axis=0)
    N_voids_all = np.sum(c_v_all[1:, :], axis=0)
    N_voids_tem = np.sum(c_v_all[mask_v, :], axis=0)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(dose, np.maximum(N_loops_all, 1), color='steelblue',
              alpha=0.3, ls=':', label=r'SIA loops (all $n\geq2$)')
    ax.loglog(dose, np.maximum(N_loops_tem, 1), color='steelblue',
              lw=2, label=rf'SIA loops (TEM, $n\geq{_N_MIN_TEM}$)')
    ax.loglog(dose, np.maximum(N_voids_all, 1), color='tomato',
              alpha=0.3, ls=':', label=r'Voids (all $m\geq2$)')
    ax.loglog(dose, np.maximum(N_voids_tem, 1), color='tomato',
              lw=2, ls='--', label=rf'Voids (TEM, $m\geq{_N_MIN_TEM}$)')
    full_dose = results['dose']
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel(r'Number density (m$^{-3}$)')
    ax.set_title(f'TEM-Visible Number Densities {title}')
    ax.legend(fontsize=8)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(left=_dose_xlim(full_dose), right=full_dose[-1])
    _apply_axis_config(ax, 'concentration')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_mean_sizes_tem(results, input_data, rate_eq_obj,
                        out_path=None, title=''):
    """
    Mean cluster sizes vs dose — TEM-visible only (n,m ≥ 10).
    Shows both atoms (left axis) and physical diameter in nm (right axis).
    """
    _check_mpl()

    c_n_all, c_v_all, dose = _reconstruct_distributions(
        results, input_data, rate_eq_obj)
    if c_n_all is None:
        return None

    N, M = c_n_all.shape[0], c_v_all.shape[0]
    n_t = c_n_all.shape[1]
    ns = np.arange(1, N + 1, dtype=float)
    ms = np.arange(1, M + 1, dtype=float)
    mask_i = ns >= _N_MIN_TEM
    mask_v = ms >= _N_MIN_TEM

    Omega = input_data.derived.get('Omega', _OMEGA)
    b_111 = input_data.derived.get('b_111', _B_111)

    mean_n_i_tem = np.zeros(n_t)
    mean_n_v_tem = np.zeros(n_t)
    d_i_nm = np.zeros(n_t)
    d_v_nm = np.zeros(n_t)

    for j in range(n_t):
        cnt_i = np.sum(c_n_all[mask_i, j])
        cnt_v = np.sum(c_v_all[mask_v, j])
        mean_n_i_tem[j] = np.dot(ns[mask_i], c_n_all[mask_i, j]) / max(cnt_i, 1e-30)
        mean_n_v_tem[j] = np.dot(ms[mask_v], c_v_all[mask_v, j]) / max(cnt_v, 1e-30)
        d_i_nm[j] = 2 * np.sqrt(mean_n_i_tem[j] * Omega / (np.pi * b_111)) * 1e9
        d_v_nm[j] = 2 * (3 * mean_n_v_tem[j] * Omega / (4 * np.pi))**(1./3.) * 1e9

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # SIA loops
    ax1.semilogx(dose, d_i_nm, color='steelblue', lw=2,
                 label=rf'TEM-visible ($n\geq{_N_MIN_TEM}$)')
    ax1.axhline(10.0, color='gray', ls='--', alpha=0.5, label='Exp. target (10 nm)')
    ax1.set_xlabel('Dose (dpa)')
    full_dose = results['dose']
    ax1.set_ylabel('Mean SIA loop diameter (nm)')
    ax1.set_title(f'TEM-Visible SIA Loop Size {title}')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(left=_dose_xlim(full_dose), right=full_dose[-1])
    ax1.set_ylim(bottom=0)

    # Voids
    ax2.semilogx(dose, d_v_nm, color='tomato', lw=2,
                 label=rf'TEM-visible ($m\geq{_N_MIN_TEM}$)')
    ax2.axhline(3.0, color='gray', ls='--', alpha=0.5, label='Exp. target (3 nm)')
    ax2.set_xlabel('Dose (dpa)')
    ax2.set_ylabel('Mean void diameter (nm)')
    ax2.set_title(f'TEM-Visible Void Size {title}')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(left=_dose_xlim(full_dose), right=full_dose[-1])
    ax2.set_ylim(bottom=0)

    _apply_axis_config_to_fig(fig, 'scalar')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_sia_distribution_tem(results, input_data, rate_eq_obj,
                              n_times=8, out_path=None, title=''):
    """
    SIA loop size distribution — TEM-visible range only.
    Top: concentration vs cluster size (log-log, n ≥ N_MIN_TEM).
    Bottom: concentration vs loop diameter in nm.
    """
    _check_mpl()

    c_n_all, _, dose = _reconstruct_distributions(
        results, input_data, rate_eq_obj)
    if c_n_all is None:
        return None

    N = c_n_all.shape[0]
    Omega = input_data.derived.get('Omega', _OMEGA)
    b_111 = input_data.derived.get('b_111', _B_111)
    ns = np.arange(1, N + 1, dtype=float)
    d_loop = 2.0 * np.sqrt(ns * Omega / (np.pi * b_111)) * 1e9
    mask = ns >= _N_MIN_TEM

    indices = _log_snapshot_indices(dose, n_times)
    cmap = plt.cm.viridis(np.linspace(0.15, 0.95, len(indices)))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 9))

    for idx, color in zip(indices, cmap):
        c_n = c_n_all[:, idx]
        ax1.loglog(ns[mask], np.maximum(c_n[mask], 1e-10), color=color,
                   marker='.', ms=3, label=f'{dose[idx]:.2e} dpa')
        ax2.semilogy(d_loop[mask], np.maximum(c_n[mask], 1e-10), color=color,
                     marker='.', ms=3, label=f'{dose[idx]:.2e} dpa')

    ax1.set_xlabel('SIA cluster size n')
    ax1.set_ylabel(_CONC_LABEL)
    ax1.set_title(rf'TEM-Visible SIA Distribution ($n\geq{_N_MIN_TEM}$) {title}')
    ax1.legend(fontsize=7, ncol=2, title='Dose')
    ax1.grid(True, which='both', alpha=0.3)

    ax2.set_xlabel('SIA loop diameter (nm)')
    ax2.set_ylabel(_CONC_LABEL)
    ax2.set_title(rf'TEM-Visible SIA Loop Diameters ($n\geq{_N_MIN_TEM}$) {title}')
    ax2.legend(fontsize=7, ncol=2, title='Dose')
    ax2.grid(True, which='both', alpha=0.3)

    _apply_axis_config_to_fig(fig, 'size_dist')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def _reconstruct_distributions(results, input_data, rate_eq_obj):
    """
    Reconstruct per-size c_n [N, n_t] and c_v [M, n_t] in m^-3 from
    ODE state, handling both full_CD and bin-moment modes.

    Returns (c_n, c_v, dose) or (None, None, None) on failure.
    """
    from .bin_moment_rates import reconstruct_distribution

    N     = input_data.I
    M     = input_data.V
    y, dose = _align_dose_to_y(results)
    Omega = results.get('Omega', input_data.derived['Omega'])
    inv_O = 1.0 / Omega
    n_t   = y.shape[1]

    is_bin = hasattr(rate_eq_obj, 'bins')
    I_bin  = getattr(rate_eq_obj, 'I_bin', getattr(rate_eq_obj, 'K_i', 0))
    V_bin  = getattr(rate_eq_obj, 'V_bin', getattr(rate_eq_obj, 'K_v', 0))
    i_d    = getattr(rate_eq_obj, 'i_discrete', 0)
    v_d    = getattr(rate_eq_obj, 'v_discrete', 0)
    iv     = getattr(rate_eq_obj, 'i_VAC', N)
    P      = getattr(rate_eq_obj, 'n_mom', 2)
    sf     = getattr(rate_eq_obj, 'shape_function', 'linear')

    c_n_all = np.zeros((N, n_t))
    c_v_all = np.zeros((M, n_t))

    for j in range(n_t):
        yj = np.maximum(y[:, j], 0.0)
        if is_bin:
            # Discrete SIA sizes
            c_n_all[:i_d, j] = yj[:i_d] * inv_O
            # Binned SIA sizes (reconstruct from moments)
            if I_bin > 0:
                mom = yj[i_d:i_d + P * I_bin]
                mu0 = mom[0::P][:I_bin]
                mu1 = mom[1::P][:I_bin] if P >= 2 else None
                mu2 = mom[2::P][:I_bin] if P >= 3 else None
                c_binned = reconstruct_distribution(
                    sf, mu0, mu1, mu2, rate_eq_obj.bins, N) * inv_O
                c_n_all[i_d:, j] = c_binned[i_d:]
        else:
            c_n_all[:, j] = yj[:N] * inv_O

        if is_bin:
            # Discrete vacancy sizes
            c_v_all[:v_d, j] = yj[iv:iv + v_d] * inv_O
            # Binned vacancy sizes
            if V_bin > 0:
                vac_start = iv + v_d
                vmom = yj[vac_start:vac_start + P * V_bin]
                vmu0 = vmom[0::P][:V_bin]
                vmu1 = vmom[1::P][:V_bin] if P >= 2 else None
                vmu2 = vmom[2::P][:V_bin] if P >= 3 else None
                c_v_binned = reconstruct_distribution(
                    sf, vmu0, vmu1, vmu2,
                    rate_eq_obj.vac_bins, M) * inv_O
                c_v_all[v_d:, j] = c_v_binned[v_d:]
        else:
            c_v_all[:, j] = yj[iv:iv + M] * inv_O

    return c_n_all, c_v_all, dose


def _reconstruct_smooth_distributions(results, input_data, rate_eq_obj):
    """
    Reconstruct smooth (midpoint-based) size distributions for plotting.

    Uses the chosen shape function closure evaluated at geometric
    bin midpoints.

    Returns
    -------
    sia_mid : ndarray [K_i] — SIA bin midpoints
    sia_conc: ndarray [K_i, n_t] — SIA concentration at midpoints (m^-3)
    vac_mid : ndarray [K_v] — vacancy bin midpoints
    vac_conc: ndarray [K_v, n_t] — vacancy concentration at midpoints (m^-3)
    dose    : ndarray [n_t]
    """
    from .bin_moment_rates import midpoint_distribution_from_moments

    y, dose = _align_dose_to_y(results)
    Omega = results.get('Omega', input_data.derived['Omega'])
    inv_O = 1.0 / Omega
    n_t   = y.shape[1]

    is_bin = hasattr(rate_eq_obj, 'bins')
    I_bin  = getattr(rate_eq_obj, 'I_bin', getattr(rate_eq_obj, 'K_i', 0))
    V_bin  = getattr(rate_eq_obj, 'V_bin', getattr(rate_eq_obj, 'K_v', 0))
    i_d    = getattr(rate_eq_obj, 'i_discrete', 0)
    v_d    = getattr(rate_eq_obj, 'v_discrete', 0)
    iv     = getattr(rate_eq_obj, 'i_VAC', input_data.I)
    P      = getattr(rate_eq_obj, 'n_mom', 2)
    sf     = getattr(rate_eq_obj, 'shape_function', 'linear')

    sia_mid = np.zeros(I_bin)
    sia_conc = np.zeros((I_bin, n_t))
    vac_mid = np.zeros(V_bin)
    vac_conc = np.zeros((V_bin, n_t))

    for j in range(n_t):
        yj = np.maximum(y[:, j], 0.0)

        if is_bin and I_bin > 0:
            mom = yj[i_d:i_d + P * I_bin]
            mu0 = mom[0::P][:I_bin]
            mu1 = mom[1::P][:I_bin] if P >= 2 else None
            mu2 = mom[2::P][:I_bin] if P >= 3 else None
            mid, c = midpoint_distribution_from_moments(
                mu0, mu1, rate_eq_obj.bins,
                mu2=mu2, shape_function=sf)
            sia_mid = mid
            sia_conc[:, j] = c * inv_O

        if is_bin and V_bin > 0:
            vac_start = iv + v_d
            vmom = yj[vac_start:vac_start + P * V_bin]
            vmu0 = vmom[0::P][:V_bin]
            vmu1 = vmom[1::P][:V_bin] if P >= 2 else None
            vmu2 = vmom[2::P][:V_bin] if P >= 3 else None
            mid, c = midpoint_distribution_from_moments(
                vmu0, vmu1, rate_eq_obj.vac_bins,
                mu2=vmu2, shape_function=sf)
            vac_mid = mid
            vac_conc[:, j] = c * inv_O

    return sia_mid, sia_conc, vac_mid, vac_conc, dose


# ── Bin-moment-aware plots ────────────────────────────────────────────────────

def plot_sia_conc_vs_dose(results, input_data, rate_eq_obj,
                          out_path=None, title=''):
    """
    (1) SIA concentrations vs. dose — discrete sizes + bin mu0.

    Discrete sizes plotted individually; each bin as one curve.
    """
    _check_mpl()
    is_bin = hasattr(rate_eq_obj, 'bins')
    I_bin  = getattr(rate_eq_obj, 'I_bin', getattr(rate_eq_obj, 'K_i', 0))
    i_d    = getattr(rate_eq_obj, 'i_discrete', 0)
    Omega  = results.get('Omega', input_data.derived['Omega'])
    inv_O  = 1.0 / Omega
    y, dose = _align_dose_to_y(results)

    fig, ax = plt.subplots(figsize=(8, 6))
    n_curves = i_d + I_bin if is_bin else min(input_data.I, 20)
    cmap = plt.cm.viridis(np.linspace(0.1, 0.9, max(n_curves, 1)))

    if is_bin:
        ci = 0
        # Discrete SIA sizes (1..i_discrete)
        for n in range(1, i_d + 1):
            cn = np.maximum(y[n - 1, :], 0.0) * inv_O
            ax.loglog(dose, np.maximum(cn, 1e-10), color=cmap[ci],
                      label=f'n={n}')
            ci += 1
        # Binned SIA (mu0 for each bin)
        P = getattr(rate_eq_obj, 'n_mom', 2)
        for k in range(I_bin):
            nlo, nhi = rate_eq_obj.bins[k]
            mu0 = np.maximum(y[i_d + P*k, :], 0.0) * inv_O
            ax.loglog(dose, np.maximum(mu0, 1e-10), color=cmap[ci],
                      label=f'n={nlo}-{nhi-1}')
            ci += 1
    else:
        I = input_data.I
        ns_show = [1, 2, 5, 10, 20, 50, 100]
        cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(ns_show)))
        for i, n in enumerate(ns_show):
            if n > I:
                break
            ax.loglog(dose, np.maximum(y[n-1, :] * inv_O, 1e-10),
                      color=cmap[i], label=f'n={n}')

    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel(_CONC_LABEL)
    ax.set_title(f'SIA Cluster Concentrations {title}')
    ax.grid(True, which='both', alpha=0.3)
    full_dose = results['dose']
    ax.set_xlim(left=_dose_xlim(full_dose), right=full_dose[-1])
    _apply_axis_config(ax, 'concentration')
    fig.tight_layout()
    n_entries = len(ax.get_lines())
    ncol = min(8, max(2, (n_entries + 5) // 6))
    _legend_below(fig, ax, ncol=ncol, bottom=0.30)
    if out_path:
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
    return fig


def plot_vac_conc_vs_dose(results, input_data, rate_eq_obj,
                          out_path=None, title=''):
    """
    (2) Cavity concentrations vs. dose — discrete sizes + bin mu0.

    Left panel: zeroth moment mu0 (number density per bin).
    Right panel: vacancy content per bin (mu0 x midpoint).
    """
    _check_mpl()
    is_bin = hasattr(rate_eq_obj, 'bins')
    V_bin  = getattr(rate_eq_obj, 'V_bin', getattr(rate_eq_obj, 'K_v', 0))
    v_d    = getattr(rate_eq_obj, 'v_discrete', 0)
    iv     = getattr(rate_eq_obj, 'i_VAC', input_data.I)
    Omega  = results.get('Omega', input_data.derived['Omega'])
    inv_O  = 1.0 / Omega
    y, dose = _align_dose_to_y(results)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    n_curves = v_d + V_bin if is_bin else min(input_data.V, 20)
    cmap = plt.cm.plasma(np.linspace(0.1, 0.85, max(n_curves, 1)))

    if is_bin:
        ci = 0
        # Discrete vacancy sizes (1..v_discrete)
        for m in range(1, v_d + 1):
            cv = np.maximum(y[iv + m - 1, :], 0.0) * inv_O
            axes[0].loglog(dose, np.maximum(cv, 1e-10), color=cmap[ci],
                           label=f'm={m}')
            axes[1].loglog(dose, np.maximum(cv * m, 1e-10), color=cmap[ci],
                           label=f'm={m}')
            ci += 1
        # Binned vacancy (mu0 for each bin)
        if V_bin > 0:
            vac_mid = rate_eq_obj.vac_mid
            for k in range(V_bin):
                mlo, mhi = rate_eq_obj.vac_bins[k]
                mu0 = np.maximum(y[iv + v_d + k, :], 0.0) * inv_O
                axes[0].loglog(dose, np.maximum(mu0, 1e-10), color=cmap[ci],
                               label=f'm={mlo}-{mhi-1}')
                content = mu0 * vac_mid[k]
                axes[1].loglog(dose, np.maximum(content, 1e-10), color=cmap[ci],
                               label=f'm={mlo}-{mhi-1}')
                ci += 1
    else:
        V = input_data.V
        ms_show = [1, 2, 5, 10, 20, 50, 100]
        cmap = plt.cm.plasma(np.linspace(0.1, 0.85, len(ms_show)))
        for i, m in enumerate(ms_show):
            if m > V:
                break
            cv = np.maximum(y[iv + m - 1, :] * inv_O, 1e-10)
            axes[0].loglog(dose, cv, color=cmap[i], label=f'm={m}')
            axes[1].loglog(dose, cv * m, color=cmap[i], label=f'm={m}')

    full_dose = results['dose']
    axes[0].set_xlabel('Dose (dpa)')
    axes[0].set_ylabel(_CONC_LABEL)
    axes[0].set_title(r'Cavity $\mu_k^{(0)}$ (number density) ' + title)
    axes[0].grid(True, which='both', alpha=0.3)
    axes[0].set_xlim(left=_dose_xlim(full_dose), right=full_dose[-1])

    axes[1].set_xlabel('Dose (dpa)')
    axes[1].set_ylabel(r'Content $\mu_k^{(0)} \times \bar{m}_k$ (m$^{-3}$)')
    axes[1].set_title(r'Cavity content per bin ' + title)
    axes[1].grid(True, which='both', alpha=0.3)
    axes[1].set_xlim(left=_dose_xlim(full_dose), right=full_dose[-1])

    _apply_axis_config_to_fig(fig, 'concentration')
    fig.tight_layout()
    n_entries = len(axes[0].get_lines())
    ncol = min(10, max(2, (n_entries + 5) // 6))
    _legend_below(fig, axes, ncol=ncol, bottom=0.30)
    if out_path:
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
    return fig


def _quarter_indices(dose):
    """Return 4 time indices at ~25%, 50%, 75%, 100% of the log-dose range."""
    return _logspaced_indices(dose, 4)


def _logspaced_indices(dose, n_snap=8):
    """Return *n_snap* time indices log-spaced over the dose range."""
    d = np.asarray(dose)
    pos = d[d > 0]
    if len(pos) < n_snap:
        return list(range(min(n_snap, len(d))))
    log_min, log_max = np.log10(pos[0]), np.log10(d[-1])
    fracs = np.linspace(1.0 / n_snap, 1.0, n_snap)
    targets = [10**(log_min + f * (log_max - log_min)) for f in fracs]
    return [int(np.argmin(np.abs(d - t))) for t in targets]


def plot_sia_size_dist_quarters(results, input_data, rate_eq_obj,
                                out_path=None, title='', n_snap=8):
    """
    (3) SIA size distribution at *n_snap* log-spaced dose snapshots.

    Top panel: log-scale concentration vs cluster size n.
    Bottom panel: linear-scale concentration vs SIA loop diameter (nm).
    """
    _check_mpl()
    is_bin = hasattr(rate_eq_obj, 'bins') and getattr(rate_eq_obj, 'K_i', 0) > 0

    if is_bin:
        sia_mid, sia_conc, _, _, dose = _reconstruct_smooth_distributions(
            results, input_data, rate_eq_obj)
    else:
        c_n, _, dose = _reconstruct_distributions(results, input_data, rate_eq_obj)
        if c_n is None:
            return None

    Omega = results.get('Omega', input_data.derived['Omega'])
    b_111 = results.get('b_111', input_data.derived['b_111'])

    indices = _logspaced_indices(dose, n_snap)
    cmap = plt.cm.viridis(np.linspace(0.15, 0.90, len(indices)))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 9))
    for i, idx in enumerate(indices):
        if is_bin:
            ns = sia_mid
            cn = np.maximum(sia_conc[:, idx], 1e-10)
            ax1.loglog(ns, cn, '-o', color=cmap[i], lw=1.5, ms=3,
                       label=f'{dose[idx]:.2e} dpa')
        else:
            ns = np.arange(1, input_data.I + 1)
            cn = np.maximum(c_n[:, idx], 1e-10)
            ax1.semilogy(ns, cn, color=cmap[i], lw=1.5,
                         label=f'{dose[idx]:.2e} dpa')

        # Loop diameter: d = 2*sqrt(n*Omega / (pi*b))
        d_loop = 2.0 * np.sqrt(ns * Omega / (np.pi * b_111)) * 1e9  # nm
        ax2.plot(d_loop, cn, color=cmap[i], lw=1.5,
                 label=f'{dose[idx]:.2e} dpa')

    ax1.set_xlabel('SIA cluster size n')
    ax1.set_ylabel(_CONC_LABEL)
    ax1.set_title(f'SIA Size Distribution {title}')
    ax1.legend(title='Dose', fontsize=7)
    ax1.grid(True, which='both', alpha=0.3)

    ax2.set_xlabel('SIA loop diameter (nm)')
    ax2.set_ylabel(_CONC_LABEL)
    ax2.set_title(f'SIA Loop Size Distribution {title}')
    ax2.legend(title='Dose', fontsize=7)
    ax2.grid(True, alpha=0.3)

    _apply_axis_config_to_fig(fig, 'size_dist')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_vac_size_dist_quarters(results, input_data, rate_eq_obj,
                                out_path=None, title='', n_snap=8):
    """
    (4) Cavity size distribution at *n_snap* log-spaced dose snapshots.

    Top panel: log-scale concentration vs vacancy cluster size m.
    Bottom panel: linear-scale concentration vs void diameter (nm).
    """
    _check_mpl()
    is_bin = hasattr(rate_eq_obj, 'bins') and getattr(rate_eq_obj, 'K_v', 0) > 0

    if is_bin:
        _, _, vac_mid, vac_conc, dose = _reconstruct_smooth_distributions(
            results, input_data, rate_eq_obj)
    else:
        _, c_v, dose = _reconstruct_distributions(results, input_data, rate_eq_obj)
        if c_v is None:
            return None

    r0 = results.get('r0', input_data.derived['r0'])

    indices = _logspaced_indices(dose, n_snap)
    cmap = plt.cm.plasma(np.linspace(0.15, 0.90, len(indices)))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 9))
    for i, idx in enumerate(indices):
        if is_bin:
            ms = vac_mid
            cv = np.maximum(vac_conc[:, idx], 1e-10)
        else:
            ms = np.arange(1, input_data.V + 1)
            cv = np.maximum(c_v[:, idx], 1e-10)

        ax1.loglog(ms, cv, '-o', color=cmap[i], lw=1.5, ms=3,
                   label=f'{dose[idx]:.2e} dpa')

        # Void diameter: d = 2*r0*m^(1/3)
        d_void = 2.0 * r0 * ms ** (1.0 / 3.0) * 1e9  # nm
        ax2.plot(d_void, cv, color=cmap[i], lw=1.5,
                 label=f'{dose[idx]:.2e} dpa')

    ax1.set_xlabel('Cavity cluster size m')
    ax1.set_ylabel(_CONC_LABEL)
    ax1.set_title(f'Cavity Size Distribution {title}')
    ax1.legend(title='Dose', fontsize=7)
    ax1.grid(True, which='both', alpha=0.3)

    ax2.set_xlabel('Void diameter (nm)')
    ax2.set_ylabel(_CONC_LABEL)
    ax2.set_title(f'Cavity Diameter Distribution {title}')
    ax2.legend(title='Dose', fontsize=7)
    ax2.grid(True, alpha=0.3)

    _apply_axis_config_to_fig(fig, 'size_dist')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_sia_fraction_breakdown(results, out_path=None, title=''):
    """
    SIA content fraction breakdown vs. dose.

    Three curves that should sum to 1.0 at every dose:
      - In clusters (content in system / total produced)
      - At fixed sinks (cumulative fixed-sink loss / total produced)
      - Recombined + cavity (mutual annihilation / total produced)
    """
    _check_mpl()
    dose = results['dose']
    eta_G = results.get('eta_G', 0)
    t = results['t']
    prod = eta_G * t

    # All in at.frac — convert C_SIA_tot back from m^-3
    Omega = results.get('Omega', 1.18e-29)
    C_SIA = results['C_SIA_tot'] * Omega
    J_fixed = results.get('J_SIA_fixed', np.zeros_like(t))
    J_mutual = results.get('J_SIA_mutual', np.zeros_like(t))

    fig, ax = plt.subplots(figsize=(8, 5))
    mask = prod > 1e-300
    f_clusters = np.where(mask, C_SIA / prod, 0.0)
    f_fixed    = np.where(mask, J_fixed / prod, 0.0)
    f_mutual   = np.where(mask, J_mutual / prod, 0.0)
    f_total    = f_clusters + f_fixed + f_mutual

    ax.semilogx(dose[mask], f_clusters[mask], color='steelblue', lw=2,
                label='In clusters')
    ax.semilogx(dose[mask], f_fixed[mask], color='firebrick', lw=2,
                label='Fixed sinks')
    ax.semilogx(dose[mask], f_mutual[mask], color='forestgreen', lw=2,
                label='Recombination + cavity')
    ax.semilogx(dose[mask], f_total[mask], color='black', lw=1, ls='--',
                label='Sum (should = 1)')

    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel('Fraction of total SIA produced')
    ax.set_title(f'SIA Content Balance {title}')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(left=_dose_xlim(dose))
    _apply_axis_config(ax, 'scalar')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def plot_vac_fraction_breakdown(results, out_path=None, title=''):
    """
    Vacancy content fraction breakdown vs. dose.

    Three curves that should sum to 1.0:
      - In clusters (swelling content / total produced)
      - At fixed sinks (cumulative fixed-sink loss / total produced)
      - Annihilated by SIA (SIA-cavity + recombination / total produced)
    """
    _check_mpl()
    dose = results['dose']
    eta_G = results.get('eta_G', 0)
    t = results['t']
    prod = eta_G * t

    Omega = results.get('Omega', 1.18e-29)
    C_VAC = results['C_VAC_tot'] * Omega
    J_fixed = results.get('J_VAC_fixed', np.zeros_like(t))
    J_mutual = results.get('J_VAC_mutual', np.zeros_like(t))  # VAC-specific mutual

    fig, ax = plt.subplots(figsize=(8, 5))
    mask = prod > 1e-300
    f_clusters = np.where(mask, C_VAC / prod, 0.0)
    f_fixed    = np.where(mask, J_fixed / prod, 0.0)
    f_mutual   = np.where(mask, J_mutual / prod, 0.0)
    f_total    = f_clusters + f_fixed + f_mutual

    ax.semilogx(dose[mask], f_clusters[mask], color='tomato', lw=2,
                label='In clusters (swelling)')
    ax.semilogx(dose[mask], f_fixed[mask], color='firebrick', lw=2,
                label='Fixed sinks')
    ax.semilogx(dose[mask], f_mutual[mask], color='forestgreen', lw=2,
                label='Annihilated by SIA')
    ax.semilogx(dose[mask], f_total[mask], color='black', lw=1, ls='--',
                label='Sum (should = 1)')

    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel('Fraction of total VAC produced')
    ax.set_title(f'Vacancy Content Balance {title}')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(left=_dose_xlim(dose))
    _apply_axis_config(ax, 'scalar')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def save_all_plots(results, input_data, out_dir, label='',
                   rate_eq_obj=None):
    """Save all standard plots to out_dir/."""
    from pathlib import Path
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    opts = dict(title=label)
    plot_sia_fraction_breakdown(results, out_path=f"{out_dir}/sia_balance.png", **opts)
    plot_vac_fraction_breakdown(results, out_path=f"{out_dir}/vac_balance.png", **opts)
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

    # (1)-(4) Bin-moment-aware plots (work for both full_CD and bin_moment)
    if rate_eq_obj is not None:
        plot_sia_conc_vs_dose(results, input_data, rate_eq_obj,
                              out_path=f"{out_dir}/sia_bin_conc.png", **opts)
        plot_vac_conc_vs_dose(results, input_data, rate_eq_obj,
                              out_path=f"{out_dir}/vac_bin_conc.png", **opts)
        plot_sia_size_dist_quarters(results, input_data, rate_eq_obj,
                                    out_path=f"{out_dir}/sia_dist_quarters.png", **opts)
        plot_vac_size_dist_quarters(results, input_data, rate_eq_obj,
                                    out_path=f"{out_dir}/vac_dist_quarters.png", **opts)

        # TEM-visible plots (n,m ≥ 10)
        plot_number_densities_tem(results, input_data, rate_eq_obj,
                                  out_path=f"{out_dir}/number_densities_tem.png", **opts)
        plot_mean_sizes_tem(results, input_data, rate_eq_obj,
                            out_path=f"{out_dir}/mean_sizes_tem.png", **opts)
        plot_sia_distribution_tem(results, input_data, rate_eq_obj,
                                  out_path=f"{out_dir}/sia_dist_tem.png", **opts)

    print(f"Saved plots to {out_dir}/")
