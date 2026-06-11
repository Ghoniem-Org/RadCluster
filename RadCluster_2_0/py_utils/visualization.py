"""
visualization.py — Plotting routines for RadCluster_2_0.

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
_INTERIOR_LEGEND_FONTSIZE = 12
if _HAS_MPL:
    plt.rcParams.update({
        'axes.titlesize':       _PLOT_FONTSIZE,
        'axes.labelsize':       _PLOT_FONTSIZE,
        'xtick.labelsize':      _PLOT_FONTSIZE,
        'ytick.labelsize':      _PLOT_FONTSIZE,
        'legend.fontsize':      _PLOT_FONTSIZE,
        'legend.title_fontsize': _PLOT_FONTSIZE,
        'figure.titlesize':     _PLOT_FONTSIZE,
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
    'legend_fontsize': _PLOT_FONTSIZE,
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
    ax.legend(loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
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
    ax.legend(loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
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
    ax.legend(loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
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
    ax.legend(loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(left=_dose_xlim(dose))
    _apply_axis_config(ax, 'concentration')
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
    ax1.legend(ncol=2, title='Dose', loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
    ax1.grid(True, which='both', alpha=0.3)

    ax2.set_xlabel('SIA loop diameter (nm)')
    ax2.set_ylabel(_CONC_LABEL)
    ax2.set_title(f'SIA Loop Size Distribution {title}')
    ax2.legend(ncol=2, title='Dose', loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
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
    ax1.legend(ncol=2, title='Dose', loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
    ax1.grid(True, which='both', alpha=0.3)

    ax2.set_xlabel('Void diameter (nm)')
    ax2.set_ylabel(_CONC_LABEL)
    ax2.set_title(f'Void Diameter Distribution {title}')
    ax2.legend(ncol=2, title='Dose', loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
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
    ax.legend(loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
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
    ax.legend(loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
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
    C_floor = float(input_data.reactions.get('C_floor', 1e-15)) / Omega

    mean_n_i_tem = np.zeros(n_t)
    mean_n_v_tem = np.zeros(n_t)
    d_i_nm = np.zeros(n_t)
    d_v_nm = np.zeros(n_t)

    for j in range(n_t):
        ci_eff = np.maximum(c_n_all[mask_i, j] - C_floor, 0.0)
        cv_eff = np.maximum(c_v_all[mask_v, j] - C_floor, 0.0)
        cnt_i = np.sum(ci_eff)
        cnt_v = np.sum(cv_eff)
        mean_n_i_tem[j] = np.dot(ns[mask_i], ci_eff) / cnt_i if cnt_i > 0 else 0.0
        mean_n_v_tem[j] = np.dot(ms[mask_v], cv_eff) / cnt_v if cnt_v > 0 else 0.0
        d_i_nm[j] = 2 * np.sqrt(mean_n_i_tem[j] * Omega / (np.pi * b_111)) * 1e9
        d_v_nm[j] = 2 * (3 * mean_n_v_tem[j] * Omega / (4 * np.pi))**(1./3.) * 1e9

    fig, ax = plt.subplots(figsize=(8, 6))
    full_dose = results['dose']

    ax.plot(dose, d_i_nm, color='steelblue', lw=2,
            label=rf'SIA loops ($n\geq{_N_MIN_TEM}$)')
    ax.plot(dose, d_v_nm, color='tomato', lw=2,
            label=rf'Voids ($m\geq{_N_MIN_TEM}$)')
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel('Mean diameter (nm)')
    ax.set_title(f'TEM-Visible Mean Cluster Sizes {title}')
    ax.legend(loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0, right=full_dose[-1])
    ax.set_ylim(bottom=0)

    _apply_axis_config_to_fig(fig, 'scalar')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


def _visualization_bins(rate_eq_obj, max_size, kind):
    """
    Build the visualization bin list covering sizes 1..max_size.

    Each entry is (n_lo, n_hi) with n_hi exclusive (matching the
    bin_moment_rates convention).  For full_CD modes every integer is
    its own width-1 bin.  For bin_moment modes the discrete portion is
    width-1 integers and the grouped portion is taken from
    rate_eq_obj.bins (kind='sia') or rate_eq_obj.vac_bins (kind='vac').
    """
    is_bin = hasattr(rate_eq_obj, 'bins')
    if not is_bin:
        return [(n, n + 1) for n in range(1, max_size + 1)]
    if kind == 'sia':
        d_disc = getattr(rate_eq_obj, 'i_discrete', 0)
        grouped = list(rate_eq_obj.bins)
    else:
        d_disc = getattr(rate_eq_obj, 'v_discrete', 0)
        grouped = list(rate_eq_obj.vac_bins)
    out = [(n, n + 1) for n in range(1, d_disc + 1)]
    out.extend(grouped)
    return out


def _aggregate_to_bins(c_per_size, bins_full):
    """
    Sum per-size concentrations over each visualization bin.

    Returns mu0 of shape [K, n_t] where K = len(bins_full).  For width-1
    bins this is just c_n; for grouped bins it is the bin's zeroth
    moment (the lognormal/linear closure preserves this exactly).
    """
    n_t = c_per_size.shape[1]
    K   = len(bins_full)
    mu0 = np.zeros((K, n_t))
    for k, (nlo, nhi) in enumerate(bins_full):
        mu0[k] = c_per_size[nlo - 1:nhi - 1].sum(axis=0)
    return mu0


def _plot_tem_density_panel(bins_full, mu0_arr, dose, *, n_times,
                            xform, n_min, xlabel, ylabel, plot_title,
                            cmap_name, use_logx, out_path):
    """
    TEM-visible per-bin density plot.

    Each visualization bin is rendered as a stair segment of height
    rho_k = mu_0(k) / (x_hi - x_lo) where the x-edges come from
    `xform` (identity for size axis, diameter formula for diameter
    axes).  This is the standard microstructure-distribution
    representation: integral under the stairs equals total
    concentration in the displayed range.  Per-bin density is well
    defined regardless of the intra-bin shape function, so the comb
    artefact of integer-resolved lognormal reconstruction is avoided
    by construction.
    """
    keep = [k for k, (nlo, _nhi) in enumerate(bins_full) if nlo >= n_min]
    if not keep:
        return None

    edges_x = np.array(
        [xform(bins_full[keep[0]][0])]
        + [xform(bins_full[k][1]) for k in keep],
        dtype=float)
    widths = np.diff(edges_x)
    widths = np.where(widths > 0.0, widths, 1.0)

    indices = _log_snapshot_indices(dose, n_times)
    cmap = getattr(plt.cm, cmap_name)(np.linspace(0.15, 0.95, len(indices)))

    fig, ax = plt.subplots(figsize=(8, 6))
    for idx, color in zip(indices, cmap):
        rho = mu0_arr[keep, idx] / widths
        rho = np.maximum(rho, 1e-30)
        ax.stairs(rho, edges_x, color=color, lw=1.5, baseline=None,
                  label=f'{dose[idx]:.2e} dpa')

    if use_logx:
        ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(plot_title)
    ax.legend(ncol=2, title='Dose', loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
    ax.grid(True, which='both', alpha=0.3)

    _apply_axis_config_to_fig(fig, 'size_dist')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


_DENSITY_LABEL_SIZE = r'$dc/dn$ (m$^{-3}$)'
_DENSITY_LABEL_DIAM = r'$dc/dD$ (m$^{-3}$ nm$^{-1}$)'


def plot_sia_distribution_tem_size(results, input_data, rate_eq_obj,
                                   n_times=8, out_path=None, title=''):
    """SIA distribution density vs cluster size (TEM-visible, n ≥ N_MIN_TEM)."""
    _check_mpl()
    c_n_all, _, dose = _reconstruct_distributions(
        results, input_data, rate_eq_obj)
    if c_n_all is None:
        return None
    N         = c_n_all.shape[0]
    bins_full = _visualization_bins(rate_eq_obj, N, kind='sia')
    mu0_arr   = _aggregate_to_bins(c_n_all, bins_full)
    return _plot_tem_density_panel(
        bins_full, mu0_arr, dose,
        n_times=n_times,
        xform=lambda n: float(n),
        n_min=_N_MIN_TEM,
        xlabel='SIA cluster size n',
        ylabel=_DENSITY_LABEL_SIZE,
        plot_title=rf'TEM-Visible SIA Distribution ($n\geq{_N_MIN_TEM}$) {title}',
        cmap_name='viridis', use_logx=True, out_path=out_path)


def plot_sia_distribution_tem_diameter(results, input_data, rate_eq_obj,
                                       n_times=8, out_path=None, title=''):
    """SIA distribution density vs loop diameter (TEM-visible, n ≥ N_MIN_TEM)."""
    _check_mpl()
    c_n_all, _, dose = _reconstruct_distributions(
        results, input_data, rate_eq_obj)
    if c_n_all is None:
        return None
    N         = c_n_all.shape[0]
    Omega     = input_data.derived.get('Omega', _OMEGA)
    b_111     = input_data.derived.get('b_111', _B_111)
    pref      = 2.0 * np.sqrt(Omega / (np.pi * b_111)) * 1e9
    bins_full = _visualization_bins(rate_eq_obj, N, kind='sia')
    mu0_arr   = _aggregate_to_bins(c_n_all, bins_full)
    return _plot_tem_density_panel(
        bins_full, mu0_arr, dose,
        n_times=n_times,
        xform=lambda n: pref * np.sqrt(float(n)),
        n_min=_N_MIN_TEM,
        xlabel='SIA loop diameter (nm)',
        ylabel=_DENSITY_LABEL_DIAM,
        plot_title=rf'TEM-Visible SIA Loop Diameters ($n\geq{_N_MIN_TEM}$) {title}',
        cmap_name='viridis', use_logx=False, out_path=out_path)


def plot_vac_distribution_tem_size(results, input_data, rate_eq_obj,
                                   n_times=8, out_path=None, title=''):
    """Cavity distribution density vs cluster size (TEM-visible, m ≥ N_MIN_TEM)."""
    _check_mpl()
    _, c_v_all, dose = _reconstruct_distributions(
        results, input_data, rate_eq_obj)
    if c_v_all is None:
        return None
    M         = c_v_all.shape[0]
    bins_full = _visualization_bins(rate_eq_obj, M, kind='vac')
    mu0_arr   = _aggregate_to_bins(c_v_all, bins_full)
    return _plot_tem_density_panel(
        bins_full, mu0_arr, dose,
        n_times=n_times,
        xform=lambda m: float(m),
        n_min=_N_MIN_TEM,
        xlabel='Cavity cluster size m',
        ylabel=_DENSITY_LABEL_SIZE,
        plot_title=rf'TEM-Visible Cavity Distribution ($m\geq{_N_MIN_TEM}$) {title}',
        cmap_name='plasma', use_logx=True, out_path=out_path)


def plot_vac_distribution_tem_diameter(results, input_data, rate_eq_obj,
                                       n_times=8, out_path=None, title=''):
    """Cavity distribution density vs void diameter (TEM-visible, m ≥ N_MIN_TEM)."""
    _check_mpl()
    _, c_v_all, dose = _reconstruct_distributions(
        results, input_data, rate_eq_obj)
    if c_v_all is None:
        return None
    M         = c_v_all.shape[0]
    r0        = results.get('r0', input_data.derived['r0'])
    pref      = 2.0 * r0 * 1e9
    bins_full = _visualization_bins(rate_eq_obj, M, kind='vac')
    mu0_arr   = _aggregate_to_bins(c_v_all, bins_full)
    return _plot_tem_density_panel(
        bins_full, mu0_arr, dose,
        n_times=n_times,
        xform=lambda m: pref * float(m) ** (1.0 / 3.0),
        n_min=_N_MIN_TEM,
        xlabel='Void diameter (nm)',
        ylabel=_DENSITY_LABEL_DIAM,
        plot_title=rf'TEM-Visible Cavity Diameters ($m\geq{_N_MIN_TEM}$) {title}',
        cmap_name='plasma', use_logx=False, out_path=out_path)


def _reconstruct_distributions(results, input_data, rate_eq_obj):
    """
    Reconstruct per-size c_n [N, n_t] and c_v [M, n_t] in m^-3 from
    ODE state, handling both full_CD and bin-moment modes.

    Returns (c_n, c_v, dose) or (None, None, None) on failure.
    """
    from .bin_moment_rates import (reconstruct_distribution,
                                    distribution_from_moments_continuous)

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

    def _bridge_discrete_to_bin(c_arr, j_idx, i_disc, bins_list):
        # Visualization-only smoothing across the discrete↔bin transition.
        # Adjusts only the leftmost-bin's first value (closure side) and
        # rescales that bin so its zeroth moment is preserved.
        if i_disc < 1 or not bins_list:
            return
        nlo, nhi = bins_list[0]
        lo, hi = nlo - 1, nhi - 1   # 0-indexed bin extent (hi exclusive)
        if hi <= lo:
            return
        v_disc = c_arr[i_disc - 1, j_idx]
        v_bin  = c_arr[lo, j_idx]
        if v_disc <= 0.0 or v_bin <= 0.0:
            return
        pre_sum = c_arr[lo:hi, j_idx].sum()
        if pre_sum <= 0.0:
            return
        c_arr[lo, j_idx] = np.sqrt(v_disc * v_bin)
        new_sum = c_arr[lo:hi, j_idx].sum()
        if new_sum > 0.0:
            c_arr[lo:hi, j_idx] *= pre_sum / new_sum

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
                    sf, mu0, mu1, mu2, rate_eq_obj.bins, N,
                    smooth_edges=True) * inv_O
                c_n_all[i_d:, j] = c_binned[i_d:]
                _bridge_discrete_to_bin(c_n_all, j, i_d, rate_eq_obj.bins)
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
                # Vacancy distribution decays steeply, so the linear hat
                # function clamps at bin edges and produces sawtooth cliffs.
                # For viz only, use the C⁰-continuous reconstruction when
                # the linear closure is active.
                if sf == 'linear':
                    c_v_binned = distribution_from_moments_continuous(
                        vmu0, vmu1, rate_eq_obj.vac_bins, M) * inv_O
                else:
                    c_v_binned = reconstruct_distribution(
                        sf, vmu0, vmu1, vmu2,
                        rate_eq_obj.vac_bins, M, smooth_edges=True) * inv_O
                c_v_all[v_d:, j] = c_v_binned[v_d:]
                _bridge_discrete_to_bin(c_v_all, j, v_d, rate_eq_obj.vac_bins)
        else:
            c_v_all[:, j] = yj[iv:iv + M] * inv_O

    return c_n_all, c_v_all, dose


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


def _plot_vac_moment_vs_dose(results, input_data, rate_eq_obj,
                             *, weight, ylabel, plot_title, out_path):
    """
    Cavity moment vs. dose — discrete sizes + bin mu0, weighted by `weight`.

    weight='mu0'    → plot mu0 (number density per bin).
    weight='content'→ plot mu0 * (m or midpoint) (vacancy content per bin).
    """
    _check_mpl()
    is_bin = hasattr(rate_eq_obj, 'bins')
    V_bin  = getattr(rate_eq_obj, 'V_bin', getattr(rate_eq_obj, 'K_v', 0))
    v_d    = getattr(rate_eq_obj, 'v_discrete', 0)
    iv     = getattr(rate_eq_obj, 'i_VAC', input_data.I)
    Omega  = results.get('Omega', input_data.derived['Omega'])
    inv_O  = 1.0 / Omega
    y, dose = _align_dose_to_y(results)
    use_content = (weight == 'content')

    fig, ax = plt.subplots(figsize=(8, 6))
    n_curves = v_d + V_bin if is_bin else min(input_data.V, 20)
    cmap = plt.cm.plasma(np.linspace(0.1, 0.85, max(n_curves, 1)))

    if is_bin:
        ci = 0
        for m in range(1, v_d + 1):
            cv = np.maximum(y[iv + m - 1, :], 0.0) * inv_O
            vals = cv * m if use_content else cv
            ax.loglog(dose, np.maximum(vals, 1e-10), color=cmap[ci],
                      label=f'm={m}')
            ci += 1
        if V_bin > 0:
            vac_mid = rate_eq_obj.vac_mid
            for k in range(V_bin):
                mlo, mhi = rate_eq_obj.vac_bins[k]
                mu0 = np.maximum(y[iv + v_d + k, :], 0.0) * inv_O
                vals = mu0 * vac_mid[k] if use_content else mu0
                ax.loglog(dose, np.maximum(vals, 1e-10), color=cmap[ci],
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
            vals = cv * m if use_content else cv
            ax.loglog(dose, vals, color=cmap[i], label=f'm={m}')

    full_dose = results['dose']
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel(ylabel)
    ax.set_title(plot_title)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim(left=_dose_xlim(full_dose), right=full_dose[-1])

    _apply_axis_config(ax, 'concentration')
    fig.tight_layout()
    n_entries = len(ax.get_lines())
    ncol = min(10, max(2, (n_entries + 5) // 6))
    _legend_below(fig, ax, ncol=ncol, bottom=0.30)
    if out_path:
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
    return fig


def plot_vac_conc_vs_dose(results, input_data, rate_eq_obj,
                          out_path=None, title=''):
    """Cavity zeroth moment μ_k^(0) (number density per bin) vs. dose."""
    return _plot_vac_moment_vs_dose(
        results, input_data, rate_eq_obj,
        weight='mu0',
        ylabel=_CONC_LABEL,
        plot_title=f'Cavity Cluster Concentrations {title}',
        out_path=out_path)


def plot_vac_content_vs_dose(results, input_data, rate_eq_obj,
                             out_path=None, title=''):
    """Cavity content per bin (μ_k^(0) × m̄_k) vs. dose."""
    return _plot_vac_moment_vs_dose(
        results, input_data, rate_eq_obj,
        weight='content',
        ylabel=r'Content $\mu_k^{(0)} \times \bar{m}_k$ (m$^{-3}$)',
        plot_title=r'Cavity content per bin ' + title,
        out_path=out_path)


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
    ax.legend(loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
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
    ax.legend(loc='best', fontsize=_INTERIOR_LEGEND_FONTSIZE)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(left=_dose_xlim(dose))
    _apply_axis_config(ax, 'scalar')
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=150)
    return fig


PLOT_DATA_FILENAME = 'plot_data.pkl'


def _atomic_pickle_dump(payload, out_path):
    """Pickle to a sibling .tmp file, fsync, then atomic-rename to out_path.

    A crash, Ctrl+C, or OOM during the dump leaves the .tmp behind (cleaned
    up here on its way out) and never touches an existing good file at
    out_path.  Callers can therefore re-load the previous run's pickle even
    if the latest save was interrupted.
    """
    import os, pickle
    out_path = str(out_path)
    tmp_path = out_path + '.tmp'
    try:
        with open(tmp_path, 'wb') as f:
            pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp_path, out_path)
    except BaseException:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise


def save_plot_data(results, input_data, out_path, *, rate_eq_obj=None, label=''):
    """
    Pickle everything needed to regenerate the plots later.

    The dump contains the full `results` dict, the `InputData` object, the
    rate-equation object (if provided), and the title `label` so a notebook
    can call any plot function later without re-running the simulation.

    Writes are atomic (temp file + os.replace) so an interrupt mid-dump
    never leaves a truncated pickle at `out_path`.
    """
    payload = {
        'results':     results,
        'input_data':  input_data,
        'rate_eq_obj': rate_eq_obj,
        'label':       label,
    }
    _atomic_pickle_dump(payload, out_path)


def load_plot_data(path):
    """Inverse of save_plot_data — returns the dict that was pickled."""
    import pickle
    with open(path, 'rb') as f:
        return pickle.load(f)


def _safe_save(artifact, fn, *args, **kwargs):
    """Run a save call, catching/printing any failure with a traceback.

    KeyboardInterrupt is NOT caught here — it propagates so the outer
    interrupt shield (see simulation._InterruptDeferred) can decide whether
    to defer it through the rest of the save or escalate.
    """
    try:
        fn(*args, **kwargs)
        return True
    except KeyboardInterrupt:
        raise
    except Exception as exc:
        import traceback
        print(f"  ✗ {artifact}: {type(exc).__name__}: {exc}")
        traceback.print_exc()
        return False


def save_all_plots(results, input_data, out_dir, label='',
                   rate_eq_obj=None):
    """Save all standard plots to out_dir/.

    Each artifact (the pickle and every PNG) is wrapped in its own
    try/except so a single failure cannot abort the rest.  The pickle uses
    atomic writes; if it fails, any pre-existing good pickle is preserved.
    """
    from pathlib import Path
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Persist plotting inputs (atomic; previous good file is preserved on failure).
    _safe_save(PLOT_DATA_FILENAME,
               save_plot_data, results, input_data,
               f"{out_dir}/{PLOT_DATA_FILENAME}",
               rate_eq_obj=rate_eq_obj, label=label)

    opts = dict(title=label)
    common = [
        ('sia_balance.png',         plot_sia_fraction_breakdown,      (results,)),
        ('vac_balance.png',         plot_vac_fraction_breakdown,      (results,)),
        ('point_defects.png',       plot_point_defects,               (results,)),
        ('totals.png',              plot_totals,                      (results,)),
        ('swelling.png',            plot_swelling,                    (results,)),
        ('mean_sizes.png',          plot_mean_sizes,                  (results,)),
        ('he_content.png',          plot_he_content,                  (results,)),
        ('number_densities.png',    plot_number_densities,            (results,)),
        ('size_distributions.png',  plot_size_distribution,           (results, input_data)),
        ('sia_dist_evolution.png',  plot_sia_distribution_evolution,  (results, input_data)),
        ('void_dist_evolution.png', plot_void_distribution_evolution, (results, input_data)),
        # Per-cluster time-series (full_CD modes only; skipped for bin_moment)
        ('sia_small.png',           plot_sia_clusters,                (results, input_data)),
        ('vac_small.png',           plot_vac_clusters,                (results, input_data)),
    ]
    for name, fn, args in common:
        _safe_save(name, fn, *args, out_path=f"{out_dir}/{name}", **opts)

    # (1)-(4) Bin-moment-aware plots (work for both full_CD and bin_moment)
    if rate_eq_obj is not None:
        bin_aware = [
            ('sia_bin_conc.png',          plot_sia_conc_vs_dose),
            ('vac_bin_conc.png',          plot_vac_conc_vs_dose),
            ('vac_bin_content.png',       plot_vac_content_vs_dose),
            # TEM-visible plots (n,m ≥ 10)
            ('number_densities_tem.png',  plot_number_densities_tem),
            ('mean_sizes_tem.png',        plot_mean_sizes_tem),
            ('sia_dist_tem_size.png',     plot_sia_distribution_tem_size),
            ('sia_dist_tem_diameter.png', plot_sia_distribution_tem_diameter),
            ('vac_dist_tem_size.png',     plot_vac_distribution_tem_size),
            ('vac_dist_tem_diameter.png', plot_vac_distribution_tem_diameter),
        ]
        for name, fn in bin_aware:
            _safe_save(name, fn, results, input_data, rate_eq_obj,
                       out_path=f"{out_dir}/{name}", **opts)

    print(f"Saved plots to {out_dir}/")
