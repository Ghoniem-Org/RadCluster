"""
visualization.py – Plotting and provenance for ClusterDynamics.

Mirrors the ZrMicroVisualizer pattern: one PNG per figure, saved to a
unique timestamped run directory under output/.

Figures produced
----------------
1. point_defects.png               – Cv1 & Ci1 vs time
2. interstitial_clusters_small.png – small  interstitial clusters vs time  (log-log)
3. interstitial_clusters_mid.png   – mid    interstitial clusters vs time  (log-log)
4. interstitial_clusters_large.png – large  interstitial clusters vs time  (log-log)
5. vacancy_clusters_small.png      – small  vacancy    clusters vs time  (log-log)
6. vacancy_clusters_mid.png        – mid    vacancy    clusters vs time  (log-log)
7. vacancy_clusters_large.png      – large  vacancy    clusters vs time  (log-log)
 8. i_cluster_size_early.png        – interstitial density vs size, early times (log-x, linear-y)
 9. i_cluster_size_mid.png          – interstitial density vs size, mid   times (log-x, linear-y)
10. i_cluster_size_late.png         – interstitial density vs size, late  times (log-x, linear-y)
11. v_cluster_size_early.png        – vacancy    density vs size, early times (log-x, linear-y)
12. v_cluster_size_mid.png          – vacancy    density vs size, mid   times (log-x, linear-y)
13. v_cluster_size_late.png         – vacancy    density vs size, late  times (log-x, linear-y)
14. i_cluster_radius_small.png      – interstitial density vs radius R [nm], small size range  (linear-x)
15. i_cluster_radius_mid.png        – interstitial density vs radius R [nm], mid   size range  (linear-x)
16. i_cluster_radius_large.png      – interstitial density vs radius R [nm], large size range  (linear-x)

Size-regime split (Figures 2–7)
--------------------------------
The full log-range [2, N] is divided into 3 equal-decade bands (small / mid /
large). ~10 log-spaced representative sizes are picked per band.  y-min is fixed
at C_floor (default 1e-20); y-max is auto-scaled to the data maximum in each band.

Time-regime split (Figures 8–13)
---------------------------------
The simulated time span is divided into 3 equal-decade bands (early / mid / late).
~10 log-spaced representative snapshots are picked per band.  Mid and late regimes
start 2 decades below their respective band edge.  x-axis is log-scale; y-axis is
linear (mask out zero density; xlim/ylim overridable via plot_config).

Radius size-distribution plots (Figures 14–16)
-----------------------------------------------
Same 3 size-regime bands (small / mid / large) as Figures 2–7, but x-axis is the
cluster radius R [nm] on a **linear** scale instead of cluster size.  Radius is
computed assuming a spherical cluster and BCC alpha-iron atomic volume:
    Omega_Fe = a_Fe^3 / 2   (a_Fe = 2.87 Å)
    R = (3 n Omega_Fe / 4 pi)^(1/3)
Multiple time snapshots (combined from early / mid / late regimes) are overlaid on
each plot.  Spline interpolation is performed in log-y / linear-x space for smooth
curves.

Public entry point:
    plot_results(sim, results, output_dir, save_plots=True, plot_config=None)

plot_config structure (all keys optional)
-----------------------------------------------------------------
{
    'point_defects':         {'xlim': (1e-6, 1e6), 'ylim': (1e-20, 1e-3)},
    'interstitial_clusters': {'xlim': (1e-6, 1e6), 'C_floor': 1e-20},
    'vacancy_clusters':      {'xlim': (1e-6, 1e6), 'C_floor': 1e-20},
    'i_cluster_size':        {'xlim': None, 'ylim': None},
    'v_cluster_size':        {'xlim': None, 'ylim': None},
    'i_cluster_radius':      {'ylim': None},   # linear-x radius plots (Figures 14–16)
}
"""

import datetime
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.interpolate import make_interp_spline


# ── Window-mode name mapping ──────────────────────────────────────────────────

_WINDOW_MODE_NAMES = {
    0: 'cpp_full',
    1: 'cpp_expand_front',
    2: 'cpp_sliding_window',
    3: 'cpp_const_width',
    4: 'cpp_openmp',
}


def _window_mode_label(sim_config):
    """Return the descriptive window-mode name for the given sim_config dict."""
    cfg = sim_config or {}
    # Prefer the explicit 'mode' key set by the notebook (e.g. 'py_segments').
    if 'mode' in cfg:
        return cfg['mode']
    if 'solver_method' not in cfg:
        return 'py_segments'
    method = cfg['solver_method']
    mode = method.get('window_mode', 0) if isinstance(method, dict) else 0
    return _WINDOW_MODE_NAMES.get(int(mode), f'cpp_mode{mode}')


# ── Run-directory management ──────────────────────────────────────────────────

def create_run_directory(output_dir, sim_config=None, tags=None):
    """
    Create and return a unique output directory:
      <output_dir>/<timestamp>_<mode>[_tag1_tag2...]/

    Falls back to timestamp-only if mode cannot be determined.
    tags : list of str, optional – appended to the directory name with '_' separators.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    mode_name = _window_mode_label(sim_config)
    run_name = f"{ts}_{mode_name}"
    if tags:
        run_name += '_' + '_'.join(str(t).strip() for t in tags)

    run_dir = output_dir / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    print(f"✓ Run directory: {run_dir}")
    return run_dir


# ── Main entry point ──────────────────────────────────────────────────────────

def plot_results(sim, results, output_dir=None, save_plots=True, sim_config=None,
                 label="", plot_config=None, wall_time=None, tags=None):
    """
    Generate all figures and save provenance.

    Parameters
    ----------
    sim        : ClusterDynamicsSimulation
    results    : dict from post_process.calculate_derived_quantities
    output_dir : str or Path  (auto-creates run sub-directory)
    save_plots : bool
    sim_config : dict, optional  (echoed into provenance)
    label      : str  (tag for plot titles, e.g. "Python LSODA" or "C++ CVODE BDF")
    plot_config : dict, optional  (axis limits per figure; see module docstring)
    wall_time  : float, optional  (solver wall-clock time in seconds)
    """
    viz = CDVisualizer(sim.input_data, plot_config=plot_config)

    if save_plots and output_dir is not None:
        run_dir = create_run_directory(output_dir, sim_config=sim_config, tags=tags)
    else:
        run_dir = None

    print("Generating figures…")
    viz.plot_point_defects(results, run_dir, label=label)
    viz.plot_interstitial_clusters(results, run_dir, label=label)
    viz.plot_vacancy_clusters(results, run_dir, label=label)
    viz.plot_i_cluster_size(results, run_dir, label=label)
    viz.plot_v_cluster_size(results, run_dir, label=label)
    viz.plot_i_cluster_size_radius(results, run_dir, label=label)

    if run_dir is not None:
        viz.save_provenance(results, run_dir, sim_config=sim_config,
                            label=label, wall_time=wall_time)

    return run_dir


# ── Helpers ───────────────────────────────────────────────────────────────────

def _selected_sizes(max_size):
    """
    Return cluster sizes to show in the time-evolution plots:
      2, 6, 10, 20, 30, 40, … (steps of 10) up to max_size.
    """
    sizes = []
    for s in [2, 6]:
        if s <= max_size:
            sizes.append(s)
    x = 10
    while x <= max_size:
        sizes.append(x)
        x += 10
    return sizes


def _size_regimes(max_size, n_regimes=3, n_per_regime=10):
    """
    Split [2, max_size] log-uniformly into n_regimes bands and pick
    ~n_per_regime log-spaced integer sizes per band.

    Returns a list of lists, e.g. [[2,4,8,16,32], [45,100,…], […]].
    """
    log_lo = np.log10(2)
    log_hi = np.log10(max(max_size, 3))
    edges  = np.linspace(log_lo, log_hi, n_regimes + 1)

    regimes = []
    for i in range(n_regimes):
        lo = max(2, int(np.round(10 ** edges[i])))
        hi = max(lo + 1, int(np.round(10 ** edges[i + 1])))
        raw   = np.geomspace(lo, hi, n_per_regime)
        sizes = sorted(set(int(np.round(s)) for s in raw))
        sizes = [s for s in sizes if 2 <= s <= max_size]
        if sizes:
            regimes.append(sizes)
    return regimes


def _time_regimes(t_arr, n_regimes=3, n_per_regime=5):
    """
    Split [t_arr[0], t_arr[-1]] log-uniformly into n_regimes bands and pick
    ~n_per_regime log-spaced time points per band.

    Returns a list of lists of indices into t_arr.
    """
    log_lo = np.log10(t_arr[0])
    log_hi = np.log10(t_arr[-1])
    edges  = np.linspace(log_lo, log_hi, n_regimes + 1)

    regimes = []
    for i in range(n_regimes):
        # mid and late regimes start 2 decades below their respective edge minimum
        if i == 0:
            t_lo = 10 ** edges[i]
        else:
            t_lo = max(t_arr[0], 10 ** (edges[i] - 2))
        t_targets = np.geomspace(t_lo, 10 ** edges[i + 1], n_per_regime)
        indices   = _nearest_indices(t_arr, t_targets)
        if indices:
            regimes.append(indices)
    return regimes


def _nearest_indices(t_arr, t_targets):
    """
    For each target time, return the index of the nearest available time point.
    Skips targets outside [t_arr[0], t_arr[-1]].
    """
    indices = []
    seen    = set()
    for t_tgt in t_targets:
        if t_tgt < t_arr[0] or t_tgt > t_arr[-1]:
            continue
        idx = int(np.argmin(np.abs(t_arr - t_tgt)))
        if idx not in seen:
            indices.append(idx)
            seen.add(idx)
    return indices


# ── CDVisualizer class ────────────────────────────────────────────────────────

class CDVisualizer:
    """All plotting methods for the ClusterDynamics simulation."""

    def __init__(self, input_data, plot_config=None):
        self.inp         = input_data
        self.plot_config = plot_config or {}
        plt.rcParams.update({'figure.dpi': 150, 'font.size': 11})

        # Atomic number density [m⁻³] for unit conversion (BCC, 2 atoms/cell)
        a_m = input_data.derived['a']                      # m (already SI)
        self._N_at = 2.0 / a_m**3                          # atoms m⁻³

    def _apply_lims(self, ax, plot_name):
        """Apply xlim/ylim from plot_config if present."""
        cfg  = self.plot_config.get(plot_name, {})
        xlim = cfg.get('xlim', None)
        ylim = cfg.get('ylim', None)
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)

    # ── Figure 1: point defects ───────────────────────────────────────────────

    def plot_point_defects(self, results, run_dir=None, label=""):
        """Cv1 (vacancy) and Ci1 (interstitial) concentrations vs time."""
        t    = results['time']
        conc = results['concentrations']

        fig, ax = plt.subplots(figsize=(8, 5.5))
        ax.loglog(t, np.maximum(conc['Cv1'], 1e-100), 'b-',  lw=2.5,
                  label=r'$C_v$  (vacancy)')
        ax.loglog(t, np.maximum(conc['Ci1'], 1e-100), 'r--', lw=2.5,
                  label=r'$C_i$  (interstitial)')

        ax.set_xlabel('Time  (s)')
        ax.set_ylabel('Concentration  (at/at)')
        ax.set_title(self._title('Point defects  $C_v$ & $C_i$', label))
        ax.legend(fontsize=10, loc='best')
        ax.grid(True, which='both', alpha=0.3)
        self._apply_lims(ax, 'point_defects')
        plt.tight_layout()
        self._save_or_show(fig, run_dir, 'point_defects.png')

    # ── Figures 2–4: interstitial clusters vs time (small / mid / large) ────────

    def plot_interstitial_clusters(self, results, run_dir=None, label=""):
        """Three log-log time-evolution plots for small / mid / large interstitial clusters."""
        t    = results['time']
        conc = results['concentrations']
        Ni   = results['Ni']

        cfg     = self.plot_config.get('interstitial_clusters', {})
        xlim    = cfg.get('xlim', None)
        C_floor = cfg.get('C_floor', 1e-20)
        n_curves = cfg.get('n_curves', 10)

        _xmin = {'small': None, 'mid': 1e2, 'large': 1e3}
        for regime_sizes, name in zip(_size_regimes(Ni, n_per_regime=n_curves), ('small', 'mid', 'large')):
            regime_xlim = xlim.get(name, None) if isinstance(xlim, dict) else xlim
            self._plot_cluster_regime(
                t, conc,
                prefix       = 'Ci',
                sub          = 'i',
                sizes        = regime_sizes,
                title        = self._title(
                    f'Interstitial clusters – {name}  (Ni={Ni})', label),
                xlim         = regime_xlim,
                xmin_default = _xmin[name],
                C_floor      = C_floor,
                cmap         = cm.plasma,
                run_dir      = run_dir,
                filename     = f'interstitial_clusters_{name}.png',
            )

    # ── Figures 5–7: vacancy clusters vs time (small / mid / large) ──────────

    def plot_vacancy_clusters(self, results, run_dir=None, label=""):
        """Three log-log time-evolution plots for small / mid / large vacancy clusters."""
        t    = results['time']
        conc = results['concentrations']
        Nv   = results['Nv']

        cfg     = self.plot_config.get('vacancy_clusters', {})
        xlim    = cfg.get('xlim', None)
        C_floor = cfg.get('C_floor', 1e-20)
        n_curves = cfg.get('n_curves', 10)

        _xmin = {'small': None, 'mid': 1e2, 'large': 1e3}
        for regime_sizes, name in zip(_size_regimes(Nv, n_per_regime=n_curves), ('small', 'mid', 'large')):
            regime_xlim = xlim.get(name, None) if isinstance(xlim, dict) else xlim
            self._plot_cluster_regime(
                t, conc,
                prefix       = 'Cv',
                sub          = 'v',
                sizes        = regime_sizes,
                title        = self._title(
                    f'Vacancy clusters – {name}  (Nv={Nv})', label),
                xlim         = regime_xlim,
                xmin_default = _xmin[name],
                C_floor      = C_floor,
                cmap         = cm.winter,
                run_dir      = run_dir,
                filename     = f'vacancy_clusters_{name}.png',
            )

    # ── Private: one size-regime time-evolution plot ──────────────────────────

    def _plot_cluster_regime(self, t, conc, prefix, sub, sizes, title,
                             xlim, C_floor, cmap, run_dir, filename,
                             xmin_default=None):
        """
        Log-log plot of concentration vs time for a list of cluster sizes.

        y-min  = C_floor (hard floor, default 1e-20)
        y-max  = data maximum across all plotted sizes, rounded up to the
                 next log decade (with one decade of headroom).
        """
        valid = [s for s in sizes if f'{prefix}{s}' in conc]
        if not valid:
            return

        colors = cmap(np.linspace(0.1, 0.9, len(valid)))

        data_max = C_floor
        fig, ax  = plt.subplots(figsize=(9, 6))
        for sz, col in zip(valid, colors):
            vals = np.asarray(conc[f'{prefix}{sz}'])
            data_max = max(data_max, float(vals.max()))
            ax.loglog(t, np.maximum(vals, C_floor),
                      lw=1.8, color=col, label=f'$C_{{{sz}{sub}}}$')

        # y-max: round up to next decade above the data maximum
        ymax = 10 ** np.ceil(np.log10(max(data_max, C_floor * 10)))

        ax.set_ylim(C_floor, ymax)
        if xlim is not None:
            ax.set_xlim(xlim)
        elif xmin_default is not None:
            ax.set_xlim(left=xmin_default)
        ax.set_xlabel('Time  (s)')
        ax.set_ylabel('Concentration  (at/at)')
        ax.set_title(title)
        ax.legend(ncol=2, fontsize=9, loc='best')
        ax.grid(True, which='both', alpha=0.3)
        plt.tight_layout()
        self._save_or_show(fig, run_dir, filename)

    # ── Figures 8–10: interstitial size distribution (early / mid / late) ───────

    def plot_i_cluster_size(self, results, run_dir=None, label=""):
        """
        Three log-x / linear-y size-distribution plots (density vs cluster size) for
        early / mid / late time regimes, each with ~10 representative snapshots.
        """
        t    = results['time']
        conc = results['concentrations']
        Ni   = results['Ni']

        sizes   = np.arange(1, Ni + 1)
        density = np.array([conc[f'Ci{x}'] for x in sizes]) * self._N_at  # (Ni, n_t)

        cfg  = self.plot_config.get('i_cluster_size', {})
        xlim = cfg.get('xlim', None)
        ylim = cfg.get('ylim', None)
        n_snapshots = cfg.get('n_snapshots', 10)

        for snap_indices, name in zip(_time_regimes(t, n_per_regime=n_snapshots), ('early', 'mid', 'late')):
            self._plot_size_dist_regime(
                sizes, density, t, snap_indices,
                xlabel   = 'Cluster size  (number of interstitials)',
                title    = self._title(
                    f'Interstitial size distribution – {name}  (Ni={Ni})', label),
                xlim     = xlim,
                ylim     = ylim,
                cmap     = cm.plasma,
                run_dir  = run_dir,
                filename = f'i_cluster_size_{name}.png',
            )

    # ── Figures 11–13: vacancy size distribution (early / mid / late) ─────────

    def plot_v_cluster_size(self, results, run_dir=None, label=""):
        """
        Three log-x / linear-y size-distribution plots (density vs cluster size) for
        early / mid / late time regimes, each with ~10 representative snapshots.
        """
        t    = results['time']
        conc = results['concentrations']
        Nv   = results['Nv']

        sizes   = np.arange(1, Nv + 1)
        density = np.array([conc[f'Cv{x}'] for x in sizes]) * self._N_at  # (Nv, n_t)

        cfg  = self.plot_config.get('v_cluster_size', {})
        xlim = cfg.get('xlim', None)
        ylim = cfg.get('ylim', None)
        n_snapshots = cfg.get('n_snapshots', 10)

        for snap_indices, name in zip(_time_regimes(t, n_per_regime=n_snapshots), ('early', 'mid', 'late')):
            self._plot_size_dist_regime(
                sizes, density, t, snap_indices,
                xlabel   = 'Cluster size  (number of vacancies)',
                title    = self._title(
                    f'Vacancy size distribution – {name}  (Nv={Nv})', label),
                xlim     = xlim,
                ylim     = ylim,
                cmap     = cm.winter,
                run_dir  = run_dir,
                filename = f'v_cluster_size_{name}.png',
            )

    # ── Figures 14–16: interstitial size distribution vs radius (small / mid / large) ──

    def plot_i_cluster_size_radius(self, results, run_dir=None, label=""):
        """
        Three linear-x size-distribution plots (density vs cluster radius R in nm) for
        small / mid / large size regimes, with representative time snapshots overlaid.

        Radius is computed from BCC alpha-iron atomic volume (spherical cluster):
            Omega_Fe = (2.87 Å)^3 / 2
            R [nm]   = (3 n Omega_Fe / 4 pi)^(1/3) × 10^9
        """
        t    = results['time']
        conc = results['concentrations']
        Ni   = results['Ni']

        # BCC alpha-iron constants
        a_Fe_bcc = 2.87e-10          # m  (lattice parameter)
        Omega_Fe = a_Fe_bcc**3 / 2.0  # m³ per atom

        sizes    = np.arange(1, Ni + 1)
        radii_nm = (3.0 * sizes * Omega_Fe / (4.0 * np.pi))**(1.0 / 3.0) * 1e9

        density = np.array([conc[f'Ci{x}'] for x in sizes]) * self._N_at  # (Ni, n_t)

        cfg  = self.plot_config.get('i_cluster_radius', {})
        ylim = cfg.get('ylim', None)
        n_snapshots = cfg.get('n_snapshots', 10)

        # Pick n_snapshots log-spaced snapshots from the full time range
        t_targets = np.geomspace(t[0], t[-1], n_snapshots)
        all_snap_indices = _nearest_indices(t, t_targets)

        # Compute the same 3 size-regime edges used by plot_interstitial_clusters
        log_lo = np.log10(2)
        log_hi = np.log10(max(Ni, 3))
        edges  = np.linspace(log_lo, log_hi, 4)  # 3 bands → 4 edges

        for i, name in enumerate(('small', 'mid', 'large')):
            s_lo = max(1, int(np.round(10 ** edges[i])))
            s_hi = max(s_lo + 1, int(np.round(10 ** edges[i + 1])))
            mask = (sizes >= s_lo) & (sizes <= s_hi)

            self._plot_size_dist_linear_regime(
                radii_nm[mask], density[mask, :],
                t, all_snap_indices,
                xlabel   = r'Radius  $R$  (nm)',
                title    = self._title(
                    f'Interstitial size distribution – {name}  (radius, Ni={Ni})', label),
                ylim     = ylim,
                cmap     = cm.plasma,
                run_dir  = run_dir,
                filename = f'i_cluster_radius_{name}.png',
            )

    # ── Private: one size-regime radius plot (linear x-axis) ─────────────────

    def _plot_size_dist_linear_regime(self, radii_nm, density, t, snap_indices,
                                      xlabel, title, ylim, cmap, run_dir, filename):
        """
        Linear-x density vs cluster-radius plot for a set of time snapshots.

        Spline interpolation is done in log-y / linear-x space so the density
        (which spans many decades) interpolates smoothly on a linear radius axis.
        """
        colors = cmap(np.linspace(0.05, 0.95, len(snap_indices)))

        fig, ax = plt.subplots(figsize=(10, 6))
        for color, idx in zip(colors, snap_indices):
            d    = density[:, idx]
            mask = d > 0
            if not mask.any():
                continue
            xs = radii_nm[mask]
            ys = d[mask]
            # Sort by radius (should already be sorted, but be safe)
            order = np.argsort(xs)
            xs, ys = xs[order], ys[order]

            if len(xs) >= 4:
                # Interpolate in log-y / linear-x space for smooth curves
                log_ys      = np.log10(np.maximum(ys, 1e-300))
                fine_xs     = np.linspace(xs[0], xs[-1], 10 * len(xs))
                spl         = make_interp_spline(xs, log_ys, k=3)
                fine_ys     = np.maximum(10 ** spl(fine_xs), 0)
                ax.plot(fine_xs, fine_ys, lw=1.5, color=color,
                        label=self._time_label(t[idx]))
            else:
                ax.plot(xs, ys, lw=1.5, color=color,
                        label=self._time_label(t[idx]))

        if ylim is not None:
            ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'Cluster density  (m$^{-3}$)')
        ax.set_title(title)
        ax.legend(ncol=2, fontsize=8, loc='best')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        self._save_or_show(fig, run_dir, filename)

    # ── Private: one time-regime size-distribution plot ───────────────────────

    def _plot_size_dist_regime(self, sizes, density, t, snap_indices,
                               xlabel, title, xlim, ylim, cmap, run_dir, filename):
        """
        Log-log density vs cluster-size plot for a set of time snapshots.

        xlim / ylim : from plot_config (None = matplotlib auto)
        """
        colors = cmap(np.linspace(0.05, 0.95, len(snap_indices)))

        fig, ax = plt.subplots(figsize=(10, 6))
        for color, idx in zip(colors, snap_indices):
            d    = density[:, idx]
            mask = d > 0
            if mask.any():
                xs = sizes[mask].astype(float)
                ys = d[mask]
                # Interpolate in log-log space to a fine grid for smooth curves
                log_xs = np.log10(xs)
                log_ys = np.log10(np.maximum(ys, 1e-300))
                if len(xs) >= 4:
                    fine_log_xs = np.linspace(log_xs[0], log_xs[-1], 10 * len(xs))
                    spl = make_interp_spline(log_xs, log_ys, k=3)
                    fine_ys = 10 ** spl(fine_log_xs)
                    fine_ys = np.maximum(fine_ys, 0)
                    ax.semilogx(10 ** fine_log_xs, fine_ys, lw=1.5, color=color,
                                label=self._time_label(t[idx]))
                else:
                    ax.semilogx(xs, ys, lw=1.5, color=color,
                                label=self._time_label(t[idx]))

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'Cluster density  (m$^{-3}$)')
        ax.set_title(title)
        ax.legend(ncol=2, fontsize=8, loc='best')
        ax.grid(True, which='both', alpha=0.3)
        plt.tight_layout()
        self._save_or_show(fig, run_dir, filename)

    # ── Provenance ─────────────────────────────────────────────────────────────

    def save_provenance(self, results, run_dir, sim_config=None,
                        label="", wall_time=None):
        """Write a Markdown provenance file with all parameters and solver stats."""
        run_dir = Path(run_dir)
        p = self.inp.material_params
        d = self.inp.derived
        m = self.inp.model_params
        meta  = results.get('metadata', {})
        stats = meta.get('solver_stats', {})

        lines = [
            "# ClusterDynamics Provenance\n",
            f"Run directory: `{run_dir.name}`\n",
            "\n## Run Info\n",
            "| Field | Value |\n|-------|-------|\n",
            f"| solver | {label if label else 'unknown'} |\n",
        ]
        if wall_time is not None:
            lines.append(f"| wall_time_s | {wall_time:.3f} |\n")
        if sim_config is not None:
            for k, v in sim_config.items():
                lines.append(f"| {k} | {v} |\n")

        lines += [
            "\n## Material Parameters\n",
            "| Parameter | Value |\n|-----------|-------|\n",
        ]
        for k, v in p.items():
            lines.append(f"| {k} | {v} |\n")
        lines += [
            "\n## Derived Parameters\n",
            "| Parameter | Value |\n|-----------|-------|\n",
        ]
        for k, v in d.items():
            val_str = f"{v:.4e}" if isinstance(v, float) else str(v)
            lines.append(f"| {k} | {val_str} |\n")
        lines += [
            "\n## Model Parameters\n",
            "| Parameter | Value |\n|-----------|-------|\n",
        ]
        for k, v in m.items():
            lines.append(f"| {k} | {v} |\n")
        lines += [
            "\n## Solver Statistics\n",
            "| Stat | Value |\n|------|-------|\n",
        ]
        for k, v in stats.items():
            lines.append(f"| {k} | {v} |\n")

        # Build filename: <YYYYMMDD_HHMMSS>_<window_mode_name>.md
        ts_part = run_dir.name[:15]  # first 15 chars = YYYYMMDD_HHMMSS
        mode_name = _window_mode_label(sim_config)
        prov_filename = f"{ts_part}_{mode_name}.md"
        prov_path = run_dir / prov_filename
        prov_path.write_text(''.join(lines), encoding='utf-8')
        print(f"  ✓ {prov_filename}")

    # ── Internal helpers ────────────────────────────────────────────────────────

    def _title(self, base, label=""):
        T = self.inp.material_params['T']
        P = self.inp.material_params['P']
        suffix = f'  [{label}]' if label else ''
        return f'316 SS  T={T - 273.15:.0f} °C  P={P:.0e} dpa/s  —  {base}{suffix}'

    @staticmethod
    def _time_label(t_val):
        """Human-readable time label for legend entries."""
        return f'{t_val:.2e} s'

    @staticmethod
    def _save_or_show(fig, run_dir, filename):
        if run_dir is not None:
            path = Path(run_dir) / filename
            fig.savefig(path, dpi=150, bbox_inches='tight')
            print(f"  ✓ {filename}")
            plt.close(fig)
        else:
            plt.show()
