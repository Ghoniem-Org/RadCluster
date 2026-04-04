"""
visualization.py — Plotting utilities for Eurofer_CD simulation results.

Produces a standard set of figures from the results dict returned by
post_process.calculate_derived_quantities().

Adapted from Full_CD/py_utils/visualization.py for the Eurofer_CD state vector
(SIA + vacancy + He species).
"""

import datetime
import subprocess
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    'font.size':        11,
    'axes.labelsize':   12,
    'axes.titlesize':   12,
    'legend.fontsize':  9,
    'lines.linewidth':  1.5,
})

OUTPUT_DIR = Path(__file__).parent.parent / 'output'


# ── Directory helpers ─────────────────────────────────────────────────────────

def create_run_directory(label='run'):
    """Create a timestamped output subdirectory and return its path."""
    ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    try:
        sha = subprocess.check_output(
            ['git', 'rev-parse', '--short', 'HEAD'],
            cwd=str(Path(__file__).parent.parent),
            stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        sha = 'nogit'
    run_dir = OUTPUT_DIR / f'{ts}_{label}_{sha}'
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / 'plots').mkdir(exist_ok=True)
    return run_dir


def write_provenance(run_dir, results, input_data):
    """Write a provenance.md file documenting the run parameters."""
    md = run_dir / 'provenance.md'
    p  = input_data.material_params
    d  = input_data.derived

    lines = [
        '# Eurofer_CD Run Provenance',
        '',
        f'- **Date**: {datetime.datetime.now().isoformat()}',
        f'- **T**: {float(p["T"])} K',
        f'- **G (dose rate)**: {float(p["G"]):.2e} dpa/s',
        f'- **rho_d**: {float(p.get("rho_d", 0)):.2e} m⁻²',
        f'- **he_mode**: {d["he_mode"]}',
        f'- **Ni**: {input_data.Ni}  **Nv**: {input_data.Nv}',
        f'- **Di**: {d["Di"]:.3e} m²/s  **Dv**: {d["Dv"]:.3e} m²/s',
        f'- **Cv_eq**: {d["Cv_eq"]:.3e}',
    ]
    if 'metadata' in results:
        meta = results['metadata'].get('solver_stats', {})
        lines += [
            '',
            '## Solver',
            f'- {meta.get("message", "")}',
            f'- wall time: {meta.get("wall_time", 0):.2f} s',
        ]
    md.write_text('\n'.join(lines), encoding='utf-8')
    print(f"Provenance written: {md}")


# ── Main visualization class ───────────────────────────────────────────────────

class EuroferCDVisualizer:
    """
    Generates the standard figure set for a Eurofer_CD run.

    Parameters
    ----------
    results    : dict  from post_process.calculate_derived_quantities
    input_data : InputData
    run_dir    : Path   output directory (plots/ subdirectory used)
    """

    def __init__(self, results, input_data, run_dir=None):
        self.res     = results
        self.inp     = input_data
        self.run_dir = Path(run_dir) if run_dir else OUTPUT_DIR / 'last_run'
        self.run_dir.mkdir(parents=True, exist_ok=True)
        (self.run_dir / 'plots').mkdir(exist_ok=True)

        self.t      = results['time']
        self.names  = list(results['concentrations'].keys())
        re_Ni = results['Ni']
        re_Nv = results['Nv']
        self.Ni = re_Ni
        self.Nv = re_Nv

    def plot_all(self):
        """Generate and save all standard figures."""
        self._plot_point_defects()
        self._plot_totals()
        self._plot_sia_clusters()
        self._plot_vac_clusters()
        self._plot_mean_sizes()
        self._plot_swelling()
        self._plot_he_content()
        self._plot_size_distributions()
        print(f"All figures saved to: {self.run_dir / 'plots'}")

    def _savefig(self, name):
        path = self.run_dir / 'plots' / f'{name}.png'
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()

    def _plot_point_defects(self):
        """Free point defect concentrations: Ci1, Cv1, C_He."""
        fig, ax = plt.subplots(figsize=(7, 4.5))
        c = self.res['concentrations']

        if 'Ci1' in c:
            ax.loglog(self.t, np.maximum(c['Ci1'], 1e-40), label='$C_{i1}$ (SIA)')
        if 'Cv1' in c:
            ax.loglog(self.t, np.maximum(c['Cv1'], 1e-40), label='$C_{v1}$ (vacancy)')
        if 'C_He' in c:
            ax.loglog(self.t, np.maximum(c['C_He'], 1e-40), label='$C_{He}$ (free He)', ls='--')

        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Concentration (atom fraction)')
        ax.set_title('Point defect concentrations')
        ax.set_ylim(bottom=1e-18)
        ax.legend()
        ax.grid(True, which='both', alpha=0.3)
        plt.tight_layout()
        self._savefig('point_defects')

    def _plot_totals(self):
        """Total interstitial and vacancy content vs time."""
        fig, ax = plt.subplots(figsize=(7, 4.5))
        tot = self.res['totals']

        ax.loglog(self.t, np.maximum(tot['total_i'], 1e-40), label='Total SIA content')
        ax.loglog(self.t, np.maximum(tot['total_v'], 1e-40), label='Total vacancy content', ls='--')

        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Total defect content (atom fraction)')
        ax.set_title('Total defect inventories')
        ax.set_ylim(bottom=1e-18)
        ax.legend()
        ax.grid(True, which='both', alpha=0.3)
        plt.tight_layout()
        self._savefig('totals')

    def _plot_he_content(self):
        """Free He and trapped He estimates."""
        fig, ax = plt.subplots(figsize=(7, 4.5))
        c   = self.res['concentrations']
        tot = self.res['totals']

        ax.loglog(self.t, np.maximum(tot['C_He'], 1e-40), label='Free He')

        ax.set_xlabel('Time (s)')
        ax.set_ylabel('He concentration (atom fraction)')
        ax.set_title('Free He concentration')
        ax.set_ylim(bottom=1e-18)
        ax.legend()
        ax.grid(True, which='both', alpha=0.3)
        plt.tight_layout()
        self._savefig('he_content')

    def _plot_sia_clusters(self):
        """SIA cluster concentrations vs time (small / mid / large)."""
        c = self.res['concentrations']
        groups = [
            ('SIA clusters: n = 1–5',    range(1, min(6,  self.Ni+1)),  'sia_small'),
            ('SIA clusters: n = 6–20',   range(6, min(21, self.Ni+1)),  'sia_mid'),
            ('SIA clusters: n = 20–100', range(20, min(101,self.Ni+1)), 'sia_large'),
        ]
        for title, rng, fname in groups:
            keys = [f'Ci{n}' for n in rng if f'Ci{n}' in c]
            if not keys:
                continue
            fig, ax = plt.subplots(figsize=(7, 4.5))
            cmap = plt.cm.viridis(np.linspace(0, 1, len(keys)))
            for key, color in zip(keys, cmap):
                n = int(key[2:])
                ax.loglog(self.t, np.maximum(c[key], 1e-40), color=color, label=f'n={n}')
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Concentration (atom fraction)')
            ax.set_title(title)
            ax.set_ylim(bottom=1e-18)
            if len(keys) <= 10:
                ax.legend(fontsize=7, ncol=2)
            ax.grid(True, which='both', alpha=0.3)
            plt.tight_layout()
            self._savefig(fname)

    def _plot_vac_clusters(self):
        """Vacancy cluster concentrations vs time."""
        c = self.res['concentrations']
        groups = [
            ('Vacancy clusters: m = 1–5',   range(1, min(6,  self.Nv+1)), 'vac_small'),
            ('Vacancy clusters: m = 6–20',  range(6, min(21, self.Nv+1)), 'vac_mid'),
            ('Vacancy clusters: m = 20–100',range(20, min(101,self.Nv+1)),'vac_large'),
        ]
        for title, rng, fname in groups:
            keys = [f'Cv{m}' for m in rng if f'Cv{m}' in c]
            if not keys:
                continue
            fig, ax = plt.subplots(figsize=(7, 4.5))
            cmap = plt.cm.plasma(np.linspace(0, 1, len(keys)))
            for key, color in zip(keys, cmap):
                m = int(key[2:])
                ax.loglog(self.t, np.maximum(c[key], 1e-40), color=color, label=f'm={m}')
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Concentration (atom fraction)')
            ax.set_title(title)
            ax.set_ylim(bottom=1e-18)
            if len(keys) <= 10:
                ax.legend(fontsize=7, ncol=2)
            ax.grid(True, which='both', alpha=0.3)
            plt.tight_layout()
            self._savefig(fname)

    def _plot_mean_sizes(self):
        """Mean cluster sizes vs time."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
        ms = self.res['mean_sizes']

        ax1.semilogx(self.t, ms['mean_n'], color='C0')
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Mean cluster size')
        ax1.set_title('Mean SIA cluster size ⟨n⟩')
        ax1.grid(True, which='both', alpha=0.3)

        ax2.semilogx(self.t, ms['mean_m'], color='C1')
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('Mean cluster size')
        ax2.set_title('Mean vacancy cluster size ⟨m⟩')
        ax2.grid(True, which='both', alpha=0.3)

        plt.tight_layout()
        self._savefig('mean_sizes')

    def _plot_swelling(self):
        """Void swelling (volumetric) vs time."""
        fig, ax = plt.subplots(figsize=(7, 4.5))
        sw = self.res['swelling']
        # Convert to % swelling
        ax.semilogx(self.t, sw * 100.0, color='C3')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Swelling (%)')
        ax.set_title('Void swelling estimate')
        ax.grid(True, which='both', alpha=0.3)
        plt.tight_layout()
        self._savefig('swelling')

    def _plot_size_distributions(self):
        """Snapshot of cluster size distributions at final time."""
        c   = self.res['concentrations']
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

        # SIA distribution at final time
        ns  = np.arange(1, self.Ni + 1)
        Ci_final = np.array([max(float(c.get(f'Ci{n}', np.array([1e-100]))[-1]), 1e-100)
                             for n in ns])
        axes[0].semilogy(ns, Ci_final, 'o-', ms=3, color='C0')
        axes[0].set_xlabel('Cluster size n')
        axes[0].set_ylabel('Concentration (atom fraction)')
        axes[0].set_title(f'SIA cluster distribution at t = {self.t[-1]:.2e} s')
        axes[0].set_ylim(bottom=1e-18)
        axes[0].grid(True, alpha=0.3)

        # Vacancy distribution at final time
        ms_arr = np.arange(1, self.Nv + 1)
        Cv_final = np.array([max(float(c.get(f'Cv{m}', np.array([1e-100]))[-1]), 1e-100)
                             for m in ms_arr])
        axes[1].semilogy(ms_arr, Cv_final, 'o-', ms=3, color='C1')
        axes[1].set_xlabel('Cluster size m')
        axes[1].set_ylabel('Concentration (atom fraction)')
        axes[1].set_title(f'Vacancy cluster distribution at t = {self.t[-1]:.2e} s')
        axes[1].set_ylim(bottom=1e-18)
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        self._savefig('size_distributions')
