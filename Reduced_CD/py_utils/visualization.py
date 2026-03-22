# -*- coding: utf-8 -*-
"""
visualization.py — ClusterDynamics plotting utilities

All functions accept a results dict (from simulation.py) and optional axes.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from py_utils.post_process import (
    mean_cluster_size, total_cluster_density,
    calculate_swelling, loop_number_density
)


def plot_monomer_concentrations(results, ax=None, x_axis='dpa'):
    """Mobile point defect concentrations vs dose or time."""
    x, xlabel = _x_axis(results, x_axis)
    if ax is None:
        _, ax = plt.subplots()
    ax.loglog(x, results['fv'][1], 'b-',  label='Vacancy ($f_v^1$)',       lw=2)
    ax.loglog(x, results['fi'][1], 'r-',  label='Interstitial ($f_i^1$)',   lw=2)
    ax.set_xlabel(xlabel); ax.set_ylabel('Concentration (at/at)')
    ax.set_title('Mobile Point Defects'); ax.legend(); ax.grid(True, alpha=0.3)
    return ax


def plot_cluster_distributions(results, ax=None, x_axis='dpa'):
    """All cluster sizes vs dose on one log–log panel."""
    x, xlabel = _x_axis(results, x_axis)
    if ax is None:
        _, ax = plt.subplots()
    cmap_v = plt.cm.Blues
    cmap_i = plt.cm.Reds
    N_v = len(results['fv'])
    N_i = len(results['fi'])
    for k, n in enumerate(range(1, N_v + 1)):
        c = cmap_v(0.4 + 0.6 * k / N_v)
        ax.loglog(x, results['fv'][n], color=c, lw=1.5, label=f'$f_v^{{{n}}}$')
    for k, n in enumerate(range(1, N_i + 1)):
        c = cmap_i(0.4 + 0.6 * k / N_i)
        ax.loglog(x, results['fi'][n], color=c, lw=1.5, linestyle='--', label=f'$f_i^{{{n}}}$')
    ax.set_xlabel(xlabel); ax.set_ylabel('Concentration (at/at)')
    ax.set_title('Cluster Size Distribution'); ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
    return ax


def plot_mean_cluster_size(results, ax=None, x_axis='dpa'):
    x, xlabel = _x_axis(results, x_axis)
    if ax is None:
        _, ax = plt.subplots()
    ax.semilogx(x, mean_cluster_size(results['fv']), 'b-',  label='Vacancy', lw=2)
    ax.semilogx(x, mean_cluster_size(results['fi']), 'r--', label='Interstitial', lw=2)
    ax.set_xlabel(xlabel); ax.set_ylabel('Mean cluster size ⟨n⟩')
    ax.set_title('Mean Cluster Size'); ax.legend(); ax.grid(True, alpha=0.3)
    return ax


def plot_loop_evolution(results, Omega, ax_density=None, ax_size=None, x_axis='dpa'):
    """Loop number density and mean radius vs dose."""
    x, xlabel = _x_axis(results, x_axis)
    nd_111, nd_100 = loop_number_density(results, Omega)

    if ax_density is None:
        _, ax_density = plt.subplots()
    ax_density.loglog(x, nd_111, 'b-',  label='1/2⟨111⟩', lw=2)
    ax_density.loglog(x, nd_100, 'r-',  label='⟨100⟩',    lw=2)
    ax_density.set_xlabel(xlabel); ax_density.set_ylabel('Number density (m$^{-3}$)')
    ax_density.set_title('Loop Density'); ax_density.legend(); ax_density.grid(True, alpha=0.3)

    if ax_size is None:
        _, ax_size = plt.subplots()
    ax_size.semilogx(x, results['r_iL_111'] * 1e9, 'b-',  label='1/2⟨111⟩', lw=2)
    ax_size.semilogx(x, results['r_iL_100'] * 1e9, 'r--', label='⟨100⟩',    lw=2)
    ax_size.set_xlabel(xlabel); ax_size.set_ylabel('Loop radius (nm)')
    ax_size.set_title('Loop Size'); ax_size.legend(); ax_size.grid(True, alpha=0.3)

    return ax_density, ax_size


def plot_void_evolution(results, Omega, ax_density=None, ax_size=None, x_axis='dpa'):
    x, xlabel = _x_axis(results, x_axis)
    nd_void = results['C_void'] / Omega

    if ax_density is None:
        _, ax_density = plt.subplots()
    ax_density.loglog(x, nd_void, 'b-', label='Void', lw=2)
    ax_density.set_xlabel(xlabel); ax_density.set_ylabel('Number density (m$^{-3}$)')
    ax_density.set_title('Void Density'); ax_density.legend(); ax_density.grid(True, alpha=0.3)

    if ax_size is None:
        _, ax_size = plt.subplots()
    ax_size.semilogx(x, results['r_void'] * 1e9, 'b-', label='Void', lw=2)
    ax_size.set_xlabel(xlabel); ax_size.set_ylabel('Void radius (nm)')
    ax_size.set_title('Void Size'); ax_size.legend(); ax_size.grid(True, alpha=0.3)

    return ax_density, ax_size


def plot_swelling(results, Omega, ax=None, x_axis='dpa'):
    x, xlabel = _x_axis(results, x_axis)
    if ax is None:
        _, ax = plt.subplots()
    ax.semilogx(x, calculate_swelling(results, Omega) * 100, 'b-', lw=2)
    ax.set_xlabel(xlabel); ax.set_ylabel('Swelling ΔV/V (%)')
    ax.set_title('Void Swelling'); ax.grid(True, alpha=0.3)
    return ax


def plot_all(results, input_data, output_dir=None, x_axis='dpa'):
    """Generate and save the full 3×3 summary figure."""
    Omega = input_data.physical_props['Omega']
    T     = int(input_data.material_params['T'] - 273)
    G     = input_data.material_params['G']

    fig, axes = plt.subplots(3, 3, figsize=(18, 12))
    fig.suptitle(
        f"ClusterDynamics — EUROFER97  T={T} °C  G={G:.1e} dpa/s",
        fontsize=14, fontweight='bold'
    )

    plot_monomer_concentrations(results, ax=axes[0, 0], x_axis=x_axis)
    plot_cluster_distributions( results, ax=axes[0, 1], x_axis=x_axis)
    plot_mean_cluster_size(     results, ax=axes[0, 2], x_axis=x_axis)

    plot_loop_evolution(results, Omega,
                        ax_density=axes[1, 0], ax_size=axes[1, 1],
                        x_axis=x_axis)
    plot_void_evolution(results, Omega,
                        ax_density=axes[1, 2], ax_size=axes[2, 0],
                        x_axis=x_axis)
    plot_swelling(results, Omega, ax=axes[2, 1], x_axis=x_axis)

    # Hide unused panel
    axes[2, 2].set_visible(False)

    plt.tight_layout()

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        tag = 'Ion' if G > 1e-5 else 'Neutron'
        fig.savefig(output_dir / f'ClusterDynamics_{T}C_{tag}.png', dpi=300, bbox_inches='tight')
        print(f"Figure saved → {output_dir}")

    plt.show()
    return fig


# ------------------------------------------------------------------
def _x_axis(results, mode):
    if mode == 'dpa':
        return results['dpa'], 'Dose (dpa)'
    return results['time'], 'Time (s)'
