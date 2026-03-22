# -*- coding: utf-8 -*-
"""
post_process.py — ClusterDynamics post-processing

Derives macroscopic microstructural quantities from the raw cluster-size
distribution stored in results dict.
"""

import numpy as np


def mean_cluster_size(f_dict):
    """
    Compute mean cluster size <n>(t) from a size distribution dict.

    Parameters
    ----------
    f_dict : dict  {n: array(t)}  — cluster number densities

    Returns
    -------
    mean_n : array(t)
    """
    n_vals = sorted(f_dict.keys())
    num = sum(n * f_dict[n] for n in n_vals)
    den = sum(    f_dict[n] for n in n_vals)
    den = np.where(den > 0, den, 1.0)
    return num / den


def total_cluster_density(f_dict, Omega=1.0):
    """
    Total number density [m^-3] = sum_n f(n) / Omega.

    Omega = 1 returns density in at/at units (same as f_dict values).
    """
    return sum(f_dict[n] for n in f_dict) / Omega


def size_distribution_at_dose(results, dpa_target):
    """
    Extract the vacancy and interstitial size distributions at a target dose.

    Returns
    -------
    fv_snap : dict {n: float}
    fi_snap : dict {n: float}
    dpa_actual : float
    """
    dpa = results['dpa']
    idx = np.argmin(np.abs(dpa - dpa_target))
    fv_snap = {n: results['fv'][n][idx] for n in results['fv']}
    fi_snap = {n: results['fi'][n][idx] for n in results['fi']}
    return fv_snap, fi_snap, dpa[idx]


def calculate_swelling(results, Omega):
    """
    Void swelling fraction ΔV/V = (4/3)π r³ C_void / 1 [dimensionless].

    Parameters
    ----------
    Omega : float  atomic volume [m^3]
    """
    r   = results['r_void']
    N   = results['C_void'] / Omega   # number density [m^-3]
    return (4/3) * np.pi * r**3 * N


def loop_number_density(results, Omega):
    """Return loop number densities [m^-3] for both loop types."""
    return (results['CiL_111'] / Omega,
            results['CiL_100'] / Omega)


def summary_table(results, Omega):
    """
    Assemble a summary dict of key macroscopic quantities vs dose.
    """
    import pandas as pd

    N_v = len(results['fv'])
    N_i = len(results['fi'])

    nd_111, nd_100 = loop_number_density(results, Omega)
    swelling = calculate_swelling(results, Omega)
    mean_fv  = mean_cluster_size(results['fv'])
    mean_fi  = mean_cluster_size(results['fi'])

    return pd.DataFrame({
        'dpa':             results['dpa'],
        'fv1':             results['fv'][1],
        'fi1':             results['fi'][1],
        'mean_vac_cluster_size': mean_fv,
        'mean_int_cluster_size': mean_fi,
        'N_iL_111_m3':     nd_111,
        'N_iL_100_m3':     nd_100,
        'r_iL_111_nm':     results['r_iL_111'] * 1e9,
        'r_iL_100_nm':     results['r_iL_100'] * 1e9,
        'C_void_m3':       results['C_void'] / Omega,
        'r_void_nm':       results['r_void'] * 1e9,
        'swelling':        swelling,
    })
