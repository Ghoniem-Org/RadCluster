#!/usr/bin/env python3
"""
Compute defect cluster production fractions epsilon_m^(i) and epsilon_n^(v)
for fission and fusion neutron spectra in alpha-iron, based on MD cascade
simulation data from the literature.

Key sources:
- Stoller (2000): clustering fractions and size distributions
- Nordlund et al. (2018): arc-dpa survival efficiency
- De Backer et al. (2016): subcascade fragmentation and scaling
- Malerba et al. (2021): recommended modeling parameters
"""

import numpy as np

# =============================================================================
# Parameters from literature
# =============================================================================

# --- Fission spectrum (typical average PKA ~ 10-20 keV, T = 100K-600K) ---
# From Stoller (2000), Malerba et al. (2021)
fission = {
    'label': 'Fission',
    'T_K': 573,            # typical operating temperature (K)
    'E_PKA_avg_keV': 15,   # spectrum-weighted average PKA energy
    'eta': 0.30,           # defect survival efficiency (N_d/N_NRT)
    'f_i_cl': 0.58,        # fraction of surviving SIA in clusters
    'f_v_cl': 0.15,        # fraction of surviving vacancies in clusters
    's_i': 1.6,            # power-law exponent for SIA cluster number distribution
    's_v': 2.5,            # power-law exponent for vacancy cluster number distribution
    'm1': 20,              # maximum SIA cluster size
    'n1': 10,              # maximum vacancy cluster size
    'He_rate': '0.5--1',   # He production (appm He/dpa)
}

# --- Fusion spectrum (14 MeV neutrons, higher PKA energies) ---
# From Stoller (2000), De Backer et al. (2016), Nordlund et al. (2018)
fusion = {
    'label': 'Fusion',
    'T_K': 573,
    'E_PKA_avg_keV': 40,   # higher average PKA from 14 MeV neutrons
    'eta': 0.28,           # slightly lower survival at higher energies
    'f_i_cl': 0.65,        # higher clustering fraction at higher energies
    'f_v_cl': 0.20,        # more vacancy clustering at higher energies
    's_i': 1.5,            # shallower exponent -> more large clusters
    's_v': 2.3,            # shallower than fission
    'm1': 50,              # larger max SIA cluster size
    'n1': 20,              # larger max vacancy cluster size
    'He_rate': '10',       # He production (appm He/dpa)
}


def compute_epsilon(f_cl, s, m_max):
    """
    Compute the production fractions epsilon_m for cluster sizes m = 2, ..., m_max.

    The cluster number distribution follows a power law:
        N(m) ~ m^{-s}

    The fraction of all surviving defects of a given type in clusters of size m:
        m * epsilon_m

    Normalization:
        sum_{m=2}^{m_max} m * epsilon_m = f_cl

    With epsilon_m = C * m^{-s}, we get:
        C * sum_{m=2}^{m_max} m^{1-s} = f_cl
        C = f_cl / sum_{m=2}^{m_max} m^{1-s}

    Parameters
    ----------
    f_cl : float
        Total fraction of surviving defects in clusters
    s : float
        Power-law exponent for cluster number distribution
    m_max : int
        Maximum cluster size

    Returns
    -------
    m_arr : array of cluster sizes (2, 3, ..., m_max)
    eps_arr : array of epsilon_m values
    """
    m_arr = np.arange(2, m_max + 1)
    norm_sum = np.sum(m_arr ** (1.0 - s))
    C = f_cl / norm_sum
    eps_arr = C * m_arr ** (-s)
    return m_arr, eps_arr


def compute_defect_survival(E_keV):
    """
    Compute athermal recombination corrected dpa efficiency xi_arc(E)
    for iron, using the Nordlund et al. (2018) parameterization.

    xi_arc(E) = (1 - c_arc) / (2 E_d / (0.8 * E))^b_arc  (simplified form)

    For Fe at 100K from Nordlund et al. (2018):
    """
    # Nordlund et al. (2018) arc-dpa parameters for Fe
    E_d = 40.0  # displacement threshold energy (eV)
    # Simplified fit: eta ~ 0.84 * (E/keV)^{-0.20} for E > 1 keV
    # Asymptotes to ~0.30 above 10 keV
    if E_keV < 0.1:
        return 1.0
    eta = min(1.0, 0.84 * (E_keV) ** (-0.20))
    return eta


def get_cluster_spectra(params):
    """
    Return interstitial and vacancy cluster production spectra for a given spectrum dict.

    Parameters
    ----------
    params : dict
        One of the spectrum parameter dicts (e.g. `fission` or `fusion`).

    Returns
    -------
    m_arr : ndarray  — SIA cluster sizes (2 … m_max)
    eps_i : ndarray  — epsilon_m^(i) values
    n_arr : ndarray  — vacancy cluster sizes (2 … n_max)
    eps_v : ndarray  — epsilon_n^(v) values
    """
    m_arr, eps_i = compute_epsilon(params['f_i_cl'], params['s_i'], params['m1'])
    n_arr, eps_v = compute_epsilon(params['f_v_cl'], params['s_v'], params['n1'])
    return m_arr, eps_i, n_arr, eps_v
