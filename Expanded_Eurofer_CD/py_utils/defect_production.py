"""
defect_production.py — Cascade defect production fractions for bcc Fe / EUROFER97.

Implements the power-law cluster number distributions for SIA and vacancy
clusters produced in displacement cascades, for fission and fusion spectra.

Physics reference
-----------------
Ghoniem, N.M. (2026), Sections 1-2 (Rate_Equations.pdf):
  Eqs. 1-13, Tables 2 and 5.

Production rates [atom frac / s]
---------------------------------
For SIA clusters (n = 2..m1):
  P_n^i = η · G · ε_n^i          (Eq. 12)
  ε_n^i = C_i · n^{−s_i}         (Eq. 7)
  C_i   = f_i^cl / Σ_{n=2}^{m1} n^{1−s_i}   (Eq. 9)

For vacancies clusters (m = 2..n1):
  P_m^v = η · G · ε_m^v          (Eq. 13)
  ε_m^v = C_v · m^{−s_v}         (Eq. 8)
  C_v   = f_v^cl / Σ_{m=2}^{n1} m^{1−s_v}   (Eq. 10)

Point-defect production (monomer):
  P_1^i = η · G · (1 − f_i^cl)   (Eq. 11 — surviving free SIAs)
  P_1^v = η · G · (1 − f_v^cl)   (Eq. 11 — surviving free vacancies)
"""

import numpy as np


# ── Cascade parameters from Rate_Equations.pdf Tables 2 and 5 ────────────────

# Fission neutron spectrum (typical PWR/FFTF, average PKA ~10-20 keV)
FISSION = {
    'label':    'Fission',
    'eta':      0.30,    # survival efficiency η  (Table 2)
    'f_cl_i':   0.58,    # SIA clustering fraction  (Table 2)
    's_i':      1.6,     # SIA power-law exponent  (Table 2, Eq. 7)
    'm1':       20,      # max SIA cluster size  (Table 2)
    'f_cl_v':   0.15,    # vacancy clustering fraction  (Table 2)
    's_v':      2.5,     # vacancy power-law exponent  (Table 2, Eq. 8)
    'n1':       10,      # max vacancy cluster size  (Table 2)
    'G_He_r':   0.75,    # He production rate [appm He/dpa]  (Table 2)
    'b0_res':   0.01,    # re-solution parameter b_0 [dpa^-1/He]  (Table 29)
}

# Fusion neutron spectrum (14 MeV, higher PKA energies)
FUSION = {
    'label':    'Fusion',
    'eta':      0.28,    # survival efficiency η  (Table 2)
    'f_cl_i':   0.65,    # SIA clustering fraction  (Table 2)
    's_i':      1.5,     # SIA power-law exponent  (Table 2, Eq. 7)
    'm1':       50,      # max SIA cluster size  (Table 2)
    'f_cl_v':   0.20,    # vacancy clustering fraction  (Table 2)
    's_v':      2.3,     # vacancy power-law exponent  (Table 2, Eq. 8)
    'n1':       20,      # max vacancy cluster size  (Table 2)
    'G_He_r':   10.0,    # He production rate [appm He/dpa]  (Table 2)
    'b0_res':   0.10,    # re-solution parameter b_0 [dpa^-1/He]  (Table 29)
}


def normalisation_constant(f_cl, s, m_max):
    """
    Compute the normalisation constant C such that:
      Σ_{m=2}^{m_max} m · C · m^{−s} = f_cl

    C = f_cl / Σ_{m=2}^{m_max} m^{1−s}   (Eqs. 9-10)

    Parameters
    ----------
    f_cl  : float  — clustering fraction
    s     : float  — power-law exponent
    m_max : int    — maximum cluster size

    Returns
    -------
    C : float
    """
    ms = np.arange(2, m_max + 1, dtype=float)
    denom = np.sum(ms**(1.0 - s))
    return f_cl / denom if denom > 0 else 0.0


def compute_epsilon(f_cl, s, m_max):
    """
    Compute production fractions ε_m = C · m^{−s} for m = 2..m_max.

    The fraction of surviving defects in clusters of size m is m·ε_m,
    summing to f_cl over all cluster sizes.

    Parameters
    ----------
    f_cl  : float  — total clustering fraction
    s     : float  — power-law exponent
    m_max : int    — maximum cluster size

    Returns
    -------
    epsilon : ndarray [m_max+1]
        epsilon[0] = 0 (size 0 unused)
        epsilon[1] = 0 (monomer handled separately)
        epsilon[m] = C · m^{−s}  for m = 2..m_max

    Notes
    -----
    Array length is m_max+1 so that epsilon[m] = ε_m (1-indexed access).
    """
    C = normalisation_constant(f_cl, s, m_max)
    eps = np.zeros(m_max + 1)
    for m in range(2, m_max + 1):
        eps[m] = C * m**(-s)
    return eps


def production_rates(G, spectrum, N_max, M_max):
    """
    Compute cascade production rates P_n^i [s^-1] and P_m^v [s^-1] scaled
    by displacement rate G [dpa/s].

    Returns arrays indexed from 1:
      Pr_SIA[n]  = η·G·ε_n^i    for n=2..m1; P_1^i = η·G·(1−f_cl_i)
      Pr_VAC[m]  = η·G·ε_m^v    for m=2..n1; P_1^v = η·G·(1−f_cl_v)

    Parameters
    ----------
    G        : float [dpa/s]
    spectrum : str  'fission' | 'fusion'
    N_max    : int  — total SIA cluster sizes tracked
    M_max    : int  — total vacancy cluster sizes tracked

    Returns
    -------
    Pr_SIA : ndarray [N_max+1]   (index 0 unused; Pr_SIA[1]=monomer)
    Pr_VAC : ndarray [M_max+1]   (index 0 unused; Pr_VAC[1]=monomer)
    G_He   : float [atom frac/s]  — He transmutation production rate
    """
    spec = FISSION if 'fiss' in spectrum.lower() else FUSION

    eta    = spec['eta']
    f_cl_i = spec['f_cl_i']
    s_i    = spec['s_i']
    m1     = spec['m1']
    f_cl_v = spec['f_cl_v']
    s_v    = spec['s_v']
    n1     = spec['n1']
    G_He_r = spec['G_He_r']   # appm He/dpa

    # SIA clusters (Eq. 12)
    eps_i = compute_epsilon(f_cl_i, s_i, min(m1, N_max))
    Pr_SIA = np.zeros(N_max + 1)
    Pr_SIA[1] = eta * G * (1.0 - f_cl_i)            # monomer (Eq. 11)
    for n in range(2, min(m1, N_max) + 1):
        Pr_SIA[n] = eta * G * eps_i[n]               # cluster (Eq. 12)

    # Vacancy clusters (Eq. 13)
    eps_v = compute_epsilon(f_cl_v, s_v, min(n1, M_max))
    Pr_VAC = np.zeros(M_max + 1)
    Pr_VAC[1] = eta * G * (1.0 - f_cl_v)            # monomer (Eq. 11)
    for m in range(2, min(n1, M_max) + 1):
        Pr_VAC[m] = eta * G * eps_v[m]               # cluster (Eq. 13)

    # He transmutation production (Eq. 5)
    # G_He_r in appm He/dpa → [atom frac/s]: appm * 1e-6 * G
    G_He = G_He_r * 1.0e-6 * G

    return Pr_SIA, Pr_VAC, G_He
