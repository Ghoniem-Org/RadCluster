"""
binding_energies.py – Defect cluster binding energies for bcc Fe / EUROFER97.

Implements the binding energy formulas from:
  Ghoniem, N.M. (2024), "Formulation of Cluster Dynamics Equations for
  Irradiated Ferritic-Martensitic Steels" (Rate_Equations.pdf, Sections 3–4).

Conventions
-----------
- Vacancy cluster of size m:  (−m, 0)
- Interstitial cluster of size n:  (n, 0)
- He-vacancy cluster of size (m, ℓ):  m vacancies + ℓ He atoms
- All energies in eV.

Section 4 references (PDF)
---------------------------
- Eq. 62–63:  Void / vacancy-cluster binding (capillary + atomistic correction)
- Eq. 64–65:  He binding to vacancy clusters (Caturla atomistic data + fit)
- Eq. 83–85:  Interstitial loop binding (blended power-law / continuum)
"""

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
_kB = 8.617333262e-5   # eV K^-1


# ===========================================================================
# 1.  Vacancy cluster (void embryo) binding energies
# ===========================================================================

def E_b_void(m, E_f_v, gamma_s, Omega):
    """
    Binding energy of the m-th vacancy in a spherical vacancy cluster
    (capillary approximation, PDF eq. 62).

    E_b(m) = E_f_v  −  A_void · [m^(2/3) − (m−1)^(2/3)]

    where  A_void = 4π · γ_s · r_0^2,  r_0 = (3Ω / 4π)^(1/3).

    Special case m = 1: return E_f_v (monomer formation energy; no surface
    curvature correction applies to a single vacancy).

    Parameters
    ----------
    m       : int or array   — cluster size (≥ 1)
    E_f_v   : float [eV]    — vacancy formation energy
    gamma_s : float [J/m²]  — surface / interface energy (converted internally)
    Omega   : float [m³]    — atomic volume

    Returns
    -------
    E_b : float or array [eV]
    """
    _J_to_eV = 6.241509074e18
    r0 = (3.0 * Omega / (4.0 * np.pi)) ** (1.0 / 3.0)   # m
    A_void = 4.0 * np.pi * gamma_s * r0**2 * _J_to_eV      # eV

    m = np.asarray(m, dtype=float)
    Eb = np.where(
        m <= 1.0,
        E_f_v,
        E_f_v - A_void * (m**(2.0/3.0) - (m - 1.0)**(2.0/3.0))
    )
    return float(Eb) if Eb.ndim == 0 else Eb


def C_v_eq_surf(m, E_f_v, gamma_s, Omega, T):
    """
    Thermal equilibrium vacancy concentration at the surface of a void of size m.

    Uses the capillary (Kelvin) equation:
        C_surf(m) = C_v^0_eq · exp(2·γ_s·Ω / (r_m · kT))
    where r_m = r_0 · m^(1/3).

    For m = 1 returns C_v^0_eq (no curvature correction for monomer).

    Parameters
    ----------
    m       : int or array
    E_f_v   : float [eV]
    gamma_s : float [J/m²]
    Omega   : float [m³]
    T       : float [K]

    Returns
    -------
    C_surf : float or array  (dimensionless atom fraction)
    """
    _J_to_eV = 6.241509074e18
    kBT = _kB * T
    r0 = (3.0 * Omega / (4.0 * np.pi)) ** (1.0 / 3.0)
    # 2·γ·Ω / (r0·kT)  — curvature correction exponent for m=1 sphere
    A_curv = 2.0 * gamma_s * Omega / r0 * _J_to_eV   # eV

    Cv0 = np.exp(-E_f_v / kBT)

    m = np.asarray(m, dtype=float)
    curv_factor = np.where(
        m <= 1.0,
        1.0,
        np.exp(A_curv / (m**(1.0/3.0) * kBT))
    )
    return Cv0 * curv_factor


# ===========================================================================
# 2.  Interstitial loop binding energies
# ===========================================================================

def E_b_loop(n, E_b_2i, E_b_inf, n_trans=8.0, alpha=0.6):
    """
    Binding energy of the n-th SIA to an interstitial loop of size n−1.
    Uses a blended exponential fit from small-cluster atomistic data to
    the continuum limit (PDF eqs. 83–85).

    Formula:
        E_b(n) = E_b_2i  for n = 2  (di-interstitial)
        E_b(n) = E_b_inf − (E_b_inf − E_b_2i) · exp(−(n−2) / n_trans)  for n ≥ 2

    This saturates smoothly at E_b_inf for large n (bulk stacking-fault limit).

    Parameters
    ----------
    n        : int or array  — cluster size (≥ 2)
    E_b_2i   : float [eV]   — di-interstitial binding (atomistic, ~0.8 eV in Fe)
    E_b_inf  : float [eV]   — large-loop binding limit (~1.8 eV for ½⟨111⟩ in Fe)
    n_trans  : float         — crossover cluster size
    alpha    : float         — not used in current form (retained for compatibility)

    Returns
    -------
    E_b : float or array [eV]
    """
    n = np.asarray(n, dtype=float)
    x = np.clip(n - 2.0, 0.0, None)
    Eb = E_b_inf - (E_b_inf - E_b_2i) * np.exp(-x / n_trans)
    return float(Eb) if Eb.ndim == 0 else Eb


# ===========================================================================
# 3.  He-vacancy cluster binding energies
# ===========================================================================

# Caturla et al. (2005) atomistic data for small (m, ℓ) He-vacancy clusters
# in bcc Fe.  Each entry: E_b_He[m][ell] = He binding energy (eV) for
# removing one He from a cluster with m vacancies and ell He atoms.
# Index convention: dict keyed by (m, ell).
_caturla_Eb_He = {
    (1, 1): 2.60,   # He in mono-vacancy
    (1, 2): 2.20,
    (2, 1): 2.80,
    (2, 2): 2.50,
    (2, 3): 2.25,
    (3, 1): 2.90,
    (3, 2): 2.65,
    (3, 3): 2.40,
    (3, 4): 2.20,
}

# Fitting parameters for large clusters (capillary + He pressure)
# E_b_He(m, ell) = E_b_He_v + delta_He * (ell/m)^beta_He
_He_fit = {
    'E_b_He_v':   2.60,   # eV  — He binding to a single vacancy site
    'delta_He':  -0.80,   # eV  — pressure coefficient
    'beta_He':    0.70,   # power-law exponent for He/V ratio
}


def E_b_He(m, ell, E_b_He_v=None, delta_He=None, beta_He=None):
    """
    He binding energy: energy to remove one He atom from cluster (m, ell).

    For small (m, ell) ≤ 3, uses Caturla atomistic values directly.
    For larger clusters, uses the fitted empirical model:
        E_b_He(m, ℓ) = E_b_He_v + delta_He · (ℓ/m)^beta_He

    Parameters
    ----------
    m         : int   — number of vacancies
    ell       : int   — number of He atoms
    E_b_He_v  : float [eV], optional  — override fit parameter
    delta_He  : float [eV], optional  — override fit parameter
    beta_He   : float,      optional  — override fit parameter

    Returns
    -------
    E_b : float [eV]
    """
    fit = _He_fit.copy()
    if E_b_He_v is not None:
        fit['E_b_He_v'] = E_b_He_v
    if delta_He is not None:
        fit['delta_He'] = delta_He
    if beta_He is not None:
        fit['beta_He'] = beta_He

    if (m, ell) in _caturla_Eb_He:
        return _caturla_Eb_He[(m, ell)]
    ratio = float(ell) / float(max(m, 1))
    return fit['E_b_He_v'] + fit['delta_He'] * ratio**fit['beta_He']


def dE_b_He_dell(m, ell_mean, params=None):
    """
    Partial derivative ∂E_b^He/∂ℓ evaluated at ℓ = ℓ_mean.

    Used in the 'decoupled' He-reduction mode (PDF Section 5.6.5) to compute
    the He-pressure correction to the effective void binding energy:
        E_b_eff(m) = E_b_void(m) + ℓ_mean · ∂E_b^He/∂ℓ

    For the empirical model:
        ∂E_b^He/∂ℓ ≈ delta_He · beta_He / (m · (ℓ/m)^(1−beta_He))

    Parameters
    ----------
    m        : int/float   — number of vacancies in cluster
    ell_mean : float       — mean He loading per cluster
    params   : dict or None  — override _He_fit values

    Returns
    -------
    dE : float [eV per He atom]
    """
    fit = _He_fit.copy()
    if params is not None:
        fit.update(params)

    delta = fit['delta_He']
    beta  = fit['beta_He']
    m     = float(max(m, 1))
    ell   = float(max(ell_mean, 1e-6))

    # d/d_ell [delta * (ell/m)^beta] = delta * beta / m * (ell/m)^(beta-1)
    return delta * beta / m * (ell / m)**(beta - 1.0)


# ===========================================================================
# 4.  He-vacancy cluster void emission barrier (capillary + He pressure)
# ===========================================================================

def E_b_bubble(m, ell, E_f_v, gamma_s, Omega):
    """
    Effective vacancy binding energy for a He-vacancy cluster (bubble) (m, ℓ).

    In the capillary model, internal He pressure P_He reduces the effective
    surface tension, lowering the emission barrier:
        E_b_bubble(m, ℓ) = E_b_void(m) + ΔE_He_pressure(m, ℓ)

    where ΔE_He_pressure = E_b_He(m, ℓ) contribution changes with ℓ loading.

    For now uses the simplified formula from PDF Section 5.6.5:
        E_b_eff(m, ℓ) = E_f_v − A_void · [m^(2/3) − (m−1)^(2/3)]
                         − ℓ · |dE_b_He/dℓ|

    Parameters
    ----------
    m       : int
    ell     : float  — He content (can be fractional mean loading)
    E_f_v   : float [eV]
    gamma_s : float [J/m²]
    Omega   : float [m³]

    Returns
    -------
    E_b : float [eV]
    """
    Eb_void = E_b_void(m, E_f_v, gamma_s, Omega)
    dE_dHe  = dE_b_He_dell(m, ell)
    # He pressure lowers effective emission barrier (reduces E_b_void)
    return Eb_void + ell * dE_dHe


# ===========================================================================
# 5.  Capture radius
# ===========================================================================

def capture_radius(n, Omega):
    """
    Capture radius r_n = r_0 · n^(1/3)  [m]
    where r_0 = (3Ω / 4π)^(1/3).

    From PDF eq. 86.

    Parameters
    ----------
    n     : int or array  — cluster size
    Omega : float [m³]    — atomic volume

    Returns
    -------
    r : float or array [m]
    """
    r0 = (3.0 * Omega / (4.0 * np.pi)) ** (1.0 / 3.0)
    return r0 * np.asarray(n, dtype=float)**(1.0 / 3.0)
