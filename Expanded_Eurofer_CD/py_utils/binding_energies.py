"""
binding_energies.py — Defect cluster binding energies for bcc Fe / EUROFER97.

All equations and tables cite:
  Ghoniem, N.M. (2026), "A Cluster Dynamics Model for Radiation Damage
  Evolution in Ferritic-Martensitic Steels" (Rate_Equations.pdf).

Conventions
-----------
- Vacancy cluster of size m (pure void):  c_{m,0}
- He-vacancy cluster:  c_{m,ℓ}  — m vacancies + ℓ He atoms
- Interstitial cluster of size n:  c_n
- All energies in eV; lengths in m unless stated.
"""

import numpy as np

_kB     = 8.617333262e-5    # eV K^-1
_J_eV   = 6.241509074e18    # J → eV


# ── Helper ────────────────────────────────────────────────────────────────────

def atomic_radius(Omega):
    """Wigner-Seitz radius r_0 = (3Ω/4π)^{1/3} [m]."""
    return (3.0 * Omega / (4.0 * np.pi)) ** (1.0 / 3.0)


def capture_radius(n, Omega):
    """r_n = r_0 · n^{1/3} [m].  Eq. 86 (equivalent form)."""
    return atomic_radius(Omega) * np.asarray(n, dtype=float) ** (1.0 / 3.0)


# ── 1.  Vacancy cluster (void) binding energies ───────────────────────────────

# Atomistic fitting amplitudes A(m) for the void capillary correction.
# Table 18 of Rate_Equations.pdf.  Index = He content ℓ = 0..4 (pure void case).
_A_void = {0: 1.2353, 1: 2.9064, 2: 3.4147, 3: 2.1504, 4: -0.1590}  # eV
_lambda_void = 0.5756   # decay constant [vac^-1], = ln(100)/8


def E_b_void(m, E_f_v, gamma_s, Omega):
    """
    Binding energy of the m-th vacancy in a spherical void cluster.

    Capillary formula with atomistic correction (Eqs. 66-67):
      E_b(m) = E_f_v − A_cap · [m^{2/3} − (m−1)^{2/3}]
               + A_void(0) · exp(−λ·(m−1))     [Eq. 67 correction]

    where  A_cap = 4π·γ_s·r_0²  [eV],  r_0 = (3Ω/4π)^{1/3}.

    Special case m=1: return E_f_v (monomer formation energy).

    Parameters
    ----------
    m       : int or array
    E_f_v   : float [eV]
    gamma_s : float [J/m²]
    Omega   : float [m³]

    Returns
    -------
    E_b : float or array [eV]
    """
    r0     = atomic_radius(Omega)
    A_cap  = 4.0 * np.pi * gamma_s * r0**2 * _J_eV   # eV  — Eq. 66
    A_atm  = _A_void[0]

    m = np.asarray(m, dtype=float)
    capillary = A_cap * (m**(2.0/3.0) - (m - 1.0)**(2.0/3.0))
    correction = A_atm * np.exp(-_lambda_void * np.maximum(m - 1.0, 0.0))

    Eb = np.where(m <= 1.0,
                  E_f_v,
                  E_f_v - capillary + correction)
    return float(Eb) if Eb.ndim == 0 else Eb


# ── 2.  He–vacancy (bubble) binding energies ─────────────────────────────────

# He virial EOS coefficients (Table 8, Eqs. 64-65)
_B2 = 1.67e-29    # m^3/atom  — 2nd virial coefficient
_B3 = 1.84e-58    # m^6/atom^2 — 3rd virial coefficient


def He_pressure(m, ell, Omega, T):
    """
    He pressure in a bubble with m vacancies and ℓ He atoms.
    Virial EOS (Eq. 64):
      P·V = N_He·k_B·T · [1 + B2·(N_He/V) + B3·(N_He/V)²]

    V = m·Ω  (bubble volume).

    Returns
    -------
    P : float [Pa]
    """
    if ell <= 0:
        return 0.0
    V  = float(m) * Omega
    rho = float(ell) / V        # He number density [m^-3]
    kB_SI = _kB / _J_eV         # eV/K → J/K
    P = float(ell) * kB_SI * T / V * (1.0 + _B2 * rho + _B3 * rho**2)
    return P


def E_b_bubble(m, ell, E_f_v, gamma_s, Omega, T):
    """
    Effective vacancy binding energy for a He-vacancy bubble (m, ℓ).

    Capillary model with He-pressure correction (Eqs. 70-73):
      E_b_bub(m, ℓ) = E_f_v − A_cap·[m^{2/3}−(m−1)^{2/3}]
                       − P_He(m,ℓ)·Ω  [He pressure work term]

    Parameters
    ----------
    m       : int
    ell     : int or float  (can be fractional mean loading)
    E_f_v   : float [eV]
    gamma_s : float [J/m²]
    Omega   : float [m³]
    T       : float [K]

    Returns
    -------
    E_b : float [eV]
    """
    r0    = atomic_radius(Omega)
    A_cap = 4.0 * np.pi * gamma_s * r0**2 * _J_eV
    m_f   = float(m)
    capillary = A_cap * (m_f**(2.0/3.0) - max(m_f - 1.0, 0.0)**(2.0/3.0))
    P     = He_pressure(m, ell, Omega, T)
    dE_P  = P * Omega * _J_eV   # Pa·m³ → eV
    return E_f_v - capillary - dE_P


# ── 3.  He binding energy to bubbles ─────────────────────────────────────────

# Atomistic fitting amplitudes A^He(α) for ℓ ≤ α.
# Table 19 of Rate_Equations.pdf.
_A_He = {2: 0.55, 3: 0.40, 4: 0.75}   # eV
_mu_He = 0.658                          # decay constant [He^-1]

def E_b_He(m, ell, E_s_He, Omega, T, B2=_B2, B3=_B3):
    """
    He binding energy: energy to remove one He from bubble (m, ℓ).

    Combined continuum + atomistic formula (Eqs. 76-77, Table 19):

      E_b^He(m, ℓ) = E_s^He + P_He(m,ℓ)·Ω   [continuum, Eq. 76]
                     + f_blend · A^He(α) · exp(−μ·(ℓ−α))  [atomistic, Eq. 77]

    where α = ℓ_max(m) = floor(alpha_He · m^{2/3}) is the trap-mutation limit,
    and f_blend fades the atomistic correction for large ℓ.

    Parameters
    ----------
    m      : int
    ell    : int  — current He content
    E_s_He : float [eV] — He interstitial solution energy (= 2.35 eV, Table 5)
    Omega  : float [m³]
    T      : float [K]

    Returns
    -------
    E_b : float [eV]
    """
    P = He_pressure(m, ell, Omega, T)
    E_cont = E_s_He + P * Omega * _J_eV   # continuum term

    # atomistic correction for small ℓ
    alpha = int(round(ell))   # closest integer loading
    A_atm = _A_He.get(alpha, 0.0)
    if A_atm != 0.0:
        E_atm = A_atm * np.exp(-_mu_He * max(ell - alpha, 0.0))
    else:
        E_atm = 0.0

    return E_cont + E_atm


def E_b_He_first(E_b_hV_1=2.30):
    """He binding energy to a monovacancy (first He, Table 5). Eq. 76."""
    return E_b_hV_1


def E_b_He_second(E_b_hV_2=2.00):
    """He binding energy to a monovacancy (second He, Table 5). Eq. 76."""
    return E_b_hV_2


# ── 4.  Interstitial loop binding energies ────────────────────────────────────

def E_b_loop_i(n, A_111=0.7501, B_111=0.3873, A_100=0.7160, B_100=0.3581,
               n_tr=25.0, sigma_tr=5.0):
    """
    Binding energy of the n-th SIA to an interstitial loop of size n−1.

    Power-law fit blended to continuum limit (Eqs. 106-108):
      E_b_111(n) = A_111 · n^{−B_111}     [Eq. 106]
      E_b_100(n) = A_100 · n^{−B_100}     [Eq. 107]
      blend      = 0.5·(1 + tanh((n − n_tr)/σ_tr))
      E_b(n)     = (1−blend)·E_b_111 + blend·E_b_100   [Eq. 108]

    Parameters
    ----------
    n        : int or array  (≥ 2)
    A_111, B_111 : float  — ½⟨111⟩ loop fit parameters (Table 18)
    A_100, B_100 : float  — ⟨100⟩ loop fit parameters  (Table 18)
    n_tr     : float      — blend center
    sigma_tr : float      — blend width

    Returns
    -------
    E_b : float or array [eV]
    """
    n = np.asarray(n, dtype=float)
    n_safe = np.maximum(n, 2.0)

    Eb_111 = A_111 * n_safe**(-B_111)
    Eb_100 = A_100 * n_safe**(-B_100)
    blend  = 0.5 * (1.0 + np.tanh((n_safe - n_tr) / sigma_tr))
    Eb     = (1.0 - blend) * Eb_111 + blend * Eb_100

    return float(Eb) if Eb.ndim == 0 else Eb


def E_b_loop_v(m, E_f_v, gamma_sf, Omega, b_111):
    """
    Vacancy loop binding energy (Eqs. 104-105).

    Continuum (Frank loop) with DFT anchor at m=2:
      E_b^vloop(m) = E_f_v − A_sf · [m^{1/2} − (m−1)^{1/2}]

    where A_sf = 2·b_111 · sqrt(π·γ_sf·Ω/√3) · _J_eV.

    Parameters
    ----------
    m       : int or array (≥ 1)
    E_f_v   : float [eV]
    gamma_sf: float [J/m²]  — stacking fault energy
    Omega   : float [m³]
    b_111   : float [m]     — ½⟨111⟩ Burgers vector

    Returns
    -------
    E_b : float or array [eV]
    """
    A_sf = 2.0 * b_111 * np.sqrt(np.pi * gamma_sf * Omega / np.sqrt(3.0)) * _J_eV
    m = np.asarray(m, dtype=float)
    Eb = np.where(
        m <= 1.0,
        E_f_v,
        E_f_v - A_sf * (m**(1.0/2.0) - (m - 1.0)**(1.0/2.0))
    )
    return float(Eb) if Eb.ndim == 0 else Eb


# ── 5.  Trap mutation barrier ─────────────────────────────────────────────────

# Table 27 of Rate_Equations.pdf
# Key: (m, ℓ)  — (vacancies, He atoms)
_E_TM_table = {
    (1, 5): 1.00,   # He5V1 — Morishita 2003
    (1, 6): 0.50,   # He6V1 — Morishita 2003
    (1, 7): 0.00,   # He7V1 — spontaneous (Gao 2011)
    (0, 4): 0.30,   # He4   self-trapping (Gao 2011)
    (0, 5): 0.10,   # He5   self-trapping (Gao 2011)
    (0, 6): 0.00,   # He6   spontaneous (Gao 2011)
}


def E_TM(m, ell, default=1e10):
    """
    Trap mutation barrier (Table 27, Eq. 83).

    Returns the barrier E_TM(m, ℓ) [eV].
    If (m, ℓ) not in the table, returns `default` (effectively disables TM).

    Parameters
    ----------
    m   : int  — vacancies
    ell : int  — He atoms
    default : float [eV]  — barrier for unknown (m, ℓ)
    """
    return _E_TM_table.get((int(m), int(ell)), default)


def Gamma_TM(m, ell, T, nu0_TM=1e12):
    """
    Trap mutation rate from bubble (m, ℓ) → (m+1, ℓ-1) [s^-1].

    Γ_TM(m, ℓ) = ν_0_TM · exp(−E_TM(m,ℓ) / k_B T)   (Eq. 142)

    Parameters
    ----------
    m, ell : int
    T      : float [K]
    nu0_TM : float [s^-1]  — attempt frequency (Table 27)
    """
    E = E_TM(m, ell)
    if E >= 1e9:
        return 0.0
    return nu0_TM * np.exp(-E / (_kB * T))


# ── 6.  Radiation re-solution rate ───────────────────────────────────────────

def Gamma_res(ell, phi_dot, b0):
    """
    Radiation re-solution rate for He from bubble (m, ℓ) [s^-1].

    Γ_res(m, ℓ) = b_0 · ℓ · φ̇   (Eq. 143)

    Parameters
    ----------
    ell     : int or float  — He content
    phi_dot : float [dpa/s] — displacement rate
    b0      : float [dpa^-1 He^-1]  — re-solution parameter (Table 29)
    """
    return b0 * float(ell) * phi_dot


# ── 7.  Maximum He loading per vacancy cluster ────────────────────────────────

def ell_max(m, alpha_He=1.7):
    """
    Maximum He loading for cluster with m vacancies.
    Trap-mutation limit: ℓ_max = floor(alpha_He · m^{2/3}).  Table 29.
    """
    return int(np.floor(alpha_He * float(m)**(2.0 / 3.0)))
