"""
reaction_rates.py – Pre-computed rate constant arrays for cluster dynamics.

Implements the capture and emission rate formulas from Ghoniem & Cho (1979)
as used in CD.ipynb.  All arrays are computed once at initialisation and
stored for fast use inside the ODE right-hand side.

Rate constant naming mirrors CD.ipynb:
  KCV[x-1] = Kc_v(x)  — spherical capture of vacancies by a size-x vacancy cluster
  KCI[x-1] = Kc_i(x)  — spherical capture of interstitials by a size-x vacancy cluster
  KLV[x-1] = Kl_v(x)  — circular-loop capture of vacancies by a size-x int. cluster
  KLI[x-1] = Kl_i(x)  — circular-loop capture of interstitials by a size-x int. cluster
  GCV[x-1] = gamma_cv(x) — thermal emission from a size-x vacancy cluster
  GLV[x-1] = gamma_lv(x) — thermal emission from a size-x interstitial loop
"""

import numpy as np


class ReactionRates:
    """
    Pre-computed rate constant arrays for the Ghoniem & Cho ODE system.

    All arrays are 0-indexed; index k corresponds to cluster size k+1.
    """

    def __init__(self, input_data):
        self.inp = input_data
        self._precompute_rate_constants()

    # ── Pre-computation ──────────────────────────────────────────────────────

    def _precompute_rate_constants(self):
        p   = self.inp.material_params
        d   = self.inp.derived
        Nv  = self.inp.model_params['Nv']
        Ni  = self.inp.model_params['Ni']
        kBT = d['kB'] * p['T']

        # Cache frequently used quantities
        nu_v = p['nu_v']; nu_i = p['nu_i']
        Z_v  = p['Z_v'];  Z_i  = p['Z_i']
        exp_v = np.exp(-p['E_m_v'] / kBT)
        exp_i = np.exp(-p['E_m_i'] / kBT)
        a    = d['a']
        g    = p['g']

        # ── Spherical-capture rates (vacancy clusters) ───────────────────────
        # Kc_v(x) — vacancy capture by size-x spherical cluster
        def Kc_v(x):
            return (2.216 * Z_v**2 * x**(2/3) / (1 + 0.1128 * Z_v * x**(1/3))
                    * nu_v * exp_v)

        # Kc_i(x) — interstitial capture by size-x spherical cluster (annihilation)
        def Kc_i(x):
            return (2.216 * Z_i**2 * x**(2/3) / (1 + 0.1128 * Z_i * x**(1/3))
                    * nu_i * exp_i)

        # ── Circular-loop capture rates (interstitial loops) ─────────────────
        # Kl_v(x) — vacancy capture by size-x interstitial loop
        def Kl_v(x):
            return 1.555 * Z_v * x**0.5 * nu_v * exp_v

        # Kl_i(x) — interstitial capture by size-x interstitial loop
        def Kl_i(x):
            return 1.555 * Z_i * x**0.5 * nu_i * exp_i

        # ── Thermal emission rates ────────────────────────────────────────────
        # gamma_cv(x) — vacancy emission from size-x spherical cluster
        def gamma_cv(x):
            return Kc_v(x) * d['Cv_eq'] * np.exp(
                (1.28 * g * a**2 / kBT) * x**(-1/3))

        # gamma_lv(x) — vacancy emission from size-x interstitial loop
        # Non-zero only for x=2 (diinterstitial dissociation)
        def gamma_lv(x):
            return Kl_v(2) * np.exp(-p['E_b_2i'] / kBT) if int(x) == 2 else 0.0

        # ── Build arrays ──────────────────────────────────────────────────────
        xs_v = np.arange(1, Nv + 1, dtype=float)
        xs_i = np.arange(1, Ni + 1, dtype=float)

        self.KCV = np.array([Kc_v(x)    for x in xs_v])
        self.KCI = np.array([Kc_i(x)    for x in xs_v])
        self.KLV = np.array([Kl_v(x)    for x in xs_i])
        self.KLI = np.array([Kl_i(x)    for x in xs_i])
        self.GCV = np.array([gamma_cv(x) for x in xs_v])
        self.GLV = np.array([gamma_lv(x) for x in xs_i])

        print(f"ReactionRates: KCV[0]={self.KCV[0]:.3e}  KLI[0]={self.KLI[0]:.3e}"
              f"  GCV[0]={self.GCV[0]:.3e}")
