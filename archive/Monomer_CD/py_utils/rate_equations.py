"""
rate_equations.py – ODE right-hand side for the Ghoniem & Cho cluster dynamics system.

Faithfully translates rhs_full() from CD.ipynb into the ZrMicro module pattern.

State vector y[N],  N = Nv + Ni:
  y[0 .. Nv-1]      — vacancy cluster concentrations Cv1, Cv2, ..., Cv_Nv
  y[Nv .. Nv+Ni-1]  — interstitial cluster concentrations Ci1, Ci2, ..., Ci_Ni

References
----------
Ghoniem & Cho, J. Nucl. Mater. 84 (1979) 202–215.
"""

import numpy as np


class RateEquations:
    """
    150-equation ODE system (Nv=50 vacancy + Ni=100 interstitial cluster sizes).

    Mirrors the structure of ZrMicro's RateEquations class so that
    simulation.py / cpp_bridge.py / post_process.py work identically.
    """

    def __init__(self, input_data, reaction_rates):
        self.inp = input_data
        self.rr  = reaction_rates

        Nv = input_data.model_params['Nv']
        Ni = input_data.model_params['Ni']
        self.Nv = Nv
        self.Ni = Ni
        self.N  = Nv + Ni

        # Names used by post_process / visualization
        self.concentration_names = (
            [f'Cv{x}' for x in range(1, Nv + 1)] +
            [f'Ci{x}' for x in range(1, Ni + 1)]
        )

        print(f"RateEquations: N={self.N} equations  (Nv={Nv}, Ni={Ni})")

    # ── Public interface (matches ZrMicro convention) ────────────────────────

    def ode_system(self, t, y):
        """Full ODE right-hand side — thin wrapper around _rhs_full."""
        return self._rhs_full(y)

    def gated_ode_system(self, t, y, x_max):
        """
        ODE right-hand side with x_max gating.

        Freezes dC/dt = 0 for all interstitial cluster sizes > x_max.
        Used by the segmented Python solver (simulation.py).
        """
        dydt = self._rhs_full(y)
        dydt[self.Nv + x_max:] = 0.0
        return dydt

    def get_initial_conditions(self):
        """Build the initial condition vector matching CD.ipynb."""
        Nv = self.Nv; Ni = self.Ni
        d  = self.inp.derived
        y0 = np.zeros(self.N)
        y0[0]      = d['Cv_eq']            # Cv1 = thermal equilibrium
        y0[1]      = 1e-6  * d['Cv_eq']   # Cv2
        y0[2]      = 1e-8  * d['Cv_eq']   # Cv3
        y0[Nv]     = 1e-20                 # Ci1
        y0[Nv + 1] = 1e-30                 # Ci2
        y0[Nv + 2] = 1e-35                 # Ci3
        return y0

    # ── ODE right-hand side ──────────────────────────────────────────────────

    def _rhs_full(self, y):
        """
        Full ODE right-hand side — direct translation of rhs_full() in CD.ipynb.

        All physics from Ghoniem & Cho (1979), Equations 6–8.
        """
        Nv  = self.Nv; Ni = self.Ni; N = self.N
        p   = self.inp.material_params
        d   = self.inp.derived
        rr  = self.rr

        # Floor concentrations (guard against solver taking negative steps)
        Cv_arr = np.maximum(y[:Nv], 1e-100)
        Ci_arr = np.maximum(y[Nv:], 1e-100)
        dydt   = np.zeros(N)

        Cv  = Cv_arr[0]  # single vacancy
        Ci  = Ci_arr[0]  # single interstitial
        C2v = Cv_arr[1]  # divacancy
        C2i = Ci_arr[1]  # diinterstitial

        # ── dCv1/dt ──────────────────────────────────────────────────────────
        # Sources: production, annihilation of C2v by interstitials, emission from C2v
        # Sinks:   self-clustering, capture by all clusters, dislocation sinks, recombination
        dCv = p['P'] + rr.KCI[1]*Ci*C2v + (2*rr.GCV[1] - rr.KCV[1]*Cv)*C2v
        dCv += np.dot(rr.GCV[2:Nv], Cv_arr[2:Nv])           # emission from Cv3..CvNv
        dCv -= np.dot(rr.KCV[2:Nv], Cv_arr[2:Nv]) * Cv       # capture by  Cv3..CvNv
        dCv -= np.dot(rr.KLV[2:Ni], Ci_arr[2:Ni]) * Cv       # capture by  Ci3..CiNi
        dCv += p['Z_v'] * p['rho_d'] * d['Dv'] * (d['Cv_eq'] - Cv)  # dislocation sink
        dCv -= d['alpha'] * Cv * Ci                            # recombination
        dCv -= rr.KCV[0] * Cv**2                               # Cv + Cv → C2v
        dCv -= rr.KLV[1] * Cv * C2i                           # C2i absorbs Cv → Ci
        dydt[0] = dCv

        # ── dCv2/dt ──────────────────────────────────────────────────────────
        dC2v = (0.5 * rr.KCV[0] * Cv**2             # Cv + Cv → C2v
                + rr.GCV[2] * Cv_arr[2]              # emission: Cv3 → C2v + Cv
                + rr.KCI[2] * Ci * Cv_arr[2]         # Ci annihilates one vacancy in Cv3
                + p['rho_d'] * d['D2v'] * (d['C2v_eq'] - C2v))  # dislocation sink
        dC2v -= (rr.KCV[1]*Cv + rr.KCI[1]*Ci + rr.GCV[1]) * C2v
        dydt[1] = dC2v

        # ── dCvx/dt  x = 3..Nv ───────────────────────────────────────────────
        for x in range(3, Nv + 1):
            i    = x - 1
            Cxm1 = Cv_arr[x - 2]   # Cv(x-1)
            Cx   = Cv_arr[x - 1]   # Cvx
            dCx  = (rr.KCV[x-2]*Cv*Cxm1           # Cv(x-1) absorbs Cv → Cvx
                    - rr.KCI[x-1]*Ci*Cx             # Cvx annihilated by Ci → Cv(x-1)
                    - rr.KCV[x-1]*Cv*Cx             # Cvx absorbs Cv → Cv(x+1)
                    - rr.GCV[x-1]*Cx)               # Cvx emits Cv → Cv(x-1)
            if x < Nv:
                dCx += (rr.KCI[x]*Ci*Cv_arr[x]     # Cv(x+1) annihilated by Ci → Cvx
                        + rr.GCV[x]*Cv_arr[x])      # Cv(x+1) emits Cv → Cvx
            dydt[i] = dCx

        # ── dCi1/dt ──────────────────────────────────────────────────────────
        dCi = (p['P']
               + rr.KLV[1]*Cv*C2i                   # C2i absorbs Cv → Ci
               - d['K_nuc_i']*Ci**2                 # nucleation: Ci + Ci → C2i
               - d['alpha']*Cv*Ci                   # recombination with vacancies
               - rr.KLI[1]*Ci*C2i                   # C2i captures Ci → C3i
               - rr.KCI[1]*Ci*C2v)                  # C2v captures Ci
        dCi -= np.dot(rr.KLI[2:Ni], Ci_arr[2:Ni]) * Ci   # capture by Ci3..CiNi
        dCi -= np.dot(rr.KCV[2:Nv], Cv_arr[2:Nv]) * Ci   # capture by vacancy clusters
        dCi -= p['Z_i'] * p['rho_d'] * d['Di'] * Ci       # dislocation sink
        dydt[Nv] = dCi

        # ── dCi2/dt ──────────────────────────────────────────────────────────
        # GLV[1] = gamma_lv(2): thermal SIA emission from C2i → 2×Ci1 (dissociation)
        dC2i = (0.5 * d['K_nuc_i'] * Ci**2          # Ci + Ci → C2i
                + rr.KLV[2]*Cv*Ci_arr[2]            # C3i absorbs Cv → C2i
                - (rr.KLI[1]*Ci + rr.KLV[1]*Cv + rr.GLV[1]) * C2i)  # losses from C2i
        dydt[Nv + 1] = dC2i

        # Thermal dissociation of C2i → 2×Ci1 feeds back into Ci1
        dCi += 2.0 * rr.GLV[1] * C2i
        dydt[Nv] = dCi

        # ── dCix/dt  x = 3..Ni ───────────────────────────────────────────────
        for x in range(3, Ni + 1):
            i    = Nv + x - 1
            Cxm1 = Ci_arr[x - 2]   # Ci(x-1)
            Cx   = Ci_arr[x - 1]   # Cix
            dCx  = (rr.KLI[x-2]*Ci*Cxm1            # Ci(x-1) captures Ci → Cix
                    - (rr.KLI[x-1]*Ci + rr.KLV[x-1]*Cv)*Cx)  # losses from Cix
            if x < Ni:
                dCx += rr.KLV[x]*Cv*Ci_arr[x]      # Ci(x+1) absorbs Cv → Cix
            dydt[i] = dCx

        return dydt
