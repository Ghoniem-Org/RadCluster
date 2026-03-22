# -*- coding: utf-8 -*-
"""
rate_equations.py — ClusterDynamics ODE system

State vector layout:
  y[0 : N_v]          → fv[1 .. N_v]   vacancy cluster size distribution
  y[N_v : N_v+N_i]    → fi[1 .. N_i]   interstitial cluster size distribution
  y[N_v+N_i]          → CiL_111        1/2<111> loop number density
  y[N_v+N_i+1]        → CiL_100        <100>   loop number density
  y[N_v+N_i+2]        → CiL_i_111      interstitials in 1/2<111> loops
  y[N_v+N_i+3]        → CiL_i_100      interstitials in <100> loops
  y[N_v+N_i+4]        → C_void         void number density
  y[N_v+N_i+5]        → r_void         mean void radius

Total length: N_v + N_i + 6
"""

import numpy as np


class ClusterDynamicsODE:
    """
    Builds and evaluates the right-hand side of the cluster-dynamics ODE system.
    """

    def __init__(self, input_data, rates):
        self.input_data = input_data
        self.rates = rates

        self.N_v    = input_data.derived['N_v']
        self.N_i    = input_data.derived['N_i']
        self.N_loop = input_data.derived['N_loop']

        # Index offsets into the state vector
        self.idx_fv      = 0
        self.idx_fi      = self.N_v
        self.idx_iL_111  = self.N_v + self.N_i
        self.idx_iL_100  = self.N_v + self.N_i + 1
        self.idx_iLi_111 = self.N_v + self.N_i + 2
        self.idx_iLi_100 = self.N_v + self.N_i + 3
        self.idx_void    = self.N_v + self.N_i + 4
        self.idx_rvoid   = self.N_v + self.N_i + 5

        self.n_equations = self.N_v + self.N_i + 6
        self._build_name_list()

        print(f"ClusterDynamicsODE initialised: {self.n_equations} equations "
              f"(N_v={self.N_v}, N_i={self.N_i}, N_loop={self.N_loop})")

    # ------------------------------------------------------------------
    def _build_name_list(self):
        self.names = (
            [f'fv_{n}' for n in range(1, self.N_v + 1)]
            + [f'fi_{n}' for n in range(1, self.N_i + 1)]
            + ['CiL_111', 'CiL_100', 'CiL_i_111', 'CiL_i_100', 'C_void', 'r_void']
        )

    # ------------------------------------------------------------------
    def _unpack(self, y):
        y = np.maximum(y, 0.0)
        fv      = y[self.idx_fv  : self.idx_fv + self.N_v]
        fi      = y[self.idx_fi  : self.idx_fi + self.N_i]
        CiL_111  = y[self.idx_iL_111]
        CiL_100  = y[self.idx_iL_100]
        CiLi_111 = y[self.idx_iLi_111]
        CiLi_100 = y[self.idx_iLi_100]
        C_void   = y[self.idx_void]
        r_void   = max(y[self.idx_rvoid], 0.5e-9)
        return fv, fi, CiL_111, CiL_100, CiLi_111, CiLi_100, C_void, r_void

    # ------------------------------------------------------------------
    def _loop_radii(self, CiL_111, CiL_100, CiLi_111, CiLi_100):
        d = self.input_data.derived
        r_111 = d['l_111'] * np.sqrt(max(CiLi_111, 1e-40) / max(CiL_111, 1e-40)) if CiL_111 > 1e-40 else 5e-9
        r_100 = d['l_100'] * np.sqrt(max(CiLi_100, 1e-40) / max(CiL_100, 1e-40)) if CiL_100 > 1e-40 else 5e-9
        return r_111, r_100

    # ------------------------------------------------------------------
    def get_initial_conditions(self):
        y0 = np.zeros(self.n_equations)
        C_v_eq = self.input_data.derived['C_v_eq']
        s1, s2 = 1e-6, 1e-10

        # Vacancy clusters
        y0[self.idx_fv]     = C_v_eq          # fv[1] = thermal equilibrium
        for n in range(2, self.N_v + 1):
            y0[self.idx_fv + n - 1] = s2 * C_v_eq

        # Interstitial clusters
        y0[self.idx_fi]     = s1 * C_v_eq     # fi[1]
        for n in range(2, self.N_i + 1):
            y0[self.idx_fi + n - 1] = s2 * C_v_eq

        y0[self.idx_iL_111]  = s2 * C_v_eq
        y0[self.idx_iL_100]  = s2 * C_v_eq
        y0[self.idx_iLi_111] = s2 * C_v_eq
        y0[self.idx_iLi_100] = s2 * C_v_eq
        y0[self.idx_void]    = s1 * C_v_eq
        y0[self.idx_rvoid]   = 2e-9

        print("Initial conditions set.")
        for i, name in enumerate(self.names):
            print(f"  {name}: {y0[i]:.2e}")
        return y0

    # ------------------------------------------------------------------
    def ode_system(self, t, y):
        """ODE right-hand side — called by scipy solve_ivp."""
        (fv, fi, CiL_111, CiL_100, CiLi_111, CiLi_100, C_void, r_void) = self._unpack(y)
        r  = self.rates
        G  = self.input_data.material_params['G']
        dpa = G * t

        # Loop radii
        r_111, r_100 = self._loop_radii(CiL_111, CiL_100, CiLi_111, CiLi_100)

        # Sink strengths
        C_v_s, C_i_s = r.sink_strengths(CiL_111, CiL_100, r_111, r_100)

        # Mobile monomer fluxes
        fv1, fi1 = fv[0], fi[0]
        phi_v = r.omega_v * fv1
        phi_i = r.omega_i * (fi1 + 2.0 * (fi[1] if self.N_i >= 2 else 0.0)
                             + 3.0 * (fi[2] if self.N_i >= 3 else 0.0))

        dydt = np.zeros(self.n_equations)

        # ----------------------------------------------------------------
        # Vacancy cluster distribution  fv[n],  n = 1 .. N_v
        # ----------------------------------------------------------------
        for n in range(1, self.N_v + 1):
            idx = self.idx_fv + n - 1
            fn  = fv[n - 1]

            # Production from cascade
            gen = r.G_v(n)

            # Absorption from left: J(n-1 → n) = beta_v(n-1)*fv[n-1]
            abs_in  = r.beta_v(n - 1) * fv[n - 2] if n > 1 else 0.0

            # Emission to left: alpha_v(n)*fv[n]
            em_out  = r.alpha_v(n) * fn

            # Absorption to right: beta_v(n)*fv[n]  (loss from size n)
            abs_out = r.beta_v(n) * fn if n < self.N_v else 0.0

            # Emission from right: alpha_v(n+1)*fv[n+1]  (gain at size n)
            em_in   = r.alpha_v(n + 1) * fv[n] if n < self.N_v else 0.0

            # i-v recombination (only affects fv[1])
            recomb = r.R_iv(fv1, fi1) if n == 1 else 0.0

            # Sink absorption (only fv[1] is mobile)
            sink = r.omega_v * fn * C_v_s if n == 1 else 0.0

            dydt[idx] = gen + abs_in - em_out - abs_out + em_in - recomb - sink

        # ----------------------------------------------------------------
        # Interstitial cluster distribution  fi[n],  n = 1 .. N_i
        # ----------------------------------------------------------------
        for n in range(1, self.N_i + 1):
            idx = self.idx_fi + n - 1
            fn  = fi[n - 1]

            gen = r.G_i(n)

            abs_in  = r.beta_i(n - 1) * fi[n - 2] if n > 1 else 0.0
            em_out  = r.alpha_i(n) * fn
            abs_out = r.beta_i(n) * fn if n < self.N_i else 0.0
            em_in   = r.alpha_i(n + 1) * fi[n] if n < self.N_i else 0.0

            # i-v recombination (only fi[1])
            recomb = r.R_iv(fv1, fi1) if n == 1 else 0.0

            # Cluster–vacancy annihilation (all mobile vacancies can absorb i-clusters)
            iv_ann = 0.0
            if n >= 2:
                iv_ann = (r.omega_v + r.omega_i) * fi[n - 1] * fv1

            # Sink (only fi[1] mobile)
            sink = r.omega_i * fn * C_i_s if n == 1 else 0.0

            # Clusters larger than N_loop feed loops (loss term)
            to_loops = 0.0
            if n == self.N_loop:
                nuc_rate = r.loop_nucleation_rate(fi, self.N_loop)
                to_loops = nuc_rate  # contributes to loop nucleation, not to fi

            dydt[idx] = gen + abs_in - em_out - abs_out + em_in - recomb - iv_ann - sink

        # ----------------------------------------------------------------
        # Loop number densities  CiL_111, CiL_100
        # ----------------------------------------------------------------
        nuc_total = r.loop_nucleation_rate(fi, self.N_loop)
        f100      = r.frac_100(dpa)
        nuc_111   = (1.0 - f100) * nuc_total
        nuc_100   = f100         * nuc_total
        trans_111 = r.loop_transition_rate(CiL_111)

        dydt[self.idx_iL_111] = nuc_111 - trans_111
        dydt[self.idx_iL_100] = nuc_100 + trans_111

        # ----------------------------------------------------------------
        # Interstitials in loops  CiL_i_111, CiL_i_100
        # ----------------------------------------------------------------
        dydt[self.idx_iLi_111] = r.loop_growth_rate('111', CiL_111, CiLi_111, phi_i, phi_v)
        dydt[self.idx_iLi_100] = r.loop_growth_rate('100', CiL_100, CiLi_100, phi_i, phi_v)

        # ----------------------------------------------------------------
        # Void number density
        # ----------------------------------------------------------------
        dydt[self.idx_void] = (r.void_nucleation_rate(C_void, r_void)
                               - r.void_dissolution_rate(C_void, r_void))

        # ----------------------------------------------------------------
        # Void radius
        # ----------------------------------------------------------------
        dydt[self.idx_rvoid] = (r.void_vacancy_absorption(fv1, r_void)
                                - r.void_interstitial_absorption(fi1, r_void)
                                - r.void_emission(r_void))

        return dydt

    # ------------------------------------------------------------------
    def loop_radii(self, y):
        """Public helper: compute loop radii from state vector."""
        (_, _, CiL_111, CiL_100, CiLi_111, CiLi_100, _, _) = self._unpack(y)
        return self._loop_radii(CiL_111, CiL_100, CiLi_111, CiLi_100)
