# -*- coding: utf-8 -*-
"""
reaction_rates.py — ClusterDynamics rate coefficients

Pre-computes all scalar rate constants (independent of the state vector)
and provides per-call functions for state-dependent rates.

Design: one instance per simulation; called from ClusterDynamicsODE.
"""

import numpy as np


class ClusterRates:
    """
    All absorption, emission, and generation rates for the cluster-size
    resolved system.

    Notation follows CLAUDE.md:
      beta_v(n)  — vacancy monomer capture rate by vacancy cluster of size n
      beta_i(n)  — interstitial monomer capture rate by interstitial cluster of size n
      alpha_v(n) — thermal emission rate from vacancy cluster of size n
      alpha_i(n) — thermal emission rate from interstitial cluster of size n
    """

    k_B = 8.617e-5  # eV/K

    def __init__(self, input_data):
        self.input_data = input_data
        self.T = input_data.material_params['T']
        self._precompute()

    # ------------------------------------------------------------------
    # Pre-computation
    # ------------------------------------------------------------------
    def _precompute(self):
        p  = self.input_data.physical_props
        mp = self.input_data.model_params
        d  = self.input_data.derived

        z_c   = p['z_c']
        nu_i  = p['nu_i']
        nu_v  = p['nu_v']
        a     = p['a']

        # --- Jump frequencies ---
        self.omega_i  = z_c * nu_i * np.exp(-p['E_m_i']  / (self.k_B * self.T))
        self.omega_2i = z_c * nu_i * np.exp(-p['E_m_2i'] / (self.k_B * self.T))
        self.omega_3i = self.omega_2i
        self.omega_v  = z_c * nu_v * np.exp(-p['E_m_v']  / (self.k_B * self.T))

        # --- Diffusion coefficients ---
        self.D_i = (a**2 / z_c) * self.omega_i
        self.D_v = (a**2 / z_c) * self.omega_v
        self.D_2i = (a**2 / z_c) * self.omega_2i

        print(f"D_v = {self.D_v:.3e}, D_i = {self.D_i:.3e} m^2/s")

        # --- Thermal emission from small clusters ---
        self.e_v  = np.exp(-p['E_F_v']  / (self.k_B * self.T))   # monomer eq. conc.
        self.e_2i = np.exp(-p['E_b_2i'] / (self.k_B * self.T))   # di-i binding
        self.e_3i = np.exp(-p['E_b_3i'] / (self.k_B * self.T))   # tri-i binding

        # --- Recombination factor ---
        self.recom = p.get('recom', 1.0)

        # --- Monomer radius ---
        self.Omega = p['Omega']
        self.r0    = d['r0']

        # --- Vacancy binding energies for clusters of size n ---
        # Simple model: E_b(n) = E_b_bulk * (1 - (n0/n)^(1/3)) + E_F_v
        # Parametrised by E_b_bulk (cohesive contribution) from Model_Parameters.
        self.E_b_bulk_v = mp.get('E_b_bulk_v', 1.6)  # eV; Fe vacancy cluster default
        self.E_b_bulk_i = mp.get('E_b_bulk_i', 3.5)  # eV; Fe SIA cluster default

        # --- Loop bias factors ---
        self.Z_iL_i = mp.get('Z_iL_i', 1.5)
        self.Z_iL_v = mp.get('Z_iL_v', 0.7)

        # --- Surface energy for void emission ---
        self.gamma = p['gamma']

    # ------------------------------------------------------------------
    # Cluster geometry
    # ------------------------------------------------------------------
    def cluster_radius(self, n):
        """Radius of a cluster of size n (spherical approximation)."""
        return self.r0 * n**(1/3)

    # ------------------------------------------------------------------
    # Binding energies  (capillary / bulk approximation)
    # ------------------------------------------------------------------
    def E_bind_v(self, n):
        """Binding energy [eV] of a vacancy to a vacancy cluster of size n."""
        if n <= 1:
            return self.input_data.physical_props['E_F_v']
        return self.E_b_bulk_v * (1.0 - (1.0 / n)**(1/3))

    def E_bind_i(self, n):
        """Binding energy [eV] of an interstitial to an interstitial cluster of size n."""
        if n == 2:
            return self.input_data.physical_props['E_b_2i']
        if n == 3:
            return self.input_data.physical_props['E_b_3i']
        return self.E_b_bulk_i * (1.0 - (1.0 / n)**(1/3))

    # ------------------------------------------------------------------
    # Absorption rate coefficients  [m^3/s] / Omega  → [at fraction / s]
    # (returns rate per unit monomer concentration, per unit cluster conc.)
    # ------------------------------------------------------------------
    def beta_v(self, n):
        """
        Absorption rate coefficient: a vacancy monomer is captured by a
        vacancy cluster of size n.
        beta_v(n) = 4*pi*r_n * D_v / Omega
        """
        r_n = self.cluster_radius(n)
        return 4 * np.pi * r_n * self.D_v / self.Omega

    def beta_i(self, n):
        """
        Absorption rate coefficient: an interstitial monomer is captured
        by an interstitial cluster of size n.
        For n=1 use D_i; for n=2 use D_2i; for n>=3 clusters are immobile
        so only the monomer diffuses.
        """
        r_n = self.cluster_radius(n)
        # Cluster mobility: only monomer diffuses for n>=3
        D_eff = self.D_i  # monomer diffusion drives capture
        return 4 * np.pi * r_n * D_eff / self.Omega

    # ------------------------------------------------------------------
    # Emission rate coefficients  [1/s]
    # ------------------------------------------------------------------
    def alpha_v(self, n):
        """
        Thermal emission rate from a vacancy cluster of size n.
        alpha_v(n) = beta_v(n-1) * exp(-E_bind_v(n) / kT)
        """
        if n <= 1:
            return 0.0
        return self.beta_v(n - 1) * np.exp(-self.E_bind_v(n) / (self.k_B * self.T))

    def alpha_i(self, n):
        """
        Thermal emission rate from an interstitial cluster of size n.
        """
        if n <= 1:
            return 0.0
        return self.beta_i(n - 1) * np.exp(-self.E_bind_i(n) / (self.k_B * self.T))

    # ------------------------------------------------------------------
    # Generation rates  [at/at/s]
    # ------------------------------------------------------------------
    def G_v(self, n):
        """Production rate into vacancy cluster of size n."""
        d = self.input_data.derived
        if n == 1:
            return d['G_v']
        elif n == 2:
            return d['G_2i'] / 2   # cascade-produced di-vacancies
        else:
            return 0.0

    def G_i(self, n):
        """Production rate into interstitial cluster of size n."""
        d = self.input_data.derived
        if n == 1:
            return d['G_i']
        elif n == 2:
            return d['G_2i'] / 2
        elif n == 3:
            return d['G_3i'] / 3
        else:
            return 0.0

    # ------------------------------------------------------------------
    # Recombination  [at/at/s]
    # ------------------------------------------------------------------
    def R_iv(self, fv1, fi1):
        """
        i–v recombination sink on monomer vacancies.
        R_iv = recom * (omega_i + omega_v) * fi1 * fv1
        """
        return self.recom * (self.omega_i + self.omega_v) * fi1 * fv1

    # ------------------------------------------------------------------
    # Sink strengths (dislocation network + loops)
    # ------------------------------------------------------------------
    def sink_strengths(self, CiL_111, CiL_100, r_iL_111, r_iL_100):
        """
        Returns (C_v_s, C_i_s) effective sink concentrations [at/at]
        for vacancy and interstitial monomers.
        """
        p   = self.input_data.physical_props
        mp  = self.input_data.model_params
        a   = p['a']
        z_c = p['z_c']
        Omega = p['Omega']
        rho_N = self.input_data.material_params.get('rho', 1e14)
        Z_N   = mp.get('Z_N', 1.05)

        rho_iL_111 = 2 * np.pi * r_iL_111 * CiL_111 / Omega
        rho_iL_100 = 2 * np.pi * r_iL_100 * CiL_100 / Omega
        rho_loops  = rho_iL_111 + rho_iL_100

        C_v_s = (a**2 / z_c) * (rho_N + rho_loops)
        C_i_s = (a**2 / z_c) * (Z_N * rho_N + Z_N * rho_loops)
        return C_v_s, C_i_s

    # ------------------------------------------------------------------
    # Void growth / emission
    # ------------------------------------------------------------------
    def void_vacancy_absorption(self, fv1, r_void):
        """dr/dt contribution from vacancy absorption."""
        a   = self.input_data.physical_props['a']
        z_c = self.input_data.physical_props['z_c']
        return (a**2 / (z_c * r_void)) * self.omega_v * fv1

    def void_interstitial_absorption(self, fi1, r_void):
        """dr/dt contribution from interstitial absorption (shrinkage)."""
        a   = self.input_data.physical_props['a']
        z_c = self.input_data.physical_props['z_c']
        return (a**2 / (z_c * r_void)) * self.omega_i * fi1

    def void_emission(self, r_void):
        """dr/dt contribution from thermal vacancy emission."""
        a     = self.input_data.physical_props['a']
        z_c   = self.input_data.physical_props['z_c']
        Omega = self.Omega
        curvature_factor = np.exp(2.0 * self.gamma * Omega / (r_void * self.k_B * self.T)) - 1.0
        return (a**2 / (z_c * r_void)) * self.omega_v * self.e_v * curvature_factor

    # ------------------------------------------------------------------
    # Loop growth
    # ------------------------------------------------------------------
    def loop_growth_rate(self, loop_type, CiL, CiL_i, phi_i, phi_v):
        """
        Growth rate dC_{iL,i}/dt for loops of a given type.
        Compact form from rate-equation derivation.
        """
        d = self.input_data.derived
        l     = d['l']
        l_hkl = d[f'l_{loop_type}']
        return (l_hkl / l) * np.sqrt(CiL_i * CiL) * (
            self.Z_iL_i * phi_i - self.Z_iL_v * phi_v
        )

    # ------------------------------------------------------------------
    # Loop nucleation / transformation
    # ------------------------------------------------------------------
    def loop_nucleation_rate(self, fi_vec, N_loop):
        """
        Rate at which f_i(N_loop) clusters nucleate new loops.
        Simplified: nucleation = omega_i * fi[N_loop]^2 (di-nucleation at cutoff).
        """
        if len(fi_vec) < N_loop:
            return 0.0
        fn = fi_vec[N_loop - 1]  # 0-indexed
        return self.omega_i * fn**2

    def loop_transition_rate(self, CiL_111):
        """1/2<111> → <100> transformation rate."""
        nu_tr = 1e11
        E_tr  = 2.4  # eV
        k_tr  = nu_tr * np.exp(-E_tr / (self.k_B * self.T))
        return k_tr * CiL_111

    def frac_100(self, dpa):
        """Temperature-dependent fraction of nucleation into <100> loops (RBF fit)."""
        from scipy.interpolate import Rbf
        points = [
            (300, 15, 77), (330, 15, 72), (330, 15, 45), (330, 32, 27),
            (330, 32, 31), (250, 16.2, 0), (350, 16.2, 50), (450, 16.2, 100),
            (250, 13.4, 10), (300, 14.6, 27), (350, 17.4, 79), (400, 17.2, 87),
            (415, 18.1, 80),
        ]
        T_pts = np.array([p[0] for p in points])
        D_pts = np.array([p[1] for p in points])
        Z_pts = np.array([p[2] for p in points])
        rbf   = Rbf(T_pts, D_pts, Z_pts, function='multiquadric', smooth=1)
        dpa   = np.clip(dpa, 1e-10, None)
        return np.clip(rbf(self.T - 273.15, dpa) / 100, 0.0, 1.0)

    # ------------------------------------------------------------------
    # Void nucleation
    # ------------------------------------------------------------------
    def void_nucleation_rate(self, C_void, r_void):
        """Cascade-induced void nucleation with capture-volume saturation."""
        G     = self.input_data.material_params['G']
        Omega = self.Omega
        n_cap = self.input_data.derived['n_cap']
        V_void = (4/3) * np.pi * r_void**3
        return G * (Omega / V_void) * (1 - n_cap * C_void)

    def void_dissolution_rate(self, C_void, r_void):
        """Void dissolution (diffusive escape of vacancies)."""
        return C_void / (r_void**2 / self.D_v)
