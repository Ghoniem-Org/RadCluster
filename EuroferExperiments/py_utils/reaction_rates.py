# -*- coding: utf-8 -*-
import numpy as np
from scipy.interpolate import Rbf

class ReactionRates:
    """
    Class to calculate all reaction rates for the rate equations system
    Contains all the frequency calculations and reaction rate member functions
    """
    def __init__(self, input_data, rate_equations=None):
        """
        Initialize with input data and optional rate_equations reference
        """
        self.input_data = input_data
        self.rate_equations = rate_equations  # Store reference
        self.k_B = 8.617e-5  # Boltzmann constant [eV/K]
        self.T = input_data.material_params['T']
        self.current_concentrations = None
        self.current_time = 0.0

        # Calculate basic frequencies
        self.calculate_basic_frequencies()

        # Calculate thermal emission probabilities
        self.calculate_thermal_emissions()

        # Calculate diffusion coefficients
        self.calculate_diffusion_coefficients()

    def calculate_basic_frequencies(self):
        """
        Calculate basic frequencies (Section 2.1)
        """
        z_c = self.input_data.physical_props['z_c']
        nu_i = self.input_data.physical_props['nu_i']
        nu_v = self.input_data.physical_props['nu_v']
        E_m_i = self.input_data.physical_props['E_m_i']
        E_m_v = self.input_data.physical_props['E_m_v']
        E_m_2i = self.input_data.physical_props['E_m_2i']
        G = self.input_data.material_params['G']

        # Basic frequencies (Equations 1-4)
        self.omega_i = z_c * nu_i * np.exp(-E_m_i / (self.k_B * self.T))
        self.omega_2i = z_c * nu_i * np.exp(-E_m_2i / (self.k_B * self.T))
        self.omega_3i = self.omega_2i  # Assuming same as di-interstitials
        self.omega_v = z_c * nu_v * np.exp(-E_m_v / (self.k_B * self.T))
        self.omega_irr = G  # Simplified - could include b_r factor

        print("omega_v, omega_i", self.omega_v, self.omega_i)

    def calculate_thermal_emissions(self):
        """
        Calculate thermal emission probabilities (Section 2.2)
        """
        E_F_v = self.input_data.physical_props['E_F_v']
        E_B_2i = self.input_data.physical_props['E_b_2i']
        E_B_3i = self.input_data.physical_props['E_b_3i']

        # Basic thermal emission (Equation 12)
        self.e_v = np.exp(-E_F_v / (self.k_B * self.T))
        self.e_2i = np.exp(-E_B_2i / (self.k_B * self.T))
        self.e_3i = np.exp(-E_B_3i / (self.k_B * self.T))

    def update_state(self, concentrations, time):
        """
        Update current state for calculations
        """
        self.current_concentrations = concentrations
        self.current_time = time

    def calculate_diffusion_coefficients(self):
        """
        Calculate diffusion coefficients (Equations 23-24)
        """
        a = self.input_data.physical_props['a']
        z_c = self.input_data.physical_props['z_c']

        self.D_v = (a ** 2 / z_c) * self.omega_v
        self.D_i = (a ** 2 / z_c) * self.omega_i
        print("D_v, D_i",self.D_v, self.D_i,)
        return self.D_v, self.D_i

    def calculate_sink_concentrations(self, concentrations):
        """
        Calculate equivalent sink concentrations (Equations 27-31)
        """
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        # Extract concentrations
        CiL_111, Cil_100, Cvoid = concentrations[4:7]

        # Dislocation densities
        rho_N = self.input_data.material_params['rho']

        # Loop radii
        riL_111 = self.rate_equations.calculate_ril_111(concentrations)
        ril_100 = self.rate_equations.calculate_ril_100(concentrations)

        atomic_volume = self.input_data.physical_props['Omega']
        rho_iL_111 = 2 * np.pi * riL_111 * CiL_111 / atomic_volume
        # print(riL_111)
        rho_iL_100 = 2 * np.pi * ril_100 * Cil_100 / atomic_volume

        # Bias factors
        Z_i = self.input_data.model_params.get('Z_N')
        Z_iL = self.input_data.model_params.get('Z_iL')

        a = self.input_data.physical_props['a']
        z_c = self.input_data.physical_props['z_c']

        # Vacancy / Interstitial sink
        self.C_v_s = (a ** 2 / z_c) * (rho_N + (rho_iL_111 + rho_iL_100))
        self.C_i_s = (a ** 2 / z_c) * (1.05 * rho_N + 1.05 * (rho_iL_111 + rho_iL_100))

        return self.C_v_s, self.C_i_s

    def R_v_s(self, concentrations=None):
        """Reaction rate between vacancies and sinks"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        Cv = concentrations[0]

        C_v_s, _ = self.calculate_sink_concentrations(concentrations)
        print("C_v_s", C_v_s)
        return self.omega_v * Cv * (C_v_s+6e-4)
        # return self.omega_v * Cv * C_v_s

    def R_i_s(self, concentrations=None):
        """Reaction rate between interstitials and sinks"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        Ci = concentrations[1]

        _, C_i_s = self.calculate_sink_concentrations(concentrations)
        print("C_i_s", C_i_s)
        return self.omega_i * Ci * (C_i_s+2e-4)
        # return self.omega_i * Ci * C_i_s

    def R_2i_s(self, concentrations=None):
        """Reaction rate between di-interstitials and sinks"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        C2i = concentrations[2]

        _, C_i_s = self.calculate_sink_concentrations(concentrations)
        return self.omega_i * C2i * C_i_s

    def R_3i_s(self, concentrations=None):
        """Reaction rate between tri-interstitials and sinks"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        C3i = concentrations[3]

        _, C_i_s = self.calculate_sink_concentrations(concentrations)

        return self.omega_i * C3i * C_i_s

    def R_i_v(self, concentrations=None):
        """Reaction rate between interstitials and vacancies"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations=np.maximum(concentrations, 0e-20)
        recom = self.input_data.physical_props['recom']
        Ci, Cv = concentrations[1], concentrations[0]
        omega_iv= recom *(self.omega_i + self.omega_v)
        return omega_iv * Ci * Cv

    def R_2i_v(self, concentrations=None):
        """Reaction rate between di-interstitials and vacancies"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        C2i, Cv = concentrations[2], concentrations[0]
        omega_2iv = self.omega_2i + self.omega_v

        return omega_2iv * C2i * Cv

    def R_3i_v(self, concentrations=None):
        """Reaction rate between tri-interstitials and vacancies"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        C3i, Cv = concentrations[3], concentrations[0]
        omega_3iv = self.omega_3i + self.omega_v

        return omega_3iv * C3i * Cv

    def R_i_i(self, concentrations=None):
        """Reaction rate between interstitials"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        Ci = concentrations[1]

        return 2 * self.omega_i * Ci ** 2

    def R_i_2i(self, concentrations=None):
        """Reaction rate between interstitials and di-interstitials"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        Ci, C2i = concentrations[1], concentrations[2]
        omega_i2i = self.omega_i + self.omega_2i

        return omega_i2i * Ci * C2i

    def R_i_3i(self, concentrations=None):
        """Reaction rate between interstitials and tri-interstitials"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        Ci, C3i = concentrations[1], concentrations[3]
        omega_i3i = self.omega_i + self.omega_3i

        return omega_i3i * Ci * C3i

    def R_2i_2i(self, concentrations=None):
        """Reaction rate between di-interstitials"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        C2i = concentrations[2]
        omega_2i2i = self.omega_2i + self.omega_2i

        return 2 * omega_2i2i * C2i ** 2

    def R_2i_3i(self, concentrations=None):
        """Reaction rate between di-interstitials"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        C2i, C3i = concentrations[2], concentrations[3]
        omega_2i3i = self.omega_2i + self.omega_3i

        return omega_2i3i * C2i * C3i

    def G_v(self):
        """Vacancy generation rate"""
        return self.input_data.derived['G_v']

    def G_i(self):
        """Interstitial generation rate"""
        return self.input_data.derived['G_i']

    def G_2i(self):
        """Di-interstitial generation rate"""
        return self.input_data.derived['G_2i'] / 2  # Divided by 2 for pairs

    def G_3i(self):
        """Tri-interstitial generation rate"""
        return self.input_data.derived['G_3i'] / 3  # Divided by 3 for pairs


    def nucleation_rate_iL(self, concentrations=None):
        """Nucleation rate for interstitial loops"""
        if concentrations is None:
            concentrations = self.current_concentrations

        concentrations = np.maximum(concentrations, 0e-20)

        # Nucleation from i+3i and 2i+2i reactions
        rate = self.R_i_3i(concentrations) + self.R_2i_2i(concentrations)
        return rate

    def T_il(self, concentrations=None):
        """Transition rate for interstitial loops from 1/2<111> to <100>."""
        if concentrations is None:
            concentrations = self.current_concentrations

        concentrations = np.maximum(concentrations, 0.0)
        CiL_111 = concentrations[4]

        nu_tr = 1e11  # attempt frequency [1/s]
        E_tr = 2.4  # activation barrier [eV] (tunable)
        tran = nu_tr * np.exp(-E_tr / (self.k_B * self.T))
        return tran * CiL_111

    def get_frac_100(self):
        dpa = self.input_data.material_params['G'] * self.current_time
        dpa = np.clip(dpa, 1e-10, None)  # avoid log(0)

        points = [(300, 15, 77), (330, 15, 72), (330, 15, 45), (330, 32, 27), (330, 32, 31), (250, 16.2, 0),
                  (350, 16.2, 50), (450, 16.2, 100), (250, 13.4, 10), (300, 14.6, 27), (350, 17.4, 79), (400, 17.2, 87),
                  (415, 18.1, 80)]

        T = np.array([p[0] for p in points])
        D = np.array([p[1] for p in points])
        Z = np.array([p[2] for p in points])
        rbf = Rbf(T, D, Z, function='multiquadric', smooth=1)

        frac_100 = rbf(self.T-273.15, dpa)/100
        return frac_100

    def f100(self, concentrations=None):
        """Fraction of nucleation rate for <100> interstitial loops."""
        if concentrations is None:
            concentrations = self.current_concentrations

        concentrations = np.maximum(concentrations, 0.0)
        frac_100 = self.get_frac_100()
        return frac_100 * self.nucleation_rate_iL(concentrations)

    def f111(self, concentrations=None):
        """Fraction of nucleation rate for 1/2<111> interstitial loops."""
        if concentrations is None:
            concentrations = self.current_concentrations

        concentrations = np.maximum(concentrations, 0.0)
        frac_100 = self.get_frac_100()

        return (1 - frac_100) * self.nucleation_rate_iL(concentrations)

    def emission_2i(self, concentrations=None):
        """Thermal emission rate from di-interstitials"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        C2i = concentrations[2]
        return self.omega_2i * self.e_2i * C2i

    def emission_3i(self, concentrations=None):
        """Thermal emission rate from tri-interstitials"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        C3i = concentrations[3]
        return self.omega_3i * self.e_3i * C3i

    def flux_v(self, concentrations=None):
        """Vacancy flux to sinks"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        Cv = concentrations[0]
        return self.omega_v * Cv

    def flux_i(self, concentrations=None):
        """Interstitial flux to sinks"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        Ci, C2i, C3i = concentrations[1], concentrations[2], concentrations[3]

        return self.omega_i * (Ci + 2 * C2i + 3 * C3i)

    def loop_growth_rate_iL(self, loop_type, concentrations=None):
        """Compute the growth rate for interstitial loops (⟨111⟩ or ⟨100⟩)."""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        # --- Extract concentrations and parameters ---
        CiL_111, CiL_100 = concentrations[4], concentrations[5]
        CiL_i_111, CiL_i_100 = concentrations[8], concentrations[9]

        l = self.input_data.derived['l']
        l_111 = self.input_data.derived['l_111']
        l_100 = self.input_data.derived['l_100']
        Z_iL_i = self.input_data.model_params.get('Z_iL_i', 1.5)
        Z_iL_v = self.input_data.model_params.get('Z_iL_v', 0.7)
        # --- Compute rate for each loop type ---
        if loop_type == "111":
            rate = (l_111 / l) * np.sqrt(CiL_i_111 * CiL_111) * (Z_iL_i * self.flux_i(concentrations) - Z_iL_v * self.flux_v(concentrations))
        elif loop_type == "100":
            rate = (l_100 / l) * np.sqrt(CiL_i_100 * CiL_100) * (Z_iL_i * self.flux_i(concentrations) - Z_iL_v * self.flux_v(concentrations))
        else:
            raise ValueError("loop_type must be '111' or '100'.")
        # print(f"flux_i={self.flux_i(concentrations):.3e}, flux_v={self.flux_v(concentrations):.3e}, rate={rate:.3e}")
        return rate

    def nucleation_rate_void(self, concentrations=None):
        """Compute the void nucleation term"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        Cvoid = concentrations[6]
        r_void = concentrations[7]
        n_cap = self.input_data.derived['n_cap']

        # From cascade generation with capture volume limitation
        Omega = self.input_data.physical_props["Omega"]
        V_void = 4 / 3 * np.pi * (r_void ** 3)
        # print((Omega/V_void))
        rate = self.input_data.material_params['G'] * (Omega/V_void) * (1 - n_cap * Cvoid)
        # rate = self.input_data.derived['G_void'] * (1 - n_cap * Cvoid)

        return rate

    def tau_diss(self, concentrations=None):
        """Compute the void dissociation term"""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        N_void = concentrations[6]
        r_void = concentrations[7]
        # vacancies diffuse away from the void.
        rate = N_void/((r_void**2)/self.D_v)
        return rate

    def Void_v(self, concentrations=None):
        """Vacancy absorption rate contribution to void growth."""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        Cv = concentrations[0]
        r_void = concentrations[7]

        a = self.input_data.physical_props['a']
        z_a = self.input_data.physical_props['z_c']

        return a**2/(z_a*r_void) * self.omega_v * Cv

    def Void_i(self, concentrations=None):
        """Interstitial absorption rate contribution to void shrinkage."""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        Ci, C2i, C3i = concentrations[1], concentrations[2], concentrations[3]
        r_void = concentrations[7]

        a = self.input_data.physical_props['a']
        z_a = self.input_data.physical_props['z_c']

        return a**2/(z_a*r_void) * self.omega_i * Ci

    def Void_emission(self, concentrations=None):
        """Vacancy emission rate from void surface."""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        gamma = self.input_data.physical_props["gamma"]
        Omega = self.input_data.physical_props["Omega"]
        r_void = concentrations[7]

        a = self.input_data.physical_props['a']
        z_a = self.input_data.physical_props['z_c']
        return a**2/(z_a*r_void) * self.omega_v * self.e_v * (np.exp((2.0*gamma*Omega)/(r_void*self.k_B*self.T)) - 1.0)

    def Cal_Cveq(self, concentrations=None):
        """Vacancy emission rate from void surface."""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        gamma = self.input_data.physical_props["gamma"]
        Omega = self.input_data.physical_props["Omega"]
        r_void = concentrations[7]

        Cveq = self.e_v * (np.exp((2.0*gamma*Omega)/(r_void*self.k_B*self.T)) - 1.0)
        return Cveq

    def Trapping(self, type, concentrations=None):
        """Trapping rate from external impurities."""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        Ci, Cv = concentrations[1], concentrations[0]
        Ctrap_i, Ctrap_v = concentrations[10], concentrations[11]
        CT0_v, CT0_i = self.input_data.derived["CT0_v"], self.input_data.derived["CT0_i"]
        Omega = self.input_data.physical_props["Omega"]
        r_trap = self.input_data.derived["r_trap"]

        Ctrap_i = np.clip(Ctrap_i, 0.0, CT0_i)
        Ctrap_v = np.clip(Ctrap_v, 0.0, CT0_v)

        if type == "interstitial":
            k_trap = 4*np.pi*r_trap*self.D_i/Omega
            # print("k_trap",k_trap)
            return k_trap * Ci * (CT0_i - Ctrap_i)
        if type == "vacancy":
            k_trap = 4*np.pi*r_trap*self.D_v/Omega
            return k_trap * Cv * (CT0_i - Ctrap_v)

    def Release(self, type, concentrations=None):
        """Release rate from external impurities."""
        if concentrations is None:
            concentrations = self.current_concentrations
        concentrations = np.maximum(concentrations, 0e-20)

        Ctrap_i, Ctrap_v = concentrations[10], concentrations[11]
        CT0_v, CT0_i = self.input_data.derived["CT0_v"], self.input_data.derived["CT0_i"]

        Ctrap_i = np.clip(Ctrap_i, 0.0, CT0_i)
        Ctrap_v = np.clip(Ctrap_v, 0.0, CT0_v)
        nu_i, nu_v = self.input_data.physical_props['nu_i'], self.input_data.physical_props['nu_v']
        if type == "interstitial":
            k_rel = nu_i * np.exp(-1.8 / (self.k_B*self.T))
            return k_rel * Ctrap_i
        if type == "vacancy":
            k_rel = nu_v * np.exp(-1 / (self.k_B*self.T))
            return k_rel * Ctrap_v

