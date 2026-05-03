# -*- coding: utf-8 -*-
import numpy as np

class RateEquations:
    def __init__(self, input_data, reaction_rates):
        """
        Initialize with input data and reaction rates
        """
        self.input_data = input_data
        self.reaction_rates = reaction_rates
        self.reaction_rates.rate_equations = self

        self.concentration_names = [
            'Cv',  # 0 - Vacancy concentration
            'Ci',  # 1 - Interstitial concentration
            'C2i',  # 2 - Di-interstitial concentration
            'C3i',  # 3 - Tri-interstitial concentration
            'CiL_111',  # 4 - 1/2<111> interstitial loop concentration
            'CiL_100',  # 5 - <100> interstitial loop concentration
            'C_void',  # 6 - Void number density
            'r_void',  # 7 - Void radius
            'CiL_i_111',  # 8 - Interstitial atoms in 1/2<111> loops
            'CiL_i_100',  # 9 - Interstitial atoms in <100> loops
            'Ctrap_i', # 10 - Trap concentration for interstitial
            'Ctrap_v', # 11 - Trap concentration for Vacancy
        ]
        self.n_equations = len(self.concentration_names)
        print(f"Initialized rate equations system with {self.n_equations} equations")

    def ode_system(self, t, y):
        """
        Main ODE system function

        Parameters:
        t: float - current time
        y: array - current concentrations

        Returns: dydt: array - time derivatives
        """

        self.reaction_rates.update_state(y, t)
        dydt = np.zeros(self.n_equations)

        dydt[0] = self.dCv_dt(y)  # Vacancy rate equation
        dydt[1] = self.dCi_dt(y)  # Interstitial rate equation
        dydt[2] = self.dC2i_dt(y)  # Di-n_capinterstitials
        dydt[3] = self.dC3i_dt(y)  # Tri-interstitials
        dydt[4] = self.dCiL_111_dt(y)  # 1/2<111> interstitial loops
        dydt[5] = self.dCiL_100_dt(y)  # <100> interstitial loops
        dydt[6] = self.dC_void_dt(y)  # Void nucleation rate
        dydt[7] = self.dr_void_dt(y)  # Void growth rate
        dydt[8] = self.dCiL_i_111_dt(y)  # Interstitial atoms in 1/2<111> loops
        dydt[9] = self.dCiL_i_100_dt(y)  # Interstitial atoms in <100> loops
        dydt[10] = self.dC_trap_i_dt(y)
        dydt[11] = self.dC_trap_v_dt(y)
        return dydt

    def set_concentration_by_name(self, y, name, value):
        """
        Set concentration by name in the y array
        """
        if name in self.concentration_names:
            index = self.concentration_names.index(name)
            y[index] = value
        else:
            raise ValueError(f"Unknown concentration name: {name}")

    def get_initial_conditions(self, custom_values=None):
        """
        Generate physically reasonable initial conditions
        """
        y0 = np.zeros(self.n_equations)

        # Default initial conditions
        C_v_eq = self.input_data.derived['C_v_eq']
        scale_1 = 1e-6
        scale_2 = 1e-10

        # --- Default initial conditions dictionary ---
        defaults = {
            'Cv': C_v_eq,  # thermal equilibrium vacancy conc.
            'Ci': scale_1 * C_v_eq,                 # small interstitial supersaturation
            'C2i': scale_2 * C_v_eq,                # di-interstitials
            'C3i': scale_2 * C_v_eq,                # tri-interstitials
            'CiL_111': scale_2 * C_v_eq,            # <111> interstitial loops
            'CiL_100': scale_2 * C_v_eq,            # <100> interstitial loops
            'C_void': scale_1 * C_v_eq,             # void number density
            'r_void': 2e-9,                         # void radius
            'CiL_i_111': scale_2 * C_v_eq,          # interstitials in <111> loops
            'CiL_i_100': scale_2 * C_v_eq,          # interstitials in <100> loops
            "Ctrap_i": 0.01*self.input_data.derived['CT0_i'],            # Trap from impurities
            "Ctrap_v": 0.01*self.input_data.derived['CT0_v'],
        }

        # Apply custom values if provided
        if custom_values:
            defaults.update(custom_values)

        # Set initial conditions
        for name, value in defaults.items():
            self.set_concentration_by_name(y0, name, value)

        print("Initial conditions set:")
        for i, name in enumerate(self.concentration_names):
            print(f"  {name}: {y0[i]:.2e}")

        # Additional validation
        if np.any(y0 < 0):
            raise ValueError("Negative initial concentrations detected")

        if np.any(np.isnan(y0)) or np.any(np.isinf(y0)):
            raise ValueError("Invalid initial concentrations (NaN or Inf)")
        return y0

    def dCv_dt(self, y):
        """Vacancy concentration rate equation"""
        y = np.maximum(y, 0e-20)

        # Generation term
        generation = self.reaction_rates.G_v()

        # Consumption terms
        reactions = self.reaction_rates.R_i_v(y) + self.reaction_rates.R_2i_v(y) + self.reaction_rates.R_3i_v(y) + self.reaction_rates.R_v_s(y)

        trap = self.reaction_rates.Release(type="vacancy", concentrations=y) - self.reaction_rates.Trapping(type="vacancy", concentrations=y)
        return generation - reactions
        # print(self.reaction_rates.Release(type="vacancy", concentrations=y),self.reaction_rates.Trapping(type="vacancy", concentrations=y))
        # return generation - reactions + trap

    def dCi_dt(self, y):
        """Interstitial concentration rate equation"""
        y = np.maximum(y, 0e-20)

        # Generation term
        generation = self.reaction_rates.G_i() + self.reaction_rates.R_2i_v(y)

        # Emission terms
        emission = 2. * self.reaction_rates.emission_2i(y) + 3. * self.reaction_rates.emission_3i(y)

        # Consumption terms
        consumption = self.reaction_rates.R_i_v(y) + self.reaction_rates.R_i_i(y) + self.reaction_rates.R_i_2i(y) + self.reaction_rates.R_i_3i(y) + self.reaction_rates.R_i_s(y)

        trap = self.reaction_rates.Release(type="interstitial", concentrations=y) - self.reaction_rates.Trapping(type="interstitial", concentrations=y)
        # return generation + emission - consumption + trap
        return generation + emission - consumption

    def dC_trap_i_dt(self, y):
        "Trap concentration rate equation"
        y = np.maximum(y, 0e-20)

        # Trapping terms
        trap = self.reaction_rates.Trapping(type="interstitial", concentrations=y)
        # Release terms
        release = self.reaction_rates.Release(type="interstitial", concentrations=y)
        return trap - release

    def dC_trap_v_dt(self, y):
        "Trap concentration rate equation"
        y = np.maximum(y, 0e-20)

        # Trapping terms
        trap = self.reaction_rates.Trapping(type="vacancy", concentrations=y)
        # Release terms
        release = self.reaction_rates.Release(type="vacancy", concentrations=y)
        return trap - release

    def dC2i_dt(self, y):
        """Di-interstitial concentration rate equation"""
        y = np.maximum(y, 0e-20)  # Element-wise maximum

        # Generation from cascades
        generation = self.reaction_rates.G_2i()

        # Formation from i+i
        formation = self.reaction_rates.R_i_i(y) + self.reaction_rates.R_3i_v(y) \
            # +self.reaction_rates.emission_2i(y)

        # Consumption terms
        consumption = self.reaction_rates.R_2i_v(y) + self.reaction_rates.R_i_2i(y) + self.reaction_rates.R_2i_2i(y) + self.reaction_rates.R_2i_s(y)
        return generation + formation - consumption

    def dC3i_dt(self, y):
        """Tri-interstitial concentration rate equation"""
        y = np.maximum(y, 0e-20)  # Element-wise maximum

        # Generation from cascades
        generation = self.reaction_rates.G_3i()

        # Formation from i+2i
        formation = self.reaction_rates.R_i_2i(y)

        # Consumption terms
        consumption = self.reaction_rates.R_3i_v(y) + self.reaction_rates.R_i_3i(y) + self.reaction_rates.R_3i_s(y) + self.reaction_rates.R_2i_3i(y)
        return generation + formation - consumption

    def dCiL_111_dt(self, y):
        """1/2<111> interstitial loop concentration rate equation"""
        y = np.maximum(y, 0e-20)

        # Formation
        formation = self.reaction_rates.f111(y)

        # Transformation
        transformation = self.reaction_rates.T_il(y)
        return formation - transformation

    def dCiL_100_dt(self, y):
        """<100> interstitial loop concentration rate equation"""
        y = np.maximum(y, 0e-20)

        # Formation
        formation = self.reaction_rates.f100(y)

        # Transformation
        transformation = self.reaction_rates.T_il(y)
        return formation + transformation

    def dC_void_dt(self, y):
        """Void number density rate equation"""
        y = np.maximum(y, 0e-20)

        # Formation
        formation = self.reaction_rates.nucleation_rate_void(y)

        # Dissolution
        dissolution = self.reaction_rates.tau_diss(y)
        # print(f"{formation:.2e}, {dissolution:.2e}, {formation - dissolution:.2e}")
        return formation - dissolution


    def dr_void_dt(self, y):
        """Void radius growth rate equation"""
        y = np.maximum(y, 0e-20)

        # Growth
        growth = self.reaction_rates.Void_v(y)

        # Annihilating
        annihilating = self.reaction_rates.Void_i(y)

        # Emission
        emission = self.reaction_rates.Void_emission(y)
        return growth - annihilating - emission

    def dCiL_i_111_dt(self, y):
        """Interstitials in interstitial loops rate equation"""
        y = np.maximum(y, 0e-20)
        return self.reaction_rates.loop_growth_rate_iL(loop_type='111', concentrations=y)

    def dCiL_i_100_dt(self, y):
        """Interstitials in interstitial loops rate equation"""
        y = np.maximum(y, 0e-20)
        return self.reaction_rates.loop_growth_rate_iL(loop_type='100', concentrations=y)

    def calculate_ril_111(self, y):
        """Calculate 1/2<111> interstitial loop radius"""
        y = np.maximum(y, 1e-20)
        CiL_111, CiL_i_111 = y[4], y[8]
        l_111 = self.input_data.derived['l_111']
        if CiL_111 > 0:
            return l_111 * np.sqrt(CiL_i_111 / CiL_111)
        else:
            return 5e-9

    def calculate_ril_100(self, y):
        """Calculate <100> interstitial loop radius"""
        y = np.maximum(y, 1e-20)
        CiL_100, CiL_i_100 = y[5], y[9]
        l_100 = self.input_data.derived['l_100']
        if CiL_100 > 0:
            return l_100 * np.sqrt(CiL_i_100 / CiL_100)
        else:
            return 5e-9
