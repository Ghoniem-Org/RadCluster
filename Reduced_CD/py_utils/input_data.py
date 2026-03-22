# -*- coding: utf-8 -*-
"""
input_data.py — ClusterDynamics input loader

Reads the 3-sheet Excel workbook and computes all derived parameters
needed by the cluster-dynamics ODE system.

Sheets expected:
  - Material_Environment  : T (K), G (dpa/s), rho (m^-2)
  - Physical_Properties   : Omega, a, b_111, b_100, migration/formation energies, ...
  - Model_Parameters      : N_v, N_i, N_loop, Z-factors, capture radii, ...
"""

import warnings
import numpy as np
import pandas as pd
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent
INPUT_FILE = BASE_DIR / 'input' / 'input_parameters.xlsx'


class InputData:
    """
    Load and validate simulation inputs from the standard 3-sheet Excel workbook.
    Extends EuroferMicrostructure.InputData with cluster-size cutoffs N_v, N_i, N_loop.
    """

    def __init__(self, excel_file=INPUT_FILE):
        self.excel_file = str(excel_file)
        file_path = Path(self.excel_file)
        if not file_path.is_file():
            raise FileNotFoundError(
                f"Excel file not found: {file_path}\n"
                f"Current working directory: {Path.cwd()}"
            )
        print(f"Found Excel file: {file_path.resolve()}")
        self.load_data()
        self.calculate_derived_parameters()
        self.validate_setup()

    # ------------------------------------------------------------------
    def _df_to_dict(self, df):
        if 'Notation' in df.columns and 'Value' in df.columns:
            return dict(zip(df['Notation'], df['Value']))
        return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))

    def load_data(self):
        try:
            self.material_env_df   = pd.read_excel(self.excel_file, sheet_name='Material_Environment')
            self.physical_props_df = pd.read_excel(self.excel_file, sheet_name='Physical_Properties')
            self.model_params_df   = pd.read_excel(self.excel_file, sheet_name='Model_Parameters')

            self.material_params = self._df_to_dict(self.material_env_df)
            self.physical_props  = self._df_to_dict(self.physical_props_df)
            self.model_params    = self._df_to_dict(self.model_params_df)
            print("Successfully loaded input data from Excel file")
        except Exception as e:
            raise RuntimeError(f"Error loading Excel file: {e}")

    # ------------------------------------------------------------------
    def calculate_derived_parameters(self):
        k_B = 8.617e-5  # eV/K
        T   = self.material_params['T']
        self.derived = {}

        # Equilibrium vacancy concentration
        self.derived['C_v_eq'] = np.exp(-self.physical_props['E_F_v'] / (k_B * T))

        # Compact loop-growth length scales (same as EuroferMicrostructure)
        Omega = self.physical_props['Omega']
        a     = self.physical_props['a']
        z_c   = self.physical_props['z_c']
        b_111 = self.physical_props['b_111']
        b_100 = self.physical_props['b_100']
        self.derived['l']     = z_c * Omega / (2 * np.pi * a**2)
        self.derived['l_111'] = np.sqrt(Omega / (np.pi * b_111))
        self.derived['l_100'] = np.sqrt(Omega / (np.pi * b_100))

        # Generation rates
        eps_2i   = self.physical_props['epsilon_2i']
        eps_3i   = self.physical_props['epsilon_3i']
        eps_void = self.physical_props['epsilon_void']
        G = self.material_params['G']
        self.derived['G_v']    = G * (1 - eps_void)
        self.derived['G_i']    = G * (1 - 2*eps_2i - 3*eps_3i)
        self.derived['G_2i']   = G * eps_2i
        self.derived['G_3i']   = G * eps_3i
        self.derived['G_void'] = G * eps_void

        # Void capture volume
        r_cap = self.model_params.get('r_cap_void', 0.5e-9)
        self.derived['n_cap'] = 4 * np.pi * r_cap**3 / (3 * Omega)

        # Trap parameters (inherited from EuroferMicrostructure defaults)
        self.derived['r_trap'] = 1.5 * 0.249e-9
        CT0 = 200e-6
        self.derived['CT0']   = CT0
        self.derived['CT0_i'] = 0.5  * CT0
        self.derived['CT0_v'] = 0.5  * CT0

        # --- Cluster-dynamics specific ---
        # Cluster size cutoffs (default to EuroferMicrostructure-compatible trimers if absent)
        self.derived['N_v']    = int(self.model_params.get('N_v',    5))
        self.derived['N_i']    = int(self.model_params.get('N_i',    5))
        self.derived['N_loop'] = int(self.model_params.get('N_loop', 4))

        # Monomer radius (Wigner-Seitz)
        self.derived['r0'] = (3 * Omega / (4 * np.pi))**(1/3)

        # Dislocation network sink strength
        rho = self.material_params.get('rho', 1e14)
        Z_N = self.model_params.get('Z_N', 1.05)
        self.derived['k2_N_v'] = 2 * np.pi * rho
        self.derived['k2_N_i'] = 2 * np.pi * rho * Z_N

        print("Derived parameters calculated.")
        print(f"  C_v_eq = {self.derived['C_v_eq']:.2e}")
        print(f"  N_v = {self.derived['N_v']},  N_i = {self.derived['N_i']},  N_loop = {self.derived['N_loop']}")

    # ------------------------------------------------------------------
    def validate_setup(self):
        T   = self.material_params['T']
        G   = self.material_params['G']
        rho = self.material_params.get('rho', 1e14)
        if not (200 <= T <= 1000):
            warnings.warn(f"Temperature {T} K outside typical range [200–1000 K]")
        if not (1e-8 <= G <= 1e-2):
            warnings.warn(f"Dose rate {G} dpa/s outside typical range")
        if not (1e12 <= rho <= 1e16):
            warnings.warn(f"Dislocation density {rho} m^-2 outside typical range")
        N_loop = self.derived['N_loop']
        N_i    = self.derived['N_i']
        if N_loop >= N_i:
            warnings.warn(f"N_loop ({N_loop}) >= N_i ({N_i}): no clusters will feed loops")
        print("Validation complete.")

    # ------------------------------------------------------------------
    def display_parameters(self):
        for section, data in [
            ("MATERIAL & ENVIRONMENT", self.material_params),
            ("PHYSICAL PROPERTIES",    self.physical_props),
            ("MODEL PARAMETERS",       self.model_params),
            ("DERIVED",                self.derived),
        ]:
            print(f"\n{'='*55}\n{section}\n{'='*55}")
            for k, v in data.items():
                print(f"  {k}: {v}")
