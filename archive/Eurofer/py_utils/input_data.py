# -*- coding: utf-8 -*-
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent
INPUT_FILE = BASE_DIR / 'input' / 'input_parameters_Neutron.xlsx'

class InputData:
    """Class to handle input data from Excel workbook with three worksheets"""
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

    def _df_to_dict(self, df):
        """Convert DataFrame to dictionary using 'Notation' column as keys"""
        if 'Notation' in df.columns and 'Value' in df.columns:
            return dict(zip(df['Notation'], df['Value']))
        else:
            return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))

    def load_data(self):
        """Load data from the three Excel worksheets"""
        try:
            # Read the three worksheets
            self.material_env_df = pd.read_excel(self.excel_file, sheet_name='Material_Environment')
            self.physical_props_df = pd.read_excel(self.excel_file, sheet_name='Physical_Properties')
            self.model_params_df = pd.read_excel(self.excel_file, sheet_name='Model_Parameters')

            # Convert to dictionaries for easier access
            self.material_params = self._df_to_dict(self.material_env_df)
            self.physical_props = self._df_to_dict(self.physical_props_df)
            self.model_params = self._df_to_dict(self.model_params_df)

            print("Successfully loaded input data from Excel file")

        except Exception as e:
            raise RuntimeError(f"Error loading Excel file: {e}")

    def validate_setup(self):
        """Validate the simulation setup and warn about potential issues"""
        print("Validating simulation setup...")

        # Check for reasonable parameter values
        T = self.material_params['T']
        G = self.material_params['G']
        rho = self.material_params['rho']

        # Temperature validation
        if not (200 <= T <= 1000):
            warnings.warn(f"Temperature {T} K may be outside typical range [200-1000 K]")

        # Dose rate validation
        if not (1e-8 <= G <= 1e-3):
            warnings.warn(f"Dose rate {G} dpa/s may be outside typical range [1e-8 - 1e-3 dpa/s]")

        # Dislocation density validation
        if not (1e12 <= rho <= 1e16):
            warnings.warn(f"Dislocation density {rho} $m^-2$ may be outside typical range [1e12 - 1e16 $m^-2$]")

        # Check derived parameters
        C_v_eq = self.derived['C_v_eq']
        if C_v_eq > 1e-6:
            warnings.warn(f"Equilibrium vacancy concentration {C_v_eq:.2e} seems high")

        print("Setup validation complete.")

    def calculate_derived_parameters(self):
        """Calculate derived parameters from input data"""
        k_B = 8.617e-5  # Boltzmann constant [eV/K]
        self.derived = {}

        # Capture atoms around a void
        self.derived['n_cap'] = 4 * np.pi * (self.model_params['r_cap_void'] ** 3) / (3 * self.physical_props['Omega'])

        # Thermal equilibrium vacancy concentration
        self.derived['C_v_eq'] = np.exp(-self.physical_props['E_F_v'] / (k_B * self.material_params['T']))

        # Compact notation parameters (Section 3.1)
        self.derived['l'] = (self.physical_props['z_c'] * self.physical_props['Omega'] / (2 * np.pi * self.physical_props['a'] ** 2))
        self.derived['l_111'] = np.sqrt(self.physical_props['Omega'] / (np.pi * self.physical_props['b_111']))
        self.derived['l_100'] = np.sqrt(self.physical_props['Omega'] / (np.pi * self.physical_props['b_100']))

        # Generation rates
        self.derived['G_v'] = self.material_params['G'] * (1 - self.physical_props['epsilon_void'])
        self.derived['G_i'] = self.material_params['G'] * (1 - 2 * self.physical_props['epsilon_2i'] - 3 * self.physical_props['epsilon_3i'])
        self.derived['G_2i'] = self.material_params['G'] * self.physical_props['epsilon_2i']
        self.derived['G_3i'] = self.material_params['G'] * self.physical_props['epsilon_3i']
        self.derived['G_void'] = self.material_params['G'] * self.physical_props['epsilon_void']

        # Trapping
        self.derived["r_trap"] = 1.5*0.249e-9
        self.derived["CT0"] = 200 * 1e-6
        self.derived["CT0_i"], self.derived["CT0_v"] = 0.5*self.derived["CT0"], (1-0.5)*self.derived["CT0"]


        # Network dislocation parameters
        if 'rho' in self.material_params:
            # Calculate network sink strength
            self.derived['k_N_v'] = 4 * np.pi * self.material_params['rho']  # Vacancy sink strength
            self.derived['k_N_i'] = 4 * np.pi * self.material_params['rho'] * self.model_params.get('Z_N_i', 1.1)  # Interstitial sink strength

        print("Derived parameters calculated successfully")
        print(f"  C_v_eq = {self.derived['C_v_eq']:.2e}")
        print(f"  G_v = {self.derived['G_v']:.2e}")
        print(f"  G_i = {self.derived['G_i']:.2e}")


    def display_parameters(self):
        """
        Display all parameters in a formatted way
        """
        print("\n" + "="*60)
        print("MATERIAL & ENVIRONMENT PARAMETERS")
        print("="*60)
        for key, value in self.material_params.items():
            print(f"{key}: {value}")

        print("\n" + "="*60)
        print("PHYSICAL PROPERTIES")
        print("="*60)
        for key, value in self.physical_props.items():
            print(f"{key}: {value}")

        print("\n" + "="*60)
        print("MODEL PARAMETERS")
        print("="*60)
        for key, value in self.model_params.items():
            print(f"{key}: {value}")

        print("\n" + "="*60)
        print("DERIVED PARAMETERS")
        print("="*60)
        for key, value in self.derived.items():
            print(f"{key}: {value}")


# Test function
def test_input_data():
    """Test the InputData class"""
    try:
        # Test with default file
        print("Testing InputData class...")
        input_data = InputData()
        input_data.display_parameters()
        print("Test completed successfully")

    except Exception as e:
        raise RuntimeError(f"Test failed: {e}")


if __name__ == "__main__":
    test_input_data()