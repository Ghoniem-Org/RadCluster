# -*- coding: utf-8 -*-
"""
input_data.py – ClusterDynamics material and irradiation parameters.

Reads material, physical, and model parameters from a 3-sheet Excel workbook
(Full_CD/input/input_parameters.xlsx) and computes all derived quantities
needed by the ODE system and the C++ solver.

Sheet layout
------------
Material_Environment  – T, P, rho_d
Physical_Properties   – lattice, migration/formation energies, attempt
                        frequencies, bias factors, surface energy
Model_Parameters      – cluster sizes, ODE tolerances, time span, Python
                        segments, and all C++ / window-solver parameters
                        (usable by Python LSODA, C++ CVODE/ARKODE, Phase I–IV
                        OpenMP window modes)

All units follow SI:
  lengths in m, concentrations dimensionless (atom fraction),
  energies in eV, time in s.
"""

import warnings
import numpy as np
import pandas as pd
from pathlib import Path

_kB = 8.617333262e-5          # Boltzmann constant [eV K^-1]

BASE_DIR   = Path(__file__).parent.parent
INPUT_FILE = BASE_DIR / 'input' / 'input_parameters.xlsx'


class InputData:
    """
    Material, irradiation, and model parameters for cluster dynamics.

    Reads from a 3-sheet Excel workbook.  Call ``display_parameters()``
    to inspect all loaded values.

    Parameters
    ----------
    excel_file : path-like, optional
        Path to ``input_parameters.xlsx``.  Defaults to
        ``Full_CD/input/input_parameters.xlsx``.
    Nv : int, optional
        Override ``Model_Parameters.Nv`` (max vacancy cluster size).
    Ni : int, optional
        Override ``Model_Parameters.Ni`` (max interstitial cluster size).
    """

    def __init__(self, excel_file=INPUT_FILE, Nv=None, Ni=None):
        self.excel_file = Path(excel_file)
        if not self.excel_file.is_file():
            raise FileNotFoundError(
                f"Excel file not found: {self.excel_file}\n"
                f"Working directory: {Path.cwd()}\n"
                f"Run Full_CD/create_excel.py to regenerate it."
            )
        print(f"Loading parameters from: {self.excel_file.resolve()}")
        self._load_data()
        # Apply caller overrides before deriving anything
        if Nv is not None:
            self.model_params['Nv'] = int(Nv)
        if Ni is not None:
            self.model_params['Ni'] = int(Ni)
        self.calculate_derived_parameters()
        self._validate()

    # ── I/O ──────────────────────────────────────────────────────────────────

    @staticmethod
    def _sheet_to_dict(df):
        """Convert a Notation/Value sheet to a plain Python dict."""
        if 'Notation' in df.columns and 'Value' in df.columns:
            return {str(k): v for k, v in zip(df['Notation'], df['Value'])
                    if pd.notna(k) and str(k).strip()}
        # Fallback: first two columns
        return {str(k): v for k, v in zip(df.iloc[:, 0], df.iloc[:, 1])
                if pd.notna(k) and str(k).strip()}

    def _load_data(self):
        """Read the three Excel worksheets into attribute dicts."""
        try:
            mat_df   = pd.read_excel(self.excel_file, sheet_name='Material_Environment')
            phys_df  = pd.read_excel(self.excel_file, sheet_name='Physical_Properties')
            model_df = pd.read_excel(self.excel_file, sheet_name='Model_Parameters')
        except Exception as exc:
            raise RuntimeError(f"Failed to read Excel file: {exc}") from exc

        # Per-sheet dicts (for structured access and new code)
        self.material_env   = self._sheet_to_dict(mat_df)
        self.physical_props = self._sheet_to_dict(phys_df)
        self.model_params   = self._sheet_to_dict(model_df)

        # Merged dict: reaction_rates.py / rate_equations.py use material_params
        # for all physics scalars (legacy interface from the hardcoded version).
        self.material_params = {**self.material_env, **self.physical_props}

        # Cast cluster sizes and integer model parameters to int
        _int_keys = ('Nv', 'Ni', 'n_segments', 'n_points', 'log_time',
                     'backend', 'lmm', 'linsol', 'mu', 'ml', 'max_order',
                     'ark_table', 'window_mode', 'window_check_every',
                     'window_w0_v', 'window_w0_i', 'window_expand_pad',
                     'window_min_active_i', 'window_prec', 'window_width',
                     'window_N_thresh', 'window_omp_threads',
                     'Ni_max', 'Ni_extend_margin')
        for k in _int_keys:
            if k in self.model_params:
                self.model_params[k] = int(self.model_params[k])

        print("Successfully loaded all three parameter sheets.")

    # ── Derived quantities ────────────────────────────────────────────────────

    def calculate_derived_parameters(self):
        """Compute all derived physics quantities from the loaded parameters."""
        p   = self.physical_props
        m   = self.material_params
        T   = float(m['T'])
        kBT = _kB * T

        a_m = float(p['a_m'])                   # m (read directly)

        self.derived = {
            'kB':   _kB,
            'a':    a_m,
            # Diffusion coefficients [m^2 s^-1]
            'Di':   float(p['nu_i']) * a_m**2 * np.exp(-float(p['E_m_i']) / kBT),
            'Dv':   float(p['nu_v']) * a_m**2 * np.exp(-float(p['E_m_v']) / kBT),
            'D2v':  float(p['nu_v']) * a_m**2 * np.exp(-float(p['E_m_2v']) / kBT),
            # Recombination coefficient [s^-1]  (48 * nu_i * exp(-E_m_i/kT))
            'alpha': 48 * float(p['nu_i']) * np.exp(-float(p['E_m_i']) / kBT),
            # Equilibrium concentrations (atom fraction)
            'Cv_eq':  np.exp(-float(p['E_f_v']) / kBT),
            'C2v_eq': 6 * np.exp(-(2 * float(p['E_f_v']) - float(p['E_b_2v'])) / kBT),
            # Interstitial-pair nucleation rate constant [s^-1]
            'K_nuc_i': float(p['C_i1']) * float(p['nu_i']) * np.exp(-float(p['E_m_i']) / kBT),
        }

        print(f"Derived:  T={T} K  Cv_eq={self.derived['Cv_eq']:.3e}"
              f"  alpha={self.derived['alpha']:.3e}"
              f"  Di={self.derived['Di']:.3e} m^2/s"
              f"  Dv={self.derived['Dv']:.3e} m^2/s")

    # ── Validation ────────────────────────────────────────────────────────────

    def _validate(self):
        T     = float(self.material_params['T'])
        P     = float(self.material_params['P'])
        rho_d = float(self.material_params['rho_d'])
        if not (300 <= T <= 1200):
            warnings.warn(f"Temperature {T} K is outside typical range [300–1200 K]")
        if not (1e-9 <= P <= 1e-3):
            warnings.warn(f"Dose rate {P} dpa/s is outside typical range [1e-9–1e-3]")
        if not (1e12 <= rho_d <= 1e15):
            warnings.warn(f"Dislocation density {rho_d} m^-2 outside typical range")
        if self.derived['Cv_eq'] > 1e-6:
            warnings.warn(f"Cv_eq={self.derived['Cv_eq']:.2e} seems high – check T and E_f_v")

    # ── Convenience accessors ─────────────────────────────────────────────────

    @property
    def Nv(self):
        return int(self.model_params['Nv'])

    @property
    def Ni(self):
        return int(self.model_params['Ni'])

    # ── Display ───────────────────────────────────────────────────────────────

    def display_parameters(self):
        """Print all loaded and derived parameters."""
        sections = [
            ('MATERIAL & ENVIRONMENT', self.material_env),
            ('PHYSICAL PROPERTIES',    self.physical_props),
            ('MODEL PARAMETERS',       self.model_params),
            ('DERIVED PARAMETERS',     self.derived),
        ]
        for title, d in sections:
            print(f"\n{'='*60}")
            print(title)
            print('='*60)
            for k, v in d.items():
                fmt = f"  {k}: {v:.4e}" if isinstance(v, float) else f"  {k}: {v}"
                print(fmt)


# ── Self-test ─────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    inp = InputData()
    inp.display_parameters()
    print("\nSelf-test passed.")
