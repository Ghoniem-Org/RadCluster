"""
input_data.py — InputData class for Eurofer_CD cluster dynamics.

Reads material, irradiation, and model parameters from a 3-sheet Excel
workbook (Eurofer_CD/input/input_parameters.xlsx) and computes all derived
quantities needed by the ODE system.

Sheet layout
------------
Material_Environment  — T, G, rho_d, neutron_spectrum
Physical_Properties   — lattice, migration/formation energies, He parameters,
                        surface energy, bias factors, cascade spectrum
Model_Parameters      — Ni, Nv, L_He_max, he_mode, ODE tolerances, time span

Physics reference
-----------------
Ghoniem, N.M. (2024), "Formulation of Cluster Dynamics Equations for
Irradiated Ferritic-Martensitic Steels" (Sections 1–5).

All units: SI (lengths in m, concentrations dimensionless as atom fractions,
energies in eV, time in s).
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
    Material, irradiation, and model parameters for Eurofer_CD.

    Reads from a 3-sheet Excel workbook.  Call ``display_parameters()``
    to inspect all loaded values.

    Parameters
    ----------
    excel_file : path-like, optional
        Defaults to ``Eurofer_CD/input/input_parameters.xlsx``.
    Nv : int, optional
        Override Model_Parameters.Nv (max vacancy cluster size).
    Ni : int, optional
        Override Model_Parameters.Ni (max interstitial cluster size).
    he_mode : str, optional
        Override Model_Parameters.he_mode.
        Choices: 'decoupled' | 'fast_eq' | 'full'
    """

    def __init__(self, excel_file=INPUT_FILE, Nv=None, Ni=None, he_mode=None):
        self.excel_file = Path(excel_file)
        if not self.excel_file.is_file():
            raise FileNotFoundError(
                f"Excel file not found: {self.excel_file}\n"
                f"Run Eurofer_CD/create_excel.py to generate it."
            )
        print(f"Loading parameters from: {self.excel_file.resolve()}")
        self._load_data()

        # Apply caller overrides before deriving anything
        if Nv is not None:
            self.model_params['Nv'] = int(Nv)
        if Ni is not None:
            self.model_params['Ni'] = int(Ni)
        if he_mode is not None:
            self.model_params['he_mode'] = str(he_mode)

        self.calculate_derived_parameters()
        self._validate()

    # ── I/O ──────────────────────────────────────────────────────────────────

    @staticmethod
    def _sheet_to_dict(df):
        """Convert a Notation/Value sheet to a plain Python dict."""
        if 'Notation' in df.columns and 'Value' in df.columns:
            return {str(k): v for k, v in zip(df['Notation'], df['Value'])
                    if pd.notna(k) and str(k).strip()}
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

        self.material_env   = self._sheet_to_dict(mat_df)
        self.physical_props = self._sheet_to_dict(phys_df)
        self.model_params   = self._sheet_to_dict(model_df)

        # Merged dict for physics modules (legacy interface)
        self.material_params = {**self.material_env, **self.physical_props}

        # Cast integer parameters
        _int_keys = ('Nv', 'Ni', 'L_He_max', 'n_segments', 'n_points',
                     'log_time', 'n1D', 'm1_spec', 'n1_spec')
        for k in _int_keys:
            if k in self.model_params:
                self.model_params[k] = int(float(self.model_params[k]))
            if k in self.material_params:
                self.material_params[k] = int(float(self.material_params[k]))

        print("Successfully loaded all three parameter sheets.")

    # ── Derived quantities ────────────────────────────────────────────────────

    def calculate_derived_parameters(self):
        """Compute all derived physics quantities from the loaded parameters."""
        p    = self.physical_props
        m    = self.material_params
        T    = float(m['T'])
        kBT  = _kB * T
        a_m  = float(p['a_m'])
        Omega = float(p['Omega'])

        # Capture radius r_0 = (3Ω/4π)^(1/3)  [m]
        r0 = (3.0 * Omega / (4.0 * np.pi)) ** (1.0 / 3.0)

        # Diffusion coefficients [m² s⁻¹]
        Di  = float(p['nu_i']) * a_m**2 * np.exp(-float(p['E_m_i']) / kBT)
        Dv  = float(p['nu_v']) * a_m**2 * np.exp(-float(p['E_m_v']) / kBT)
        DHe = float(p['nu_He']) * a_m**2 * np.exp(-float(p['E_m_He']) / kBT)

        # Equilibrium concentrations (atom fractions)
        Cv_eq  = np.exp(-float(p['E_f_v']) / kBT)
        C2v_eq = 6.0 * np.exp(-(2.0*float(p['E_f_v']) - float(p['E_b_2v'])) / kBT)

        # Recombination coefficient α = 48·ν_i·exp(−E_m_i/kT)  [s⁻¹]
        alpha = 48.0 * float(p['nu_i']) * np.exp(-float(p['E_m_i']) / kBT)

        # Nucleation rate constant for di-interstitial formation
        K_nuc_i = float(p.get('C_i1', 1.0)) * Di

        # He-vacancy reduction mode
        he_mode = str(self.model_params.get('he_mode', 'decoupled')).strip().lower()

        # Surface energy capillary parameter A_v = 4π·γ_s·r_0²  [eV]
        _J_to_eV = 6.241509074e18
        A_void = 4.0 * np.pi * float(p['gamma_s']) * r0**2 * _J_to_eV

        self.derived = {
            'kB':       _kB,
            'a':        a_m,
            'Omega':    Omega,
            'r0':       r0,
            'Di':       Di,
            'Dv':       Dv,
            'DHe':      DHe,
            'alpha':    alpha,
            'Cv_eq':    Cv_eq,
            'C2v_eq':   C2v_eq,
            'K_nuc_i':  K_nuc_i,
            'A_void':   A_void,
            'he_mode':  he_mode,
        }

        print(f"Derived:  T={T} K  Cv_eq={Cv_eq:.3e}"
              f"  Di={Di:.3e} m2/s  Dv={Dv:.3e} m2/s"
              f"  he_mode='{he_mode}'")

    # ── Validation ────────────────────────────────────────────────────────────

    def _validate(self):
        T     = float(self.material_params['T'])
        G     = float(self.material_params.get('G', 1e-6))
        rho_d = float(self.material_params.get('rho_d', 5e14))

        if not (300 <= T <= 1200):
            warnings.warn(f"Temperature {T} K is outside typical range [300–1200 K]")
        if not (1e-9 <= G <= 1e-3):
            warnings.warn(f"Dose rate {G} dpa/s is outside typical range [1e-9–1e-3]")
        if not (1e12 <= rho_d <= 1e16):
            warnings.warn(f"Dislocation density {rho_d} m^-2 outside typical range")
        if self.derived['Cv_eq'] > 1e-6:
            warnings.warn(f"Cv_eq={self.derived['Cv_eq']:.2e} seems high – check T and E_f_v")

        he_mode = self.derived['he_mode']
        if he_mode not in ('decoupled', 'fast_eq', 'full'):
            warnings.warn(f"Unknown he_mode='{he_mode}'; using 'decoupled'")
            self.derived['he_mode'] = 'decoupled'

    # ── Convenience accessors ─────────────────────────────────────────────────

    @property
    def Nv(self):
        return int(self.model_params['Nv'])

    @property
    def Ni(self):
        return int(self.model_params['Ni'])

    @property
    def L_He_max(self):
        return int(self.model_params.get('L_He_max', 5))

    @property
    def he_mode(self):
        return str(self.model_params.get('he_mode', 'decoupled')).strip().lower()

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
