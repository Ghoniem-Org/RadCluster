#!/usr/bin/env python3
"""
create_excel.py — Generate input_parameters.xlsx for Eurofer_CD.

Creates the 3-sheet workbook with bcc Fe / EUROFER97 default parameters
based on the formulation in:
  Ghoniem, N.M. (2024), "Formulation of Cluster Dynamics Equations for
  Irradiated Ferritic-Martensitic Steels" (Sections 1–5).

Run from the Eurofer_CD/ directory:
    python create_excel.py

Overwrites input/input_parameters.xlsx.

Sheet layout
------------
Material_Environment  — irradiation conditions: T, G, rho_d, neutron_spectrum
Physical_Properties   — bcc Fe / Eurofer97 material parameters
Model_Parameters      — cluster sizes, ODE settings, He mode
"""

from pathlib import Path
import pandas as pd

OUTPUT_FILE = Path(__file__).parent / 'input' / 'input_parameters.xlsx'


# ---------------------------------------------------------------------------
# Parameter tables
# ---------------------------------------------------------------------------

MATERIAL_ENV = [
    # (Notation, Value, Units, Description)
    ('T',               600.0,     'K',        'Irradiation temperature'),
    ('G',               1.0e-6,    'dpa/s',    'Displacement damage rate'),
    ('rho_d',           5.0e14,    'm^-2',     'Network dislocation density'),
    ('neutron_spectrum','fission',  '-',        'Neutron spectrum: fission or fusion'),
]

PHYSICAL_PROPS = [
    # bcc Fe / Eurofer97 lattice and thermodynamic parameters
    # Section 2, Table 5 of Rate_Equations.pdf; Malerba et al. (2021)
    ('a_m',             2.87e-10,   'm',        'BCC Fe lattice parameter'),
    ('Omega',           1.178e-29,  'm^3',      'Atomic volume (a^3/2 for bcc)'),

    # Vacancy properties
    ('E_f_v',           1.73,       'eV',       'Vacancy formation energy'),
    ('E_m_v',           0.55,       'eV',       'Vacancy migration energy'),
    ('nu_v',            6.25e12,    'Hz',       'Vacancy attempt frequency'),
    ('E_b_2v',          0.30,       'eV',       'Di-vacancy binding energy'),

    # Interstitial (SIA) properties
    ('E_m_i',           0.013,      'eV',       'SIA migration energy (1D glide, bcc Fe)'),
    ('nu_i',            6.25e12,    'Hz',       'SIA attempt frequency'),
    ('E_b_2i',          0.80,       'eV',       'Di-interstitial binding energy'),
    ('E_b_inf_loop',    1.80,       'eV',       'Large-loop SIA binding limit (½⟨111⟩)'),
    ('n_trans_loop',    8.0,        '-',        'Loop binding saturation cluster size'),
    ('n_1D',            4,          '-',        'Min. SIA cluster size for 1D glide'),
    ('D_1D_factor',     1.0,        '-',        '1D-to-3D effective diffusion ratio'),

    # Helium properties
    ('E_f_He',          4.37,       'eV',       'He interstitial formation energy in Fe'),
    ('E_m_He',          0.06,       'eV',       'He migration energy in bcc Fe'),
    ('nu_He',           6.25e12,    'Hz',       'He attempt frequency'),
    ('E_b_HeV',         2.60,       'eV',       'He binding energy in mono-vacancy (Caturla)'),
    ('E_b_He_delta',   -0.80,       'eV',       'He pressure fit parameter delta_He'),
    ('E_b_He_beta',     0.70,       '-',        'He pressure fit exponent beta_He'),

    # Surface / interface energy for capillary model
    ('gamma_s',         1.50,       'J/m^2',    'Fe surface energy for capillary model'),

    # Dislocation sink bias factors
    ('Z_i',             1.05,       '-',        'Interstitial bias factor at dislocations'),
    ('Z_v',             1.00,       '-',        'Vacancy capture factor at dislocations'),

    # Cascade production parameters (fission spectrum, Stoller 2000)
    ('eta',             0.30,       '-',        'Cascade survival efficiency (N_d/N_NRT)'),
    ('f_i_cl',          0.58,       '-',        'Fraction of SIA born in clusters'),
    ('f_v_cl',          0.15,       '-',        'Fraction of vacancies born in clusters'),
    ('s_i',             1.6,        '-',        'SIA cluster power-law exponent'),
    ('s_v',             2.5,        '-',        'Vacancy cluster power-law exponent'),
    ('m1_spec',         20,         '-',        'Max SIA cluster size in production spectrum'),
    ('n1_spec',         10,         '-',        'Max vacancy cluster size in production spectrum'),

    # He transmutation production
    ('G_He_per_dpa',    0.5,        'appm/dpa', 'He production rate (transmutation)'),
]

MODEL_PARAMS = [
    # Cluster size truncation
    ('Ni',              100,        '-',        'Max SIA cluster size tracked'),
    ('Nv',              100,        '-',        'Max vacancy cluster size tracked'),
    ('L_He_max',        5,          '-',        'Max He per void (full/fast_eq modes)'),

    # He-vacancy state-space reduction mode (PDF Section 5.6.5):
    #   decoupled : He modifies void binding energy only (for fission; low He/dpa)
    #   fast_eq   : Fast He equilibration within each void class (general)
    #   full      : Explicit (m, ell) pairs  [large state space]
    ('he_mode',         'decoupled', '-',       'He-vacancy reduction mode'),

    # ODE solver settings
    ('t_begin',         1.0e-8,     's',        'Integration start time'),
    ('t_end',           1.0e6,      's',        'Integration end time'),
    ('n_points',        200,        '-',        'Number of output time points'),
    ('n_segments',      60,         '-',        'Segments for segmented LSODA'),
    ('rtol',            1.0e-8,     '-',        'Relative ODE tolerance'),
    ('atol',            1.0e-50,    '-',        'Absolute ODE tolerance'),
    ('method',          'LSODA',    '-',        'ODE solver method'),
    ('log_time',        1,          '-',        'Use log-spaced output times (1=yes)'),

    # Floor concentration (prevents log(0) in Jacobian)
    ('C_floor',         1.0e-100,   '-',        'Minimum allowed concentration'),
]


# ---------------------------------------------------------------------------
# Build and write
# ---------------------------------------------------------------------------

def _make_df(rows, has_units=True):
    if has_units:
        cols = ['Notation', 'Value', 'Units', 'Description']
    else:
        cols = ['Notation', 'Value', 'Description']
    return pd.DataFrame(rows, columns=cols)


def create_excel(output_file=OUTPUT_FILE):
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    df_mat   = _make_df(MATERIAL_ENV)
    df_phys  = _make_df(PHYSICAL_PROPS)
    df_model = _make_df(MODEL_PARAMS)

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df_mat.to_excel(writer,   sheet_name='Material_Environment', index=False)
        df_phys.to_excel(writer,  sheet_name='Physical_Properties',  index=False)
        df_model.to_excel(writer, sheet_name='Model_Parameters',     index=False)

    print(f"Written: {output_file.resolve()}")
    print(f"  Material_Environment : {len(df_mat)} rows")
    print(f"  Physical_Properties  : {len(df_phys)} rows")
    print(f"  Model_Parameters     : {len(df_model)} rows")


if __name__ == '__main__':
    create_excel()
