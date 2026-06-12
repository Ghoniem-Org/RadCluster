#!/usr/bin/env python3
"""
create_excel.py — Generate input_parameters.xlsx for RadCluster_2_0.

Creates a 5-sheet workbook with exact bcc Fe / EUROFER97 parameters from:
  Ghoniem, N.M. (2026), "A Cluster Dynamics Model for Radiation Damage
  Evolution in Ferritic-Martensitic Steels" (Rate_Equations.pdf).

Sheet layout
------------
Production    — cascade production parameters (Tables 2, 5): fission vs. fusion
Energetics    — lattice, migration/formation energies, He EOS (Tables 5, 8)
Diffusion     — diffusion pre-factors, solute trapping (Tables 5, 16)
Dissociation  — binding energy parameters, trap mutation (Tables 18-19, 27)
Reactions     — rate constant prefactors, sink parameters, solver settings (Tables 25-26, 29)

Run from RadCluster_2_0/:
    python create_excel.py

Overwrites input/input_parameters.xlsx.
"""

from pathlib import Path
import pandas as pd

OUTPUT_FILE = Path(__file__).parent / 'input' / 'input_parameters.xlsx'


# ── Sheet 1: Production ───────────────────────────────────────────────────────
# Tables 2 and 5 of Rate_Equations.pdf

PRODUCTION = [
    # (Parameter, Symbol, Fission, Fusion, Units)
    ('Survival and Dose Rate', None, None, None, None),
    ('Survival efficiency (spectrum-avg.)', 'eta',      0.30,   0.28,   '−'),
    ('Displacement rate',                   'phi_dot',  'user', 'user', 'dpa/s'),
    ('He production rate',                  'G_He/G',   '0.5−1','~10',  'appm He/dpa'),

    ('Interstitial Cluster Production', None, None, None, None),
    ('SIA clustering fraction',    'f_cl_i', 0.58,   0.65,   '−'),     # Table 2
    ('SIA power-law exponent',     's_i',    1.6,    1.5,    '−'),     # Table 2, Eq. 7
    ('Max SIA cascade size',       'i_cascade', 20,   50,     'SIAs'),  # Table 2
    ('SIA normalization constant', 'C_i',    0.1093, 0.0553, '−'),     # Eq. 9

    ('Vacancy Cluster Production', None, None, None, None),
    ('Vacancy clustering fraction',    'f_cl_v', 0.15,   0.20,   '−'),     # Table 2
    ('Vacancy power-law exponent',     's_v',    2.5,    2.3,    '−'),     # Table 2, Eq. 8
    ('Max vacancy cascade size',       'v_cascade', 10,  20,     'vacancies'),  # Table 2
    ('Vacancy normalization constant', 'C_v',    0.1506, 0.1296, '−'),     # Eq. 10
]

# ── Sheet 2: Energetics ───────────────────────────────────────────────────────
# Tables 5 and 8 of Rate_Equations.pdf

ENERGETICS = [
    # (Parameter, Symbol, Value, Units, Notes)
    ('Lattice and Elastic Constants', None, None, None, None),
    ('Lattice constant',           'a',      0.2867,     'nm',    'bcc Fe'),       # Table 5
    ('Atomic volume',              'Omega',  1.18e-29,   'm^3',   'a^3/2'),        # Table 5
    ('Burgers vector (1/2<111>)',  'b_111',  0.2482,     'nm',    None),           # Table 5
    ('Burgers vector (<100>)',     'b_100',  0.2867,     'nm',    None),           # Table 5
    ('Shear modulus (RT)',         'mu',     82,         'GPa',   None),           # Table 5
    ('Poisson ratio',              'nu',     0.29,       '−',     None),           # Table 5
    ('Surface energy',             'gamma_s',2.0,        'J/m^2', 'DFT average'),  # Table 5

    ('Vacancy Energetics', None, None, None, None),
    ('Formation energy',           'E_f_v',  2.0,        'eV',   'DFT consensus'),  # Table 5
    ('Migration energy',           'E_m_v',  0.67,       'eV',   'DFT, Stage III'), # Table 5
    ('Divacancy binding (2NN)',    'E_B_2v', 0.22,       'eV',   'DFT'),            # Table 5
    ('Divacancy migration',        'E_m_2v', 0.70,       'eV',   'DFT'),            # Table 5

    ('SIA Energetics', None, None, None, None),
    ('Formation energy (<110> dumbbell)', 'E_f_i',  3.64, 'eV', 'DFT'),           # Table 5
    ('Migration energy (<110> T-R)',      'E_m_i',  0.34, 'eV', 'DFT, Stage I_E'),# Table 5
    ('Di-SIA binding energy',             'E_B_2i', 0.80, 'eV', 'DFT'),           # Table 5
    ('Di-SIA migration energy',           'E_m_2i', 0.42, 'eV', 'DFT'),           # Table 5

    ('Helium Energetics', None, None, None, None),
    ('He migration energy (tetrahedral)', 'E_m_h',    0.06, 'eV', 'DFT'),                 # Table 5
    ('He-vacancy binding (1st He)',       'E_b_hV_1', 2.30, 'eV', 'Fu & Willaime 2005'),  # Table 5
    ('He-vacancy binding (2nd He)',       'E_b_hV_2', 2.00, 'eV', 'Fu & Willaime 2005'),  # Table 5
    ('He solution energy (interstitial)', 'E_s_He',   2.35, 'eV', 'Fu & Willaime 2005'),  # Table 5

    ('Attempt Frequencies', None, None, None, None),
    ('Vacancy',        'nu_v', 1.0e13, 's^-1', None),   # Table 5
    ('SIA (<110>)',    'nu_i', 1.0e13, 's^-1', None),   # Table 5
    ('He interstitial','nu_h', 3.0e12, 's^-1', None),   # Table 5

    ('He Equation of State (Virial, 600-1000 K)', None, None, None, None),
    ('Second virial coefficient', 'B2', 1.67e-29, 'm^3/atom', 'Eq. 64, Table 8'),  # Table 8
    ('Third virial coefficient',  'B3', 1.84e-58, 'm^6/atom^2','Eq. 65, Table 8'), # Table 8
]

# ── Sheet 3: Diffusion ────────────────────────────────────────────────────────
# Tables 5 and 16 of Rate_Equations.pdf

DIFFUSION = [
    # (Parameter, Symbol, Value, Units, Notes)
    ('Pure Fe Baseline Diffusion', None, None, None, None),
    ('Vacancy D0',  'D0_v',   8.2e-8,  'm^2/s', 'Eq. 17'),   # Table 5
    ('Vacancy E_m', 'E_m_v',  0.67,    'eV',    None),
    ('SIA D0',      'D0_i',   8.2e-8,  'm^2/s', 'Eq. 17'),   # Table 5
    ('SIA E_m',     'E_m_i',  0.34,    'eV',    None),
    ('He D0',       'D0_h',   2.5e-8,  'm^2/s', 'Eq. 17'),   # Table 5
    ('He E_m',      'E_m_h',  0.06,    'eV',    None),

    ('SIA Cluster 1D Glide (n >= 4)', None, None, None, None),
    ('Attempt frequency',  'nu0_1D',  6.0e12, 's^-1',  'Eq. 33'),           # Table 5
    ('1D migration energy','E_m_1D',  0.03,   'eV',    'mid-range estimate'),# Table 5
    ('Size exponent',      's_1D',    0.7,    '−',     'MD range 0.5-1.0'), # Table 5
    ('Max mobile SIA size', 'i_mobile', 100,   'SIAs',  None),               # Table 5

    ('Vacancy Cluster Diffusion', None, None, None, None),
    ('Size exponent',      's_vc',     1.0, '−', None),            # Table 5
    ('Max mobile vacancy size', 'v_mobile', 5,  'vacancies', None),    # Table 5

    ('EUROFER Dissolved Solute Concentrations', None, None, None, None),
    ('Cr', 'c_Cr', 0.094,  'at/at', '9 wt% Cr'),           # Table 16
    ('W',  'c_W',  0.0033, 'at/at', '1.1 wt% W'),          # Table 16
    ('Mn', 'c_Mn', 0.0047, 'at/at', '0.47 wt% Mn'),        # Table 16
    ('C (dissolved)', 'c_C', 5.0e-4, 'at/at', 'after tempering'),  # Table 16
    ('N (dissolved)', 'c_N', 2.0e-4, 'at/at', 'after tempering'),  # Table 16

    ('SIA Solute Trapping: z_s * exp(E_b / kT)  (Eq. 42)', None, None, None, None),
    ('C-SIA binding',  'E_b_C_SIA',  0.45, 'eV', 'z=4'),  # Table 16
    ('N-SIA binding',  'E_b_N_SIA',  0.40, 'eV', 'z=4'),  # Table 16
    ('Cr-SIA binding', 'E_b_Cr_SIA', 0.10, 'eV', 'z=8'),  # Table 16
    ('Mn-SIA binding', 'E_b_Mn_SIA', 0.20, 'eV', 'z=6'),  # Table 16

    ('Vacancy Solute Trapping: z_s * exp(E_b / kT)  (Eq. 48)', None, None, None, None),
    ('C-V binding',  'E_b_C_V',   0.45, 'eV', 'z=3'),  # Table 16
    ('N-V binding',  'E_b_N_V',   0.40, 'eV', 'z=3'),  # Table 16
    ('W-V binding',  'E_b_W_V',   0.27, 'eV', 'z=8'),  # Table 16
    ('Mn-V binding', 'E_b_Mn_V',  0.10, 'eV', 'z=8'),  # Table 16
    ('Cr-V binding', 'E_b_Cr_V',  0.05, 'eV', 'z=8'),  # Table 16

    ('SIA Cluster Loop Trapping (1D, Eq. 52)', None, None, None, None),
    ('C-loop binding',  'E_b_C_loop',  0.50, 'eV', 'z=2'),  # Table 16
    ('N-loop binding',  'E_b_N_loop',  0.40, 'eV', 'z=2'),  # Table 16
    ('Cr-loop binding', 'E_b_Cr_loop', 0.10, 'eV', 'z=4'),  # Table 16

    ('Mixed 1D/3D Parameters', None, None, None, None),
    ('Mean free path (EUROFER)', 'L_hat', 50.0,   '−',  'L/a, alloy-dependent'),  # Table 5
    ('1D/3D interpolation const.','B_rot', 2.627,  '−',  'Eq. 128'),              # Eq. 128
]

# ── Sheet 4: Dissociation ─────────────────────────────────────────────────────
# Tables 18-19, 27 of Rate_Equations.pdf

DISSOCIATION = [
    # (Parameter, Symbol, Value, Units, Notes)
    ('Void Binding Energy Model (Eqs. 66-67)', None, None, None, None),
    ('Vacancy formation energy', 'E_f_v',   2.0,    'eV',    None),
    ('Surface energy',           'gamma_s', 2.0,    'J/m^2', None),
    ('Divacancy binding (DFT)',  'E_b_v2',  0.22,   'eV',    'anchor point'),
    ('Decay constant',           'lambda',  0.5756, 'vac^-1','ln(100)/8'),

    ('Atomistic Fitting Amplitudes A(m) — Table 18', None, None, None, None),
    ('A(m=0)', 'A_void_0', 1.2353, 'eV', 'pure void correction'),  # Table 18
    ('A(m=1)', 'A_void_1', 2.9064, 'eV', None),                    # Table 18
    ('A(m=2)', 'A_void_2', 3.4147, 'eV', None),                    # Table 18
    ('A(m=3)', 'A_void_3', 2.1504, 'eV', None),                    # Table 18
    ('A(m=4)', 'A_void_4',-0.1590, 'eV', None),                    # Table 18

    ('Vacancy Loop Binding (Eqs. 104-105)', None, None, None, None),
    ('Stacking fault energy', 'gamma_sf', 0.6,    'J/m^2', None),
    ('Core radius',           'r_c',      0.333,  'b_111', None),
    ('Divacancy binding',     'E_b_v2',   0.22,   'eV',    None),
    ('Decay constant',        'lambda_v', 0.575,  'vac^-1',None),

    ('Interstitial Loop Binding (power-law + blend, Eqs. 106-108)', None, None, None, None),
    ('1/2<111> amplitude',  'A_111',   0.7501, 'eV',   None),   # Table 18 / Eq. 106
    ('1/2<111> exponent',   'B_111',   0.3873, '−',    None),   # Table 18 / Eq. 106
    ('<100> amplitude',     'A_100',   0.7160, 'eV',   None),   # Table 18 / Eq. 107
    ('<100> exponent',      'B_100',   0.3581, '−',    None),   # Table 18 / Eq. 107
    ('Blend center',        'n_tr',    25,     'SIAs', None),   # Eq. 108
    ('Blend width',         'sigma_tr',5,      'SIAs', None),   # Eq. 108

    ('He Maximum Loading', None, None, None, None),
    ('He/V ratio parameter', 'alpha_He', 1.7, '−', 'l_max = alpha * m^(2/3)'),  # Table 29

    ('He Binding Energy to Bubbles (Eqs. 76-77)', None, None, None, None),
    ('He atomic mass',        'm_He',  4.0026,  'amu',       '4-He'),
    ('He virial B2',          'B2',    1.67e-29,'m^3/atom',  'Table 8'),
    ('He virial B3',          'B3',    1.84e-58,'m^6/atom^2','Table 8'),

    ('He Binding: Atomistic-Continuum Blending (Table 19)', None, None, None, None),
    ('Decay constant',             'mu_He',    0.658, 'He^-1', 'ln(100)/7'),    # Table 19
    ('Fitting amplitude A^He(2)',  'A_He_2',   0.55,  'eV',    'Caturla/Fu'),   # Table 19
    ('Fitting amplitude A^He(3)',  'A_He_3',   0.40,  'eV',    'Caturla/Fu'),   # Table 19
    ('Fitting amplitude A^He(4)',  'A_He_4',   0.75,  'eV',    'Caturla/Fu'),   # Table 19

    ('Trap Mutation (Eq. 83, Eq. 142, Table 27)', None, None, None, None),
    ('Attempt frequency',   'nu0_TM',      1.0e12, 's^-1', None),           # Table 27
    ('E_TM: He5V1',         'E_TM_1_5',    1.00,   'eV',   'Morishita 2003'),# Table 27
    ('E_TM: He6V1',         'E_TM_1_6',    0.50,   'eV',   'Morishita 2003'),# Table 27
    ('E_TM: He7V1',         'E_TM_1_7',    0.00,   'eV',   'spontaneous'),   # Table 27
    ('E_TM: He4 (alpha=0)', 'E_TM_0_4',    0.30,   'eV',   'Gao 2011'),      # Table 27
    ('E_TM: He5 (alpha=0)', 'E_TM_0_5',    0.10,   'eV',   'Gao 2011'),      # Table 27
    ('E_TM: He6 (alpha=0)', 'E_TM_0_6',    0.00,   'eV',   'spontaneous'),   # Table 27
]

# ── Sheet 5: Reactions ────────────────────────────────────────────────────────
# Tables 25-26, 29 of Rate_Equations.pdf

REACTIONS = [
    # (Parameter, Symbol, Value, Units, Notes)
    ('Geometric Rate Constant Prefactors (Eq. 128, Table 25)', None, None, None, None),
    ('Spherical cluster',     'A_sph',  7.818, '−', '(48*pi^2)^(1/3); Eq. 128'),  # Table 25
    ('Dislocation loop',      'A_loop', 10.78, '−', '8*sqrt(pi/sqrt(3)); Eq. 128'),# Table 25
    ('Pure 1D sink',          'A_1D',   2.632, '−', '9/(8*pi^(2/3)); Eq. 128'),    # Table 25
    ('1D/3D interpolation',   'B_rot',  2.627, '−', '(4/pi)*(8*pi/3)^(1/3)'),      # Table 25
    ('V-SIA recombination',   'K_iv_pf',21.77, '−', '4*sqrt(3)*pi'),               # Eq. 130

    ('Dislocation Network Sink (Table 26)', None, None, None, None),
    ('Network dislocation density', 'rho_d',  1.0e14, 'm^-2',  'EUROFER tempered'),  # Table 26
    ('Interstitial bias factor',    'Z_i',    1.10,   '−',     'Table 26'),          # Table 26
    ('Vacancy bias factor',         'Z_v',    1.00,   '−',     'Table 26'),          # Table 26
    ('He bias factor',              'Z_He',   1.00,   '−',     'Table 26'),          # Table 26

    ('Grain Boundary Sink (Eq. 134-137)', None, None, None, None),
    ('Grain diameter', 'd_g', 5.0e-6, 'm', 'EUROFER typical'),   # Table 26

    ('Precipitate Sink (Eq. 134-137)', None, None, None, None),
    ('Precipitate density',  'rho_p', 1.0e21, 'm^-3', 'M23C6 precipitates'),  # Table 26
    ('Precipitate radius',   'r_p',   5.0e-9, 'm',    'EUROFER typical'),     # Table 26
    ('Precipitate bias SIA', 'Z_p_i', 1.0,    '−',    None),                  # Table 26
    ('Precipitate bias vac', 'Z_p_v', 1.0,    '−',    None),                  # Table 26

    ('Radiation Re-solution (Eq. 143)', None, None, None, None),
    ('Re-solution parameter (fission)', 'b0_fission', 0.01, 'dpa^-1/He', None),  # Table 29
    ('Re-solution parameter (fusion)',  'b0_fusion',  0.10, 'dpa^-1/He', None),  # Table 29

    ('He State-Space Parameters (Table 29)', None, None, None, None),
    ('He max loading parameter',  'alpha_He',  1.7,    '−',      None),          # Table 29
    ('Max SIA cluster tracked',    'I',         500,    '−',      'user: 500-5000'),# Table 29
    ('Max vacancy cluster tracked','V',         500,    '−',      'user: 500-5000'),# Table 29
    ('Max He per cluster',        'L_He_max',  'mf',   '−',      'mf=mean-field'), # Table 29

    ('Size-Bin Moment Parameters (Chapter 9)', None, None, None, None),
    ('Moments per bin',         'n_moments',   2,      '−',      'zeroth + first'),
    ('Max discrete SIA size',   'i_discrete',  1,      '−',      'sizes 1..i_discrete resolved individually'),
    ('Max discrete vacancy size','v_discrete', 5,      '−',      'sizes 1..v_discrete resolved individually; default = v_mobile'),
    ('Number of SIA bin-moment eqs',  'I_bin', 14,     '−',      'grouped bins for SIA clusters'),
    ('Number of VAC bin-moment eqs',  'V_bin', 14,     '−',      'grouped bins for vacancy clusters'),
    ('Grouping transition',     'n_group',     50,     '−',      'Eq. 190'),

    ('Loop 1/2<111> -> <100> Conversion (Marian 2002 + Dudarev 2008)', None, None, None, None),
    ('Unary barrier offset',      'E_a0_conv',    0.8,    'eV',  'Dudarev unary E_a^0; Marian dH2 ~ 1.0 eV'),
    ('Unary barrier size slope',  'gamma_a_conv', 0.03,   'eV',  'per perimeter-segment; calibrate to Arakawa onset + size window'),
    ('Unary attempt frequency',   'nu0_conv',     1.0e13, '1/s', 'Debye'),
    ('Conversion crossover temp', 'T_star_conv_C',450.0,  'C',   'dF(n_ref,T*)=0; in [350,550] (Dudarev Fig. 4)'),
    ('Conversion reference size', 'n_ref_conv',   50,     '-',   'calibration anchor size for LoopEnergetics'),
    ('Junction peak yield',       'phi_max_junc', 0.5,    '-',   'Marian; yield at n=n_prime (0-1)'),
    ('Junction log-size tol.',    'sigma_s_junc', 0.35,   '-',   'Marian comparable-size width in ln(n/n_prime)'),
    ('Junction min size',         'n_j_min_junc', 30,     '-',   'Marian; junctions from n ~ 34-37'),
    ('<100> loop-onset size',     'n_loop_min',   4,      '-',   'bulk-100 n_min; below this no <100> loop exists'),

    ('Solver Settings (Table 29)', None, None, None, None),
    ('Solver mode',      'solver_mode',   'full_system',       '−', 'full_system | active_window'),
    ('Physics option',   'physics_option','full_CD_fission',   '−', 'full_CD_fission | full_CD_fusion | bin_moment_CD_fission | bin_moment_CD_fusion'),
    ('Free He kinetics', 'he_kinetics',   'dynamic',           '−', 'dynamic | quasi_steady_state'),
    ('Start time',       't_begin',       1.0e-8,              's', None),
    ('End time',         't_end',         1.0e7,               's', None),
    ('Output points',    'n_points',      200,                 '−', None),
    ('Log time spacing', 'log_time',      1,                   '−', '1=yes, 0=no'),
    ('Rel. tolerance',   'rtol',          1.0e-8,              '−', None),
    ('Abs. tolerance',   'atol',          1.0e-50,             '−', None),
    ('Linear solver',    'linsol',        'dense',             '−', 'dense | band | gmres | klu (sparse, full_CD only)'),
    ('Band upper bw',    'mu_band',       0,                   '−', '0 = auto'),
    ('Band lower bw',    'ml_band',       0,                   '−', '0 = auto'),
    ('Window width',     'window_width',  100,                 '−', 'sliding window initial width (shared SIA/VAC)'),
    ('Window expand thr','window_C_exp',  1.0e-18,             '−', 'expand when C > this'),
    ('Concentration floor','C_floor',     1.0e-100,            '−', None),

    ('Irradiation Conditions', None, None, None, None),
    ('Temperature',           'T',      600.0,  'K',      'irradiation temperature'),
    ('Displacement rate',     'G',      1.0e-6, 'dpa/s',  None),
    ('He production rate',    'G_He_r', 1.0,    'appm/dpa','for fusion use ~10'),
    ('Neutron spectrum',      'spectrum','fission','−',   'fission | fusion'),
]


# ── Writer ───────────────────────────────────────────────────────────────────

def build_df(rows, columns=('Parameter', 'Symbol', 'Value', 'Units', 'Notes')):
    return pd.DataFrame(rows, columns=list(columns))


def write_excel():
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    prod_df  = build_df(PRODUCTION,
                        ('Parameter', 'Symbol', 'Fission', 'Fusion', 'Units'))
    ener_df  = build_df(ENERGETICS)
    diff_df  = build_df(DIFFUSION)
    diss_df  = build_df(DISSOCIATION)
    reac_df  = build_df(REACTIONS)

    with pd.ExcelWriter(OUTPUT_FILE, engine='openpyxl') as writer:
        prod_df.to_excel(writer, sheet_name='Production',   index=False)
        ener_df.to_excel(writer, sheet_name='Energetics',   index=False)
        diff_df.to_excel(writer, sheet_name='Diffusion',    index=False)
        diss_df.to_excel(writer, sheet_name='Dissociation', index=False)
        reac_df.to_excel(writer, sheet_name='Reactions',    index=False)

    print(f"Written: {OUTPUT_FILE.resolve()}")
    print(f"  Sheets: Production, Energetics, Diffusion, Dissociation, Reactions")


if __name__ == '__main__':
    write_excel()
