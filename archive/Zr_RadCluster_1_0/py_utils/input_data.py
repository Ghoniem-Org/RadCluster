"""
input_data.py — InputData class for RadCluster_1_0.

Reads material, irradiation, and model parameters from a 5-sheet Excel
workbook (RadCluster_1_0/input/input_parameters.xlsx) and computes
all derived quantities needed by the ODE system.

Sheet layout
------------
Production    — cascade production parameters (fission vs. fusion)
Energetics    — lattice constants, migration/formation energies, He EOS
Diffusion     — diffusion pre-factors, EUROFER solute trapping
Dissociation  — binding energy parameters, trap mutation barriers
Reactions     — rate constant prefactors, sink parameters, solver settings

Physics reference
-----------------
Ghoniem, N.M. (2026), "A Cluster Dynamics Model for Radiation Damage
Evolution in Ferritic-Martensitic Steels" (Rate_Equations.pdf).

Naming convention
-----------------
I           — max SIA cluster size
V           — max vacancy cluster size
i_mobile    — max mobile SIA cluster size
v_mobile    — max mobile vacancy cluster size
i_discrete  — max discrete SIA size (individually tracked; default = i_mobile)
v_discrete  — max discrete vacancy size (individually tracked; default = v_mobile)
I_bin       — number of SIA bin-moment equations beyond i_discrete
V_bin       — number of VAC bin-moment equations beyond v_discrete
i_cascade   — max SIA cluster size from cascade
v_cascade   — max vacancy cluster size from cascade

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

# Valid options
_SOLVER_MODES   = ('full_system', 'active_window')

# Backward-compat aliases for the previous solver-mode names.  Old configs,
# Excel files, and saved output dirs still reference the legacy strings; this
# map silently translates them to the new names with no warning, since the
# physics is identical.
_SOLVER_MODE_ALIASES = {
    'cpp_full':       'full_system',
    'sliding_OpenMP': 'active_window',
}
_PHYSICS_OPTIONS = ('full_CD_fission', 'full_CD_fusion',
                    'bin_moment_CD_fission', 'bin_moment_CD_fusion')

# Two-axis decomposition of physics_option.
_EQUATIONS = ('discrete', 'bin_moment')
_CASCADES  = ('fission', 'fusion')

# Legacy equations-axis aliases (silently accepted).  The internal
# physics_option strings still use the historical "full_CD" token; only the
# user-facing equations selector has been renamed.
_EQUATIONS_ALIASES = {
    'full_CD': 'discrete',
}


def make_physics_option(equations, cascade):
    """Combine an (equations, cascade) pair into the canonical physics_option string.

    equations : 'discrete' | 'bin_moment'   (legacy 'full_CD' also accepted)
    cascade   : 'fission'  | 'fusion'
    """
    equations = _EQUATIONS_ALIASES.get(equations, equations)
    if equations not in _EQUATIONS:
        raise ValueError(f"equations must be one of {_EQUATIONS}, got '{equations}'")
    if cascade not in _CASCADES:
        raise ValueError(f"cascade must be one of {_CASCADES}, got '{cascade}'")
    if equations == 'bin_moment':
        return f'bin_moment_CD_{cascade}'
    return f'full_CD_{cascade}'


def split_physics_option(po):
    """Reverse of make_physics_option: return (equations, cascade).

    The 'full_CD_*' physics_option prefix maps back to the user-facing
    equations name 'discrete'.
    """
    if po.startswith('bin_moment_CD_'):
        return ('bin_moment', po[len('bin_moment_CD_'):])
    if po.startswith('full_CD_'):
        return ('discrete', po[len('full_CD_'):])
    raise ValueError(f"Cannot split physics_option='{po}'")
_SHAPE_FUNCTIONS = ('constant', 'linear', 'lognormal')

# Backward-compat key aliases (old → new)
# NOTE: n1_bin is NOT aliased to i_discrete — different semantics.
# n1_bin was the starting bin edge; i_discrete is the number of discrete sizes.
_KEY_ALIASES = {
    'N':       'I',
    'M':       'V',
    'n_max_i': 'i_mobile',
    'm_max_v': 'v_mobile',
    'm1':      'i_cascade',
    'n1':      'v_cascade',
}


class InputData:
    """
    Material, irradiation, and model parameters for RadCluster_1_0.

    Reads from the 5-sheet Excel workbook.  Call display_parameters() to
    inspect all loaded values.

    Parameters
    ----------
    excel_file     : path-like, optional
    I              : int, optional  — override max SIA cluster size
    V              : int, optional  — override max vacancy cluster size
    solver_mode    : str, optional  — override solver mode
    physics_option : str, optional  — override physics option
    """

    def __init__(self, excel_file=INPUT_FILE, I=None, V=None,
                 solver_mode=None, physics_option=None,
                 # Backward-compat kwargs
                 N=None, M=None):
        # Support old N/M kwargs
        if I is None and N is not None:
            I = N
        if V is None and M is not None:
            V = M

        self.excel_file = Path(excel_file)
        if not self.excel_file.is_file():
            raise FileNotFoundError(
                f"Excel file not found: {self.excel_file}\n"
                f"Run RadCluster_1_0/create_excel.py to generate it."
            )
        print(f"Loading parameters from: {self.excel_file.resolve()}")
        self._load_data()

        # Apply caller overrides
        if I is not None:
            self.reactions['I'] = int(I)
        if V is not None:
            self.reactions['V'] = int(V)
        if solver_mode is not None:
            self.reactions['solver_mode'] = str(solver_mode)
        if physics_option is not None:
            self.reactions['physics_option'] = str(physics_option)

        self._calculate_derived()
        self._validate()

    # ── Sheet loading ─────────────────────────────────────────────────────────

    @staticmethod
    def _sheet_to_dict(df):
        """Convert Symbol/Value sheet to dict.  Skips header rows (no symbol)."""
        d = {}
        sym_col = None
        val_col = None
        for c in df.columns:
            cl = str(c).lower()
            if 'symbol' in cl:
                sym_col = c
            if 'value' in cl or 'fission' in cl:
                val_col = c
        if sym_col is None:
            sym_col = df.columns[1]
        if val_col is None:
            val_col = df.columns[2]
        for _, row in df.iterrows():
            k = row.get(sym_col, None)
            v = row.get(val_col, None)
            if pd.notna(k) and str(k).strip() and str(k).strip() != 'nan':
                d[str(k).strip()] = v
        return d

    @staticmethod
    def _production_to_dict(df):
        """Read Production sheet — returns fission and fusion sub-dicts."""
        fis, fus = {}, {}
        for _, row in df.iterrows():
            sym = row.get('Symbol', None)
            if not pd.notna(sym) or not str(sym).strip():
                continue
            k = str(sym).strip()
            f_val = row.get('Fission', None)
            u_val = row.get('Fusion', None)
            if pd.notna(f_val):
                fis[k] = f_val
            if pd.notna(u_val):
                fus[k] = u_val
        return fis, fus

    @staticmethod
    def _apply_aliases(d):
        """Translate old key names to new convention (in-place)."""
        for old, new in _KEY_ALIASES.items():
            if old in d and new not in d:
                d[new] = d.pop(old)

    def _load_data(self):
        """Read all five Excel worksheets."""
        try:
            prod_df  = pd.read_excel(self.excel_file, sheet_name='Production')
            ener_df  = pd.read_excel(self.excel_file, sheet_name='Energetics')
            diff_df  = pd.read_excel(self.excel_file, sheet_name='Diffusion')
            diss_df  = pd.read_excel(self.excel_file, sheet_name='Dissociation')
            reac_df  = pd.read_excel(self.excel_file, sheet_name='Reactions')
        except Exception as exc:
            raise RuntimeError(f"Failed to read Excel file: {exc}") from exc

        self.production_fission, self.production_fusion = \
            self._production_to_dict(prod_df)

        self.energetics  = self._sheet_to_dict(ener_df)
        self.diffusion   = self._sheet_to_dict(diff_df)
        self.dissociation= self._sheet_to_dict(diss_df)
        self.reactions   = self._sheet_to_dict(reac_df)

        # Apply backward-compat aliases (old Excel → new names)
        for d in (self.reactions, self.diffusion, self.production_fission,
                  self.production_fusion):
            self._apply_aliases(d)

        # Cast integer fields
        _int_keys = ('I', 'V', 'L_He_max', 'n_points', 'log_time',
                     'i_discrete', 'v_discrete', 'I_bin', 'V_bin',
                     'n_moments', 'n_group', 'window_width',
                     'v_mobile', 'i_mobile',
                     'i_cascade', 'v_cascade')
        for k in _int_keys:
            for d in (self.reactions, self.diffusion, self.production_fission,
                      self.production_fusion):
                if k in d:
                    try:
                        d[k] = int(float(d[k]))
                    except (TypeError, ValueError):
                        pass

        print("Successfully loaded all five parameter sheets.")

    # ── Derived quantities ────────────────────────────────────────────────────

    def _calculate_derived(self):
        """Compute all derived physics quantities."""
        e   = self.energetics
        d   = self.diffusion
        re  = self.reactions

        # Temperature and irradiation conditions
        T     = float(re.get('T', 600.0))
        G     = float(re.get('G', 1.0e-6))
        kBT   = _kB * T

        # Lattice
        a_nm  = float(e.get('a',     0.2867))
        a_m   = a_nm * 1.0e-9            # nm → m
        Omega = float(e.get('Omega',  1.18e-29))
        r0    = (3.0 * Omega / (4.0 * np.pi)) ** (1.0 / 3.0)
        b_111 = float(e.get('b_111', 0.2482)) * 1.0e-9   # nm → m

        # Energetics
        E_f_v  = float(e.get('E_f_v',  2.0))
        E_m_v  = float(e.get('E_m_v',  0.67))
        E_m_i  = float(e.get('E_m_i',  0.34))
        E_m_h  = float(e.get('E_m_h',  0.06))
        E_s_He = float(e.get('E_s_He', 2.35))
        gamma_s = float(e.get('gamma_s', 2.0))   # J/m^2

        # Attempt frequencies [s^-1]
        nu_v = float(e.get('nu_v', 1.0e13))
        nu_i = float(e.get('nu_i', 1.0e13))
        nu_h = float(e.get('nu_h', 3.0e12))

        # Pure Fe diffusivities [m^2/s]  — Eq. 17
        Dv_Fe = a_m**2 * nu_v * np.exp(-E_m_v / kBT)
        Di_Fe = a_m**2 * nu_i * np.exp(-E_m_i / kBT)
        Dh_Fe = a_m**2 * nu_h * np.exp(-E_m_h / kBT)

        # Jump frequencies [s^-1]  — Eq. 22-24
        omega_v_Fe = Dv_Fe / a_m**2
        omega_i_Fe = Di_Fe / a_m**2
        omega_h_Fe = Dh_Fe / a_m**2

        # EUROFER solute trapping (Eq. 42, 48)
        c_Cr  = float(d.get('c_Cr',  0.094))
        c_W   = float(d.get('c_W',   0.0033))
        c_Mn  = float(d.get('c_Mn',  0.0047))
        c_C   = float(d.get('c_C',   5.0e-4))
        c_N   = float(d.get('c_N',   2.0e-4))

        # SIA trapping (Eq. 42)
        def _trap_sum_SIA():
            # z_s · c_s · exp(E_b^{s,i} / kBT)
            E_b_C_SIA  = float(d.get('E_b_C_SIA',  0.45)); z_C_SIA  = 4
            E_b_N_SIA  = float(d.get('E_b_N_SIA',  0.40)); z_N_SIA  = 4
            E_b_Cr_SIA = float(d.get('E_b_Cr_SIA', 0.10)); z_Cr_SIA = 8
            E_b_Mn_SIA = float(d.get('E_b_Mn_SIA', 0.20)); z_Mn_SIA = 6
            return (z_C_SIA  * c_C  * np.exp(E_b_C_SIA  / kBT) +
                    z_N_SIA  * c_N  * np.exp(E_b_N_SIA  / kBT) +
                    z_Cr_SIA * c_Cr * np.exp(E_b_Cr_SIA / kBT) +
                    z_Mn_SIA * c_Mn * np.exp(E_b_Mn_SIA / kBT))

        # Vacancy trapping (Eq. 48)
        def _trap_sum_VAC():
            E_b_C_V   = float(d.get('E_b_C_V',   0.45)); z_C_V   = 3
            E_b_N_V   = float(d.get('E_b_N_V',   0.40)); z_N_V   = 3
            E_b_W_V   = float(d.get('E_b_W_V',   0.27)); z_W_V   = 8
            E_b_Mn_V  = float(d.get('E_b_Mn_V',  0.10)); z_Mn_V  = 8
            E_b_Cr_V  = float(d.get('E_b_Cr_V',  0.05)); z_Cr_V  = 8
            return (z_C_V   * c_C  * np.exp(E_b_C_V   / kBT) +
                    z_N_V   * c_N  * np.exp(E_b_N_V   / kBT) +
                    z_W_V   * c_W  * np.exp(E_b_W_V   / kBT) +
                    z_Mn_V  * c_Mn * np.exp(E_b_Mn_V  / kBT) +
                    z_Cr_V  * c_Cr * np.exp(E_b_Cr_V  / kBT))

        # SIA cluster loop trapping (Eq. 52)
        def _trap_sum_loop():
            E_b_C_loop  = float(d.get('E_b_C_loop',  0.50)); z_C_loop  = 2
            E_b_N_loop  = float(d.get('E_b_N_loop',  0.40)); z_N_loop  = 2
            E_b_Cr_loop = float(d.get('E_b_Cr_loop', 0.10)); z_Cr_loop = 4
            return (z_C_loop  * c_C  * np.exp(E_b_C_loop  / kBT) +
                    z_N_loop  * c_N  * np.exp(E_b_N_loop  / kBT) +
                    z_Cr_loop * c_Cr * np.exp(E_b_Cr_loop / kBT))

        trap_SIA  = _trap_sum_SIA()
        trap_VAC  = _trap_sum_VAC()
        trap_loop = _trap_sum_loop()

        # Effective diffusivities and jump frequencies (Eqs. 42, 48, 52)
        omega_i_eff = omega_i_Fe / (1.0 + trap_SIA)
        omega_v_eff = omega_v_Fe / (1.0 + trap_VAC)
        omega_h_eff = omega_h_Fe                       # He not trapped by solutes

        Di_eff = omega_i_eff * a_m**2
        Dv_eff = omega_v_eff * a_m**2
        Dh_eff = omega_h_eff * a_m**2

        # SIA cluster 1D glide (Eq. 33)
        nu0_1D  = float(d.get('nu0_1D', 6.0e12))
        E_m_1D  = float(d.get('E_m_1D', 0.03))
        s_1D    = float(d.get('s_1D',   0.7))
        i_mobile = int(float(d.get('i_mobile', 1)))
        v_mobile = int(float(d.get('v_mobile', 1)))

        # Boundary flux option for coalescence at upper size limit
        # 'absorption' (default): product lost at boundary (open boundary)
        # 'reflection': product folded back into largest tracked size (closed)
        boundary_flux = str(re.get('boundary_flux',
                            d.get('boundary_flux', 'absorption'))).lower().strip()
        if boundary_flux not in ('absorption', 'reflection'):
            import warnings
            warnings.warn(f"Unknown boundary_flux='{boundary_flux}', using 'absorption'")
            boundary_flux = 'absorption'

        # D_n^{1D}(n) = (3a²ν_0^{1D}) / (2n^{s_1D}) · exp(−E_m^{1D}/k_BT)  (Eq. 33)
        D1D_base = (3.0 * a_m**2 * nu0_1D / 2.0) * np.exp(-E_m_1D / kBT)

        # Loop trapping correction for 1D glide (Eq. 52)
        def D1D(n):
            return D1D_base / float(n)**s_1D / (1.0 + trap_loop)

        # Mean free path for 1D/3D mixed (Eq. 121)
        L_hat = float(d.get('L_hat', 50.0))   # L/a (dimensionless)
        B_rot = float(d.get('B_rot', 2.627))

        # Equilibrium vacancy concentration
        Cv_eq = np.exp(-E_f_v / kBT)

        # Geometric rate constant prefactors (Eq. 128)
        A_sph  = (48.0 * np.pi**2)**(1.0/3.0)                   # ≈ 7.818
        A_loop = 8.0 * np.sqrt(np.pi / np.sqrt(3.0))            # ≈ 10.78
        A_1D   = 9.0 / (8.0 * np.pi**(2.0/3.0))                 # ≈ 2.632

        # He production
        spectrum = str(re.get('spectrum', 'fission')).lower()
        from .defect_production import FISSION, FUSION
        spec = FISSION if 'fiss' in spectrum else FUSION
        G_He_r = float(re.get('G_He_r', spec['G_He_r']))
        G_He   = G_He_r * 1.0e-6 * G    # appm/dpa * dpa/s → atom frac/s

        self.derived = {
            'T':           T,
            'G':           G,
            'kBT':         kBT,
            'a_m':         a_m,
            'Omega':       Omega,
            'r0':          r0,
            'b_111':       b_111,
            'E_f_v':       E_f_v,
            'E_m_v':       E_m_v,
            'E_m_i':       E_m_i,
            'E_m_h':       E_m_h,
            'E_s_He':      E_s_He,
            'gamma_s':     gamma_s,
            'nu_v':        nu_v,
            'nu_i':        nu_i,
            'nu_h':        nu_h,
            'Di_Fe':       Di_Fe,
            'Dv_Fe':       Dv_Fe,
            'Dh_Fe':       Dh_Fe,
            'Di_eff':      Di_eff,
            'Dv_eff':      Dv_eff,
            'Dh_eff':      Dh_eff,
            'omega_i_eff': omega_i_eff,
            'omega_v_eff': omega_v_eff,
            'omega_h_eff': omega_h_eff,
            'trap_SIA':    trap_SIA,
            'trap_VAC':    trap_VAC,
            'trap_loop':   trap_loop,
            'D1D_base':    D1D_base,
            'D1D':         D1D,
            's_1D':        s_1D,
            'i_mobile':    i_mobile,
            'v_mobile':    v_mobile,
            'L_hat':       L_hat,
            'B_rot':       B_rot,
            'Cv_eq':       Cv_eq,
            'A_sph':       A_sph,
            'A_loop':      A_loop,
            'A_1D':        A_1D,
            'G_He':        G_He,
            'G_He_r':      G_He_r,
            'spectrum':    spectrum,
            'boundary_flux': boundary_flux,
        }

        print(f"Derived: T={T} K  Cv_eq={Cv_eq:.3e}"
              f"  Di_eff={Di_eff:.3e}  Dv_eff={Dv_eff:.3e} m2/s"
              f"  spectrum='{spectrum}'")

    # ── Validation ────────────────────────────────────────────────────────────

    def _validate(self):
        T   = self.derived['T']
        G   = self.derived['G']
        rho = float(self.reactions.get('rho_d', 1.0e14))

        if not (300 <= T <= 1200):
            warnings.warn(f"Temperature {T} K is outside typical range [300–1200 K]")
        if not (1e-9 <= G <= 1e-3):
            warnings.warn(f"Dose rate {G} dpa/s is outside [1e-9–1e-3]")
        if not (1e12 <= rho <= 1e16):
            warnings.warn(f"Dislocation density {rho} m^-2 outside typical range")

        sm = self.solver_mode
        if sm in _SOLVER_MODE_ALIASES:
            self.reactions['solver_mode'] = _SOLVER_MODE_ALIASES[sm]
            sm = self.reactions['solver_mode']
        if sm not in _SOLVER_MODES:
            warnings.warn(f"Unknown solver_mode='{sm}'. Using 'full_system'.")
            self.reactions['solver_mode'] = 'full_system'

        po = self.physics_option
        if po not in _PHYSICS_OPTIONS:
            warnings.warn(f"Unknown physics_option='{po}'. Using 'full_CD_fission'.")
            self.reactions['physics_option'] = 'full_CD_fission'

        sf = self.shape_function
        if sf not in _SHAPE_FUNCTIONS:
            warnings.warn(f"Unknown shape_function='{sf}'. Using 'linear'.")
            self.reactions['shape_function'] = 'linear'

    # ── Convenience properties ────────────────────────────────────────────────

    @property
    def I(self):
        """Max SIA cluster size."""
        return int(float(self.reactions.get('I', 500)))

    @property
    def V(self):
        """Max vacancy cluster size."""
        return int(float(self.reactions.get('V', 500)))

    # Backward-compat aliases
    @property
    def N(self):
        return self.I

    @N.setter
    def N(self, val):
        self.reactions['I'] = int(val)

    @property
    def M(self):
        return self.V

    @M.setter
    def M(self, val):
        self.reactions['V'] = int(val)

    @property
    def i_mobile(self):
        """Max mobile SIA cluster size."""
        return self.derived.get('i_mobile', 1)

    @property
    def v_mobile(self):
        """Max mobile vacancy cluster size."""
        return self.derived.get('v_mobile', 1)

    @property
    def i_discrete(self):
        """Max discrete SIA size (individually tracked). Defaults to i_mobile."""
        val = self.reactions.get('i_discrete', None)
        if val is not None:
            return int(float(val))
        return self.i_mobile

    @property
    def v_discrete(self):
        """Max discrete vacancy size (individually tracked). Defaults to v_mobile."""
        val = self.reactions.get('v_discrete', None)
        if val is not None:
            return int(float(val))
        return self.v_mobile

    @property
    def I_bin(self):
        """Number of SIA bin-moment equations beyond i_discrete.
        Defaults to ceil(log(I/i_discrete)/log(2)) if not specified."""
        val = self.reactions.get('I_bin', None)
        if val is not None:
            return int(float(val))
        # Default: auto-compute from r=2.0
        id = self.i_discrete
        if id >= self.I:
            return 0
        return int(np.ceil(np.log(self.I / max(id, 1)) / np.log(2.0)))

    @property
    def V_bin(self):
        """Number of VAC bin-moment equations beyond v_discrete.
        Defaults to ceil(log(V/v_discrete)/log(2)) if not specified."""
        val = self.reactions.get('V_bin', None)
        if val is not None:
            return int(float(val))
        # Default: auto-compute from r=2.0
        vd = self.v_discrete
        if vd >= self.V:
            return 0
        return int(np.ceil(np.log(self.V / max(vd, 1)) / np.log(2.0)))

    @property
    def L_He_max(self):
        val = self.reactions.get('L_He_max', 'mf')
        if str(val).lower() in ('mf', 'mean-field', 'nan', 'none', ''):
            return None   # signal to use mean-field reduction
        return int(float(val))

    @property
    def solver_mode(self):
        sm = str(self.reactions.get('solver_mode', 'full_system')).strip()
        return _SOLVER_MODE_ALIASES.get(sm, sm)

    @property
    def physics_option(self):
        return str(self.reactions.get('physics_option', 'full_CD_fission')).strip()

    @property
    def equations(self):
        """Equation system: 'full_CD' or 'bin_moment' (derived from physics_option)."""
        return split_physics_option(self.physics_option)[0]

    @property
    def cascade(self):
        """Cascade type: 'fission' or 'fusion' (derived from physics_option)."""
        return split_physics_option(self.physics_option)[1]

    @property
    def alpha_He(self):
        return float(self.reactions.get('alpha_He', 1.7))

    @property
    def shape_function(self):
        """Intra-bin shape function: 'constant', 'linear', or 'lognormal'."""
        return str(self.reactions.get('shape_function', 'linear')).strip().lower()

    # ── Display ───────────────────────────────────────────────────────────────

    def display_parameters(self):
        sections = [
            ('ENERGETICS',  self.energetics),
            ('DIFFUSION',   self.diffusion),
            ('REACTIONS',   self.reactions),
            ('DERIVED',     {k: v for k, v in self.derived.items()
                             if not callable(v)}),
        ]
        for title, d in sections:
            print(f"\n{'='*60}\n{title}\n{'='*60}")
            for k, v in d.items():
                fmt = f"  {k}: {v:.4e}" if isinstance(v, float) else f"  {k}: {v}"
                print(fmt)
