/**
 * parameters.h — RadCluster_2_0 solver parameter struct.
 *
 * All quantities that the ODE right-hand side needs are pre-computed on the
 * Python side (InputData + ReactionRates) and forwarded via a parameter file
 * (--param_file=<path>). build_parameters() unpacks them into this struct so
 * the RHS performs only arithmetic.
 *
 * Physics reference: Ghoniem (2026), Rate_Equations.pdf.
 *
 * State vector y[N_eq]
 * --------------------
 * full_CD_fission  (Case 2, he_mode=0):
 *   y[0..I-1]       — SIA clusters c_n, n=1..I
 *   y[I..I+V-1]     — vacancy/bubble c_m (marginal), m=1..V
 *   y[I+V]          — Q_tot (total He in voids)
 *   y[I+V+1]        — c_h (free He)  [omitted when he_kinetics=1 (QSS)]
 *   N_eq = I + V + 2  (dynamic)  or  I + V + 1  (quasi_steady_state)
 *
 * full_CD_fusion  (Case 1, he_mode=1):
 *   y[0..I-1]       — SIA clusters c_n
 *   y[I..I+V-1]     — c_m^tot
 *   y[I+V..I+2V-1]  — Q_m (He per class)
 *   y[I+2V]         — c_h  [omitted when he_kinetics=1 (QSS)]
 *   N_eq = I + 2V + 1  (dynamic)  or  I + 2V  (quasi_steady_state)
 *
 * bin_moment_CD_fission/fusion (physics_option 2/3):
 *   y[0..2Ib-1]     — SIA bin moments [μ_0^(0), μ_0^(1), ..., μ_{Ib-1}^(1)]
 *   y[2Ib..2Ib+V-1] — c_m
 *   y[2Ib+V..]      — He variables (same as Case 2 or Case 1)
 *   N_eq = 2Ib + V + 2  or  2Ib + 2V + 1  (dynamic)
 *   N_eq = 2Ib + V + 1  or  2Ib + 2V      (quasi_steady_state)
 *
 * he_kinetics: 0 = dynamic (c_h is an ODE state, Eq. 157)
 *              1 = quasi_steady_state (c_h computed algebraically from
 *                  dc_h/dt = 0;  valid because E_m_h = 0.06 eV is small)
 *
 * C_floor: concentration floor.  In the RHS, any state variable below C_floor
 *   is clamped to C_floor for rate computation, and any derivative of a
 *   floor-clamped variable is clamped to >= 0 to prevent negative excursions.
 */
#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

struct Parameters {
    // ── Cluster size limits ─────────────────────────────────────────────────
    int I, V, N_eq;
    int Ni_max;   // = I (kept for legacy compatibility)

    // ── Physics and solver mode ─────────────────────────────────────────────
    // physics_option: 0=full_CD_fission, 1=full_CD_fusion,
    //                  2=bin_moment_CD_fission, 3=bin_moment_CD_fusion
    // he_mode:         0=Case 2 (decoupled/fission), 1=Case 1 (mean-field/fusion)
    // he_kinetics:     0=dynamic (ODE), 1=quasi_steady_state (algebraic c_h)
    int physics_option;
    int he_mode;
    int he_kinetics;  // 0=dynamic, 1=quasi_steady_state

    // ── Geometric rate constant prefactors (Eq. 128) ─────────────────────────
    double A_sph;   // (48π²)^{1/3} ≈ 7.818
    double A_loop;  // 8√(π/√3) ≈ 10.78
    double B_rot;   // (4/π)(8π/3)^{1/3} ≈ 2.627
    double L_hat;   // mean free path L/a (dimensionless)

    // ── Solute trapping sums (Eq. 42, 48, 52) ────────────────────────────────
    double trap_SIA;   // Σ z_s·c_s·exp(E_b^{s,i}/kBT)
    double trap_VAC;   // Σ z_s·c_s·exp(E_b^{s,v}/kBT)
    double trap_loop;  // Σ z_s·c_s·exp(E_b^{s,loop}/kBT)

    // ── Mobility cutoffs ──────────────────────────────────────────────────────
    int i_mobile;   // max mobile SIA cluster size (1D glide cutoff)
    int v_mobile;   // max mobile vacancy cluster size

    // ── Boundary flux option ─────────────────────────────────────────────────
    // 0 = absorption (open boundary): product at I+1 or V+1 is lost
    // 1 = reflection (closed boundary): product folded back into I or V
    int boundary_flux;

    // ── Vacancy cluster arrays [V] ────────────────────────────────────────────
    std::vector<double> KVV;      // Cv capture by void m  (K_VAC_grow)
    std::vector<double> KVI;      // Ci annihilation at void m  (K_VAC_shrink)
    std::vector<double> GVV;      // thermal vacancy emission from void m
    std::vector<double> KHeV;     // He capture by void m
    std::vector<double> Pr_VAC;   // cascade vacancy production rate [V]
    std::vector<double> m13;      // m^{1/3} factors [V]

    // ── SIA cluster arrays [I] ────────────────────────────────────────────────
    std::vector<double> KII;      // Ci capture by SIA loop n  (K_SIA_grow)
    std::vector<double> KIV;      // Cv capture by SIA loop n  (K_SIA_shrink)
    std::vector<double> GII;      // thermal SIA emission from loop n
    std::vector<double> k2_SIA;   // total sink rate for SIA cluster n
    std::vector<double> Pr_SIA;   // cascade SIA production rate [I]

    // ── 1D glide prefactors (Eq. 141) ─────────────────────────────────────────
    // K_1D_eff(n, m) = K_1D_pref[n-1] · m^{1/3} / (1 + B_rot·L̂²·m^{-1/3})
    std::vector<double> K_1D_pref;   // [I]: A_sph · ω_n^{1D} / Ω

    // Legacy separable cross-term (for backward compatibility with solver.cpp)
    std::vector<double> K_IclV_ns;   // [I]: 4π·r0·Di·n^{-2/3} / Ω
    std::vector<double> K_IclV_ni;   // [I]: 4π·r0·Di / (n·Ω)

    // ── Mobile cluster effective 3D diffusivities (for coalescence) ─────────
    std::vector<double> D_SIA_eff;   // [I]: effective 3D D for SIA cluster n
    std::vector<double> D_VAC_eff;   // [V]: effective 3D D for vac cluster m
    double A_sph_inv_O23;            // A_sph / Ω^{2/3}  [m^-2]
    double A_loop_inv_O23;           // A_loop / Ω^{2/3} [m^-2]  (loop geometry for n≥4)
    double Z_i_loop;                 // SIA dislocation bias at loops (Table 26, Eq. P3_i)
    double Z_ii;                     // SIA-SIA coalescence bias factor (elastic interaction)

    // ── Scalar physics ────────────────────────────────────────────────────────
    double G_He;          // He transmutation rate [at.frac/s]
    double k2_disl_v;     // vacancy fixed sink [s^-1]
    double k2_disl_i;     // SIA fixed sink [s^-1]
    double k2_disl_He;    // He fixed sink [s^-1]
    double Cv_eq;         // thermal vacancy equilibrium concentration
    double K_iv;          // V–SIA recombination prefactor (P1, mutual diffusivity)
    double K_3D_cav_pref; // 3D cavity absorption prefactor: A_sph · Di_eff / Ω^{2/3}

    // ── He parameters ─────────────────────────────────────────────────────────
    double beta_He;       // He emission: ν_h·exp(−(E_b_hV+E_m_h)/kBT) [s^-1]
    double delta_He;      // He pressure coefficient [eV]
    double beta_He_exp;   // He pressure power-law exponent
    double kBT;           // k_B·T [eV]

    // ── Bin-moment parameters (Chapter 9) ────────────────────────────────────
    int    I_bin;         // number of SIA logarithmic bins
    int    V_bin;         // number of vacancy logarithmic bins
    double r_ratio;       // bin ratio r > 1 (Eq. 188) — diagnostic only;
                          //   bin edges are transmitted explicitly (see below)
    int    i_discrete;    // max discrete SIA size (minimum cluster size in binned region)
    int    v_discrete;    // max discrete vacancy size (default = v_mobile)
    int    shape_function; // intra-bin closure: 0=constant, 1=linear, 2=lognormal
    int    n_mom;          // moments per bin: 1 (constant), 2 (linear), 3 (lognormal)
                           //   DERIVED as shape_function + 1 (not read independently)

    // Explicit integer bin edges transmitted from Python (cpp_bridge.py).
    // Bin k covers integer sizes [bin_lo[k], bin_hi[k]-1].  These are the
    // exact partition Python computed; re-deriving them in C++ from r_ratio
    // via std::floor can diverge from numpy.floor over many bins.
    std::vector<int> sia_bin_lo;   // [I_bin]
    std::vector<int> sia_bin_hi;   // [I_bin]
    std::vector<int> vac_bin_lo;   // [V_bin]
    std::vector<int> vac_bin_hi;   // [V_bin]

    // ── Initial conditions ─────────────────────────────────────────────────────
    std::vector<double> y0;   // [N_eq]

    // ── Concentration floor ────────────────────────────────────────────────────
    double C_floor;

    // ── Solver settings ────────────────────────────────────────────────────────
    double t_begin;
    double t_end;
    int    n_points;
    bool   log_time;
    double rtol;
    double atol;

    // ── Integration method ─────────────────────────────────────────────────────
    // CVODE BDF is the only integrator wired up; the linear solver is selectable.
    int linsol;     // 0=dense, 1=band, 2=gmres, 3=klu (sparse direct, full_CD only)
    int mu, ml;     // band solver bandwidths
    int max_order;  // 0 = solver default
    double hmin;    // minimum step size (0 = no limit)
    double hmax;    // maximum step size (0 = no limit)

    // ── Dynamic window (active_window) ────────────────────────────────────────
    // window_mode: 0=full_system, 4=active_window (two independent sliding windows + OpenMP-parallel RHS)
    // Two independent windows: one for SIA cluster indices, one for VAC cluster indices.
    int    window_mode;
    int    window_width;             // initial window width (shared by SIA and VAC)
    double concentration_threshold;  // window expansion threshold (shared by SIA and VAC)
    int    window_pad;               // SIA window expansion pad
    int    window_pad_v;             // VAC window expansion pad (default = window_pad)
    double window_expand_factor;
    int    window_check_every;
    double window_C_contract;
    int    window_min_active_i;
    double window_nuc_guard;
    double window_t_start;
    double Ni_extend_tol;
    int    Ni_extend_margin;

    // ── Preconditioner storage ───────────────────────────────────────────────
    // Jacobi diagonal (used when prec_type==0)
    std::vector<double> prec_diag;

    // Woodbury bordered-banded preconditioner (used when prec_type==1)
    //
    // The Jacobian has the structure  J = T + U·V^T  where:
    //   T is banded (half-bandwidth prec_bw)
    //   U is N_eq × prec_rank  (dense columns from mobile species)
    //   V = [e_{j1}, ..., e_{jr}]  (selector for mobile indices)
    //
    // The preconditioner solves  (I - γJ)x = r  via SMW:
    //   M = T̂ - γ U V^T,  T̂ = I - γT  (banded)
    //   M^{-1} = T̂^{-1} + T̂^{-1} U S^{-1} V^T T̂^{-1}
    //   S = -I/γ_scale + V^T T̂^{-1} U   (r × r Schur complement)
    int prec_type;       // 0=Jacobi (legacy), 1=Woodbury (default when coalescence)
    int prec_bw;         // half-bandwidth of T  (auto: max(2*i_mobile, 2*v_mobile) + 1)
    int prec_rank;       // rank of dense border  (auto: i_mobile + v_mobile)

    // Storage (allocated in prec_setup):
    //   prec_band:      banded LU of T̂ in LAPACK dgbtrf layout
    //                   rows: 2*prec_bw + prec_bw + 1 = 3*prec_bw + 1
    //                   cols: N_eq
    //   prec_Tinv_U:    T̂^{-1} U  [N_eq × prec_rank]
    //   prec_schur:     S factored  [prec_rank × prec_rank]
    //   prec_ipiv_band: pivot array for dgbtrf [N_eq]
    //   prec_ipiv_schur:pivot array for dgetrf [prec_rank]
    //   prec_mobile_idx:indices of mobile species in state vector [prec_rank]
    //   prec_f0:        base RHS evaluation [N_eq]
    std::vector<double> prec_band;
    std::vector<double> prec_Tinv_U;
    std::vector<double> prec_schur;
    std::vector<int>    prec_ipiv_band;
    std::vector<int>    prec_ipiv_schur;
    std::vector<int>    prec_mobile_idx;
    std::vector<double> prec_f0;
    std::vector<double> prec_work;     // scratch [N_eq]
    std::vector<double> prec_y_save;   // saved y during FD probing [N_eq]
    std::vector<double> prec_deltas;   // FD perturbation sizes [N_eq]
    std::vector<double> prec_f_pert;   // perturbed RHS [N_eq]
    double              prec_gamma;    // cached γ from last setup

    // ── KLU sparse direct solver (linsol == 3) ───────────────────────────────
    // CSC sparsity pattern of the Jacobian (built once at startup). Values
    // are filled by sparse_fd_jac() via colored finite differences.
    //   jac_col_ptr      [N_eq + 1]: CSC column pointers
    //   jac_row_idx      [nnz]:      CSC row indices, sorted within each column
    //   jac_colors       [N_eq]:     CPR color of each column
    //   jac_n_colors:               number of color groups
    //   jac_color_groups [N_eq]:    columns flattened by color
    //   jac_color_offsets[n_colors+1]: offsets into jac_color_groups
    std::vector<int> jac_col_ptr;
    std::vector<int> jac_row_idx;
    std::vector<int> jac_colors;
    std::vector<int> jac_color_groups;
    std::vector<int> jac_color_offsets;
    int              jac_n_colors;

    // ── Diagnostics ────────────────────────────────────────────────────────────
    bool verbose;   // if false, suppress per-timestep [diag] / [ci5_rates] output
};

// ── CLI / file parsing helpers ───────────────────────────────────────────────

inline double require_param(const std::map<std::string, double>& m,
                             const std::string& key) {
    auto it = m.find(key);
    if (it == m.end()) {
        std::cerr << "Missing required parameter: " << key << "\n";
        exit(1);
    }
    return it->second;
}

inline double optional_param(const std::map<std::string, double>& m,
                              const std::string& key, double def) {
    auto it = m.find(key);
    return (it != m.end()) ? it->second : def;
}

/**
 * Read a parameter file written by cpp_bridge.py.
 * Format: key=value per line; blank lines and '#' lines ignored.
 */
inline std::map<std::string, double> parse_param_file(const std::string& path) {
    std::map<std::string, double> props;
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "Cannot open parameter file: " << path << "\n";
        exit(1);
    }
    std::string line;
    int lineno = 0;
    while (std::getline(f, line)) {
        ++lineno;
        if (line.empty() || line[0] == '#') continue;
        auto pos = line.find('=');
        if (pos == std::string::npos) continue;
        std::string key = line.substr(0, pos);
        double val = 0.0;
        try {
            val = std::stod(line.substr(pos + 1));
        } catch (...) {
            std::cerr << "Line " << lineno << ": bad value for \"" << key << "\"\n";
            exit(1);
        }
        props[key] = val;
    }
    return props;
}

/**
 * Unpack a parsed parameter map into a Parameters struct.
 */
inline Parameters build_parameters(const std::map<std::string, double>& p) {
    Parameters P{};

    // Cluster sizes — support old (Ni/Nv, N/M) and new (I/V) key names
    P.I = static_cast<int>(optional_param(p, "I",
               optional_param(p, "N",
                   optional_param(p, "Ni", -1))));
    if (P.I < 0) { std::cerr << "Missing required parameter: I (or N or Ni)\n"; exit(1); }
    P.V = static_cast<int>(optional_param(p, "V",
               optional_param(p, "M",
                   optional_param(p, "Nv", -1))));
    if (P.V < 0) { std::cerr << "Missing required parameter: V (or M or Nv)\n"; exit(1); }

    P.Ni_max = static_cast<int>(optional_param(p, "Ni_max", static_cast<double>(P.I)));

    // Physics mode
    P.physics_option = static_cast<int>(optional_param(p, "physics_option_int", 0));
    P.he_mode        = static_cast<int>(optional_param(p, "he_mode",             0));
    P.he_kinetics    = static_cast<int>(optional_param(p, "qss_He",              0));
    // he_kinetics: 0=dynamic (c_h is ODE state), 1=quasi_steady_state (c_h algebraic)

    // he_mode duplicates information already encoded in physics_option:
    //   0=full_CD_fission(Case2), 1=full_CD_fusion(Case1),
    //   2=bin_moment_fission(Case2), 3=bin_moment_fusion(Case1)
    // so the He-reduction case is physics_option % 2.  Fail loudly on mismatch
    // rather than silently running the wrong He physics.
    if (P.he_mode != (P.physics_option % 2)) {
        std::cerr << "Inconsistent contract: he_mode=" << P.he_mode
                  << " but physics_option=" << P.physics_option
                  << " requires he_mode=" << (P.physics_option % 2) << "\n";
        exit(1);
    }

    // State vector size depends on he_mode and he_kinetics:
    //   Case 2 dynamic:           I + V + 2
    //   Case 2 quasi_steady_state: I + V + 1  (c_h removed from state)
    //   Case 1 dynamic:           I + 2V + 1
    //   Case 1 quasi_steady_state: I + 2V     (c_h removed from state)
    const bool qss = (P.he_kinetics == 1);
    int n_he_extra;
    if (P.he_mode == 1)
        n_he_extra = qss ? P.V     : P.V + 1;
    else
        n_he_extra = qss ? 1       : 2;
    P.N_eq = P.I + P.V + n_he_extra + 5;  // +5 for conservation accounting ODEs

    // Geometric prefactors
    P.A_sph  = optional_param(p, "A_sph",  7.818);
    P.A_loop = optional_param(p, "A_loop", 10.78);
    P.B_rot  = optional_param(p, "B_rot",  2.627);
    P.L_hat  = optional_param(p, "L_hat",  50.0);

    // Solute trapping
    P.trap_SIA  = optional_param(p, "trap_SIA",  0.0);
    P.trap_VAC  = optional_param(p, "trap_VAC",  0.0);
    P.trap_loop = optional_param(p, "trap_loop", 0.0);

    // Mobility cutoffs — support both old and new key names
    P.i_mobile = static_cast<int>(optional_param(p, "i_mobile",
                     optional_param(p, "n_max_i", 100.0)));
    P.v_mobile = static_cast<int>(optional_param(p, "v_mobile",
                     optional_param(p, "m_max_v", 5.0)));

    // Boundary flux: 0=absorption (default), 1=reflection
    P.boundary_flux = static_cast<int>(optional_param(p, "boundary_flux", 0.0));

    // Resize arrays
    P.KVV.resize(P.V);   P.KVI.resize(P.V);   P.GVV.resize(P.V);
    P.KHeV.resize(P.V);  P.Pr_VAC.resize(P.V); P.m13.resize(P.V);
    P.KII.resize(P.I);   P.KIV.resize(P.I);    P.GII.resize(P.I);
    P.k2_SIA.resize(P.I); P.Pr_SIA.resize(P.I);
    P.K_1D_pref.resize(P.I);
    P.K_IclV_ns.resize(P.I); P.K_IclV_ni.resize(P.I);
    P.D_SIA_eff.resize(P.I); P.D_VAC_eff.resize(P.V);
    P.y0.resize(P.N_eq);

    // Vacancy arrays
    for (int k = 0; k < P.V; ++k) {
        P.KVV[k]    = require_param(p, "KVV_"    + std::to_string(k));
        P.KVI[k]    = require_param(p, "KVI_"    + std::to_string(k));
        P.GVV[k]    = require_param(p, "GVV_"    + std::to_string(k));
        P.KHeV[k]   = require_param(p, "KHeV_"   + std::to_string(k));
        P.Pr_VAC[k] = require_param(p, "Pr_VAC_" + std::to_string(k));
        P.m13[k]    = require_param(p, "m13_"    + std::to_string(k));
    }

    // SIA arrays
    for (int k = 0; k < P.I; ++k) {
        P.KII[k]         = require_param(p, "KII_"        + std::to_string(k));
        P.KIV[k]         = require_param(p, "KIV_"        + std::to_string(k));
        P.GII[k]         = require_param(p, "GII_"        + std::to_string(k));
        P.k2_SIA[k]      = require_param(p, "k2_SIA_"     + std::to_string(k));
        P.Pr_SIA[k]      = require_param(p, "Pr_SIA_"     + std::to_string(k));
        P.K_1D_pref[k]   = optional_param(p, "K_1D_pref_" + std::to_string(k), 0.0);
        P.K_IclV_ns[k]   = optional_param(p, "K_IclV_ns_" + std::to_string(k), 0.0);
        P.K_IclV_ni[k]   = optional_param(p, "K_IclV_ni_" + std::to_string(k), 0.0);
    }

    // Mobile cluster diffusivities (for coalescence)
    for (int k = 0; k < P.I; ++k)
        P.D_SIA_eff[k] = optional_param(p, "D_SIA_eff_" + std::to_string(k), 0.0);
    for (int k = 0; k < P.V; ++k)
        P.D_VAC_eff[k] = optional_param(p, "D_VAC_eff_" + std::to_string(k), 0.0);
    P.A_sph_inv_O23  = optional_param(p, "A_sph_inv_O23", 0.0);
    P.A_loop_inv_O23 = optional_param(p, "A_loop_inv_O23", 0.0);
    P.Z_i_loop       = optional_param(p, "Z_i_loop", 1.05);
    P.Z_ii           = optional_param(p, "Z_ii", 1.0);

    // Scalar physics
    P.G_He       = require_param(p, "G_He");
    P.k2_disl_v  = require_param(p, "k2_disl_v");
    P.k2_disl_i  = require_param(p, "k2_disl_i");
    P.k2_disl_He = require_param(p, "k2_disl_He");
    P.Cv_eq      = require_param(p, "Cv_eq");
    // K_iv (recombination) and K_3D_cav_pref (cavity absorption) are required:
    // a silent default of 0.0 would disable those physics entirely.
    P.K_iv           = require_param(p, "K_iv");
    P.K_3D_cav_pref  = require_param(p, "K_3D_cav_pref");
    P.beta_He    = require_param(p, "beta_He");
    P.delta_He   = require_param(p, "delta_He");
    P.beta_He_exp= require_param(p, "beta_He_exp");
    P.kBT        = require_param(p, "kBT");

    // Bin-moment parameters — support both old and new key names
    P.I_bin    = static_cast<int>(optional_param(p, "I_bin",
                     optional_param(p, "K_bins", 0.0)));
    P.V_bin    = static_cast<int>(optional_param(p, "V_bin",
                     optional_param(p, "K_v_bins", 0.0)));
    P.r_ratio  = optional_param(p, "r_ratio", 2.0);
    P.i_discrete = static_cast<int>(optional_param(p, "i_discrete",
                       static_cast<double>(P.i_mobile)));
    P.v_discrete = static_cast<int>(optional_param(p, "v_discrete",
                       static_cast<double>(P.v_mobile)));
    P.shape_function = static_cast<int>(optional_param(p, "shape_function", 1.0));
    // n_mom is DERIVED from shape_function (it duplicates that information):
    //   constant(0)->1 moment, linear(1)->2 moments, lognormal(2)->3 moments.
    P.n_mom = P.shape_function + 1;

    // Explicit integer bin edges transmitted from Python.  Bin k covers
    // integer sizes [bin_lo[k], bin_hi[k]-1].  Consumed directly so the C++
    // never re-derives edges from r_ratio (numpy.floor vs std::floor can
    // diverge over many bins, corrupting the last bins).
    P.sia_bin_lo.resize(P.I_bin);
    P.sia_bin_hi.resize(P.I_bin);
    for (int k = 0; k < P.I_bin; ++k) {
        P.sia_bin_lo[k] = static_cast<int>(
            require_param(p, "sia_bin_lo_" + std::to_string(k)));
        P.sia_bin_hi[k] = static_cast<int>(
            require_param(p, "sia_bin_hi_" + std::to_string(k)));
    }
    P.vac_bin_lo.resize(P.V_bin);
    P.vac_bin_hi.resize(P.V_bin);
    for (int k = 0; k < P.V_bin; ++k) {
        P.vac_bin_lo[k] = static_cast<int>(
            require_param(p, "vac_bin_lo_" + std::to_string(k)));
        P.vac_bin_hi[k] = static_cast<int>(
            require_param(p, "vac_bin_hi_" + std::to_string(k)));
    }

    // If bin_moment mode (physics_option 2 or 3), override N_eq
    // Layout: [discrete SIA | binned SIA moments | discrete VAC |
    //          binned VAC moments | He state]
    // Always override when physics_option >= 2, even if I_bin=V_bin=0
    // (pure-discrete subset of the full domain).
    if (P.physics_option >= 2) {
        int n_sia = P.i_discrete + P.n_mom * P.I_bin;
        int n_vac = P.v_discrete + P.n_mom * P.V_bin;
        int n_he;
        if (P.he_mode == 1)
            n_he = qss ? P.V_bin : P.V_bin + 1;  // Q_k per vac bin
        else
            n_he = qss ? 1       : 2;             // Q_tot scalar
        P.N_eq = n_sia + n_vac + n_he + 5;  // +5 for conservation accounting
        P.y0.resize(P.N_eq);
    }

    // Initial conditions
    for (int k = 0; k < P.N_eq; ++k)
        P.y0[k] = optional_param(p, "y0_" + std::to_string(k), 1e-100);

    P.C_floor = optional_param(p, "C_floor", 1e-15);

    // Solver settings
    P.t_begin  = require_param(p, "t_begin");
    P.t_end    = require_param(p, "t_end");
    P.n_points = static_cast<int>(require_param(p, "n_points"));
    P.log_time = (optional_param(p, "log_time", 1.0) > 0.5);
    P.rtol     = optional_param(p, "rtol",  1e-8);
    P.atol     = optional_param(p, "atol",  1e-20);

    // Integration method
    P.linsol    = static_cast<int>(optional_param(p, "linsol",    0.0));
    P.mu        = static_cast<int>(optional_param(p, "mu",
                                   static_cast<double>(P.N_eq - 1)));
    P.ml        = static_cast<int>(optional_param(p, "ml",
                                   static_cast<double>(P.N_eq - 1)));
    P.max_order = static_cast<int>(optional_param(p, "max_order", 4.0));
    P.hmin      = optional_param(p, "hmin", 0.0);
    P.hmax      = optional_param(p, "hmax", 0.0);

    // Window parameters
    P.window_mode          = static_cast<int>(optional_param(p, "window_mode",       0.0));
    // Single initial width shared by SIA and VAC; defaults to the larger domain.
    P.window_width         = static_cast<int>(optional_param(p, "window_width",
                                 static_cast<double>(std::max(P.I, P.V))));
    P.concentration_threshold = optional_param(p, "concentration_threshold", 1e-18);
    P.window_pad           = static_cast<int>(optional_param(p, "window_pad", 10.0));
    // VAC expansion pad defaults to the SIA pad if not supplied.
    P.window_pad_v         = static_cast<int>(optional_param(p, "window_pad_v",
                                 static_cast<double>(P.window_pad)));
    P.window_expand_factor = optional_param(p, "window_expand_factor", 0.0);
    P.window_check_every   = static_cast<int>(optional_param(p, "window_check_every", 1.0));
    P.window_C_contract    = optional_param(p, "window_C_contract",    0.0);
    P.window_min_active_i  = static_cast<int>(optional_param(p, "window_min_active_i", 5.0));
    P.window_nuc_guard     = optional_param(p, "window_nuc_guard",     0.0);
    P.window_t_start       = optional_param(p, "window_t_start",  10.0);
    P.Ni_extend_tol        = optional_param(p, "Ni_extend_tol",    0.0);
    P.Ni_extend_margin     = static_cast<int>(optional_param(p, "Ni_extend_margin", 0.0));

    // Preconditioner auto-selection: Jacobi by default, Woodbury when the
    // Jacobian has a true bordered-arrow structure from coalescence reactions.
    //
    //   Bordered-arrow ⇔ multiple mobile cluster sizes coalesce
    //   (i_mobile ≥ 2 or v_mobile ≥ 2).  With only monomer mobility
    //   (i_mobile=1, v_mobile=1) the dense border collapses to two columns
    //   and Jacobi+GMRES converges adequately on its own.
    //
    // Additional gating:
    //   - linsol == 2 (GMRES) — preconditioners are only consulted by the
    //     iterative linear solver.
    //   - window_mode == 0 (full domain) — sliding-window modes keep the
    //     active system at 50–200 unknowns, where Woodbury's ~58-RHS setup
    //     cost outweighs its convergence advantage.
    //
    // The user can always override via the "prec_type" key in the param file.
    {
        const bool gmres          = (P.linsol == 2);
        const bool full_domain    = (P.window_mode == 0);
        const bool has_coalescence = (P.i_mobile >= 2) || (P.v_mobile >= 2);
        const bool use_woodbury   = gmres && full_domain && has_coalescence;
        P.prec_type = static_cast<int>(optional_param(p, "prec_type",
                          use_woodbury ? 1.0 : 0.0));
    }
    // Bandwidth auto-computed from mobility cutoffs
    P.prec_bw   = static_cast<int>(optional_param(p, "prec_bw",
                      static_cast<double>(std::max(2 * P.i_mobile, 2 * P.v_mobile) + 1)));
    // Rank = number of mobile species (dense columns in Jacobian)
    P.prec_rank = static_cast<int>(optional_param(p, "prec_rank",
                      static_cast<double>(P.i_mobile + P.v_mobile)));
    P.prec_gamma = 0.0;

    // Build mobile species index list for Woodbury preconditioner
    if (P.prec_type == 1) {
        P.prec_mobile_idx.resize(P.prec_rank);
        // First i_mobile entries: SIA mobile indices 0..i_mobile-1
        for (int k = 0; k < P.i_mobile; ++k)
            P.prec_mobile_idx[k] = k;
        // Next v_mobile entries: VAC mobile indices I..I+v_mobile-1
        // (or i_discrete offset for bin_moment modes)
        const int vac_off = (P.physics_option >= 2)
                            ? (P.i_discrete + P.n_mom * P.I_bin)
                            : P.I;
        for (int k = 0; k < P.v_mobile; ++k)
            P.prec_mobile_idx[P.i_mobile + k] = vac_off + k;
    }

    P.verbose = (optional_param(p, "verbose", 0.0) > 0.5);

    return P;
}
