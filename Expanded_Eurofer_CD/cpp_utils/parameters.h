/**
 * parameters.h — Expanded_Eurofer_CD solver parameter struct.
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
 *   y[0..N-1]       — SIA clusters c_n, n=1..N
 *   y[N..N+M-1]     — vacancy/bubble c_m (marginal), m=1..M
 *   y[N+M]          — Q_tot (total He in voids)
 *   y[N+M+1]        — c_h (free He)  [omitted when he_options=1 (QSS)]
 *   N_eq = N + M + 2  (dynamic)  or  N + M + 1  (quasi_steady_state)
 *
 * full_CD_fusion  (Case 1, he_mode=1):
 *   y[0..N-1]       — SIA clusters c_n
 *   y[N..N+M-1]     — c_m^tot
 *   y[N+M..N+2M-1]  — Q_m (He per class)
 *   y[N+2M]         — c_h  [omitted when he_options=1 (QSS)]
 *   N_eq = N + 2M + 1  (dynamic)  or  N + 2M  (quasi_steady_state)
 *
 * bin_moment_CD_fission/fusion (physics_option 2/3):
 *   y[0..2K-1]      — SIA bin moments [μ_0^(0), μ_0^(1), ..., μ_{K-1}^(1)]
 *   y[2K..2K+M-1]   — c_m
 *   y[2K+M..]       — He variables (same as Case 2 or Case 1)
 *   N_eq = 2K + M + 2  or  2K + 2M + 1  (dynamic)
 *   N_eq = 2K + M + 1  or  2K + 2M      (quasi_steady_state)
 *
 * he_options:  0 = dynamic (c_h is an ODE state, Eq. 157)
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
    int N, M, N_eq;
    int Ni_max;   // = N (kept for legacy compatibility)

    // ── Physics and solver mode ─────────────────────────────────────────────
    // physics_option: 0=full_CD_fission, 1=full_CD_fusion,
    //                  2=bin_moment_CD_fission, 3=bin_moment_CD_fusion
    // he_mode:         0=Case 2 (decoupled/fission), 1=Case 1 (mean-field/fusion)
    // he_options:      0=dynamic (ODE), 1=quasi_steady_state (algebraic c_h)
    int physics_option;
    int he_mode;
    int he_options;   // 0=dynamic, 1=quasi_steady_state

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
    int n_max_i;   // max mobile SIA cluster size (1D glide cutoff)
    int m_max_v;   // max mobile vacancy cluster size

    // ── Vacancy cluster arrays [M] ────────────────────────────────────────────
    std::vector<double> KVV;      // Cv capture by void m  (K_VAC_grow)
    std::vector<double> KVI;      // Ci annihilation at void m  (K_VAC_shrink)
    std::vector<double> GVV;      // thermal vacancy emission from void m
    std::vector<double> KHeV;     // He capture by void m
    std::vector<double> Pr_VAC;   // cascade vacancy production rate [M]
    std::vector<double> m13;      // m^{1/3} factors [M]

    // ── SIA cluster arrays [N] ────────────────────────────────────────────────
    std::vector<double> KII;      // Ci capture by SIA loop n  (K_SIA_grow)
    std::vector<double> KIV;      // Cv capture by SIA loop n  (K_SIA_shrink)
    std::vector<double> GII;      // thermal SIA emission from loop n
    std::vector<double> k2_SIA;   // total sink rate for SIA cluster n
    std::vector<double> Pr_SIA;   // cascade SIA production rate [N]

    // ── 1D glide prefactors (Eq. 141) ─────────────────────────────────────────
    // K_1D_eff(n, m) = K_1D_pref[n-1] · m^{1/3} / (1 + B_rot·L̂²·m^{-1/3})
    std::vector<double> K_1D_pref;   // [N]: A_sph · ω_n^{1D} / Ω

    // Legacy separable cross-term (for backward compatibility with solver.cpp)
    std::vector<double> K_IclV_ns;   // [N]: 4π·r0·Di·n^{-2/3} / Ω
    std::vector<double> K_IclV_ni;   // [N]: 4π·r0·Di / (n·Ω)

    // ── Scalar physics ────────────────────────────────────────────────────────
    double G_He;          // He transmutation rate [at.frac/s]
    double k2_disl_v;     // vacancy fixed sink [s^-1]
    double k2_disl_i;     // SIA fixed sink [s^-1]
    double k2_disl_He;    // He fixed sink [s^-1]
    double Cv_eq;         // thermal vacancy equilibrium concentration
    double K_iv;          // V–SIA recombination prefactor

    // ── He parameters ─────────────────────────────────────────────────────────
    double beta_He;       // He emission: ν_h·exp(−(E_b_hV+E_m_h)/kBT) [s^-1]
    double delta_He;      // He pressure coefficient [eV]
    double beta_He_exp;   // He pressure power-law exponent
    double kBT;           // k_B·T [eV]
    double L_He_max;      // max He loading per vacancy cluster

    // ── Bin-moment parameters (Chapter 9) ────────────────────────────────────
    int    K_bins;        // number of logarithmic bins
    double r_ratio;       // bin ratio r > 1 (Eq. 188)
    int    n1_bin;        // minimum cluster size in binned region

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
    int backend;    // 0=CVODE, 1=ARKODE ARKStep DIRK
    int lmm;        // CVODE: 2=CV_BDF (default), 1=CV_ADAMS
    int linsol;     // 0=dense, 1=band, 2=gmres
    int mu, ml;     // band solver bandwidths
    int max_order;  // 0 = solver default
    int ark_table;  // ARKODE_DIRKTableID (default 111)

    // ── Dynamic window (cpp_sliding_win / sliding_OpenMP) ─────────────────────
    // window_mode: 0=full, 3=Phase III (constant width), 4=Phase IV (OpenMP)
    int    window_mode;
    int    window_w0_v;
    int    window_w0_i;
    double window_C_expand;
    int    window_expand_pad;
    double window_expand_factor;
    int    window_check_every;
    double window_C_contract;
    int    window_min_active_i;
    int    window_prec;
    double window_nuc_guard;
    int    window_width;
    double window_t_start;
    int    window_N_thresh;
    int    window_omp_threads;
    int    window_gmres_maxl;
    double Ni_extend_tol;
    int    Ni_extend_margin;

    // ── Jacobi preconditioner storage ─────────────────────────────────────────
    std::vector<double> prec_diag;
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

    // Cluster sizes — support both old (Ni/Nv) and new (N/M) key names
    if (p.count("N"))  P.N = static_cast<int>(require_param(p, "N"));
    else               P.N = static_cast<int>(require_param(p, "Ni"));
    if (p.count("M"))  P.M = static_cast<int>(require_param(p, "M"));
    else               P.M = static_cast<int>(require_param(p, "Nv"));

    P.Ni_max = static_cast<int>(optional_param(p, "Ni_max", static_cast<double>(P.N)));

    // Physics mode
    P.physics_option = static_cast<int>(optional_param(p, "physics_option_int", 0));
    P.he_mode        = static_cast<int>(optional_param(p, "he_mode",             0));
    P.he_options     = static_cast<int>(optional_param(p, "qss_He",              0));
    // he_options: 0=dynamic (c_h is ODE state), 1=quasi_steady_state (c_h algebraic)

    // State vector size depends on he_mode and he_options:
    //   Case 2 dynamic:           N + M + 2
    //   Case 2 quasi_steady_state: N + M + 1  (c_h removed from state)
    //   Case 1 dynamic:           N + 2M + 1
    //   Case 1 quasi_steady_state: N + 2M     (c_h removed from state)
    const bool qss = (P.he_options == 1);
    int n_he_extra;
    if (P.he_mode == 1)
        n_he_extra = qss ? P.M     : P.M + 1;
    else
        n_he_extra = qss ? 1       : 2;
    P.N_eq = P.N + P.M + n_he_extra;

    // Geometric prefactors
    P.A_sph  = optional_param(p, "A_sph",  7.818);
    P.A_loop = optional_param(p, "A_loop", 10.78);
    P.B_rot  = optional_param(p, "B_rot",  2.627);
    P.L_hat  = optional_param(p, "L_hat",  50.0);

    // Solute trapping
    P.trap_SIA  = optional_param(p, "trap_SIA",  0.0);
    P.trap_VAC  = optional_param(p, "trap_VAC",  0.0);
    P.trap_loop = optional_param(p, "trap_loop", 0.0);

    // Mobility cutoffs
    P.n_max_i = static_cast<int>(optional_param(p, "n_max_i", 100.0));
    P.m_max_v = static_cast<int>(optional_param(p, "m_max_v", 5.0));

    // Resize arrays
    P.KVV.resize(P.M);   P.KVI.resize(P.M);   P.GVV.resize(P.M);
    P.KHeV.resize(P.M);  P.Pr_VAC.resize(P.M); P.m13.resize(P.M);
    P.KII.resize(P.N);   P.KIV.resize(P.N);    P.GII.resize(P.N);
    P.k2_SIA.resize(P.N); P.Pr_SIA.resize(P.N);
    P.K_1D_pref.resize(P.N);
    P.K_IclV_ns.resize(P.N); P.K_IclV_ni.resize(P.N);
    P.y0.resize(P.N_eq);

    // Vacancy arrays
    for (int k = 0; k < P.M; ++k) {
        P.KVV[k]    = require_param(p, "KVV_"    + std::to_string(k));
        P.KVI[k]    = require_param(p, "KVI_"    + std::to_string(k));
        P.GVV[k]    = require_param(p, "GVV_"    + std::to_string(k));
        P.KHeV[k]   = require_param(p, "KHeV_"   + std::to_string(k));
        P.Pr_VAC[k] = require_param(p, "Pr_VAC_" + std::to_string(k));
        P.m13[k]    = require_param(p, "m13_"    + std::to_string(k));
    }

    // SIA arrays
    for (int k = 0; k < P.N; ++k) {
        P.KII[k]         = require_param(p, "KII_"        + std::to_string(k));
        P.KIV[k]         = require_param(p, "KIV_"        + std::to_string(k));
        P.GII[k]         = require_param(p, "GII_"        + std::to_string(k));
        P.k2_SIA[k]      = require_param(p, "k2_SIA_"     + std::to_string(k));
        P.Pr_SIA[k]      = require_param(p, "Pr_SIA_"     + std::to_string(k));
        P.K_1D_pref[k]   = optional_param(p, "K_1D_pref_" + std::to_string(k), 0.0);
        P.K_IclV_ns[k]   = optional_param(p, "K_IclV_ns_" + std::to_string(k), 0.0);
        P.K_IclV_ni[k]   = optional_param(p, "K_IclV_ni_" + std::to_string(k), 0.0);
    }

    // Scalar physics
    P.G_He       = require_param(p, "G_He");
    P.k2_disl_v  = require_param(p, "k2_disl_v");
    P.k2_disl_i  = require_param(p, "k2_disl_i");
    P.k2_disl_He = require_param(p, "k2_disl_He");
    P.Cv_eq      = require_param(p, "Cv_eq");
    P.K_iv       = optional_param(p, "K_iv", 0.0);
    P.beta_He    = require_param(p, "beta_He");
    P.delta_He   = require_param(p, "delta_He");
    P.beta_He_exp= require_param(p, "beta_He_exp");
    P.kBT        = require_param(p, "kBT");
    P.L_He_max   = require_param(p, "L_He_max");

    // Bin-moment parameters
    P.K_bins  = static_cast<int>(optional_param(p, "K_bins", 0.0));
    P.r_ratio = optional_param(p, "r_ratio", 2.0);
    P.n1_bin  = static_cast<int>(optional_param(p, "n1_bin", 1.0));

    // If bin_moment mode, override N_eq
    if (P.K_bins > 0) {
        int n_he_extra2;
        if (P.he_mode == 1)
            n_he_extra2 = qss ? P.M     : P.M + 1;
        else
            n_he_extra2 = qss ? 1       : 2;
        P.N_eq = 2 * P.K_bins + P.M + n_he_extra2;
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
    P.backend   = static_cast<int>(optional_param(p, "backend",   0.0));
    P.lmm       = static_cast<int>(optional_param(p, "lmm",       2.0));
    P.linsol    = static_cast<int>(optional_param(p, "linsol",    0.0));
    P.mu        = static_cast<int>(optional_param(p, "mu",
                                   static_cast<double>(P.N_eq - 1)));
    P.ml        = static_cast<int>(optional_param(p, "ml",
                                   static_cast<double>(P.N_eq - 1)));
    P.max_order = static_cast<int>(optional_param(p, "max_order", 4.0));
    P.ark_table = static_cast<int>(optional_param(p, "ark_table", 111.0));

    // Window parameters
    P.window_mode          = static_cast<int>(optional_param(p, "window_mode",       0.0));
    P.window_w0_v          = static_cast<int>(optional_param(p, "window_w0_v",
                                 static_cast<double>(P.M)));
    P.window_w0_i          = static_cast<int>(optional_param(p, "window_w0_i",
                                 static_cast<double>(P.N)));
    P.window_C_expand      = optional_param(p, "window_C_expand",      1e-18);
    P.window_expand_pad    = static_cast<int>(optional_param(p, "window_expand_pad", 10.0));
    P.window_expand_factor = optional_param(p, "window_expand_factor", 0.0);
    P.window_check_every   = static_cast<int>(optional_param(p, "window_check_every", 1.0));
    P.window_C_contract    = optional_param(p, "window_C_contract",    0.0);
    P.window_min_active_i  = static_cast<int>(optional_param(p, "window_min_active_i", 5.0));
    P.window_prec          = static_cast<int>(optional_param(p, "window_prec",          0.0));
    P.window_nuc_guard     = optional_param(p, "window_nuc_guard",     0.0);
    P.window_width         = static_cast<int>(optional_param(p, "window_width",   500.0));
    P.window_t_start       = optional_param(p, "window_t_start",  10.0);
    P.window_N_thresh      = static_cast<int>(optional_param(p, "window_N_thresh",1000.0));
    P.window_omp_threads   = static_cast<int>(optional_param(p, "window_omp_threads", 0.0));
    P.window_gmres_maxl    = static_cast<int>(optional_param(p, "window_gmres_maxl",  20.0));
    P.Ni_extend_tol        = optional_param(p, "Ni_extend_tol",    0.0);
    P.Ni_extend_margin     = static_cast<int>(optional_param(p, "Ni_extend_margin", 0.0));

    return P;
}
