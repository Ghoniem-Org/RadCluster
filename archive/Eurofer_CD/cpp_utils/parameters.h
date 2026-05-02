/**
 * parameters.h – Eurofer_CD solver parameter struct.
 *
 * All quantities that the ODE right-hand side needs are pre-computed on the
 * Python side (via InputData + ReactionRates) and forwarded via a parameter
 * file (--param_file=<path>).  build_parameters() unpacks those into this
 * struct so the RHS performs only arithmetic (no exp/log except for the
 * He-pressure GVV_eff correction, which is O(Nv) per call).
 *
 * Mirrors: py_utils/input_data.py + py_utils/reaction_rates.py
 *
 * State vector y[N_EQ], N_EQ = Ni + Nv + 1  (runtime values):
 *   y[0 .. Ni-1]       – Ci1 .. Ci_Ni   (SIA clusters)
 *   y[Ni .. Ni+Nv-1]   – Cv1 .. Cv_Nv   (vacancy clusters)
 *   y[Ni+Nv]           – C_He            (free He)
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
    // ── Cluster size limits (runtime) ─────────────────────────────────────────
    int Nv, Ni, N_EQ;

    // ── Pre-computed rate constant arrays ─────────────────────────────────────
    // Vacancy cluster arrays (0-indexed; index k → cluster size k+1)
    std::vector<double> KVV;           // Cv capture by vacancy cluster [Nv]
    std::vector<double> KVI;           // Ci capture (annihilation) by vacancy cluster [Nv]
    std::vector<double> GVV;           // thermal vacancy emission from vacancy cluster [Nv]
    std::vector<double> KHeV;          // He capture by vacancy cluster [Nv]
    std::vector<double> Pr_VAC;        // cascade vacancy production rate [Nv]

    // Interstitial cluster arrays (0-indexed; index k → cluster size k+1)
    std::vector<double> KII;           // Ci capture by SIA cluster [Ni]
    std::vector<double> KIV;           // Cv capture (shrink) by SIA cluster [Ni]
    std::vector<double> GII;           // thermal SIA emission from SIA cluster [Ni]
    std::vector<double> k2_SIA;        // dislocation sink rate for SIA cluster n [Ni]
    std::vector<double> Pr_SIA;        // cascade SIA production rate [Ni]

    // ── K_IclV separable cross-term coefficients ──────────────────────────────
    // K_IclV[n-1,m-1] = K_IclV_ns[n-1] + K_IclV_ni[n-1] * m13[m-1]
    // where:
    //   K_IclV_ns[n-1] = 4π·r0·Di·n^{-2/3} / Ω   (n-scale factor)
    //   K_IclV_ni[n-1] = 4π·r0·Di / (n·Ω)          (n-inverse factor)
    //   m13[m-1]        = m^{1/3}                   (m-radius factor)
    // Index k=0 (n=1) is 0 (mono-SIA excluded from cluster cross-recombination).
    std::vector<double> K_IclV_ns;     // [Ni]: n-scale factor
    std::vector<double> K_IclV_ni;     // [Ni]: n-inverse factor
    std::vector<double> m13;           // [Nv]: m^{1/3}

    // ── Scalar physics ────────────────────────────────────────────────────────
    double G_He;          // He transmutation production rate [atom frac / s]
    double k2_disl_v;     // dislocation sink for mono-vacancy [s^-1]
    double k2_disl_i;     // dislocation sink for mono-SIA [s^-1]
    double k2_disl_He;    // dislocation sink for free He [s^-1]
    double Cv_eq;         // thermal equilibrium vacancy concentration

    // ── He-pressure GVV_eff parameters ────────────────────────────────────────
    // GVV_eff[m] = GVV[m] * exp(-ell_m * dE / kBT)
    // where ell_m = min(KHeV[m]*C_He/beta_He, L_He_max)
    //       dE    = delta_He * beta_He_exp / m * (ell_c/m)^{beta_He_exp - 1}
    double beta_He;       // He emission rate = nu_He * exp(-(E_b_HeV+E_m_He)/kBT) [s^-1]
    double delta_He;      // He pressure coefficient [eV]
    double beta_He_exp;   // He pressure exponent (power-law)
    double kBT;           // k_B * T [eV]
    double L_He_max;      // max He loading per vacancy cluster (clamping)

    // ── Initial conditions ─────────────────────────────────────────────────────
    std::vector<double> y0;   // [N_EQ]

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
    int backend;    // 0=CVODE (default), 1=ARKODE ARKStep DIRK
    int lmm;        // CVODE: 2=CV_BDF (stiff, default), 1=CV_ADAMS
    int linsol;     // 0=dense (default), 1=band, 2=gmres
    int mu;         // upper bandwidth for band solver
    int ml;         // lower bandwidth for band solver
    int max_order;  // 0 = solver default
    int ark_table;  // ARKODE_DIRKTableID integer (default 111)

    // ── Dynamic window solver ──────────────────────────────────────────────────
    // The window operates on the SIA cluster equations only.
    // Vacancy clusters and free He are always fully active.
    //
    // window_mode=0: full system (default)
    // window_mode=1: Phase I   — upper truncation on SIA: active [Ci1..Ci_{x_hi_i}]
    // window_mode=2: Phase II  — sliding window on SIA: [Ci1] + [Ci_{x_lo_i}..Ci_{x_hi_i}]
    //                            x_hi expands when top SIA cluster exceeds window_C_expand
    //                            x_lo_i advances when lowest cluster reaches QSS
    // window_mode=3: Phase III — constant-width window (window_width) sliding upward
    // window_mode=4: Phase IV  — same as Phase III + OpenMP intra-RHS parallelism
    int    window_mode;
    int    window_w0_v;           // initial vacancy window size  (unused; always Nv)
    int    window_w0_i;           // initial SIA window size (default: Ni)
    double window_C_expand;       // expand upper bound when C[x_hi] > this
    int    window_expand_pad;     // minimum additive increment per expansion step
    double window_expand_factor;  // geometric expansion factor (0 = additive only)
    int    window_check_every;    // check expansion every N output points
    // Phase II only:
    double window_C_contract;     // contract lower bound when |dC/dt|/C < this (0=off)
    int    window_min_active_i;   // minimum active SIA window size
    int    window_prec;           // 0=no preconditioner, 1=Jacobi diagonal
    double window_nuc_guard;      // nucleation guard threshold (0=disabled)
    // Phase III:
    int    window_width;          // constant window width (default 500)
    double window_t_start;        // suppress lower sliding until t > this
    int    window_N_thresh;       // activate Phase III only if Ni > this
    // Phase IV:
    int window_omp_threads;       // 0 = OMP_NUM_THREADS, >0 = explicit
    int window_gmres_maxl;        // GMRES Krylov subspace size (0 = SUNDIALS default of 5)

    // ── Dynamic Ni extension ───────────────────────────────────────────────────
    int    Ni_max;
    double Ni_extend_tol;
    int    Ni_extend_margin;

    // ── Jacobi preconditioner storage ─────────────────────────────────────────
    std::vector<double> prec_diag;
};

// ── CLI / file argument helpers ───────────────────────────────────────────────

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
 * Format: one "key=value" per line; blank lines and '#' lines are ignored.
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
        if (pos == std::string::npos) {
            std::cerr << "Parameter file line " << lineno
                      << ": missing '=' in \"" << line << "\"\n";
            continue;
        }
        std::string key = line.substr(0, pos);
        double val = 0.0;
        try {
            val = std::stod(line.substr(pos + 1));
        } catch (...) {
            std::cerr << "Parameter file line " << lineno
                      << ": invalid value for \"" << key << "\"\n";
            exit(1);
        }
        props[key] = val;
    }
    return props;
}

// Unpack a parsed parameter map into a Parameters struct.
inline Parameters build_parameters(const std::map<std::string, double>& p) {
    Parameters P{};

    // ── Cluster size limits ───────────────────────────────────────────────────
    P.Nv     = static_cast<int>(require_param(p, "Nv"));
    P.Ni     = static_cast<int>(require_param(p, "Ni"));
    P.Ni_max = static_cast<int>(optional_param(p, "Ni_max", static_cast<double>(P.Ni)));
    P.N_EQ   = P.Ni_max + P.Nv + 1;  // SIA + vacancy + He

    // ── Resize arrays ─────────────────────────────────────────────────────────
    P.KVV.resize(P.Nv);   P.KVI.resize(P.Nv);   P.GVV.resize(P.Nv);
    P.KHeV.resize(P.Nv);  P.Pr_VAC.resize(P.Nv); P.m13.resize(P.Nv);

    P.KII.resize(P.Ni_max);   P.KIV.resize(P.Ni_max);   P.GII.resize(P.Ni_max);
    P.k2_SIA.resize(P.Ni_max); P.Pr_SIA.resize(P.Ni_max);
    P.K_IclV_ns.resize(P.Ni_max); P.K_IclV_ni.resize(P.Ni_max);

    P.y0.resize(P.N_EQ);

    // ── Vacancy cluster arrays ────────────────────────────────────────────────
    for (int k = 0; k < P.Nv; ++k) {
        P.KVV[k]    = require_param(p, "KVV_"    + std::to_string(k));
        P.KVI[k]    = require_param(p, "KVI_"    + std::to_string(k));
        P.GVV[k]    = require_param(p, "GVV_"    + std::to_string(k));
        P.KHeV[k]   = require_param(p, "KHeV_"   + std::to_string(k));
        P.Pr_VAC[k] = require_param(p, "Pr_VAC_" + std::to_string(k));
        P.m13[k]    = require_param(p, "m13_"    + std::to_string(k));
    }

    // ── SIA cluster arrays ────────────────────────────────────────────────────
    for (int k = 0; k < P.Ni_max; ++k) {
        P.KII[k]        = require_param(p, "KII_"        + std::to_string(k));
        P.KIV[k]        = require_param(p, "KIV_"        + std::to_string(k));
        P.GII[k]        = require_param(p, "GII_"        + std::to_string(k));
        P.k2_SIA[k]     = require_param(p, "k2_SIA_"     + std::to_string(k));
        P.Pr_SIA[k]     = require_param(p, "Pr_SIA_"     + std::to_string(k));
        P.K_IclV_ns[k]  = require_param(p, "K_IclV_ns_"  + std::to_string(k));
        P.K_IclV_ni[k]  = require_param(p, "K_IclV_ni_"  + std::to_string(k));
    }

    // ── Scalar physics ────────────────────────────────────────────────────────
    P.G_He       = require_param(p, "G_He");
    P.k2_disl_v  = require_param(p, "k2_disl_v");
    P.k2_disl_i  = require_param(p, "k2_disl_i");
    P.k2_disl_He = require_param(p, "k2_disl_He");
    P.Cv_eq      = require_param(p, "Cv_eq");
    P.beta_He    = require_param(p, "beta_He");
    P.delta_He   = require_param(p, "delta_He");
    P.beta_He_exp= require_param(p, "beta_He_exp");
    P.kBT        = require_param(p, "kBT");
    P.L_He_max   = require_param(p, "L_He_max");

    // ── Initial conditions ─────────────────────────────────────────────────────
    for (int k = 0; k < P.N_EQ; ++k)
        P.y0[k] = require_param(p, "y0_" + std::to_string(k));

    // ── Floor ─────────────────────────────────────────────────────────────────
    P.C_floor = optional_param(p, "C_floor", 1e-100);

    // ── Solver settings ────────────────────────────────────────────────────────
    P.t_begin  = require_param(p, "t_begin");
    P.t_end    = require_param(p, "t_end");
    P.n_points = static_cast<int>(require_param(p, "n_points"));
    P.log_time = (optional_param(p, "log_time", 1.0) > 0.5);
    P.rtol     = optional_param(p, "rtol",   1e-8);
    P.atol     = optional_param(p, "atol",   1e-50);

    // ── Integration method ─────────────────────────────────────────────────────
    P.backend   = static_cast<int>(optional_param(p, "backend",   0.0));
    P.lmm       = static_cast<int>(optional_param(p, "lmm",       2.0));
    P.linsol    = static_cast<int>(optional_param(p, "linsol",    0.0));
    P.mu        = static_cast<int>(optional_param(p, "mu",  static_cast<double>(P.N_EQ - 1)));
    P.ml        = static_cast<int>(optional_param(p, "ml",  static_cast<double>(P.N_EQ - 1)));
    P.max_order = static_cast<int>(optional_param(p, "max_order", 0.0));
    P.ark_table = static_cast<int>(optional_param(p, "ark_table", 111.0));

    // ── Dynamic window parameters ──────────────────────────────────────────────
    P.window_mode          = static_cast<int>(optional_param(p, "window_mode", 0.0));
    P.window_w0_v          = static_cast<int>(optional_param(p, "window_w0_v",
                                 static_cast<double>(P.Nv)));
    P.window_w0_i          = static_cast<int>(optional_param(p, "window_w0_i",
                                 static_cast<double>(P.Ni)));
    P.window_C_expand      = optional_param(p, "window_C_expand",      1e-18);
    P.window_expand_pad    = static_cast<int>(optional_param(p, "window_expand_pad",    10.0));
    P.window_expand_factor = optional_param(p, "window_expand_factor", 0.0);
    P.window_check_every   = static_cast<int>(optional_param(p, "window_check_every",    1.0));
    P.window_C_contract    = optional_param(p, "window_C_contract",    0.0);
    P.window_min_active_i  = static_cast<int>(optional_param(p, "window_min_active_i",  5.0));
    P.window_prec          = static_cast<int>(optional_param(p, "window_prec",           0.0));
    P.window_nuc_guard     = optional_param(p, "window_nuc_guard", 0.0);
    P.window_width         = static_cast<int>(optional_param(p, "window_width",    500.0));
    P.window_t_start       = optional_param(p, "window_t_start",  10.0);
    P.window_N_thresh      = static_cast<int>(optional_param(p, "window_N_thresh", 1000.0));
    P.window_omp_threads   = static_cast<int>(optional_param(p, "window_omp_threads", 0.0));
    P.window_gmres_maxl    = static_cast<int>(optional_param(p, "window_gmres_maxl",  20.0));

    P.Ni_extend_tol    = optional_param(p, "Ni_extend_tol",    0.0);
    P.Ni_extend_margin = static_cast<int>(optional_param(p, "Ni_extend_margin", 0.0));

    return P;
}
