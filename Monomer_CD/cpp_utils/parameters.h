/**
 * parameters.h – ClusterDynamics solver parameter struct.
 *
 * All quantities that the ODE right-hand side needs are pre-computed on the
 * Python side (via InputData + ReactionRates) and forwarded via a parameter
 * file (--param_file=<path>).  build_parameters() unpacks those into this
 * struct so the RHS performs only arithmetic.
 *
 * Mirrors: py_utils/input_data.py + py_utils/reaction_rates.py
 *
 * State vector y[N_EQ], N_EQ = Nv + Ni  (runtime values):
 *   y[0 .. Nv-1]       – Cv1 .. Cv_Nv   (vacancy clusters)
 *   y[Nv .. Nv+Ni-1]   – Ci1 .. Ci_Ni   (interstitial clusters)
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

    // ── Pre-computed rate constant arrays (0-indexed, size k → cluster size k+1) ──
    std::vector<double> KCV;   // spherical capture of vacancies    by vacancy clusters  [Nv]
    std::vector<double> KCI;   // spherical capture of interstitials by vacancy clusters [Nv]
    std::vector<double> KLV;   // circular-loop capture of vacancies      by i-clusters  [Ni]
    std::vector<double> KLI;   // circular-loop capture of interstitials   by i-clusters  [Ni]
    std::vector<double> GCV;   // thermal emission (vacancy) from vacancy    clusters     [Nv]
    std::vector<double> GLV;   // thermal emission (vacancy) from interstitial clusters   [Ni]

    // ── Scalar physics ────────────────────────────────────────────────────────
    double P_prod;    // production rate [dpa s^-1]
    double alpha;     // recombination coefficient [s^-1]
    double Cv_eq;     // thermal equilibrium vacancy concentration
    double C2v_eq;    // thermal equilibrium divacancy concentration
    double Z_v;       // vacancy–dislocation bias factor
    double Z_i;       // interstitial–dislocation bias factor
    double rho_d;     // dislocation density [cm cm^-3]
    double Dv;        // vacancy diffusion coefficient [cm^2 s^-1]
    double D2v;       // divacancy diffusion coefficient [cm^2 s^-1]
    double Di;        // interstitial diffusion coefficient [cm^2 s^-1]
    double K_nuc_i;   // interstitial nucleation rate constant

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
    // window_mode=0: full system (default)
    // window_mode=1: Phase I   — upper truncation only [1..x_hi_v] x [1..x_hi_i]
    // window_mode=2: Phase II  — sliding window: upper + lower truncation
    //                            Active set: [Cv1..Cv_{x_hi_v}] + [Ci1] + [Ci_{x_lo_i}..Ci_{x_hi_i}]
    //                            x_hi expands when top cluster exceeds window_C_expand.
    //                            x_lo_i advances when lowest cluster reaches QSS.
    // window_mode=3: Phase III — constant-width window (window_width) sliding upward.
    // window_mode=4: Phase IV  — same as Phase III + OpenMP intra-RHS parallelism
    //                            + pre-allocated scratch buffers (no malloc per RHS call).
    // All window modes use GMRES (matrix-free). Phase II-IV support Jacobi preconditioner.
    // The full Nv+Ni output row is always written; truncated species are zero.
    int    window_mode;           // 0=off, 1=Phase I, 2=Phase II
    int    window_w0_v;           // initial vacancy  window size (default: Nv)
    int    window_w0_i;           // initial interstitial window size (default: Ni)
    double window_C_expand;       // expand upper bound when C[x_hi] > this
    int    window_expand_pad;     // minimum additive increment per expansion step
    double window_expand_factor;  // geometric expansion: new_hi = max(hi*factor, hi+pad)
                                  //   0.0 = additive only (default, backward-compatible)
    int    window_check_every;    // check expansion/contraction every N output points
    // Phase II only:
    double window_C_contract;     // contract lower bound when |dC/dt|/C < this (0 = off)
    int    window_min_active_i;   // minimum active interstitial window size (safety guard)
    int    window_prec;           // 0=no preconditioner, 1=Jacobi diagonal
    // Phase II contraction guard: block contracting x_lo_i from 2 when nucleation
    // rate exceeds this fraction of the Ci2 outflow rate.  0.0 = disabled (default).
    double window_nuc_guard;      // 0.0 = trust QSS criterion only

    // ── Phase III constant-width sliding window ────────────────────────────────
    // window_mode=3: x_lo_i = max(2, x_hi_i - window_width + 1) coupled to x_hi_i.
    // Lower sliding is suppressed until t > window_t_start (skip nucleation phase).
    // Automatically falls back to full solver when N_EQ <= window_N_thresh.
    int    window_width;      // constant window width (default 500)
    double window_t_start;    // suppress lower sliding until t > this (default 10.0 s)
    int    window_N_thresh;   // activate Phase III only if N_EQ > this (default 1000)

    // ── Phase IV: Multithread-OpenMP ───────────────────────────────────────────
    // window_mode=4: identical sliding-window algorithm as Phase III but all hot
    // loops inside rhs_window_omp() are parallelised with OpenMP.  Pre-allocated
    // scratch buffers (WindowDataOMP::Cv_buf / Ci_buf) eliminate per-call heap
    // allocation from the RHS hot path.
    // Requires compilation with CD_HAVE_OPENMP=1 (cmake finds libomp).
    int window_omp_threads;   // 0 = honour OMP_NUM_THREADS (default), >0 = explicit

    // ── Dynamic Ni extension ───────────────────────────────────────────────────
    // Ni_max: pre-allocated array size; Ni grows from Ni_initial up to Ni_max.
    // Triggers: proximity (x_hi_i >= Ni - Ni_extend_margin) OR
    //           conservation (KLI[Ni-1]*Ci1*Ci_Ni / P_prod > Ni_extend_tol).
    // Extension step: window_expand_pad (reused).  No CVODE reinit needed.
    int    Ni_max;           // maximum Ni (pre-allocated; default = Ni)
    double Ni_extend_tol;    // conservation trigger threshold (0 = disabled)
    int    Ni_extend_margin; // proximity trigger margin in cluster sizes (0 = disabled)

    // ── Phase I Jacobi preconditioner storage ─────────────────────────────────
    // Populated by prec_setup_win1(); sized to Nv+Ni of the active Phase I window.
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
 * Format: one "key=value" per line; blank lines and lines starting with '#'
 * are ignored.  Values are doubles (scientific notation accepted).
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
    P.Nv  = static_cast<int>(require_param(p, "Nv"));
    P.Ni  = static_cast<int>(require_param(p, "Ni"));   // Ni_initial (active upper bound)
    P.Ni_max = static_cast<int>(optional_param(p, "Ni_max", static_cast<double>(P.Ni)));
    P.N_EQ = P.Nv + P.Ni_max;   // always sized to Ni_max for fixed-width output

    // ── Resize arrays ─────────────────────────────────────────────────────────
    P.KCV.resize(P.Nv); P.KCI.resize(P.Nv); P.GCV.resize(P.Nv);
    P.KLV.resize(P.Ni_max); P.KLI.resize(P.Ni_max); P.GLV.resize(P.Ni_max);
    P.y0.resize(P.N_EQ);

    // ── Rate constant arrays ──────────────────────────────────────────────────
    for (int k = 0; k < P.Nv; ++k) {
        P.KCV[k] = require_param(p, "KCV_" + std::to_string(k));
        P.KCI[k] = require_param(p, "KCI_" + std::to_string(k));
        P.GCV[k] = require_param(p, "GCV_" + std::to_string(k));
    }
    for (int k = 0; k < P.Ni_max; ++k) {
        P.KLV[k] = require_param(p, "KLV_" + std::to_string(k));
        P.KLI[k] = require_param(p, "KLI_" + std::to_string(k));
        P.GLV[k] = require_param(p, "GLV_" + std::to_string(k));
    }

    // ── Scalar physics ────────────────────────────────────────────────────────
    P.P_prod   = require_param(p, "P_prod");
    P.alpha    = require_param(p, "alpha");
    P.Cv_eq    = require_param(p, "Cv_eq");
    P.C2v_eq   = require_param(p, "C2v_eq");
    P.Z_v      = require_param(p, "Z_v");
    P.Z_i      = require_param(p, "Z_i");
    P.rho_d    = require_param(p, "rho_d");
    P.Dv       = require_param(p, "Dv");
    P.D2v      = require_param(p, "D2v");
    P.Di       = require_param(p, "Di");
    P.K_nuc_i  = require_param(p, "K_nuc_i");

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

    // Phase III
    P.window_width    = static_cast<int>(optional_param(p, "window_width",    500.0));
    P.window_t_start  = optional_param(p, "window_t_start",  10.0);
    P.window_N_thresh = static_cast<int>(optional_param(p, "window_N_thresh", 1000.0));

    // Phase IV: Multithread-OpenMP
    P.window_omp_threads = static_cast<int>(optional_param(p, "window_omp_threads", 0.0));

    // Dynamic Ni extension
    P.Ni_extend_tol    = optional_param(p, "Ni_extend_tol",    0.0);
    P.Ni_extend_margin = static_cast<int>(optional_param(p, "Ni_extend_margin", 0.0));

    return P;
}
