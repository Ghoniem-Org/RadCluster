/**
 * solver.cpp — RadCluster_1_0 main C++ ODE solver.
 *
 * Drives the RadCluster_1_0 cluster dynamics system for bcc Fe / EUROFER97.
 * Invoked by py_utils/cpp_bridge.py via a parameter file:
 *   solver.exe --param_file=<path>
 *
 * Physics reference: Ghoniem (2026), Rate_Equations.pdf, Eqs. 152, 155, 157,
 *   Chapter 9 (bin moments), Section 8 (He-reduction).
 *
 * Solver modes (window_mode parameter):
 *   0 = full_system     — full system, CVODE BDF.
 *   4 = active_window   — two independent sliding windows (SIA + VAC) with
 *                          OpenMP-parallel RHS. Thread count is auto-selected
 *                          from N_eq (overridable via OMP_NUM_THREADS); when
 *                          OpenMP is unavailable or the auto-pick lands on 1,
 *                          the same code path simply runs serial.
 *
 * Physics options (physics_option_int):
 *   0 = full_CD_fission      — I+V+2 equations (Case 2, Eq. 175)
 *   1 = full_CD_fusion       — I+2V+1 equations (Case 1, Eq. 174)
 *   2 = bin_moment_CD_fission — 2Ib+V+2 equations (Chapter 9 + Case 2)
 *   3 = bin_moment_CD_fusion  — 2Ib+2V+1 equations (Chapter 9 + Case 1)
 *
 * Output: n_points rows × (1 + N_eq) columns written as raw float64 binary
 * to a companion .bin file, or space-separated text to stdout (fallback).
 *
 * Build:
 *   cd RadCluster_1_0/cpp_utils
 *   cmake -S . -B ../build -DCMAKE_BUILD_TYPE=Release
 *   cmake --build ../build --config Release
 */

#include "parameters.h"
#include "rate_equations.h"

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sundials/sundials_types.h>

#ifdef CD_HAVE_KLU
#  include "sparse_jacobian.h"
#  include <sunlinsol/sunlinsol_klu.h>
#  include <sunmatrix/sunmatrix_sparse.h>
#endif

#include <atomic>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#ifdef CD_HAVE_OPENMP
#  include <omp.h>
#endif

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  ifndef NOMINMAX
#    define NOMINMAX     // prevent windows.h from defining min/max macros
#  endif                 // that collide with std::min / std::max
#  include <windows.h>
#endif

// ── Graceful interrupt handling ──────────────────────────────────────────────
// A SIGINT / SIGTERM (POSIX) or CTRL_C_EVENT / CTRL_BREAK_EVENT (Windows)
// flips this flag.  The integration loop checks it at the top of every
// output-step iteration and breaks cleanly so the post-loop block (final
// CVODE stats, file close) still runs and the partial .bin output is a
// valid prefix of a normal completed run.
namespace {
std::atomic<bool> g_interrupt_requested{false};
}

#ifdef _WIN32
static BOOL WINAPI win_console_handler(DWORD ctrl_type) {
    switch (ctrl_type) {
        case CTRL_C_EVENT:
        case CTRL_BREAK_EVENT:
        case CTRL_CLOSE_EVENT:
        case CTRL_LOGOFF_EVENT:
        case CTRL_SHUTDOWN_EVENT:
            g_interrupt_requested.store(true);
            return TRUE;   // signal handled — don't terminate immediately
    }
    return FALSE;
}
#else
extern "C" void posix_interrupt_handler(int /*sig*/) {
    g_interrupt_requested.store(true);
}
#endif

static void install_interrupt_handlers() {
#ifdef _WIN32
    SetConsoleCtrlHandler(win_console_handler, TRUE);
#else
    std::signal(SIGINT,  posix_interrupt_handler);
    std::signal(SIGTERM, posix_interrupt_handler);
#endif
}

// ── CLI argument parser ────────────────────────────────────────────────────────

static std::map<std::string, double> parse_cli_args(int argc, char* argv[]) {
    std::map<std::string, double> props;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.size() < 3 || arg[0] != '-' || arg[1] != '-') continue;
        auto pos = arg.find('=');
        if (pos == std::string::npos) continue;
        std::string key = arg.substr(2, pos - 2);
        try { props[key] = std::stod(arg.substr(pos + 1)); }
        catch (...) { std::cerr << "Bad value for '" << key << "'\n"; }
    }
    return props;
}

// ── SUNDIALS error check ───────────────────────────────────────────────────────

#define CHECK_SUNDIALS(call) do {                            \
    int _r = (call);                                         \
    if (_r < 0) {                                            \
        std::cerr << "SUNDIALS error " << _r                 \
                  << " at " #call "\n";                      \
        return 1;                                            \
    }                                                        \
} while (0)

// ── Main ──────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {

    install_interrupt_handlers();

    std::map<std::string, double> args;
    std::string param_path;

    if (argc >= 2) {
        const std::string pf = "--param_file=";
        std::string a1 = argv[1];
        if (a1.rfind(pf, 0) == 0) {
            param_path = a1.substr(pf.size());
            args = parse_param_file(param_path);
        } else {
            args = parse_cli_args(argc, argv);
        }
    }
    if (args.empty() && argc > 1) return 1;

    Parameters P = build_parameters(args);
    const int N_EQ = P.N_eq;

    // ── Binary output ─────────────────────────────────────────────────────────
    FILE* fp_bin = nullptr;
    FILE* fp_win = nullptr;   // sidecar CSV: t,x_hi_i,x_hi_v per output point
    if (!param_path.empty()) {
        std::string bp = param_path;
        const std::string suf = ".txt";
        if (bp.size() >= suf.size() &&
            bp.compare(bp.size() - suf.size(), suf.size(), suf) == 0)
            bp.replace(bp.size() - suf.size(), suf.size(), ".bin");
        else bp += ".bin";
#ifdef _MSC_VER
        fopen_s(&fp_bin, bp.c_str(), "wb");
#else
        fp_bin = std::fopen(bp.c_str(), "wb");
#endif
        if (!fp_bin)
            std::cerr << "Warning: cannot open " << bp << " — using stdout\n";

        // Window-bounds sidecar: <bin_path>.window.csv
        const std::string wp = bp + ".window.csv";
#ifdef _MSC_VER
        fopen_s(&fp_win, wp.c_str(), "w");
#else
        fp_win = std::fopen(wp.c_str(), "w");
#endif
        if (fp_win)
            std::fprintf(fp_win, "t,x_hi_i,x_hi_v\n");
    }

    auto write_row = [&](double t, const double* data, int n) {
        if (fp_bin) {
            std::fwrite(&t,   sizeof(double), 1, fp_bin);
            std::fwrite(data, sizeof(double), n, fp_bin);
            std::fflush(fp_bin);  // flush so partial output survives kill/interrupt
        } else {
            std::cout << t;
            for (int k = 0; k < n; ++k) std::cout << ' ' << data[k];
            std::cout << '\n';
        }
    };

    auto write_window = [&](double t, int x_hi_i, int x_hi_v) {
        if (fp_win) {
            std::fprintf(fp_win, "%.17e,%d,%d\n", t, x_hi_i, x_hi_v);
            std::fflush(fp_win);
        }
    };

    // Window mode validation
    constexpr int WINDOW_N_THRESH = 500;
    if (P.window_mode == 4 &&
        P.I < WINDOW_N_THRESH && P.V < WINDOW_N_THRESH) {
        std::cerr << "[Window] I=" << P.I << " and V=" << P.V
                  << " both < threshold=" << WINDOW_N_THRESH
                  << " — using full solver.\n";
        P.window_mode = 0;
    }
#ifndef CD_HAVE_OPENMP
    if (P.window_mode == 4) {
        std::cerr << "[active_window] OpenMP unavailable — running serial (single thread).\n";
    }
#endif

    // Auto-select OMP thread count from problem size.  Per-row work in the
    // hot RHS loops scales like O(i_mobile + v_mobile) × O(N), and the loop
    // length is ~N_eq.  Below ~500 unknowns the parallel overhead dominates;
    // above that, useful threads ≈ N_eq / 1000 (one thread per ~1k rows of
    // work), capped at the runtime maximum.  User can still override via
    // OMP_NUM_THREADS — if it is set in the environment, we honor it as-is.
#ifdef CD_HAVE_OPENMP
    int omp_threads = 1;
#ifdef _MSC_VER
    char* omp_env = nullptr;
    size_t omp_env_len = 0;
    _dupenv_s(&omp_env, &omp_env_len, "OMP_NUM_THREADS");
#else
    const char* omp_env = std::getenv("OMP_NUM_THREADS");
#endif
    if (omp_env && *omp_env) {
        omp_threads = omp_get_max_threads();   // already honored by runtime
        std::cerr << "[OpenMP] OMP_NUM_THREADS=" << omp_env
                  << " honored — using " << omp_threads << " thread"
                  << (omp_threads == 1 ? "" : "s") << "\n";
    } else {
        const int hw_max = omp_get_max_threads();
        int picked;
        if      (N_EQ <   500) picked = 1;
        else if (N_EQ <  2000) picked = 2;
        else if (N_EQ <  5000) picked = 4;
        else if (N_EQ < 10000) picked = 8;
        else if (N_EQ < 20000) picked = 12;
        else if (N_EQ < 40000) picked = 16;
        else if (N_EQ < 80000) picked = 20;
        else                   picked = hw_max;   // largest problems get all cores
        if (picked > hw_max) picked = hw_max;     // cap at machine maximum
        omp_set_num_threads(picked);
        omp_threads = picked;
        std::cerr << "[OpenMP] auto-selected " << omp_threads
                  << " thread" << (omp_threads == 1 ? "" : "s")
                  << " for N_eq=" << N_EQ
                  << " (hw_max=" << hw_max << ")\n";
    }
    // Emit a machine-readable line so the Python provenance can capture
    // the actual thread count chosen by the solver.
    std::cerr << "[OpenMP_threads] " << omp_threads << "\n";
#ifdef _MSC_VER
    if (omp_env) std::free(omp_env);
#endif
#else
    std::cerr << "[OpenMP] not available — running serial (1 thread)\n";
    std::cerr << "[OpenMP_threads] 1\n";
#endif

    const char* he_kin_str = (P.he_kinetics == 1) ? "quasi_steady_state" : "dynamic";
    std::cerr << "RadCluster_1_0 solver: N_eq=" << N_EQ
              << "  physics_option=" << P.physics_option
              << "  he_mode=" << P.he_mode
              << "  he_kinetics=" << he_kin_str
              << "  C_floor=" << P.C_floor
              << "  window_mode=" << P.window_mode;
    if (P.window_mode != 0) {
        std::cerr << "  win_SIA=[0," << std::min(P.window_width - 1, P.I - 1) << "->" << (P.I - 1) << "]"
                  << "  win_VAC=[0," << std::min(P.window_width - 1, P.V - 1) << "->" << (P.V - 1) << "]";
    }
    std::cerr << "\n";

    // ── Time grid ─────────────────────────────────────────────────────────────
    std::vector<double> t_eval(P.n_points);
    if (P.log_time) {
        double l0 = std::log10(P.t_begin), lf = std::log10(P.t_end);
        double step = (lf - l0) / (P.n_points - 1);
        for (int i = 0; i < P.n_points; ++i)
            t_eval[i] = std::pow(10.0, l0 + i * step);
    } else {
        double step = (P.t_end - P.t_begin) / (P.n_points - 1);
        for (int i = 0; i < P.n_points; ++i)
            t_eval[i] = P.t_begin + i * step;
    }

    // ── SUNDIALS context ──────────────────────────────────────────────────────
    SUNContext sunctx;
    if (SUNContext_Create(SUN_COMM_NULL, &sunctx) != 0) {
        std::cerr << "SUNContext_Create failed\n"; return 1;
    }

    // ── State vector ──────────────────────────────────────────────────────────
    N_Vector y = N_VNew_Serial(N_EQ, sunctx);
    if (!y) { std::cerr << "N_VNew_Serial failed\n"; return 1; }
    double* ydata = N_VGetArrayPointer_Serial(y);
    for (int k = 0; k < N_EQ; ++k)
        ydata[k] = std::max(P.y0[k], P.C_floor);

    // ── Offset of the first VAC state in the full state vector ───────────────
    // full_CD:      P.I
    // bin_moment:   P.i_discrete + P.n_mom * P.I_bin
    const int i_VAC_off = (P.physics_option >= 2)
                          ? (P.i_discrete + P.n_mom * P.I_bin)
                          : P.I;
    // Number of VAC state-vector entries per window domain
    const int V_states = (P.physics_option >= 2)
                         ? (P.v_discrete + P.n_mom * P.V_bin)
                         : P.V;

    // ── Select RHS ────────────────────────────────────────────────────────────
    CVRhsFn rhs_fn;
    if (P.physics_option >= 2)
        rhs_fn = rhs_bin_moment;
    else
        rhs_fn = rhs_full_CD;

    // ── User data ─────────────────────────────────────────────────────────────
    UserData ud;
    ud.P             = &P;
    ud.rhs_fn        = rhs_fn;
    ud.x_lo_i        = 0;
    ud.x_hi_i        = P.I - 1;
    ud.x_lo_v        = 0;
    ud.x_hi_v        = V_states - 1;
    ud.window_active = (P.window_mode != 0);

    // ── Create CVODE ──────────────────────────────────────────────────────────
    void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (!cvode_mem) { std::cerr << "CVodeCreate failed\n"; return 1; }

    CHECK_SUNDIALS(CVodeInit(cvode_mem, rhs_fn, P.t_begin, y));
    CHECK_SUNDIALS(CVodeSStolerances(cvode_mem, P.rtol, P.atol));
    CHECK_SUNDIALS(CVodeSetUserData(cvode_mem, &ud));
    CHECK_SUNDIALS(CVodeSetMaxNumSteps(cvode_mem, 500000));
    if (P.max_order > 0)
        CHECK_SUNDIALS(CVodeSetMaxOrd(cvode_mem, P.max_order));
    if (P.hmin > 0.0)
        CHECK_SUNDIALS(CVodeSetMinStep(cvode_mem, P.hmin));

    // Non-negativity is enforced via C_floor clamping at output time (line 317-318).
    // CVODE constraint enforcement (CVodeSetConstraints) is intentionally disabled
    // because it causes Newton corrector failures when near-floor concentrations
    // interact with the stiff SIA–vacancy coupling at moderate doses.

    // ── Linear solver ─────────────────────────────────────────────────────────
    SUNMatrix     sunmat = nullptr;
    SUNLinearSolver sunls = nullptr;

    if (P.linsol == 2) {
        // GMRES is always preconditioned (Jacobi or Woodbury, selected by prec_type).
        constexpr int GMRES_MAXL = 50;
        sunls = SUNLinSol_SPGMR(y, SUN_PREC_RIGHT, GMRES_MAXL, sunctx);
        if (!sunls) { std::cerr << "SPGMR create failed\n"; return 1; }
        CHECK_SUNDIALS(CVodeSetLinearSolver(cvode_mem, sunls, nullptr));
        CHECK_SUNDIALS(CVodeSetPreconditioner(cvode_mem, prec_setup, prec_solve));
        std::cout << "[solver] preconditioner: "
                  << (P.prec_type == 1 ? "Woodbury (bordered-arrow, rank "
                                         + std::to_string(P.prec_rank)
                                         + ", bw " + std::to_string(P.prec_bw) + ")"
                                       : "Jacobi (diagonal)")
                  << "\n";
    } else if (P.linsol == 1) {
        // Band
        int mu = (P.mu > 0) ? P.mu : N_EQ - 1;
        int ml = (P.ml > 0) ? P.ml : N_EQ - 1;
        sunmat = SUNBandMatrix(N_EQ, mu, ml, sunctx);
        sunls  = SUNLinSol_Band(y, sunmat, sunctx);
        if (!sunmat || !sunls) { std::cerr << "Band solver create failed\n"; return 1; }
        CHECK_SUNDIALS(CVodeSetLinearSolver(cvode_mem, sunls, sunmat));
    } else if (P.linsol == 3) {
#ifdef CD_HAVE_KLU
        // KLU sparse direct solver with colored-FD Jacobian.
        // Restricted to full_CD modes (physics_option 0 or 1).
        if (P.physics_option != 0 && P.physics_option != 1) {
            std::cerr << "[KLU] linsol=klu only supports full_CD modes "
                         "(physics_option 0 or 1).  Use gmres for bin_moment.\n";
            return 1;
        }
        int nnz = build_sparsity_pattern_full_CD(P);
        if (nnz <= 0) {
            std::cerr << "[KLU] sparsity pattern build failed\n"; return 1;
        }
        color_columns_greedy(P);
        build_color_groups(P);
        std::cout << "[KLU] N_eq=" << N_EQ
                  << "  nnz=" << nnz
                  << "  density=" << std::scientific << std::setprecision(2)
                  << double(nnz) / (double(N_EQ) * double(N_EQ))
                  << "  colors=" << P.jac_n_colors
                  << "  (~" << P.jac_n_colors + 1
                  << " RHS evals per Jacobian)\n";

        sunmat = SUNSparseMatrix(N_EQ, N_EQ, nnz, CSC_MAT, sunctx);
        if (!sunmat) { std::cerr << "Sparse matrix create failed\n"; return 1; }
        // Install the precomputed CSC pattern (sorted row indices per column).
        sunindextype* col_ptr = SUNSparseMatrix_IndexPointers(sunmat);
        sunindextype* row_idx = SUNSparseMatrix_IndexValues  (sunmat);
        for (int j = 0; j <= N_EQ; ++j) col_ptr[j] = P.jac_col_ptr[j];
        for (int k = 0; k < nnz;     ++k) row_idx[k] = P.jac_row_idx[k];

        sunls = SUNLinSol_KLU(y, sunmat, sunctx);
        if (!sunls) { std::cerr << "KLU create failed\n"; return 1; }
        CHECK_SUNDIALS(CVodeSetLinearSolver(cvode_mem, sunls, sunmat));
        CHECK_SUNDIALS(CVodeSetJacFn(cvode_mem, sparse_fd_jac));
#else
        std::cerr << "[KLU] linsol=klu requested but solver was built without "
                     "KLU support (CD_HAVE_KLU not defined). Rebuild SUNDIALS "
                     "with -DENABLE_KLU=ON and link sunlinsolklu.\n";
        return 1;
#endif
    } else {
        // Dense (default)
        sunmat = SUNDenseMatrix(N_EQ, N_EQ, sunctx);
        sunls  = SUNLinSol_Dense(y, sunmat, sunctx);
        if (!sunmat || !sunls) { std::cerr << "Dense solver create failed\n"; return 1; }
        CHECK_SUNDIALS(CVodeSetLinearSolver(cvode_mem, sunls, sunmat));
    }

    // ── Integration loop ──────────────────────────────────────────────────────
    // Two independent sliding windows:
    //   SIA window: active SIA state indices  x_lo_i .. x_hi_i  (0-based, ≤ I-1)
    //   VAC window: active VAC state indices  x_lo_v .. x_hi_v  (0-based, ≤ V_states-1)
    // window_mode==0: windows span the full domain (no truncation).
    // window_mode==4: start from a user-specified initial width (window_width,
    //   shared by SIA and VAC) and expand independently as the leading
    //   concentration in each axis exceeds concentration_threshold.
    int x_lo_i = 0;
    int x_hi_i = (P.window_mode == 0) ? P.I - 1
                                       : std::min(P.window_width - 1, P.I - 1);
    int x_lo_v = 0;
    int x_hi_v = (P.window_mode == 0) ? V_states - 1
                                       : std::min(P.window_width - 1, V_states - 1);

    // Output buffer
    std::vector<double> out_row(N_EQ);

    // Write initial condition
    for (int k = 0; k < N_EQ; ++k)
        out_row[k] = std::max(ydata[k], 0.0);
    write_row(P.t_begin, out_row.data(), N_EQ);
    write_window(P.t_begin, x_hi_i, x_hi_v);

    int n_written = 1;
    int check_every = std::max(P.window_check_every, 1);
    bool was_interrupted = false;

    for (int i = 1; i < P.n_points; ++i) {
        if (g_interrupt_requested.load()) {
            std::cerr << "[interrupt] received at pt=" << i << "/" << P.n_points
                      << " — finalising " << n_written
                      << " saved time points and exiting cleanly.\n";
            was_interrupted = true;
            break;
        }
        double t_out = t_eval[i];

        // Update both window upper bounds (active_window)
        if (P.window_mode != 0) {
            if ((i - 1) % check_every == 0) {
                // ── SIA window expansion ───────────────────────────────────
                // Expand when the leading SIA cluster concentration exceeds
                // concentration_threshold (absolute index x_hi_i in state vector).
                if (x_hi_i < P.I - 1 && x_hi_i < N_EQ - 1 &&
                    ydata[x_hi_i] > P.concentration_threshold) {
                    x_hi_i = std::min(x_hi_i + P.window_pad, P.I - 1);
                }
                // ── VAC window expansion ───────────────────────────────────
                // ydata[i_VAC_off + x_hi_v] is the leading VAC concentration
                // in the state vector (absolute index i_VAC_off + x_hi_v).
                if (x_hi_v < V_states - 1) {
                    const int vac_abs = i_VAC_off + x_hi_v;
                    if (vac_abs < N_EQ && ydata[vac_abs] > P.concentration_threshold) {
                        x_hi_v = std::min(x_hi_v + P.window_pad_v, V_states - 1);
                    }
                }
            }
            ud.x_hi_i = x_hi_i;
            ud.x_lo_i = x_lo_i;
            ud.x_hi_v = x_hi_v;
            ud.x_lo_v = x_lo_v;
        }

        double t_now = P.t_begin;
        if (i > 0) CVodeGetCurrentTime(cvode_mem, &t_now);

        int retval = CVode(cvode_mem, t_out, y, &t_now, CV_NORMAL);

        // ── Progress diagnostics (always, to stderr) ──────────────────────
        {
            long int nst = 0, nfe = 0, nni = 0, ncfn = 0, netf = 0;
            double hlast = 0.0;
            CVodeGetNumSteps(cvode_mem, &nst);
            CVodeGetNumRhsEvals(cvode_mem, &nfe);
            CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
            CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
            CVodeGetNumErrTestFails(cvode_mem, &netf);
            CVodeGetLastStep(cvode_mem, &hlast);
            std::cerr << "[cvode] pt=" << i << "/" << P.n_points
                      << "  t=" << t_now
                      << "  steps=" << nst << "  rhs=" << nfe
                      << "  nlcf=" << ncfn << "  etf=" << netf
                      << "  h=" << hlast
                      << "  ret=" << retval << "\n";
        }

        if (retval < 0) {
            std::cerr << "CVode failed at t=" << t_out << "  retval=" << retval
                      << " — stopping integration (no reinit to preserve conservation)\n";
            // Do NOT reinit: CVodeReInit resets BDF history, which corrupts the
            // cumulative conservation-accounting ODEs and produces the "Sum > 1"
            // artifacts and solution branch jumps visible in post-processed plots.
            // Instead, stop cleanly so the partial output up to t_now is valid.
            break;
        }

        // Only write output for time points where the solver succeeded.
        // Writing a post-ReInit failure row produces the spike artefacts in plots.
        if (retval >= 0) {
            // Apply C_floor post-step: clamp the CVODE state before writing so
            // output concentrations never fall below the prescribed minimum.
            // This keeps the RHS smooth (no kink above zero) while still
            // enforcing the user-specified floor at every output point.
            for (int k = 0; k < N_EQ; ++k)
                ydata[k] = std::max(ydata[k], P.C_floor);
            for (int k = 0; k < N_EQ; ++k)
                out_row[k] = ydata[k];
            write_row(t_out, out_row.data(), N_EQ);
            write_window(t_out, x_hi_i, x_hi_v);
            ++n_written;

            // ── Progress diagnostics ───────────────────────────────────────
            // Enabled only when verbose=1 is passed in the parameter file.
            if (!P.verbose) continue;

            const int I = P.I;
            const int V = P.V;

            // State vector layout depends on physics_option:
            //   full_CD:      y[0..I-1] = c_i, y[I..I+V-1] = c_v, ...
            //   bin_moment:   y[0..i_d-1] = discrete SIA,
            //                 y[i_d..i_d+PM*Ib-1] = binned SIA moments,
            //                 y[i_VAC..i_VAC+v_d-1] = discrete VAC,
            //                 y[i_VAC+v_d..i_VAC+v_d+PM*Kv-1] = binned VAC moments,
            //                 y[i_Q_base] = Q_tot (case2) or Q_k (case1)
            const bool is_bin = (P.physics_option >= 2);
            const int Ib      = is_bin ? P.I_bin : 0;
            const int Kv      = is_bin ? P.V_bin : 0;
            const int PM      = is_bin ? P.n_mom : 0;
            const int i_d     = is_bin ? P.i_discrete : 0;
            const int v_d     = is_bin ? P.v_discrete : 0;
            const int i_VAC   = is_bin ? (i_d + PM * Ib) : I;
            const int i_Q_base = is_bin ? (i_VAC + v_d + PM * Kv)
                                        : (I + V);

            // Reconstruct c_i1 from bin-moment state (bin 0 = monomer only)
            const double c_i1 = ydata[0];
            const double c_v1 = ydata[i_VAC];

            // Total SIA content Sigma n*c_n
            double SIA_content = 0.0;
            if (is_bin) {
                // Discrete SIA sizes: Sigma n*c_n for n=1..i_discrete
                for (int n = 0; n < i_d; ++n)
                    SIA_content += (n + 1.0) * ydata[n];
                // Binned SIA: first moments mu_k^{(1)}
                for (int k = 0; k < Ib; ++k)
                    SIA_content += ydata[i_d + PM*k + 1];
            } else {
                for (int n = 0; n < I; ++n)
                    SIA_content += (n + 1.0) * ydata[n];
            }

            // Total vacancy content Sigma m*c_m
            double VAC_content = 0.0;
            if (is_bin) {
                // Discrete VAC sizes: Sigma m*c_m for m=1..v_discrete
                for (int m = 0; m < v_d; ++m)
                    VAC_content += (m + 1.0) * ydata[i_VAC + m];
                // Binned VAC: first moments (if PM>=2) or midpoint approx
                if (Kv > 0) {
                    const int vac_mom = i_VAC + v_d;
                    if (PM >= 2) {
                        for (int k = 0; k < Kv; ++k)
                            VAC_content += ydata[vac_mom + PM*k + 1];
                    } else {
                        int edge = v_d + 1;
                        for (int k = 0; k < Kv; ++k) {
                            int mlo = edge;
                            int next = std::max(static_cast<int>(std::floor(edge * P.r_ratio)), edge + 1);
                            int mhi  = std::min(next, V + 1);
                            double mid = 0.5 * (mlo + mhi - 1);
                            VAC_content += ydata[vac_mom + k] * mid;
                            edge = mhi;
                        }
                    }
                }
            } else {
                for (int m = 0; m < V; ++m)
                    VAC_content += (m + 1.0) * ydata[I + m];
            }

            // Representative cluster sizes
            const double c_i2  = is_bin ? ((i_d > 1) ? ydata[1] : 0.0)
                                        : ((I > 1) ? ydata[1] : 0.0);
            const double c_i5  = is_bin ? ((i_d > 4) ? ydata[4] : 0.0)
                                        : ((I > 4) ? ydata[4] : 0.0);
            const double c_v2  = is_bin ? ((v_d > 1) ? ydata[i_VAC + 1] : 0.0)
                                        : ((V > 1) ? ydata[I + 1] : 0.0);
            const double c_v5  = is_bin ? ((v_d > 4) ? ydata[i_VAC + 4] : 0.0)
                                        : ((V > 4) ? ydata[I + 4] : 0.0);

            // He
            const double Q_tot = ydata[i_Q_base];
            const double c_h   = (P.he_kinetics == 0 && N_EQ > i_Q_base + 1)
                                 ? ydata[i_Q_base + 1] : -1.0;

            std::cerr << std::scientific << std::setprecision(3)
                      << "  [diag] t=" << t_out
                      << "  c_i1=" << c_i1
                      << "  c_v1=" << c_v1
                      << "  c_i2=" << c_i2
                      << "  c_i5=" << c_i5
                      << "  c_v2=" << c_v2
                      << "  c_v5=" << c_v5
                      << "  Q_tot=" << Q_tot;
            if (c_h >= 0.0)
                std::cerr << "  c_h=" << c_h;
            std::cerr << "  SIA_tot=" << SIA_content
                      << "  VAC_tot=" << VAC_content
                      << "\n";

            // ── C_i5 reaction rate breakdown (n=5, index 4) ───────────────
            // For bin_moment mode, c_i/c_v indices differ from full_CD
            if (!is_bin && I > 4) {
                const double ci4 = std::max(ydata[3], 0.0);
                const double ci5_ = c_i5;
                const double ci6 = (I > 5) ? std::max(ydata[5], 0.0) : 0.0;

                const double rate_prod       =  P.Pr_SIA[4];
                const double rate_emit_gain  =  (I > 5) ? P.GII[5] * ci6 : 0.0;
                const double rate_emit_loss  = -P.GII[4] * ci5_;
                const double rate_grow_gain  =  P.KII[3] * c_i1 * ci4;
                const double rate_grow_loss  = -P.KII[4] * c_i1 * ci5_;
                const double rate_shrink_gain = (I > 5) ? P.KIV[5] * c_v1 * ci6 : 0.0;
                const double rate_shrink_loss = -P.KIV[4] * c_v1 * ci5_;

                double rate_1D_loss = 0.0;
                if (4 < P.i_mobile && P.K_1D_pref[4] > 1e-300) {
                    for (int m = 0; m < V; ++m) {
                        const double m_f  = static_cast<double>(m + 1);
                        const double m13  = std::cbrt(m_f);
                        const double denom = 1.0 + P.B_rot * P.L_hat * P.L_hat / m13;
                        const double k1d  = P.K_1D_pref[4] * m13 / denom;
                        rate_1D_loss -= k1d * ci5_ * std::max(ydata[i_VAC + m], 0.0);
                    }
                }

                const double rate_sink = -P.k2_SIA[4] * ci5_;
                const double dc_i5_total = rate_prod + rate_emit_gain + rate_emit_loss
                                         + rate_grow_gain + rate_grow_loss
                                         + rate_shrink_gain + rate_shrink_loss
                                         + rate_1D_loss + rate_sink;

                std::cerr << std::scientific << std::setprecision(3)
                          << "  [ci5_rates] t=" << t_out
                          << "  prod=" << rate_prod
                          << "  emit_in=" << rate_emit_gain
                          << "  emit_out=" << rate_emit_loss
                          << "  grow_in=" << rate_grow_gain
                          << "  grow_out=" << rate_grow_loss
                          << "  shrink_in=" << rate_shrink_gain
                          << "  shrink_out=" << rate_shrink_loss
                          << "  1D_loss=" << rate_1D_loss
                          << "  sink=" << rate_sink
                          << "  dc_i5=" << dc_i5_total
                          << "\n";
            }

            // ── C_v5 reaction rate breakdown (m=5, index 4) ───────────────
            if (!is_bin && V > 4) {
                const double cv4_ = (V > 3) ? std::max(ydata[i_VAC+3], 0.0) : 0.0;
                const double cv5_ = c_v5;
                const double cv6_ = (V > 5) ? std::max(ydata[i_VAC+5], 0.0) : 0.0;

                const double rate_prod        =  P.Pr_VAC[4];
                const double rate_emit_gain   =  (V > 5) ? P.GVV[5] * cv6_ : 0.0;
                const double rate_emit_loss   = -P.GVV[4] * cv5_;
                const double rate_grow_gain   =  P.KVV[3] * c_v1 * cv4_;
                const double rate_grow_loss   = -P.KVV[4] * c_v1 * cv5_;
                const double rate_shrink_gain =  (V > 5) ? P.KVI[5] * c_i1 * cv6_ : 0.0;
                const double rate_shrink_loss = -P.KVI[4] * c_i1 * cv5_;

                double rate_1D_loss = 0.0;
                if (!is_bin) {
                    const double m5f  = 5.0;
                    const double m513 = std::cbrt(m5f);
                    for (int n = 4; n < std::min(I, P.i_mobile); ++n) {
                        if (P.K_1D_pref[n] < 1e-300) continue;
                        const double denom = 1.0 + P.B_rot * P.L_hat * P.L_hat / m513;
                        const double k1d   = P.K_1D_pref[n] * m513 / denom;
                        rate_1D_loss -= k1d * std::max(ydata[n], 0.0) * cv5_;
                    }
                }

                const double rate_sink    = -P.k2_disl_v * cv5_;
                const double dc_v5_total  = rate_prod + rate_emit_gain + rate_emit_loss
                                          + rate_grow_gain + rate_grow_loss
                                          + rate_shrink_gain + rate_shrink_loss
                                          + rate_1D_loss + rate_sink;

                std::cerr << std::scientific << std::setprecision(3)
                          << "  [cv5_rates] t=" << t_out
                          << "  prod=" << rate_prod
                          << "  emit_in=" << rate_emit_gain
                          << "  emit_out=" << rate_emit_loss
                          << "  grow_in=" << rate_grow_gain
                          << "  grow_out=" << rate_grow_loss
                          << "  shrink_in=" << rate_shrink_gain
                          << "  shrink_out=" << rate_shrink_loss
                          << "  1D_loss=" << rate_1D_loss
                          << "  sink=" << rate_sink
                          << "  dc_v5=" << dc_v5_total
                          << "\n";
            }
        }
    }

    // ── Final integrator/linsol statistics ───────────────────────────────────
    // Emit a single machine-readable line so the Python bridge can capture
    // GMRES iteration count + preconditioner activity into provenance.
    {
        long nst = 0, nfe = 0, nni = 0, ncfn = 0, netf = 0;
        long nli = 0, npe = 0, nps = 0, nlsetup = 0;
        CVodeGetNumSteps              (cvode_mem, &nst);
        CVodeGetNumRhsEvals           (cvode_mem, &nfe);
        CVodeGetNumNonlinSolvIters    (cvode_mem, &nni);
        CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
        CVodeGetNumErrTestFails       (cvode_mem, &netf);
        CVodeGetNumLinSolvSetups      (cvode_mem, &nlsetup);
        // GMRES-only stats — return 0 / no-op for dense/band linsols
        CVodeGetNumLinIters (cvode_mem, &nli);
        CVodeGetNumPrecEvals(cvode_mem, &npe);
        CVodeGetNumPrecSolves(cvode_mem, &nps);
        const double nli_per_nni = (nni > 0) ? double(nli) / double(nni) : 0.0;
        std::cerr << "[stats] steps=" << nst
                  << " nfe="  << nfe
                  << " nni="  << nni
                  << " nli="  << nli
                  << " nli_per_nni=" << std::fixed << std::setprecision(2)
                  << nli_per_nni
                  << " npe="  << npe
                  << " nps="  << nps
                  << " ncfn=" << ncfn
                  << " netf=" << netf
                  << " nlsetup=" << nlsetup
                  << " interrupted=" << (was_interrupted ? 1 : 0)
                  << "\n";
    }

    std::cerr << "Done: " << n_written << " time points written"
              << (was_interrupted ? " (interrupted)" : "") << ".\n";

    // ── Cleanup ───────────────────────────────────────────────────────────────
    if (fp_bin) std::fclose(fp_bin);
    if (fp_win) std::fclose(fp_win);
    N_VDestroy_Serial(y);
    CVodeFree(&cvode_mem);
    if (sunls)  SUNLinSolFree(sunls);
    if (sunmat) SUNMatDestroy(sunmat);
    SUNContext_Free(&sunctx);

    return 0;
}
