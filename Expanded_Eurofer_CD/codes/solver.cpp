/**
 * solver.cpp — Expanded_Eurofer_CD main C++ ODE solver.
 *
 * Drives the Expanded_Eurofer_CD cluster dynamics system for bcc Fe / EUROFER97.
 * Invoked by py_utils/cpp_bridge.py via a parameter file:
 *   solver.exe --param_file=<path>
 *
 * Physics reference: Ghoniem (2026), Rate_Equations.pdf, Eqs. 152, 155, 157,
 *   Chapter 9 (bin moments), Section 8 (He-reduction).
 *
 * Solver modes (window_mode parameter):
 *   0 = cpp_full        — full system, CVODE BDF
 *   3 = cpp_sliding_win — Phase III: constant-width sliding window on SIA
 *   4 = sliding_OpenMP  — Phase III + OpenMP intra-RHS parallelism
 *
 * Physics options (physics_option_int):
 *   0 = full_CD_fission      — N+M+2 equations (Case 2, Eq. 175)
 *   1 = full_CD_fusion       — N+2M+1 equations (Case 1, Eq. 174)
 *   2 = bin_moment_CD_fission — 2K+M+2 equations (Chapter 9 + Case 2)
 *   3 = bin_moment_CD_fusion  — 2K+2M+1 equations (Chapter 9 + Case 1)
 *
 * Output: n_points rows × (1 + N_eq) columns written as raw float64 binary
 * to a companion .bin file, or space-separated text to stdout (fallback).
 *
 * Build:
 *   cd Expanded_Eurofer_CD/cpp_utils
 *   cmake -S . -B ../build -DCMAKE_BUILD_TYPE=Release
 *   cmake --build ../build --config Release
 */

#include "parameters.h"
#include "rate_equations.h"

#include <cvode/cvode.h>
#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sundials/sundials_types.h>

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#ifdef CD_HAVE_OPENMP
#  include <omp.h>
#endif

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
    }

    auto write_row = [&](double t, const double* data, int n) {
        if (fp_bin) {
            std::fwrite(&t,   sizeof(double), 1, fp_bin);
            std::fwrite(data, sizeof(double), n, fp_bin);
        } else {
            std::cout << t;
            for (int k = 0; k < n; ++k) std::cout << ' ' << data[k];
            std::cout << '\n';
        }
    };

    // Window mode validation
    if ((P.window_mode == 3 || P.window_mode == 4) &&
        P.N < P.window_N_thresh) {
        std::cerr << "[Window] N=" << P.N << " < threshold=" << P.window_N_thresh
                  << " — using full solver.\n";
        P.window_mode = 0;
    }
#ifndef CD_HAVE_OPENMP
    if (P.window_mode == 4) {
        std::cerr << "[sliding_OpenMP] OpenMP unavailable — using cpp_sliding_win.\n";
        P.window_mode = 3;
    }
#endif

    const char* he_opts_str = (P.he_options == 1) ? "quasi_steady_state" : "dynamic";
    std::cerr << "Expanded_Eurofer_CD solver: N_eq=" << N_EQ
              << "  physics_option=" << P.physics_option
              << "  he_mode=" << P.he_mode
              << "  he_options=" << he_opts_str
              << "  C_floor=" << P.C_floor
              << "  window_mode=" << P.window_mode << "\n";

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

    // ── User data ─────────────────────────────────────────────────────────────
    UserData ud;
    ud.P            = &P;
    ud.x_lo_i       = 0;
    ud.x_hi_i       = P.N - 1;
    ud.x_hi_v       = P.M - 1;
    ud.window_active = (P.window_mode != 0);

    // ── Select RHS ────────────────────────────────────────────────────────────
    CVRhsFn rhs_fn;
    if (P.physics_option >= 2)
        rhs_fn = rhs_bin_moment;
    else
        rhs_fn = rhs_full_CD;

    // ── Create CVODE ──────────────────────────────────────────────────────────
    void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (!cvode_mem) { std::cerr << "CVodeCreate failed\n"; return 1; }

    CHECK_SUNDIALS(CVodeInit(cvode_mem, rhs_fn, P.t_begin, y));
    CHECK_SUNDIALS(CVodeSStolerances(cvode_mem, P.rtol, P.atol));
    CHECK_SUNDIALS(CVodeSetUserData(cvode_mem, &ud));
    CHECK_SUNDIALS(CVodeSetMaxNumSteps(cvode_mem, 100000));
    if (P.max_order > 0)
        CHECK_SUNDIALS(CVodeSetMaxOrd(cvode_mem, P.max_order));

    // ── Non-negativity constraints ─────────────────────────────────────────────
    // Enforce c_n, c_m ≥ 0 natively in CVODE's nonlinear solver so it never
    // proposes negative concentrations. This eliminates the primary trigger for
    // the step-failure / CVodeReInit spike pattern seen in the solution output.
    {
        N_Vector constraints = N_VNew_Serial(N_EQ, sunctx);
        N_VConst(1.0, constraints);   // 1.0 = non-negative for every component
        int cret = CVodeSetConstraints(cvode_mem, constraints);
        N_VDestroy_Serial(constraints);
        if (cret < 0)
            std::cerr << "Warning: CVodeSetConstraints failed (retval=" << cret
                      << ") — continuing without constraint enforcement\n";
    }

    // ── Linear solver ─────────────────────────────────────────────────────────
    SUNMatrix     sunmat = nullptr;
    SUNLinearSolver sunls = nullptr;

    if (P.linsol == 2) {
        // GMRES
        sunls = SUNLinSol_SPGMR(y, SUN_PREC_RIGHT,
                                  P.window_gmres_maxl > 0 ? P.window_gmres_maxl : 20,
                                  sunctx);
        if (!sunls) { std::cerr << "SPGMR create failed\n"; return 1; }
        CHECK_SUNDIALS(CVodeSetLinearSolver(cvode_mem, sunls, nullptr));
        if (P.window_prec) {
            CHECK_SUNDIALS(CVodeSetPreconditioner(cvode_mem, prec_setup, prec_solve));
        }
    } else if (P.linsol == 1) {
        // Band
        int mu = (P.mu > 0) ? P.mu : N_EQ - 1;
        int ml = (P.ml > 0) ? P.ml : N_EQ - 1;
        sunmat = SUNBandMatrix(N_EQ, mu, ml, sunctx);
        sunls  = SUNLinSol_Band(y, sunmat, sunctx);
        if (!sunmat || !sunls) { std::cerr << "Band solver create failed\n"; return 1; }
        CHECK_SUNDIALS(CVodeSetLinearSolver(cvode_mem, sunls, sunmat));
    } else {
        // Dense (default)
        sunmat = SUNDenseMatrix(N_EQ, N_EQ, sunctx);
        sunls  = SUNLinSol_Dense(y, sunmat, sunctx);
        if (!sunmat || !sunls) { std::cerr << "Dense solver create failed\n"; return 1; }
        CHECK_SUNDIALS(CVodeSetLinearSolver(cvode_mem, sunls, sunmat));
    }

    // ── Integration loop ──────────────────────────────────────────────────────
    // Window state (applies to SIA cluster indices 0..N-1)
    int x_lo_i = 0;
    int x_hi_i = (P.window_mode == 0) ? P.N - 1 : std::min(P.window_w0_i - 1, P.N - 1);

    // Output buffer
    std::vector<double> out_row(N_EQ);

    // Write initial condition
    for (int k = 0; k < N_EQ; ++k)
        out_row[k] = std::max(ydata[k], 0.0);
    write_row(P.t_begin, out_row.data(), N_EQ);

    int n_written = 1;
    int check_every = std::max(P.window_check_every, 1);

    for (int i = 1; i < P.n_points; ++i) {
        double t_out = t_eval[i];

        // Update window upper bound (cpp_sliding_win / sliding_OpenMP)
        if (P.window_mode != 0) {
            // Check if leading SIA concentration exceeds expand threshold
            if ((i - 1) % check_every == 0 && x_hi_i < P.N - 1) {
                if (x_hi_i < N_EQ - 1 && ydata[x_hi_i] > P.window_C_expand) {
                    x_hi_i = std::min(x_hi_i + P.window_expand_pad, P.N - 1);
                }
            }
            ud.x_hi_i = x_hi_i;
            ud.x_lo_i = x_lo_i;
        }

        double t_now = P.t_begin;
        if (i > 0) CVodeGetCurrentTime(cvode_mem, &t_now);

        int retval = CVode(cvode_mem, t_out, y, &t_now, CV_NORMAL);
        if (retval < 0) {
            std::cerr << "CVode failed at t=" << t_out << "  retval=" << retval
                      << " — reinitialising at t=" << t_now << "\n";
            // Clamp negatives before reinit so the RHS starts from a valid state
            for (int k = 0; k < N_EQ; ++k)
                ydata[k] = std::max(ydata[k], P.C_floor);
            // CVodeReInit resets BDF order to 1 and clears step history,
            // allowing the solver to bootstrap through a sharp transient.
            int ri = CVodeReInit(cvode_mem, t_now, y);
            if (ri < 0) {
                std::cerr << "CVodeReInit failed (retval=" << ri << ") — skipping point\n";
            } else {
                // Reset initial step size: heuristic h0 = (t_out - t_now) / 1000
                // so BDF restarts at order 1 with a small but not hmin step.
                double h0 = std::max((t_out - t_now) * 1e-3, P.rtol * t_now * 1e-3);
                CVodeSetInitStep(cvode_mem, h0);
                retval = CVode(cvode_mem, t_out, y, &t_now, CV_NORMAL);
                if (retval < 0)
                    std::cerr << "CVode still failed after reinit at t=" << t_out
                              << "  retval=" << retval << "\n";
            }
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
            ++n_written;

            // ── Progress diagnostics ───────────────────────────────────────
            // Print key concentrations at every output point so the user can
            // monitor evolution without waiting for post-processing.
            //
            // Layout (Case 2 / fission):
            //   y[0..N-1]   = c_i (SIA clusters, n=1..N)
            //   y[N..N+M-1] = c_v (vacancy clusters, m=1..M)
            //   y[N+M]      = Q_tot (total He in voids)
            //   y[N+M+1]    = c_h  (free He, dynamic mode only)
            const int N = P.N;
            const int M = P.M;

            // Point-defect monomers
            const double c_i1 = ydata[0];
            const double c_v1 = ydata[N];
            // Total SIA content Σ n·c_n  and void content Σ m·c_m
            double SIA_content = 0.0, VAC_content = 0.0;
            for (int n = 0; n < N; ++n) SIA_content += (n + 1.0) * ydata[n];
            for (int m = 0; m < M; ++m) VAC_content += (m + 1.0) * ydata[N + m];
            // A few representative cluster sizes
            const double c_i2  = (N > 1) ? ydata[1]   : 0.0;  // n=2 SIA cluster
            const double c_i5  = (N > 4) ? ydata[4]   : 0.0;  // n=5 SIA cluster
            const double c_v2  = (M > 1) ? ydata[N+1] : 0.0;  // m=2 void
            const double c_v5  = (M > 4) ? ydata[N+4] : 0.0;  // m=5 void
            // He
            const double Q_tot = ydata[N + M];
            const double c_h   = (P.he_options == 0 && N_EQ > N+M+1)
                                 ? ydata[N + M + 1] : -1.0;    // -1 = QSS mode

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
        }

        if (i % (P.n_points / 10 + 1) == 0)
            std::cerr << "  t=" << t_out << "  n=" << n_written
                      << "  x_hi_i=" << x_hi_i << "\n";
    }

    std::cerr << "Done: " << n_written << " time points written.\n";

    // ── Cleanup ───────────────────────────────────────────────────────────────
    if (fp_bin) std::fclose(fp_bin);
    N_VDestroy_Serial(y);
    CVodeFree(&cvode_mem);
    if (sunls)  SUNLinSolFree(sunls);
    if (sunmat) SUNMatDestroy(sunmat);
    SUNContext_Free(&sunctx);

    return 0;
}
