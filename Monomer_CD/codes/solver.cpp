/**
 * solver.cpp – ClusterDynamics main C++ ODE solver.
 *
 * Drives a runtime-sized cluster dynamics system (Nv vacancy + Ni interstitial
 * cluster sizes, both set at run time via the parameter file).
 *
 * Invoked by py_utils/cpp_bridge.py via a parameter file to avoid Windows
 * command-line length limits (WinError 206) with large Nv/Ni:
 *
 *   solver.exe --param_file=<path>
 *
 * The parameter file contains one "key=value" entry per line (written by
 * cpp_bridge.write_param_file).  See parameters.h for the full key list.
 *
 * Legacy CLI invocation (--key=value ...) is still supported for small systems.
 *
 * Backend (--backend in param file):
 *   0 = CVODE    — linear multistep BDF/Adams (default; best for stiff systems)
 *   1 = ARKODE   — implicit Runge-Kutta DIRK
 *
 * Output: n_points rows × (1 + N_EQ) columns, space-separated, scientific:
 *   t  Cv1 Cv2 ... Cv_Nv  Ci1 Ci2 ... Ci_Ni
 *
 * Build:
 *   cd ClusterDynamics/cpp_utils
 *   cmake -S . -B ../build -DCMAKE_BUILD_TYPE=Release
 *   cmake --build ../build --config Release
 */

#include "parameters.h"
#include "rate_equations.h"

#include <cvodes/cvodes.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_butcher_dirk.h>
#include <arkode/arkode_butcher_erk.h>
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

// ── CLI argument parser (legacy --key=value ... interface) ────────────────────

std::map<std::string, double> parse_args(int argc, char* argv[]) {
    std::map<std::string, double> props;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.size() < 3 || arg[0] != '-' || arg[1] != '-') {
            std::cerr << "Unexpected argument: " << arg << "\n";
            return {};
        }
        auto pos = arg.find('=');
        if (pos == std::string::npos) {
            std::cerr << "Invalid argument format (missing '='): " << arg << "\n";
            return {};
        }
        std::string key = arg.substr(2, pos - 2);
        double val = 0.0;
        try {
            val = std::stod(arg.substr(pos + 1));
        } catch (...) {
            std::cerr << "Invalid numeric value for " << key << "\n";
            return {};
        }
        props[key] = val;
    }
    return props;
}

// ── Main ───────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {

    std::map<std::string, double> args;
    std::string param_path_str;   // remembered for deriving the binary output path

    // Detect --param_file=<path> (single-argument file-based invocation)
    if (argc == 2) {
        std::string a1 = argv[1];
        const std::string prefix = "--param_file=";
        if (a1.rfind(prefix, 0) == 0) {
            param_path_str = a1.substr(prefix.size());
            args = parse_param_file(param_path_str);
        } else {
            args = parse_args(argc, argv);
        }
    } else {
        args = parse_args(argc, argv);
    }

    if (args.empty() && argc > 1) return 1;

    Parameters P = build_parameters(args);
    const int N_EQ = P.N_EQ;

    // ── Binary output file ─────────────────────────────────────────────────────
    // Derive the output path from the param file path (replace .txt → .bin).
    // If no param file was used (legacy CLI), fp_bin stays null and we fall back
    // to the original text-stdout path.
    FILE* fp_bin = nullptr;
    std::string bin_path;
    if (!param_path_str.empty()) {
        bin_path = param_path_str;
        const std::string suffix = ".txt";
        if (bin_path.size() >= suffix.size() &&
            bin_path.compare(bin_path.size() - suffix.size(), suffix.size(), suffix) == 0)
            bin_path.replace(bin_path.size() - suffix.size(), suffix.size(), ".bin");
        else
            bin_path += ".bin";
#ifdef _MSC_VER
        fopen_s(&fp_bin, bin_path.c_str(), "wb");
#else
        fp_bin = std::fopen(bin_path.c_str(), "wb");
#endif
        if (!fp_bin)
            std::cerr << "Warning: cannot open binary output " << bin_path
                      << " — falling back to text stdout\n";
    }

    // Helper: write one output row [t, data[0..n-1]] as raw doubles, or as
    // space-separated text if the binary file is unavailable.
    auto write_row = [&](double t, const double* data, int n) {
        if (fp_bin) {
            std::fwrite(&t,    sizeof(double), 1, fp_bin);
            std::fwrite(data,  sizeof(double), n, fp_bin);
        } else {
            std::cout << t;
            for (int k = 0; k < n; ++k) std::cout << ' ' << data[k];
            std::cout << '\n';
        }
    };

    // Phase III/IV: fall back to full solver when system is too small
    if ((P.window_mode == 3 || P.window_mode == 4) && P.N_EQ <= P.window_N_thresh) {
        std::cerr << "[Phase III/IV] N_EQ=" << P.N_EQ << " <= N_thresh=" << P.window_N_thresh
                  << " — using full solver\n";
        P.window_mode = 0;
    }

#ifndef CD_HAVE_OPENMP
    if (P.window_mode == 4) {
        std::cerr << "[Phase IV] OpenMP not available in this build "
                     "(rebuild with brew install libomp then re-run cmake).\n"
                     "Falling back to Phase III (window_mode=3).\n";
        P.window_mode = 3;
    }
#endif

    // ── Build time evaluation grid ─────────────────────────────────────────────
    std::vector<double> t_eval(P.n_points);
    if (P.log_time) {
        double log_t0 = std::log10(P.t_begin);
        double log_tf = std::log10(P.t_end);
        double step   = (log_tf - log_t0) / (P.n_points - 1);
        for (int i = 0; i < P.n_points; ++i)
            t_eval[i] = std::pow(10.0, log_t0 + i * step);
    } else {
        double step = (P.t_end - P.t_begin) / (P.n_points - 1);
        for (int i = 0; i < P.n_points; ++i)
            t_eval[i] = P.t_begin + i * step;
    }

    // ── SUNDIALS context ────────────────────────────────────────────────────────
    SUNContext sunctx;
    if (SUNContext_Create(SUN_COMM_NULL, &sunctx) != 0) {
        std::cerr << "Error creating SUNContext\n";
        return 1;
    }

    // ── State vector — initialised from Python-computed y0 ────────────────────
    N_Vector y = N_VNew_Serial(N_EQ, sunctx);
    if (!y) { std::cerr << "Error allocating N_Vector\n"; return 1; }
    for (int k = 0; k < N_EQ; ++k)
        NV_Ith_S(y, k) = P.y0[k];

    int flag;
    SUNMatrix       A  = nullptr;
    SUNLinearSolver LS = nullptr;
    const char* linsol_name = "dense";

    // Shared linear-solver factory (used by both CVODE and ARKODE)
    auto make_linear_solver = [&](void* solver_mem, bool is_arkode) -> bool {
        if (P.linsol == 1) {
            linsol_name = "band";
            A  = SUNBandMatrix(N_EQ, P.mu, P.ml, sunctx);
            if (!A)  { std::cerr << "Error creating SUNBandMatrix\n"; return false; }
            LS = SUNLinSol_Band(y, A, sunctx);
            if (!LS) { std::cerr << "Error creating SUNLinSol_Band\n"; return false; }
        } else if (P.linsol == 2) {
            linsol_name = "gmres";
            LS = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx);
            if (!LS) { std::cerr << "Error creating SUNLinSol_SPGMR\n"; return false; }
            A  = nullptr;
        } else {
            A  = SUNDenseMatrix(N_EQ, N_EQ, sunctx);
            if (!A)  { std::cerr << "Error creating SUNDenseMatrix\n"; return false; }
            LS = SUNLinSol_Dense(y, A, sunctx);
            if (!LS) { std::cerr << "Error creating SUNLinSol_Dense\n"; return false; }
        }
        int r = is_arkode ? ARKodeSetLinearSolver(solver_mem, LS, A)
                          : CVodeSetLinearSolver(solver_mem, LS, A);
        if (r != 0) { std::cerr << "Error in SetLinearSolver\n"; return false; }
        return true;
    };

    std::cout << std::scientific << std::setprecision(10);

    if (P.backend == 1) {
        // ── ARKODE ARKStep — implicit Runge-Kutta (DIRK) ──────────────────────
        void* ark_mem = ARKStepCreate(nullptr, rhs_cd, t_eval[0], y, sunctx);
        if (!ark_mem) { std::cerr << "Error in ARKStepCreate\n"; return 1; }

        flag = ARKStepSetImplicit(ark_mem);
        if (flag != ARK_SUCCESS) { std::cerr << "Error in ARKStepSetImplicit\n"; return 1; }

        flag = ARKStepSetTableNum(ark_mem,
                                  static_cast<ARKODE_DIRKTableID>(P.ark_table),
                                  ARKODE_ERK_NONE);
        if (flag != ARK_SUCCESS) { std::cerr << "Error in ARKStepSetTableNum\n"; return 1; }

        flag = ARKodeSetUserData(ark_mem, &P);
        if (flag != ARK_SUCCESS) { std::cerr << "Error in ARKodeSetUserData\n"; return 1; }

        flag = ARKodeSStolerances(ark_mem, P.rtol, P.atol);
        if (flag != ARK_SUCCESS) { std::cerr << "Error in ARKodeSStolerances\n"; return 1; }

        flag = ARKodeSetMaxNumSteps(ark_mem, 1000000);
        if (flag != ARK_SUCCESS) { std::cerr << "Error in ARKodeSetMaxNumSteps\n"; return 1; }

        if (P.max_order > 0) {
            flag = ARKodeSetOrder(ark_mem, P.max_order);
            if (flag != ARK_SUCCESS) { std::cerr << "Error in ARKodeSetOrder\n"; return 1; }
        }

        if (!make_linear_solver(ark_mem, true)) return 1;

        const char* tbl_name = ARKodeButcherTable_DIRKIDToName(
                                    static_cast<ARKODE_DIRKTableID>(P.ark_table));
        std::cerr << "Solver: ARKODE ARKStep DIRK"
                  << " | table: " << (tbl_name ? tbl_name : "?")
                  << " (id=" << P.ark_table << ")"
                  << " | linsol: " << linsol_name
                  << " | N_EQ=" << N_EQ << "\n";

        // Initial point
        write_row(t_eval[0], NV_DATA_S(y), N_EQ);

        sunrealtype t_current = t_eval[0];
        for (int i = 1; i < P.n_points; ++i) {
            flag = ARKodeEvolve(ark_mem, t_eval[i], y, &t_current, ARK_NORMAL);
            if (flag < 0) {
                std::cerr << "ARKodeEvolve error at step " << i
                          << " (t=" << t_eval[i] << "): flag=" << flag << "\n";
                N_VDestroy(y); SUNLinSolFree(LS); if (A) SUNMatDestroy(A);
                ARKodeFree(&ark_mem); SUNContext_Free(&sunctx);
                if (fp_bin) std::fclose(fp_bin);
                return 1;
            }
            write_row(t_eval[i], NV_DATA_S(y), N_EQ);
        }

        N_VDestroy(y); SUNLinSolFree(LS); if (A) SUNMatDestroy(A);
        ARKodeFree(&ark_mem);

    } else if (P.window_mode == 1) {
        // ── CVODE — Dynamic window solver (Phase I: upper truncation) ─────────
        //
        // Only clusters [1..x_hi_v] (vacancy) and [1..x_hi_i] (interstitial) are
        // integrated by CVODE.  When the concentration at the upper boundary
        // exceeds window_C_expand the window is expanded and CVODE is reinitialised
        // with the new (larger) state vector.  The full Nv+Ni row is always written
        // to stdout; truncated species are held at zero.
        //
        // rhs_cd is reused unchanged: a P_win copy with Nv=x_hi_v / Ni=x_hi_i
        // drives the correct truncated loops with no modifications to rate_equations.
        // GMRES (matrix-free) is used so the linear solver needs no resizing.

        int lmm_win = (P.lmm == 1) ? CV_ADAMS : CV_BDF;
        const char* lmm_win_name = (P.lmm == 1) ? "ADAMS" : "BDF";

        // Current window upper bounds (number of active clusters, 1-indexed)
        int x_hi_v = std::min(std::max(P.window_w0_v, 2), P.Nv);
        int x_hi_i = std::min(std::max(P.window_w0_i, 2), P.Ni);

        // P_win: a copy of P with reduced Nv/Ni/N_EQ; updated in-place on expansion.
        // Must outlive every CVode call (passed as user_data pointer).
        Parameters P_win = P;
        P_win.Nv   = x_hi_v;
        P_win.Ni   = x_hi_i;
        P_win.N_EQ = x_hi_v + x_hi_i;

        // Full concentration buffer (size Nv+Ni); output always uses this.
        std::vector<double> full_conc(P.N_EQ, 0.0);
        for (int k = 0; k < P.N_EQ; ++k) full_conc[k] = P.y0[k];

        // Pack active window [Cv1..Cv_{x_hi_v}, Ci1..Ci_{x_hi_i}] from full_conc.
        // Layout matches rhs_cd when called with P_win.Nv=x_hi_v, P_win.Ni=x_hi_i.
        auto pack_active = [&](N_Vector ya) {
            for (int k = 0; k < x_hi_v; ++k)
                NV_Ith_S(ya, k) = full_conc[k];
            for (int k = 0; k < x_hi_i; ++k)
                NV_Ith_S(ya, x_hi_v + k) = full_conc[P.Nv + k];
        };
        auto unpack_active = [&](N_Vector ya) {
            for (int k = 0; k < x_hi_v; ++k)
                full_conc[k] = NV_Ith_S(ya, k);
            for (int k = 0; k < x_hi_i; ++k)
                full_conc[P.Nv + k] = NV_Ith_S(ya, x_hi_v + k);
        };

        // Active state vector
        N_Vector y_win = N_VNew_Serial(P_win.N_EQ, sunctx);
        if (!y_win) { std::cerr << "Error allocating window N_Vector\n"; return 1; }
        pack_active(y_win);

        // CVODE and linear-solver handles; recreated on each window expansion.
        void*           cvode_win = nullptr;
        SUNLinearSolver LS_win    = nullptr;
        double          h_last_win = 0.0;   // last step size before reinit

        // Helper: (re)create CVODE + GMRES for current x_hi_v / x_hi_i / y_win.
        int n_reinits = 0;
        long int acc_nsteps_I=0, acc_nfevals_I=0, acc_nniters_I=0, acc_nliters_I=0;
        auto setup_cvode_win = [&](double t0) -> bool {
            if (LS_win)    { SUNLinSolFree(LS_win);  LS_win    = nullptr; }
            if (cvode_win) {
                long int _ns=0,_nf=0,_nn=0,_nl=0;
                CVodeGetNumSteps(cvode_win,&_ns); CVodeGetNumRhsEvals(cvode_win,&_nf);
                CVodeGetNumNonlinSolvIters(cvode_win,&_nn); CVodeGetNumLinIters(cvode_win,&_nl);
                acc_nsteps_I+=_ns; acc_nfevals_I+=_nf; acc_nniters_I+=_nn; acc_nliters_I+=_nl;
                CVodeFree(&cvode_win); cvode_win = nullptr;
            }

            // Sync P_win with current window bounds
            P_win.Nv   = x_hi_v;
            P_win.Ni   = x_hi_i;
            P_win.N_EQ = x_hi_v + x_hi_i;

            cvode_win = CVodeCreate(lmm_win, sunctx);
            if (!cvode_win) { std::cerr << "Window CVodeCreate failed\n"; return false; }

            if (CVodeSetUserData(cvode_win, &P_win) != CV_SUCCESS) return false;
            if (CVodeInit(cvode_win, rhs_cd, t0, y_win) != CV_SUCCESS) return false;
            if (CVodeSStolerances(cvode_win, P.rtol, P.atol) != CV_SUCCESS) return false;
            if (CVodeSetMaxNumSteps(cvode_win, 1000000) != CV_SUCCESS) return false;
            if (P.max_order > 0)
                if (CVodeSetMaxOrd(cvode_win, P.max_order) != CV_SUCCESS) return false;

            // Cap the post-reinit step to the last step of the old system.
            // CVodeSetMaxStep (not InitStep) is used so CVODE can still choose a
            // *smaller* step for the expanded (potentially stiffer) system, while
            // being prevented from jumping to a step that is too large and would
            // cause Newton convergence failures on the very first BDF iteration.
            if (h_last_win > 0.0)
                CVodeSetMaxStep(cvode_win, h_last_win);

            // GMRES: matrix-free, no matrix to resize when window expands.
            // Optional Jacobi diagonal preconditioner (window_prec=1) mirrors
            // Phase II behaviour and prevents GMRES from stalling at large N.
            const int prec_type = (P.window_prec == 1) ? SUN_PREC_LEFT : SUN_PREC_NONE;
            LS_win = SUNLinSol_SPGMR(y_win, prec_type, 0, sunctx);
            if (!LS_win) { std::cerr << "Window SUNLinSol_SPGMR failed\n"; return false; }
            if (CVodeSetLinearSolver(cvode_win, LS_win, nullptr) != CV_SUCCESS) return false;
            if (P.window_prec == 1) {
                if (CVodeSetPreconditioner(cvode_win,
                                           prec_setup_win1,
                                           prec_solve_win1) != CV_SUCCESS) return false;
            }

            ++n_reinits;
            std::cerr << "[window] reinit #" << n_reinits
                      << "  x_hi_v=" << x_hi_v << "  x_hi_i=" << x_hi_i
                      << "  N_active=" << P_win.N_EQ << "  t0=" << t0 << "\n";
            return true;
        };

        if (!setup_cvode_win(t_eval[0])) {
            N_VDestroy(y_win); N_VDestroy(y); SUNContext_Free(&sunctx); return 1;
        }

        std::cerr << "Solver: CVODE " << lmm_win_name
                  << " WINDOW | linsol: GMRES"
                  << " | w0_v=" << P.window_w0_v << " w0_i=" << P.window_w0_i
                  << " | C_expand=" << P.window_C_expand
                  << " | expand_pad=" << P.window_expand_pad
                  << " | Nv=" << P.Nv << " Ni=" << P.Ni << "\n";

        // Output initial point (full Nv+Ni row)
        write_row(t_eval[0], full_conc.data(), N_EQ);

        sunrealtype t_cur_win = t_eval[0];

        for (int i = 1; i < P.n_points; ++i) {
            flag = CVode(cvode_win, t_eval[i], y_win, &t_cur_win, CV_NORMAL);
            if (flag < 0) {
                std::cerr << "Window CVode error at step " << i
                          << " (t=" << t_eval[i] << "): flag=" << flag << "\n";
                N_VDestroy(y_win); if (LS_win) SUNLinSolFree(LS_win);
                if (cvode_win) CVodeFree(&cvode_win);
                N_VDestroy(y); SUNContext_Free(&sunctx);
                if (fp_bin) std::fclose(fp_bin);
                return 1;
            }

            unpack_active(y_win);

            // Write full Nv+Ni row (truncated species remain 0 in full_conc)
            write_row(t_eval[i], full_conc.data(), N_EQ);

            // ── Window expansion check ────────────────────────────────────────
            if (i % P.window_check_every == 0) {
                bool expanded = false;

                // Interstitial upper bound: expand if top cluster exceeds threshold
                if (x_hi_i < P.Ni &&
                    full_conc[P.Nv + x_hi_i - 1] > P.window_C_expand) {
                    int new_hi = (P.window_expand_factor > 1.0)
                        ? std::max(static_cast<int>(x_hi_i * P.window_expand_factor),
                                   x_hi_i + P.window_expand_pad)
                        : x_hi_i + P.window_expand_pad;
                    x_hi_i = std::min(new_hi, P.Ni);
                    expanded = true;
                }
                // Vacancy upper bound
                if (x_hi_v < P.Nv &&
                    full_conc[x_hi_v - 1] > P.window_C_expand) {
                    int new_hi = (P.window_expand_factor > 1.0)
                        ? std::max(static_cast<int>(x_hi_v * P.window_expand_factor),
                                   x_hi_v + P.window_expand_pad)
                        : x_hi_v + P.window_expand_pad;
                    x_hi_v = std::min(new_hi, P.Nv);
                    expanded = true;
                }

                if (expanded) {
                    // Save last step size so reinit can restore continuity.
                    CVodeGetLastStep(cvode_win, &h_last_win);
                    // Resize active vector and reinitialise CVODE
                    N_VDestroy(y_win);
                    y_win = N_VNew_Serial(x_hi_v + x_hi_i, sunctx);
                    if (!y_win) {
                        std::cerr << "Error allocating expanded window N_Vector\n";
                        if (LS_win) SUNLinSolFree(LS_win);
                        if (cvode_win) CVodeFree(&cvode_win);
                        N_VDestroy(y); SUNContext_Free(&sunctx);
                        return 1;
                    }
                    pack_active(y_win);
                    if (!setup_cvode_win(t_cur_win)) {
                        N_VDestroy(y_win); N_VDestroy(y);
                        SUNContext_Free(&sunctx); return 1;
                    }
                }
            }
        }

        { long int _ns=0,_nf=0,_nn=0,_nl=0;
          if (cvode_win) { CVodeGetNumSteps(cvode_win,&_ns); CVodeGetNumRhsEvals(cvode_win,&_nf);
            CVodeGetNumNonlinSolvIters(cvode_win,&_nn); CVodeGetNumLinIters(cvode_win,&_nl); }
          acc_nsteps_I+=_ns; acc_nfevals_I+=_nf; acc_nniters_I+=_nn; acc_nliters_I+=_nl; }
        std::cerr << "  Phase I: " << n_reinits << " CVODE reinits"
                  << "  final window: Cv[1.." << x_hi_v << "]  Ci[1.." << x_hi_i << "]\n"
                  << "  Stats: steps=" << acc_nsteps_I << " rhs_evals=" << acc_nfevals_I
                  << " newton_iters=" << acc_nniters_I << " lin_iters=" << acc_nliters_I
                  << "  [" << (acc_nniters_I>0?(double)acc_nliters_I/acc_nniters_I:0.0)
                  << " GMRES iters/Newton]\n";
        if (LS_win)    SUNLinSolFree(LS_win);
        if (cvode_win) CVodeFree(&cvode_win);
        N_VDestroy(y_win);

    } else if (P.window_mode == 2) {
        // ── CVODE — Phase II sliding-window solver (upper + lower truncation) ──
        //
        // Active window: [Cv1..Cv_{x_hi_v}] + [Ci1] + [Ci_{x_lo_i}..Ci_{x_hi_i}]
        //   x_hi grows (geometric) when the top cluster exceeds window_C_expand.
        //   x_lo_i advances when the lowest active cluster reaches QSS:
        //     |dCi_{x_lo}/dt| / Ci_{x_lo} < window_C_contract
        // Optional Jacobi diagonal preconditioner (window_prec=1).

        const int lmm2 = (P.lmm == 1) ? CV_ADAMS : CV_BDF;
        const char* lmm2_name = (P.lmm == 1) ? "ADAMS" : "BDF";

        // Initial window bounds
        int x_hi_v2 = std::min(std::max(P.window_w0_v, 2), P.Nv);
        int x_lo_i2 = 2;   // always start at Ci2; Ci1 is always a separate slot
        int x_hi_i2 = std::min(std::max(P.window_w0_i, 2), P.Ni);

        // Build WindowData
        WindowData W;
        W.P_full         = &P;
        W.x_hi_v         = x_hi_v2;
        W.x_lo_i         = x_lo_i2;
        W.x_hi_i         = x_hi_i2;
        W.N_active       = x_hi_v2 + 1 + (x_hi_i2 - x_lo_i2 + 1);
        W.frozen_KLI_sum = 0.0;
        W.frozen_KLV_sum = 0.0;
        W.Ci_frozen_top  = 0.0;
        W.full_conc.assign(P.N_EQ, 0.0);
        for (int k = 0; k < P.N_EQ; ++k) W.full_conc[k] = P.y0[k];

        // Pack/unpack between full_conc and the active state vector.
        // Layout: y[0..x_hi_v-1]=Cv1..Cv_{x_hi_v},  y[x_hi_v]=Ci1,
        //         y[x_hi_v+1+j]=Ci_{x_lo_i+j}  for j=0..n-1 (n=x_hi_i-x_lo_i+1)
        auto pack2 = [&](N_Vector ya) {
            for (int k = 0; k < W.x_hi_v; ++k)
                NV_Ith_S(ya, k) = W.full_conc[k];
            NV_Ith_S(ya, W.x_hi_v) = W.full_conc[P.Nv];   // Ci1
            const int n = W.x_hi_i - W.x_lo_i + 1;
            for (int j = 0; j < n; ++j)
                NV_Ith_S(ya, W.x_hi_v + 1 + j) =
                    W.full_conc[P.Nv + W.x_lo_i - 1 + j];  // Ci_{x_lo_i+j}
        };
        auto unpack2 = [&](N_Vector ya) {
            for (int k = 0; k < W.x_hi_v; ++k)
                W.full_conc[k] = NV_Ith_S(ya, k);
            W.full_conc[P.Nv] = NV_Ith_S(ya, W.x_hi_v);   // Ci1
            const int n = W.x_hi_i - W.x_lo_i + 1;
            for (int j = 0; j < n; ++j)
                W.full_conc[P.Nv + W.x_lo_i - 1 + j] =
                    NV_Ith_S(ya, W.x_hi_v + 1 + j);
        };

        N_Vector y_win2 = N_VNew_Serial(W.N_active, sunctx);
        if (!y_win2) { std::cerr << "Error allocating Phase II N_Vector\n"; return 1; }
        pack2(y_win2);

        void*           cvode2 = nullptr;
        SUNLinearSolver LS2    = nullptr;
        int n_reinits2         = 0;
        double h_last2         = 0.0;   // last step size before reinit
        long int acc_nsteps_II=0, acc_nfevals_II=0, acc_nniters_II=0, acc_nliters_II=0;

        // (Re)create CVODE + GMRES for the current window size.
        auto setup2 = [&](double t0) -> bool {
            if (LS2)    { SUNLinSolFree(LS2); LS2    = nullptr; }
            if (cvode2) {
                long int _ns=0,_nf=0,_nn=0,_nl=0;
                CVodeGetNumSteps(cvode2,&_ns); CVodeGetNumRhsEvals(cvode2,&_nf);
                CVodeGetNumNonlinSolvIters(cvode2,&_nn); CVodeGetNumLinIters(cvode2,&_nl);
                acc_nsteps_II+=_ns; acc_nfevals_II+=_nf; acc_nniters_II+=_nn; acc_nliters_II+=_nl;
                CVodeFree(&cvode2); cvode2 = nullptr;
            }

            W.x_hi_v  = x_hi_v2;
            W.x_lo_i  = x_lo_i2;
            W.x_hi_i  = x_hi_i2;
            W.N_active = x_hi_v2 + 1 + (x_hi_i2 - x_lo_i2 + 1);

            cvode2 = CVodeCreate(lmm2, sunctx);
            if (!cvode2) return false;
            if (CVodeSetUserData(cvode2, &W)                   != CV_SUCCESS) return false;
            if (CVodeInit(cvode2, rhs_window, t0, y_win2)      != CV_SUCCESS) return false;
            if (CVodeSStolerances(cvode2, P.rtol, P.atol)      != CV_SUCCESS) return false;
            if (CVodeSetMaxNumSteps(cvode2, 1000000)            != CV_SUCCESS) return false;
            if (P.max_order > 0)
                CVodeSetMaxOrd(cvode2, P.max_order);

            // Cap the post-reinit step to the last step of the old system.
            // CVodeSetMaxStep (not InitStep) is used so CVODE can still choose a
            // *smaller* step for the expanded (potentially stiffer) system, while
            // being prevented from jumping to a step that is too large and would
            // cause Newton convergence failures on the very first BDF iteration.
            if (h_last2 > 0.0)
                CVodeSetMaxStep(cvode2, h_last2);

            const int prec_type = (P.window_prec > 0) ? SUN_PREC_LEFT : SUN_PREC_NONE;
            LS2 = SUNLinSol_SPGMR(y_win2, prec_type, 0, sunctx);
            if (!LS2) return false;
            if (CVodeSetLinearSolver(cvode2, LS2, nullptr) != CV_SUCCESS) return false;
            if (P.window_prec > 0)
                if (CVodeSetPreconditioner(cvode2,
                        prec_setup_window, prec_solve_window) != CV_SUCCESS) return false;
            ++n_reinits2;
            return true;
        };

        if (!setup2(t_eval[0])) {
            N_VDestroy(y_win2); N_VDestroy(y); SUNContext_Free(&sunctx); return 1;
        }

        std::cerr << "Solver: CVODE " << lmm2_name << " WINDOW-II (sliding)"
                  << " | linsol: GMRES" << (P.window_prec ? " + Jacobi" : "")
                  << " | w0_v=" << P.window_w0_v << " w0_i=" << P.window_w0_i
                  << " | C_expand=" << P.window_C_expand
                  << " | expand_factor=" << P.window_expand_factor
                  << " | C_contract=" << P.window_C_contract
                  << " | Nv=" << P.Nv << " Ni=" << P.Ni << "\n";

        // Output initial point (full Nv+Ni row)
        write_row(t_eval[0], W.full_conc.data(), N_EQ);

        sunrealtype t_cur2 = t_eval[0];

        for (int i = 1; i < P.n_points; ++i) {
            flag = CVode(cvode2, t_eval[i], y_win2, &t_cur2, CV_NORMAL);
            if (flag < 0) {
                std::cerr << "Phase II CVode error at step " << i
                          << " (t=" << t_eval[i] << "): flag=" << flag << "\n";
                N_VDestroy(y_win2); if (LS2) SUNLinSolFree(LS2);
                if (cvode2) CVodeFree(&cvode2);
                N_VDestroy(y); SUNContext_Free(&sunctx);
                if (fp_bin) std::fclose(fp_bin);
                return 1;
            }

            unpack2(y_win2);

            // Write full Nv+Ni output row
            write_row(t_eval[i], W.full_conc.data(), N_EQ);

            if (i % P.window_check_every != 0) continue;

            bool changed = false;

            // ── Upper expansion (geometric or additive) ───────────────────────
            if (x_hi_i2 < P.Ni &&
                W.full_conc[P.Nv + x_hi_i2 - 1] > P.window_C_expand) {
                int new_hi = (P.window_expand_factor > 1.0)
                    ? std::max(static_cast<int>(x_hi_i2 * P.window_expand_factor),
                               x_hi_i2 + P.window_expand_pad)
                    : x_hi_i2 + P.window_expand_pad;
                x_hi_i2 = std::min(new_hi, P.Ni);
                changed = true;
            }
            if (x_hi_v2 < P.Nv &&
                W.full_conc[x_hi_v2 - 1] > P.window_C_expand) {
                int new_hi = (P.window_expand_factor > 1.0)
                    ? std::max(static_cast<int>(x_hi_v2 * P.window_expand_factor),
                               x_hi_v2 + P.window_expand_pad)
                    : x_hi_v2 + P.window_expand_pad;
                x_hi_v2 = std::min(new_hi, P.Nv);
                changed = true;
            }

            // ── Lower contraction (advance x_lo_i) ────────────────────────────
            // Criterion: |dCi_{x_lo}/dt| / Ci_{x_lo} < window_C_contract
            if (P.window_C_contract > 0.0 &&
                (x_hi_i2 - x_lo_i2 + 1) > P.window_min_active_i) {

                N_Vector ydot_tmp = N_VNew_Serial(W.N_active, sunctx);
                if (ydot_tmp) {
                    rhs_window(t_cur2, y_win2, ydot_tmp, &W);
                    const double Ci_lo  = NV_Ith_S(y_win2,   x_hi_v2 + 1);
                    const double dCi_lo = NV_Ith_S(ydot_tmp, x_hi_v2 + 1);
                    N_VDestroy(ydot_tmp);

                    const double rel = std::abs(dCi_lo)
                                       / std::max(Ci_lo, P.C_floor);

                    if (rel < P.window_C_contract) {
                        bool do_contract = true;

                        // Optional nucleation guard: block freezing Ci2 while
                        // Ci1+Ci1->Ci2 nucleation dominates its outflow.
                        // Disabled when window_nuc_guard == 0.0 (default).
                        if (do_contract && P.window_nuc_guard > 0.0 && x_lo_i2 == 2) {
                            const double Ci1_g = NV_Ith_S(y_win2, x_hi_v2);
                            const double Cv1_g = NV_Ith_S(y_win2, 0);
                            const double nuc_g = 0.5 * P.K_nuc_i * Ci1_g * Ci1_g;
                            const double out_g = (P.KLI[1]*Ci1_g + P.KLV[1]*Cv1_g)
                                                 * std::max(Ci_lo, P.C_floor);
                            if (nuc_g > P.window_nuc_guard * out_g)
                                do_contract = false;
                        }

                        if (do_contract) {
                            // Extend frozen sums (Ci2 always explicit; only Ci3+)
                            if (x_lo_i2 >= 3) {
                                W.frozen_KLI_sum += P.KLI[x_lo_i2 - 1] * Ci_lo;
                                W.frozen_KLV_sum += P.KLV[x_lo_i2 - 1] * Ci_lo;
                            }
                            W.Ci_frozen_top = Ci_lo;   // = Ci_{old x_lo_i}
                            x_lo_i2 += 1;
                            changed = true;
                        }
                    }
                }
            }

            if (changed) {
                // Save last step size so reinit can restore continuity.
                CVodeGetLastStep(cvode2, &h_last2);
                // Sync W fields from local vars BEFORE pack2 and N_VNew_Serial.
                // pack2 reads W.x_hi_v / W.x_lo_i / W.x_hi_i to determine the
                // layout; if they still hold the old values after a contraction
                // step (x_lo_i2 advanced by 1), pack2 would overflow the newly
                // allocated (smaller) vector by exactly one element.
                W.x_hi_v   = x_hi_v2;
                W.x_lo_i   = x_lo_i2;
                W.x_hi_i   = x_hi_i2;
                W.N_active = x_hi_v2 + 1 + (x_hi_i2 - x_lo_i2 + 1);
                N_VDestroy(y_win2);
                y_win2 = N_VNew_Serial(W.N_active, sunctx);
                if (!y_win2) {
                    std::cerr << "Error allocating Phase II expanded N_Vector\n";
                    if (LS2) SUNLinSolFree(LS2); if (cvode2) CVodeFree(&cvode2);
                    N_VDestroy(y); SUNContext_Free(&sunctx); return 1;
                }
                pack2(y_win2);
                if (!setup2(t_cur2)) {
                    N_VDestroy(y_win2); N_VDestroy(y);
                    SUNContext_Free(&sunctx); return 1;
                }
            }
        }

        { long int _ns=0,_nf=0,_nn=0,_nl=0;
          if (cvode2) { CVodeGetNumSteps(cvode2,&_ns); CVodeGetNumRhsEvals(cvode2,&_nf);
            CVodeGetNumNonlinSolvIters(cvode2,&_nn); CVodeGetNumLinIters(cvode2,&_nl); }
          acc_nsteps_II+=_ns; acc_nfevals_II+=_nf; acc_nniters_II+=_nn; acc_nliters_II+=_nl; }
        std::cerr << "  Phase II: " << n_reinits2 << " CVODE reinits"
                  << "  final window: Cv[1.." << x_hi_v2 << "]"
                  << "  Ci[" << x_lo_i2 << ".." << x_hi_i2 << "]\n"
                  << "  Stats: steps=" << acc_nsteps_II << " rhs_evals=" << acc_nfevals_II
                  << " newton_iters=" << acc_nniters_II << " lin_iters=" << acc_nliters_II
                  << "  [" << (acc_nniters_II>0?(double)acc_nliters_II/acc_nniters_II:0.0)
                  << " GMRES iters/Newton]\n";

        if (LS2)    SUNLinSolFree(LS2);
        if (cvode2) CVodeFree(&cvode2);
        N_VDestroy(y_win2);

    } else if (P.window_mode == 3) {
        // ── CVODE — Phase III constant-width sliding window ───────────────────
        //
        // Window of fixed width W = window_width slides upward with the cluster
        // front:  x_lo_i = max(2, x_hi_i - W + 1).
        // Upper bound expands the same way as Phase I/II (threshold-triggered).
        // Lower sliding is gated behind t > window_t_start to skip the nucleation
        // transient (during which Ci2 must remain active).
        // Reuses WindowData / rhs_window / recompute_frozen_sums from Phase II.

        const int lmm3 = (P.lmm == 1) ? CV_ADAMS : CV_BDF;
        const char* lmm3_name = (P.lmm == 1) ? "ADAMS" : "BDF";
        const int   win_width = std::max(P.window_width, 2);

        // Initial window bounds (lower stays at 2 until t_start)
        int x_hi_v3 = std::min(std::max(P.window_w0_v, 2), P.Nv);
        int x_hi_i3 = std::min(std::max(P.window_w0_i, 2), P.Ni);
        int x_lo_i3 = 2;

        WindowData W3d;
        W3d.P_full         = &P;
        W3d.x_hi_v         = x_hi_v3;
        W3d.x_lo_i         = x_lo_i3;
        W3d.x_hi_i         = x_hi_i3;
        W3d.N_active       = x_hi_v3 + 1 + (x_hi_i3 - x_lo_i3 + 1);
        W3d.frozen_KLI_sum = 0.0;
        W3d.frozen_KLV_sum = 0.0;
        W3d.Ci_frozen_top  = 0.0;
        W3d.full_conc.assign(P.N_EQ, 0.0);
        for (int k = 0; k < P.N_EQ; ++k) W3d.full_conc[k] = P.y0[k];

        // pack/unpack — identical layout to Phase II
        auto pack3 = [&](N_Vector ya) {
            for (int k = 0; k < W3d.x_hi_v; ++k)
                NV_Ith_S(ya, k) = W3d.full_conc[k];
            NV_Ith_S(ya, W3d.x_hi_v) = W3d.full_conc[P.Nv];  // Ci1
            const int n = W3d.x_hi_i - W3d.x_lo_i + 1;
            for (int j = 0; j < n; ++j)
                NV_Ith_S(ya, W3d.x_hi_v + 1 + j) =
                    W3d.full_conc[P.Nv + W3d.x_lo_i - 1 + j];
        };
        auto unpack3 = [&](N_Vector ya) {
            for (int k = 0; k < W3d.x_hi_v; ++k)
                W3d.full_conc[k] = NV_Ith_S(ya, k);
            W3d.full_conc[P.Nv] = NV_Ith_S(ya, W3d.x_hi_v);  // Ci1
            const int n = W3d.x_hi_i - W3d.x_lo_i + 1;
            for (int j = 0; j < n; ++j)
                W3d.full_conc[P.Nv + W3d.x_lo_i - 1 + j] =
                    NV_Ith_S(ya, W3d.x_hi_v + 1 + j);
        };

        N_Vector y_win3 = N_VNew_Serial(W3d.N_active, sunctx);
        if (!y_win3) { std::cerr << "Error allocating Phase III N_Vector\n"; return 1; }
        pack3(y_win3);

        void*           cvode3   = nullptr;
        SUNLinearSolver LS3      = nullptr;
        int             n_reinits3 = 0;
        double          h_last3  = 0.0;
        long int acc_nsteps_III=0, acc_nfevals_III=0, acc_nniters_III=0, acc_nliters_III=0;

        auto setup3 = [&](double t0) -> bool {
            if (LS3)    { SUNLinSolFree(LS3);  LS3    = nullptr; }
            if (cvode3) {
                long int _ns=0,_nf=0,_nn=0,_nl=0;
                CVodeGetNumSteps(cvode3,&_ns); CVodeGetNumRhsEvals(cvode3,&_nf);
                CVodeGetNumNonlinSolvIters(cvode3,&_nn); CVodeGetNumLinIters(cvode3,&_nl);
                acc_nsteps_III+=_ns; acc_nfevals_III+=_nf; acc_nniters_III+=_nn; acc_nliters_III+=_nl;
                CVodeFree(&cvode3); cvode3 = nullptr;
            }

            W3d.x_hi_v   = x_hi_v3;
            W3d.x_lo_i   = x_lo_i3;
            W3d.x_hi_i   = x_hi_i3;
            W3d.N_active = x_hi_v3 + 1 + (x_hi_i3 - x_lo_i3 + 1);

            cvode3 = CVodeCreate(lmm3, sunctx);
            if (!cvode3) return false;
            if (CVodeSetUserData(cvode3, &W3d)                  != CV_SUCCESS) return false;
            if (CVodeInit(cvode3, rhs_window, t0, y_win3)       != CV_SUCCESS) return false;
            if (CVodeSStolerances(cvode3, P.rtol, P.atol)       != CV_SUCCESS) return false;
            if (CVodeSetMaxNumSteps(cvode3, 1000000)             != CV_SUCCESS) return false;
            if (P.max_order > 0)
                CVodeSetMaxOrd(cvode3, P.max_order);
            if (h_last3 > 0.0)
                CVodeSetMaxStep(cvode3, h_last3);

            const int prec_type = (P.window_prec > 0) ? SUN_PREC_LEFT : SUN_PREC_NONE;
            LS3 = SUNLinSol_SPGMR(y_win3, prec_type, 0, sunctx);
            if (!LS3) return false;
            if (CVodeSetLinearSolver(cvode3, LS3, nullptr) != CV_SUCCESS) return false;
            if (P.window_prec > 0)
                if (CVodeSetPreconditioner(cvode3,
                        prec_setup_window, prec_solve_window) != CV_SUCCESS) return false;
            ++n_reinits3;
            return true;
        };

        if (!setup3(t_eval[0])) {
            N_VDestroy(y_win3); N_VDestroy(y); SUNContext_Free(&sunctx); return 1;
        }

        std::cerr << "Solver: CVODE " << lmm3_name << " WINDOW-III (const-width)"
                  << " | linsol: GMRES" << (P.window_prec ? " + Jacobi" : "")
                  << " | w0_v=" << P.window_w0_v << " w0_i=" << P.window_w0_i
                  << " | W=" << win_width
                  << " | C_expand=" << P.window_C_expand
                  << " | t_start=" << P.window_t_start
                  << " | Nv=" << P.Nv << " Ni=" << P.Ni << "\n";

        // Output initial point (full Nv+Ni row)
        write_row(t_eval[0], W3d.full_conc.data(), N_EQ);

        sunrealtype t_cur3 = t_eval[0];

        for (int i = 1; i < P.n_points; ++i) {
            flag = CVode(cvode3, t_eval[i], y_win3, &t_cur3, CV_NORMAL);
            if (flag < 0) {
                std::cerr << "Phase III CVode error at step " << i
                          << " (t=" << t_eval[i] << "): flag=" << flag << "\n";
                N_VDestroy(y_win3); if (LS3) SUNLinSolFree(LS3);
                if (cvode3) CVodeFree(&cvode3);
                N_VDestroy(y); SUNContext_Free(&sunctx);
                if (fp_bin) std::fclose(fp_bin);
                return 1;
            }

            unpack3(y_win3);

            // Write full Nv+Ni output row
            write_row(t_eval[i], W3d.full_conc.data(), N_EQ);

            if (i % P.window_check_every != 0) continue;

            // ── Ni domain extension ───────────────────────────────────────────
            if (P.Ni < P.Ni_max) {
                const bool near_bdy = (P.Ni_extend_margin > 0 &&
                                       x_hi_i3 >= P.Ni - P.Ni_extend_margin);
                const double Ci_top = W3d.full_conc[P.Nv + P.Ni - 1];
                const double Ci1val = W3d.full_conc[P.Nv];
                const double F_bdy  = P.KLI[P.Ni - 1] * Ci1val * Ci_top;
                const bool consv    = (P.Ni_extend_tol > 0.0 && P.P_prod > 0.0 &&
                                       F_bdy / P.P_prod > P.Ni_extend_tol);
                if (near_bdy || consv) {
                    P.Ni = std::min(P.Ni + P.window_expand_pad, P.Ni_max);
                    std::cerr << "  [Ni -> " << P.Ni << "  t=" << t_cur3 << "]\n";
                }
            }

            bool changed = false;

            // Save upper bound before any expansion (needed by lower-bound cap)
            const int x_hi_i3_pre = x_hi_i3;

            // ── Upper expansion (same as Phase I/II) ─────────────────────────
            if (x_hi_i3 < P.Ni &&
                W3d.full_conc[P.Nv + x_hi_i3 - 1] > P.window_C_expand) {
                int new_hi = (P.window_expand_factor > 1.0)
                    ? std::max(static_cast<int>(x_hi_i3 * P.window_expand_factor),
                               x_hi_i3 + P.window_expand_pad)
                    : x_hi_i3 + P.window_expand_pad;
                x_hi_i3 = std::min(new_hi, P.Ni);
                changed  = true;
            }
            if (x_hi_v3 < P.Nv &&
                W3d.full_conc[x_hi_v3 - 1] > P.window_C_expand) {
                int new_hi = (P.window_expand_factor > 1.0)
                    ? std::max(static_cast<int>(x_hi_v3 * P.window_expand_factor),
                               x_hi_v3 + P.window_expand_pad)
                    : x_hi_v3 + P.window_expand_pad;
                x_hi_v3 = std::min(new_hi, P.Nv);
                changed  = true;
            }

            // ── Lower bound: constant-width coupling ─────────────────────────
            // Slide only after t_start; x_lo_i = max(2, x_hi_i - W + 1).
            // Safety cap: never skip unactivated clusters — new_lo must not
            // exceed (x_hi_i_before_expand + 1) to avoid a frozen gap.
            if (t_cur3 > P.window_t_start) {
                const int new_lo = std::max(2, x_hi_i3 - win_width + 1);
                // x_hi_i3_pre holds the upper bound before this expansion step;
                // it was saved to x_hi_i3 at the top of this check block.
                const int safe_lo = std::min(new_lo, x_hi_i3_pre + 1);
                if (safe_lo > x_lo_i3) {
                    x_lo_i3      = safe_lo;
                    W3d.x_lo_i   = x_lo_i3;
                    recompute_frozen_sums(W3d);  // recompute from full_conc
                    changed      = true;
                }
            }

            if (changed) {
                CVodeGetLastStep(cvode3, &h_last3);
                W3d.x_hi_v   = x_hi_v3;
                W3d.x_lo_i   = x_lo_i3;
                W3d.x_hi_i   = x_hi_i3;
                W3d.N_active = x_hi_v3 + 1 + (x_hi_i3 - x_lo_i3 + 1);
                N_VDestroy(y_win3);
                y_win3 = N_VNew_Serial(W3d.N_active, sunctx);
                if (!y_win3) {
                    std::cerr << "Error allocating Phase III N_Vector\n";
                    if (LS3) SUNLinSolFree(LS3); if (cvode3) CVodeFree(&cvode3);
                    N_VDestroy(y); SUNContext_Free(&sunctx); return 1;
                }
                pack3(y_win3);
                if (!setup3(t_cur3)) {
                    N_VDestroy(y_win3); N_VDestroy(y);
                    SUNContext_Free(&sunctx); return 1;
                }
            }
        }

        { long int _ns=0,_nf=0,_nn=0,_nl=0;
          if (cvode3) { CVodeGetNumSteps(cvode3,&_ns); CVodeGetNumRhsEvals(cvode3,&_nf);
            CVodeGetNumNonlinSolvIters(cvode3,&_nn); CVodeGetNumLinIters(cvode3,&_nl); }
          acc_nsteps_III+=_ns; acc_nfevals_III+=_nf; acc_nniters_III+=_nn; acc_nliters_III+=_nl; }
        std::cerr << "  Phase III: " << n_reinits3 << " CVODE reinits"
                  << "  final window: Cv[1.." << x_hi_v3 << "]"
                  << "  Ci[" << x_lo_i3 << ".." << x_hi_i3 << "]\n"
                  << "  Stats: steps=" << acc_nsteps_III << " rhs_evals=" << acc_nfevals_III
                  << " newton_iters=" << acc_nniters_III << " lin_iters=" << acc_nliters_III
                  << "  [" << (acc_nniters_III>0?(double)acc_nliters_III/acc_nniters_III:0.0)
                  << " GMRES iters/Newton]\n";

        if (LS3)    SUNLinSolFree(LS3);
        if (cvode3) CVodeFree(&cvode3);
        N_VDestroy(y_win3);

    } else if (P.window_mode == 4) {
        // ── CVODE — Phase IV: Multithread-OpenMP constant-width sliding window ──
        //
        // Algorithm: identical to Phase III.
        // Differences from Phase III:
        //   • rhs_window_omp()  — hot loops parallelised with OpenMP
        //   • WindowDataOMP     — pre-allocated Cv_buf/Ci_buf (no malloc per RHS)
        //   • window_omp_threads controls thread count (0 = OMP_NUM_THREADS)

        const int lmm4 = (P.lmm == 1) ? CV_ADAMS : CV_BDF;
        const char* lmm4_name = (P.lmm == 1) ? "ADAMS" : "BDF";
        const int   win_width4 = std::max(P.window_width, 2);

        // Set thread count before any parallel region
#ifdef CD_HAVE_OPENMP
        if (P.window_omp_threads > 0)
            omp_set_num_threads(P.window_omp_threads);
        const int n_thr_actual = (P.window_omp_threads > 0)
                                     ? P.window_omp_threads
                                     : omp_get_max_threads();
#else
        const int n_thr_actual = 1;
#endif

        int x_hi_v4 = std::min(std::max(P.window_w0_v, 2), P.Nv);
        int x_hi_i4 = std::min(std::max(P.window_w0_i, 2), P.Ni);
        int x_lo_i4 = 2;

        WindowDataOMP W4;
        W4.P_full         = &P;
        W4.x_hi_v         = x_hi_v4;
        W4.x_lo_i         = x_lo_i4;
        W4.x_hi_i         = x_hi_i4;
        W4.N_active       = x_hi_v4 + 1 + (x_hi_i4 - x_lo_i4 + 1);
        W4.frozen_KLI_sum = 0.0;
        W4.frozen_KLV_sum = 0.0;
        W4.Ci_frozen_top  = 0.0;
        W4.n_omp_threads  = P.window_omp_threads;
        W4.full_conc.assign(P.N_EQ, 0.0);
        for (int k = 0; k < P.N_EQ; ++k) W4.full_conc[k] = P.y0[k];
        W4.resize_buffers();

        // pack/unpack — identical layout to Phase II/III
        auto pack4 = [&](N_Vector ya) {
            for (int k = 0; k < W4.x_hi_v; ++k)
                NV_Ith_S(ya, k) = W4.full_conc[k];
            NV_Ith_S(ya, W4.x_hi_v) = W4.full_conc[P.Nv];  // Ci1
            const int n = W4.x_hi_i - W4.x_lo_i + 1;
            for (int j = 0; j < n; ++j)
                NV_Ith_S(ya, W4.x_hi_v + 1 + j) =
                    W4.full_conc[P.Nv + W4.x_lo_i - 1 + j];
        };
        auto unpack4 = [&](N_Vector ya) {
            for (int k = 0; k < W4.x_hi_v; ++k)
                W4.full_conc[k] = NV_Ith_S(ya, k);
            W4.full_conc[P.Nv] = NV_Ith_S(ya, W4.x_hi_v);  // Ci1
            const int n = W4.x_hi_i - W4.x_lo_i + 1;
            for (int j = 0; j < n; ++j)
                W4.full_conc[P.Nv + W4.x_lo_i - 1 + j] =
                    NV_Ith_S(ya, W4.x_hi_v + 1 + j);
        };

        N_Vector y_win4 = N_VNew_Serial(W4.N_active, sunctx);
        if (!y_win4) { std::cerr << "Error allocating Phase IV N_Vector\n"; return 1; }
        pack4(y_win4);

        void*           cvode4   = nullptr;
        SUNLinearSolver LS4      = nullptr;
        int             n_reinits4 = 0;
        double          h_last4  = 0.0;
        long int acc_nsteps_IV=0, acc_nfevals_IV=0, acc_nniters_IV=0, acc_nliters_IV=0;

        auto setup4 = [&](double t0) -> bool {
            if (LS4)    { SUNLinSolFree(LS4);  LS4    = nullptr; }
            if (cvode4) {
                long int _ns=0,_nf=0,_nn=0,_nl=0;
                CVodeGetNumSteps(cvode4,&_ns); CVodeGetNumRhsEvals(cvode4,&_nf);
                CVodeGetNumNonlinSolvIters(cvode4,&_nn); CVodeGetNumLinIters(cvode4,&_nl);
                acc_nsteps_IV+=_ns; acc_nfevals_IV+=_nf; acc_nniters_IV+=_nn; acc_nliters_IV+=_nl;
                CVodeFree(&cvode4); cvode4 = nullptr;
            }

            W4.x_hi_v   = x_hi_v4;
            W4.x_lo_i   = x_lo_i4;
            W4.x_hi_i   = x_hi_i4;
            W4.N_active = x_hi_v4 + 1 + (x_hi_i4 - x_lo_i4 + 1);
            W4.resize_buffers();   // keep scratch buffers in sync

            cvode4 = CVodeCreate(lmm4, sunctx);
            if (!cvode4) return false;
            if (CVodeSetUserData(cvode4, &W4)                    != CV_SUCCESS) return false;
            if (CVodeInit(cvode4, rhs_window_omp, t0, y_win4)   != CV_SUCCESS) return false;
            if (CVodeSStolerances(cvode4, P.rtol, P.atol)        != CV_SUCCESS) return false;
            if (CVodeSetMaxNumSteps(cvode4, 1000000)              != CV_SUCCESS) return false;
            if (P.max_order > 0)
                CVodeSetMaxOrd(cvode4, P.max_order);
            if (h_last4 > 0.0)
                CVodeSetMaxStep(cvode4, h_last4);

            const int prec_type = (P.window_prec > 0) ? SUN_PREC_LEFT : SUN_PREC_NONE;
            LS4 = SUNLinSol_SPGMR(y_win4, prec_type, 0, sunctx);
            if (!LS4) return false;
            if (CVodeSetLinearSolver(cvode4, LS4, nullptr) != CV_SUCCESS) return false;
            if (P.window_prec > 0)
                if (CVodeSetPreconditioner(cvode4,
                        prec_setup_window_omp, prec_solve_window_omp) != CV_SUCCESS)
                    return false;
            ++n_reinits4;
            return true;
        };

        if (!setup4(t_eval[0])) {
            N_VDestroy(y_win4); N_VDestroy(y); SUNContext_Free(&sunctx); return 1;
        }

        std::cerr << "Solver: CVODE " << lmm4_name << " WINDOW-IV (OpenMP)"
                  << " | linsol: GMRES" << (P.window_prec ? " + Jacobi" : "")
                  << " | threads=" << n_thr_actual
                  << " | w0_v=" << P.window_w0_v << " w0_i=" << P.window_w0_i
                  << " | W=" << win_width4
                  << " | C_expand=" << P.window_C_expand
                  << " | t_start=" << P.window_t_start
                  << " | Nv=" << P.Nv << " Ni=" << P.Ni << "\n";

        // Output initial point (full Nv+Ni row)
        write_row(t_eval[0], W4.full_conc.data(), N_EQ);

        sunrealtype t_cur4 = t_eval[0];

        for (int i = 1; i < P.n_points; ++i) {
            flag = CVode(cvode4, t_eval[i], y_win4, &t_cur4, CV_NORMAL);
            if (flag < 0) {
                std::cerr << "Phase IV CVode error at step " << i
                          << " (t=" << t_eval[i] << "): flag=" << flag << "\n";
                N_VDestroy(y_win4); if (LS4) SUNLinSolFree(LS4);
                if (cvode4) CVodeFree(&cvode4);
                N_VDestroy(y); SUNContext_Free(&sunctx);
                if (fp_bin) std::fclose(fp_bin);
                return 1;
            }

            unpack4(y_win4);

            // Write full Nv+Ni output row
            write_row(t_eval[i], W4.full_conc.data(), N_EQ);

            // ── Progress report to stderr every 10 output points ─────────────
            if (i % 10 == 0) {
                std::cerr << "  [pt " << std::setw(4) << i << "/" << (P.n_points - 1)
                          << "  t=" << std::scientific << std::setprecision(3) << t_cur4
                          << "  N_active=" << std::setw(6) << W4.N_active
                          << "  Cv[1.." << x_hi_v4 << "]"
                          << "  Ci[" << x_lo_i4 << ".." << x_hi_i4 << "]]\n";
                std::cerr.flush();
            }

            if (i % P.window_check_every != 0) continue;

            // ── Ni domain extension ───────────────────────────────────────────
            if (P.Ni < P.Ni_max) {
                const bool near_bdy = (P.Ni_extend_margin > 0 &&
                                       x_hi_i4 >= P.Ni - P.Ni_extend_margin);
                const double Ci_top = W4.full_conc[P.Nv + P.Ni - 1];
                const double Ci1val = W4.full_conc[P.Nv];
                const double F_bdy  = P.KLI[P.Ni - 1] * Ci1val * Ci_top;
                const bool consv    = (P.Ni_extend_tol > 0.0 && P.P_prod > 0.0 &&
                                       F_bdy / P.P_prod > P.Ni_extend_tol);
                if (near_bdy || consv) {
                    P.Ni = std::min(P.Ni + P.window_expand_pad, P.Ni_max);
                    std::cerr << "  [Ni -> " << P.Ni << "  t=" << t_cur4 << "]\n";
                }
            }

            bool changed4 = false;
            const int x_hi_i4_pre = x_hi_i4;

            // ── Upper expansion (same as Phase III) ───────────────────────────
            if (x_hi_i4 < P.Ni &&
                W4.full_conc[P.Nv + x_hi_i4 - 1] > P.window_C_expand) {
                int new_hi = (P.window_expand_factor > 1.0)
                    ? std::max(static_cast<int>(x_hi_i4 * P.window_expand_factor),
                               x_hi_i4 + P.window_expand_pad)
                    : x_hi_i4 + P.window_expand_pad;
                x_hi_i4  = std::min(new_hi, P.Ni);
                changed4 = true;
            }
            if (x_hi_v4 < P.Nv &&
                W4.full_conc[x_hi_v4 - 1] > P.window_C_expand) {
                int new_hi = (P.window_expand_factor > 1.0)
                    ? std::max(static_cast<int>(x_hi_v4 * P.window_expand_factor),
                               x_hi_v4 + P.window_expand_pad)
                    : x_hi_v4 + P.window_expand_pad;
                x_hi_v4  = std::min(new_hi, P.Nv);
                changed4 = true;
            }

            // ── Lower bound: constant-width coupling ──────────────────────────
            if (t_cur4 > P.window_t_start) {
                const int new_lo  = std::max(2, x_hi_i4 - win_width4 + 1);
                const int safe_lo = std::min(new_lo, x_hi_i4_pre + 1);
                if (safe_lo > x_lo_i4) {
                    x_lo_i4    = safe_lo;
                    W4.x_lo_i  = x_lo_i4;
                    recompute_frozen_sums_omp(W4);
                    changed4   = true;
                }
            }

            if (changed4) {
                CVodeGetLastStep(cvode4, &h_last4);
                W4.x_hi_v   = x_hi_v4;
                W4.x_lo_i   = x_lo_i4;
                W4.x_hi_i   = x_hi_i4;
                W4.N_active = x_hi_v4 + 1 + (x_hi_i4 - x_lo_i4 + 1);
                N_VDestroy(y_win4);
                y_win4 = N_VNew_Serial(W4.N_active, sunctx);
                if (!y_win4) {
                    std::cerr << "Error allocating Phase IV N_Vector\n";
                    if (LS4) SUNLinSolFree(LS4); if (cvode4) CVodeFree(&cvode4);
                    N_VDestroy(y); SUNContext_Free(&sunctx); return 1;
                }
                pack4(y_win4);
                if (!setup4(t_cur4)) {
                    N_VDestroy(y_win4); N_VDestroy(y);
                    SUNContext_Free(&sunctx); return 1;
                }
            }
        }

        { long int _ns=0,_nf=0,_nn=0,_nl=0;
          if (cvode4) { CVodeGetNumSteps(cvode4,&_ns); CVodeGetNumRhsEvals(cvode4,&_nf);
            CVodeGetNumNonlinSolvIters(cvode4,&_nn); CVodeGetNumLinIters(cvode4,&_nl); }
          acc_nsteps_IV+=_ns; acc_nfevals_IV+=_nf; acc_nniters_IV+=_nn; acc_nliters_IV+=_nl; }
        std::cerr << "  Phase IV: " << n_reinits4 << " CVODE reinits"
                  << "  final window: Cv[1.." << x_hi_v4 << "]"
                  << "  Ci[" << x_lo_i4 << ".." << x_hi_i4 << "]\n"
                  << "  Stats: steps=" << acc_nsteps_IV << " rhs_evals=" << acc_nfevals_IV
                  << " newton_iters=" << acc_nniters_IV << " lin_iters=" << acc_nliters_IV
                  << "  [" << (acc_nniters_IV>0?(double)acc_nliters_IV/acc_nniters_IV:0.0)
                  << " GMRES iters/Newton]\n";

        if (LS4)    SUNLinSolFree(LS4);
        if (cvode4) CVodeFree(&cvode4);
        N_VDestroy(y_win4);

    } else {
        // ── CVODE — linear multistep BDF or Adams (default) ───────────────────
        int lmm_flag = (P.lmm == 1) ? CV_ADAMS : CV_BDF;
        const char* lmm_name = (P.lmm == 1) ? "ADAMS" : "BDF";

        void* cvode_mem = CVodeCreate(lmm_flag, sunctx);
        if (!cvode_mem) { std::cerr << "Error in CVodeCreate\n"; return 1; }

        flag = CVodeSetUserData(cvode_mem, &P);
        if (flag != CV_SUCCESS) { std::cerr << "Error in CVodeSetUserData\n"; return 1; }

        flag = CVodeInit(cvode_mem, rhs_cd, t_eval[0], y);
        if (flag != CV_SUCCESS) { std::cerr << "Error in CVodeInit\n"; return 1; }

        flag = CVodeSStolerances(cvode_mem, P.rtol, P.atol);
        if (flag != CV_SUCCESS) { std::cerr << "Error in CVodeSStolerances\n"; return 1; }

        flag = CVodeSetMaxNumSteps(cvode_mem, 1000000);
        if (flag != CV_SUCCESS) { std::cerr << "Error in CVodeSetMaxNumSteps\n"; return 1; }

        if (P.max_order > 0) {
            flag = CVodeSetMaxOrd(cvode_mem, P.max_order);
            if (flag != CV_SUCCESS) { std::cerr << "Error in CVodeSetMaxOrd\n"; return 1; }
        }

        if (!make_linear_solver(cvode_mem, false)) return 1;

        std::cerr << "Solver: CVODE " << lmm_name
                  << " | linsol: " << linsol_name
                  << " | max_order: " << (P.max_order > 0 ? P.max_order : -1) << " (default)"
                  << " | N_EQ=" << N_EQ << "\n";

        // Initial point
        write_row(t_eval[0], NV_DATA_S(y), N_EQ);

        sunrealtype t_current = t_eval[0];
        for (int i = 1; i < P.n_points; ++i) {
            flag = CVode(cvode_mem, t_eval[i], y, &t_current, CV_NORMAL);
            if (flag < 0) {
                std::cerr << "CVode error at step " << i
                          << " (t=" << t_eval[i] << "): flag=" << flag << "\n";
                N_VDestroy(y); SUNLinSolFree(LS); if (A) SUNMatDestroy(A);
                CVodeFree(&cvode_mem); SUNContext_Free(&sunctx);
                if (fp_bin) std::fclose(fp_bin);
                return 1;
            }
            write_row(t_eval[i], NV_DATA_S(y), N_EQ);
        }

        { long int ns=0,nf=0,nn=0,nl=0;
          CVodeGetNumSteps(cvode_mem,&ns); CVodeGetNumRhsEvals(cvode_mem,&nf);
          CVodeGetNumNonlinSolvIters(cvode_mem,&nn); CVodeGetNumLinIters(cvode_mem,&nl);
          std::cerr << "  Stats: steps=" << ns << " rhs_evals=" << nf
                    << " newton_iters=" << nn << " lin_iters=" << nl
                    << "  [" << (nn>0?(double)nl/nn:0.0) << " GMRES iters/Newton]\n"; }
        N_VDestroy(y); SUNLinSolFree(LS); if (A) SUNMatDestroy(A);
        CVodeFree(&cvode_mem);
    }

    if (fp_bin) std::fclose(fp_bin);
    SUNContext_Free(&sunctx);
    return 0;
}
