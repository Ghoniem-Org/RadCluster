/**
 * solver.cpp – Eurofer_CD main C++ ODE solver.
 *
 * Drives the Eurofer_CD cluster dynamics system (Ni SIA + Nv vacancy + 1 He).
 * Invoked by py_utils/cpp_bridge.py via a parameter file:
 *
 *   solver.exe --param_file=<path>
 *
 * State vector layout:
 *   y[0..Ni-1]        SIA clusters Ci1..Ci_Ni
 *   y[Ni..Ni+Nv-1]    vacancy clusters Cv1..Cv_Nv
 *   y[Ni+Nv]          free He C_He
 *
 * Window modes (window_mode):
 *   0 = full system
 *   1 = Phase I   — upper truncation on SIA: active Ci1..Ci_{x_hi_i}
 *   2 = Phase II  — sliding window: [Ci1] + [Ci_{x_lo}..Ci_{x_hi}]
 *   3 = Phase III — constant-width sliding window
 *   4 = Phase IV  — Phase III + OpenMP
 *
 * Vacancy clusters and He are always fully active in all window modes.
 * Active state vector layout for window modes:
 *   y[0..Nv-1]    Cv1..Cv_Nv
 *   y[Nv]         C_He
 *   y[Nv+1]       Ci1                  (always active)
 *   y[Nv+2..]     Ci_{x_lo}..Ci_{x_hi} (window)
 *
 * Output: n_points rows × (1 + Ni + Nv + 1) columns (t, Ci1.., Cv1.., C_He),
 * written as raw float64 binary to a companion .bin file (fast) or space-
 * separated text to stdout (fallback).
 *
 * Build:
 *   cd Eurofer_CD/cpp_utils
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

// ── CLI argument parser ───────────────────────────────────────────────────────

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
            std::cerr << "Invalid argument (missing '='): " << arg << "\n";
            return {};
        }
        std::string key = arg.substr(2, pos - 2);
        double val = 0.0;
        try { val = std::stod(arg.substr(pos + 1)); }
        catch (...) { std::cerr << "Invalid value for " << key << "\n"; return {}; }
        props[key] = val;
    }
    return props;
}

// ── Main ──────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {

    std::map<std::string, double> args;
    std::string param_path_str;

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
    const int N_EQ = P.N_EQ;   // Ni_max + Nv + 1

    // ── Binary output file ─────────────────────────────────────────────────────
    FILE* fp_bin = nullptr;
    std::string bin_path;
    if (!param_path_str.empty()) {
        bin_path = param_path_str;
        const std::string suffix = ".txt";
        if (bin_path.size() >= suffix.size() &&
            bin_path.compare(bin_path.size() - suffix.size(),
                             suffix.size(), suffix) == 0)
            bin_path.replace(bin_path.size() - suffix.size(),
                             suffix.size(), ".bin");
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

    // Phase III/IV fallback when system too small
    if ((P.window_mode == 3 || P.window_mode == 4) &&
        P.Ni < P.window_N_thresh) {
        std::cerr << "[Phase III/IV] Ni=" << P.Ni
                  << " < N_thresh=" << P.window_N_thresh
                  << " — using full solver\n";
        P.window_mode = 0;
    }

#ifndef CD_HAVE_OPENMP
    if (P.window_mode == 4) {
        std::cerr << "[Phase IV] OpenMP not available — falling back to Phase III.\n";
        P.window_mode = 3;
    }
#endif

    // ── Time grid ─────────────────────────────────────────────────────────────
    std::vector<double> t_eval(P.n_points);
    if (P.log_time) {
        double l0 = std::log10(P.t_begin), lf = std::log10(P.t_end);
        double dt = (lf - l0) / (P.n_points - 1);
        for (int i = 0; i < P.n_points; ++i)
            t_eval[i] = std::pow(10.0, l0 + i * dt);
    } else {
        double dt = (P.t_end - P.t_begin) / (P.n_points - 1);
        for (int i = 0; i < P.n_points; ++i)
            t_eval[i] = P.t_begin + i * dt;
    }

    // ── SUNDIALS context ────────────────────────────────────────────────────────
    SUNContext sunctx;
    if (SUNContext_Create(SUN_COMM_NULL, &sunctx) != 0) {
        std::cerr << "Error creating SUNContext\n"; return 1;
    }

    // ── State vector ───────────────────────────────────────────────────────────
    N_Vector y = N_VNew_Serial(N_EQ, sunctx);
    if (!y) { std::cerr << "Error allocating N_Vector\n"; return 1; }
    for (int k = 0; k < N_EQ; ++k) NV_Ith_S(y, k) = P.y0[k];

    int flag;
    SUNMatrix       A  = nullptr;
    SUNLinearSolver LS = nullptr;
    const char*     linsol_name = "dense";

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
        // ── ARKODE ───────────────────────────────────────────────────────────
        void* ark_mem = ARKStepCreate(nullptr, rhs_eurofer, t_eval[0], y, sunctx);
        if (!ark_mem) { std::cerr << "ARKStepCreate failed\n"; return 1; }

        ARKStepSetImplicit(ark_mem);
        ARKStepSetTableNum(ark_mem,
                           static_cast<ARKODE_DIRKTableID>(P.ark_table),
                           ARKODE_ERK_NONE);
        ARKodeSetUserData(ark_mem, &P);
        ARKodeSStolerances(ark_mem, P.rtol, P.atol);
        ARKodeSetMaxNumSteps(ark_mem, 1000000);
        if (P.max_order > 0) ARKodeSetOrder(ark_mem, P.max_order);
        if (!make_linear_solver(ark_mem, true)) return 1;

        const char* tbl_name = ARKodeButcherTable_DIRKIDToName(
                                    static_cast<ARKODE_DIRKTableID>(P.ark_table));
        std::cerr << "Solver: ARKODE ARKStep DIRK"
                  << " | table: " << (tbl_name ? tbl_name : "?")
                  << " | linsol: " << linsol_name
                  << " | N_EQ=" << N_EQ << "\n";

        write_row(t_eval[0], NV_DATA_S(y), N_EQ);
        sunrealtype t_cur = t_eval[0];
        for (int i = 1; i < P.n_points; ++i) {
            flag = ARKodeEvolve(ark_mem, t_eval[i], y, &t_cur, ARK_NORMAL);
            if (flag < 0) {
                std::cerr << "ARKodeEvolve error step " << i
                          << " t=" << t_eval[i] << " flag=" << flag << "\n";
                N_VDestroy(y); SUNLinSolFree(LS); if (A) SUNMatDestroy(A);
                ARKodeFree(&ark_mem); SUNContext_Free(&sunctx);
                if (fp_bin) std::fclose(fp_bin);
                return 1;
            }
            write_row(t_eval[i], NV_DATA_S(y), N_EQ);
        }
        N_VDestroy(y); SUNLinSolFree(LS); if (A) SUNMatDestroy(A);
        ARKodeFree(&ark_mem);

    } else if (P.window_mode == 0) {
        // ── CVODE — full system ────────────────────────────────────────────────
        int lmm0 = (P.lmm == 1) ? CV_ADAMS : CV_BDF;
        const char* lmm0_name = (P.lmm == 1) ? "ADAMS" : "BDF";

        void* cvode0 = CVodeCreate(lmm0, sunctx);
        if (!cvode0) { std::cerr << "CVodeCreate failed\n"; return 1; }
        CVodeSetUserData(cvode0, &P);
        CVodeInit(cvode0, rhs_eurofer, t_eval[0], y);
        CVodeSStolerances(cvode0, P.rtol, P.atol);
        CVodeSetMaxNumSteps(cvode0, 10000000);
        if (P.max_order > 0) CVodeSetMaxOrd(cvode0, P.max_order);
        if (!make_linear_solver(cvode0, false)) return 1;

        std::cerr << "Solver: CVODE " << lmm0_name
                  << " FULL | linsol: " << linsol_name
                  << " | N_EQ=" << N_EQ << "\n";

        write_row(t_eval[0], NV_DATA_S(y), N_EQ);
        sunrealtype t_cur = t_eval[0];
        for (int i = 1; i < P.n_points; ++i) {
            flag = CVode(cvode0, t_eval[i], y, &t_cur, CV_NORMAL);
            if (flag < 0) {
                std::cerr << "CVode error step " << i
                          << " t=" << t_eval[i] << " flag=" << flag << "\n";
                N_VDestroy(y); SUNLinSolFree(LS); if (A) SUNMatDestroy(A);
                CVodeFree(&cvode0); SUNContext_Free(&sunctx);
                if (fp_bin) std::fclose(fp_bin);
                return 1;
            }
            write_row(t_eval[i], NV_DATA_S(y), N_EQ);
        }
        {
            long int ns=0, nf=0, nn=0, nl=0;
            CVodeGetNumSteps(cvode0,&ns); CVodeGetNumRhsEvals(cvode0,&nf);
            CVodeGetNumNonlinSolvIters(cvode0,&nn); CVodeGetNumLinIters(cvode0,&nl);
            std::cerr << "  Stats: steps=" << ns << " rhs_evals=" << nf
                      << " newton_iters=" << nn << " lin_iters=" << nl << "\n";
        }
        N_VDestroy(y); SUNLinSolFree(LS); if (A) SUNMatDestroy(A);
        CVodeFree(&cvode0);

    } else {
        // ── Window modes (Phase I / II / III / IV) ────────────────────────────
        //
        // The window operates on SIA clusters only.
        // Active state vector layout:
        //   y[0..Nv-1]    Cv1..Cv_Nv        (always all active)
        //   y[Nv]         C_He               (always active)
        //   y[Nv+1]       Ci1                (always active)
        //   y[Nv+2+j]     Ci_{x_lo+j}        (active window)
        //
        // For Phase I the window goes from Ci1..Ci_{x_hi_i} (no lower contraction).
        //
        // Full concentration buffer layout (full_conc[N_EQ]):
        //   full_conc[k]       = Ci_{k+1}  for k = 0..Ni-1
        //   full_conc[Ni+k]    = Cv_{k+1}  for k = 0..Nv-1
        //   full_conc[Ni+Nv]   = C_He

        const int lmm_win = (P.lmm == 1) ? CV_ADAMS : CV_BDF;
        const char* lmm_win_name = (P.lmm == 1) ? "ADAMS" : "BDF";
        const int win_width = std::max(P.window_width, 2);

        // Initial window bounds
        int x_hi_i = std::min(std::max(P.window_w0_i, 2), P.Ni);
        int x_lo_i = 2;  // Phase I: always from Ci2; Phase II+: starts at 2

        // For Phase I (upper truncation only), x_lo_i is always 2 (= all from Ci2)
        // For Phase II/III, x_lo_i can advance via contraction

        // WindowData for Phases II–IV
        WindowData   Wd;
        WindowDataOMP Womp;

        auto init_window_data = [&](auto& W) {
            W.P_full          = &P;
            W.x_hi_i          = x_hi_i;
            W.x_lo_i          = x_lo_i;
            W.N_active        = P.Nv + 1 + 1 + (x_hi_i - x_lo_i + 1);
            W.frozen_KII_sum  = 0.0;
            W.frozen_K_IclV_A = 0.0;
            W.frozen_K_IclV_B = 0.0;
            W.frozen_GII_sum  = 0.0;
            W.Ci_frozen_top   = 0.0;
            W.full_conc.assign(P.N_EQ, 0.0);
            for (int k = 0; k < P.N_EQ; ++k) W.full_conc[k] = P.y0[k];
        };

        if (P.window_mode >= 2) {
            if (P.window_mode == 4) {
                init_window_data(Womp);
                Womp.n_omp_threads = P.window_omp_threads;
                Womp.resize_buffers();
            } else {
                init_window_data(Wd);
            }
        }

        // Full concentration buffer (used by Phase I too)
        std::vector<double> full_conc(P.N_EQ);
        for (int k = 0; k < P.N_EQ; ++k) full_conc[k] = P.y0[k];

        // Active size: Nv (vacancy) + 1 (He) + 1 (Ci1) + window_len (Ci window)
        auto active_n = [&]() -> int {
            return P.Nv + 1 + 1 + (x_hi_i - x_lo_i + 1);
        };

        // Pack: full_conc → y_win (active window vector)
        // Layout: [Cv1..Cv_Nv | C_He | Ci1 | Ci_{x_lo}..Ci_{x_hi}]
        auto pack_win = [&](N_Vector ya) {
            for (int k = 0; k < P.Nv; ++k)
                NV_Ith_S(ya, k) = full_conc[P.Ni + k];          // Cv
            NV_Ith_S(ya, P.Nv) = full_conc[P.Ni + P.Nv];        // C_He
            NV_Ith_S(ya, P.Nv + 1) = full_conc[0];               // Ci1
            int n = x_hi_i - x_lo_i + 1;
            for (int j = 0; j < n; ++j)
                NV_Ith_S(ya, P.Nv + 2 + j) = full_conc[x_lo_i - 1 + j]; // Ci window
        };

        auto unpack_win = [&](N_Vector ya) {
            for (int k = 0; k < P.Nv; ++k)
                full_conc[P.Ni + k]     = NV_Ith_S(ya, k);
            full_conc[P.Ni + P.Nv]     = NV_Ith_S(ya, P.Nv);
            full_conc[0]               = NV_Ith_S(ya, P.Nv + 1);
            int n = x_hi_i - x_lo_i + 1;
            for (int j = 0; j < n; ++j)
                full_conc[x_lo_i - 1 + j] = NV_Ith_S(ya, P.Nv + 2 + j);
        };

        // Sync pack/unpack for WindowData (uses its own full_conc)
        auto pack_wd = [&](N_Vector ya, auto& W) {
            for (int k = 0; k < P.Nv; ++k)
                NV_Ith_S(ya, k)        = W.full_conc[P.Ni + k];
            NV_Ith_S(ya, P.Nv)        = W.full_conc[P.Ni + P.Nv];
            NV_Ith_S(ya, P.Nv + 1)    = W.full_conc[0];
            int n = W.x_hi_i - W.x_lo_i + 1;
            for (int j = 0; j < n; ++j)
                NV_Ith_S(ya, P.Nv + 2 + j) = W.full_conc[W.x_lo_i - 1 + j];
        };

        auto unpack_wd = [&](N_Vector ya, auto& W) {
            for (int k = 0; k < P.Nv; ++k)
                W.full_conc[P.Ni + k]     = NV_Ith_S(ya, k);
            W.full_conc[P.Ni + P.Nv]     = NV_Ith_S(ya, P.Nv);
            W.full_conc[0]               = NV_Ith_S(ya, P.Nv + 1);
            int n = W.x_hi_i - W.x_lo_i + 1;
            for (int j = 0; j < n; ++j)
                W.full_conc[W.x_lo_i - 1 + j] = NV_Ith_S(ya, P.Nv + 2 + j);
        };

        // Allocate active N_Vector
        N_Vector y_win = N_VNew_Serial(active_n(), sunctx);
        if (!y_win) { std::cerr << "Error allocating window N_Vector\n"; return 1; }

        if (P.window_mode >= 2) {
            if (P.window_mode == 4) pack_wd(y_win, Womp);
            else                     pack_wd(y_win, Wd);
        } else {
            pack_win(y_win);
        }

        // Phase I: rhs_cd-style with reduced Ni; use a sub-Parameters copy
        Parameters P_win1 = P;   // for Phase I only
        P_win1.Ni   = x_hi_i;
        P_win1.N_EQ = x_hi_i + P.Nv + 1;

        // CVODE handle and linear solver for the window
        void*           cvode_win = nullptr;
        SUNLinearSolver LS_win    = nullptr;
        int n_reinits = 0;
        double h_last = 0.0;
        long int acc_steps=0, acc_fevals=0, acc_niters=0, acc_liters=0;

        // Prec type for GMRES
        const int prec_type = (P.window_prec > 0) ? SUN_PREC_LEFT : SUN_PREC_NONE;

        // (Re)create CVODE + GMRES for the current window size
        auto setup_cvode = [&](double t0) -> bool {
            if (LS_win)    { SUNLinSolFree(LS_win); LS_win = nullptr; }
            if (cvode_win) {
                long int _s=0,_f=0,_n=0,_l=0;
                CVodeGetNumSteps(cvode_win,&_s); CVodeGetNumRhsEvals(cvode_win,&_f);
                CVodeGetNumNonlinSolvIters(cvode_win,&_n); CVodeGetNumLinIters(cvode_win,&_l);
                acc_steps+=_s; acc_fevals+=_f; acc_niters+=_n; acc_liters+=_l;
                CVodeFree(&cvode_win); cvode_win = nullptr;
            }

            int N_act = active_n();
            cvode_win = CVodeCreate(lmm_win, sunctx);
            if (!cvode_win) { std::cerr << "Window CVodeCreate failed\n"; return false; }

            // Select RHS and user_data based on window mode
            auto rhs_fn = rhs_eurofer;
            void* ud = nullptr;

            if (P.window_mode == 1) {
                // Phase I: rhs_eurofer with reduced Ni (copy P_win1)
                P_win1.Ni   = x_hi_i;
                P_win1.N_EQ = x_hi_i + P.Nv + 1;
                P_win1.prec_diag.resize(N_act);
                rhs_fn = rhs_eurofer;
                ud     = &P_win1;
            } else if (P.window_mode == 2 || P.window_mode == 3) {
                Wd.x_hi_i  = x_hi_i;
                Wd.x_lo_i  = x_lo_i;
                Wd.N_active = N_act;
                rhs_fn = rhs_window;
                ud     = &Wd;
            } else if (P.window_mode == 4) {
                Womp.x_hi_i  = x_hi_i;
                Womp.x_lo_i  = x_lo_i;
                Womp.N_active = N_act;
                Womp.resize_buffers();
                rhs_fn = rhs_window_omp;
                ud     = &Womp;
            }

            if (CVodeSetUserData(cvode_win, ud)                   != CV_SUCCESS) return false;
            if (CVodeInit(cvode_win, rhs_fn, t0, y_win)           != CV_SUCCESS) return false;
            if (CVodeSStolerances(cvode_win, P.rtol, P.atol)      != CV_SUCCESS) return false;
            if (CVodeSetMaxNumSteps(cvode_win, 10000000)           != CV_SUCCESS) return false;
            if (P.max_order > 0)
                CVodeSetMaxOrd(cvode_win, P.max_order);
            if (h_last > 0.0) CVodeSetMaxStep(cvode_win, h_last);

            LS_win = SUNLinSol_SPGMR(y_win, prec_type, P.window_gmres_maxl, sunctx);
            if (!LS_win) { std::cerr << "Window SUNLinSol_SPGMR failed\n"; return false; }
            if (CVodeSetLinearSolver(cvode_win, LS_win, nullptr) != CV_SUCCESS) return false;

            if (P.window_prec > 0) {
                if (P.window_mode == 1) {
                    CVodeSetPreconditioner(cvode_win, prec_setup_win1, prec_solve_win1);
                } else if (P.window_mode == 4) {
                    CVodeSetPreconditioner(cvode_win, prec_setup_window_omp, prec_solve_window_omp);
                } else {
                    CVodeSetPreconditioner(cvode_win, prec_setup_window, prec_solve_window);
                }
            }

            ++n_reinits;
            std::cerr << "[window] reinit #" << n_reinits
                      << "  x_hi_i=" << x_hi_i
                      << "  x_lo_i=" << x_lo_i
                      << "  N_active=" << N_act << "  t0=" << t0 << "\n";
            return true;
        };

        if (!setup_cvode(t_eval[0])) {
            N_VDestroy(y_win); N_VDestroy(y); SUNContext_Free(&sunctx); return 1;
        }

        // Print solver configuration
        const char* phase_name =
            (P.window_mode == 1) ? "WINDOW-I (upper-trunc)" :
            (P.window_mode == 2) ? "WINDOW-II (sliding)"   :
            (P.window_mode == 3) ? "WINDOW-III (const-width)" :
                                    "WINDOW-IV (OpenMP)";
        std::cerr << "Solver: CVODE " << lmm_win_name << " " << phase_name
                  << " | linsol: GMRES(maxl=" << P.window_gmres_maxl << ")"
                  << (P.window_prec ? " + SMW-Jacobi" : "")
                  << " | w0_i=" << P.window_w0_i
                  << " | width=" << win_width
                  << " | C_expand=" << P.window_C_expand
                  << " | Nv=" << P.Nv << " Ni=" << P.Ni << "\n";

        // Output initial point (full N_EQ row)
        write_row(t_eval[0], full_conc.data(), N_EQ);

        // ── Progress printer: called after each output time step ──────────────
        // full_conc layout: [0..Ni-1]=Ci1..CiNi, [Ni..Ni+Nv-1]=Cv1..CvNv, [Ni+Nv]=C_He
        auto print_progress = [&](int step_i, double t, const std::vector<double>& fc) {
            const double floor_val = P.C_floor;
            // Point defects
            double Ci1  = fc[0];
            double Cv1  = fc[P.Ni];
            double CHe  = fc[P.Ni + P.Nv];
            // SIA cluster stats (size >= 2; skip monomer)
            double sia_n1 = 0.0, sia_n2 = 0.0;
            int    sia_max = 1;
            for (int k = 1; k < P.Ni; ++k) {
                double c = fc[k];
                if (c > floor_val) {
                    sia_n1 += c;
                    sia_n2 += (k + 1) * c;
                    sia_max = k + 1;
                }
            }
            double sia_mean = (sia_n1 > 0.0) ? sia_n2 / sia_n1 : 0.0;
            // Vacancy cluster stats (size >= 2; skip monomer)
            double vac_n1 = 0.0, vac_n2 = 0.0;
            int    vac_max = 1;
            for (int k = 1; k < P.Nv; ++k) {
                double c = fc[P.Ni + k];
                if (c > floor_val) {
                    vac_n1 += c;
                    vac_n2 += (k + 1) * c;
                    vac_max = k + 1;
                }
            }
            double vac_mean = (vac_n1 > 0.0) ? vac_n2 / vac_n1 : 0.0;

            std::cerr << std::scientific << std::setprecision(3)
                      << "[" << step_i << "/" << P.n_points - 1 << " t=" << t << "]"
                      << "  Ci1=" << Ci1 << "  Cv1=" << Cv1 << "  CHe=" << CHe
                      << "  |  <n_SIA>=" << std::fixed << std::setprecision(1) << sia_mean
                      << " max_SIA=" << sia_max
                      << "  |  <m_vac>=" << vac_mean
                      << " max_vac=" << vac_max
                      << "\n";
            std::cerr.flush();
        };

        print_progress(0, t_eval[0], full_conc);

        sunrealtype t_cur = t_eval[0];

        for (int i = 1; i < P.n_points; ++i) {
            flag = CVode(cvode_win, t_eval[i], y_win, &t_cur, CV_NORMAL);
            if (flag < 0) {
                std::cerr << "Window CVode error step " << i
                          << " t=" << t_eval[i] << " flag=" << flag << "\n";
                N_VDestroy(y_win);
                if (LS_win) SUNLinSolFree(LS_win);
                if (cvode_win) CVodeFree(&cvode_win);
                N_VDestroy(y); SUNContext_Free(&sunctx);
                if (fp_bin) std::fclose(fp_bin);
                return 1;
            }

            // Unpack into the appropriate full_conc buffer
            if (P.window_mode >= 2) {
                if (P.window_mode == 4) unpack_wd(y_win, Womp);
                else                     unpack_wd(y_win, Wd);
                // Copy to local full_conc for write_row
                const auto& wfc = (P.window_mode == 4) ? Womp.full_conc : Wd.full_conc;
                for (int k = 0; k < P.N_EQ; ++k) full_conc[k] = wfc[k];
            } else {
                unpack_win(y_win);
            }

            write_row(t_eval[i], full_conc.data(), N_EQ);
            print_progress(i, t_eval[i], full_conc);

            if (i % P.window_check_every != 0) continue;

            bool changed = false;

            // ── Upper expansion: expand x_hi_i when Ci_{x_hi} exceeds threshold ──
            if (x_hi_i < P.Ni &&
                full_conc[x_hi_i - 1] > P.window_C_expand) {
                int new_hi = (P.window_expand_factor > 1.0)
                    ? std::max(static_cast<int>(x_hi_i * P.window_expand_factor),
                               x_hi_i + P.window_expand_pad)
                    : x_hi_i + P.window_expand_pad;
                x_hi_i = std::min(new_hi, P.Ni);
                changed = true;
            }

            // ── Phase II lower contraction ──────────────────────────────────────
            if (P.window_mode == 2 && P.window_C_contract > 0.0 &&
                (x_hi_i - x_lo_i + 1) > P.window_min_active_i) {

                N_Vector ydot_tmp = N_VNew_Serial(active_n(), sunctx);
                if (ydot_tmp) {
                    rhs_window(t_cur, y_win, ydot_tmp, &Wd);
                    // x_lo is at position Nv+2 in y_win
                    const double Ci_lo  = NV_Ith_S(y_win, P.Nv + 2);
                    const double dCi_lo = NV_Ith_S(ydot_tmp, P.Nv + 2);
                    N_VDestroy(ydot_tmp);

                    double rel = std::abs(dCi_lo) / std::max(Ci_lo, P.C_floor);
                    if (rel < P.window_C_contract) {
                        // Extend frozen sums and advance lower bound
                        Wd.frozen_KII_sum  += P.KII[x_lo_i - 1]       * Ci_lo;
                        Wd.frozen_K_IclV_A += P.K_IclV_ns[x_lo_i - 1] * Ci_lo;
                        Wd.frozen_K_IclV_B += P.K_IclV_ni[x_lo_i - 1] * Ci_lo;
                        Wd.frozen_GII_sum  += P.GII[x_lo_i - 1]        * Ci_lo;
                        Wd.Ci_frozen_top    = Ci_lo;
                        x_lo_i += 1;
                        changed = true;
                    }
                }
            }

            // ── Phase III / IV lower bound: constant-width ────────────────────
            if ((P.window_mode == 3 || P.window_mode == 4) &&
                t_cur > P.window_t_start) {
                int new_lo = std::max(2, x_hi_i - win_width + 1);
                if (new_lo > x_lo_i) {
                    // Advance lower bound to new_lo and recompute frozen sums
                    auto& Wref = (P.window_mode == 4) ? (WindowData&)(Womp) : Wd;
                    // We need to build a temporary WindowData-compatible object.
                    // For simplicity, use the existing Wd or Womp directly.
                    // Recompute frozen sums after advancing x_lo_i.
                    x_lo_i = new_lo;
                    if (P.window_mode == 4) {
                        Womp.x_lo_i = x_lo_i;
                        recompute_frozen_sums_omp(Womp);
                    } else {
                        Wd.x_lo_i = x_lo_i;
                        recompute_frozen_sums(Wd);
                    }
                    changed = true;
                }
            }

            if (changed) {
                CVodeGetLastStep(cvode_win, &h_last);
                N_VDestroy(y_win);
                y_win = N_VNew_Serial(active_n(), sunctx);
                if (!y_win) {
                    std::cerr << "Error allocating expanded window N_Vector\n";
                    if (LS_win) SUNLinSolFree(LS_win);
                    if (cvode_win) CVodeFree(&cvode_win);
                    N_VDestroy(y); SUNContext_Free(&sunctx); return 1;
                }

                // Re-pack from appropriate full_conc
                if (P.window_mode >= 2) {
                    if (P.window_mode == 4) {
                        Womp.x_hi_i = x_hi_i;
                        Womp.x_lo_i = x_lo_i;
                        Womp.N_active = active_n();
                        pack_wd(y_win, Womp);
                    } else {
                        Wd.x_hi_i = x_hi_i;
                        Wd.x_lo_i = x_lo_i;
                        Wd.N_active = active_n();
                        pack_wd(y_win, Wd);
                    }
                } else {
                    pack_win(y_win);
                }

                if (!setup_cvode(t_cur)) {
                    N_VDestroy(y_win); N_VDestroy(y);
                    SUNContext_Free(&sunctx); return 1;
                }
            }
        }

        {
            long int _s=0,_f=0,_n=0,_l=0;
            if (cvode_win) {
                CVodeGetNumSteps(cvode_win,&_s); CVodeGetNumRhsEvals(cvode_win,&_f);
                CVodeGetNumNonlinSolvIters(cvode_win,&_n); CVodeGetNumLinIters(cvode_win,&_l);
            }
            acc_steps+=_s; acc_fevals+=_f; acc_niters+=_n; acc_liters+=_l;
        }
        std::cerr << "  Window solver: " << n_reinits << " CVODE reinits"
                  << "  final: Ci[" << x_lo_i << ".." << x_hi_i << "]\n"
                  << "  Stats: steps=" << acc_steps << " rhs_evals=" << acc_fevals
                  << " newton_iters=" << acc_niters << " lin_iters=" << acc_liters
                  << "  [" << (acc_niters>0?(double)acc_liters/acc_niters:0.0)
                  << " GMRES/Newton]\n";

        if (LS_win)    SUNLinSolFree(LS_win);
        if (cvode_win) CVodeFree(&cvode_win);
        N_VDestroy(y_win);
        N_VDestroy(y);
    }

    SUNContext_Free(&sunctx);
    if (fp_bin) std::fclose(fp_bin);
    return 0;
}
