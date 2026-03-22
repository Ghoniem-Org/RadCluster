/**
 * rate_equations.h – ODE right-hand side for the 150-equation cluster dynamics system.
 *
 * Mirrors: py_utils/rate_equations.py (RateEquations._rhs_full)
 *
 * State vector y[N_EQ = 150]:
 *   y[0 .. Nv-1]       – Cv1 .. Cv_Nv   (vacancy clusters, Nv=50)
 *   y[Nv .. Nv+Ni-1]   – Ci1 .. Ci_Ni   (interstitial clusters, Ni=100)
 */
#pragma once

#include "parameters.h"
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>

/**
 * CVODE/ARKODE-compatible RHS callback.
 * user_data must point to a Parameters struct.
 */
int rhs_cd(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

// ── Phase II sliding-window data and callbacks ────────────────────────────────

/**
 * WindowData – passed as user_data in the Phase II sliding-window CVODE session.
 *
 * State vector layout:
 *   y[0 .. x_hi_v-1]          : Cv1 .. Cv_{x_hi_v}
 *   y[x_hi_v]                  : Ci1  (always active)
 *   y[x_hi_v+1 .. x_hi_v+n-1] : Ci_{x_lo_i} .. Ci_{x_hi_i}  (n = x_hi_i-x_lo_i+1)
 *
 * Frozen clusters Ci2 .. Ci_{x_lo_i-1} are held in full_conc and contribute
 * precomputed correction sums (frozen_KLI_sum / frozen_KLV_sum) to the Ci1/Cv1
 * equations.  Ci2 is always handled explicitly regardless of x_lo_i.
 */
struct WindowData {
    const Parameters* P_full;   // immutable full parameter set (all Nv/Ni rate arrays)

    // Current window bounds (1-indexed cluster sizes)
    int x_hi_v;   // active vacancy upper bound  (<= Nv)
    int x_lo_i;   // active interstitial lower bound (>= 2; Ci1 always separate)
    int x_hi_i;   // active interstitial upper bound (<= Ni)
    int N_active; // total equations in state vector = x_hi_v + 1 + (x_hi_i - x_lo_i + 1)

    // Full concentration buffer (size Nv+Ni); always up-to-date after unpack.
    std::vector<double> full_conc;

    // Precomputed frozen correction sums (updated after each contraction event).
    // Cover clusters Ci3 .. Ci_{x_lo_i-1} only (Ci2 is always handled explicitly).
    //   frozen_KLI_sum = sum_{x=3}^{x_lo_i-1} KLI[x-1] * Ci_x
    //   frozen_KLV_sum = sum_{x=3}^{x_lo_i-1} KLV[x-1] * Ci_x
    double frozen_KLI_sum;
    double frozen_KLV_sum;

    // Ci_{x_lo_i - 1}: left-neighbour ghost for the lowest active cluster.
    // Used only when x_lo_i >= 3 (when x_lo_i==2 the neighbour is live Ci1).
    double Ci_frozen_top;

    // Diagonal preconditioner storage (size N_active); populated by prec_setup_window.
    std::vector<double> prec_diag;
};

/**
 * Phase II RHS callback.  user_data must point to a WindowData struct.
 */
int rhs_window(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/**
 * Recompute frozen_KLI_sum, frozen_KLV_sum, and Ci_frozen_top from full_conc.
 * Call after every contraction event (x_lo_i advance).
 */
void recompute_frozen_sums(WindowData& W);

/**
 * CVODE GMRES Jacobi-diagonal preconditioner callbacks.
 */
int prec_setup_window(sunrealtype t, N_Vector y, N_Vector fy,
                      sunbooleantype jok, sunbooleantype* jcurPtr,
                      sunrealtype gamma, void* user_data);

int prec_solve_window(sunrealtype t, N_Vector y, N_Vector fy,
                      N_Vector r, N_Vector z,
                      sunrealtype gamma, sunrealtype delta,
                      int lr, void* user_data);

// ── Phase I Jacobi-diagonal preconditioner callbacks ─────────────────────────
// user_data must point to a Parameters struct (the active P_win used by Phase I).
// prec_diag is stored in Parameters::prec_diag (sized to P.Nv + P.Ni).

int prec_setup_win1(sunrealtype t, N_Vector y, N_Vector fy,
                    sunbooleantype jok, sunbooleantype* jcurPtr,
                    sunrealtype gamma, void* user_data);

int prec_solve_win1(sunrealtype t, N_Vector y, N_Vector fy,
                    N_Vector r, N_Vector z,
                    sunrealtype gamma, sunrealtype delta,
                    int lr, void* user_data);

// ── Phase IV: Multithread-OpenMP window data and callbacks ────────────────────

/**
 * WindowDataOMP – Phase IV extension of the sliding-window data.
 *
 * Adds to WindowData:
 *   Cv_buf / Ci_buf : pre-allocated scratch vectors that replace the per-call
 *                     std::vector allocations in rhs_window_omp(), eliminating
 *                     heap traffic from the hot RHS path.
 *   n_omp_threads   : number of OpenMP threads (0 → honour OMP_NUM_THREADS).
 *
 * Call resize_buffers() after every CVODE reinitialisation (whenever x_hi_v
 * or the interstitial window width changes) to keep the buffers in sync.
 */
struct WindowDataOMP {
    const Parameters* P_full;

    int x_hi_v;
    int x_lo_i;
    int x_hi_i;
    int N_active;

    std::vector<double> full_conc;

    double frozen_KLI_sum;
    double frozen_KLV_sum;
    double Ci_frozen_top;

    std::vector<double> prec_diag;

    // Phase IV additions ──────────────────────────────────────────────────────
    std::vector<double> Cv_buf;   // pre-allocated [x_hi_v]
    std::vector<double> Ci_buf;   // pre-allocated [x_hi_i - x_lo_i + 1]
    int n_omp_threads;            // 0 = OMP_NUM_THREADS env var

    void resize_buffers() {
        Cv_buf.resize(static_cast<size_t>(x_hi_v));
        Ci_buf.resize(static_cast<size_t>(x_hi_i - x_lo_i + 1));
    }
};

/**
 * Phase IV RHS: same physics as rhs_window but with OpenMP-parallelised loops
 * and pre-allocated scratch buffers (no heap allocation per call).
 * user_data must point to a WindowDataOMP struct.
 */
int rhs_window_omp(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/**
 * Recompute frozen sums for WindowDataOMP (mirrors recompute_frozen_sums).
 */
void recompute_frozen_sums_omp(WindowDataOMP& W);

/**
 * Phase IV Jacobi-diagonal preconditioner callbacks (same logic as Phase II/III,
 * adapted for WindowDataOMP).
 */
int prec_setup_window_omp(sunrealtype t, N_Vector y, N_Vector fy,
                           sunbooleantype jok, sunbooleantype* jcurPtr,
                           sunrealtype gamma, void* user_data);

int prec_solve_window_omp(sunrealtype t, N_Vector y, N_Vector fy,
                           N_Vector r, N_Vector z,
                           sunrealtype gamma, sunrealtype delta,
                           int lr, void* user_data);
