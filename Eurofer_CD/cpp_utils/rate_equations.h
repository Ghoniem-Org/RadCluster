/**
 * rate_equations.h – ODE right-hand side for the Eurofer_CD cluster dynamics system.
 *
 * Mirrors: py_utils/rate_equations.py (RateEquations._rhs_decoupled / _rhs_shared)
 *
 * State vector y[N_EQ = Ni + Nv + 1]:
 *   y[0 .. Ni-1]       – Ci1 .. Ci_Ni   (SIA clusters)
 *   y[Ni .. Ni+Nv-1]   – Cv1 .. Cv_Nv   (vacancy clusters)
 *   y[Ni+Nv]           – C_He            (free He)
 *
 * Window modes operate on the SIA dimension only.
 * Vacancy clusters and free He are always fully active.
 *
 * WindowData active state vector layout (Phases II–IV):
 *   y[0 .. Nv-1]              : Cv1 .. Cv_Nv           (all vacancy, always active)
 *   y[Nv]                     : C_He                   (always active)
 *   y[Nv+1]                   : Ci1                    (always active)
 *   y[Nv+2 .. Nv+1+n]         : Ci_{x_lo_i}..Ci_{x_hi_i}  (n = x_hi_i - x_lo_i + 1)
 *
 * For Phase I the layout is:
 *   y[0 .. Nv-1]              : Cv1 .. Cv_Nv
 *   y[Nv]                     : C_He
 *   y[Nv+1 .. Nv+x_hi_i]     : Ci1 .. Ci_{x_hi_i}
 */
#pragma once

#include "parameters.h"
#include <array>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>

/**
 * CVODE/ARKODE-compatible RHS callback for the full system.
 * user_data must point to a Parameters struct.
 */
int rhs_eurofer(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

// ── Phase I / Phase II window data ───────────────────────────────────────────

/**
 * WindowData – passed as user_data for Phase II/III sliding-window CVODE sessions.
 *
 * Active state vector layout (see file header above).
 * Frozen SIA clusters Ci_{2}..Ci_{x_lo_i-1} are held in full_conc and contribute
 * precomputed correction sums (frozen_KII_sum, frozen_GII_sum, etc.) to the Ci1,
 * Cv1, and K_IclV cross-term equations.
 */
struct WindowData {
    const Parameters* P_full;   // immutable full parameter set

    // Current window bounds (1-indexed SIA cluster sizes)
    int x_hi_i;    // active SIA upper bound (<= Ni)
    int x_lo_i;    // active SIA lower bound (>= 2; Ci1 always separate)
    int N_active;  // Nv + 1 (He) + 1 (Ci1) + (x_hi_i - x_lo_i + 1)

    // Full concentration buffer (size Ni + Nv + 1); always up-to-date after unpack.
    std::vector<double> full_conc;

    // Precomputed frozen correction sums.
    // Cover frozen SIA clusters Ci_{3}..Ci_{x_lo_i-1} (Ci2 is always explicit).
    //
    //   frozen_KII_sum  = Σ_{n=3}^{x_lo_i-1} KII[n-1] * Ci_n
    //                     (Ci1 sink term due to frozen SIA clusters absorbing Ci1)
    //   frozen_K_IclV_A = Σ_{n=3}^{x_lo_i-1} K_IclV_ns[n-1] * Ci_n
    //                     (frozen contribution to Vac_recom part A)
    //   frozen_K_IclV_B = Σ_{n=3}^{x_lo_i-1} K_IclV_ni[n-1] * Ci_n
    //                     (frozen contribution to Vac_recom part B)
    //   frozen_GII_sum  = Σ_{n=3}^{x_lo_i-1} GII[n-1] * Ci_n
    //                     (frozen GII terms returning SIA to pool; Ci1 source)
    double frozen_KII_sum;
    double frozen_K_IclV_A;
    double frozen_K_IclV_B;
    double frozen_GII_sum;

    // Ci_{x_lo_i - 1}: left-neighbour ghost for the lowest active SIA cluster.
    // Used only when x_lo_i >= 3 (when x_lo_i==2 the neighbour is live Ci1).
    double Ci_frozen_top;

    // Diagonal preconditioner storage (size N_active)
    std::vector<double> prec_diag;

    // Sherman-Morrison-Woodbury rank-4 preconditioner data.
    //
    // The K_IclV separable cross-coupling adds a rank-4 correction to the
    // Jacobian that the Jacobi (diagonal) preconditioner ignores.  Storing the
    // four D^{-1}*u vectors and the inverted 4×4 Woodbury matrix allows the
    // preconditioner solve to apply the full correction in O(N) time:
    //
    //   M^{-1} r = D^{-1}r - D^{-1}U · S^{-1} · V^T D^{-1}r
    //
    // where S = I_4 + V^T D^{-1} U  (built and inverted in prec_setup).
    //
    // Vector layout (DinvU* are nonzero only in their respective subspaces):
    //   smw_DinvU1[j] = prec_diag[Nv+2+j] * γ * K_IclV_ns[x_lo-1+j] * Ci_win[j]
    //   smw_DinvU2[j] = prec_diag[Nv+2+j] * γ * K_IclV_ni[x_lo-1+j] * Ci_win[j]
    //   smw_DinvU3[k] = prec_diag[k]       * γ * Cv[k]
    //   smw_DinvU4[k] = prec_diag[k]       * γ * Cv[k] * m13[k]
    std::vector<double> smw_DinvU1;   // length n_ci_win (SIA window)
    std::vector<double> smw_DinvU2;   // length n_ci_win
    std::vector<double> smw_DinvU3;   // length Nv (vacancy)
    std::vector<double> smw_DinvU4;   // length Nv
    std::array<std::array<double,4>,4> smw_Sinv = {};  // inverse of 4×4 Woodbury matrix
};

/**
 * Phase II RHS callback.  user_data must point to a WindowData struct.
 */
int rhs_window(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/**
 * Recompute frozen sums from full_conc.
 * Call after every lower-bound contraction event.
 */
void recompute_frozen_sums(WindowData& W);

/**
 * CVODE GMRES Jacobi-diagonal preconditioner callbacks (Phase II/III).
 */
int prec_setup_window(sunrealtype t, N_Vector y, N_Vector fy,
                      sunbooleantype jok, sunbooleantype* jcurPtr,
                      sunrealtype gamma, void* user_data);

int prec_solve_window(sunrealtype t, N_Vector y, N_Vector fy,
                      N_Vector r, N_Vector z,
                      sunrealtype gamma, sunrealtype delta,
                      int lr, void* user_data);

/**
 * Phase I Jacobi-diagonal preconditioner callbacks.
 * user_data must point to a Parameters struct (with prec_diag sized to N_active).
 */
int prec_setup_win1(sunrealtype t, N_Vector y, N_Vector fy,
                    sunbooleantype jok, sunbooleantype* jcurPtr,
                    sunrealtype gamma, void* user_data);

int prec_solve_win1(sunrealtype t, N_Vector y, N_Vector fy,
                    N_Vector r, N_Vector z,
                    sunrealtype gamma, sunrealtype delta,
                    int lr, void* user_data);

// ── Phase IV: Multithread-OpenMP ──────────────────────────────────────────────

/**
 * WindowDataOMP – Phase IV extension of WindowData.
 * Pre-allocated scratch vectors eliminate per-call heap allocation.
 */
struct WindowDataOMP {
    const Parameters* P_full;

    int x_hi_i;
    int x_lo_i;
    int N_active;

    std::vector<double> full_conc;

    double frozen_KII_sum;
    double frozen_K_IclV_A;
    double frozen_K_IclV_B;
    double frozen_GII_sum;
    double Ci_frozen_top;

    std::vector<double> prec_diag;

    // Sherman-Morrison-Woodbury rank-4 preconditioner (same semantics as WindowData)
    std::vector<double> smw_DinvU1;
    std::vector<double> smw_DinvU2;
    std::vector<double> smw_DinvU3;
    std::vector<double> smw_DinvU4;
    std::array<std::array<double,4>,4> smw_Sinv = {};

    // Phase IV scratch buffers (no malloc on hot RHS path)
    std::vector<double> Ci_buf;     // pre-allocated [x_hi_i - x_lo_i + 1]
    std::vector<double> Cv_buf;     // pre-allocated [Nv]
    int n_omp_threads;

    void resize_buffers() {
        Ci_buf.resize(static_cast<size_t>(x_hi_i - x_lo_i + 1));
        Cv_buf.resize(static_cast<size_t>(P_full->Nv));
    }
};

int rhs_window_omp(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
void recompute_frozen_sums_omp(WindowDataOMP& W);

int prec_setup_window_omp(sunrealtype t, N_Vector y, N_Vector fy,
                           sunbooleantype jok, sunbooleantype* jcurPtr,
                           sunrealtype gamma, void* user_data);

int prec_solve_window_omp(sunrealtype t, N_Vector y, N_Vector fy,
                           N_Vector r, N_Vector z,
                           sunrealtype gamma, sunrealtype delta,
                           int lr, void* user_data);
