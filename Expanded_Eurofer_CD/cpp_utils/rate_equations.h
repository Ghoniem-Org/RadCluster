/**
 * rate_equations.h — ODE RHS declarations for Expanded_Eurofer_CD.
 *
 * Two RHS callbacks are provided:
 *   rhs_full_CD()    — full per-size master equations (Eqs. 152, 155, 157)
 *                      dispatches between Case 1 (fusion) and Case 2 (fission)
 *                      depending on Parameters.he_mode.
 *   rhs_bin_moment() — Chapter 9 size-bin moment RHS (Eqs. 193-208)
 *
 * Window-solver data structures are shared with the solver.cpp infrastructure
 * inherited from Eurofer_CD.
 *
 * Physics reference: Ghoniem (2026), Rate_Equations.pdf.
 */
#pragma once

#include "parameters.h"

#ifdef __cplusplus
extern "C" {
#endif

/* SUNDIALS N_Vector */
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus
}
#endif

// ── User data struct passed to SUNDIALS RHS ───────────────────────────────────

struct UserData {
    const Parameters* P;
    // Dynamic window state (used by cpp_sliding_win and sliding_OpenMP)
    int x_lo_i;    // lower active SIA index (0-based)
    int x_hi_i;    // upper active SIA index (inclusive, 0-based)
    int x_hi_v;    // upper active vacancy index (always = M-1)
    bool window_active;
};

// ── RHS callbacks ─────────────────────────────────────────────────────────────

/**
 * Full per-size RHS (Eqs. 152, 155, 157).
 *
 * Dispatches between:
 *   he_mode == 0 (Case 2, fission/decoupled, Eq. 175)
 *   he_mode == 1 (Case 1, fusion/mean-field, Eq. 174)
 */
int rhs_full_CD(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/**
 * Size-bin moment RHS (Chapter 9, Eqs. 193-208).
 * Piecewise-constant closure (Eq. 198-200).
 */
int rhs_bin_moment(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/**
 * Jacobi diagonal preconditioner setup (for GMRES, window mode).
 */
int prec_setup(sunrealtype t, N_Vector y, N_Vector fy,
               sunbooleantype jok, sunbooleantype* jcurPtr,
               sunrealtype gamma, void* user_data);

int prec_solve(sunrealtype t, N_Vector y, N_Vector fy,
               N_Vector r, N_Vector z,
               sunrealtype gamma, sunrealtype delta,
               int lr, void* user_data);
