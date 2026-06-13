/**
 * sparse_jacobian.h — Colored finite-difference sparse Jacobian for KLU.
 *
 * Builds a conservative CSC sparsity pattern for the RadCluster_2_1 RHS,
 * applies Curtis–Powell–Reid (CPR) greedy column coloring, and fills the
 * SUNSparseMatrix one color group at a time using forward finite differences.
 *
 * Restricted to full_CD modes (physics_option 0 or 1). For bin_moment the
 * solver should fall back to GMRES.
 *
 * Reference: Curtis, Powell & Reid (1974), "On the estimation of sparse
 * Jacobian matrices."
 */
#pragma once

#include "parameters.h"
#include "rate_equations.h"

#include <sunmatrix/sunmatrix_sparse.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>

/**
 * Build the CSC sparsity pattern for full_CD modes (physics_option 0 or 1).
 *
 * Fills P.jac_col_ptr (size N_eq+1) and P.jac_row_idx (size nnz) with the
 * structural-nonzero indices. Pattern is conservative: it includes every
 * entry that may become nonzero across the whole integration. Returns nnz.
 *
 * Returns -1 for unsupported physics options.
 */
int build_sparsity_pattern_full_CD(Parameters& P);

/**
 * Greedy column coloring (CPR). Two columns share a color iff they have no
 * common nonzero row. Fills P.jac_colors (size N_eq) and P.jac_n_colors.
 */
void color_columns_greedy(Parameters& P);

/**
 * Group columns by color into P.jac_color_groups (flattened, indexed by
 * P.jac_color_offsets). Used by the colored FD evaluator.
 */
void build_color_groups(Parameters& P);

/**
 * Sparse Jacobian function for CVODE (CVLsJacFn signature).
 *
 * Computes J = ∂f/∂y by evaluating the RHS once per color group with all
 * columns in that group perturbed simultaneously. Total RHS evaluations =
 * P.jac_n_colors + 1 (one base + one per color group).
 */
int sparse_fd_jac(sunrealtype t, N_Vector y, N_Vector fy,
                  SUNMatrix Jac, void* user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
