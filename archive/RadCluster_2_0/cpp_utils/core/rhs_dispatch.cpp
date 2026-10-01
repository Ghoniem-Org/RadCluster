/**
 * rhs_dispatch.cpp — RHS-dispatch scaffolding and GMRES preconditioners.
 *
 * This is the host-INDEPENDENT half of the former rate_equations.cpp.  It holds
 * the solver-side scaffolding that is shared by every material model:
 *
 *   • prec_setup() / prec_solve()  — preconditioner dispatch entry points
 *                                    invoked by SUNDIALS (CVODE + SPGMR).
 *   • Jacobi preconditioner        — diagonal scaling (legacy, prec_type==0).
 *   • Woodbury preconditioner      — bordered-banded SMW (prec_type==1),
 *                                    LAPACK dgbtrf/dgbtrs/dgetrf/dgetrs.
 *
 * The EUROFER-specific P1-P8 rate-equation kernels (rhs_full_CD / rhs_case1 /
 * rhs_case2 / rhs_bin_moment) live in
 *   cpp_utils/materials/eurofer97/rate_kernels.cpp
 * and are reached here only through the RHS function pointer carried in
 * UserData.rhs_fn — the preconditioner FD-probes whatever RHS the active
 * material model installed, so this file has no compile-time dependency on
 * the EUROFER kernel arithmetic.
 *
 * NOTE ON THE REORGANIZATION (RadCluster_2_0 core/materials split):
 *   The byte content below is relocated verbatim from rate_equations.cpp
 *   lines 1827-2190.  No arithmetic, no LAPACK call, and no control flow was
 *   altered — this is a pure file-level relocation.
 *
 * Physics reference: Ghoniem (2026), Rate_Equations.pdf.
 * Preconditioner derivation: Docs/Formulation/Jacobian_Preconditioner.tex.
 */

#include "rate_equations.h"
#include <nvector/nvector_serial.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

// ══════════════════════════════════════════════════════════════════════════════
// Preconditioner: Jacobi (legacy) and Woodbury (bordered-banded)
//
// The Woodbury preconditioner exploits the fact that the Jacobian has the
// structure  J = T + U·V^T  where:
//   T = banded part (half-bandwidth prec_bw)
//   U = N_eq × r dense columns from mobile species coupling
//   V = selector matrix [e_{j1}, ..., e_{jr}]
//
// The Newton correction system (I - γJ)x = r is solved via SMW:
//   M = T̂ - γ U V^T,  T̂ = I - γT
//   M^{-1} = T̂^{-1} + T̂^{-1} U S^{-1} V^T T̂^{-1}
//   where S = -I/γ_scale + V^T T̂^{-1} U   (r × r Schur complement)
//
// LAPACK band storage (dgbtrf):
//   A(i,j) is stored at band[kl + ku + i - j][j]  (0-indexed)
//   Total rows = 2*kl + ku + 1 = 3*bw + 1  (kl = ku = bw)
// ══════════════════════════════════════════════════════════════════════════════

// Forward declarations for LAPACK band and dense solvers.
// The Woodbury preconditioner needs LAPACK (dgbtrf/dgbtrs/dgetrf/dgetrs).
// When LAPACK is unavailable (CD_HAVE_LAPACK undefined), the entire Woodbury
// section is excluded and the dispatcher falls back to Jacobi.
#ifdef CD_HAVE_LAPACK
extern "C" {
    void dgbtrf_(const int* m, const int* n, const int* kl, const int* ku,
                 double* ab, const int* ldab, int* ipiv, int* info);
    void dgbtrs_(const char* trans, const int* n, const int* kl, const int* ku,
                 const int* nrhs, const double* ab, const int* ldab,
                 const int* ipiv, double* b, const int* ldb, int* info);
    void dgetrf_(const int* m, const int* n, double* a, const int* lda,
                 int* ipiv, int* info);
    void dgetrs_(const char* trans, const int* n, const int* nrhs,
                 const double* a, const int* lda, const int* ipiv,
                 double* b, const int* ldb, int* info);
}
#endif

// ── Jacobi preconditioner (legacy, prec_type==0) ─────────────────────────────

static int prec_setup_jacobi(sunrealtype /*t*/, N_Vector /*y*/, N_Vector /*fy*/,
                              sunbooleantype /*jok*/, sunbooleantype* jcurPtr,
                              sunrealtype /*gamma*/, void* /*user_data*/) {
    *jcurPtr = SUNTRUE;
    return 0;
}

static int prec_solve_jacobi(sunrealtype /*t*/, N_Vector /*y*/, N_Vector /*fy*/,
                              N_Vector r, N_Vector z,
                              sunrealtype gamma, sunrealtype /*delta*/,
                              int /*lr*/, void* user_data) {
    const UserData*   ud = static_cast<const UserData*>(user_data);
    const Parameters& P  = *ud->P;
    const double* rv = N_VGetArrayPointer_Serial(r);
    double*       zv = N_VGetArrayPointer_Serial(z);
    if (P.prec_diag.empty()) {
        for (int k = 0; k < P.N_eq; ++k) zv[k] = rv[k];
        return 0;
    }
    for (int k = 0; k < P.N_eq; ++k) {
        double d = 1.0 - gamma * P.prec_diag[k];
        zv[k] = rv[k] / (std::abs(d) > 1e-20 ? d : 1e-20);
    }
    return 0;
}

// ── Woodbury preconditioner (prec_type==1) ───────────────────────────────────
#ifdef CD_HAVE_LAPACK

// Finite-difference perturbation size (same as CVODE internal formula)
static inline double fd_delta(double yj, double atol) {
    const double srur = std::sqrt(2.2e-16);  // sqrt(machine epsilon)
    double dy = srur * std::max(std::abs(yj), std::max(atol, 1e-30));
    // Ensure (yj + dy) - yj is exactly representable
    volatile double temp = yj + dy;
    return temp - yj;
}

static int prec_setup_woodbury(sunrealtype t, N_Vector y, N_Vector fy,
                                sunbooleantype jok, sunbooleantype* jcurPtr,
                                sunrealtype gamma, void* user_data) {
    UserData*   ud = static_cast<UserData*>(user_data);
    Parameters& P  = *ud->P;
    const int N    = P.N_eq;
    const int bw   = P.prec_bw;
    const int r    = P.prec_rank;

    // ── Reuse check: if CVODE says the old Jacobian is still OK (jok=TRUE)
    // and the band storage has been initialised, skip the expensive rebuild.
    // CVODE will pass jok=FALSE when it wants a fresh Jacobian (e.g. after
    // a convergence failure or at the start of a new step).
    if (jok && !P.prec_band.empty()) {
        *jcurPtr = SUNFALSE;   // tell CVODE we kept the old Jacobian
        P.prec_gamma = gamma;  // update γ (step size may have changed)
        return 0;
    }

    *jcurPtr = SUNTRUE;
    P.prec_gamma = gamma;

    double* ydata = N_VGetArrayPointer_Serial(y);

    // ── Allocate storage on first call ───────────────────────────────────────
    const int ldab = 3 * bw + 1;   // LAPACK dgbtrf: 2*kl + ku + 1 with kl=ku=bw
    if (static_cast<int>(P.prec_band.size()) != ldab * N) {
        P.prec_band.resize(ldab * N, 0.0);
        P.prec_Tinv_U.resize(N * r, 0.0);
        P.prec_schur.resize(r * r, 0.0);
        P.prec_ipiv_band.resize(N, 0);
        P.prec_ipiv_schur.resize(r, 0);
        P.prec_f0.resize(N, 0.0);
        P.prec_work.resize(N, 0.0);
        P.prec_y_save.resize(N, 0.0);
        P.prec_deltas.resize(N, 0.0);
        P.prec_f_pert.resize(N, 0.0);
    }

    // ── Base RHS evaluation ──────────────────────────────────────────────────
    // fy already contains f(t,y); copy it so we can reuse fy as scratch
    const double* f0_ptr = N_VGetArrayPointer_Serial(fy);
    for (int k = 0; k < N; ++k) P.prec_f0[k] = f0_ptr[k];

    // Save the original y for restoration after each probe
    for (int k = 0; k < N; ++k) P.prec_y_save[k] = ydata[k];

    // ── Step 1: Build banded part T via Curtis-Powell-Reid coloring ──────────
    // Columns j with j mod (2*bw+1) == c can be probed simultaneously
    // because they are more than bw apart → only (2*bw+1) RHS evaluations.
    const int n_colors = 2 * bw + 1;

    // Zero the band storage
    std::fill(P.prec_band.begin(), P.prec_band.end(), 0.0);
    std::fill(P.prec_deltas.begin(), P.prec_deltas.end(), 0.0);

    for (int c = 0; c < n_colors; ++c) {
        // Perturb all columns in this color group
        for (int j = c; j < N; j += n_colors) {
            double dy = fd_delta(ydata[j], P.atol);
            P.prec_deltas[j] = dy;
            ydata[j] += dy;
        }

        // Evaluate perturbed RHS (writes into fy, which we treat as scratch)
        ud->rhs_fn(t, y, fy, user_data);
        const double* f_tmp = N_VGetArrayPointer_Serial(fy);
        for (int k = 0; k < N; ++k) P.prec_f_pert[k] = f_tmp[k];

        // Extract banded Jacobian entries and form T̂ = I - γT
        for (int j = c; j < N; j += n_colors) {
            double inv_dy = 1.0 / P.prec_deltas[j];
            int klo = std::max(0, j - bw);
            int khi = std::min(N - 1, j + bw);
            for (int k = klo; k <= khi; ++k) {
                double Jkj = (P.prec_f_pert[k] - P.prec_f0[k]) * inv_dy;
                // LAPACK column-major band storage:
                //   A(i,j) is at band[kl + ku + i - j, j]
                //   = band[(kl + ku + i - j) + j * ldab]
                int band_row = 2 * bw + (k - j);   // kl + ku + i - j
                P.prec_band[band_row + j * ldab] = (k == j)
                    ? (1.0 - gamma * Jkj)
                    : (-gamma * Jkj);
            }
        }

        // Restore y for this color group
        for (int j = c; j < N; j += n_colors) {
            ydata[j] = P.prec_y_save[j];
            P.prec_deltas[j] = 0.0;
        }
    }

    // ── Step 2: Factor banded T̂ via LAPACK dgbtrf ────────────────────────────
    {
        int info = 0;
        int kl = bw, ku = bw;
        dgbtrf_(&N, &N, &kl, &ku, P.prec_band.data(), &ldab,
                P.prec_ipiv_band.data(), &info);
        if (info != 0) {
            std::cerr << "[prec_setup_woodbury] dgbtrf failed, info=" << info << "\n";
            // Restore fy before returning
            double* fy_data = N_VGetArrayPointer_Serial(fy);
            for (int k = 0; k < N; ++k) fy_data[k] = P.prec_f0[k];
            return -1;
        }
    }

    // ── Step 3: Probe dense columns for mobile species ──────────────────────
    // For each mobile index j, compute full column J[:,j] via one FD probe,
    // subtract the banded part already in T, giving U[:,j] (the correction).
    // Then solve T̂^{-1} U[:,j] in-place.
    for (int jj = 0; jj < r; ++jj) {
        int j = P.prec_mobile_idx[jj];
        double dy = fd_delta(ydata[j], P.atol);
        ydata[j] += dy;

        ud->rhs_fn(t, y, fy, user_data);
        const double* f_tmp = N_VGetArrayPointer_Serial(fy);

        double inv_dy = 1.0 / dy;
        double* col = &P.prec_Tinv_U[jj * N];  // column jj of Tinv_U

        for (int k = 0; k < N; ++k) {
            double Jkj = (f_tmp[k] - P.prec_f0[k]) * inv_dy;
            // U[:,jj] = -γ * (J[k,j] - T[k,j]) for rows outside the band.
            // Within the band, T already captures J, so the correction is zero.
            int dist = std::abs(k - j);
            if (dist <= bw) {
                col[k] = 0.0;
            } else {
                col[k] = -gamma * Jkj;
            }
        }

        ydata[j] = P.prec_y_save[j];

        // Solve T̂^{-1} col in-place via dgbtrs
        {
            int info = 0;
            int kl = bw, ku = bw;
            int nrhs = 1;
            char trans = 'N';
            dgbtrs_(&trans, &N, &kl, &ku, &nrhs, P.prec_band.data(), &ldab,
                    P.prec_ipiv_band.data(), col, &N, &info);
            if (info != 0) {
                std::cerr << "[prec_setup_woodbury] dgbtrs for col " << jj
                          << " failed, info=" << info << "\n";
                double* fy_data = N_VGetArrayPointer_Serial(fy);
                for (int k = 0; k < N; ++k) fy_data[k] = P.prec_f0[k];
                return -1;
            }
        }
    }

    // ── Step 4: Form and factor Schur complement S = I + V^T T̂^{-1} U ──────
    // V = [e_{j1}, ..., e_{jr}], so V^T picks rows j1..jr from T̂^{-1} U.
    // S[ii][jj] = δ_{ii,jj} + (T̂^{-1} U)[mobile_idx[ii], jj]
    for (int ii = 0; ii < r; ++ii) {
        int row_idx = P.prec_mobile_idx[ii];
        for (int jj = 0; jj < r; ++jj) {
            P.prec_schur[ii + jj * r] = (ii == jj ? 1.0 : 0.0)
                                        + P.prec_Tinv_U[jj * N + row_idx];
        }
    }

    // Factor S via LAPACK dgetrf
    {
        int info = 0;
        dgetrf_(&r, &r, P.prec_schur.data(), &r,
                P.prec_ipiv_schur.data(), &info);
        if (info != 0) {
            std::cerr << "[prec_setup_woodbury] dgetrf (Schur) failed, info="
                      << info << "\n";
            double* fy_data = N_VGetArrayPointer_Serial(fy);
            for (int k = 0; k < N; ++k) fy_data[k] = P.prec_f0[k];
            return -1;
        }
    }

    // ── Restore fy to original f(t,y) ────────────────────────────────────────
    // CVODE expects fy to remain unchanged after prec_setup returns.
    double* fy_data = N_VGetArrayPointer_Serial(fy);
    for (int k = 0; k < N; ++k) fy_data[k] = P.prec_f0[k];

    return 0;
}

static int prec_solve_woodbury(sunrealtype t, N_Vector y, N_Vector fy,
                                N_Vector rv, N_Vector zv,
                                sunrealtype gamma, sunrealtype /*delta*/,
                                int /*lr*/, void* user_data) {
    const UserData*   ud = static_cast<const UserData*>(user_data);
    const Parameters& P  = *ud->P;
    const int N = P.N_eq;
    const int bw = P.prec_bw;
    const int r  = P.prec_rank;
    const int ldab = 3 * bw + 1;

    double* rdata = N_VGetArrayPointer_Serial(rv);
    double* zdata = N_VGetArrayPointer_Serial(zv);

    (void)t; (void)y; (void)fy; (void)gamma;

    // Step 1: w1 = T̂^{-1} r  (banded back-solve)
    for (int k = 0; k < N; ++k) zdata[k] = rdata[k];
    {
        int info = 0;
        int kl = bw, ku = bw;
        int nrhs = 1;
        char trans = 'N';
        dgbtrs_(&trans, &N, &kl, &ku, &nrhs, P.prec_band.data(), &ldab,
                P.prec_ipiv_band.data(), zdata, &N, &info);
        if (info != 0) {
            std::cerr << "[prec_solve_woodbury] dgbtrs failed, info=" << info << "\n";
            return -1;
        }
    }

    if (r == 0) return 0;  // No dense border — pure banded solve

    // Step 2: w2 = V^T w1  (extract r components at mobile indices)
    std::vector<double> w2(r);
    for (int jj = 0; jj < r; ++jj)
        w2[jj] = zdata[P.prec_mobile_idx[jj]];

    // Step 3: w3 = S^{-1} w2  (dense solve)
    {
        int info = 0;
        int nrhs = 1;
        char trans = 'N';
        dgetrs_(&trans, &r, &nrhs, P.prec_schur.data(), &r,
                P.prec_ipiv_schur.data(), w2.data(), &r, &info);
        if (info != 0) {
            std::cerr << "[prec_solve_woodbury] dgetrs (Schur) failed, info="
                      << info << "\n";
            return -1;
        }
    }

    // Step 4: z = w1 - (T̂^{-1} U) · w3
    for (int k = 0; k < N; ++k) {
        double corr = 0.0;
        for (int jj = 0; jj < r; ++jj)
            corr += P.prec_Tinv_U[jj * N + k] * w2[jj];
        zdata[k] -= corr;
    }

    return 0;
}

#endif  // CD_HAVE_LAPACK

// ── Dispatch wrappers ────────────────────────────────────────────────────────

int prec_setup(sunrealtype t, N_Vector y, N_Vector fy,
               sunbooleantype jok, sunbooleantype* jcurPtr,
               sunrealtype gamma, void* user_data) {
    const UserData* ud = static_cast<const UserData*>(user_data);
#ifdef CD_HAVE_LAPACK
    if (ud->P->prec_type == 1)
        return prec_setup_woodbury(t, y, fy, jok, jcurPtr, gamma, user_data);
#else
    if (ud->P->prec_type == 1) {
        static bool warned = false;
        if (!warned) {
            std::cerr << "[prec_setup] Woodbury requested but binary was built "
                         "without LAPACK -- falling back to Jacobi.\n";
            warned = true;
        }
    }
#endif
    return prec_setup_jacobi(t, y, fy, jok, jcurPtr, gamma, user_data);
}

int prec_solve(sunrealtype t, N_Vector y, N_Vector fy,
               N_Vector r, N_Vector z,
               sunrealtype gamma, sunrealtype delta,
               int lr, void* user_data) {
    const UserData* ud = static_cast<const UserData*>(user_data);
#ifdef CD_HAVE_LAPACK
    if (ud->P->prec_type == 1)
        return prec_solve_woodbury(t, y, fy, r, z, gamma, delta, lr, user_data);
#endif
    return prec_solve_jacobi(t, y, fy, r, z, gamma, delta, lr, user_data);
}
