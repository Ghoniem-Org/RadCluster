/**
 * rate_equations.cpp — ODE right-hand side for Expanded_Eurofer_CD.
 *
 * Implements:
 *   rhs_full_CD()    — full per-size master equations, Eqs. 152, 155, 157
 *   rhs_bin_moment() — size-bin moment equations, Chapter 9, Eqs. 193-208
 *
 * He-vacancy reduction (he_mode):
 *   he_mode == 0  — Case 2, fission/decoupled (Eq. 175)
 *   he_mode == 1  — Case 1, fusion/mean-field  (Eq. 174)
 *
 * Free He treatment (he_options):
 *   he_options == 0  — dynamic: c_h integrated as a full ODE (Eq. 157)
 *   he_options == 1  — quasi_steady_state: c_h computed algebraically from
 *                      dc_h/dt = 0 at each RHS call (E_m_h = 0.06 eV → fast)
 *
 * Concentration floor (C_floor):
 *   Enforced post-step in solver.cpp (after each successful CVode() call),
 *   NOT inside the RHS.  Clamping inside the RHS via max(y, C_floor) creates
 *   a Jacobian kink at y = C_floor — tested and confirmed to break CVODE's
 *   BDF corrector regardless of IC choice or tolerance.  The RHS uses
 *   max(y, 0) throughout: the kink is only at zero, where CVodeSetConstraints
 *   prevents the solver from probing, keeping the effective Jacobian smooth.
 *
 * He-balance corrections vs. pre-2026-04 version:
 *   • Removed erroneous KHeV term from dcv (He capture does not change the
 *     vacancy cluster size class m; marginal Σ_ℓ c_{m,ℓ} is unchanged).
 *   • Added missing k2_disl_v · Q_tot sink in the Q_tot equation (Case 2)
 *     and k2_disl_v · Q_m[m] sink in Q_m equations (Case 1).
 *
 * Physics reference: Ghoniem (2026), Rate_Equations.pdf.
 */

#include "rate_equations.h"
#include <nvector/nvector_serial.h>
#include <cmath>
#include <algorithm>

#ifdef CD_HAVE_OPENMP
#include <omp.h>
#endif

// ── Helpers ───────────────────────────────────────────────────────────────────

static inline double K_1D_eff(const Parameters& P, int n_idx, int m_idx) {
    // Mixed 1D/3D effective rate (Eq. 141)
    // K_n,m^eff = K_1D_pref[n] · m^{1/3} / (1 + B_rot · L̂² · m^{-1/3})
    int n = n_idx + 1;
    int m = m_idx + 1;
    if (n < 4 || n > P.n_max_i) return 0.0;
    double k_pref = P.K_1D_pref[n_idx];
    if (k_pref < 1e-300) return 0.0;
    double m_f   = static_cast<double>(m);
    double m13   = std::cbrt(m_f);
    double denom = 1.0 + P.B_rot * P.L_hat * P.L_hat / m13;
    return k_pref * m13 / denom;
}

/**
 * QSS free He concentration from dc_h/dt = 0 (Case 2):
 *   c_h = (G_He + beta_He · Q_tot) / (Σ KHeV[m] · c_v[m] + k2_He)
 *
 * All c_v values are clamped to 0 before summing to avoid amplifying
 * transiently negative concentrations (CVODE's Newton solver can probe
 * slightly negative intermediate states).
 */
static inline double c_h_qss_case2(const Parameters& P,
                                    const double* c_v, double Q_tot) {
    double sink = P.k2_disl_He;
    for (int m = 0; m < P.M; ++m)
        sink += P.KHeV[m] * std::max(c_v[m], 0.0);
    double source = P.G_He + P.beta_He * std::max(Q_tot, 0.0);
    return source / (sink > 1e-300 ? sink : 1e-300);
}

/**
 * QSS free He concentration from dc_h/dt = 0 (Case 1):
 *   c_h = (G_He + beta_He · Σ Q_m) / (Σ KHeV[m] · c_v[m] + k2_He)
 */
static inline double c_h_qss_case1(const Parameters& P,
                                    const double* c_v, const double* Q_m) {
    double sink    = P.k2_disl_He;
    double He_emit = 0.0;
    for (int m = 0; m < P.M; ++m) {
        sink    += P.KHeV[m] * std::max(c_v[m], 0.0);
        He_emit += P.beta_He * std::max(Q_m[m], 0.0);
    }
    double source = P.G_He + He_emit;
    return source / (sink > 1e-300 ? sink : 1e-300);
}

// ── Case 2 — fission/decoupled (Eq. 175) ─────────────────────────────────────
//
// State (dynamic):           [c_i(N) | c_v(M) | Q_tot | c_h]  N_eq = N+M+2
// State (quasi_steady_state):[c_i(N) | c_v(M) | Q_tot]        N_eq = N+M+1
//
// RHS design notes:
//   • All state variables clamped to max(y, 0) at access time.
//     This prevents negative amplification without introducing kinks above 0
//     (which would break CVODE's Jacobian estimation).
//   • C_floor is enforced post-step in solver.cpp, NOT inside the RHS.
//   • KHeV term removed from dcv: He capture preserves marginal Σ_ℓ c_{m,ℓ}.
//   • Q_tot equation includes -k2_disl_v · Q_tot sink.

static int rhs_case2(sunrealtype /*t*/, N_Vector yv, N_Vector ydotv,
                      const Parameters& P) {
    const double* y    = N_VGetArrayPointer_Serial(yv);
    double*       dydt = N_VGetArrayPointer_Serial(ydotv);
    const int  N   = P.N;
    const int  M   = P.M;
    const bool qss = (P.he_options == 1);

    // Unpack and clamp to 0 (CVODE can probe transiently negative values)
    const double* c_i = y;        // [N] — use max(y[n], 0) at access time
    const double* c_v = y + N;    // [M]
    const double  Q_tot = std::max(y[N + M], 0.0);

    // Free He: QSS algebraic or ODE state
    const double c_h = qss ? c_h_qss_case2(P, c_v, Q_tot)
                           : std::max(y[N + M + 1], 0.0);

    const double ci1 = std::max(c_i[0], 0.0);
    const double cv1 = std::max(c_v[0], 0.0);

    for (int k = 0; k < P.N_eq; ++k) dydt[k] = 0.0;

    double* dci = dydt;        // [N]
    double* dcv = dydt + N;    // [M]

    // Mean He per void (scalar, Eq. 175) — for He-pressure correction of GVV
    double C_vac_tot = 0.0;
    for (int m = 0; m < M; ++m) C_vac_tot += std::max(c_v[m], 0.0);
    const double ell_bar = (C_vac_tot > 1e-300) ? Q_tot / C_vac_tot : 0.0;

    // GVV corrected for He pressure (simplified expression)
    auto GVV_eff = [&](int m_idx) -> double {
        const int    m     = m_idx + 1;
        const double ell_m = ell_bar * std::pow(static_cast<double>(m), 2.0/3.0);
        if (ell_m < 1e-6) return P.GVV[m_idx];
        const double ratio = ell_m / static_cast<double>(m);
        if (ratio < 1e-10) return P.GVV[m_idx];
        const double dE = P.delta_He * P.beta_He_exp / static_cast<double>(m)
                          * std::pow(ratio, P.beta_He_exp - 1.0);
        return P.GVV[m_idx] * std::exp(std::min(-ell_m * dE / P.kBT, 0.0));
    };

    // ── SIA cluster equations (Eq. 152) ──────────────────────────────────────
    for (int n = 0; n < N; ++n) {
        const double cn = std::max(c_i[n], 0.0);

        dci[n] += P.Pr_SIA[n];

        // Thermal SIA emission: gain from n+1 → n
        if (n + 1 < N) dci[n] += P.GII[n + 1] * std::max(c_i[n + 1], 0.0);
        dci[n] -= P.GII[n] * cn;

        // SIA capture: gain from n-1 → n
        if (n > 0) dci[n] += P.KII[n - 1] * ci1 * std::max(c_i[n - 1], 0.0);
        dci[n] -= P.KII[n] * ci1 * cn;

        // Vacancy capture (SIA loop shrink): gain from n+1 → n
        if (n + 1 < N) dci[n] += P.KIV[n + 1] * cv1 * std::max(c_i[n + 1], 0.0);
        // Loss (n=0: monomer annihilates with vacancy)
        dci[n] -= P.KIV[n] * cv1 * cn;

        // 1D glide recombination (Eq. 141): n ≥ 4 (n_idx ≥ 3)
        if (n >= 3 && n < P.n_max_i) {
            for (int m = 0; m < M; ++m)
                dci[n] -= K_1D_eff(P, n, m) * cn * std::max(c_v[m], 0.0);
        }

        dci[n] -= P.k2_SIA[n] * cn;
    }

    // ── Vacancy cluster equations (Eq. 155, Case 2) ──────────────────────────
    // He capture does NOT change void size class m (only He occupancy ℓ changes),
    // so marginal c_m = Σ_ℓ c_{m,ℓ} is conserved by He capture.
    // No KHeV term here — He balance is handled entirely by Q_tot below.
    for (int m = 0; m < M; ++m) {
        const double cm    = std::max(c_v[m], 0.0);
        const double gvv_m = GVV_eff(m);

        dcv[m] += P.Pr_VAC[m];

        if (m + 1 < M) dcv[m] += GVV_eff(m + 1) * std::max(c_v[m + 1], 0.0);
        dcv[m] -= gvv_m * cm;

        if (m > 0) dcv[m] += P.KVV[m - 1] * cv1 * std::max(c_v[m - 1], 0.0);
        dcv[m] -= P.KVV[m] * cv1 * cm;

        if (m + 1 < M) dcv[m] += P.KVI[m + 1] * ci1 * std::max(c_v[m + 1], 0.0);
        if (m > 0)     dcv[m] -= P.KVI[m]     * ci1 * cm;

        for (int n = 3; n < std::min(N, P.n_max_i); ++n)
            dcv[m] -= K_1D_eff(P, n, m) * std::max(c_i[n], 0.0) * cm;

        dcv[m] -= P.k2_disl_v * cm;
    }

    // ── Q_tot equation (total He in voids) ───────────────────────────────────
    double He_uptake = 0.0;
    for (int m = 0; m < M; ++m)
        He_uptake += P.KHeV[m] * c_h * std::max(c_v[m], 0.0);
    const double He_emit = P.beta_He * Q_tot;
    const double He_sink = P.k2_disl_v * Q_tot;   // He removed when voids reach sinks
    dydt[N + M] = He_uptake - He_emit - He_sink;

    // ── Free He (Eq. 157) — dynamic mode only ────────────────────────────────
    if (!qss)
        dydt[N + M + 1] = P.G_He - He_uptake - P.k2_disl_He * c_h + He_emit;

    return 0;
}

// ── Case 1 — fusion/mean-field (Eq. 174) ─────────────────────────────────────
//
// State (dynamic):            [c_i(N) | c_v(M) | Q_m(M) | c_h]  N_eq = N+2M+1
// State (quasi_steady_state): [c_i(N) | c_v(M) | Q_m(M)]        N_eq = N+2M

static int rhs_case1(sunrealtype /*t*/, N_Vector yv, N_Vector ydotv,
                      const Parameters& P) {
    const double* y    = N_VGetArrayPointer_Serial(yv);
    double*       dydt = N_VGetArrayPointer_Serial(ydotv);
    const int  N   = P.N;
    const int  M   = P.M;
    const bool qss = (P.he_options == 1);

    const double* c_i = y;
    const double* c_v = y + N;
    const double* Q_m = y + N + M;

    const double c_h = qss ? c_h_qss_case1(P, c_v, Q_m)
                           : std::max(y[N + 2*M], 0.0);

    const double ci1 = std::max(c_i[0], 0.0);
    const double cv1 = std::max(c_v[0], 0.0);

    for (int k = 0; k < P.N_eq; ++k) dydt[k] = 0.0;

    double* dci = dydt;
    double* dcv = dydt + N;
    double* dQ  = dydt + N + M;

    // Per-class He loading ℓ̄_m = Q_m / c_m (Eq. 174)
    auto ell_bar_m = [&](int m_idx) -> double {
        const double cm = std::max(c_v[m_idx], 1e-300);
        return std::max(Q_m[m_idx], 0.0) / cm;
    };

    auto GVV_eff_m = [&](int m_idx) -> double {
        const int    m   = m_idx + 1;
        const double ell = ell_bar_m(m_idx);
        if (ell < 1e-6) return P.GVV[m_idx];
        const double ratio = ell / static_cast<double>(m);
        if (ratio < 1e-10) return P.GVV[m_idx];
        const double dE = P.delta_He * P.beta_He_exp / static_cast<double>(m)
                          * std::pow(ratio, P.beta_He_exp - 1.0);
        return P.GVV[m_idx] * std::exp(std::min(-ell * dE / P.kBT, 0.0));
    };

    // ── SIA clusters (identical to Case 2) ───────────────────────────────────
    for (int n = 0; n < N; ++n) {
        const double cn = std::max(c_i[n], 0.0);
        dci[n] += P.Pr_SIA[n];
        if (n + 1 < N) dci[n] += P.GII[n+1] * std::max(c_i[n+1], 0.0);
        dci[n] -= P.GII[n] * cn;
        if (n > 0)     dci[n] += P.KII[n-1] * ci1 * std::max(c_i[n-1], 0.0);
        dci[n] -= P.KII[n] * ci1 * cn;
        if (n + 1 < N) dci[n] += P.KIV[n+1] * cv1 * std::max(c_i[n+1], 0.0);
        dci[n] -= P.KIV[n] * cv1 * cn;
        if (n >= 3 && n < P.n_max_i) {
            for (int m = 0; m < M; ++m)
                dci[n] -= K_1D_eff(P, n, m) * cn * std::max(c_v[m], 0.0);
        }
        dci[n] -= P.k2_SIA[n] * cn;
    }

    // ── Vacancy clusters — no KHeV term (see Case 2 comment) ─────────────────
    for (int m = 0; m < M; ++m) {
        const double cm  = std::max(c_v[m], 0.0);
        const double gvv = GVV_eff_m(m);
        dcv[m] += P.Pr_VAC[m];
        if (m + 1 < M) dcv[m] += GVV_eff_m(m+1) * std::max(c_v[m+1], 0.0);
        dcv[m] -= gvv * cm;
        if (m > 0)     dcv[m] += P.KVV[m-1] * cv1 * std::max(c_v[m-1], 0.0);
        dcv[m] -= P.KVV[m] * cv1 * cm;
        if (m + 1 < M) dcv[m] += P.KVI[m+1] * ci1 * std::max(c_v[m+1], 0.0);
        if (m > 0)     dcv[m] -= P.KVI[m]   * ci1 * cm;
        for (int n = 3; n < std::min(N, P.n_max_i); ++n)
            dcv[m] -= K_1D_eff(P, n, m) * std::max(c_i[n], 0.0) * cm;
        dcv[m] -= P.k2_disl_v * cm;
    }

    // ── Q_m equations (He content per void class, Eq. 174) ───────────────────
    double He_cap_total  = 0.0;
    double He_emit_total = 0.0;
    for (int m = 0; m < M; ++m) {
        const double cm       = std::max(c_v[m], 0.0);
        const double qm       = std::max(Q_m[m], 0.0);
        const double he_cap_m = P.KHeV[m] * c_h * cm;
        const double he_emit_m = P.beta_He * qm;
        dQ[m]        += he_cap_m;
        dQ[m]        -= he_emit_m;
        dQ[m]        -= P.k2_disl_v * qm;   // He lost when voids absorbed at sinks
        He_cap_total  += he_cap_m;
        He_emit_total += he_emit_m;
    }

    // ── Free He (Eq. 157) — dynamic mode only ────────────────────────────────
    if (!qss)
        dydt[N + 2*M] = P.G_He - He_cap_total - P.k2_disl_He * c_h + He_emit_total;

    return 0;
}

// ── Dispatch ─────────────────────────────────────────────────────────────────

int rhs_full_CD(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data) {
    const UserData*   ud = static_cast<const UserData*>(user_data);
    const Parameters& P  = *ud->P;
    if (P.he_mode == 1)
        return rhs_case1(t, y, ydot, P);
    else
        return rhs_case2(t, y, ydot, P);
}

// ── Size-bin moment RHS (Chapter 9, Eqs. 193-208) ────────────────────────────

int rhs_bin_moment(sunrealtype t, N_Vector yv, N_Vector ydotv, void* user_data) {
    const UserData*   ud = static_cast<const UserData*>(user_data);
    const Parameters& P  = *ud->P;

    const double* y    = N_VGetArrayPointer_Serial(yv);
    double*       dydt = N_VGetArrayPointer_Serial(ydotv);

    const int  K   = P.K_bins;
    const int  M   = P.M;
    const int  N   = P.N;
    const bool qss = (P.he_options == 1);

    for (int k = 0; k < P.N_eq; ++k) dydt[k] = 0.0;

    // Build bin edges
    std::vector<int> n_lo(K), n_hi(K);
    {
        int edge = P.n1_bin;
        for (int k = 0; k < K; ++k) {
            n_lo[k] = edge;
            int next = std::max(static_cast<int>(std::floor(edge * P.r_ratio)), edge + 1);
            n_hi[k]  = std::min(next, N + 1);
            edge      = n_hi[k];
        }
    }

    // Piecewise-constant c_n from bin moments
    std::vector<double> c_n(N, 0.0);
    for (int k = 0; k < K; ++k) {
        const double mu0_k = std::max(y[2*k], 0.0);
        const double bw    = static_cast<double>(n_hi[k] - n_lo[k]);
        const double val   = (bw > 0) ? mu0_k / bw : 0.0;
        for (int n = n_lo[k]; n < n_hi[k]; ++n)
            if (n - 1 >= 0 && n - 1 < N) c_n[n - 1] = val;
    }

    const int     i_VAC = 2 * K;
    const double* c_v   = y + i_VAC;
    int           i_He_idx = -1;
    double        Q_tot    = 0.0;

    if (P.he_mode == 1) {
        const double* Q_m_ptr = y + i_VAC + M;
        for (int m = 0; m < M; ++m) Q_tot += std::max(Q_m_ptr[m], 0.0);
        i_He_idx = qss ? -1 : i_VAC + 2*M;
    } else {
        Q_tot    = std::max(y[i_VAC + M], 0.0);
        i_He_idx = qss ? -1 : i_VAC + M + 1;
    }

    const double c_h = qss
        ? (P.he_mode == 1 ? c_h_qss_case1(P, c_v, y + i_VAC + M)
                          : c_h_qss_case2(P, c_v, Q_tot))
        : std::max(y[i_He_idx], 0.0);

    const double ci1 = c_n.empty() ? 0.0 : std::max(c_n[0], 0.0);
    const double cv1 = std::max(c_v[0], 0.0);

    // dc_n/dt per-size
    std::vector<double> dc_n(N, 0.0);
    for (int n_idx = 0; n_idx < N; ++n_idx) {
        const double cn = std::max(c_n[n_idx], 0.0);
        dc_n[n_idx] += P.Pr_SIA[n_idx];
        if (n_idx + 1 < N) dc_n[n_idx] += P.GII[n_idx+1] * std::max(c_n[n_idx+1], 0.0);
        dc_n[n_idx] -= P.GII[n_idx] * cn;
        if (n_idx > 0)     dc_n[n_idx] += P.KII[n_idx-1] * ci1 * std::max(c_n[n_idx-1], 0.0);
        dc_n[n_idx] -= P.KII[n_idx] * ci1 * cn;
        if (n_idx + 1 < N) dc_n[n_idx] += P.KIV[n_idx+1] * cv1 * std::max(c_n[n_idx+1], 0.0);
        dc_n[n_idx] -= P.KIV[n_idx] * cv1 * cn;
        if (n_idx >= 3 && n_idx < P.n_max_i) {
            for (int m = 0; m < M; ++m)
                dc_n[n_idx] -= K_1D_eff(P, n_idx, m) * cn * std::max(c_v[m], 0.0);
        }
        dc_n[n_idx] -= P.k2_SIA[n_idx] * cn;
    }

    // Project onto bin moments
    for (int k = 0; k < K; ++k) {
        double dmu0 = 0.0, dmu1 = 0.0;
        for (int n = n_lo[k]; n < n_hi[k]; ++n) {
            if (n - 1 >= 0 && n - 1 < N) {
                dmu0 += dc_n[n - 1];
                dmu1 += static_cast<double>(n) * dc_n[n - 1];
            }
        }
        dydt[2*k]     = dmu0;
        dydt[2*k + 1] = dmu1;
    }

    // Inter-bin upwind flux
    for (int k = 0; k < K - 1; ++k) {
        const int    n_edge  = n_hi[k];
        if (n_edge < 1 || n_edge > N) continue;
        const int    n_idx   = n_edge - 1;
        const double c_edge  = (n_idx < N) ? c_n[n_idx] : 0.0;
        const double flux_fw = (n_idx < N) ? P.KII[n_idx] * ci1 * c_edge : 0.0;
        const double flux_bk = (n_idx < N) ? P.GII[n_idx] * c_edge : 0.0;
        dydt[2*k]         -= flux_fw;
        dydt[2*k + 1]     -= static_cast<double>(n_edge) * flux_fw;
        dydt[2*(k+1)]     += flux_fw;
        dydt[2*(k+1) + 1] += static_cast<double>(n_edge) * flux_fw;
        dydt[2*k]         += flux_bk;
        dydt[2*k + 1]     += static_cast<double>(n_edge) * flux_bk;
        dydt[2*(k+1)]     -= flux_bk;
        dydt[2*(k+1) + 1] -= static_cast<double>(n_edge) * flux_bk;
    }

    // Vacancy equations (no KHeV term)
    double C_vac_tot = 0.0;
    for (int m = 0; m < M; ++m) C_vac_tot += std::max(c_v[m], 0.0);
    const double ell_bar = (C_vac_tot > 1e-300) ? Q_tot / C_vac_tot : 0.0;

    double* dcv = dydt + i_VAC;
    for (int m = 0; m < M; ++m) {
        const double cm    = std::max(c_v[m], 0.0);
        const double ell_m = ell_bar * std::pow(static_cast<double>(m+1), 2.0/3.0);
        const double ratio = (ell_m > 1e-6) ? ell_m / (m+1) : 0.0;
        const double dE    = (ratio > 1e-10) ? P.delta_He * P.beta_He_exp / (m+1)
                                               * std::pow(ratio, P.beta_He_exp - 1.0) : 0.0;
        const double gvv_eff = P.GVV[m] * std::exp(std::min(-ell_m * dE / P.kBT, 0.0));
        dcv[m] += P.Pr_VAC[m];
        if (m + 1 < M) dcv[m] += P.GVV[m+1] * std::max(c_v[m+1], 0.0);
        dcv[m] -= gvv_eff * cm;
        if (m > 0)     dcv[m] += P.KVV[m-1] * cv1 * std::max(c_v[m-1], 0.0);
        dcv[m] -= P.KVV[m] * cv1 * cm;
        if (m + 1 < M) dcv[m] += P.KVI[m+1] * ci1 * std::max(c_v[m+1], 0.0);
        if (m > 0)     dcv[m] -= P.KVI[m]   * ci1 * cm;
        dcv[m] -= P.k2_disl_v * cm;
    }

    // He equations
    const int i_Qtot = i_VAC + M;
    if (P.he_mode != 1) {
        double He_up = 0.0;
        for (int m = 0; m < M; ++m)
            He_up += P.KHeV[m] * c_h * std::max(c_v[m], 0.0);
        const double He_emit = P.beta_He * Q_tot;
        dydt[i_Qtot] = He_up - He_emit - P.k2_disl_v * Q_tot;
        if (!qss)
            dydt[i_He_idx] = P.G_He - He_up - P.k2_disl_He * c_h + He_emit;
    } else {
        const double* Q_m_ptr = y + i_VAC + M;
        double*       dQ      = dydt + i_VAC + M;
        double He_cap_total = 0.0, He_emit_total = 0.0;
        for (int m = 0; m < M; ++m) {
            const double cm        = std::max(c_v[m], 0.0);
            const double qm        = std::max(Q_m_ptr[m], 0.0);
            const double he_cap_m  = P.KHeV[m] * c_h * cm;
            const double he_emit_m = P.beta_He * qm;
            dQ[m]        += he_cap_m;
            dQ[m]        -= he_emit_m;
            dQ[m]        -= P.k2_disl_v * qm;
            He_cap_total  += he_cap_m;
            He_emit_total += he_emit_m;
        }
        if (!qss)
            dydt[i_He_idx] = P.G_He - He_cap_total - P.k2_disl_He * c_h + He_emit_total;
    }

    return 0;
}

// ── Jacobi preconditioner ─────────────────────────────────────────────────────

int prec_setup(sunrealtype /*t*/, N_Vector /*y*/, N_Vector /*fy*/,
               sunbooleantype /*jok*/, sunbooleantype* jcurPtr,
               sunrealtype /*gamma*/, void* user_data) {
    *jcurPtr = SUNTRUE;
    (void)user_data;
    return 0;
}

int prec_solve(sunrealtype /*t*/, N_Vector /*y*/, N_Vector /*fy*/,
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
