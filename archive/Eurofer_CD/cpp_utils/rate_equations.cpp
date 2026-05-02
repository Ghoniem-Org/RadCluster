/**
 * rate_equations.cpp – ODE right-hand side for Eurofer_CD cluster dynamics.
 *
 * Faithfully translates RateEquations._rhs_decoupled / _rhs_shared from
 *   py_utils/rate_equations.py
 *
 * Physics: Ghoniem (2024), bcc Fe / EUROFER97.
 * All rate-constant arrays arrive pre-computed from the Python side.
 * He-pressure correction (GVV_eff) requires O(Nv) exp calls per RHS call
 * (unavoidable, as it depends on C_He which changes with time).
 *
 * K_IclV cross-term (SIA cluster ↔ vacancy cluster recombination) is computed
 * using a separable decomposition — O(Ni + Nv) instead of O(Ni × Nv):
 *
 *   K_IclV[n-1,m-1] = K_IclV_ns[n-1] + K_IclV_ni[n-1] * m13[m-1]
 *
 *   Vac_recom[m-1] = A + m13[m-1] * B
 *     where A = Σ_{n≥2} K_IclV_ns[n-1] * Ci_n
 *           B = Σ_{n≥2} K_IclV_ni[n-1] * Ci_n
 *
 *   SIA_recom[n-1] = K_IclV_ns[n-1]*C + K_IclV_ni[n-1]*D
 *     where C = Σ_m Cv_m
 *           D = Σ_m m13[m-1] * Cv_m
 */
#include "rate_equations.h"

#ifdef CD_HAVE_OPENMP
#  include <omp.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

// Floor a concentration to avoid underflow
static inline double fl(double c, double floor_val) {
    return c > floor_val ? c : floor_val;
}

// ── Helper: compute per-m GVV_eff array in-place ─────────────────────────────
// GVV_eff[m] = GVV[m] * exp(max(-ell_m * dE / kBT, -100))
// where:
//   ell_m = min(KHeV[m] * C_He / beta_He, L_He_max)
//   ell_c = max(ell_m, 0.1)
//   dE    = delta_He * beta_He_exp / (m+1) * (ell_c/(m+1))^{beta_He_exp-1}
// If C_He ≈ 0 the correction vanishes (GVV_eff = GVV).
static void compute_GVV_eff(const Parameters& P, double C_He,
                             std::vector<double>& GVV_eff) {
    const int Nv = P.Nv;
    GVV_eff.resize(Nv);
    if (C_He < 1e-200) {
        for (int k = 0; k < Nv; ++k) GVV_eff[k] = P.GVV[k];
        return;
    }
    const double inv_bHe  = 1.0 / std::max(P.beta_He, 1e-200);
    const double dHe      = P.delta_He;
    const double bexp     = P.beta_He_exp;
    const double L_max    = P.L_He_max;
    const double inv_kBT  = 1.0 / P.kBT;

    for (int k = 0; k < Nv; ++k) {
        double m     = static_cast<double>(k + 1);
        double ell_m = std::min(P.KHeV[k] * C_He * inv_bHe, L_max);
        if (ell_m < 1e-10) {
            GVV_eff[k] = P.GVV[k];
            continue;
        }
        double ell_c    = std::max(ell_m, 0.1);
        double dE       = dHe * bexp / m * std::pow(ell_c / m, bexp - 1.0);
        double exp_arg  = std::max(-ell_m * dE * inv_kBT, -100.0);
        GVV_eff[k]      = P.GVV[k] * std::exp(exp_arg);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Full-system RHS: rhs_eurofer
// ─────────────────────────────────────────────────────────────────────────────

int rhs_eurofer(sunrealtype /*t*/, N_Vector y, N_Vector ydot, void* user_data) {
    const Parameters& P  = *static_cast<const Parameters*>(user_data);
    const int Ni          = P.Ni;
    const int Nv          = P.Nv;
    const double C_floor  = P.C_floor;

    // ── Extract concentrations ────────────────────────────────────────────────
    // SIA: y[0..Ni-1], Vacancy: y[Ni..Ni+Nv-1], He: y[Ni+Nv]
    std::vector<double> Ci(Ni), Cv(Nv);
    for (int k = 0; k < Ni; ++k) Ci[k] = fl(NV_Ith_S(y, k),      C_floor);
    for (int k = 0; k < Nv; ++k) Cv[k] = fl(NV_Ith_S(y, Ni + k), C_floor);
    double C_He = fl(NV_Ith_S(y, Ni + Nv), C_floor);

    const double Ci1 = Ci[0];
    const double Cv1 = Cv[0];

    // ── He-pressure correction to GVV ────────────────────────────────────────
    std::vector<double> GVV_eff;
    compute_GVV_eff(P, C_He, GVV_eff);

    // ── K_IclV cross-term partial sums ────────────────────────────────────────
    // A = Σ_{n=2}^{Ni} K_IclV_ns[n-1] * Ci_n   (n=1 excluded: K_IclV_ns[0]=0)
    // B = Σ_{n=2}^{Ni} K_IclV_ni[n-1] * Ci_n
    // C_sum = Σ_{m=1}^{Nv} Cv_m
    // D     = Σ_{m=1}^{Nv} m13[m-1] * Cv_m
    double A = 0.0, B = 0.0, C_sum = 0.0, D = 0.0;
    for (int k = 1; k < Ni; ++k) {  // k=1 → n=2
        A += P.K_IclV_ns[k] * Ci[k];
        B += P.K_IclV_ni[k] * Ci[k];
    }
    for (int k = 0; k < Nv; ++k) {
        C_sum += Cv[k];
        D     += P.m13[k] * Cv[k];
    }

    // ── dCi equations ─────────────────────────────────────────────────────────
    // --- Ci1 (n=1) ---
    {
        double dCi1 = P.Pr_SIA[0]
                      - P.KII[0] * Ci1 * Ci1         // Ci+Ci → C2i (nucleation)
                      + 2.0 * P.GII[1] * Ci[1]        // C2i → 2×Ci1
                      - P.KIV[0] * Cv1 * Ci1          // Ci1 annihilated by Cv1
                      - P.k2_disl_i * Ci1;             // dislocation sink

        // Ci1 absorbed by vacancy clusters Cv2..Cv_Nv
        for (int k = 1; k < Nv; ++k)
            dCi1 -= P.KVI[k] * Cv[k] * Ci1;

        // Ci1 absorbed by SIA clusters Ci2..Ci_Ni
        for (int k = 1; k < Ni; ++k)
            dCi1 -= P.KII[k] * Ci[k] * Ci1;

        // Emission from Ci3..Ci_Ni returns one SIA to pool
        for (int k = 2; k < Ni; ++k)
            dCi1 += P.GII[k] * Ci[k];

        NV_Ith_S(ydot, 0) = dCi1;
    }

    // --- Ci2 (n=2) ---
    if (Ni >= 2) {
        double SIA_recom_2 = P.K_IclV_ns[1] * C_sum + P.K_IclV_ni[1] * D;
        double dCi2 = P.Pr_SIA[1]
                      + 0.5 * P.KII[0] * Ci1 * Ci1    // Ci+Ci → Ci2
                      - (P.KII[1] * Ci1                // Ci2 captures Ci → Ci3
                         + P.KIV[1] * Cv1              // Ci2 shrinks by Cv
                         + P.GII[1]                    // Ci2 emits SIA → 2×Ci1
                         + P.k2_SIA[1]                 // dislocation sink
                         + SIA_recom_2) * Ci[1];        // K_IclV recombination
        if (Ni >= 3)
            dCi2 += P.GII[2] * Ci[2] + P.KIV[2] * Cv1 * Ci[2];
        NV_Ith_S(ydot, 1) = dCi2;
    }

    // --- Ci_{n} for n = 3..Ni ---
    for (int k = 2; k < Ni; ++k) {  // k = n-1
        double SIA_recom_n = P.K_IclV_ns[k] * C_sum + P.K_IclV_ni[k] * D;
        double dCin = P.Pr_SIA[k]
                      + P.KII[k-1] * Ci1 * Ci[k-1]   // Ci_{n-1} + Ci1 → Ci_n
                      - P.KII[k]   * Ci1 * Ci[k]      // Ci_n + Ci1 → Ci_{n+1}
                      - P.KIV[k]   * Cv1 * Ci[k]      // Ci_n shrinks by Cv
                      - P.GII[k]   * Ci[k]             // Ci_n emits SIA
                      - P.k2_SIA[k] * Ci[k]            // dislocation sink
                      - SIA_recom_n * Ci[k];            // K_IclV recombination
        if (k < Ni - 1)
            dCin += P.GII[k+1]   * Ci[k+1]
                 +  P.KIV[k+1]  * Cv1 * Ci[k+1];
        NV_Ith_S(ydot, k) = dCin;
    }

    // ── dCv equations ─────────────────────────────────────────────────────────
    // --- Cv1 (m=1) ---
    {
        double Vac_recom_1 = A + P.m13[0] * B;   // m=1, m13[0]=1^{1/3}=1
        double dCv1 = P.Pr_VAC[0]
                      - P.KVV[0] * Cv1 * Cv1       // Cv+Cv → Cv2
                      + 2.0 * GVV_eff[1] * Cv[1]   // Cv2 emits Cv → 2×Cv1
                      - P.KVI[0] * Ci1 * Cv1        // Cv1 annihilated by Ci1
                      - P.KHeV[0] * C_He * Cv1      // He trapping at monovacancy
                      + P.k2_disl_v * (P.Cv_eq - Cv1)  // dislocation sink/source
                      - Vac_recom_1 * Cv1;           // SIA-cluster recombination

        // Emission from Cv3..Cv_Nv
        for (int k = 2; k < Nv; ++k)
            dCv1 += GVV_eff[k] * Cv[k];
        // Cv1 absorbed by larger vacancy clusters
        for (int k = 2; k < Nv; ++k)
            dCv1 -= P.KVV[k] * Cv[k] * Cv1;
        // Cv1 annihilated by SIA clusters Ci2..Ci_Ni
        for (int k = 1; k < Ni; ++k)
            dCv1 -= P.KIV[k] * Ci[k] * Cv1;

        NV_Ith_S(ydot, Ni) = dCv1;
    }

    // --- Cv2 (m=2) ---
    if (Nv >= 2) {
        double Vac_recom_2 = A + P.m13[1] * B;
        double dCv2 = P.Pr_VAC[1]
                      + 0.5 * P.KVV[0] * Cv1 * Cv1    // Cv+Cv → Cv2
                      - (P.KVV[1] * Cv1                 // Cv2 captures Cv → Cv3
                         + P.KVI[1] * Ci1               // Ci1 annihilates Cv2
                         + GVV_eff[1]                   // Cv2 emits Cv
                         + Vac_recom_2) * Cv[1];         // SIA-cluster recombination
        if (Nv >= 3)
            dCv2 += GVV_eff[2] * Cv[2] + P.KVI[2] * Ci1 * Cv[2];
        dCv2 -= P.KHeV[1] * C_He * Cv[1];
        NV_Ith_S(ydot, Ni + 1) = dCv2;
    }

    // --- Cv_{m} for m = 3..Nv ---
    for (int k = 2; k < Nv; ++k) {  // k = m-1
        double Vac_recom_m = A + P.m13[k] * B;
        double dCvm = P.Pr_VAC[k]
                      + P.KVV[k-1] * Cv1 * Cv[k-1]  // Cv_{m-1} + Cv → Cv_m
                      - P.KVV[k]   * Cv1 * Cv[k]     // Cv_m + Cv → Cv_{m+1}
                      - P.KVI[k]   * Ci1 * Cv[k]     // Ci1 annihilates Cv_m
                      - GVV_eff[k] * Cv[k]            // Cv_m emits Cv
                      - P.KHeV[k]  * C_He * Cv[k]    // He capture
                      - Vac_recom_m * Cv[k];           // SIA-cluster recombination
        if (k < Nv - 1)
            dCvm += GVV_eff[k+1]  * Cv[k+1]
                 +  P.KVI[k+1] * Ci1 * Cv[k+1];
        NV_Ith_S(ydot, Ni + k) = dCvm;
    }

    // ── dC_He ─────────────────────────────────────────────────────────────────
    {
        double He_capture = 0.0;
        for (int k = 0; k < Nv; ++k)
            He_capture += P.KHeV[k] * Cv[k];
        He_capture *= C_He;
        NV_Ith_S(ydot, Ni + Nv) = P.G_He - He_capture - P.k2_disl_He * C_He;
    }

    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase I Jacobi preconditioner  (user_data → Parameters, prec_diag sized)
// ─────────────────────────────────────────────────────────────────────────────

int prec_setup_win1(sunrealtype /*t*/, N_Vector y, N_Vector /*fy*/,
                    sunbooleantype /*jok*/, sunbooleantype* jcurPtr,
                    sunrealtype gamma, void* user_data) {
    Parameters& P = *static_cast<Parameters*>(user_data);
    const int N = P.Ni + P.Nv + 1;   // active window size
    P.prec_diag.resize(N);

    // Approximate diagonal of the Jacobian: d(dydt_k)/d(y_k)
    // Use a simple diagonal approximation based on the dominant loss terms.
    double Ci1 = std::max(NV_Ith_S(y, 0),      P.C_floor);
    double Cv1 = std::max(NV_Ith_S(y, P.Ni),   P.C_floor);
    double C_He= std::max(NV_Ith_S(y, P.Ni + P.Nv), P.C_floor);

    // SIA clusters
    P.prec_diag[0] = -(2.0 * P.KII[0] * Ci1 + P.k2_disl_i);
    for (int k = 1; k < P.Ni; ++k)
        P.prec_diag[k] = -(P.KII[k] * Ci1 + P.KIV[k] * Cv1 + P.GII[k] + P.k2_SIA[k]);

    // Vacancy clusters
    P.prec_diag[P.Ni] = -(2.0 * P.KVV[0] * Cv1 + P.KVI[0] * Ci1 + P.k2_disl_v);
    for (int k = 1; k < P.Nv; ++k)
        P.prec_diag[P.Ni + k] = -(P.KVV[k] * Cv1 + P.KVI[k] * Ci1 + P.GVV[k]);

    // Free He
    P.prec_diag[P.Ni + P.Nv] = -P.k2_disl_He;

    // Apply: prec_diag[k] = 1 / (1 - gamma * J_kk)
    for (int k = 0; k < N; ++k) {
        double d = 1.0 - gamma * P.prec_diag[k];
        P.prec_diag[k] = (std::abs(d) > 1e-100) ? 1.0 / d : 1.0;
    }
    *jcurPtr = SUNTRUE;
    return 0;
}

int prec_solve_win1(sunrealtype /*t*/, N_Vector /*y*/, N_Vector /*fy*/,
                    N_Vector r, N_Vector z,
                    sunrealtype /*gamma*/, sunrealtype /*delta*/,
                    int /*lr*/, void* user_data) {
    const Parameters& P = *static_cast<const Parameters*>(user_data);
    const int N = static_cast<int>(P.prec_diag.size());
    for (int k = 0; k < N; ++k)
        NV_Ith_S(z, k) = P.prec_diag[k] * NV_Ith_S(r, k);
    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase II frozen-sum maintenance
// ─────────────────────────────────────────────────────────────────────────────

void recompute_frozen_sums(WindowData& W) {
    const Parameters& P = *W.P_full;
    W.frozen_KII_sum   = 0.0;
    W.frozen_K_IclV_A  = 0.0;
    W.frozen_K_IclV_B  = 0.0;
    W.frozen_GII_sum   = 0.0;
    // Frozen range: Ci3 .. Ci_{x_lo_i-1}  (Ci2 always explicit)
    for (int n = 3; n <= W.x_lo_i - 1; ++n) {
        double Ci_n = W.full_conc[n - 1];   // Ci_n is at full_conc[n-1] in SIA layout
        int k = n - 1;  // 0-indexed
        W.frozen_KII_sum  += P.KII[k]       * Ci_n;
        W.frozen_K_IclV_A += P.K_IclV_ns[k] * Ci_n;
        W.frozen_K_IclV_B += P.K_IclV_ni[k] * Ci_n;
        W.frozen_GII_sum  += P.GII[k]        * Ci_n;
    }
    // Ghost for lowest active cluster (Ci_{x_lo_i-1})
    if (W.x_lo_i >= 3)
        W.Ci_frozen_top = W.full_conc[W.x_lo_i - 2];
    else
        W.Ci_frozen_top = 0.0;
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase II RHS: rhs_window
// Active layout: y[0..Nv-1]=Cv, y[Nv]=C_He, y[Nv+1]=Ci1, y[Nv+2..]=Ci_active
// ─────────────────────────────────────────────────────────────────────────────

int rhs_window(sunrealtype /*t*/, N_Vector y, N_Vector ydot, void* user_data) {
    const WindowData& W  = *static_cast<const WindowData*>(user_data);
    const Parameters& P  = *W.P_full;
    const int Nv          = P.Nv;
    const int x_lo        = W.x_lo_i;   // lowest active SIA cluster (1-indexed)
    const int x_hi        = W.x_hi_i;   // highest active SIA cluster
    const int n_ci_win    = x_hi - x_lo + 1;
    const double C_floor  = P.C_floor;

    // ── Unpack state ─────────────────────────────────────────────────────────
    std::vector<double> Cv(Nv);
    for (int k = 0; k < Nv; ++k)
        Cv[k] = fl(NV_Ith_S(y, k), C_floor);
    double C_He = fl(NV_Ith_S(y, Nv),     C_floor);
    double Ci1  = fl(NV_Ith_S(y, Nv + 1), C_floor);

    // Active SIA window (0-indexed within active vector)
    // active[j] = Ci_{x_lo + j}  for j = 0..n_ci_win-1
    std::vector<double> Ci_win(n_ci_win);
    for (int j = 0; j < n_ci_win; ++j)
        Ci_win[j] = fl(NV_Ith_S(y, Nv + 2 + j), C_floor);

    const double Cv1 = Cv[0];

    // ── He-pressure GVV_eff ───────────────────────────────────────────────────
    std::vector<double> GVV_eff;
    compute_GVV_eff(P, C_He, GVV_eff);

    // ── K_IclV partial sums ───────────────────────────────────────────────────
    // Active contributions from the window [x_lo..x_hi], plus frozen sums.
    double A = W.frozen_K_IclV_A;
    double B = W.frozen_K_IclV_B;
    for (int j = 0; j < n_ci_win; ++j) {
        int k = x_lo - 1 + j;  // 0-indexed
        A += P.K_IclV_ns[k] * Ci_win[j];
        B += P.K_IclV_ni[k] * Ci_win[j];
    }
    // Also include Ci2 (always explicit when x_lo_i == 2)
    // When x_lo_i > 2 Ci2 is in the frozen range (already in W.frozen_K_IclV_A/B)
    if (x_lo > 2) {
        // Ci2 is frozen: contribution already in frozen_K_IclV_A/B via recompute
    }
    // n=1 is always excluded (K_IclV_ns[0] = 0 from Python)

    double C_sum = 0.0, D = 0.0;
    for (int k = 0; k < Nv; ++k) {
        C_sum += Cv[k];
        D     += P.m13[k] * Cv[k];
    }

    // ── dCi1/dt ───────────────────────────────────────────────────────────────
    {
        double dCi1 = P.Pr_SIA[0]
                      - P.KII[0] * Ci1 * Ci1
                      - P.KIV[0] * Cv1 * Ci1
                      - P.k2_disl_i * Ci1
                      + W.frozen_GII_sum;    // emission from frozen Ci3..Ci_{x_lo-1}

        // C2i contribution (always handled explicitly)
        dCi1 += 2.0 * P.GII[1] * (x_lo == 2 ? Ci_win[0] : W.full_conc[1]);

        // Ci1 absorbed by active SIA clusters
        double active_KII_sum = P.KII[x_lo - 1] * (x_lo == 2 ? Ci_win[0] : W.full_conc[1]);
        if (x_lo == 2) {
            for (int j = 0; j < n_ci_win; ++j)
                dCi1 -= P.KII[x_lo - 1 + j] * Ci_win[j] * Ci1;
        } else {
            // Ci2 is frozen
            dCi1 -= P.KII[1] * W.full_conc[1] * Ci1;
            for (int j = 0; j < n_ci_win; ++j)
                dCi1 -= P.KII[x_lo - 1 + j] * Ci_win[j] * Ci1;
        }
        dCi1 -= W.frozen_KII_sum * Ci1;

        // Emission from active SIA window
        for (int j = 0; j < n_ci_win; ++j)
            dCi1 += P.GII[x_lo - 1 + j] * Ci_win[j];

        // Ci1 absorbed by vacancy clusters Cv2..Cv_Nv
        for (int k = 1; k < Nv; ++k)
            dCi1 -= P.KVI[k] * Cv[k] * Ci1;

        NV_Ith_S(ydot, Nv + 1) = dCi1;
    }

    // ── dCi for active window [x_lo..x_hi] ───────────────────────────────────
    for (int j = 0; j < n_ci_win; ++j) {
        int n = x_lo + j;      // cluster size (1-indexed)
        int k = n - 1;         // 0-indexed into arrays
        double Cin = Ci_win[j];

        double SIA_recom_n = P.K_IclV_ns[k] * C_sum + P.K_IclV_ni[k] * D;

        double dCin = P.Pr_SIA[k]
                      - P.KII[k]  * Ci1 * Cin
                      - P.KIV[k]  * Cv1 * Cin
                      - P.GII[k]  * Cin
                      - P.k2_SIA[k] * Cin
                      - SIA_recom_n * Cin;

        // Growth from n-1 → n:
        if (n == 2) {
            // Ci2 forms from Ci1+Ci1
            dCin += 0.5 * P.KII[0] * Ci1 * Ci1;
        } else {
            // Ci_{n-1} captures Ci1 → Ci_n
            double Ci_nm1 = (j > 0) ? Ci_win[j - 1] : W.Ci_frozen_top;
            dCin += P.KII[k - 1] * Ci1 * Ci_nm1;
        }

        // Shrink from n+1 → n (emission and vacancy absorption):
        if (j < n_ci_win - 1) {
            dCin += P.GII[k + 1]   * Ci_win[j + 1]
                 +  P.KIV[k + 1] * Cv1 * Ci_win[j + 1];
        }
        // At upper boundary (j == n_ci_win-1): Ci_{x_hi+1} = 0 (frozen above)

        NV_Ith_S(ydot, Nv + 2 + j) = dCin;
    }

    // ── dCv and dC_He (same as full system) ───────────────────────────────────
    // Vacancy clusters are always fully active.
    {
        double Vac_recom_1 = A + P.m13[0] * B;
        double dCv1 = P.Pr_VAC[0]
                      - P.KVV[0] * Cv1 * Cv1
                      + 2.0 * GVV_eff[1] * Cv[1]
                      - P.KVI[0] * Ci1 * Cv1
                      - P.KHeV[0] * C_He * Cv1
                      + P.k2_disl_v * (P.Cv_eq - Cv1)
                      - Vac_recom_1 * Cv1;

        for (int k = 2; k < Nv; ++k) dCv1 += GVV_eff[k] * Cv[k];
        for (int k = 2; k < Nv; ++k) dCv1 -= P.KVV[k] * Cv[k] * Cv1;

        // Cv1 annihilated by active SIA clusters
        if (x_lo == 2) {
            for (int j = 0; j < n_ci_win; ++j)
                dCv1 -= P.KIV[x_lo - 1 + j] * Ci_win[j] * Cv1;
        } else {
            dCv1 -= P.KIV[1] * W.full_conc[1] * Cv1;
            for (int j = 0; j < n_ci_win; ++j)
                dCv1 -= P.KIV[x_lo - 1 + j] * Ci_win[j] * Cv1;
        }
        NV_Ith_S(ydot, 0) = dCv1;
    }

    for (int k = 1; k < Nv; ++k) {
        double Vac_recom_m = A + P.m13[k] * B;
        double dCvm = P.Pr_VAC[k]
                      + P.KVV[k-1] * Cv1 * Cv[k-1]
                      - P.KVV[k]   * Cv1 * Cv[k]
                      - P.KVI[k]   * Ci1 * Cv[k]
                      - GVV_eff[k] * Cv[k]
                      - P.KHeV[k]  * C_He * Cv[k]
                      - Vac_recom_m * Cv[k];
        if (k < Nv - 1)
            dCvm += GVV_eff[k+1] * Cv[k+1] + P.KVI[k+1] * Ci1 * Cv[k+1];
        NV_Ith_S(ydot, k) = dCvm;
    }

    {
        double He_capture = 0.0;
        for (int k = 0; k < Nv; ++k) He_capture += P.KHeV[k] * Cv[k];
        NV_Ith_S(ydot, Nv) = P.G_He - He_capture * C_He - P.k2_disl_He * C_He;
    }

    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// Sherman-Morrison-Woodbury (SMW) preconditioner helpers
// ─────────────────────────────────────────────────────────────────────────────
//
// The K_IclV separable cross-coupling contributes a rank-4 off-diagonal block
// to the Jacobian.  Written in (I - γJ) form the correction is:
//
//   +γ * (u1⊗v1ᵀ + u2⊗v2ᵀ + u3⊗v3ᵀ + u4⊗v4ᵀ)
//
// where (using window-layout indices: 0..Nv-1 = vacancy, Nv+2+j = SIA window):
//   u1[Nv+2+j] = K_IclV_ns[n] * Ci_win[j]   v1[k∈vac] = 1
//   u2[Nv+2+j] = K_IclV_ni[n] * Ci_win[j]   v2[k∈vac] = m13[k]
//   u3[k∈vac]  = Cv[k]                        v3[Nv+2+j] = K_IclV_ns[n]
//   u4[k∈vac]  = Cv[k]*m13[k]                 v4[Nv+2+j] = K_IclV_ni[n]
//
// The preconditioner approximates M = D + γ·U·Vᵀ where D is diagonal.
// By Woodbury:  M⁻¹r = D⁻¹r − D⁻¹U · S⁻¹ · Vᵀ D⁻¹r
//               S = I₄ + Vᵀ D⁻¹ U  (4×4, inverted once per setup)
//
// Because u1,u2 and v3,v4 live in the SIA subspace while u3,u4 and v1,v2 live
// in the vacancy subspace, eight of the sixteen S entries are zero, leaving an
// anti-diagonal 2×2 block structure that is still inverted generically below.

// Invert a 4×4 matrix via Gauss-Jordan with partial pivoting.
// Returns false if singular (Ainv is untouched); true on success.
static bool invert_4x4(const double A[4][4],
                        std::array<std::array<double,4>,4>& Ainv)
{
    double M[4][8];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) M[i][j]   = A[i][j];
        for (int j = 0; j < 4; ++j) M[i][4+j] = (i == j) ? 1.0 : 0.0;
    }
    for (int col = 0; col < 4; ++col) {
        int piv = col;
        for (int row = col+1; row < 4; ++row)
            if (std::fabs(M[row][col]) > std::fabs(M[piv][col])) piv = row;
        if (std::fabs(M[piv][col]) < 1e-200) return false;
        if (piv != col)
            for (int j = 0; j < 8; ++j) std::swap(M[col][j], M[piv][j]);
        double sc = 1.0 / M[col][col];
        for (int j = 0; j < 8; ++j) M[col][j] *= sc;
        for (int row = 0; row < 4; ++row) {
            if (row == col) continue;
            double f = M[row][col];
            for (int j = 0; j < 8; ++j) M[row][j] -= f * M[col][j];
        }
    }
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            Ainv[i][j] = M[i][4+j];
    return true;
}

// Build the improved diagonal + SMW rank-4 correction for any WindowData-like
// struct W.  Called from prec_setup_window and prec_setup_window_omp.
//
// Improvements over the plain Jacobi preconditioner:
//   1. Diagonal entries now include the K_IclV contribution (SIA_recom_n and
//      Vac_recom_m), which the original code omitted.
//   2. The rank-4 off-diagonal K_IclV coupling is captured exactly via the
//      Woodbury formula, shrinking GMRES iteration counts from O(N) to O(10).
template<typename WType>
static void smw_setup(WType& W, N_Vector y, sunrealtype gamma,
                      sunbooleantype* jcurPtr)
{
    const Parameters& P  = *W.P_full;
    const int Nv         = P.Nv;
    const int x_lo       = W.x_lo_i;
    const int x_hi       = W.x_hi_i;
    const int n_ci_win   = x_hi - x_lo + 1;
    const int N          = W.N_active;
    const double C_floor = P.C_floor;

    // ── Extract key point-defect concentrations ──────────────────────────────
    double Cv1  = std::max(NV_Ith_S(y, 0),      C_floor);
    double C_He = std::max(NV_Ith_S(y, Nv),     C_floor);
    double Ci1  = std::max(NV_Ith_S(y, Nv + 1), C_floor);

    // ── K_IclV scalar sums needed for the improved diagonal ──────────────────
    // A_ki = Σ_{active n} K_IclV_ns[n-1]*Ci_n  (includes frozen contribution)
    // B_ki = Σ_{active n} K_IclV_ni[n-1]*Ci_n
    double A_ki = W.frozen_K_IclV_A;
    double B_ki = W.frozen_K_IclV_B;
    for (int j = 0; j < n_ci_win; ++j) {
        int k = x_lo - 1 + j;
        double Ci_n = std::max(NV_Ith_S(y, Nv + 2 + j), C_floor);
        A_ki += P.K_IclV_ns[k] * Ci_n;
        B_ki += P.K_IclV_ni[k] * Ci_n;
    }
    // C_sum = Σ_m Cv_m,  D_m13 = Σ_m m^{1/3}*Cv_m
    double C_sum = 0.0, D_m13 = 0.0;
    for (int k = 0; k < Nv; ++k) {
        double cv = std::max(NV_Ith_S(y, k), C_floor);
        C_sum += cv;
        D_m13 += P.m13[k] * cv;
    }

    // ── Build improved diagonal ──────────────────────────────────────────────
    W.prec_diag.resize(N);

    // Vacancy clusters: Cv1 (index 0), Cv2..Cv_Nv (index 1..Nv-1)
    W.prec_diag[0] = -(2.0 * P.KVV[0] * Cv1 + P.KVI[0] * Ci1 + P.k2_disl_v
                       + (A_ki + P.m13[0] * B_ki));
    for (int k = 1; k < Nv; ++k)
        W.prec_diag[k] = -(P.KVV[k] * Cv1 + P.KVI[k] * Ci1 + P.GVV[k]
                           + P.KHeV[k] * C_He + (A_ki + P.m13[k] * B_ki));

    // Free He (index Nv)
    W.prec_diag[Nv] = -P.k2_disl_He;

    // Ci1 (index Nv+1)
    W.prec_diag[Nv + 1] = -(2.0 * P.KII[0] * Ci1 + P.k2_disl_i);

    // SIA window (index Nv+2+j for j=0..n_ci_win-1)
    for (int j = 0; j < n_ci_win; ++j) {
        int k = x_lo - 1 + j;
        W.prec_diag[Nv + 2 + j] = -(P.KII[k] * Ci1 + P.KIV[k] * Cv1 + P.GII[k]
                                     + P.k2_SIA[k]
                                     + P.K_IclV_ns[k] * C_sum
                                     + P.K_IclV_ni[k] * D_m13);
    }

    // Apply Jacobi inversion: prec_diag[i] ← 1 / (1 − γ·J_ii)
    for (int i = 0; i < N; ++i) {
        double d = 1.0 - gamma * W.prec_diag[i];
        W.prec_diag[i] = (std::fabs(d) > 1e-100) ? 1.0 / d : 1.0;
    }

    // ── Build SMW rank-4 vectors: D⁻¹·(γ·u_k) ───────────────────────────────
    W.smw_DinvU1.resize(n_ci_win);
    W.smw_DinvU2.resize(n_ci_win);
    W.smw_DinvU3.resize(Nv);
    W.smw_DinvU4.resize(Nv);

    for (int j = 0; j < n_ci_win; ++j) {
        int k = x_lo - 1 + j;
        double Ci_n = std::max(NV_Ith_S(y, Nv + 2 + j), C_floor);
        W.smw_DinvU1[j] = W.prec_diag[Nv + 2 + j] * gamma * P.K_IclV_ns[k] * Ci_n;
        W.smw_DinvU2[j] = W.prec_diag[Nv + 2 + j] * gamma * P.K_IclV_ni[k] * Ci_n;
    }
    for (int k = 0; k < Nv; ++k) {
        double cv = std::max(NV_Ith_S(y, k), C_floor);
        W.smw_DinvU3[k] = W.prec_diag[k] * gamma * cv;
        W.smw_DinvU4[k] = W.prec_diag[k] * gamma * cv * P.m13[k];
    }

    // ── Build the 4×4 Woodbury matrix S = I₄ + Vᵀ D⁻¹ U ────────────────────
    // The v1/v2 rows (vacancy DOF) are orthogonal to the u1/u2 columns (SIA),
    // and v3/v4 rows (SIA) are orthogonal to u3/u4 columns (vacancy).
    // So only the two off-diagonal 2×2 blocks are nonzero.
    double S[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

    double v1u3=0, v1u4=0, v2u3=0, v2u4=0;
    for (int k = 0; k < Nv; ++k) {
        v1u3 += W.smw_DinvU3[k];
        v1u4 += W.smw_DinvU4[k];
        v2u3 += P.m13[k] * W.smw_DinvU3[k];
        v2u4 += P.m13[k] * W.smw_DinvU4[k];
    }
    S[0][2] += v1u3;  S[0][3] += v1u4;
    S[1][2] += v2u3;  S[1][3] += v2u4;

    double v3u1=0, v3u2=0, v4u1=0, v4u2=0;
    for (int j = 0; j < n_ci_win; ++j) {
        int k = x_lo - 1 + j;
        v3u1 += P.K_IclV_ns[k] * W.smw_DinvU1[j];
        v3u2 += P.K_IclV_ns[k] * W.smw_DinvU2[j];
        v4u1 += P.K_IclV_ni[k] * W.smw_DinvU1[j];
        v4u2 += P.K_IclV_ni[k] * W.smw_DinvU2[j];
    }
    S[2][0] += v3u1;  S[2][1] += v3u2;
    S[3][0] += v4u1;  S[3][1] += v4u2;

    if (!invert_4x4(S, W.smw_Sinv)) {
        // Singular: fall back to pure Jacobi (Sinv = 0 ⟹ no Woodbury correction)
        for (auto& row : W.smw_Sinv) row.fill(0.0);
    }

    *jcurPtr = SUNTRUE;
}

// Apply the SMW preconditioner solve: z = M⁻¹ r
template<typename WType>
static void smw_solve(const WType& W, N_Vector r, N_Vector z)
{
    const Parameters& P  = *W.P_full;
    const int Nv         = P.Nv;
    const int x_lo       = W.x_lo_i;
    const int x_hi       = W.x_hi_i;
    const int n_ci_win   = x_hi - x_lo + 1;
    const int N          = W.N_active;

    // Step 1: w = D⁻¹ r  (Jacobi apply)
    for (int i = 0; i < N; ++i)
        NV_Ith_S(z, i) = W.prec_diag[i] * NV_Ith_S(r, i);

    // Step 2: t = Vᵀ w  (4-vector; v1/v2 live in vacancy subspace, v3/v4 in SIA)
    double t[4] = {};
    for (int k = 0; k < Nv; ++k) {
        double wk = NV_Ith_S(z, k);
        t[0] += wk;
        t[1] += P.m13[k] * wk;
    }
    for (int j = 0; j < n_ci_win; ++j) {
        int  k  = x_lo - 1 + j;
        double wj = NV_Ith_S(z, Nv + 2 + j);
        t[2] += P.K_IclV_ns[k] * wj;
        t[3] += P.K_IclV_ni[k] * wj;
    }

    // Step 3: c = S⁻¹ t  (4×4 mat-vec)
    double c[4] = {};
    for (int i = 0; i < 4; ++i)
        for (int p = 0; p < 4; ++p)
            c[i] += W.smw_Sinv[i][p] * t[p];

    // Step 4: z -= D⁻¹ U c
    for (int j = 0; j < n_ci_win; ++j)
        NV_Ith_S(z, Nv + 2 + j) -= W.smw_DinvU1[j] * c[0]
                                  + W.smw_DinvU2[j] * c[1];
    for (int k = 0; k < Nv; ++k)
        NV_Ith_S(z, k) -= W.smw_DinvU3[k] * c[2]
                        + W.smw_DinvU4[k] * c[3];
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase II Jacobi preconditioner
// ─────────────────────────────────────────────────────────────────────────────

int prec_setup_window(sunrealtype /*t*/, N_Vector y, N_Vector /*fy*/,
                      sunbooleantype /*jok*/, sunbooleantype* jcurPtr,
                      sunrealtype gamma, void* user_data) {
    WindowData& W = *static_cast<WindowData*>(user_data);
    smw_setup(W, y, gamma, jcurPtr);
    return 0;
}

int prec_solve_window(sunrealtype /*t*/, N_Vector /*y*/, N_Vector /*fy*/,
                      N_Vector r, N_Vector z,
                      sunrealtype /*gamma*/, sunrealtype /*delta*/,
                      int /*lr*/, void* user_data) {
    const WindowData& W = *static_cast<const WindowData*>(user_data);
    smw_solve(W, r, z);
    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// Phase IV: OpenMP-parallelised RHS
// ─────────────────────────────────────────────────────────────────────────────

void recompute_frozen_sums_omp(WindowDataOMP& W) {
    // Same logic as recompute_frozen_sums
    const Parameters& P = *W.P_full;
    W.frozen_KII_sum   = 0.0;
    W.frozen_K_IclV_A  = 0.0;
    W.frozen_K_IclV_B  = 0.0;
    W.frozen_GII_sum   = 0.0;
    for (int n = 3; n <= W.x_lo_i - 1; ++n) {
        double Ci_n = W.full_conc[n - 1];
        int k = n - 1;
        W.frozen_KII_sum  += P.KII[k]        * Ci_n;
        W.frozen_K_IclV_A += P.K_IclV_ns[k]  * Ci_n;
        W.frozen_K_IclV_B += P.K_IclV_ni[k]  * Ci_n;
        W.frozen_GII_sum  += P.GII[k]         * Ci_n;
    }
    if (W.x_lo_i >= 3)
        W.Ci_frozen_top = W.full_conc[W.x_lo_i - 2];
    else
        W.Ci_frozen_top = 0.0;
}

int rhs_window_omp(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data) {
#ifdef CD_HAVE_OPENMP
    WindowDataOMP& W     = *static_cast<WindowDataOMP*>(user_data);
    const Parameters& P  = *W.P_full;
    const int Nv          = P.Nv;
    const int x_lo        = W.x_lo_i;
    const int x_hi        = W.x_hi_i;
    const int n_ci_win    = x_hi - x_lo + 1;
    const double C_floor  = P.C_floor;
    constexpr int OMP_MIN_WORK = 20000;

    // Unpack into pre-allocated buffers
    for (int k = 0; k < Nv; ++k)
        W.Cv_buf[k] = fl(NV_Ith_S(y, k), C_floor);
    double C_He = fl(NV_Ith_S(y, Nv),     C_floor);
    double Ci1  = fl(NV_Ith_S(y, Nv + 1), C_floor);
    for (int j = 0; j < n_ci_win; ++j)
        W.Ci_buf[j] = fl(NV_Ith_S(y, Nv + 2 + j), C_floor);

    if (n_ci_win < OMP_MIN_WORK) {
        // Serial path: reuse the Phase II logic (no parallel overhead for small windows)
        // Build a temporary WindowData from WindowDataOMP fields and call rhs_window
        WindowData Ws;
        Ws.P_full          = W.P_full;
        Ws.x_hi_i          = W.x_hi_i;
        Ws.x_lo_i          = W.x_lo_i;
        Ws.N_active        = W.N_active;
        Ws.full_conc       = W.full_conc;
        Ws.frozen_KII_sum  = W.frozen_KII_sum;
        Ws.frozen_K_IclV_A = W.frozen_K_IclV_A;
        Ws.frozen_K_IclV_B = W.frozen_K_IclV_B;
        Ws.frozen_GII_sum  = W.frozen_GII_sum;
        Ws.Ci_frozen_top   = W.Ci_frozen_top;
        return rhs_window(t, y, ydot, &Ws);
    }

    const double Cv1 = W.Cv_buf[0];

    // He-pressure GVV_eff
    std::vector<double> GVV_eff;
    compute_GVV_eff(P, C_He, GVV_eff);

    // K_IclV partial sums (parallel reduction for large windows)
    double A = W.frozen_K_IclV_A, B = W.frozen_K_IclV_B;
    double C_sum = 0.0, D = 0.0;

    int n_thr = (W.n_omp_threads > 0) ? W.n_omp_threads : 0;
    if (n_thr > 0) omp_set_num_threads(n_thr);

    #pragma omp parallel reduction(+:A,B,C_sum,D)
    {
        #pragma omp for schedule(static) nowait
        for (int j = 0; j < n_ci_win; ++j) {
            int k = x_lo - 1 + j;
            A += P.K_IclV_ns[k] * W.Ci_buf[j];
            B += P.K_IclV_ni[k] * W.Ci_buf[j];
        }
        #pragma omp for schedule(static)
        for (int k = 0; k < Nv; ++k) {
            C_sum += W.Cv_buf[k];
            D     += P.m13[k] * W.Cv_buf[k];
        }
    }

    // dCi1
    {
        double dCi1 = P.Pr_SIA[0]
                      - P.KII[0] * Ci1 * Ci1
                      - P.KIV[0] * Cv1 * Ci1
                      - P.k2_disl_i * Ci1
                      + W.frozen_GII_sum
                      + 2.0 * P.GII[1] * (x_lo == 2 ? W.Ci_buf[0] : W.full_conc[1]);
        dCi1 -= W.frozen_KII_sum * Ci1;
        if (x_lo > 2) dCi1 -= P.KII[1] * W.full_conc[1] * Ci1;
        for (int k = 1; k < Nv; ++k) dCi1 -= P.KVI[k] * W.Cv_buf[k] * Ci1;
        double kii_active = 0.0, gii_active = 0.0;
        #pragma omp parallel reduction(+:kii_active,gii_active)
        {
            #pragma omp for schedule(static)
            for (int j = 0; j < n_ci_win; ++j) {
                kii_active += P.KII[x_lo - 1 + j] * W.Ci_buf[j];
                gii_active += P.GII[x_lo - 1 + j] * W.Ci_buf[j];
            }
        }
        dCi1 -= kii_active * Ci1;
        dCi1 += gii_active;
        NV_Ith_S(ydot, Nv + 1) = dCi1;
    }

    // dCi window (parallel)
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < n_ci_win; ++j) {
        int n = x_lo + j;
        int k = n - 1;
        double Cin = W.Ci_buf[j];
        double SIA_recom = P.K_IclV_ns[k] * C_sum + P.K_IclV_ni[k] * D;
        double dCin = P.Pr_SIA[k]
                      - P.KII[k]   * Ci1 * Cin
                      - P.KIV[k]   * Cv1 * Cin
                      - P.GII[k]   * Cin
                      - P.k2_SIA[k]* Cin
                      - SIA_recom  * Cin;
        if (n == 2)
            dCin += 0.5 * P.KII[0] * Ci1 * Ci1;
        else {
            double Ci_nm1 = (j > 0) ? W.Ci_buf[j-1] : W.Ci_frozen_top;
            dCin += P.KII[k-1] * Ci1 * Ci_nm1;
        }
        if (j < n_ci_win - 1)
            dCin += P.GII[k+1] * W.Ci_buf[j+1] + P.KIV[k+1] * Cv1 * W.Ci_buf[j+1];
        NV_Ith_S(ydot, Nv + 2 + j) = dCin;
    }

    // dCv and dC_He (serial: Nv usually small)
    {
        double Vac_recom_1 = A + P.m13[0] * B;
        double dCv1 = P.Pr_VAC[0]
                      - P.KVV[0] * Cv1 * Cv1
                      + 2.0 * GVV_eff[1] * W.Cv_buf[1]
                      - P.KVI[0] * Ci1 * Cv1
                      - P.KHeV[0] * C_He * Cv1
                      + P.k2_disl_v * (P.Cv_eq - Cv1)
                      - Vac_recom_1 * Cv1;
        for (int k = 2; k < Nv; ++k) dCv1 += GVV_eff[k] * W.Cv_buf[k];
        for (int k = 2; k < Nv; ++k) dCv1 -= P.KVV[k] * W.Cv_buf[k] * Cv1;
        if (x_lo > 2) dCv1 -= P.KIV[1] * W.full_conc[1] * Cv1;
        for (int j = 0; j < n_ci_win; ++j)
            dCv1 -= P.KIV[x_lo - 1 + j] * W.Ci_buf[j] * Cv1;
        NV_Ith_S(ydot, 0) = dCv1;
    }
    for (int k = 1; k < Nv; ++k) {
        double Vac_recom_m = A + P.m13[k] * B;
        double dCvm = P.Pr_VAC[k]
                      + P.KVV[k-1] * Cv1 * W.Cv_buf[k-1]
                      - P.KVV[k]   * Cv1 * W.Cv_buf[k]
                      - P.KVI[k]   * Ci1 * W.Cv_buf[k]
                      - GVV_eff[k] * W.Cv_buf[k]
                      - P.KHeV[k]  * C_He * W.Cv_buf[k]
                      - Vac_recom_m * W.Cv_buf[k];
        if (k < Nv - 1)
            dCvm += GVV_eff[k+1] * W.Cv_buf[k+1] + P.KVI[k+1] * Ci1 * W.Cv_buf[k+1];
        NV_Ith_S(ydot, k) = dCvm;
    }
    {
        double He_cap = 0.0;
        for (int k = 0; k < Nv; ++k) He_cap += P.KHeV[k] * W.Cv_buf[k];
        NV_Ith_S(ydot, Nv) = P.G_He - He_cap * C_He - P.k2_disl_He * C_He;
    }
    return 0;
#else
    // OpenMP not available: fall back to Phase II
    WindowData Ws;
    const WindowDataOMP& W = *static_cast<const WindowDataOMP*>(user_data);
    Ws.P_full          = W.P_full;
    Ws.x_hi_i          = W.x_hi_i;
    Ws.x_lo_i          = W.x_lo_i;
    Ws.N_active        = W.N_active;
    Ws.full_conc       = W.full_conc;
    Ws.frozen_KII_sum  = W.frozen_KII_sum;
    Ws.frozen_K_IclV_A = W.frozen_K_IclV_A;
    Ws.frozen_K_IclV_B = W.frozen_K_IclV_B;
    Ws.frozen_GII_sum  = W.frozen_GII_sum;
    Ws.Ci_frozen_top   = W.Ci_frozen_top;
    return rhs_window(t, y, ydot, &Ws);
#endif
}

int prec_setup_window_omp(sunrealtype /*t*/, N_Vector y, N_Vector /*fy*/,
                           sunbooleantype /*jok*/, sunbooleantype* jcurPtr,
                           sunrealtype gamma, void* user_data) {
    WindowDataOMP& W = *static_cast<WindowDataOMP*>(user_data);
    smw_setup(W, y, gamma, jcurPtr);
    return 0;
}

int prec_solve_window_omp(sunrealtype /*t*/, N_Vector /*y*/, N_Vector /*fy*/,
                           N_Vector r, N_Vector z,
                           sunrealtype /*gamma*/, sunrealtype /*delta*/,
                           int /*lr*/, void* user_data) {
    const WindowDataOMP& W = *static_cast<const WindowDataOMP*>(user_data);
    smw_solve(W, r, z);
    return 0;
}
