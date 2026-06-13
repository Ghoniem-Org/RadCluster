/**
 * rate_equations.cpp – ODE right-hand side implementation for cluster dynamics.
 *
 * Faithfully translates RateEquations._rhs_full() from
 *   py_utils/rate_equations.py
 * which in turn translates rhs_full() from ClusterDynamics/codes/CD.ipynb.
 *
 * Physics: Ghoniem & Cho (1979), 316 stainless steel at 450 °C.
 * All parameters arrive pre-computed from the Python side via Parameters.
 * No std::exp inside the RHS — only arithmetic.
 *
 * Nv and Ni are runtime values stored in Parameters (no compile-time constants).
 *
 * Phase IV (window_mode=4, Multithread-OpenMP):
 *   rhs_window_omp()        — hot loops parallelised with OpenMP
 *   recompute_frozen_sums_omp()
 *   prec_setup_window_omp() / prec_solve_window_omp()
 * Compiled only when CD_HAVE_OPENMP is defined (cmake finds libomp).
 * Guards (#ifdef CD_HAVE_OPENMP) keep the code valid without OpenMP too.
 */
#include "rate_equations.h"

#ifdef CD_HAVE_OPENMP
#  include <omp.h>
#endif

#include <algorithm>   // std::max
#include <cmath>
#include <vector>

// Floor a concentration to avoid underflow
static inline double fl(double c, double floor) {
    return c > floor ? c : floor;
}

int rhs_cd(sunrealtype /*t*/, N_Vector y, N_Vector ydot, void* user_data) {
    const Parameters& P = *static_cast<const Parameters*>(user_data);
    const int    Nv      = P.Nv;
    const int    Ni      = P.Ni;
    const double C_floor = P.C_floor;

    // ── Extract and floor vacancy concentrations ──────────────────────────────
    // Cv_arr[k] = concentration of (k+1)-vacancy cluster, k = 0..Nv-1
    std::vector<double> Cv_arr(Nv);
    for (int k = 0; k < Nv; ++k)
        Cv_arr[k] = fl(NV_Ith_S(y, k), C_floor);

    // ── Extract and floor interstitial concentrations ─────────────────────────
    // Ci_arr[k] = concentration of (k+1)-interstitial cluster, k = 0..Ni-1
    std::vector<double> Ci_arr(Ni);
    for (int k = 0; k < Ni; ++k)
        Ci_arr[k] = fl(NV_Ith_S(y, Nv + k), C_floor);

    // Convenience aliases for single defects and dimer
    const double Cv  = Cv_arr[0];  // Cv1
    const double Ci  = Ci_arr[0];  // Ci1
    const double C2v = Cv_arr[1];  // Cv2
    const double C2i = Ci_arr[1];  // Ci2

    // ── dCv1/dt ───────────────────────────────────────────────────────────────
    {
        double dCv = (P.P_prod
                      + P.KCI[1]*Ci*C2v
                      + (2.0*P.GCV[1] - P.KCV[1]*Cv)*C2v);

        // emission from Cv3..CvNv  and  capture by Cv3..CvNv
        for (int k = 2; k < Nv; ++k) {
            dCv += P.GCV[k] * Cv_arr[k];
            dCv -= P.KCV[k] * Cv_arr[k] * Cv;
        }
        // capture by Ci3..CiNi  (vacancy absorbed by interstitial loop)
        for (int k = 2; k < Ni; ++k)
            dCv -= P.KLV[k] * Ci_arr[k] * Cv;

        dCv += P.Z_v * P.rho_d * P.Dv * (P.Cv_eq - Cv);  // dislocation sink
        dCv -= P.alpha  * Cv * Ci;                          // recombination
        dCv -= P.KCV[0] * Cv * Cv;                          // Cv + Cv → Cv2
        dCv -= P.KLV[1] * Cv * C2i;                         // C2i absorbs Cv → Ci

        NV_Ith_S(ydot, 0) = dCv;
    }

    // ── dCv2/dt ───────────────────────────────────────────────────────────────
    {
        double dC2v = (0.5 * P.KCV[0] * Cv * Cv
                       + P.GCV[2] * Cv_arr[2]
                       + P.KCI[2] * Ci * Cv_arr[2]
                       + P.rho_d * P.D2v * (P.C2v_eq - C2v));
        dC2v -= (P.KCV[1]*Cv + P.KCI[1]*Ci + P.GCV[1]) * C2v;
        NV_Ith_S(ydot, 1) = dC2v;
    }

    // ── dCvx/dt  for x = 3..Nv (k = x-1 = 2..Nv-1) ─────────────────────────
    for (int k = 2; k < Nv; ++k) {
        // k = x-1, so x = k+1
        // Cxm1 = Cv_arr[k-1], Cx = Cv_arr[k]
        double Cxm1 = Cv_arr[k - 1];
        double Cx   = Cv_arr[k];

        double dCx = (P.KCV[k-1]*Cv*Cxm1       // Cv(x-1) + Cv → Cvx
                      - P.KCI[k]*Ci*Cx           // Cvx + Ci → Cv(x-1)
                      - P.KCV[k]*Cv*Cx           // Cvx + Cv → Cv(x+1)
                      - P.GCV[k]*Cx);            // Cvx emits Cv → Cv(x-1)

        if (k < Nv - 1) {
            // gain from Cv(x+1) → Cvx + Cv  (emission and annihilation by Ci)
            dCx += P.KCI[k+1]*Ci*Cv_arr[k+1] + P.GCV[k+1]*Cv_arr[k+1];
        }
        NV_Ith_S(ydot, k) = dCx;
    }

    // ── dCi1/dt ───────────────────────────────────────────────────────────────
    {
        double dCi = (P.P_prod
                      + P.KLV[1]*Cv*C2i          // C2i absorbs Cv → Ci
                      - P.K_nuc_i*Ci*Ci           // Ci + Ci → C2i
                      - P.alpha*Cv*Ci             // recombination
                      - P.KLI[1]*Ci*C2i           // C2i captures Ci → C3i
                      - P.KCI[1]*Ci*C2v           // C2v captures Ci
                      + 2.0*P.GLV[1]*C2i);        // C2i dissociates → 2×Ci1

        // capture by Ci3..CiNi (interstitial clusters grow by absorbing Ci)
        for (int k = 2; k < Ni; ++k)
            dCi -= P.KLI[k] * Ci_arr[k] * Ci;

        // capture by vacancy clusters Cv3..CvNv (Ci absorbed into Cv cluster)
        for (int k = 2; k < Nv; ++k)
            dCi -= P.KCV[k] * Cv_arr[k] * Ci;

        dCi -= P.Z_i * P.rho_d * P.Di * Ci;     // dislocation sink

        NV_Ith_S(ydot, Nv) = dCi;
    }

    // ── dCi2/dt ───────────────────────────────────────────────────────────────
    {
        // GLV[1] = gamma_lv(2): thermal SIA emission C2i → 2×Ci1 (dissociation)
        double dC2i = (0.5 * P.K_nuc_i * Ci * Ci      // Ci + Ci → C2i
                       + P.KLV[2]*Cv*Ci_arr[2]         // C3i absorbs Cv → C2i
                       - (P.KLI[1]*Ci + P.KLV[1]*Cv + P.GLV[1]) * C2i);
        NV_Ith_S(ydot, Nv + 1) = dC2i;
    }

    // ── dCix/dt  for x = 3..Ni (k = x-1 = 2..Ni-1) ──────────────────────────
    for (int k = 2; k < Ni; ++k) {
        // k = x-1, so x = k+1
        // Cxm1 = Ci_arr[k-1], Cx = Ci_arr[k]
        double Cxm1 = Ci_arr[k - 1];
        double Cx   = Ci_arr[k];

        double dCx = (P.KLI[k-1]*Ci*Cxm1         // Ci(x-1) + Ci → Cix
                      - (P.KLI[k]*Ci + P.KLV[k]*Cv)*Cx);  // losses from Cix

        if (k < Ni - 1)
            dCx += P.KLV[k+1]*Cv*Ci_arr[k+1];    // Ci(x+1) absorbs Cv → Cix

        NV_Ith_S(ydot, Nv + k) = dCx;
    }

    // ── Floor enforcement ─────────────────────────────────────────────────────
    // Zero out any derivative that would drive a floored concentration further negative.
    for (int k = 0; k < P.N_EQ; ++k) {
        if (NV_Ith_S(y, k) <= C_floor && NV_Ith_S(ydot, k) < 0.0)
            NV_Ith_S(ydot, k) = 0.0;
    }

    return 0;
}

// ── Phase II: recompute frozen correction sums ────────────────────────────────

void recompute_frozen_sums(WindowData& W) {
    const Parameters& P = *W.P_full;
    W.frozen_KLI_sum = 0.0;
    W.frozen_KLV_sum = 0.0;
    // Ci2 is always handled explicitly; sums cover Ci3 .. Ci_{x_lo_i-1}.
    for (int x = 3; x < W.x_lo_i; ++x) {
        const double Ci_x = W.full_conc[P.Nv + x - 1];  // full_conc[Nv+x-1] = Ci_x
        W.frozen_KLI_sum += P.KLI[x - 1] * Ci_x;         // KLI[x-1] for cluster size x
        W.frozen_KLV_sum += P.KLV[x - 1] * Ci_x;
    }
    // Left-neighbour ghost for the lowest active cluster Ci_{x_lo_i}.
    W.Ci_frozen_top = (W.x_lo_i >= 3) ? W.full_conc[P.Nv + W.x_lo_i - 2] : 0.0;
    //   full_conc[Nv + (x_lo_i-1) - 1] = full_conc[Nv + x_lo_i - 2] = Ci_{x_lo_i - 1}
}

// ── Phase II: sliding-window RHS ─────────────────────────────────────────────
//
// State vector (size N_active = x_hi_v + 1 + (x_hi_i - x_lo_i + 1)):
//   y[0 .. x_hi_v-1]          : Cv1 .. Cv_{x_hi_v}
//   y[x_hi_v]                  : Ci1  (always active)
//   y[x_hi_v+1 .. x_hi_v+n-1] : Ci_{x_lo_i} .. Ci_{x_hi_i}  (n = x_hi_i-x_lo_i+1)
//
// Frozen clusters Ci2 .. Ci_{x_lo_i-1} live in W.full_conc.
// Ci2 is always accessed explicitly (never via frozen_KLI/KLV_sum).

int rhs_window(sunrealtype /*t*/, N_Vector y, N_Vector ydot, void* user_data) {
    const WindowData& W   = *static_cast<const WindowData*>(user_data);
    const Parameters& P   = *W.P_full;
    const int  x_hi_v     = W.x_hi_v;
    const int  x_lo_i     = W.x_lo_i;   // >= 2
    const int  x_hi_i     = W.x_hi_i;
    const int  n_ci_win   = x_hi_i - x_lo_i + 1;
    const double C_floor  = P.C_floor;

    // ── Unpack vacancy concentrations (Cv1..Cv_{x_hi_v}) ─────────────────────
    std::vector<double> Cv_arr(x_hi_v);
    for (int k = 0; k < x_hi_v; ++k)
        Cv_arr[k] = fl(NV_Ith_S(y, k), C_floor);
    const double Cv  = Cv_arr[0];   // Cv1

    // ── Ci1 (always at index x_hi_v) ─────────────────────────────────────────
    const double Ci = fl(NV_Ith_S(y, x_hi_v), C_floor);

    // ── Active interstitial window Ci_{x_lo_i}..Ci_{x_hi_i} ─────────────────
    std::vector<double> Ci_win(n_ci_win);
    for (int j = 0; j < n_ci_win; ++j)
        Ci_win[j] = fl(NV_Ith_S(y, x_hi_v + 1 + j), C_floor);

    // C2i: active when x_lo_i==2 (Ci_win[0]), frozen otherwise (from full_conc).
    const double C2i = (x_lo_i == 2) ? Ci_win[0]
                                      : fl(W.full_conc[P.Nv + 1], C_floor);
    // Convenience: C2v
    const double C2v = (x_hi_v >= 2) ? Cv_arr[1] : 0.0;

    // Sink loop start index in Ci_win: skip j=0 (Ci_{x_lo_i}) when x_lo_i==2
    // because Ci2 is handled explicitly in dCv1 and dCi1.
    const int j_sink_start = (x_lo_i == 2) ? 1 : 0;

    // ── dCv1/dt ───────────────────────────────────────────────────────────────
    {
        double dCv = P.P_prod;

        if (x_hi_v >= 2)
            dCv += P.KCI[1]*Ci*C2v + (2.0*P.GCV[1] - P.KCV[1]*Cv)*C2v;

        for (int k = 2; k < x_hi_v; ++k) {
            dCv += P.GCV[k] * Cv_arr[k];
            dCv -= P.KCV[k] * Cv_arr[k] * Cv;
        }

        // Interstitial sinks: Ci2 explicit, frozen Ci3..Ci_{x_lo_i-1}, active window.
        dCv -= P.KLV[1] * Cv * C2i;                  // Ci2 absorbs Cv (always explicit)
        dCv -= W.frozen_KLV_sum * Cv;                 // Ci3..Ci_{x_lo_i-1} frozen
        for (int j = j_sink_start; j < n_ci_win; ++j)
            dCv -= P.KLV[x_lo_i + j - 1] * Ci_win[j] * Cv;  // KLV[x-1] for cluster x

        dCv += P.Z_v * P.rho_d * P.Dv * (P.Cv_eq - Cv);
        dCv -= P.alpha * Cv * Ci;
        dCv -= P.KCV[0] * Cv * Cv;

        NV_Ith_S(ydot, 0) = dCv;
    }

    // ── dCv2/dt ───────────────────────────────────────────────────────────────
    if (x_hi_v >= 2) {
        const double C3v = (x_hi_v >= 3) ? Cv_arr[2] : 0.0;
        double dC2v = (0.5 * P.KCV[0] * Cv * Cv
                       + P.GCV[2] * C3v
                       + P.KCI[2] * Ci * C3v
                       + P.rho_d * P.D2v * (P.C2v_eq - C2v));
        dC2v -= (P.KCV[1]*Cv + P.KCI[1]*Ci + P.GCV[1]) * C2v;
        NV_Ith_S(ydot, 1) = dC2v;
    }

    // ── dCvx/dt for x = 3..x_hi_v ────────────────────────────────────────────
    for (int k = 2; k < x_hi_v; ++k) {
        const double Cxm1 = Cv_arr[k - 1];
        const double Cx   = Cv_arr[k];
        double dCx = (P.KCV[k-1]*Cv*Cxm1
                      - P.KCI[k]*Ci*Cx
                      - P.KCV[k]*Cv*Cx
                      - P.GCV[k]*Cx);
        if (k < x_hi_v - 1)
            dCx += P.KCI[k+1]*Ci*Cv_arr[k+1] + P.GCV[k+1]*Cv_arr[k+1];
        NV_Ith_S(ydot, k) = dCx;
    }

    // ── dCi1/dt ───────────────────────────────────────────────────────────────
    {
        double dCi = P.P_prod;
        dCi += P.KLV[1]*Cv*C2i;             // C2i absorbs Cv → Ci1 (gain)
        dCi -= P.K_nuc_i*Ci*Ci;
        dCi -= P.alpha*Cv*Ci;
        dCi -= P.KLI[1]*Ci*C2i;             // C2i captures Ci → C3i (loss)
        dCi -= P.KCI[1]*Ci*C2v;             // C2v captures Ci
        dCi += 2.0*P.GLV[1]*C2i;            // C2i dissociates → 2×Ci1

        dCi -= W.frozen_KLI_sum * Ci;       // Ci3..Ci_{x_lo_i-1} frozen sinks
        for (int j = j_sink_start; j < n_ci_win; ++j)
            dCi -= P.KLI[x_lo_i + j - 1] * Ci_win[j] * Ci;  // KLI[x-1] for cluster x

        for (int k = 2; k < x_hi_v; ++k)   // vacancy cluster sinks
            dCi -= P.KCV[k] * Cv_arr[k] * Ci;

        dCi -= P.Z_i * P.rho_d * P.Di * Ci;

        NV_Ith_S(ydot, x_hi_v) = dCi;
    }

    // ── dCi2/dt (only when x_lo_i == 2, i.e. Ci2 is active) ─────────────────
    if (x_lo_i == 2) {
        const double C3i = (n_ci_win > 1) ? Ci_win[1] : 0.0;
        // GLV[1] = gamma_lv(2): thermal SIA emission C2i → 2×Ci1 (dissociation)
        const double dC2i = (0.5 * P.K_nuc_i * Ci * Ci
                             + P.KLV[2] * Cv * C3i
                             - (P.KLI[1]*Ci + P.KLV[1]*Cv + P.GLV[1]) * C2i);
        NV_Ith_S(ydot, x_hi_v + 1) = dC2i;
    }

    // ── dCix/dt for the active interstitial window ────────────────────────────
    // Loop over Ci_{max(x_lo_i,3)} .. Ci_{x_hi_i}  (j_start skips Ci2 above).
    // At j == j_sink_start (first entry):
    //   x_lo_i==2: j_start=1, cluster=Ci3, left neighbour=Ci2=Ci_win[0]  ✓
    //   x_lo_i>=3: j_start=0, cluster=Ci_{x_lo_i}, left neighbour=Ci_frozen_top ✓
    for (int j = j_sink_start; j < n_ci_win; ++j) {
        const int    x  = x_lo_i + j;      // cluster size (1-indexed)
        const int    k  = x - 1;            // 0-indexed array position
        const double Cx = Ci_win[j];

        double Cxm1;
        if (j == j_sink_start && x_lo_i >= 3)
            Cxm1 = W.Ci_frozen_top;         // Ci_{x_lo_i-1} is frozen
        else
            Cxm1 = Ci_win[j - 1];           // previous active cluster

        double dCx = (P.KLI[k-1]*Ci*Cxm1
                      - (P.KLI[k]*Ci + P.KLV[k]*Cv)*Cx);
        if (j < n_ci_win - 1)
            dCx += P.KLV[k+1]*Cv*Ci_win[j+1];

        NV_Ith_S(ydot, x_hi_v + 1 + j) = dCx;
    }

    // ── Floor enforcement ─────────────────────────────────────────────────────
    for (int k = 0; k < W.N_active; ++k) {
        if (NV_Ith_S(y, k) <= C_floor && NV_Ith_S(ydot, k) < 0.0)
            NV_Ith_S(ydot, k) = 0.0;
    }

    return 0;
}

// ── Jacobi diagonal preconditioner ───────────────────────────────────────────

int prec_setup_window(sunrealtype /*t*/, N_Vector y, N_Vector /*fy*/,
                      sunbooleantype /*jok*/, sunbooleantype* jcurPtr,
                      sunrealtype gamma, void* user_data) {
    WindowData& W        = *static_cast<WindowData*>(user_data);
    const Parameters& P  = *W.P_full;
    const int x_hi_v     = W.x_hi_v;
    const int x_lo_i     = W.x_lo_i;
    const int x_hi_i     = W.x_hi_i;
    const double C_floor = P.C_floor;

    const double Cv = fl(NV_Ith_S(y, 0),      C_floor);
    const double Ci = fl(NV_Ith_S(y, x_hi_v), C_floor);

    W.prec_diag.resize(W.N_active);

    // Vacancy block: diagonal ≈ -(KCV[k]*Cv + KCI[k]*Ci + GCV[k])
    for (int k = 0; k < x_hi_v; ++k)
        W.prec_diag[k] = 1.0 - gamma * (-(P.KCV[k]*Cv + P.KCI[k]*Ci + P.GCV[k]));

    // Ci1: diagonal dominated by loss terms
    {
        double J_ii = -(P.K_nuc_i*2.0*Ci + P.alpha*Cv + P.Z_i*P.rho_d*P.Di
                        + P.KLI[1]*fl(NV_Ith_S(y, x_hi_v+1), C_floor)   // C2i loss
                        + W.frozen_KLI_sum);
        const int n = x_hi_i - x_lo_i + 1;
        const int j0 = (x_lo_i == 2) ? 1 : 0;
        for (int j = j0; j < n; ++j)
            J_ii -= P.KLI[x_lo_i + j - 1] * fl(NV_Ith_S(y, x_hi_v + 1 + j), C_floor);
        W.prec_diag[x_hi_v] = 1.0 - gamma * J_ii;
    }

    // Active interstitial window: diagonal = -(KLI[x-1]*Ci + KLV[x-1]*Cv)
    const int n_ci_win = x_hi_i - x_lo_i + 1;
    for (int j = 0; j < n_ci_win; ++j) {
        const int x = x_lo_i + j;
        W.prec_diag[x_hi_v + 1 + j] =
            1.0 - gamma * (-(P.KLI[x-1]*Ci + P.KLV[x-1]*Cv));
    }

    *jcurPtr = SUNTRUE;
    return 0;
}

int prec_solve_window(sunrealtype /*t*/, N_Vector /*y*/, N_Vector /*fy*/,
                      N_Vector r, N_Vector z,
                      sunrealtype /*gamma*/, sunrealtype /*delta*/,
                      int /*lr*/, void* user_data) {
    const WindowData& W = *static_cast<const WindowData*>(user_data);
    for (int i = 0; i < W.N_active; ++i) {
        const double m = W.prec_diag[i];
        NV_Ith_S(z, i) = NV_Ith_S(r, i) / (std::abs(m) > 1e-300 ? m : 1.0);
    }
    return 0;
}

// ── Phase I Jacobi-diagonal preconditioner ───────────────────────────────────
//
// user_data is a Parameters* (the active P_win in Phase I, with Nv=x_hi_v, Ni=x_hi_i).
// The state vector layout is flat: y[0..Nv-1] = Cv, y[Nv..Nv+Ni-1] = Ci.
//
// Diagonal approximations (same approach as Phase II):
//   Cv block (k=0..Nv-1): diag ≈ 1 - γ*(-(KCV[k]*Cv1 + KCI[k]*Ci1 + GCV[k]))
//   Ci1      (k=Nv)     : diag ≈ 1 - γ*(-(KLI[0]*Ci1 + KLV[0]*Cv1))
//   Ci block (k=Nv+j, j=1..Ni-1): diag ≈ 1 - γ*(-(KLI[j]*Ci1 + KLV[j]*Cv1))

int prec_setup_win1(sunrealtype /*t*/, N_Vector y, N_Vector /*fy*/,
                    sunbooleantype /*jok*/, sunbooleantype* jcurPtr,
                    sunrealtype gamma, void* user_data) {
    Parameters& P       = *static_cast<Parameters*>(user_data);
    const int    Nv     = P.Nv;
    const int    Ni     = P.Ni;
    const double C_floor = P.C_floor;

    const double Cv = fl(NV_Ith_S(y, 0),  C_floor);  // Cv1
    const double Ci = fl(NV_Ith_S(y, Nv), C_floor);  // Ci1

    P.prec_diag.resize(static_cast<size_t>(Nv + Ni));

    // Vacancy cluster block
    for (int k = 0; k < Nv; ++k)
        P.prec_diag[k] = 1.0 - gamma * (-(P.KCV[k]*Cv + P.KCI[k]*Ci + P.GCV[k]));

    // Interstitial cluster block (Ci1 and larger)
    for (int j = 0; j < Ni; ++j)
        P.prec_diag[Nv + j] = 1.0 - gamma * (-(P.KLI[j]*Ci + P.KLV[j]*Cv));

    *jcurPtr = SUNTRUE;
    return 0;
}

int prec_solve_win1(sunrealtype /*t*/, N_Vector /*y*/, N_Vector /*fy*/,
                    N_Vector r, N_Vector z,
                    sunrealtype /*gamma*/, sunrealtype /*delta*/,
                    int /*lr*/, void* user_data) {
    const Parameters& P = *static_cast<const Parameters*>(user_data);
    const int n = P.Nv + P.Ni;
    for (int i = 0; i < n; ++i) {
        const double m = P.prec_diag[i];
        NV_Ith_S(z, i) = NV_Ith_S(r, i) / (std::abs(m) > 1e-300 ? m : 1.0);
    }
    return 0;
}

// ── Phase IV: Multithread-OpenMP ──────────────────────────────────────────────
//
// rhs_window_omp() mirrors rhs_window() exactly in physics but:
//   1. Uses pre-allocated W.Cv_buf / W.Ci_buf instead of local std::vector,
//      eliminating heap allocation from every RHS call.
//   2. Parallelises the dominant loops with OpenMP:
//        - Buffer fill (Cv, Ci)
//        - dCv1 KLV sink accumulation  (reduction)
//        - dCi1 KLI sink accumulation  (reduction)
//        - dCvx loop (independent writes)
//        - dCix active-window loop (independent writes)
//        - Floor enforcement
// All loops use schedule(static) for deterministic assignment and minimal
// synchronisation overhead on the M3 Max's 12 P-cores.

void recompute_frozen_sums_omp(WindowDataOMP& W) {
    const Parameters& P = *W.P_full;
    W.frozen_KLI_sum = 0.0;
    W.frozen_KLV_sum = 0.0;
    for (int x = 3; x < W.x_lo_i; ++x) {
        const double Ci_x = W.full_conc[P.Nv + x - 1];
        W.frozen_KLI_sum += P.KLI[x - 1] * Ci_x;
        W.frozen_KLV_sum += P.KLV[x - 1] * Ci_x;
    }
    W.Ci_frozen_top = (W.x_lo_i >= 3) ? W.full_conc[P.Nv + W.x_lo_i - 2] : 0.0;
}

int rhs_window_omp(sunrealtype /*t*/, N_Vector y, N_Vector ydot, void* user_data) {
    WindowDataOMP&    W        = *static_cast<WindowDataOMP*>(user_data);
    const Parameters& P        = *W.P_full;
    const int  x_hi_v          = W.x_hi_v;
    const int  x_lo_i          = W.x_lo_i;
    const int  x_hi_i          = W.x_hi_i;
    const int  n_ci_win        = x_hi_i - x_lo_i + 1;
    const double C_floor       = P.C_floor;
    const int  n_thr           = W.n_omp_threads;  // 0 = OMP default

    double* Cv_arr = W.Cv_buf.data();
    double* Ci_win = W.Ci_buf.data();

    // ── Minimum-work threshold ─────────────────────────────────────────────────
    // OpenMP fork-join overhead (~5–50 µs/call depending on OS/runtime) dominates
    // for small windows.  Below OMP_MIN_WORK iterations the parallel overhead
    // exceeds the compute benefit; fall through to serial execution while still
    // using the pre-allocated buffers (which eliminate per-call heap allocation).
    //
    // Threshold derivation (MSVC OpenMP on Xeon Silver 4116 @ 2.1 GHz, 12 cores):
    //   fork cost ≈ 10 µs; compute per iteration ≈ 6 FLOPs @ ~25 GFLOPS/thread
    //   → break-even: n_ci_win / 12 * 6 / 25e9 = 10e-6 → n_ci_win ≈ 500 000
    //   → 20 000 is conservative but safe; at W=1000 serial path always wins.
    // At W ≥ 20 000 (larger Ni problems) the parallel path activates and provides
    // additional speedup beyond the AVX2/SIMD gains from /arch:AVX2.

    static constexpr int OMP_MIN_WORK = 20000;

#ifdef CD_HAVE_OPENMP
    const int n_thr_eff = (n_thr > 0) ? n_thr : omp_get_max_threads();
    const bool use_omp  = (n_ci_win >= OMP_MIN_WORK);
#else
    const bool use_omp  = false;
#endif

    // =========================================================================
    // PARALLEL PATH — single persistent parallel region (one fork-join per call)
    // =========================================================================
    //
    // Thread structure:
    //   1. omp for (nowait) — fill Cv_arr
    //   2. omp for          — fill Ci_win  [barrier → both buffers complete]
    //   3. omp single       — set scalars Cv, Ci, C2i, C2v  [barrier]
    //   4. omp for          — combined sink_lv + sink_li reduction  [barrier]
    //   5. omp single nowait— scalar equations: dCv1, dCv2, dCi1, dCi2
    //                         (writes ydot[0], ydot[1], ydot[x_hi_v])
    //   6. omp for nowait   — dCvx  (writes ydot[2..x_hi_v-1])
    //   7. omp for nowait   — dCix  (writes ydot[x_hi_v+1..N_active-1])
    //   Steps 5–7 run concurrently; they write disjoint ydot entries → no race.
    //   8. omp for          — floor enforcement  [barrier → end of parallel]

    if (use_omp) {
#ifdef CD_HAVE_OPENMP
        // Shared scalars set inside the parallel region (after buffer barrier).
        double Cv = 0.0, Ci = 0.0, C2i = 0.0, C2v = 0.0;
        double sink_lv = 0.0, sink_li = 0.0;
        const int j_ss = (x_lo_i == 2) ? 1 : 0;  // j_sink_start (firstprivate below)

        #pragma omp parallel num_threads(n_thr_eff) \
            default(none) \
            shared(W, P, y, ydot, Cv_arr, Ci_win, Cv, Ci, C2i, C2v, sink_lv, sink_li) \
            firstprivate(x_hi_v, x_lo_i, x_hi_i, n_ci_win, C_floor, j_ss)
        {
            // 1 ── Fill Cv_arr (nowait: Ci fill will provide the needed barrier)
            #pragma omp for schedule(static) nowait
            for (int k = 0; k < x_hi_v; ++k)
                Cv_arr[k] = fl(NV_Ith_S(y, k), C_floor);

            // 2 ── Fill Ci_win [implicit barrier: both buffers ready after this]
            #pragma omp for schedule(static)
            for (int j = 0; j < n_ci_win; ++j)
                Ci_win[j] = fl(NV_Ith_S(y, x_hi_v + 1 + j), C_floor);

            // 3 ── Set scalars (single thread, implicit barrier → all threads wait)
            #pragma omp single
            {
                Cv  = Cv_arr[0];
                Ci  = fl(NV_Ith_S(y, x_hi_v), C_floor);
                C2i = (x_lo_i == 2) ? Ci_win[0]
                                     : fl(W.full_conc[P.Nv + 1], C_floor);
                C2v = (x_hi_v >= 2) ? Cv_arr[1] : 0.0;
            }

            // 4 ── Combined KLV + KLI reductions in one loop pass [barrier]
            #pragma omp for schedule(static) reduction(+:sink_lv) reduction(+:sink_li)
            for (int j = j_ss; j < n_ci_win; ++j) {
                sink_lv += P.KLV[x_lo_i + j - 1] * Ci_win[j];
                sink_li += P.KLI[x_lo_i + j - 1] * Ci_win[j];
            }

            // 5 ── Scalar equations (thread 0, nowait → others proceed to 6 & 7)
            #pragma omp single nowait
            {
                // dCv1
                double dCv = P.P_prod;
                if (x_hi_v >= 2)
                    dCv += P.KCI[1]*Ci*C2v + (2.0*P.GCV[1] - P.KCV[1]*Cv)*C2v;
                for (int k = 2; k < x_hi_v; ++k) {
                    dCv += P.GCV[k] * Cv_arr[k];
                    dCv -= P.KCV[k] * Cv_arr[k] * Cv;
                }
                dCv -= P.KLV[1]*Cv*C2i + W.frozen_KLV_sum*Cv + sink_lv*Cv;
                dCv += P.Z_v*P.rho_d*P.Dv*(P.Cv_eq - Cv);
                dCv -= P.alpha*Cv*Ci + P.KCV[0]*Cv*Cv;
                NV_Ith_S(ydot, 0) = dCv;

                // dCv2
                if (x_hi_v >= 2) {
                    const double C3v = (x_hi_v >= 3) ? Cv_arr[2] : 0.0;
                    double dC2v = (0.5*P.KCV[0]*Cv*Cv + P.GCV[2]*C3v
                                   + P.KCI[2]*Ci*C3v
                                   + P.rho_d*P.D2v*(P.C2v_eq - C2v));
                    dC2v -= (P.KCV[1]*Cv + P.KCI[1]*Ci + P.GCV[1]) * C2v;
                    NV_Ith_S(ydot, 1) = dC2v;
                }

                // dCi1
                double dCi = P.P_prod;
                dCi += P.KLV[1]*Cv*C2i;
                dCi -= P.K_nuc_i*Ci*Ci + P.alpha*Cv*Ci;
                dCi -= P.KLI[1]*Ci*C2i + P.KCI[1]*Ci*C2v;
                dCi += 2.0*P.GLV[1]*C2i;            // C2i dissociates → 2×Ci1
                for (int k = 2; k < x_hi_v; ++k)
                    dCi -= P.KCV[k] * Cv_arr[k] * Ci;
                dCi -= W.frozen_KLI_sum*Ci + sink_li*Ci;
                dCi -= P.Z_i*P.rho_d*P.Di*Ci;
                NV_Ith_S(ydot, x_hi_v) = dCi;

                // dCi2
                if (x_lo_i == 2) {
                    const double C3i = (n_ci_win > 1) ? Ci_win[1] : 0.0;
                    // GLV[1]: thermal SIA emission C2i → 2×Ci1
                    NV_Ith_S(ydot, x_hi_v + 1) =
                        (0.5*P.K_nuc_i*Ci*Ci + P.KLV[2]*Cv*C3i
                         - (P.KLI[1]*Ci + P.KLV[1]*Cv + P.GLV[1])*C2i);
                }
            }
            // nowait: other threads have already started 6 & 7 (disjoint writes)

            // 6 ── dCvx (parallel, independent writes to ydot[2..x_hi_v-1])
            #pragma omp for schedule(static) nowait
            for (int k = 2; k < x_hi_v; ++k) {
                const double Cxm1 = Cv_arr[k - 1];
                const double Cx   = Cv_arr[k];
                double dCx = (P.KCV[k-1]*Cv*Cxm1 - P.KCI[k]*Ci*Cx
                              - P.KCV[k]*Cv*Cx - P.GCV[k]*Cx);
                if (k < x_hi_v - 1)
                    dCx += P.KCI[k+1]*Ci*Cv_arr[k+1] + P.GCV[k+1]*Cv_arr[k+1];
                NV_Ith_S(ydot, k) = dCx;
            }

            // 7 ── dCix (parallel, independent writes to ydot[x_hi_v+1..N-1])
            // Ci_win[] is fully filled (barrier after step 2).
            // Reads Ci_win[j-1] — read-only from pre-filled buffer → no race.
            #pragma omp for schedule(static) nowait
            for (int j = j_ss; j < n_ci_win; ++j) {
                const int    x  = x_lo_i + j;
                const int    k  = x - 1;
                const double Cx = Ci_win[j];
                const double Cxm1 = (j == j_ss && x_lo_i >= 3)
                                    ? W.Ci_frozen_top : Ci_win[j - 1];
                double dCx = (P.KLI[k-1]*Ci*Cxm1
                              - (P.KLI[k]*Ci + P.KLV[k]*Cv)*Cx);
                if (j < n_ci_win - 1)
                    dCx += P.KLV[k+1]*Cv*Ci_win[j+1];
                NV_Ith_S(ydot, x_hi_v + 1 + j) = dCx;
            }

            // 8 ── Floor enforcement [implicit barrier at end of omp for]
            #pragma omp for schedule(static)
            for (int k = 0; k < W.N_active; ++k) {
                if (NV_Ith_S(y, k) <= C_floor && NV_Ith_S(ydot, k) < 0.0)
                    NV_Ith_S(ydot, k) = 0.0;
            }
        } // end omp parallel  — exactly one fork-join per RHS call
#endif
        return 0;
    }

    // =========================================================================
    // SERIAL PATH — pre-allocated buffers, no thread overhead
    // =========================================================================
    // Used when n_ci_win < OMP_MIN_WORK (thread overhead > compute benefit)
    // or when OpenMP is not available.  Identical physics to rhs_window().

    for (int k = 0; k < x_hi_v; ++k)
        Cv_arr[k] = fl(NV_Ith_S(y, k), C_floor);
    for (int j = 0; j < n_ci_win; ++j)
        Ci_win[j] = fl(NV_Ith_S(y, x_hi_v + 1 + j), C_floor);

    const double Cv  = Cv_arr[0];
    const double Ci  = fl(NV_Ith_S(y, x_hi_v), C_floor);
    const double C2i = (x_lo_i == 2) ? Ci_win[0]
                                      : fl(W.full_conc[P.Nv + 1], C_floor);
    const double C2v = (x_hi_v >= 2) ? Cv_arr[1] : 0.0;
    const int j_sink_start = (x_lo_i == 2) ? 1 : 0;

    // dCv1
    {
        double dCv = P.P_prod;
        if (x_hi_v >= 2)
            dCv += P.KCI[1]*Ci*C2v + (2.0*P.GCV[1] - P.KCV[1]*Cv)*C2v;
        for (int k = 2; k < x_hi_v; ++k) {
            dCv += P.GCV[k] * Cv_arr[k];
            dCv -= P.KCV[k] * Cv_arr[k] * Cv;
        }
        double sink_lv = 0.0;
        for (int j = j_sink_start; j < n_ci_win; ++j)
            sink_lv += P.KLV[x_lo_i + j - 1] * Ci_win[j];
        dCv -= P.KLV[1]*Cv*C2i + W.frozen_KLV_sum*Cv + sink_lv*Cv;
        dCv += P.Z_v*P.rho_d*P.Dv*(P.Cv_eq - Cv);
        dCv -= P.alpha*Cv*Ci + P.KCV[0]*Cv*Cv;
        NV_Ith_S(ydot, 0) = dCv;
    }
    // dCv2
    if (x_hi_v >= 2) {
        const double C3v = (x_hi_v >= 3) ? Cv_arr[2] : 0.0;
        double dC2v = (0.5*P.KCV[0]*Cv*Cv + P.GCV[2]*C3v + P.KCI[2]*Ci*C3v
                       + P.rho_d*P.D2v*(P.C2v_eq - C2v));
        dC2v -= (P.KCV[1]*Cv + P.KCI[1]*Ci + P.GCV[1]) * C2v;
        NV_Ith_S(ydot, 1) = dC2v;
    }
    // dCvx
    for (int k = 2; k < x_hi_v; ++k) {
        const double Cxm1 = Cv_arr[k - 1];
        const double Cx   = Cv_arr[k];
        double dCx = (P.KCV[k-1]*Cv*Cxm1 - P.KCI[k]*Ci*Cx
                      - P.KCV[k]*Cv*Cx - P.GCV[k]*Cx);
        if (k < x_hi_v - 1)
            dCx += P.KCI[k+1]*Ci*Cv_arr[k+1] + P.GCV[k+1]*Cv_arr[k+1];
        NV_Ith_S(ydot, k) = dCx;
    }
    // dCi1
    {
        double dCi = P.P_prod;
        dCi += P.KLV[1]*Cv*C2i;
        dCi -= P.K_nuc_i*Ci*Ci + P.alpha*Cv*Ci;
        dCi -= P.KLI[1]*Ci*C2i + P.KCI[1]*Ci*C2v;
        dCi += 2.0*P.GLV[1]*C2i;            // C2i dissociates → 2×Ci1
        for (int k = 2; k < x_hi_v; ++k)
            dCi -= P.KCV[k] * Cv_arr[k] * Ci;
        double sink_li = 0.0;
        for (int j = j_sink_start; j < n_ci_win; ++j)
            sink_li += P.KLI[x_lo_i + j - 1] * Ci_win[j];
        dCi -= W.frozen_KLI_sum*Ci + sink_li*Ci;
        dCi -= P.Z_i*P.rho_d*P.Di*Ci;
        NV_Ith_S(ydot, x_hi_v) = dCi;
    }
    // dCi2
    if (x_lo_i == 2) {
        const double C3i = (n_ci_win > 1) ? Ci_win[1] : 0.0;
        // GLV[1]: thermal SIA emission C2i → 2×Ci1
        NV_Ith_S(ydot, x_hi_v + 1) =
            (0.5*P.K_nuc_i*Ci*Ci + P.KLV[2]*Cv*C3i
             - (P.KLI[1]*Ci + P.KLV[1]*Cv + P.GLV[1])*C2i);
    }
    // dCix
    for (int j = j_sink_start; j < n_ci_win; ++j) {
        const int    x  = x_lo_i + j;
        const int    k  = x - 1;
        const double Cx = Ci_win[j];
        const double Cxm1 = (j == j_sink_start && x_lo_i >= 3)
                            ? W.Ci_frozen_top : Ci_win[j - 1];
        double dCx = (P.KLI[k-1]*Ci*Cxm1
                      - (P.KLI[k]*Ci + P.KLV[k]*Cv)*Cx);
        if (j < n_ci_win - 1)
            dCx += P.KLV[k+1]*Cv*Ci_win[j+1];
        NV_Ith_S(ydot, x_hi_v + 1 + j) = dCx;
    }
    // floor
    for (int k = 0; k < W.N_active; ++k) {
        if (NV_Ith_S(y, k) <= C_floor && NV_Ith_S(ydot, k) < 0.0)
            NV_Ith_S(ydot, k) = 0.0;
    }

    return 0;
}

// ── Phase IV: Jacobi preconditioner (mirrors Phase II, adapted for WindowDataOMP)

int prec_setup_window_omp(sunrealtype /*t*/, N_Vector y, N_Vector /*fy*/,
                           sunbooleantype /*jok*/, sunbooleantype* jcurPtr,
                           sunrealtype gamma, void* user_data) {
    WindowDataOMP&    W        = *static_cast<WindowDataOMP*>(user_data);
    const Parameters& P        = *W.P_full;
    const int x_hi_v           = W.x_hi_v;
    const int x_lo_i           = W.x_lo_i;
    const int x_hi_i           = W.x_hi_i;
    const double C_floor       = P.C_floor;

    const double Cv = fl(NV_Ith_S(y, 0),      C_floor);
    const double Ci = fl(NV_Ith_S(y, x_hi_v), C_floor);

    W.prec_diag.resize(W.N_active);

    // Vacancy block
    for (int k = 0; k < x_hi_v; ++k)
        W.prec_diag[k] = 1.0 - gamma * (-(P.KCV[k]*Cv + P.KCI[k]*Ci + P.GCV[k]));

    // Ci1
    {
        double J_ii = -(P.K_nuc_i*2.0*Ci + P.alpha*Cv + P.Z_i*P.rho_d*P.Di
                        + P.KLI[1]*fl(NV_Ith_S(y, x_hi_v+1), C_floor)
                        + W.frozen_KLI_sum);
        const int n  = x_hi_i - x_lo_i + 1;
        const int j0 = (x_lo_i == 2) ? 1 : 0;
        for (int j = j0; j < n; ++j)
            J_ii -= P.KLI[x_lo_i + j - 1] * fl(NV_Ith_S(y, x_hi_v + 1 + j), C_floor);
        W.prec_diag[x_hi_v] = 1.0 - gamma * J_ii;
    }

    // Active interstitial window
    const int n_ci_win = x_hi_i - x_lo_i + 1;
    for (int j = 0; j < n_ci_win; ++j) {
        const int x = x_lo_i + j;
        W.prec_diag[x_hi_v + 1 + j] =
            1.0 - gamma * (-(P.KLI[x-1]*Ci + P.KLV[x-1]*Cv));
    }

    *jcurPtr = SUNTRUE;
    return 0;
}

int prec_solve_window_omp(sunrealtype /*t*/, N_Vector /*y*/, N_Vector /*fy*/,
                           N_Vector r, N_Vector z,
                           sunrealtype /*gamma*/, sunrealtype /*delta*/,
                           int /*lr*/, void* user_data) {
    const WindowDataOMP& W = *static_cast<const WindowDataOMP*>(user_data);
    for (int i = 0; i < W.N_active; ++i) {
        const double m = W.prec_diag[i];
        NV_Ith_S(z, i) = NV_Ith_S(r, i) / (std::abs(m) > 1e-300 ? m : 1.0);
    }
    return 0;
}
