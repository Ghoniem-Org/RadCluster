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
 *   • Fixed vacancy and He fixed-sink terms to be size-dependent:
 *     only mobile vacancy clusters (m <= m_max_v) diffuse to fixed sinks.
 *     Previously k2_disl_v was applied to ALL sizes, unphysically draining
 *     He from immobile voids and breaking He conservation (δ_He → 1).
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
    if (n < 4 || n > P.i_mobile) return 0.0;
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
    for (int m = 0; m < P.V; ++m)
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
    for (int m = 0; m < P.V; ++m) {
        sink    += P.KHeV[m] * std::max(c_v[m], 0.0);
        He_emit += P.beta_He * std::max(Q_m[m], 0.0);
    }
    double source = P.G_He + He_emit;
    return source / (sink > 1e-300 ? sink : 1e-300);
}

// ── Coalescence rate constants ───────────────────────────────────────────────
//
// K_ii(n, n') = A_sph/Ω^{2/3} · (n^{1/3} + n'^{1/3}) · D_{n'}^{eff}
// Rate for mobile SIA cluster n' encountering SIA cluster n.
// Only the mobile partner's diffusivity enters; target can be sessile.

static inline double K_ii_coal(const Parameters& P, int n, int np) {
    // n, np are 1-indexed sizes
    // Z_ii bias factor accounts for elastic interaction between SIA loops.
    // Prismatic loops have strong stress fields; mutual attraction can
    // significantly enhance the encounter rate beyond geometric capture.
    //
    // Geometry: SIA clusters n ≥ 4 form dislocation loops with capture
    // cross-section ∝ circumference (n^{1/2}), not surface area (n^{1/3}).
    // Small clusters n < 4 (dumbbells) remain spherical.
    // The Z_i_loop bias applies only to the loop target absorbing SIA.
    double D_np = P.D_SIA_eff[np - 1];
    if (D_np < 1e-300) return 0.0;

    // Target (n) geometry: loop for n ≥ 4, spherical for n < 4
    double target_factor, target_pref, target_bias;
    if (n >= 4) {
        target_factor = std::sqrt(static_cast<double>(n));
        target_pref   = P.A_loop_inv_O23;
        target_bias   = P.Z_i_loop;
    } else {
        target_factor = std::cbrt(static_cast<double>(n));
        target_pref   = P.A_sph_inv_O23;
        target_bias   = 1.0;
    }
    // Projectile (np) geometry: loop for np ≥ 4, spherical for np < 4
    double proj_factor;
    if (np >= 4) {
        proj_factor = std::sqrt(static_cast<double>(np));
    } else {
        proj_factor = std::cbrt(static_cast<double>(np));
    }
    // Combined size factor uses average of target and projectile prefactors
    // weighted by their respective geometry.  For simplicity and consistency
    // with the additive cross-section convention, use target geometry for
    // the full rate (target is the absorber; projectile diffuses to it).
    return P.Z_ii * target_bias * target_pref * (target_factor + proj_factor) * D_np;
}

// K_vv(m, m') = A_sph/Ω^{2/3} · (m^{1/3} + m'^{1/3}) · D_{m'}^{eff}
// Rate for mobile vacancy cluster m' encountering vacancy cluster m.

static inline double K_vv_coal(const Parameters& P, int m, int mp) {
    double D_mp = P.D_VAC_eff[mp - 1];
    if (D_mp < 1e-300) return 0.0;
    double size_factor = std::cbrt(static_cast<double>(m))
                       + std::cbrt(static_cast<double>(mp));
    return P.A_sph_inv_O23 * size_factor * D_mp;
}

// K_vi(n, m') = A_sph/Ω^{2/3} · n^{1/3} · D_{m'}^{eff}
// Rate for mobile vacancy cluster m' shrinking SIA cluster n.

static inline double K_vi_coal(const Parameters& P, int n, int mp) {
    double D_mp = P.D_VAC_eff[mp - 1];
    if (D_mp < 1e-300) return 0.0;
    return P.A_sph_inv_O23 * std::cbrt(static_cast<double>(n)) * D_mp;
}

// ── Case 2 — fission/decoupled (Eq. 175) ─────────────────────────────────────
//
// State (dynamic):           [c_i(I) | c_v(V) | Q_tot | c_h]  N_eq = I+V+2
// State (quasi_steady_state):[c_i(I) | c_v(V) | Q_tot]        N_eq = I+V+1
//
// RHS design notes:
//   • All state variables clamped to max(y, 0) at access time.
//     This prevents negative amplification without introducing kinks above 0
//     (which would break CVODE's Jacobian estimation).
//   • C_floor is enforced post-step in solver.cpp, NOT inside the RHS.
//   • KHeV term removed from dcv: He capture preserves marginal Σ_ℓ c_{m,ℓ}.
//   • Fixed-sink loss limited to mobile voids (m+1 ≤ m_max_v).

// x_hi_i_win / x_hi_v_win: inclusive upper bounds for active SIA / VAC state
// indices (0-based).  Full solver passes I-1 / V-1.  Sliding-window modes pass
// the current window frontier.
// use_omp: true only for sliding_OpenMP (window_mode==4) with CD_HAVE_OPENMP.
static int rhs_case2(sunrealtype /*t*/, N_Vector yv, N_Vector ydotv,
                      const Parameters& P,
                      int x_hi_i_win, int x_hi_v_win, bool use_omp) {
    const double* y    = N_VGetArrayPointer_Serial(yv);
    double*       dydt = N_VGetArrayPointer_Serial(ydotv);
    const int  I   = P.I;
    const int  V   = P.V;
    const bool qss = (P.he_options == 1);

    // Unpack and clamp to 0 (CVODE can probe transiently negative values)
    const double* c_i = y;        // [I] — use max(y[n], 0) at access time
    const double* c_v = y + I;    // [V]
    const double  Q_tot = std::max(y[I + V], 0.0);

    // Free He: QSS algebraic or ODE state
    const double c_h = qss ? c_h_qss_case2(P, c_v, Q_tot)
                           : std::max(y[I + V + 1], 0.0);

    const double ci1 = std::max(c_i[0], 0.0);
    const double cv1 = std::max(c_v[0], 0.0);

    for (int k = 0; k < P.N_eq; ++k) dydt[k] = 0.0;

    double* dci = dydt;        // [I]
    double* dcv = dydt + I;    // [V]

    // Mean He per void (scalar, Eq. 175) — for He-pressure correction of GVV
    double C_vac_tot = 0.0;
    for (int m = 0; m < V; ++m) C_vac_tot += std::max(c_v[m], 0.0);
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

    // OpenMP thread count for sliding_OpenMP (mode 4)
    const int omp_threads = (use_omp && P.window_omp_threads > 0)
                            ? P.window_omp_threads : 1;

    // ── SIA cluster equations (Eq. ME_SIA) ─────────────────────────────────
    // Loop is restricted to [0, x_hi_i_win] — clusters beyond the window
    // keep dydt[n]=0 (initialised above), preventing CVODE from evolving them.
    // Each iteration writes only to dci[n], so the loop is race-free under OMP.
#ifdef CD_HAVE_OPENMP
#pragma omp parallel for if(use_omp) schedule(dynamic) num_threads(omp_threads)
#endif
    for (int n = 0; n <= x_hi_i_win; ++n) {
        const int    sn = n + 1;   // 1-indexed size
        const double cn = std::max(c_i[n], 0.0);

        // Production
        dci[n] += P.Pr_SIA[n];

        // Thermal SIA emission (gain + loss)
        if (n + 1 < I) dci[n] += P.GII[n + 1] * std::max(c_i[n + 1], 0.0);
        dci[n] -= P.GII[n] * cn;

        // i–i coalescence with ALL SIA clusters
        // Gain: pairs (np, sn-np) where np is mobile (projectile).
        // With single-D rate constant K(target, projectile) = geo·sf·D_proj,
        // both orderings fire for mobile-mobile pairs, summing D_a + D_b = K_full.
        // No symmetry factor needed (each ordering contributes its own D).
        for (int np = 1; np <= std::min(sn - 1, P.i_mobile); ++np) {
            const int npp = sn - np;   // partner size (1-indexed)
            if (npp < 1 || npp > I) continue;
            const double c_np  = std::max(c_i[np  - 1], 0.0);
            const double c_npp = std::max(c_i[npp - 1], 0.0);
            dci[n] += K_ii_coal(P, npp, np) * c_np * c_npp;
        }
        // Loss (D_np contribution): any mobile np hits this cluster
        for (int np = 1; np <= P.i_mobile; ++np) {
            const double c_np = std::max(c_i[np - 1], 0.0);
            dci[n] -= K_ii_coal(P, sn, np) * cn * c_np;
        }
        // Loss (D_sn contribution): this mobile cluster hits ALL targets
        if (sn <= P.i_mobile) {
            for (int np = 1; np <= I; ++np) {
                const double c_np = std::max(c_i[np - 1], 0.0);
                dci[n] -= K_ii_coal(P, np, sn) * cn * c_np;
            }
        }

        // V–I annihilation: all mobile vacancy clusters m' = 1..v_mobile
        if (n == 0) {
            // n=1, m'=1: P1 recombination (V_1 + I_1 → nothing)
            dci[0] -= P.K_iv * cv1 * ci1;
            // m'=1 gain: V_1 + I_2 → I_1 (was missing — Bug fix)
            if (1 < I)
                dci[0] += P.KIV[1] * cv1 * std::max(c_i[1], 0.0);
            // m'=2..v_mobile: mobile vac cluster shrinks I_{1+m'} → I_1
            for (int mp = 2; mp <= P.v_mobile; ++mp) {
                const double c_mp = std::max(c_v[mp - 1], 0.0);
                if (sn + mp - 1 < I) {
                    // c_i index: target I_{sn+mp} has 0-index (sn+mp-1)
                    const double c_target = std::max(c_i[sn + mp - 1], 0.0);
                    dci[0] += K_vi_coal(P, sn + mp, mp) * c_mp * c_target;
                }
                dci[0] -= K_vi_coal(P, sn, mp) * c_mp * ci1;
            }
        } else {
            // n>=2: all mobile vacancy clusters m' = 1..v_mobile
            for (int mp = 1; mp <= P.v_mobile; ++mp) {
                const double c_mp = std::max(c_v[mp - 1], 0.0);
                double K_shrink;
                if (mp == 1) {
                    K_shrink = P.KIV[n];   // existing rate for mono-vacancy
                } else {
                    K_shrink = K_vi_coal(P, sn, mp);
                }
                // gain from larger cluster I_{n+m'} shrunk by V_{m'}
                if (n + mp < I) {
                    double K_gain = (mp == 1) ? P.KIV[n + mp] :
                                                K_vi_coal(P, sn + mp, mp);
                    dci[n] += K_gain * c_mp * std::max(c_i[n + mp], 0.0);
                }
                dci[n] -= K_shrink * c_mp * cn;
            }
        }

        // SIA cluster–cavity absorption: all mobile n = 1..i_mobile
        if (n == 0) {
            for (int m = 1; m < V; ++m)
                dci[0] -= P.K_3D_cav_pref * P.m13[m] * cn * std::max(c_v[m], 0.0);
        } else if (n < 3 && n < P.i_mobile) {
            for (int m = 0; m < V; ++m)
                dci[n] -= P.K_3D_cav_pref * P.m13[m] * cn * std::max(c_v[m], 0.0);
        } else if (n >= 3 && n < P.i_mobile) {
            for (int m = 0; m < V; ++m)
                dci[n] -= K_1D_eff(P, n, m) * cn * std::max(c_v[m], 0.0);
        }
        // Partial SIA survival (channel a): when I_{n+1} > V_{m+1} (n > m),
        // the SIA is only partially absorbed: I_{n+1} + V_{m+1} → I_{n-m}.
        if (n >= 1 && n < P.i_mobile) {
            for (int m = 0; m < std::min(n, V); ++m) {
                const double cvm = std::max(c_v[m], 0.0);
                const double K_cav = (n < 3) ? P.K_3D_cav_pref * P.m13[m]
                                              : K_1D_eff(P, n, m);
                dci[n - m - 1] += K_cav * cn * cvm;
            }
        }

        // Fixed sinks
        dci[n] -= P.k2_SIA[n] * cn;
    }

    // Window boundary: suppress SIA coalescence reactions whose product
    // exceeds the current window frontier but stays within I.  When
    // I_{k} + I_{np} → I_{k+np} with k+np > wlim (1-indexed window
    // limit), the gain at the product is never computed because the
    // outer loop stops at x_hi_i_win.  Undo the corresponding target
    // and projectile loss so SIA content is conserved across the window.
    if (x_hi_i_win < I - 1) {
        const int wlim = x_hi_i_win + 1;  // 1-indexed window size limit
        for (int np = 1; np <= P.i_mobile; ++np) {
            const double c_np = std::max(c_i[np - 1], 0.0);
            if (c_np < 1e-300) continue;
            // Target sizes k (1-indexed) in window where product k+np
            // overflows the window (k+np > wlim) but not the domain
            // (k+np ≤ I; domain overflow handled by boundary_flux).
            const int k_lo = std::max(wlim - np + 1, 1);
            const int k_hi = std::min(wlim, I - np);
            for (int k = k_lo; k <= k_hi; ++k) {
                const double ck = std::max(c_i[k - 1], 0.0);
                const double rate = K_ii_coal(P, k, np) * c_np * ck;
                dci[k - 1]  += rate;  // undo target loss
                dci[np - 1] += rate;  // undo projectile depletion
            }
        }
    }

    // Reflection boundary: suppress reactions whose product exceeds I
    // Undo target loss and monomer/projectile depletion for overflow reactions.
    if (P.boundary_flux == 1) {
        for (int np = 1; np <= P.i_mobile; ++np) {
            const double c_np = std::max(c_i[np - 1], 0.0);
            for (int k = std::max(I - np + 1, 1); k <= I; ++k) {
                const double rate = K_ii_coal(P, k, np) * c_np * std::max(c_i[k - 1], 0.0);
                dci[k - 1]  += rate;  // undo target loss
                dci[np - 1] += rate;  // undo projectile depletion
            }
        }
    }

    // ── Vacancy cluster equations (Eq. ME_vac, Case 2) ──────────────────────
    // He capture does NOT change void size class m — handled by Q_tot below.
    // Pre-accumulate emitted monomers from thermal vacancy emission within
    // the active VAC window [1, x_hi_v_win].  Out-of-window clusters have
    // near-zero concentration so their contribution is negligible.
    {
        double emit_mono = 0.0;
        for (int m = 1; m <= x_hi_v_win; ++m)
            emit_mono += GVV_eff(m) * std::max(c_v[m], 0.0);
        dcv[0] += emit_mono;
    }
    // Loop restricted to [0, x_hi_v_win]; each iteration writes only to dcv[m].
#ifdef CD_HAVE_OPENMP
#pragma omp parallel for if(use_omp) schedule(dynamic) num_threads(omp_threads)
#endif
    for (int m = 0; m <= x_hi_v_win; ++m) {
        const double cm    = std::max(c_v[m], 0.0);
        const double gvv_m = GVV_eff(m);

        // Production
        dcv[m] += P.Pr_VAC[m];

        // Thermal vacancy emission (gain + loss)
        if (m + 1 < V) dcv[m] += GVV_eff(m + 1) * std::max(c_v[m + 1], 0.0);
        dcv[m] -= gvv_m * cm;

        // V–V coalescence with ALL vacancy clusters
        const int sm = m + 1;   // 1-indexed size
        // Gain: pairs (mp, sm-mp) where mp is mobile (projectile)
        for (int mp = 1; mp <= std::min(sm - 1, P.v_mobile); ++mp) {
            const int mpp = sm - mp;
            if (mpp < 1 || mpp > V) continue;
            const double c_mp  = std::max(c_v[mp  - 1], 0.0);
            const double c_mpp = std::max(c_v[mpp - 1], 0.0);
            dcv[m] += K_vv_coal(P, mpp, mp) * c_mp * c_mpp;
        }
        // Loss (D_mp contribution): any mobile mp hits this cluster
        for (int mp = 1; mp <= P.v_mobile; ++mp) {
            const double c_mp = std::max(c_v[mp - 1], 0.0);
            dcv[m] -= K_vv_coal(P, sm, mp) * cm * c_mp;
        }
        // Loss (D_sm contribution): this mobile cluster hits ALL targets
        if (sm <= P.v_mobile) {
            for (int mp = 1; mp <= V; ++mp) {
                const double c_mp = std::max(c_v[mp - 1], 0.0);
                dcv[m] -= K_vv_coal(P, mp, sm) * cm * c_mp;
            }
        }

        // SIA-induced cavity shrinkage: all mobile SIA n = 1..i_mobile
        // n=1 (monomer):
        if (m == 0) {
            dcv[0] -= P.K_iv * ci1 * cv1;
            if (V >= 2) dcv[0] += P.KVI[1] * ci1 * std::max(c_v[1], 0.0);
            // Vacancy monomer consumed by SIA loop shrinkage:
            // V_1 + I_n → I_{n-1} for n>=2 (n=1 already in K_iv above)
            {
                double sia_shrink_sink = 0.0;
                for (int np = 1; np < I; ++np)
                    sia_shrink_sink += P.KIV[np] * std::max(c_i[np], 0.0);
                dcv[0] -= cv1 * sia_shrink_sink;
            }
        } else {
            dcv[m] -= P.KVI[m] * ci1 * cm;
            if (m + 1 < V) dcv[m] += P.KVI[m + 1] * ci1 * std::max(c_v[m + 1], 0.0);
        }

        // n=2,3 (3D mobile SIA clusters): absorb into all cavities
        for (int n = 1; n < std::min(3, P.i_mobile); ++n) {
            const double cn = std::max(c_i[n], 0.0);
            dcv[m] -= P.KVI[m] * cn * cm;
            if (m + n + 1 < V)
                dcv[m] += P.KVI[m + n + 1] * cn * std::max(c_v[m + n + 1], 0.0);
        }

        // n=4..i_mobile (1D/3D mixed): gain + loss
        for (int n = 3; n < std::min(I, P.i_mobile); ++n) {
            const double cn = std::max(c_i[n], 0.0);
            dcv[m] -= K_1D_eff(P, n, m) * cn * cm;
            if (m + n + 1 < V) {
                const int mp  = m + n + 1;
                const double k_gain = K_1D_eff(P, n, mp);
                dcv[m] += k_gain * cn * std::max(c_v[mp], 0.0);
            }
        }

        // V–I annihilation channel (b): mobile V_{sm} diffuses to SIA.
        // V_1 loss handled above (lines 386–392); add loss for sm = 2..v_mobile
        // and gain at all m from V_{sm+sn} + I_{sn} → V_{sm} (source mobile).
        if (sm >= 2 && sm <= P.v_mobile) {
            for (int sn = 1; sn <= I; ++sn) {
                const double c_sn = std::max(c_i[sn - 1], 0.0);
                dcv[m] -= K_vi_coal(P, sn, sm) * cm * c_sn;
            }
        }
        for (int sn = 1; sm + sn <= P.v_mobile && sn <= I; ++sn) {
            if (sm + sn > V) break;
            dcv[m] += K_vi_coal(P, sn, sm + sn)
                      * std::max(c_v[sm + sn - 1], 0.0)
                      * std::max(c_i[sn - 1], 0.0);
        }

        // Fixed sinks — only mobile vacancy clusters diffuse to sinks
        if (m + 1 <= P.v_mobile)
            dcv[m] -= P.k2_disl_v * cm;
    }

    // Window boundary: suppress VAC coalescence reactions whose product
    // exceeds the current VAC window frontier but stays within V.
    if (x_hi_v_win < V - 1) {
        const int wlim_v = x_hi_v_win + 1;
        for (int mp = 1; mp <= P.v_mobile; ++mp) {
            const double c_mp = std::max(c_v[mp - 1], 0.0);
            if (c_mp < 1e-300) continue;
            const int k_lo = std::max(wlim_v - mp + 1, 1);
            const int k_hi = std::min(wlim_v, V - mp);
            for (int k = k_lo; k <= k_hi; ++k) {
                const double ck = std::max(c_v[k - 1], 0.0);
                const double rate = K_vv_coal(P, k, mp) * c_mp * ck;
                dcv[k - 1]  += rate;  // undo target loss
                dcv[mp - 1] += rate;  // undo projectile depletion
            }
        }
    }

    // Reflection boundary: suppress reactions whose product exceeds V
    if (P.boundary_flux == 1) {
        for (int mp = 1; mp <= P.v_mobile; ++mp) {
            const double c_mp = std::max(c_v[mp - 1], 0.0);
            for (int k = std::max(V - mp + 1, 1); k <= V; ++k) {
                const double rate = K_vv_coal(P, k, mp) * c_mp * std::max(c_v[k - 1], 0.0);
                dcv[k - 1]  += rate;  // undo target loss
                dcv[mp - 1] += rate;  // undo projectile depletion
            }
        }
    }

    // ── Q_tot equation (total He in voids) ───────────────────────────────────
    double He_uptake = 0.0;
    for (int m = 0; m < V; ++m)
        He_uptake += P.KHeV[m] * c_h * std::max(c_v[m], 0.0);
    const double He_emit = P.beta_He * Q_tot;

    // He lost when mobile voids (m <= v_mobile) are absorbed at fixed sinks.
    // Immobile voids do not diffuse to sinks, so their He is retained.
    double He_sink = 0.0;
    for (int m = 0; m < std::min(V, P.v_mobile); ++m) {
        const double ell_m = ell_bar * std::pow(static_cast<double>(m + 1), 2.0/3.0);
        He_sink += P.k2_disl_v * ell_m * std::max(c_v[m], 0.0);
    }
    dydt[I + V] = He_uptake - He_emit - He_sink;

    // ── Free He (Eq. 157) — dynamic mode only ────────────────────────────────
    if (!qss)
        dydt[I + V + 1] = P.G_He - He_uptake - P.k2_disl_He * c_h + He_emit;

    // ── Conservation accounting ODEs ─────────────────────────────────────────
    // J_SIA_fixed: SIA content lost to fixed sinks
    {
        double sia_fixed = 0.0;
        for (int n = 0; n < I; ++n)
            sia_fixed += static_cast<double>(n + 1) * P.k2_SIA[n] * std::max(c_i[n], 0.0);
        dydt[P.N_eq - 5] = sia_fixed;
    }

    // J_SIA_mutual: ALL SIA content lost to SIA-vacancy annihilation.
    // Two channels: (a) mobile SIA hitting voids, (b) mobile vacancies hitting SIA.
    // Weight = min(sn, m+1): when I_n hits V_m with n > m, only m defects
    // are annihilated; the remainder forms a smaller SIA cluster.
    {
        double mutual = 0.0;
        for (int n = 0; n < I; ++n) {
            const int sn = n + 1;
            const double cn = std::max(c_i[n], 0.0);
            if (cn < 1e-300) continue;
            // (a) Mobile SIA cluster → cavity absorption
            if (n == 0) {
                for (int m = 1; m < V; ++m)
                    mutual += P.K_3D_cav_pref * P.m13[m] * cn * std::max(c_v[m], 0.0);
            } else if (n < 3 && n < P.i_mobile) {
                for (int m = 0; m < V; ++m)
                    mutual += std::min(sn, m + 1) * P.K_3D_cav_pref * P.m13[m] * cn * std::max(c_v[m], 0.0);
            } else if (n >= 3 && n < P.i_mobile) {
                for (int m = 0; m < V; ++m)
                    mutual += std::min(sn, m + 1) * K_1D_eff(P, n, m) * cn * std::max(c_v[m], 0.0);
            }
            // (b) Mobile vacancy hitting this SIA cluster
            for (int mp = 1; mp <= P.v_mobile; ++mp) {
                const double c_mp = std::max(c_v[mp - 1], 0.0);
                if (c_mp < 1e-300) continue;
                double K_s;
                if (mp == 1 && n > 0)      K_s = P.KIV[n];
                else if (mp == 1 && n == 0) K_s = P.K_iv;
                else                        K_s = K_vi_coal(P, sn, mp);
                mutual += std::min(mp, sn) * K_s * c_mp * cn;
            }
        }
        dydt[P.N_eq - 4] = mutual;
    }

    // J_VAC_fixed: VAC content lost to fixed sinks
    {
        double vac_fixed = 0.0;
        for (int m = 0; m < std::min(P.v_mobile, P.V); ++m)
            vac_fixed += static_cast<double>(m + 1) * P.k2_disl_v * std::max(c_v[m], 0.0);
        dydt[P.N_eq - 3] = vac_fixed;
    }

    // J_VAC_mutual: VAC content lost to mutual annihilation.
    // For both channels (a) and (b), the vacancy content destroyed per
    // reaction equals min(m', n) — the same as the SIA content destroyed.
    // When V_{m'} hits I_n with m'>n, the vacancy cluster shrinks to
    // V_{m'-n}, losing only n vacancies (not m').
    dydt[P.N_eq - 2] = dydt[P.N_eq - 4];  // J_VAC_mutual = J_SIA_mutual

    // J_He_sink: He lost to sinks
    {
        double he_sink = P.k2_disl_He * c_h;
        if (C_vac_tot > 1e-300 && Q_tot > 0.0) {
            double ell_bar_loc = Q_tot / C_vac_tot;
            for (int m = 0; m < std::min(P.v_mobile, P.V); ++m) {
                double ell_m = ell_bar_loc * std::pow(static_cast<double>(m + 1), 2.0/3.0);
                he_sink += P.k2_disl_v * ell_m * std::max(c_v[m], 0.0);
            }
        }
        dydt[P.N_eq - 1] = he_sink;
    }

    return 0;
}

// ── Case 1 — fusion/mean-field (Eq. 174) ─────────────────────────────────────
//
// State (dynamic):            [c_i(I) | c_v(V) | Q_m(V) | c_h]  N_eq = I+2V+1
// State (quasi_steady_state): [c_i(I) | c_v(V) | Q_m(V)]        N_eq = I+2V

static int rhs_case1(sunrealtype /*t*/, N_Vector yv, N_Vector ydotv,
                      const Parameters& P,
                      int x_hi_i_win, int x_hi_v_win, bool use_omp) {
    const double* y    = N_VGetArrayPointer_Serial(yv);
    double*       dydt = N_VGetArrayPointer_Serial(ydotv);
    const int  I   = P.I;
    const int  V   = P.V;
    const bool qss = (P.he_options == 1);

    const double* c_i = y;
    const double* c_v = y + I;
    const double* Q_m = y + I + V;

    const double c_h = qss ? c_h_qss_case1(P, c_v, Q_m)
                           : std::max(y[I + 2*V], 0.0);

    const double ci1 = std::max(c_i[0], 0.0);
    const double cv1 = std::max(c_v[0], 0.0);

    for (int k = 0; k < P.N_eq; ++k) dydt[k] = 0.0;

    double* dci = dydt;
    double* dcv = dydt + I;
    double* dQ  = dydt + I + V;

    // Per-class He loading ell_bar_m = Q_m / c_m (Eq. 174)
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

    const int omp_threads = (use_omp && P.window_omp_threads > 0)
                            ? P.window_omp_threads : 1;

    // ── SIA clusters (general coalescence, same structure as Case 2) ─────────
    // Loop restricted to [0, x_hi_i_win]; each iteration writes only to dci[n].
#ifdef CD_HAVE_OPENMP
#pragma omp parallel for if(use_omp) schedule(dynamic) num_threads(omp_threads)
#endif
    for (int n = 0; n <= x_hi_i_win; ++n) {
        const int    sn = n + 1;
        const double cn = std::max(c_i[n], 0.0);
        dci[n] += P.Pr_SIA[n];
        if (n + 1 < I) dci[n] += P.GII[n+1] * std::max(c_i[n+1], 0.0);
        dci[n] -= P.GII[n] * cn;

        // i–i coalescence with ALL SIA clusters
        for (int np = 1; np <= std::min(sn - 1, P.i_mobile); ++np) {
            const int npp = sn - np;
            if (npp < 1 || npp > I) continue;
            const double c_np  = std::max(c_i[np  - 1], 0.0);
            const double c_npp = std::max(c_i[npp - 1], 0.0);
            dci[n] += K_ii_coal(P, npp, np) * c_np * c_npp;
        }
        // Loss (D_np contribution): any mobile np hits this cluster
        for (int np = 1; np <= P.i_mobile; ++np) {
            dci[n] -= K_ii_coal(P, sn, np) * cn * std::max(c_i[np - 1], 0.0);
        }
        // Loss (D_sn contribution): this mobile cluster hits ALL targets
        if (sn <= P.i_mobile) {
            for (int np = 1; np <= I; ++np) {
                dci[n] -= K_ii_coal(P, np, sn) * cn * std::max(c_i[np - 1], 0.0);
            }
        }

        // V–I annihilation: all mobile vacancy clusters m' = 1..v_mobile
        if (n == 0) {
            // P1 recombination: V_1 + I_1 → nothing
            dci[0] -= P.K_iv * cv1 * ci1;
            // m'=1 gain: V_1 + I_2 → I_1 (Bug fix — was missing)
            if (1 < I)
                dci[0] += P.KIV[1] * cv1 * std::max(c_i[1], 0.0);
            // m'=2..v_mobile: V_{m'} + I_{1+m'} → I_1
            for (int mp = 2; mp <= P.v_mobile; ++mp) {
                const double c_mp = std::max(c_v[mp - 1], 0.0);
                if (sn + mp - 1 < I) {
                    // target I_{sn+mp} at 0-index (sn+mp-1) — fixed off-by-one
                    dci[0] += K_vi_coal(P, sn + mp, mp) * c_mp * std::max(c_i[sn + mp - 1], 0.0);
                }
                dci[0] -= K_vi_coal(P, sn, mp) * c_mp * ci1;
            }
        } else {
            for (int mp = 1; mp <= P.v_mobile; ++mp) {
                const double c_mp = std::max(c_v[mp - 1], 0.0);
                double K_s = (mp == 1) ? P.KIV[n] : K_vi_coal(P, sn, mp);
                if (n + mp < I) {
                    double K_g = (mp == 1) ? P.KIV[n + mp] : K_vi_coal(P, sn + mp, mp);
                    dci[n] += K_g * c_mp * std::max(c_i[n + mp], 0.0);
                }
                dci[n] -= K_s * c_mp * cn;
            }
        }

        // Cavity absorption
        if (n == 0) {
            for (int m = 1; m < V; ++m)
                dci[0] -= P.K_3D_cav_pref * P.m13[m] * cn * std::max(c_v[m], 0.0);
        } else if (n < 3 && n < P.i_mobile) {
            for (int m = 0; m < V; ++m)
                dci[n] -= P.K_3D_cav_pref * P.m13[m] * cn * std::max(c_v[m], 0.0);
        } else if (n >= 3 && n < P.i_mobile) {
            for (int m = 0; m < V; ++m)
                dci[n] -= K_1D_eff(P, n, m) * cn * std::max(c_v[m], 0.0);
        }
        // Partial SIA survival (channel a): I_{n+1} + V_{m+1} → I_{n-m}
        if (n >= 1 && n < P.i_mobile) {
            for (int m = 0; m < std::min(n, V); ++m) {
                const double cvm = std::max(c_v[m], 0.0);
                const double K_cav = (n < 3) ? P.K_3D_cav_pref * P.m13[m]
                                              : K_1D_eff(P, n, m);
                dci[n - m - 1] += K_cav * cn * cvm;
            }
        }

        dci[n] -= P.k2_SIA[n] * cn;
    }

    // Window boundary: suppress SIA coalescence beyond window frontier.
    if (x_hi_i_win < I - 1) {
        const int wlim = x_hi_i_win + 1;
        for (int np = 1; np <= P.i_mobile; ++np) {
            const double c_np = std::max(c_i[np - 1], 0.0);
            if (c_np < 1e-300) continue;
            const int k_lo = std::max(wlim - np + 1, 1);
            const int k_hi = std::min(wlim, I - np);
            for (int k = k_lo; k <= k_hi; ++k) {
                const double ck = std::max(c_i[k - 1], 0.0);
                const double rate = K_ii_coal(P, k, np) * c_np * ck;
                dci[k - 1]  += rate;
                dci[np - 1] += rate;
            }
        }
    }

    // Reflection boundary: suppress reactions whose product exceeds I
    if (P.boundary_flux == 1) {
        for (int np = 1; np <= P.i_mobile; ++np) {
            const double c_np = std::max(c_i[np - 1], 0.0);
            for (int k = std::max(I - np + 1, 1); k <= I; ++k) {
                const double rate = K_ii_coal(P, k, np) * c_np * std::max(c_i[k - 1], 0.0);
                dci[k - 1]  += rate;  // undo target loss
                dci[np - 1] += rate;  // undo projectile depletion
            }
        }
    }

    // ── Vacancy clusters (general coalescence) ───────────────────────────────
    // Pre-accumulate emitted monomers within active VAC window.
    {
        double emit_mono = 0.0;
        for (int m = 1; m <= x_hi_v_win; ++m)
            emit_mono += GVV_eff_m(m) * std::max(c_v[m], 0.0);
        dcv[0] += emit_mono;
    }
    // Loop restricted to [0, x_hi_v_win]; each iteration writes only to dcv[m].
#ifdef CD_HAVE_OPENMP
#pragma omp parallel for if(use_omp) schedule(dynamic) num_threads(omp_threads)
#endif
    for (int m = 0; m <= x_hi_v_win; ++m) {
        const int    sm  = m + 1;
        const double cm  = std::max(c_v[m], 0.0);
        const double gvv = GVV_eff_m(m);
        dcv[m] += P.Pr_VAC[m];
        if (m + 1 < V) dcv[m] += GVV_eff_m(m+1) * std::max(c_v[m+1], 0.0);
        dcv[m] -= gvv * cm;

        // V–V coalescence with ALL vacancy clusters
        for (int mp = 1; mp <= std::min(sm - 1, P.v_mobile); ++mp) {
            const int mpp = sm - mp;
            if (mpp < 1 || mpp > V) continue;
            const double c_mp  = std::max(c_v[mp  - 1], 0.0);
            const double c_mpp = std::max(c_v[mpp - 1], 0.0);
            dcv[m] += K_vv_coal(P, mpp, mp) * c_mp * c_mpp;
        }
        // Loss (D_mp contribution): any mobile mp hits this cluster
        for (int mp = 1; mp <= P.v_mobile; ++mp) {
            dcv[m] -= K_vv_coal(P, sm, mp) * cm * std::max(c_v[mp - 1], 0.0);
        }
        // Loss (D_sm contribution): this mobile cluster hits ALL targets
        if (sm <= P.v_mobile) {
            for (int mp = 1; mp <= V; ++mp) {
                dcv[m] -= K_vv_coal(P, mp, sm) * cm * std::max(c_v[mp - 1], 0.0);
            }
        }

        // SIA-induced cavity shrinkage
        if (m == 0) {
            dcv[0] -= P.K_iv * ci1 * cv1;
            if (V >= 2) dcv[0] += P.KVI[1] * ci1 * std::max(c_v[1], 0.0);
            // Vacancy monomer consumed by SIA loop shrinkage:
            // V_1 + I_n → I_{n-1} for n>=2 (n=1 already in K_iv above)
            {
                double sia_shrink_sink = 0.0;
                for (int np = 1; np < I; ++np)
                    sia_shrink_sink += P.KIV[np] * std::max(c_i[np], 0.0);
                dcv[0] -= cv1 * sia_shrink_sink;
            }
        } else {
            dcv[m] -= P.KVI[m] * ci1 * cm;
            if (m + 1 < V) dcv[m] += P.KVI[m+1] * ci1 * std::max(c_v[m+1], 0.0);
        }
        for (int n = 1; n < std::min(3, P.i_mobile); ++n) {
            const double cn = std::max(c_i[n], 0.0);
            dcv[m] -= P.KVI[m] * cn * cm;
            if (m + n + 1 < V)
                dcv[m] += P.KVI[m + n + 1] * cn * std::max(c_v[m + n + 1], 0.0);
        }
        for (int n = 3; n < std::min(I, P.i_mobile); ++n) {
            const double cn = std::max(c_i[n], 0.0);
            dcv[m] -= K_1D_eff(P, n, m) * cn * cm;
            if (m + n + 1 < V) {
                const int mp = m + n + 1;
                dcv[m] += K_1D_eff(P, n, mp) * cn * std::max(c_v[mp], 0.0);
            }
        }

        // V–I annihilation channel (b): mobile V_{sm} diffuses to SIA.
        if (sm >= 2 && sm <= P.v_mobile) {
            for (int sn = 1; sn <= I; ++sn) {
                const double c_sn = std::max(c_i[sn - 1], 0.0);
                dcv[m] -= K_vi_coal(P, sn, sm) * cm * c_sn;
            }
        }
        for (int sn = 1; sm + sn <= P.v_mobile && sn <= I; ++sn) {
            if (sm + sn > V) break;
            dcv[m] += K_vi_coal(P, sn, sm + sn)
                      * std::max(c_v[sm + sn - 1], 0.0)
                      * std::max(c_i[sn - 1], 0.0);
        }

        // Fixed sinks — only mobile vacancy clusters diffuse to sinks
        if (m + 1 <= P.v_mobile)
            dcv[m] -= P.k2_disl_v * cm;
    }

    // Window boundary: suppress VAC coalescence beyond window frontier.
    if (x_hi_v_win < V - 1) {
        const int wlim_v = x_hi_v_win + 1;
        for (int mp = 1; mp <= P.v_mobile; ++mp) {
            const double c_mp = std::max(c_v[mp - 1], 0.0);
            if (c_mp < 1e-300) continue;
            const int k_lo = std::max(wlim_v - mp + 1, 1);
            const int k_hi = std::min(wlim_v, V - mp);
            for (int k = k_lo; k <= k_hi; ++k) {
                const double ck = std::max(c_v[k - 1], 0.0);
                const double rate = K_vv_coal(P, k, mp) * c_mp * ck;
                dcv[k - 1]  += rate;
                dcv[mp - 1] += rate;
            }
        }
    }

    // Reflection boundary: suppress reactions whose product exceeds V
    if (P.boundary_flux == 1) {
        for (int mp = 1; mp <= P.v_mobile; ++mp) {
            const double c_mp = std::max(c_v[mp - 1], 0.0);
            for (int k = std::max(V - mp + 1, 1); k <= V; ++k) {
                const double rate = K_vv_coal(P, k, mp) * c_mp * std::max(c_v[k - 1], 0.0);
                dcv[k - 1]  += rate;  // undo target loss
                dcv[mp - 1] += rate;  // undo projectile depletion
            }
        }
    }

    // ── Q_m equations (He content per void class, Eq. 174) ───────────────────
    // Also restricted to the VAC window; Q_m for out-of-window voids stays zero.
    double He_cap_total  = 0.0;
    double He_emit_total = 0.0;
    for (int m = 0; m <= x_hi_v_win; ++m) {
        const double cm       = std::max(c_v[m], 0.0);
        const double qm       = std::max(Q_m[m], 0.0);
        const double he_cap_m = P.KHeV[m] * c_h * cm;
        const double he_emit_m = P.beta_He * qm;
        dQ[m]        += he_cap_m;
        dQ[m]        -= he_emit_m;
        // He lost when mobile voids are absorbed at fixed sinks
        if (m + 1 <= P.v_mobile)
            dQ[m] -= P.k2_disl_v * qm;
        He_cap_total  += he_cap_m;
        He_emit_total += he_emit_m;
    }

    // ── Free He (Eq. 157) — dynamic mode only ────────────────────────────────
    if (!qss)
        dydt[I + 2*V] = P.G_He - He_cap_total - P.k2_disl_He * c_h + He_emit_total;

    // ── Conservation accounting ODEs ─────────────────────────────────────────
    // J_SIA_fixed: SIA content lost to fixed sinks
    {
        double sia_fixed = 0.0;
        for (int n = 0; n < I; ++n)
            sia_fixed += static_cast<double>(n + 1) * P.k2_SIA[n] * std::max(c_i[n], 0.0);
        dydt[P.N_eq - 5] = sia_fixed;
    }

    // J_SIA_mutual: ALL SIA content lost to SIA-vacancy annihilation.
    // Two channels: (a) mobile SIA hitting voids, (b) mobile vacancies hitting SIA.
    // Weight = min(sn, m+1): when I_n hits V_m with n > m, only m defects
    // are annihilated; the remainder forms a smaller SIA cluster.
    {
        double mutual = 0.0;
        for (int n = 0; n < I; ++n) {
            const int sn = n + 1;
            const double cn = std::max(c_i[n], 0.0);
            if (cn < 1e-300) continue;
            // (a) Mobile SIA cluster → cavity absorption
            if (n == 0) {
                for (int m = 1; m < V; ++m)
                    mutual += P.K_3D_cav_pref * P.m13[m] * cn * std::max(c_v[m], 0.0);
            } else if (n < 3 && n < P.i_mobile) {
                for (int m = 0; m < V; ++m)
                    mutual += std::min(sn, m + 1) * P.K_3D_cav_pref * P.m13[m] * cn * std::max(c_v[m], 0.0);
            } else if (n >= 3 && n < P.i_mobile) {
                for (int m = 0; m < V; ++m)
                    mutual += std::min(sn, m + 1) * K_1D_eff(P, n, m) * cn * std::max(c_v[m], 0.0);
            }
            // (b) Mobile vacancy hitting this SIA cluster
            for (int mp = 1; mp <= P.v_mobile; ++mp) {
                const double c_mp = std::max(c_v[mp - 1], 0.0);
                if (c_mp < 1e-300) continue;
                double K_s;
                if (mp == 1 && n > 0)      K_s = P.KIV[n];
                else if (mp == 1 && n == 0) K_s = P.K_iv;
                else                        K_s = K_vi_coal(P, sn, mp);
                mutual += std::min(mp, sn) * K_s * c_mp * cn;
            }
        }
        dydt[P.N_eq - 4] = mutual;
    }

    // J_VAC_fixed: VAC content lost to fixed sinks
    {
        double vac_fixed = 0.0;
        for (int m = 0; m < std::min(P.v_mobile, P.V); ++m)
            vac_fixed += static_cast<double>(m + 1) * P.k2_disl_v * std::max(c_v[m], 0.0);
        dydt[P.N_eq - 3] = vac_fixed;
    }

    // J_VAC_mutual: VAC content lost to mutual annihilation.
    // Vacancy content destroyed per reaction = min(m', n) = SIA content destroyed.
    dydt[P.N_eq - 2] = dydt[P.N_eq - 4];  // J_VAC_mutual = J_SIA_mutual

    // J_He_sink: He lost to sinks
    {
        double C_vac_tot = 0.0;
        for (int m = 0; m < V; ++m) C_vac_tot += std::max(c_v[m], 0.0);
        double Q_tot_loc = 0.0;
        for (int m = 0; m < V; ++m) Q_tot_loc += std::max(Q_m[m], 0.0);
        double he_sink = P.k2_disl_He * c_h;
        if (C_vac_tot > 1e-300 && Q_tot_loc > 0.0) {
            double ell_bar_loc = Q_tot_loc / C_vac_tot;
            for (int m = 0; m < std::min(P.v_mobile, P.V); ++m) {
                double ell_m = ell_bar_loc * std::pow(static_cast<double>(m + 1), 2.0/3.0);
                he_sink += P.k2_disl_v * ell_m * std::max(c_v[m], 0.0);
            }
        }
        dydt[P.N_eq - 1] = he_sink;
    }

    return 0;
}

// ── Dispatch ─────────────────────────────────────────────────────────────────

int rhs_full_CD(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data) {
    const UserData*   ud = static_cast<const UserData*>(user_data);
    const Parameters& P  = *ud->P;

    // Resolve window bounds:
    //   window_mode == 0 (full): entire SIA / VAC domains are active.
    //   window_mode == 3 (cpp_sliding_win): two independent windows, serial.
    //   window_mode == 4 (sliding_OpenMP): two independent windows + OMP.
    const int  x_hi_i = ud->window_active ? ud->x_hi_i : P.I - 1;
    const int  x_hi_v = ud->window_active ? ud->x_hi_v : P.V - 1;
    const bool use_omp = (P.window_mode == 4);

    if (P.he_mode == 1)
        return rhs_case1(t, y, ydot, P, x_hi_i, x_hi_v, use_omp);
    else
        return rhs_case2(t, y, ydot, P, x_hi_i, x_hi_v, use_omp);
}

// ── Size-bin moment RHS (Chapter 9, Eqs. 193-208) ────────────────────────────

int rhs_bin_moment(sunrealtype t, N_Vector yv, N_Vector ydotv, void* user_data) {
    const UserData*   ud = static_cast<const UserData*>(user_data);
    const Parameters& P  = *ud->P;

    const double* y    = N_VGetArrayPointer_Serial(yv);
    double*       dydt = N_VGetArrayPointer_Serial(ydotv);

    const int  Ib  = P.I_bin;
    const int  V   = P.V;
    const int  I   = P.I;
    const bool qss = (P.he_options == 1);

    for (int k = 0; k < P.N_eq; ++k) dydt[k] = 0.0;

    // ── Shape function and moment stride ────────────────────────────────────
    const int  PM  = P.n_mom;  // moments per bin: 1 (constant), 2 (linear), 3 (lognormal)

    // ── Hybrid discrete + binned SIA reconstruction ────────────────────────
    const int i_d = P.i_discrete;
    const int v_d = P.v_discrete;
    const int i_VAC = i_d + PM * Ib;   // first vacancy index in state vector

    // Build SIA bin edges (bins cover sizes i_discrete+1 .. I)
    std::vector<int> n_lo(Ib), n_hi(Ib);
    if (Ib > 0) {
        int edge = i_d + 1;   // bins start after discrete region
        for (int k = 0; k < Ib; ++k) {
            n_lo[k] = edge;
            int next = std::max(static_cast<int>(std::floor(edge * P.r_ratio)), edge + 1);
            n_hi[k]  = std::min(next, I + 1);
            edge      = n_hi[k];
        }
    }

    // Reconstruct full c_n[0..I-1] from discrete + binned
    std::vector<double> c_n(I, 0.0);
    // Discrete sizes: y[0..i_d-1] = c_1..c_{i_discrete}
    for (int n = 0; n < i_d; ++n)
        c_n[n] = std::max(y[n], 0.0);
    // Binned sizes: closure from moments (shape_function selects method)
    for (int k = 0; k < Ib; ++k) {
        const double mu0_k = std::max(y[i_d + PM*k], 0.0);
        const double mu1_k = (PM >= 2) ? y[i_d + PM*k + 1] : 0.0;
        const double mu2_k = (PM >= 3) ? y[i_d + PM*k + 2] : 0.0;
        const double bw    = static_cast<double>(n_hi[k] - n_lo[k]);
        if (bw <= 0 || mu0_k <= 0.0) continue;
        if (bw == 1.0) {
            int ni = n_lo[k] - 1;
            if (ni >= 0 && ni < I) c_n[ni] = mu0_k;
            continue;
        }

        // Select reconstruction method
        bool use_lognormal = false;
        if (P.shape_function == 2 && PM >= 3) {
            const double n_bar  = mu1_k / std::max(mu0_k, 1e-300);
            const double n2_bar = mu2_k / std::max(mu0_k, 1e-300);
            const double ratio  = n2_bar / std::max(n_bar * n_bar, 1e-300);
            if (ratio > 1.5) {
                use_lognormal = true;
                const double sig2 = std::log(ratio);
                const double m_k  = std::log(std::max(n_bar, 1e-300)) - 0.5 * sig2;
                double f_sum = 0.0;
                double max_log_f = -1e300;
                std::vector<double> log_f_arr(static_cast<int>(bw));
                for (int n = n_lo[k]; n < n_hi[k]; ++n) {
                    double nf = static_cast<double>(n);
                    double ln_n = std::log(nf);
                    double lf = std::max(-(ln_n - m_k) * (ln_n - m_k) / (2.0 * sig2), -500.0) - ln_n;
                    log_f_arr[n - n_lo[k]] = lf;
                    if (lf > max_log_f) max_log_f = lf;
                }
                for (int j = 0; j < static_cast<int>(bw); ++j) {
                    log_f_arr[j] -= max_log_f;
                    f_sum += std::exp(log_f_arr[j]);
                }
                if (f_sum < 1e-300) {
                    use_lognormal = false;  // degenerate — fall back
                } else {
                    for (int n = n_lo[k]; n < n_hi[k]; ++n) {
                        if (n - 1 >= 0 && n - 1 < I)
                            c_n[n - 1] = std::max(
                                mu0_k * std::exp(log_f_arr[n - n_lo[k]]) / f_sum, 0.0);
                    }
                }
            }
        }

        if (!use_lognormal) {
            if (P.shape_function == 0) {
                // Piecewise-constant
                const double val = mu0_k / bw;
                for (int n = n_lo[k]; n < n_hi[k]; ++n)
                    if (n - 1 >= 0 && n - 1 < I) c_n[n - 1] = val;
            } else {
                // Linear (hat-function / dual-basis) — default and fallback
                double S1 = 0.0, S2 = 0.0;
                for (int n = n_lo[k]; n < n_hi[k]; ++n) {
                    double nf = static_cast<double>(n);
                    S1 += nf; S2 += nf * nf;
                }
                double det = bw * S2 - S1 * S1;
                if (std::abs(det) < 1e-30) {
                    double val = mu0_k / bw;
                    for (int n = n_lo[k]; n < n_hi[k]; ++n)
                        if (n - 1 >= 0 && n - 1 < I) c_n[n - 1] = val;
                } else {
                    for (int n = n_lo[k]; n < n_hi[k]; ++n) {
                        if (n - 1 < 0 || n - 1 >= I) continue;
                        double nf = static_cast<double>(n);
                        double phi0 = (S2 - S1 * nf) / det;
                        double phi1 = (bw * nf - S1) / det;
                        c_n[n - 1] = std::max(phi0 * mu0_k + phi1 * mu1_k, 0.0);
                    }
                }
            }
        }
    }

    // ── Hybrid discrete + binned vacancy reconstruction ──────────────────
    const int  Kv = P.V_bin;

    // Build vacancy bin edges (bins cover sizes v_discrete+1 .. V)
    std::vector<int> m_lo(Kv), m_hi(Kv);
    if (Kv > 0) {
        int edge = v_d + 1;   // bins start after discrete region
        for (int k = 0; k < Kv; ++k) {
            m_lo[k] = edge;
            int next = std::max(static_cast<int>(std::floor(edge * P.r_ratio)), edge + 1);
            m_hi[k] = std::min(next, V + 1);
            edge     = m_hi[k];
        }
    }

    // Reconstruct full c_v[0..V-1] from discrete + binned
    std::vector<double> c_v_vec(V, 0.0);
    // Discrete sizes: y[i_VAC..i_VAC+v_d-1] = c_1..c_{v_discrete}
    for (int m = 0; m < v_d; ++m)
        c_v_vec[m] = std::max(y[i_VAC + m], 0.0);
    // Binned sizes: closure from (mu0, mu1[, mu2]) using shape_function
    if (Kv > 0) {
        const int vac_mom_start = i_VAC + v_d;
        for (int k = 0; k < Kv; ++k) {
            const double mu0_k = std::max(y[vac_mom_start + PM*k], 0.0);
            const double mu1_k = (PM >= 2) ? y[vac_mom_start + PM*k + 1] : 0.0;
            const double mu2_k = (PM >= 3) ? y[vac_mom_start + PM*k + 2] : 0.0;
            const double bw    = static_cast<double>(m_hi[k] - m_lo[k]);
            if (bw <= 0 || mu0_k <= 0.0) continue;
            if (bw == 1.0) {
                int mi = m_lo[k] - 1;
                if (mi >= 0 && mi < V) c_v_vec[mi] = mu0_k;
                continue;
            }

            // Select reconstruction method
            bool use_lognormal_v = false;
            if (P.shape_function == 2 && PM >= 3) {
                const double n_bar  = mu1_k / std::max(mu0_k, 1e-300);
                const double n2_bar = mu2_k / std::max(mu0_k, 1e-300);
                const double ratio  = n2_bar / std::max(n_bar * n_bar, 1e-300);
                if (ratio > 1.5) {
                    use_lognormal_v = true;
                    const double sig2 = std::log(ratio);
                    const double m_k  = std::log(std::max(n_bar, 1e-300)) - 0.5 * sig2;
                    double f_sum = 0.0;
                    double max_log_f = -1e300;
                    std::vector<double> log_f_arr(static_cast<int>(bw));
                    for (int m = m_lo[k]; m < m_hi[k]; ++m) {
                        double mf = static_cast<double>(m);
                        double ln_m = std::log(mf);
                        double lf = std::max(-(ln_m - m_k) * (ln_m - m_k) / (2.0 * sig2), -500.0) - ln_m;
                        log_f_arr[m - m_lo[k]] = lf;
                        if (lf > max_log_f) max_log_f = lf;
                    }
                    for (int j = 0; j < static_cast<int>(bw); ++j) {
                        log_f_arr[j] -= max_log_f;
                        f_sum += std::exp(log_f_arr[j]);
                    }
                    if (f_sum < 1e-300) {
                        use_lognormal_v = false;
                    } else {
                        for (int m = m_lo[k]; m < m_hi[k]; ++m) {
                            if (m - 1 >= 0 && m - 1 < V)
                                c_v_vec[m - 1] = std::max(
                                    mu0_k * std::exp(log_f_arr[m - m_lo[k]]) / f_sum, 0.0);
                        }
                    }
                }
            }

            if (!use_lognormal_v) {
                if (P.shape_function == 0) {
                    // Piecewise-constant
                    double val = mu0_k / bw;
                    for (int m = m_lo[k]; m < m_hi[k]; ++m)
                        if (m - 1 >= 0 && m - 1 < V) c_v_vec[m - 1] = std::max(val, 0.0);
                } else {
                    // Linear (hat-function) — default and fallback
                    double S1 = 0.0, S2 = 0.0;
                    for (int m = m_lo[k]; m < m_hi[k]; ++m) {
                        double mf = static_cast<double>(m);
                        S1 += mf; S2 += mf * mf;
                    }
                    double det = bw * S2 - S1 * S1;
                    if (std::abs(det) < 1e-30) {
                        double val = mu0_k / bw;
                        for (int m = m_lo[k]; m < m_hi[k]; ++m)
                            if (m - 1 >= 0 && m - 1 < V) c_v_vec[m - 1] = std::max(val, 0.0);
                    } else {
                        for (int m = m_lo[k]; m < m_hi[k]; ++m) {
                            if (m - 1 < 0 || m - 1 >= V) continue;
                            double mf = static_cast<double>(m);
                            double phi0 = (S2 - S1 * mf) / det;
                            double phi1 = (bw * mf - S1) / det;
                            c_v_vec[m - 1] = std::max(phi0 * mu0_k + phi1 * mu1_k, 0.0);
                        }
                    }
                }
            }
        }
    }
    const double* c_v = c_v_vec.data();

    // He state
    int    i_He_idx = -1;
    double Q_tot    = 0.0;
    const int i_Q_base = i_VAC + v_d + PM * Kv;  // after discrete VAC + binned VAC moments

    if (P.he_mode == 1) {
        int n_Q = Kv;  // Q_k per vacancy bin
        for (int k = 0; k < n_Q; ++k) Q_tot += std::max(y[i_Q_base + k], 0.0);
        i_He_idx = qss ? -1 : i_Q_base + n_Q;
    } else {
        Q_tot    = std::max(y[i_Q_base], 0.0);
        i_He_idx = qss ? -1 : i_Q_base + 1;
    }

    const double c_h = qss
        ? (P.he_mode == 1 ? c_h_qss_case1(P, c_v, y + i_Q_base)
                          : c_h_qss_case2(P, c_v, Q_tot))
        : std::max(y[i_He_idx], 0.0);

    const double ci1 = c_n.empty() ? 0.0 : std::max(c_n[0], 0.0);
    const double cv1 = std::max(c_v[0], 0.0);

    // dc_n/dt per-size
    // IMPORTANT: For bin-moment mode, general SIA-SIA coalescence is
    // restricted to DISCRETE sizes only (n, n' <= i_discrete).  Binned
    // sizes use piecewise-constant reconstruction, which amplifies
    // coalescence rates by O(bw^2) in wide bins.  Monomer-driven growth
    // (monomer absorbed by any cluster) is still computed for all sizes.
    std::vector<double> dc_n(I, 0.0);
    for (int n = 0; n < I; ++n) {
        const int    sn = n + 1;
        const double cn = std::max(c_n[n], 0.0);
        dc_n[n] += P.Pr_SIA[n];
        if (n + 1 < I) dc_n[n] += P.GII[n+1] * std::max(c_n[n+1], 0.0);
        dc_n[n] -= P.GII[n] * cn;

        // i–i coalescence: discrete–discrete (n, n' <= i_discrete, np >= 2)
        // Monomer growth (np=1) is handled separately below for ALL sizes.
        // Discrete–binned coalescence handled in a separate block after this loop.
        if (sn <= i_d) {
            for (int np = 2; np <= std::min(sn - 1, i_d); ++np) {
                const int npp = sn - np;
                if (npp < 1 || npp > i_d) continue;
                const double c_np  = std::max(c_n[np  - 1], 0.0);
                const double c_npp = std::max(c_n[npp - 1], 0.0);
                dc_n[n] += K_ii_coal(P, npp, np) * c_np * c_npp;
            }
            // Loss (D_np contribution): mobile discrete np>=2 hits this cluster
            for (int np = 2; np <= std::min(i_d, P.i_mobile); ++np) {
                dc_n[n] -= K_ii_coal(P, sn, np) * cn * std::max(c_n[np - 1], 0.0);
            }
            // Loss (D_sn contribution): this mobile cluster hits ALL discrete targets
            // Completes the full rate K_full = K(sn,np;D_np) + K(np,sn;D_sn)
            // for discrete–discrete pairs.  Binned targets handled separately below.
            if (sn >= 2 && sn <= P.i_mobile) {
                for (int np = 1; np <= i_d; ++np) {
                    dc_n[n] -= K_ii_coal(P, np, sn) * cn * std::max(c_n[np - 1], 0.0);
                }
            }
        }

        // Monomer growth: I_1 + I_n → I_{n+1}  (all sizes, including binned)
        // Gain from I_1 + I_{n-1} → I_n
        if (n > 0)
            dc_n[n] += K_ii_coal(P, n, 1) * ci1 * std::max(c_n[n - 1], 0.0);
        // Loss: I_n + I_1 → I_{n+1}
        dc_n[n] -= K_ii_coal(P, sn, 1) * cn * ci1;
        // Monomer projectile depletion (n=0 only): I_1 consumed by all clusters
        if (n == 0) {
            double mono_sink = 0.0;
            for (int np = 1; np <= I; ++np)
                mono_sink += K_ii_coal(P, np, 1) * std::max(c_n[np - 1], 0.0);
            dc_n[0] -= ci1 * mono_sink;
        }

        // V–I annihilation: all mobile vacancy clusters m' = 1..v_mobile
        if (n == 0) {
            // P1 recombination: V_1 + I_1 → nothing
            dc_n[0] -= P.K_iv * cv1 * ci1;
            // m'=1 gain: V_1 + I_2 → I_1 (Bug fix — was missing)
            if (1 < I)
                dc_n[0] += P.KIV[1] * cv1 * std::max(c_n[1], 0.0);
            // m'=2..v_mobile: V_{m'} + I_{1+m'} → I_1
            for (int mp = 2; mp <= P.v_mobile; ++mp) {
                const double c_mp = std::max(c_v[mp - 1], 0.0);
                if (sn + mp - 1 < I) {
                    // target I_{sn+mp} at 0-index (sn+mp-1) — fixed off-by-one
                    dc_n[0] += K_vi_coal(P, sn + mp, mp) * c_mp * std::max(c_n[sn + mp - 1], 0.0);
                }
                dc_n[0] -= K_vi_coal(P, sn, mp) * c_mp * ci1;
            }
        } else {
            for (int mp = 1; mp <= P.v_mobile; ++mp) {
                const double c_mp = std::max(c_v[mp - 1], 0.0);
                double K_s = (mp == 1) ? P.KIV[n] : K_vi_coal(P, sn, mp);
                if (n + mp < I) {
                    double K_g = (mp == 1) ? P.KIV[n + mp] : K_vi_coal(P, sn + mp, mp);
                    dc_n[n] += K_g * c_mp * std::max(c_n[n + mp], 0.0);
                }
                dc_n[n] -= K_s * c_mp * cn;
            }
        }

        if (n == 0) {
            for (int m = 1; m < V; ++m)
                dc_n[0] -= P.K_3D_cav_pref * P.m13[m] * cn * std::max(c_v[m], 0.0);
        } else if (n < 3 && n < P.i_mobile) {
            for (int m = 0; m < V; ++m)
                dc_n[n] -= P.K_3D_cav_pref * P.m13[m] * cn * std::max(c_v[m], 0.0);
        } else if (n >= 3 && n < P.i_mobile) {
            for (int m = 0; m < V; ++m)
                dc_n[n] -= K_1D_eff(P, n, m) * cn * std::max(c_v[m], 0.0);
        }
        // Partial SIA survival (channel a): I_{n+1} + V_{m+1} → I_{n-m}
        if (n >= 1 && n < P.i_mobile) {
            for (int m = 0; m < std::min(n, V); ++m) {
                const double cvm = std::max(c_v[m], 0.0);
                const double K_cav = (n < 3) ? P.K_3D_cav_pref * P.m13[m]
                                              : K_1D_eff(P, n, m);
                dc_n[n - m - 1] += K_cav * cn * cvm;
            }
        }

        dc_n[n] -= P.k2_SIA[n] * cn;
    }

    // Reflection boundary: suppress overflow reactions whose product exceeds I
    if (P.boundary_flux == 1) {
        // Monomer overflow: I_1 + I_I → I_{I+1}
        const double rate_i = P.KII[I - 1] * ci1 * std::max(c_n[I - 1], 0.0);
        dc_n[I - 1] += rate_i;  // undo target loss
        dc_n[0]     += rate_i;  // undo monomer depletion
        // Mobile cluster overflow: I_sn + I_np → I_{sn+np} where sn+np > I
        // Undo both D_np and D_sn loss contributions for discrete–discrete pairs
        for (int sn = 2; sn <= std::min(i_d, P.i_mobile); ++sn) {
            const double c_sn = std::max(c_n[sn - 1], 0.0);
            if (c_sn < 1e-300) continue;
            for (int np = std::max(I - sn + 1, 1); np <= i_d; ++np) {
                const double c_np = std::max(c_n[np - 1], 0.0);
                if (c_np < 1e-300) continue;
                // Undo D_np loss (np hits sn): applied in D_np loop
                if (np >= 2 && np <= P.i_mobile) {
                    const double rate_np = K_ii_coal(P, sn, np) * c_sn * c_np;
                    dc_n[sn - 1] += rate_np;
                }
                // Undo D_sn loss (sn hits np): applied in new D_sn loop
                const double rate_sn = K_ii_coal(P, np, sn) * c_sn * c_np;
                dc_n[sn - 1] += rate_sn;
            }
        }
    }

    // ── Mobile discrete SIA (sn=2..i_mobile): D_sn coalescence ────────────
    // Adds the "this cluster diffuses to ALL targets" loss channel and
    // the compensating moment corrections for target/product bins.
    // (sn=1 monomer is fully handled by monomer growth + monomer depletion.)
    //
    // Corrections are accumulated at the MOMENT level (per bin), NOT at the
    // per-size level, to avoid dense cross-bin Jacobian coupling that would
    // collapse CVODE step size.  This is exact for piecewise-constant bins
    // and O(bw) accurate for linear/lognormal.
    //
    // For reaction: mobile sn + immobile np → product(sn+np)
    //   rate = K(np, sn) * c_sn * c_np = Z_ii * A * (np^{1/3}+sn^{1/3}) * D_sn * c_sn * c_np
    //
    //   dc_n[n]     -= rate      (projectile loss — into per-size dc_n)
    //   dmu0[k_tgt] -= rate      (target bin: one cluster removed)
    //   dmu1[k_tgt] -= np*rate   (target bin: np atoms removed)
    //   dmu0[k_prd] += rate      (product bin: one cluster added)
    //   dmu1[k_prd] += (sn+np)*rate (product bin: sn+np atoms added)
    //
    // Net inventory change: -sn*rate (projectile) + sn*rate (product - target) = 0. ✓

    // Accumulate moment corrections per SIA bin [indexed 0..Ib-1]
    std::vector<double> coal_dmu0(Ib, 0.0), coal_dmu1(Ib, 0.0), coal_dmu2(Ib, 0.0);

    for (int n = 1; n < std::min(i_d, P.i_mobile); ++n) {
        const int sn = n + 1;
        const double cn = std::max(c_n[n], 0.0);
        if (cn < 1e-300) continue;
        const double D_sn = P.D_SIA_eff[sn - 1];
        if (D_sn < 1e-300) continue;
        const double pref = P.Z_ii * P.A_sph_inv_O23 * D_sn * cn;
        const double sn_cbrt = std::cbrt(static_cast<double>(sn));

        // Single pass over all SIA bins: target loss + product gain + projectile sink
        int kp = 0;  // product bin index (monotonic scan)
        for (int k = 0; k < Ib; ++k) {
            for (int np = n_lo[k]; np < n_hi[k]; ++np) {
                if (np - 1 < 0 || np - 1 >= I) continue;
                const double c_np = std::max(c_n[np - 1], 0.0);
                if (c_np < P.C_floor) continue;  // skip negligible sizes
                const double nf = static_cast<double>(np);
                const double sf = std::cbrt(nf) + sn_cbrt;
                const double rate = pref * sf * c_np;

                // Projectile loss (per-size dc_n)
                dc_n[n] -= rate;

                // Target bin: one cluster and np atoms removed
                coal_dmu0[k] -= rate;
                coal_dmu1[k] -= nf * rate;
                if (PM >= 3) coal_dmu2[k] -= nf * nf * rate;

                // Product bin: one cluster and (sn+np) atoms added
                const int prod = sn + np;
                if (prod <= I) {
                    if (prod <= i_d) {
                        dc_n[prod - 1] += rate;
                    } else {
                        while (kp < Ib && prod >= n_hi[kp]) ++kp;
                        if (kp < Ib && prod >= n_lo[kp]) {
                            const double pf = static_cast<double>(prod);
                            coal_dmu0[kp] += rate;
                            coal_dmu1[kp] += pf * rate;
                            if (PM >= 3) coal_dmu2[kp] += pf * pf * rate;
                        }
                    }
                } else if (P.boundary_flux == 1) {
                    // Reflection: suppress the overflow reaction entirely
                    dc_n[n] += rate;              // undo projectile loss
                    coal_dmu0[k] += rate;          // undo target bin mu0 loss
                    coal_dmu1[k] += nf * rate;     // undo target bin mu1 loss
                    if (PM >= 3) coal_dmu2[k] += nf * nf * rate;
                }
            }
        }
    }

    // ── Project dc_n into dydt ────────────────────────────────────────────
    // Discrete sizes: direct copy
    for (int n = 0; n < i_d; ++n)
        dydt[n] = dc_n[n];
    // Binned sizes: project onto PM moments per bin + add coalescence corrections
    for (int k = 0; k < Ib; ++k) {
        double dmu0 = 0.0, dmu1 = 0.0, dmu2 = 0.0;
        for (int n = n_lo[k]; n < n_hi[k]; ++n) {
            if (n - 1 >= 0 && n - 1 < I) {
                const double nf = static_cast<double>(n);
                dmu0 += dc_n[n - 1];
                if (PM >= 2) dmu1 += nf * dc_n[n - 1];
                if (PM >= 3) dmu2 += nf * nf * dc_n[n - 1];
            }
        }
        dydt[i_d + PM*k] = dmu0 + coal_dmu0[k];
        if (PM >= 2) dydt[i_d + PM*k + 1] = dmu1 + coal_dmu1[k];
        if (PM >= 3) dydt[i_d + PM*k + 2] = dmu2 + coal_dmu2[k];
    }

    // Note: no inter-bin flux needed — the per-size dc_n projection
    // already accounts for growth/emission across bin boundaries exactly.

    // ── Vacancy per-size dc_v/dt → project onto vacancy bin moments ─────────
    double C_vac_tot = 0.0;
    for (int m = 0; m < V; ++m) C_vac_tot += std::max(c_v[m], 0.0);
    const double ell_bar = (C_vac_tot > 1e-300) ? Q_tot / C_vac_tot : 0.0;

    // Precompute He-modified vacancy emission rates for all sizes
    std::vector<double> GVV_eff(V);
    for (int m = 0; m < V; ++m) {
        const double ell_m = ell_bar * std::pow(static_cast<double>(m+1), 2.0/3.0);
        const double ratio = (ell_m > 1e-6) ? ell_m / (m+1) : 0.0;
        const double dE    = (ratio > 1e-10) ? P.delta_He * P.beta_He_exp / (m+1)
                                               * std::pow(ratio, P.beta_He_exp - 1.0) : 0.0;
        GVV_eff[m] = P.GVV[m] * std::exp(std::min(-ell_m * dE / P.kBT, 0.0));
    }

    std::vector<double> dc_v(V, 0.0);
    // Pre-accumulate emitted monomers from thermal vacancy emission:
    // V_m → V_{m-1} + V_1 — the emitted V_1 monomer must be added to dc_v[0].
    {
        double emit_mono = 0.0;
        for (int m = 1; m < V; ++m)
            emit_mono += GVV_eff[m] * std::max(c_v[m], 0.0);
        dc_v[0] += emit_mono;
    }
    for (int m = 0; m < V; ++m) {
        const double cm    = std::max(c_v[m], 0.0);
        const int sm = m + 1;
        dc_v[m] += P.Pr_VAC[m];
        // Thermal emission: gain from m+1 and loss from m — same He-modified rate
        if (m + 1 < V) dc_v[m] += GVV_eff[m+1] * std::max(c_v[m+1], 0.0);
        dc_v[m] -= GVV_eff[m] * cm;

        // V–V coalescence with ALL vacancy clusters
        for (int mp = 1; mp <= std::min(sm - 1, P.v_mobile); ++mp) {
            const int mpp = sm - mp;
            if (mpp < 1 || mpp > V) continue;
            const double c_mp  = std::max(c_v[mp  - 1], 0.0);
            const double c_mpp = std::max(c_v[mpp - 1], 0.0);
            dc_v[m] += K_vv_coal(P, mpp, mp) * c_mp * c_mpp;
        }
        // Loss (D_mp contribution): any mobile mp hits this cluster
        for (int mp = 1; mp <= P.v_mobile; ++mp)
            dc_v[m] -= K_vv_coal(P, sm, mp) * cm * std::max(c_v[mp - 1], 0.0);
        // Loss (D_sm contribution): this mobile cluster hits ALL targets
        if (sm <= P.v_mobile) {
            for (int mp = 1; mp <= V; ++mp)
                dc_v[m] -= K_vv_coal(P, mp, sm) * cm * std::max(c_v[mp - 1], 0.0);
        }

        // SIA-induced cavity shrinkage
        if (m == 0) {
            dc_v[0] -= P.K_iv * ci1 * cv1;
            if (V >= 2) dc_v[0] += P.KVI[1] * ci1 * std::max(c_v[1], 0.0);
            // Vacancy monomer consumed by SIA loop shrinkage:
            // I_n + V_1 -> I_{n-1} for n>=2 (n=1 already in K_iv above)
            double sia_shrink_sink = 0.0;
            for (int np = 1; np < I; ++np)
                sia_shrink_sink += P.KIV[np] * std::max(c_n[np], 0.0);
            dc_v[0] -= cv1 * sia_shrink_sink;
        } else {
            dc_v[m] -= P.KVI[m] * ci1 * cm;
            if (m + 1 < V) dc_v[m] += P.KVI[m+1] * ci1 * std::max(c_v[m+1], 0.0);
        }

        // n=2,3 (3D mobile SIA clusters): absorb into all cavities
        for (int n = 1; n < std::min(3, P.i_mobile); ++n) {
            const double cn = std::max(c_n[n], 0.0);
            dc_v[m] -= P.KVI[m] * cn * cm;
            if (m + n + 1 < V)
                dc_v[m] += P.KVI[m + n + 1] * cn * std::max(c_v[m + n + 1], 0.0);
        }

        // n=4..i_mobile (1D/3D mixed): gain + loss
        for (int n = 3; n < std::min(I, P.i_mobile); ++n) {
            const double cn = std::max(c_n[n], 0.0);
            dc_v[m] -= K_1D_eff(P, n, m) * cn * cm;
            if (m + n + 1 < V) {
                const int mp  = m + n + 1;
                const double k_gain = K_1D_eff(P, n, mp);
                dc_v[m] += k_gain * cn * std::max(c_v[mp], 0.0);
            }
        }

        // V–I annihilation channel (b): mobile V_{sm} diffuses to SIA.
        if (sm >= 2 && sm <= P.v_mobile) {
            for (int sn = 1; sn <= I; ++sn) {
                const double c_sn = std::max(c_n[sn - 1], 0.0);
                dc_v[m] -= K_vi_coal(P, sn, sm) * cm * c_sn;
            }
        }
        for (int sn = 1; sm + sn <= P.v_mobile && sn <= I; ++sn) {
            if (sm + sn > V) break;
            dc_v[m] += K_vi_coal(P, sn, sm + sn)
                       * std::max(c_v[sm + sn - 1], 0.0)
                       * std::max(c_n[sn - 1], 0.0);
        }

        // Fixed sinks — only mobile vacancy clusters diffuse to sinks
        if (m + 1 <= P.v_mobile)
            dc_v[m] -= P.k2_disl_v * cm;
    }

    // Reflection boundary: suppress reactions whose product exceeds V
    if (P.boundary_flux == 1) {
        for (int mp = 1; mp <= P.v_mobile; ++mp) {
            const double c_mp = std::max(c_v[mp - 1], 0.0);
            for (int k = std::max(V - mp + 1, 1); k <= V; ++k) {
                const double rate = K_vv_coal(P, k, mp) * c_mp * std::max(c_v[k - 1], 0.0);
                dc_v[k - 1]  += rate;  // undo target loss
                dc_v[mp - 1] += rate;  // undo projectile depletion
            }
        }
    }

    // ── Project dc_v into dydt ────────────────────────────────────────────
    // Discrete sizes: direct copy
    for (int m = 0; m < v_d; ++m)
        dydt[i_VAC + m] = dc_v[m];
    // Binned sizes: project onto PM moments per bin
    const int vac_mom_start = i_VAC + v_d;
    if (Kv > 0) {
        for (int k = 0; k < Kv; ++k) {
            double dmu0 = 0.0, dmu1 = 0.0, dmu2 = 0.0;
            for (int m = m_lo[k]; m < m_hi[k]; ++m) {
                if (m - 1 >= 0 && m - 1 < V) {
                    const double mf = static_cast<double>(m);
                    dmu0 += dc_v[m - 1];
                    if (PM >= 2) dmu1 += mf * dc_v[m - 1];
                    if (PM >= 3) dmu2 += mf * mf * dc_v[m - 1];
                }
            }
            dydt[vac_mom_start + PM*k] = dmu0;
            if (PM >= 2) dydt[vac_mom_start + PM*k + 1] = dmu1;
            if (PM >= 3) dydt[vac_mom_start + PM*k + 2] = dmu2;
        }
    }

    // ── He equations ─────────────────────────────────────────────────────────
    if (P.he_mode != 1) {
        // Case 2: scalar Q_tot
        double He_up = 0.0;
        for (int m = 0; m < V; ++m)
            He_up += P.KHeV[m] * c_h * std::max(c_v[m], 0.0);
        const double He_emit = P.beta_He * Q_tot;
        // He lost only from mobile voids reaching fixed sinks
        double He_sink_2 = 0.0;
        for (int m = 0; m < std::min(V, P.v_mobile); ++m) {
            const double ell_m = ell_bar * std::pow(static_cast<double>(m + 1), 2.0/3.0);
            He_sink_2 += P.k2_disl_v * ell_m * std::max(c_v[m], 0.0);
        }
        dydt[i_Q_base] = He_up - He_emit - He_sink_2;
        if (!qss)
            dydt[i_He_idx] = P.G_He - He_up - P.k2_disl_He * c_h + He_emit;
    } else {
        // Case 1: Q per bin (Kv>0) or Q per size (Kv==0)
        int n_Q = Kv;
        double He_cap_total = 0.0, He_emit_total = 0.0;
        if (Kv > 0) {
            for (int k = 0; k < Kv; ++k) {
                const double mu0_k = std::max(y[vac_mom_start + PM*k], 0.0);
                const double qk    = std::max(y[i_Q_base + k], 0.0);
                const double ell_k = (mu0_k > 1e-200) ? qk / mu0_k : 0.0;
                double dqk = 0.0;
                for (int m = m_lo[k]; m < m_hi[k]; ++m) {
                    int mi = m - 1;
                    if (mi < 0 || mi >= V) continue;
                    const double cm = std::max(c_v[mi], 0.0);
                    const double he_cap  = P.KHeV[mi] * c_h * cm;
                    const double q_approx = ell_k * cm;
                    const double he_emit  = P.beta_He * q_approx;
                    dqk += he_cap - he_emit;
                    if (m <= P.v_mobile)   // only mobile voids
                        dqk -= P.k2_disl_v * q_approx;
                    He_cap_total  += he_cap;
                    He_emit_total += he_emit;
                }
                dydt[i_Q_base + k] = dqk;
            }
        } else {
            for (int m = 0; m < V; ++m) {
                const double cm = std::max(c_v[m], 0.0);
                const double qm = std::max(y[i_Q_base + m], 0.0);
                const double he_cap_m  = P.KHeV[m] * c_h * cm;
                const double he_emit_m = P.beta_He * qm;
                dydt[i_Q_base + m] += he_cap_m - he_emit_m;
                if (m + 1 <= P.v_mobile)   // only mobile voids
                    dydt[i_Q_base + m] -= P.k2_disl_v * qm;
                He_cap_total  += he_cap_m;
                He_emit_total += he_emit_m;
            }
        }
        if (!qss)
            dydt[i_He_idx] = P.G_He - He_cap_total - P.k2_disl_He * c_h + He_emit_total;
    }

    // ── Conservation accounting ODEs ─────────────────────────────────────────
    // J_SIA_fixed: SIA content lost to fixed sinks
    {
        double sia_fixed = 0.0;
        for (int n = 0; n < I; ++n)
            sia_fixed += static_cast<double>(n + 1) * P.k2_SIA[n] * std::max(c_n[n], 0.0);
        dydt[P.N_eq - 5] = sia_fixed;
    }

    // J_SIA_mutual: ALL SIA content lost to SIA-vacancy annihilation (bin_moment).
    // Weight = min(sn, m+1): when I_n hits V_m with n > m, only m defects
    // are annihilated; the remainder forms a smaller SIA cluster.
    {
        double mutual = 0.0;
        for (int n = 0; n < I; ++n) {
            const int sn = n + 1;
            const double cn_val = std::max(c_n[n], 0.0);
            if (cn_val < 1e-300) continue;
            // (a) Mobile SIA cluster → cavity absorption
            if (n == 0) {
                for (int m = 1; m < V; ++m)
                    mutual += P.K_3D_cav_pref * P.m13[m] * cn_val * std::max(c_v[m], 0.0);
            } else if (n < 3 && n < P.i_mobile) {
                for (int m = 0; m < V; ++m)
                    mutual += std::min(sn, m + 1) * P.K_3D_cav_pref * P.m13[m] * cn_val * std::max(c_v[m], 0.0);
            } else if (n >= 3 && n < P.i_mobile) {
                for (int m = 0; m < V; ++m)
                    mutual += std::min(sn, m + 1) * K_1D_eff(P, n, m) * cn_val * std::max(c_v[m], 0.0);
            }
            // (b) Mobile vacancy hitting this SIA cluster
            for (int mp = 1; mp <= P.v_mobile; ++mp) {
                const double c_mp = std::max(c_v[mp - 1], 0.0);
                if (c_mp < 1e-300) continue;
                double K_s;
                if (mp == 1 && n > 0)      K_s = P.KIV[n];
                else if (mp == 1 && n == 0) K_s = P.K_iv;
                else                        K_s = K_vi_coal(P, sn, mp);
                mutual += std::min(mp, sn) * K_s * c_mp * cn_val;
            }
        }
        dydt[P.N_eq - 4] = mutual;
    }

    // J_VAC_fixed: VAC content lost to fixed sinks
    {
        double vac_fixed = 0.0;
        for (int m = 0; m < std::min(P.v_mobile, P.V); ++m)
            vac_fixed += static_cast<double>(m + 1) * P.k2_disl_v * std::max(c_v[m], 0.0);
        dydt[P.N_eq - 3] = vac_fixed;
    }

    // J_VAC_mutual: VAC content lost to mutual annihilation.
    // For both channels (a) and (b), the vacancy content destroyed per
    // reaction equals min(m', n) — the same as the SIA content destroyed.
    // When V_{m'} hits I_n with m'>n, the vacancy cluster shrinks to
    // V_{m'-n}, losing only n vacancies (not m').
    dydt[P.N_eq - 2] = dydt[P.N_eq - 4];  // J_VAC_mutual = J_SIA_mutual

    // J_He_sink: He lost to sinks
    {
        double he_sink = P.k2_disl_He * c_h;
        if (C_vac_tot > 1e-300 && Q_tot > 0.0) {
            double ell_bar_loc = Q_tot / C_vac_tot;
            for (int m = 0; m < std::min(P.v_mobile, P.V); ++m) {
                double ell_m = ell_bar_loc * std::pow(static_cast<double>(m + 1), 2.0/3.0);
                he_sink += P.k2_disl_v * ell_m * std::max(c_v[m], 0.0);
            }
        }
        dydt[P.N_eq - 1] = he_sink;
    }

    // ── Sliding-window masking for bin_moment mode ───────────────────────────
    // For the bin_moment RHS the outer cluster loops are not restructured here
    // (they are governed by discrete index + bin index, not a simple 0..N range).
    // Instead, any derivative for a state index beyond the current window
    // frontier is zeroed post-hoc.  Out-of-window concentrations are near zero
    // so their contributions to in-window derivatives are negligible.
    if (ud->window_active) {
        const int n_sia = i_d + PM * Ib;
        const int n_vac = v_d + PM * P.V_bin;

        // Zero SIA state derivatives beyond the SIA window
        const int hi_i = std::min(ud->x_hi_i, n_sia - 1);
        for (int k = hi_i + 1; k < n_sia; ++k) dydt[k] = 0.0;

        // Zero VAC state derivatives beyond the VAC window
        const int hi_v = std::min(ud->x_hi_v, n_vac - 1);
        for (int k = hi_v + 1; k < n_vac; ++k) dydt[i_VAC + k] = 0.0;
    }

    return 0;
}

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

// Forward declarations for LAPACK band and dense solvers
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

// ── Dispatch wrappers ────────────────────────────────────────────────────────

int prec_setup(sunrealtype t, N_Vector y, N_Vector fy,
               sunbooleantype jok, sunbooleantype* jcurPtr,
               sunrealtype gamma, void* user_data) {
    const UserData* ud = static_cast<const UserData*>(user_data);
    if (ud->P->prec_type == 1)
        return prec_setup_woodbury(t, y, fy, jok, jcurPtr, gamma, user_data);
    else
        return prec_setup_jacobi(t, y, fy, jok, jcurPtr, gamma, user_data);
}

int prec_solve(sunrealtype t, N_Vector y, N_Vector fy,
               N_Vector r, N_Vector z,
               sunrealtype gamma, sunrealtype delta,
               int lr, void* user_data) {
    const UserData* ud = static_cast<const UserData*>(user_data);
    if (ud->P->prec_type == 1)
        return prec_solve_woodbury(t, y, fy, r, z, gamma, delta, lr, user_data);
    else
        return prec_solve_jacobi(t, y, fy, r, z, gamma, delta, lr, user_data);
}
