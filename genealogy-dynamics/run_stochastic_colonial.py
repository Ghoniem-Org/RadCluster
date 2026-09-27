#!/usr/bin/env python3
"""Step 7 — Stochastic colonial era (1650-1800): tau-leaping ensemble.

Method choice: TAU-LEAPING ensemble (task option (a)), preferred for direct
comparability with the deterministic ODE. Moment equations were not needed:
the 30-state jump process vectorizes cleanly over the trajectory axis, so a
few hundred trajectories x 150 years run in minutes.

Jump process (density-dependent, mean = the deterministic ODE):
  state X(t) in N^{2(n+1)} (PEOPLE, block order [pat, mat]) with
  - immigration:   empty -> (0,pat) at sigma(t)*I(t); -> (0,mat) at (1-sigma)*I(t)
  - births:        (q,pat)+(r,mat) -> parents (CATALYTIC, not consumed)
                   + child in (g(q,r),pat) / (g(q,r),mat), g = child rule,
                   at mating flux J_{qr}(X) (people/yr); child sex thinned
                   s : (1-s)  [Poisson thinning is exact]
  - death/emig:    (d,alpha) -> empty at (mu_alpha(t)+eps)*X_{d,alpha}
                   [exact per-leap via Binomial(X, 1-exp(-rate*tau))]

Tau selection (adaptive, Cao-Gillespie-Petzold style bound): the fastest
per-capita propensity rmax(t) = max(beta, mu) + eps sets
tau = min(tau_max, LEAP_EPS / rmax), so no propensity changes by more than
~LEAP_EPS in relative terms over a leap (J_{qr} is homogeneous degree-1 in
X with scale-free shares, so bounding relative state change bounds the
propensity change). Births/immigration are Poisson (unbounded, no
negativity); losses are Binomial (bounded by X, exact for constant
per-capita rate over the leap) so populations can never go negative.

Ensembles: shallowest/deepest x alpha in {0.8, 0.0}, N=300 trajectories,
1650 -> 1800, all-C_0 IC (50,368 people, sigma-split) — the same founding
bottleneck as the deterministic colonial run. Plus a mean-field check:
lambda=1000 population-scaled run (N=30); the ODE is exactly homogeneous of
degree 1 (J(lambda X)=lambda J(X), sources/losses linear), so the scaled
ensemble mean / lambda must converge to the ODE.

Writes outputs/<stamp>_stochastic/: yearly ensemble stats, 1800 cluster CVs,
summary table, fan/CV/trajectory/establishment figures, provenance.md.
"""
import os
import time
from datetime import datetime

import numpy as np
import pandas as pd

from genealogy.populations import (
    N_CLUSTERS, MAX_DEPTH, N_STATE, RULES, RULE_DISPLAY, state_index,
    child_cluster_matrix,
)
from genealogy.kernels import SEX_RATIO_AT_BIRTH
from genealogy.model import rhs_td
from genealogy.solver import run_hindcast, cluster_totals
from genealogy.sex_rates import load_sigma_series, load_mu_pair_series

BASE = os.path.dirname(os.path.abspath(__file__))
STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
OUTDIR = os.path.join(BASE, "outputs", f"{STAMP}_stochastic")
os.makedirs(OUTDIR, exist_ok=True)

S = SEX_RATIO_AT_BIRTH
SIGMA = load_sigma_series()
MU_PAIR = load_mu_pair_series()

inputs = pd.read_csv(os.path.join(BASE, "data", "inputs.csv"))
inputs = inputs.sort_values("year").reset_index(drop=True)
calib = {}
with open(os.path.join(BASE, "data", "emigration_calibration.txt")) as f:
    for line in f:
        line = line.strip()
        if line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        try:
            calib[k.strip()] = float(v)
        except ValueError:
            pass
EPS = calib.get("eps_calibrated", 0.00075)

years = inputs["year"].to_numpy()
I_ARR = inputs["immigration_M"].to_numpy()
BETA_ARR = 2.0 * inputs["cbr_per_1000"].to_numpy() / 1000.0


def rate_series(t):
    i = int(np.clip(np.searchsorted(years, t, side="right") - 1, 0, len(years) - 1))
    bvec = np.full(N_CLUSTERS, BETA_ARR[i])
    return I_ARR[i], bvec, bvec, MU_PAIR(t), EPS


# --- deterministic baseline (same drivers, same IC) ---------------------------
C0_M = np.zeros(N_STATE)
C0_M[state_index(0, "pat")] = SIGMA(1650.0) * 0.050368
C0_M[state_index(0, "mat")] = (1.0 - SIGMA(1650.0)) * 0.050368

ODE = {}  # (rule, alpha) -> dict(t, total_M, share)


def run_ode(rule, alpha):
    t, Y = run_hindcast(C0_M, 1650.0, 1800.0, rate_series, alpha=alpha,
                        sigma=SIGMA, s=S, rule=rule, dt=1.0)
    T = cluster_totals(Y)
    total = T.sum(axis=1)
    ODE[(rule, alpha)] = dict(t=t, total_M=total, share=T / total[:, None])


# --- tau-leaping ensemble ------------------------------------------------------
N = N_CLUSTERS
IDX_PAT = np.arange(N)
IDX_MAT = np.arange(N, 2 * N)
I0_PAT = state_index(0, "pat")
I0_MAT = state_index(0, "mat")

TAU_MAX = 0.5      # yr
LEAP_EPS = 0.002   # max relative propensity change per leap (adaptive).
                   # Convergence-tested (shallowest, alpha=0.8): the tau-leaping
                   # ensemble mean == explicit-Euler mean to MC noise, and
                   # Euler -> RK45 ODE linearly in tau. The tau error is
                   # DEPTH-AMPLIFIED (compounds through the generation
                   # chain): at LEAP_EPS=0.05 the C14 1800 share is 15% low
                   # vs the ODE, at 0.01 it is 3.2% low, at 0.002 it is
                   # 0.6% low, while totals are off by 1.6e-3/4e-4/2e-5.
                   # A too-large tau therefore mimics a "stochastic deficit"
                   # of deep clusters; 0.002 keeps the residual tau bias
                   # below the demographic-noise floor (CV ~6-11% at C14).


def mating_flux_counts(cpat, cmat, bpat, bmat, alpha):
    """J_{qr} (people/yr) for the whole ensemble. cpat/cmat: (K, n) counts."""
    wpat = bpat[None, :] * cpat
    wmat = bmat[None, :] * cmat
    Wpat = wpat.sum(axis=1)
    Wmat = wmat.sum(axis=1)
    K = cpat.shape[0]
    J = np.zeros((K, N, N))
    ok = (Wpat > 0.0) & (Wmat > 0.0)
    if not np.any(ok):
        return J
    B = np.minimum(Wpat[ok], Wmat[ok])                       # (Ko,)
    mpat = wpat[ok] / Wpat[ok, None]
    mmat = wmat[ok] / Wmat[ok, None]
    J[ok] = B[:, None, None] * (
        (1.0 - alpha) * mpat[:, :, None] * mmat[:, None, :]
        + alpha * mpat[:, :, None] * np.eye(N)[None, :, :])
    return J


def run_ensemble(rule, alpha, n_traj, seed, lam=1.0, t0=1650.0, t1=1800.0):
    """Vectorized tau-leaping ensemble. lam: population scale factor
    (mean-field check). Returns dict with yearly snapshots."""
    rng = np.random.default_rng(seed)
    K = n_traj
    child = child_cluster_matrix(rule)          # (n, n): (q, r) -> g
    gidx = child.ravel()                        # flat q*n+r -> g
    masks = [(gidx == g) for g in range(N)]

    X = np.zeros((K, 2 * N), dtype=np.int64)
    n0 = int(round(50368 * lam))
    n0pat = int(round(SIGMA(t0) * n0))
    X[:, I0_PAT] = n0pat
    X[:, I0_MAT] = n0 - n0pat

    snap_years = np.arange(int(t0), int(t1) + 1)
    n_snap = len(snap_years)
    tot_snap = np.zeros((n_snap, K))
    clu_snap = np.zeros((n_snap, K, N))
    tot_snap[0] = X.sum(axis=1)
    clu_snap[0] = X[:, :N] + X[:, N:]

    t = t0
    si = 1  # next snapshot index
    n_leaps = 0
    while t < t1 - 1e-12:
        I, bpat, bmat, (mu_pat, mu_mat), eps = rate_series(t)
        I = I * lam
        sig = SIGMA(t)
        rmax = max(float(bpat.max()), float(bmat.max()), mu_pat, mu_mat) + eps
        tau = min(TAU_MAX, LEAP_EPS / rmax, snap_years[si] - t)
        tau = max(tau, 1e-6)

        cpat = X[:, :N].astype(float)
        cmat = X[:, N:].astype(float)
        J = mating_flux_counts(cpat, cmat, bpat, bmat, alpha)   # (K,n,n) /yr

        # births by (q,r) pair type, child sex via Poisson thinning
        lam_b = J * tau
        Kpat = rng.poisson(S * lam_b).reshape(K, N * N)
        Kmat = rng.poisson((1.0 - S) * lam_b).reshape(K, N * N)

        # losses: exact per-leap Binomial(X, 1-exp(-rate*tau))
        p_pat = 1.0 - np.exp(-(mu_pat + eps) * tau)
        p_mat = 1.0 - np.exp(-(mu_mat + eps) * tau)
        Lpat = rng.binomial(X[:, :N], p_pat)
        Lmat = rng.binomial(X[:, N:], p_mat)
        X[:, :N] -= Lpat
        X[:, N:] -= Lmat

        # immigration
        X[:, I0_PAT] += rng.poisson(sig * I * 1e6 * tau, size=K)
        X[:, I0_MAT] += rng.poisson((1.0 - sig) * I * 1e6 * tau, size=K)

        # birth inflow accumulated by child cluster
        for g in range(N):
            m = masks[g]
            if np.any(m):
                X[:, g] += Kpat[:, m].sum(axis=1)
                X[:, N + g] += Kmat[:, m].sum(axis=1)

        t += tau
        n_leaps += 1
        while si < n_snap and t >= snap_years[si] - 1e-12:
            tot_snap[si] = X.sum(axis=1)
            clu_snap[si] = X[:, :N] + X[:, N:]
            si += 1
    return dict(years=snap_years, total=tot_snap, clusters=clu_snap,
                n_leaps=n_leaps, rule=rule, alpha=alpha, seed=seed, lam=lam)


def summarize(rule, alpha, ens):
    """Comparison stats vs the ODE baseline."""
    ode = ODE[(rule, alpha)]
    yrs = ens["years"]
    # align the ODE baseline to the ensemble snapshot years
    oi = np.searchsorted(np.round(ode["t"]).astype(int), yrs)
    ode_tot_M = ode["total_M"][oi]
    ode_sh = ode["share"][oi]
    tot = ens["total"]                       # (T, K) people
    clu = ens["clusters"]                    # (T, K, n) people
    K = tot.shape[1]
    ode_tot_p = ode_tot_M * 1e6              # people

    ens_mean = tot.mean(axis=1)
    ens_std = tot.std(axis=1, ddof=1)
    se = ens_std / np.sqrt(K)
    rel_dev = (ens_mean - ode_tot_p) / ode_tot_p

    share = clu / tot[:, :, None]            # (T, K, n)

    # 1800 (last snapshot) per-cluster stats
    i1800 = -1
    c1800 = clu[i1800]                       # (K, n) people
    s1800 = share[i1800]                     # (K, n)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv_share = s1800.std(axis=0, ddof=1) / s1800.mean(axis=0)
        cv_count = (c1800.std(axis=0, ddof=1)
                    / np.maximum(c1800.mean(axis=0), 1e-300))
        cv_c5 = (share[:, :, 5].std(axis=1, ddof=1)
                 / share[:, :, 5].mean(axis=1))
    ode_s1800 = ode_sh[i1800]
    ens_s1800 = s1800.mean(axis=0)
    se_s1800 = s1800.std(axis=0, ddof=1) / np.sqrt(K)
    # extinction: fraction of trajectories with zero people in C_d at 1800
    extinct = (c1800 == 0).mean(axis=0)

    # establishment year: first year T_d >= 1 (median over trajectories)
    est_year = np.full(N, np.nan)
    for d in range(N):
        hit = np.argmax(clu[:, :, d] >= 1, axis=0)   # first True index
        ever = (clu[:, :, d] >= 1).any(axis=0)
        if ever.any():
            est_year[d] = np.median(yrs[hit[ever]])

    # crossover: CV_total < 1% sustained; CV of C_5 share < 10% sustained.
    # (CV starts at 0 from the deterministic IC, peaks, then decays ~1/sqrt(N).)
    cv_tot = ens_std / ens_mean
    peak_i = int(np.argmax(cv_tot))
    cross_tot = next((y for y in range(len(yrs))
                      if np.all(cv_tot[y:] < 0.01)), None)
    cross_c5 = next((y for y in range(len(yrs))
                     if np.all(cv_c5[y:] < 0.10)), None)

    return dict(
        years=yrs, ode_tot_M=ode_tot_M, ode_share=ode_sh,
        tot=tot, share=share,
        ens_mean_total=ens_mean, ens_std_total=ens_std, se_total=se,
        rel_dev_total=rel_dev, cv_total=cv_tot,
        cv_share_1800=cv_share, cv_count_1800=cv_count,
        ode_share_1800=ode_s1800, ens_share_1800=ens_s1800, se_share_1800=se_s1800,
        extinct_1800=extinct, est_year=est_year,
        peak_cv_total=float(cv_tot[peak_i]), peak_cv_total_yr=int(yrs[peak_i]),
        cross_tot_yr=yrs[cross_tot] if cross_tot is not None else None,
        cross_c5_yr=yrs[cross_c5] if cross_c5 is not None else None,
        n_leaps=ens["n_leaps"],
    )


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t_start = time.time()
    CONFIGS = [
        ("shallowest_inheritance", 0.8, 300, 20260926, 1.0),
        ("shallowest_inheritance", 0.0, 300, 20260927, 1.0),
        ("deepest_inheritance", 0.8, 300, 20260928, 1.0),
        ("deepest_inheritance", 0.0, 300, 20260929, 1.0),
        ("shallowest_inheritance", 0.8, 30, 20260930, 1000.0),   # mean-field check
    ]
    for rule, alpha, _, _, _ in CONFIGS:
        if (rule, alpha) not in ODE:
            run_ode(rule, alpha)
            print(f"ODE done: {rule} alpha={alpha}", flush=True)

    results = {}
    for rule, alpha, n_traj, seed, lam in CONFIGS:
        tag = f"{rule}_a{alpha}" + ("_meanfield" if lam != 1.0 else "")
        print(f"ensemble start: {tag} N={n_traj} lam={lam} seed={seed}", flush=True)
        ens = run_ensemble(rule, alpha, n_traj, seed, lam=lam)
        print(f"  done: {ens['n_leaps']} leaps, "
              f"{(time.time()-t_start)/60:.1f} min elapsed", flush=True)
        if lam == 1.0:
            results[(rule, alpha)] = (ens, summarize(rule, alpha, ens))

    # ---- summary table ------------------------------------------------------
    rows = []
    for (rule, alpha), (ens, sm) in results.items():
        i1800 = -1
        ode_tot_1800 = sm["ode_tot_M"][i1800]
        m = sm["ens_mean_total"][i1800] / 1e6
        se = sm["se_total"][i1800] / 1e6
        z = (sm["ens_mean_total"][i1800] - ode_tot_1800 * 1e6) / sm["se_total"][i1800]
        rows.append(dict(
            rule=rule, alpha=alpha,
            ode_total_1800_M=ode_tot_1800,
            ens_mean_total_1800_M=m, ens_se_total_1800_M=se,
            z_total_1800=z,
            max_abs_rel_dev_total=np.abs(sm["rel_dev_total"]).max(),
            CV_total_1800=sm["cv_total"][i1800],
            CV_share_C6_1800=sm["cv_share_1800"][6],
            CV_share_C10_1800=sm["cv_share_1800"][10],
            CV_share_C14_1800=sm["cv_share_1800"][14],
            extinct_C12_1800=sm["extinct_1800"][12],
            extinct_C13_1800=sm["extinct_1800"][13],
            extinct_C14_1800=sm["extinct_1800"][14],
            est_year_C6=sm["est_year"][6], est_year_C10=sm["est_year"][10],
            peak_CV_total=sm["peak_cv_total"],
            peak_CV_total_yr=sm["peak_cv_total_yr"],
            crossover_yr_CVtot_lt_1pct=sm["cross_tot_yr"],
            crossover_yr_CV_C5share_lt_10pct=sm["cross_c5_yr"],
            n_leaps=sm["n_leaps"],
        ))
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(OUTDIR, "summary_step7.csv"), index=False)
    print(summ.to_string(index=False), flush=True)

    # ---- per-cluster 1800 CV table ------------------------------------------
    cvrows = []
    for (rule, alpha), (ens, sm) in results.items():
        for d in range(N):
            cvrows.append(dict(
                rule=rule, alpha=alpha, cluster=d,
                ode_share_1800=sm["ode_share_1800"][d],
                ens_mean_share_1800=sm["ens_share_1800"][d],
                se_share_1800=sm["se_share_1800"][d],
                CV_share_1800=sm["cv_share_1800"][d],
                CV_count_1800=sm["cv_count_1800"][d],
                extinct_frac_1800=sm["extinct_1800"][d],
                est_year=sm["est_year"][d],
            ))
    pd.DataFrame(cvrows).to_csv(os.path.join(OUTDIR, "cluster_cv_1800.csv"),
                                index=False)

    # ---- mean-field check ----------------------------------------------------
    print("mean-field check: shallowest alpha=0.8 lam=1000 N=30", flush=True)
    mf = run_ensemble("shallowest_inheritance", 0.8, 30, 20260930, lam=1000.0)
    ode = ODE[("shallowest_inheritance", 0.8)]
    mf_mean = mf["total"].mean(axis=1) / 1000.0          # back to people scale
    ode_p = ode["total_M"] * 1e6
    mf_rel = (mf_mean - ode_p) / ode_p
    # shares: compare at the 1800 snapshot only (at 1650 the ODE shares of
    # C_1.. are exactly 0, so relative deviations are undefined there)
    mf_s1800 = (mf["clusters"][-1] / mf["total"][-1][:, None]).mean(axis=0)
    ode_s1800 = ode["share"][-1]
    mask = ode_s1800 > 1e-6
    mf_share_rel_1800 = (mf_s1800[mask] - ode_s1800[mask]) / ode_s1800[mask]
    print(f"  total: max|rel dev| = {np.abs(mf_rel).max():.3e}", flush=True)
    print(f"  1800 shares (C_d with ODE share > 1e-6): max|rel dev| = "
          f"{np.abs(mf_share_rel_1800).max():.3e}", flush=True)

    # ---- figures --------------------------------------------------------------
    yrs = ODE[("shallowest_inheritance", 0.8)]["t"]

    # 1. fan chart: total population 1650-1800
    fig, ax = plt.subplots(figsize=(9, 5.2))
    sm = results[("shallowest_inheritance", 0.8)][1]
    tot_M = sm["tot"] / 1e6
    ax.fill_between(yrs, np.percentile(tot_M, 5, axis=1),
                    np.percentile(tot_M, 95, axis=1),
                    color="C0", alpha=0.25, label="ensemble 5-95% (N=300)")
    ax.plot(yrs, sm["ode_tot_M"], "k-", lw=1.6, label="deterministic ODE")
    ax.plot(yrs, sm["ens_mean_total"] / 1e6, "C0--", lw=1.2,
            label="ensemble mean")
    ax.plot(yrs, mf_mean / 1e6, "C2:", lw=1.4,
            label=r"mean-field $\lambda$=1000, N=30 (scaled)")
    if sm["cross_tot_yr"] is not None:
        ax.axvline(sm["cross_tot_yr"], color="gray", ls=":", lw=1)
        ax.text(sm["cross_tot_yr"] + 1, ax.get_ylim()[1] * 0.9,
                f"crossover {int(sm['cross_tot_yr'])}", fontsize=8,
                color="gray")
    ax.set_yscale("log")
    ax.set_xlabel("year")
    ax.set_ylabel("total population (millions)")
    ax.set_title("Stochastic colonial era: total population, 1650-1800\n"
                 "tau-leaping ensemble, shallowest inheritance, α=0.8, all-C₀ IC")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "stochastic_fan_total.png"), dpi=150)
    plt.close(fig)

    # 2. CV of 1800 cluster share vs depth
    fig, ax = plt.subplots(figsize=(9, 5.2))
    for (rule, alpha), (_, sm2) in results.items():
        if alpha == 0.8:
            ax.semilogy(range(N), sm2["cv_share_1800"], "o-", ms=4,
                        label=f"{RULE_DISPLAY[rule]}, α=0.8")
    ax.set_xlabel("cluster depth d")
    ax.set_ylabel("CV of 1800 cluster share (std/mean over trajectories)")
    ax.set_title("Demographic noise vs pedigree depth, 1800\n"
                 "CV of cluster share across N=300 tau-leaping trajectories")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "stochastic_cv_depth.png"), dpi=150)
    plt.close(fig)

    # 3. mean-vs-ODE deep-share trajectories with 2-sigma bands
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for i, rule in enumerate(RULES):
        sm3 = results[(rule, 0.8)][1]
        for j, d in enumerate((6, 10)):
            ax = axes[i, j]
            mu = sm3["share"][:, :, d].mean(axis=1)
            sd = sm3["share"][:, :, d].std(axis=1, ddof=1)
            ax.plot(yrs, sm3["ode_share"][:, d] * 100, "k-", lw=1.5,
                    label="ODE")
            ax.plot(yrs, mu * 100, "C0--", lw=1.2, label="ensemble mean")
            ax.fill_between(yrs, (mu - 2 * sd) * 100, (mu + 2 * sd) * 100,
                            color="C0", alpha=0.25, label="mean ± 2σ")
            ax.set_title(f"{RULE_DISPLAY[rule]}, C$_{d}$ share")
            ax.set_ylabel("% of population")
            ax.set_xlabel("year")
            ax.legend(fontsize=8)
    fig.suptitle("Deep-cluster shares 1650-1800: ensemble mean vs ODE "
                 "(N=300, α=0.8)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "stochastic_mean_vs_ode.png"), dpi=150)
    plt.close(fig)

    # 4. establishment curves: fraction of trajectories with T_d > 0 vs year
    fig, ax = plt.subplots(figsize=(9, 5.2))
    sm4 = results[("shallowest_inheritance", 0.8)][1]
    for d in (6, 8, 10, 12, 14):
        frac = (sm4["tot"] * 0 + (results[("shallowest_inheritance", 0.8)][0]
                                  ["clusters"][:, :, d] > 0)).mean(axis=1)
        ax.plot(sm4["years"], frac, label=f"C$_{d}$")
    ax.set_xlabel("year")
    ax.set_ylabel("fraction of trajectories with cluster occupied")
    ax.set_title("Establishment of deep clusters (shallowest, α=0.8, N=300)\n"
                 "ODE has all clusters > 0 from the first birth; "
                 "stochastic trajectories establish them discretely")
    ax.legend(fontsize=9, title="cluster")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "stochastic_establishment.png"), dpi=150)
    plt.close(fig)

    # ---- yearly stats CSV (central configs) -----------------------------------
    for (rule, alpha), (ens, sm5) in results.items():
        if alpha != 0.8:
            continue
        df = pd.DataFrame({
            "year": sm5["years"],
            "ode_total_M": sm5["ode_tot_M"],
            "ens_mean_total_M": sm5["ens_mean_total"] / 1e6,
            "ens_p5_total_M": np.percentile(sm5["tot"] / 1e6, 5, axis=1),
            "ens_p95_total_M": np.percentile(sm5["tot"] / 1e6, 95, axis=1),
            "CV_total": sm5["cv_total"],
        })
        for d in (0, 1, 2, 3, 4, 5, 6):
            df[f"ode_share_C{d}"] = sm5["ode_share"][:, d]
            df[f"ens_mean_share_C{d}"] = sm5["share"][:, :, d].mean(axis=1)
        df.to_csv(os.path.join(OUTDIR, f"yearly_stats_{rule}_a08.csv"),
                  index=False)

    with open(os.path.join(OUTDIR, "provenance.md"), "w") as f:
        f.write("# Step 7 — stochastic colonial era: provenance\n\n")
        f.write(f"Generated {STAMP}. Driver: `run_stochastic_colonial.py`.\n\n")
        f.write("Rates: identical drivers as `run_colonial.py` — "
                "`data/inputs.csv` (1650-1800), `data/immigrant_sex_share.csv` "
                "sigma(t), `data/sex_specific_mortality.csv` (mu_pat, mu_mat), "
                "eps=0.00075/yr calibrated, s=0.512 (NCHS, MEASURED), "
                "all-C_0 IC 1650 (50,368 people, sigma-split).\n")
        f.write("Method: tau-leaping on the density-dependent jump process; "
                f"adaptive tau = min({TAU_MAX}, {LEAP_EPS}/rmax) yr; "
                "Poisson births/immigration, Binomial deaths+emigration.\n")
        f.write("Ensembles: " + ", ".join(
            f"{r} α={a} N={n} seed={s}" + (f" λ={l}" if l != 1 else "")
            for r, a, n, s, l in CONFIGS) + ".\n")
        f.write("Deterministic baselines re-integrated 1650-1800 with "
                "`solver.run_hindcast` (same rate series, same IC) for "
                "apples-to-apples comparison.\n")
    print("outdir:", OUTDIR, flush=True)
    print(f"total wall time: {(time.time()-t_start)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
