#!/usr/bin/env python3
"""Colonial regional run: 1650 -> 2025 with Census-region compartments.

State c_{d,s,r}: 15 clusters x 2 sexes x 4 Census regions (120 ODEs) at the
n = 14 default.

- National rates: data/inputs.csv (same as run_colonial.py); mu/eps applied
  uniformly across regions (documented choice; regional vital-rate
  differences are future work).
- Immigration allocation lam_r(t): DHS Yearbook Table 4 LPR-by-state shares
  (data/regional_immigration_shares.csv); 2014-2023 measured, pre-2014 held
  at 2014 (assumed), 2024-2025 held at 2023 (interpolated).
- Inter-regional migration M(t): Census ACS 1-year state-to-state flows
  aggregated to regions (data/regional_migration_rates.csv); 2023/2024
  measured benchmarks, pre-2023 held at 2023 (assumed), 2023-2024 linearly
  interpolated, 2025 held at 2024.
- Mating: WITHIN-REGION (region assortativity = 1, assumed); cluster
  assortativity alpha applied within each region.
- IC 1650: 50,368 people in C_0, split sigma/(1-sigma) by sex and
  55% Northeast / 45% South / 0% Midwest / 0% West by region (assumed
  spin-up split: New England + Middle colonies vs Chesapeake + Carolinas;
  first ~100 yr treated as transient, as in the national run).
- alpha swept {0.0, 0.8, 0.9} x rules {shallowest, deepest}.

Writes to outputs/<stamp>_colonial_regional/.
"""
import os
from datetime import datetime

import numpy as np
import pandas as pd

from genealogy.populations import N_CLUSTERS, MAX_DEPTH, RULES, RULE_DISPLAY
from genealogy.kernels import (
    fertility_schedule, MORTALITY_RATE,
    SEX_RATIO_AT_BIRTH,
)
from genealogy.regional import (
    REGIONS, REGION_FULL, N_REGIONS, rstate_index, rhs_regional_td,
    total_balance_regional, region_totals, cluster_totals_regional,
    region_cluster_shares,
)
from genealogy.sex_rates import load_sigma_series, load_mu_pair_series
from genealogy.solver import run_hindcast
from genealogy.model import rhs_td as rhs_td_national

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "outputs",
                      f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_colonial_regional")
os.makedirs(OUTDIR, exist_ok=True)

S = SEX_RATIO_AT_BIRTH
SIGMA = load_sigma_series()  # step 5: sigma(t)
MU_PAIR = load_mu_pair_series()  # step 5: (mu_pat(t), mu_mat(t))

# --- national inputs ----------------------------------------------------------
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

# --- regional inputs ----------------------------------------------------------
lam_df = pd.read_csv(os.path.join(BASE, "data", "regional_immigration_shares.csv"))
lam_years = lam_df["year"].to_numpy()
lam_vals = lam_df[list(REGIONS)].to_numpy()          # (years, 4)

mig_df = pd.read_csv(os.path.join(BASE, "data", "regional_migration_rates.csv"))
mig_years = sorted(mig_df["year"].unique())
mig_mats = {}                                   # year -> (4,4) rate matrix
for y in mig_years:
    M = np.zeros((4, 4))
    sub = mig_df[mig_df["year"] == y]
    for _, r in sub.iterrows():
        i = REGION_FULL.index(r["from_region"])
        j = REGION_FULL.index(r["to_region"])
        M[i, j] = r["rate_per_yr"]
    mig_mats[y] = M


def lam_at(t):
    """Regional immigration shares; pre-2014 held at 2014 (assumed)."""
    i = int(np.clip(np.searchsorted(lam_years, t, side="right") - 1,
                    0, len(lam_years) - 1))
    return lam_vals[i]


def M_at(t):
    """(4,4) migration-rate matrix; linear interp between benchmarks,
    held constant outside (pre-2023 = assumed)."""
    if t <= mig_years[0]:
        return mig_mats[mig_years[0]]
    if t >= mig_years[-1]:
        return mig_mats[mig_years[-1]]
    for y0, y1 in zip(mig_years[:-1], mig_years[1:]):
        if y0 <= t < y1:
            w = (t - y0) / (y1 - y0)
            return (1 - w) * mig_mats[y0] + w * mig_mats[y1]
    return mig_mats[mig_years[-1]]


def rate_series_factory():
    years = inputs["year"].to_numpy()
    I = inputs["immigration_M"].to_numpy()
    beta = (2.0 * inputs["cbr_per_1000"].to_numpy() / 1000.0)
    def rate_series(t):
        i = int(np.clip(np.searchsorted(years, t, side="right") - 1,
                        0, len(years) - 1))
        bvec = np.full(N_CLUSTERS, beta[i])
        return I[i], lam_at(t), bvec, bvec, MU_PAIR(t), EPS, M_at(t)
    return rate_series


rate_series = rate_series_factory()

# --- initial condition: 1650 ---------------------------------------------------
IC_REGION_SHARES = np.array([0.55, 0.0, 0.45, 0.0])  # NE, MW, S, W (assumed)
c0 = np.zeros(4 * 2 * N_CLUSTERS)
for ri, share in enumerate(IC_REGION_SHARES):
    c0[rstate_index(0, "pat", ri)] = SIGMA(1650.0) * 0.050368 * share
    c0[rstate_index(0, "mat", ri)] = (1.0 - SIGMA(1650.0)) * 0.050368 * share

# --- anchors -------------------------------------------------------------------
anchors = pd.read_csv(os.path.join(BASE, "data", "validation_long.csv"))[["year", "pop_M"]]

ALPHAS = [0.0, 0.8, 0.9]
CENTRAL_ALPHA = 0.8

from functools import partial
from scipy.integrate import solve_ivp


def run_regional(c0, alpha, rule, dt=1.0, t0=1650.0, t1=2025.0,
                 rate_series_fn=None):
    rs = rate_series_fn or rate_series
    f = partial(rhs_regional_td, rate_series=rs, alpha=alpha,
                sigma=SIGMA, s=S, rule=rule)
    t_eval = np.arange(t0, t1 + dt, dt)
    sol = solve_ivp(f, (t0, t1), np.asarray(c0, dtype=float),
                    t_eval=t_eval, method="RK45", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")
    if np.any(sol.y < -1e-6):
        raise RuntimeError("negative state detected")
    return sol.t, sol.y.T


def write_rule(rule):
    rows = []
    snap_years = list(range(1650, 2026, 25))
    for alpha in ALPHAS:
        t, Y = run_regional(c0, alpha, rule)
        T = cluster_totals_regional(Y)          # (steps, n+1), all regions+sexes
        RT = region_totals(Y)                   # (steps, 4)
        RS = region_cluster_shares(Y)           # (steps, 4, n+1)
        total = T.sum(axis=1)
        share = T / total[:, None]
        # exact balance identity with migration ON: sum(rhs) ==
        # I + B - (mu+eps)*total (migration must cancel exactly)
        from genealogy.regional import rhs_regional
        dev = 0.0
        for i in range(0, len(t), 25):
            I_i, lam_i, bpat, bmat, mu_i, eps_i, M_i = rate_series(t[i])
            rsum = rhs_regional_td(t[i], Y[i], rate_series=rate_series,
                                   alpha=alpha, sigma=SIGMA, s=S,
                                   rule=rule).sum()
            _, B_i = rhs_regional(t[i], Y[i], I=I_i, lam=lam_i,
                                  beta_pat=bpat, beta_mat=bmat, alpha=alpha,
                                  mu=mu_i, eps=eps_i, M=M_i, sigma=SIGMA,
                                  s=S, rule=rule)
            Yr = Y[i].reshape(N_REGIONS, 2, N_CLUSTERS)
            pred = total_balance_regional(I_i, B_i,
                                          (Yr[:, 0, :].sum(), Yr[:, 1, :].sum()),
                                          mu_i, eps_i)
            dev = max(dev, abs(rsum - pred))
        row = dict(alpha=alpha, total_2025_M=total[-1],
                   balance_max_dev=dev)
        for ri, rname in enumerate(REGION_FULL):
            row[f"total_{rname}_2025_M"] = RT[-1, ri]
            row[f"Cn_share_{rname}_2025"] = RS[-1, ri, MAX_DEPTH]
        row["Cn_share_2025"] = share[-1, MAX_DEPTH]
        rows.append(row)

        SUM = pd.DataFrame({"year": t, "total_M": total,
                            "Cn_share": share[:, MAX_DEPTH]})
        for ri, rname in enumerate(REGION_FULL):
            SUM[f"total_{rname}_M"] = RT[:, ri]
            for d in range(N_CLUSTERS):
                SUM[f"share_C{d}_{rname}"] = RS[:, ri, d]
        for d in range(N_CLUSTERS):
            SUM[f"share_C{d}"] = share[:, d]
        SUM.to_csv(os.path.join(
            OUTDIR, f"results_regional_{rule}_alpha{alpha}.csv"), index=False)
        print(f"[{rule} a={alpha}] 2025 total={total[-1]:.1f}M "
              f"Cn={share[-1, MAX_DEPTH]*100:.2f}% "
              + " ".join(f"{REGION_FULL[ri]} Cn={RS[-1, ri, MAX_DEPTH]*100:.2f}%"
                         for ri in range(N_REGIONS))
              + f" balance={dev:.1e}", flush=True)
    SUMM = pd.DataFrame(rows)
    SUMM.to_csv(os.path.join(OUTDIR, f"summary_regional_{rule}.csv"),
                index=False)
    return SUMM


def main():
    allrows = []
    for rule in RULES:
        SUMM = write_rule(rule)
        for _, r in SUMM.iterrows():
            allrows.append(dict(rule=rule, **r.to_dict()))
    pd.DataFrame(allrows).to_csv(
        os.path.join(OUTDIR, "summary_all_regional.csv"), index=False)

    # --- anchor validation on the NATIONAL total (central alpha) ---------------
    val_lines = []
    for rule in RULES:
        df = pd.read_csv(os.path.join(
            OUTDIR, f"results_regional_{rule}_alpha0.8.csv"))
        merged = pd.merge(anchors, df[["year", "total_M"]], on="year",
                          how="left")
        merged["rel_err"] = (merged["total_M"] - merged.pop_M) / merged.pop_M
        merged.to_csv(os.path.join(
            OUTDIR, f"anchor_validation_{rule}.csv"), index=False)
        worst = merged.loc[merged["rel_err"].abs().idxmax()]
        post1960 = merged[merged["year"] >= 1960]
        val_lines.append(
            f"{rule}: worst anchor rel err {worst['rel_err']*100:+.2f}% "
            f"({int(worst['year'])}); "
            f"max |err| since 1960: {post1960['rel_err'].abs().max()*100:.2f}%")
    print("\n".join(val_lines))

    with open(os.path.join(OUTDIR, "colonial_regional.md"), "w") as f:
        f.write("# Colonial regional run 1650-2025 (4 Census regions)\n\n")
        f.write(f"{4 * 2 * N_CLUSTERS} ODE states ({N_CLUSTERS} clusters "
                f"x 2 sexes x 4 regions); "
                f"eps={EPS}/yr (calibrated); alpha swept {ALPHAS}; s={S} "
                f"(measured NCHS), sigma(t) from "
                f"data/immigrant_sex_share.csv (step 5), mu=(mu_pat,mu_mat) "
                f"from data/sex_specific_mortality.csv (step 5).\n")
        f.write("Mating within-region (region assortativity = 1, assumed); "
                "national mu/eps applied uniformly (regional vital-rate "
                "differences are future work).\n")
        f.write("Inter-regional migration: Census ACS 1-yr state-to-state "
                "flows aggregated to regions (2023/2024 measured; pre-2023 "
                "held at 2023, assumed). Immigration allocation: DHS Yearbook "
                "Table 4 LPR-by-state shares (2014-2023 measured; pre-2014 "
                "held at 2014, assumed).\n")
        f.write("IC 1650: 50,368 people in C_0; 55% Northeast / 45% South / "
                "0% Midwest / 0% West (assumed spin-up split).\n\n")
        f.write("## 2025 deepest-cluster (C_n) share by region, alpha=0.8\n\n")
        f.write("| rule | " + " | ".join(REGION_FULL) + " | national |\n")
        f.write("|---|---|---|---|---|---|\n")
        summ = pd.read_csv(os.path.join(OUTDIR, "summary_all_regional.csv"))
        for rule in RULES:
            r = summ[(summ.rule == rule)
                     & (summ.alpha == CENTRAL_ALPHA)].iloc[0]
            cells = [f"{r[f'Cn_share_{rn}_2025']*100:.2f}%"
                     for rn in REGION_FULL]
            f.write(f"| {RULE_DISPLAY[rule]} | " + " | ".join(cells)
                    + f" | {r['Cn_share_2025']*100:.2f}% |\n")
        f.write("\n## Anchor validation (national total)\n\n"
                + "\n".join(f"- {l}" for l in val_lines) + "\n")
    print("outdir:", OUTDIR)


if __name__ == "__main__":
    main()
