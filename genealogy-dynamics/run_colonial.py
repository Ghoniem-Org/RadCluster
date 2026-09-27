#!/usr/bin/env python3
"""Colonial run: 1650 -> 2025 with anchored historical rates, both rules.

- Rates: data/inputs.csv (1650-2024; 2025 holds the 2024 values).
- IC 1650: entire colonial population (50,368 people) in cluster C_0, split
  sigma/(1-sigma) between P_pat and P_mat; C_1..C_n = 0.
  This is a SPIN-UP approximation, not a claim about 1650 pedigrees:
  treat roughly the first 100 years as transient.
- Emigration eps: calibrated (data/emigration_calibration.txt).
- Assortativity alpha = 0.8 (enforced central assumption, matching the mean of
  the religious-group endogamy estimates; no annual series exists).
- Sex ratio at birth s = 0.512 (measured, NCHS); immigrant male fraction
  sigma = 0.5 (assumed).

Writes to outputs/<stamp>_colonial/: per-rule results (cluster totals =
paternal + maternal), per-population snapshot CSVs and heatmaps, anchor
validation, and colonial.md reading guide.
"""
import os
from datetime import datetime

import numpy as np
import pandas as pd

from genealogy.populations import (
    N_CLUSTERS, MAX_DEPTH, N_STATE, RULES, RULE_DISPLAY, state_index,
)
from genealogy.kernels import (
    fertility_schedule, MORTALITY_RATE,
    SEX_RATIO_AT_BIRTH,
)
from genealogy.model import rhs_td
from genealogy.solver import run_hindcast, cluster_totals, sex_totals
from genealogy.visualization import plot_pedigree_heatmap

from genealogy.sex_rates import load_sigma_series, load_mu_pair_series

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "outputs",
                      f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_colonial")
os.makedirs(OUTDIR, exist_ok=True)

S = SEX_RATIO_AT_BIRTH
SIGMA = load_sigma_series()            # step 5: sigma(t) immigrant male share
MU_PAIR = load_mu_pair_series()        # step 5: (mu_pat(t), mu_mat(t))

# --- inputs -----------------------------------------------------------------
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


def rate_series_factory():
    years = inputs["year"].to_numpy()
    I = inputs["immigration_M"].to_numpy()
    beta = (2.0 * inputs["cbr_per_1000"].to_numpy() / 1000.0)

    def rate_series(t):
        i = int(np.clip(np.searchsorted(years, t, side="right") - 1, 0, len(years) - 1))
        bvec = np.full(N_CLUSTERS, beta[i])
        return I[i], bvec, bvec, MU_PAIR(t), EPS
    return rate_series


rate_series = rate_series_factory()

# --- initial condition: 1650 -------------------------------------------------
c0 = np.zeros(N_STATE)
c0[state_index(0, "pat")] = SIGMA(1650.0) * 0.050368  # 50,368 people, M/F split by sigma
c0[state_index(0, "mat")] = (1.0 - SIGMA(1650.0)) * 0.050368

# --- anchors -----------------------------------------------------------------
anchors = pd.read_csv(os.path.join(BASE, "data", "validation_long.csv"))
anchors = anchors.rename(columns={"pop_M": "pop_M"})[["year", "pop_M"]]

ALPHAS = [0.0, 0.8, 0.9]
CENTRAL_ALPHA = 0.8

SUMMARY = {}


def run_rule(rule, alpha):
    t, Y = run_hindcast(c0, 1650.0, 2025.0, rate_series, alpha=alpha,
                        sigma=SIGMA, s=S, rule=rule, dt=1.0)
    T = cluster_totals(Y)                    # (steps, n+1) cluster totals
    Tpat, Tmat = sex_totals(Y)
    total = T.sum(axis=1)
    share = T / total[:, None]
    return t, Y, T, Tpat, Tmat, total, share


def write_rule(rule):
    rows, snaps = [], []
    snap_years = list(range(1650, 2026, 25))
    heatmaps = {}
    for alpha in ALPHAS:
        t, Y, T, Tpat, Tmat, total, share = run_rule(rule, alpha)
        # exact balance identity: sum(rhs) == I + B - (mu+eps)*total_pat/mat
        # (verifies the mating flux conserves people in the td driver)
        from genealogy.model import total_balance as _tb
        from genealogy.kernels import mating_flux as _mf
        dev = 0.0
        for i in range(0, len(t), 25):
            I_i, bpat, bmat, mu_i, eps_i = rate_series(t[i])
            rsum = rhs_td(t[i], Y[i], rate_series=rate_series, alpha=alpha,
                          sigma=SIGMA, s=S, rule=rule).sum()
            _, B_i = _mf(Y[i], bpat, bmat, alpha, rule)
            pred = _tb(I_i, B_i, Y[i, :N_CLUSTERS].sum(),
                       Y[i, N_CLUSTERS:].sum(), mu_i, eps_i)
            dev = max(dev, abs(rsum - pred))
        rows.append(dict(alpha=alpha, total_2025_M=total[-1],
                         pat_2025_M=Tpat[-1], mat_2025_M=Tmat[-1],
                         pat_share_2025=Tpat[-1] / total[-1],
                         Cn_share_2025=share[-1, MAX_DEPTH],
                         balance_max_dev=dev))
        SUM = pd.DataFrame({"year": t, "total_M": total,
                            "total_pat_M": Tpat, "total_mat_M": Tmat})
        for d in range(N_CLUSTERS):
            SUM[f"share_C{d}_pat"] = Y[:, state_index(d, "pat")] / total
            SUM[f"share_C{d}_mat"] = Y[:, state_index(d, "mat")] / total
            SUM[f"share_C{d}"] = share[:, d]
        SUM.to_csv(os.path.join(OUTDIR, f"results_{rule}_alpha{alpha}.csv"),
                   index=False)
        if alpha == CENTRAL_ALPHA:
            SUM.to_csv(os.path.join(OUTDIR, f"results_{rule}.csv"), index=False)
            idx = [int(np.argmin(np.abs(t - yy))) for yy in snap_years]
            snap = pd.DataFrame({"year": t[idx]})
            for d in range(N_CLUSTERS):
                snap[f"share_C{d}_pat"] = Y[idx, state_index(d, "pat")] / total[idx]
                snap[f"share_C{d}_mat"] = Y[idx, state_index(d, "mat")] / total[idx]
                snap[f"share_C{d}"] = share[idx, d]
            snap.to_csv(os.path.join(OUTDIR, f"snapshots_{rule}.csv"), index=False)
            Zpat = np.array([Y[idx, state_index(d, "pat")] / total[idx]
                             for d in range(N_CLUSTERS)])
            Zmat = np.array([Y[idx, state_index(d, "mat")] / total[idx]
                             for d in range(N_CLUSTERS)])
            for pop, Z in (("pat", Zpat), ("mat", Zmat)):
                out = os.path.join(OUTDIR, f"heatmap_{rule}_{pop}.png")
                plot_pedigree_heatmap(
                    np.array(snap_years), Z, out,
                    title=(f"Population cluster distribution, 1650-2025 "
                           f"[{RULE_DISPLAY[rule]} rule]"),
                    subtitle=(f"eps={EPS}/yr, α=0.8 (central), s={S} · "
                              f"{'paternal' if pop == 'pat' else 'maternal'} population"),
                    vmin=0.0, vmax=0.40, population=pop)
                heatmaps[pop] = out
                print("wrote", out)
    SUMM = pd.DataFrame(rows)
    SUMM.to_csv(os.path.join(OUTDIR, f"summary_{rule}.csv"), index=False)
    return SUMM, heatmaps


def main():
    allrows = []
    for rule in RULES:
        SUMM, heatmaps = write_rule(rule)
        for _, r in SUMM.iterrows():
            allrows.append(dict(rule=rule, **r.to_dict()))
        print(SUMM.to_string(index=False))
    pd.DataFrame(allrows).to_csv(os.path.join(OUTDIR, "summary_all.csv"),
                                 index=False)

    # --- anchor validation (central-alpha runs) --------------------------------
    val_lines = []
    for rule in RULES:
        df = pd.read_csv(os.path.join(OUTDIR, f"results_{rule}.csv"))
        merged = pd.merge(anchors, df[["year", "total_M"]], on="year", how="left")
        merged["rel_err"] = (merged["total_M"] - merged.pop_M) / merged.pop_M
        merged.to_csv(os.path.join(OUTDIR, f"anchor_validation_{rule}.csv"),
                      index=False)
        worst = merged.loc[merged["rel_err"].abs().idxmax()]
        post1960 = merged[merged["year"] >= 1960]
        val_lines.append(
            f"{rule}: worst anchor rel err {worst['rel_err']*100:+.2f}% "
            f"({int(worst['year'])}); "
            f"max |err| since 1960: {post1960['rel_err'].abs().max()*100:.2f}%")
    print("\n".join(val_lines))

    with open(os.path.join(OUTDIR, "colonial.md"), "w") as f:
        f.write("# Colonial run 1650-2025 (two-sex cluster model)\n\n")
        f.write("Emigration eps={EPS}/yr (calibrated); alpha swept "
                f"{ALPHAS}; s={S} (measured NCHS sex ratio at birth). "
                f"sigma(t) = immigrant male share from "
                "data/immigrant_sex_share.csv (step 5). mu(t) = "
                "(mu_pat(t), mu_mat(t)) from "
                "data/sex_specific_mortality.csv (step 5).\n\n")
        f.write("results_<rule>.csv: alpha=0.8 (central), per-cluster total shares "
                "(share_C{d}) plus paternal/maternal shares (share_C{d}_pat/_mat).\n")
        f.write("results_<rule>_alpha{a}.csv: alpha=0.0/0.8/0.9 variants.\n")
        f.write("snapshots_<rule>.csv: 25-year snapshots (alpha=0.8).\n")
        f.write("heatmap_<rule>_pat.png / heatmap_<rule>_mat.png: per-population "
                "cluster distributions.\n\n")
        f.write("## 2025 deepest-cluster (C_n) total share\n\n")
        f.write("| rule | alpha=0.0 | **0.8 (central)** | 0.9 |\n|---|---|---|---|\n")
        summ = pd.read_csv(os.path.join(OUTDIR, "summary_all.csv"))
        for rule in RULES:
            r = summ[summ.rule == rule].set_index("alpha")
            f.write(f"| {RULE_DISPLAY[rule]} | "
                    + " | ".join(f"{r.loc[a, 'Cn_share_2025']*100:.2f}%"
                                 for a in ALPHAS) + " |\n")
        f.write("\n## 2025 paternal share of total population\n\n")
        for rule in RULES:
            r = summ[summ.rule == rule].set_index("alpha")
            f.write(f"| {RULE_DISPLAY[rule]} | "
                    + " | ".join(f"{r.loc[a, 'pat_share_2025']*100:.2f}%"
                                 for a in ALPHAS) + " |\n")
        f.write("\n## Anchor validation\n\n" + "\n".join(f"- {l}" for l in val_lines) + "\n")
    print("outdir:", OUTDIR)


if __name__ == "__main__":
    main()
