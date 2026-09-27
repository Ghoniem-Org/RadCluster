#!/usr/bin/env python3
"""Per-age-group cluster-distribution panels for the colonial run.

Runs the age-structured colonial model (1650->2025) at alpha=0.8 for the
shallowest and deepest rules, saves the 2025 age-resolved cluster distribution,
and draws a 2x4 panel figure: rows = shallowest/deepest, columns = ages 0-14 /
15-49 / 50+ / all ages. Each panel shows the 2025 pedigree-cluster
distribution (C_0..C_14 shares) for that age group.

Outputs:
  outputs/age_step1/cluster_by_age_2025.png
  outputs/age_step1/colonial_age/cluster_by_age_2025.csv  (tidy data)
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import run_colonial_age as R
from genealogy import solver
from genealogy.age import N_CLUSTERS, N_AGE, initial_condition_age
from genealogy.populations import RULES
from genealogy.kernels import SEX_RATIO_AT_BIRTH
from genealogy.sex_rates import load_sigma_series

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "outputs", "age_step1")
AGEDIR = os.path.join(OUTDIR, "colonial_age")
os.makedirs(OUTDIR, exist_ok=True)

# sigma(t): same time-varying immigrant male-share series as run_colonial_age.py
# (Step 5 data-driven series); a constant here would silently change results.
SIGMA = load_sigma_series()
S = SEX_RATIO_AT_BIRTH
NCAP = N_CLUSTERS - 1  # cap cluster index (C_14 at default MAX_DEPTH=14)

# age-group definitions (band indices)
GROUPS = {
    "0-14": (0, 1, 2),
    "15-49": tuple(range(3, 10)),
    "50+": tuple(range(10, 18)),
    "all ages": tuple(range(18)),
}
ALPHA = 0.8


def main():
    series = R.AgeSeries()
    _, beta0, mu_p0, mu_m0, _ = series.at_year(1650)
    c0 = initial_condition_age(0.050368, beta0, mu_p0, mu_m0,
                               sigma=SIGMA(1650.0), s=S)

    rows = []
    dist = {}
    for rule in RULES:
        print(f"running {rule} ...", flush=True)
        t, Y = solver.run_age_hindcast(c0, 1650, 2025, series, ALPHA,
                                       SIGMA, S, rule, dt=1.0)
        C = Y[-1].reshape(2, N_CLUSTERS, N_AGE)  # [sex, cluster, age]
        dist[rule] = {}
        for gname, bands in GROUPS.items():
            tot = C[:, :, list(bands)].sum(axis=(0, 2))  # per cluster
            share = tot / tot.sum()
            dist[rule][gname] = share
            for d in range(N_CLUSTERS):
                rows.append({"rule": rule, "age_group": gname,
                             "cluster": d, "share": share[d]})
    pd.DataFrame(rows).to_csv(
        os.path.join(AGEDIR, "cluster_by_age_2025.csv"), index=False)

    fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharex=True, sharey=True)
    for i, rule in enumerate(RULES):
        for j, gname in enumerate(GROUPS):
            ax = axes[i][j]
            s = dist[rule][gname]
            ax.bar(range(N_CLUSTERS), s * 100, color="steelblue")
            ax.set_title(f"{gname}", fontsize=10)
            ax.set_xticks(range(0, N_CLUSTERS, 2))
            ax.grid(True, alpha=0.3, axis="y")
            # annotate cap cluster
            ax.text(NCAP, s[NCAP] * 100 * 1.02, f"{s[NCAP]*100:.1f}%",
                    ha="center", va="bottom", fontsize=8, color="darkred")
        axes[i][0].set_ylabel(f"{rule.replace('_', ' ').title()}\nshare (%)", fontsize=10)
    for j in range(4):
        axes[1][j].set_xlabel("pedigree cluster d")
    fig.suptitle("2025 pedigree-cluster distribution by age group "
                 f"(age-structured colonial run, alpha={ALPHA})")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "cluster_by_age_2025.png"), dpi=110)
    print("wrote", os.path.join(OUTDIR, "cluster_by_age_2025.png"))

    # summary for the log
    for rule in RULES:
        ccap = {g: dist[rule][g][NCAP] * 100 for g in GROUPS}
        print(f"{rule} C_{NCAP} by age: " +
              ", ".join(f"{g}={v:.2f}%" for g, v in ccap.items()))


if __name__ == "__main__":
    main()
