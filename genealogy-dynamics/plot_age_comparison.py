#!/usr/bin/env python3
"""Step 1 comparison figures: age-structured vs aggregate colonial runs.

- pyramid_2025.png : age pyramid in 2025 (age model, shallowest alpha=0.8)
  by sex, compared with the stable pyramid at 2024 rates.
- cluster_dist_2025.png : 2025 cluster distributions, age vs aggregate,
  shallowest and deepest at alpha=0.8.
- wavefront_colonial.png : C_d share trajectories 1650-2025 (age vs agg).
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from genealogy.populations import N_CLUSTERS

BASE = os.path.dirname(os.path.abspath(__file__))
AGEDIR = os.path.join(BASE, "outputs", "age_step1", "colonial_age")
AGGDIR = os.path.join(BASE, "outputs", "20260927_161901_colonial")
OUTDIR = os.path.join(BASE, "outputs", "age_step1")
os.makedirs(OUTDIR, exist_ok=True)

NCAP = N_CLUSTERS - 1  # cap cluster index (C_14 at default MAX_DEPTH=14)
BANDS = ["0_4", "5_9", "10_14", "15_19", "20_24", "25_29", "30_34",
         "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69",
         "70_74", "75_79", "80_84", "85p"]
MIDS = [2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5,
        57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 90.0]

# --- 2025 cluster distributions: age vs aggregate ---------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
for ax, rule, title in zip(axes, ("shallowest_inheritance", "deepest_inheritance"),
                           ("Shallowest inheritance", "Deepest inheritance")):
    age = pd.read_csv(os.path.join(AGEDIR, f"results_{rule}_alpha0.8.csv"))
    agg = pd.read_csv(os.path.join(AGGDIR, f"results_{rule}_alpha0.8.csv"))
    a = age[age.year == 2025].iloc[0]
    g = agg[agg.year == 2025].iloc[0]
    d = np.arange(N_CLUSTERS)
    ax.bar(d - 0.2, [a[f"share_C{i}"] * 100 for i in d], 0.4,
           label=f"age-structured ({a['total_M']:.0f}M)")
    ax.bar(d + 0.2, [g[f"share_C{i}"] * 100 for i in d], 0.4,
           label=f"aggregate ({g['total_M']:.0f}M)")
    ax.set_xlabel("pedigree cluster d")
    ax.set_title(title)
    ax.set_xticks(d)
    ax.grid(True, alpha=0.3, axis="y")
axes[0].set_ylabel("share of 2025 population (%)")
axes[0].legend(fontsize=8)
fig.suptitle("2025 pedigree-cluster distribution at alpha=0.8: "
             "age-structured vs aggregate")
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "cluster_dist_2025.png"), dpi=110)

# --- C_n share trajectories --------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
for ax, rule, title in zip(axes, ("shallowest_inheritance", "deepest_inheritance"),
                           ("Shallowest inheritance", "Deepest inheritance")):
    age = pd.read_csv(os.path.join(AGEDIR, f"results_{rule}_alpha0.8.csv"))
    agg = pd.read_csv(os.path.join(AGGDIR, f"results_{rule}_alpha0.8.csv"))
    ax.plot(age["year"], age[f"share_C{NCAP}"] * 100, label="age-structured")
    ax.plot(agg["year"], agg[f"share_C{NCAP}"] * 100, label="aggregate")
    ax.set_xlabel("year")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
axes[0].set_ylabel(f"C_{NCAP} share (%)")
axes[0].legend(fontsize=8)
fig.suptitle(f"Deepest-cluster (C_{NCAP}) share 1650-2025 at alpha=0.8")
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "cn_trajectory.png"), dpi=110)

print("comparison figures written to", OUTDIR)
