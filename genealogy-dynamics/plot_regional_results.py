#!/usr/bin/env python3
"""Result figures for Step 2 (regional compartments), Case 1 1650-2025.

Reads outputs/<run>/results_regional_*.csv and writes:
- drive_doc/figs/regional_dist_2025.png : per-region 2025 cluster
  distributions at alpha=0.8 (shallowest + deepest), log scale.
- drive_doc/figs/regional_deepest_traj.png : regional deepest-cluster
  (C_12) share trajectories 1650-2025 at alpha=0.8, both rules.
"""

import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "drive_doc", "figs")
OUTDIR = os.path.join(
    HERE, "outputs", "20260927_161913_colonial_regional")

REGIONS = ["Northeast", "Midwest", "South", "West"]
RSHORT = {"Northeast": "NE", "Midwest": "MW", "South": "S", "West": "W"}
RCOLOR = {"Northeast": "#1f77b4", "Midwest": "#ff7f0e",
          "South": "#2ca02c", "West": "#d62728"}
RULES = ["shallowest_inheritance", "deepest_inheritance"]
RNICE = {"shallowest_inheritance": "shallowest", "deepest_inheritance": "deepest"}
N_CLUSTERS = 13


def load(rule, alpha):
    p = os.path.join(OUTDIR, f"results_regional_{rule}_alpha{alpha}.csv")
    return np.genfromtxt(p, delimiter=",", names=True)


def region_shares(data, region):
    return np.array([data[f"share_C{d}_{region}"] for d in range(N_CLUSTERS)])


def fig_dist_2025():
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.2), sharey=True)
    x = np.arange(N_CLUSTERS)
    for i, rule in enumerate(RULES):
        d = load(rule, 0.8)
        last = -1
        for j, r in enumerate(REGIONS):
            ax = axes[i][j]
            shares = np.array([d[f"share_C{d_}_{r}"][-1]
                               for d_ in range(N_CLUSTERS)])
            ax.bar(x, shares, color=RCOLOR[r], edgecolor="black",
                   linewidth=0.4)
            ax.set_yscale("log")
            ax.set_ylim(3e-7, 2.0)
            ax.set_xticks([0, 3, 6, 9, 12])
            if i == 1:
                ax.set_xlabel("cluster")
            if j == 0:
                ax.set_ylabel(f"{RNICE[rule]}\nshare (log)")
            ax.set_title(f"{r}", fontsize=10)
            ax.grid(True, which="major", axis="y", alpha=0.3)
    fig.suptitle("2025 per-region cluster distributions, Case 1 (1650-2025), "
                 "alpha=0.8", fontsize=12)
    fig.tight_layout()
    p = os.path.join(FIGDIR, "regional_dist_2025.png")
    fig.savefig(p, dpi=150)
    print("wrote", p)


def fig_deepest_traj():
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for i, rule in enumerate(RULES):
        ax = axes[i]
        d = load(rule, 0.8)
        for r in REGIONS:
            ax.plot(d["year"], d[f"share_C12_{r}"] * 100,
                    color=RCOLOR[r], lw=1.6, label=r)
        ax.plot(d["year"], d["Cn_share"] * 100, color="black", lw=1.2,
                ls="--", label="national agg.")
        ax.set_xlabel("year")
        if i == 0:
            ax.set_ylabel("C_12 share (%)")
        ax.set_title(f"{RNICE[rule]} bloodline, alpha=0.8")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("Regional deepest-cluster trajectories, Case 1 (1650-2025)",
                 fontsize=12)
    fig.tight_layout()
    p = os.path.join(FIGDIR, "regional_deepest_traj.png")
    fig.savefig(p, dpi=150)
    print("wrote", p)


if __name__ == "__main__":
    os.makedirs(FIGDIR, exist_ok=True)
    fig_dist_2025()
    fig_deepest_traj()
