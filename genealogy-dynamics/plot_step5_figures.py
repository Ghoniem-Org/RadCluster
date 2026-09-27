#!/usr/bin/env python3
"""Step-5 figures.

1. sigma(t) 1650-2025 with era labels (data/immigrant_sex_share.csv).
2. 2025 deepest-cluster share: old (sigma=0.5, uniform mu) vs new
   (sigma(t), mu_pat/mu_mat), two-sex model, shallowest/deepest x 3 alphas.

Usage: plot_step5_figures.py <new_colonial_dir> <outdir>
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NEW_DIR, OUTDIR = sys.argv[1], sys.argv[2]
BASE = os.path.dirname(os.path.abspath(__file__))
OLD_DIR = os.path.join(BASE, "outputs", "20260927_161901_colonial")
os.makedirs(OUTDIR, exist_ok=True)


def fig_sigma():
    df = pd.read_csv(os.path.join(BASE, "data", "immigrant_sex_share.csv"))
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.plot(df["year"], df["sigma"], lw=1.4, color="tab:blue")
    for y0, y1, lab in [(1650, 1819, "assumed 0.65\n(indentured-labor skew)"),
                        (1820, 1917, "estimated 0.629\n(HSUS period aggregate\napplied annually)"),
                        (1918, 1949, "estimated\nquota-era convergence"),
                        (1950, 2014, "estimated/\ninterpolated"),
                        (2015, 2024, "measured\n(DHS Table 9) +\nestimated unauth.\nstock-as-inflow proxy")]:
        ax.axvspan(y0, y1, alpha=0.08, color="tab:blue")
        ax.text((y0 + y1) / 2, 0.66, lab, ha="center", va="top", fontsize=7,
                color="tab:blue")
    ax.axhline(0.5, ls="--", lw=1, color="k", alpha=0.5, label="old sigma=0.5")
    ax.set_xlim(1650, 2025)
    ax.set_xlabel("year")
    ax.set_ylabel("immigrant male share  sigma(t)")
    ax.set_title("Step 5: data-driven immigrant male share")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "step5_sigma_series.png"), dpi=150)


def fig_compare():
    old = pd.read_csv(os.path.join(OLD_DIR, "summary_all.csv"))
    new = pd.read_csv(os.path.join(NEW_DIR, "summary_all.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, rule in zip(axes, ["shallowest_inheritance", "deepest_inheritance"]):
        o = old[old["rule"] == rule].sort_values("alpha")
        n = new[new["rule"] == rule].sort_values("alpha")
        x = np.arange(len(o))
        w = 0.35
        ax.bar(x - w / 2, o["Cn_share_2025"] * 100, w, label="old: sigma=0.5, uniform mu")
        ax.bar(x + w / 2, n["Cn_share_2025"] * 100, w, label="new: sigma(t), mu_pat/mu_mat")
        ax.set_xticks(x)
        ax.set_xticklabels([f"alpha={a}" for a in o["alpha"]])
        ax.set_title(rule.replace("_", " "))
        ax.set_ylabel("2025 deepest-cluster share (%)")
        ax.legend(fontsize=8)
    fig.suptitle("Step 5: 2025 deepest-cluster share, old vs new sex-specific rates")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "step5_old_vs_new_cn.png"), dpi=150)


if __name__ == "__main__":
    fig_sigma()
    fig_compare()
    print("wrote", OUTDIR)
