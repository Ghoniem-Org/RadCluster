#!/usr/bin/env python3
"""Step 4 figures: alpha(t) schedule + 2025 comparison vs constant-alpha.

Writes (project convention: dpi=150):
  drive_doc/figs/alpha_schedule.png        the schedule with quality labels
  drive_doc/figs/alpha_t_comparison.png   2025 C_n shares: alpha(t) vs a=0.0/0.8/0.9
  <outdir>/cn_trajectory_alpha_t.png      C_n(t) 1650-2025, alpha(t) vs a=0.8
"""
import glob
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = sorted(glob.glob(os.path.join(BASE, "outputs", "*_colonial_alpha_t")))[-1]
FIGDIR = os.path.join(BASE, "drive_doc", "figs")
os.makedirs(FIGDIR, exist_ok=True)

SCHED = pd.read_csv(os.path.join(BASE, "data", "assortativity_schedule.csv"))
SUM_T = pd.read_csv(os.path.join(OUTDIR, "summary_alpha_t.csv"))
SUM_C2 = pd.read_csv(os.path.join(OUTDIR, "summary_case2_alpha_t.csv"))
BASELINE = pd.read_csv(os.path.join(
    BASE, "outputs", "20260927_161901_colonial", "summary_all.csv"))

QCOLOR = {"assumed": "#d95f02", "estimated": "#1b9e77",
          "measured-series": "#1f78b4", "interpolated": "#999999"}
QLABEL = {"assumed": "assumed", "estimated": "estimated",
          "measured-series": "measured series\n(Pew, calibrated map)",
          "interpolated": "interpolated"}

# --- 1. schedule -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 4.2))
ys = SCHED["year"].to_numpy()
aa = SCHED["alpha"].to_numpy()
# fine grid via the actual schedule class
from genealogy.assortativity import alpha_schedule
at = alpha_schedule()
tg = np.linspace(1650, 2050, 801)
ag = at(tg)
# color by quality of the nearest anchor segment
for i in range(len(ys) - 1):
    m = (tg >= ys[i]) & (tg <= ys[i + 1])
    ax.plot(tg[m], ag[m], color=QCOLOR[SCHED["quality"].iloc[i + 1]],
            lw=2.2, solid_capstyle="round")
ax.axhline(0.8, color="k", ls="--", lw=1.2, label="constant α=0.8 (central)")
ax.plot([2013], [0.70], "o", color="k", ms=6, zorder=5,
        label="α_CPS=0.70 (2013, Step 3)")
ax.axvline(1967, color="k", ls=":", lw=1.0, alpha=0.6)
ax.text(1967.5, 0.97, "Loving v. Virginia", fontsize=8, va="top")
ax.axvline(2025, color="k", ls=":", lw=1.0, alpha=0.6)
ax.text(2025.5, 0.97, "Case 2: held flat\n(assumed)", fontsize=8, va="top")
ax.set_xlim(1650, 2050)
ax.set_ylim(0.60, 1.0)
ax.set_xlabel("year")
ax.set_ylabel("α(t): pedigree-cluster assortativity")
ax.set_title("Time-varying assortativity α(t): era anchors with quality labels")
for q, lab in QLABEL.items():
    ax.plot([], [], color=QCOLOR[q], lw=2.5, label=lab)
ax.legend(loc="lower left", fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "alpha_schedule.png"), dpi=150)
print("wrote", os.path.join(FIGDIR, "alpha_schedule.png"))

# --- 2. 2025 comparison -------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
for ax, rule in zip(axes, ["shallowest_inheritance", "deepest_inheritance"]):
    b = BASELINE[BASELINE["rule"] == rule].set_index("alpha")
    r = SUM_T[SUM_T["rule"] == rule].iloc[0]
    xs = np.arange(4)
    vals = [b.loc[0.0, "Cn_share_2025"] * 100, b.loc[0.8, "Cn_share_2025"] * 100,
            b.loc[0.9, "Cn_share_2025"] * 100, r["Cn_share_2025"] * 100]
    colors = ["#999999", "#1f78b4", "#6a3d9a", "#d95f02"]
    bars = ax.bar(xs, vals, color=colors, edgecolor="k", lw=0.6)
    for x, v in zip(xs, vals):
        ax.text(x, v + 0.08, f"{v:.2f}%", ha="center", fontsize=9)
    ax.set_xticks(xs)
    ax.set_xticklabels(["α=0.0\n(null)", "α=0.8\n(central)",
                        "α=0.9\n(bound)", "α(t)\nschedule"])
    ax.set_title(f"{'Shallowest' if rule.startswith('shallow') else 'Deepest'} inheritance")
    ax.set_ylabel("2025 deepest-cluster (C₁₂) share of population")
ax = axes[0]
ax.set_ylim(0, max(ax.get_ylim()[1], axes[1].get_ylim()[1]) * 1.18)
fig.suptitle("2025 C₁₂ shares: time-varying α(t) vs constant-α baselines "
             "(colonial 1650→2025)")
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "alpha_t_comparison.png"), dpi=150)
print("wrote", os.path.join(FIGDIR, "alpha_t_comparison.png"))

# --- 3. C_n trajectories -------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, rule in zip(axes, ["shallowest_inheritance", "deepest_inheritance"]):
    dft = pd.read_csv(os.path.join(OUTDIR, f"results_alpha_t_{rule}.csv"))
    dfb = pd.read_csv(os.path.join(BASE, "outputs", "20260927_161901_colonial",
                                   f"results_{rule}.csv"))
    ax.plot(dft["year"], dft["share_C12"] * 100, color="#d95f02", lw=1.8,
            label="α(t) schedule")
    ax.plot(dfb["year"], dfb["share_C12"] * 100, color="#1f78b4", lw=1.4,
            ls="--", label="α=0.8 constant")
    ax.set_xlim(1650, 2025)
    ax.set_xlabel("year")
    ax.set_title(f"{'Shallowest' if rule.startswith('shallow') else 'Deepest'} inheritance")
    ax.legend(fontsize=9)
axes[0].set_ylabel("C₁₂ share of population (%)")
fig.suptitle("Deepest-cluster share trajectory: α(t) vs constant α=0.8")
fig.tight_layout()
out = os.path.join(OUTDIR, "cn_trajectory_alpha_t.png")
fig.savefig(out, dpi=150)
print("wrote", out)
