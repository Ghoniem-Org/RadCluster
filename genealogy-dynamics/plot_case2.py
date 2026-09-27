#!/usr/bin/env python3
"""Case-2 (2025-2050) publication figures -> drive_doc/figs/ (proj2050_* names).

Reads the newest outputs/*_case2 and outputs/*_case2_religious directories.
Figures (dpi=150, matching project convention):
  proj2050_cluster_evolution_{shallowest,deepest}.png -- 2025-2050 per-cluster
      share evolution, two-sex global, alpha=0.8 (central)
  proj2050_cn_by_rule.png -- 2050 deepest-cluster shares by rule x alpha
      (table-ready numbers; the CSV is summary_case2.csv)
  proj2050_religious_shares.png -- religious composition 2025-2050 (global),
      deepest rule, alpha=0.8 (central)
  proj2050_religious_groups.png -- per-group panels: each group's share of
      total population 2025-2050, deepest rule, alpha=0.8 (central)
"""
import glob
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from genealogy.populations import N_CLUSTERS, MAX_DEPTH, RULE_DISPLAY

BASE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(BASE, "drive_doc", "figs")


def _newest(pattern):
    cands = sorted(glob.glob(os.path.join(BASE, "outputs", pattern)))
    if not cands:
        raise RuntimeError(f"no outputs/{pattern} directory found")
    return cands[-1]


CASE2 = _newest("*_case2")
CASE2R = _newest("*_case2_religious")
print("two-sex dir:", CASE2)
print("religious dir:", CASE2R)

CLASS_COLORS = plt.get_cmap("tab20").colors
GROUPS = ["pro", "cat", "mor", "jew", "mus", "dhr", "oth", "una"]
GD = {"pro": "Protestant", "cat": "Catholic", "mor": "Mormon (LDS)",
      "jew": "Jewish", "mus": "Muslim", "dhr": "Hindu/Buddhist",
      "oth": "Orthodox/other", "una": "Unaffiliated"}
COLORS = {"pro": "#1f77b4", "cat": "#ff7f0e", "mor": "#2ca02c",
          "jew": "#d62728", "mus": "#9467bd", "dhr": "#8c564b",
          "oth": "#e377c2", "una": "#7f7f7f"}

# --- 1. two-sex cluster-share evolution 2025-2050 (alpha=0.8 central) --------
for rule in ("shallowest_inheritance", "deepest_inheritance"):
    df = pd.read_csv(
        os.path.join(CASE2, f"results_case2_{rule}_alpha0.8.csv"))
    years = df["year"].to_numpy()
    fig, ax = plt.subplots(figsize=(13, 6.8))
    for d in range(N_CLUSTERS):
        y = (df[f"share_C{d}_pat"].to_numpy()
             + df[f"share_C{d}_mat"].to_numpy()) * 100.0
        ax.plot(years, y, color=CLASS_COLORS[d % 20],
                lw=2.4 if d == 0 else 1.4,
                label="C0 foreign-born" if d == 0
                else (f"C{d} ({MAX_DEPTH}+)" if d == MAX_DEPTH else f"C{d}"),
                zorder=3)
    ax.axvline(2025, color="black", lw=1.2, ls="--", alpha=0.6)
    ax.text(2025.2, 96.5, "Case-1 2025 handoff", fontsize=8, va="top",
            color="#555555")
    ax.set_xlim(2025, 2050)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Year")
    ax.set_ylabel("Share of population (%)")
    disp = RULE_DISPLAY[rule]
    ax.set_title(f"Case 2 projection 2025-2050: cluster shares [{disp} rule]\n"
                 "projected inputs (see outputs/*_case2/case2_inputs.md), "
                 "\u03b1=0.8 (enforced central endogamy)", fontsize=11)
    ax.grid(True, alpha=0.2)
    ax.legend(title="Cluster (total population)", fontsize=8,
              title_fontsize=9, loc="upper left",
              bbox_to_anchor=(1.01, 1.0), frameon=True)
    fig.tight_layout()
    out = os.path.join(FIGDIR, f"proj2050_cluster_evolution_{rule}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)

# --- 2. 2050 deepest-cluster shares by rule x alpha ---------------------------
summ = pd.read_csv(os.path.join(CASE2, "summary_case2.csv"))
rules = ["shallowest_inheritance", "deepest_inheritance"]
alphas = [0.0, 0.8, 0.9]
x = np.arange(len(alphas))
w = 0.30
fig, ax = plt.subplots(figsize=(10.5, 5.2))
for i, rule in enumerate(rules):
    vals = [100 * float(summ[(summ.rule == rule)
                             & (summ.alpha == a)]["Cn_share_2050"].iloc[0])
            for a in alphas]
    ax.bar(x + (i - 0.5) * w, vals, w, label=RULE_DISPLAY[rule].capitalize(),
           color=plt.cm.Set2(i))
    for j, v in enumerate(vals):
        ax.text(x[j] + (i - 0.5) * w, v + 0.15, f"{v:.2f}%", ha="center",
                fontsize=8)
ax.set_xticks(x)
ax.set_xticklabels(["0.0 (null)", "0.8 (central)", "0.9 (bound)"])
ax.set_xlabel("Assortativity \u03b1")
ax.set_ylabel(f"Deepest-cluster (C{MAX_DEPTH}) share of total 2050 population (%)")
ax.set_title("Case 2: deepest-cluster share in 2050, by bloodline rule and "
             "\u03b1\n2025-2050 projection from the Case-1 2025 state")
ax.legend()
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
out = os.path.join(FIGDIR, "proj2050_cn_by_rule.png")
fig.savefig(out, dpi=150)
plt.close(fig)
print("wrote", out)

# --- 3. religious composition 2025-2050 (global; deepest, alpha=0.8) ----------
df = pd.read_csv(os.path.join(
    CASE2R, "results_case2_religious_deepest_inheritance_alpha0.8.csv"))
t = df["year"].to_numpy()
fig, ax = plt.subplots(figsize=(10, 5.5))
bottom = np.zeros(len(t))
for g in GROUPS:
    y = 100 * df[f"group_{g}_share"].to_numpy()
    ax.fill_between(t, bottom, bottom + y, color=COLORS[g], label=GD[g],
                    alpha=0.85)
    bottom += y
ax.axvline(2025, color="black", lw=1.2, ls="--", alpha=0.6)
ax.text(2025.2, 96, "Case-1 2025 handoff", fontsize=8, va="top",
        color="#222222")
ax.set_xlim(2025, 2050)
ax.set_ylim(0, 100)
ax.set_xlabel("Year")
ax.set_ylabel("Share of modeled population (%)")
ax.set_title("Case 2 projection 2025-2050: religious composition "
             "[deepest rule, \u03b1=0.8 (central)]\ndisaffiliation continued "
             "at calibrated \u03b4=8\u00d710\u207b\u00b3/yr")
ax.legend(ncol=4, fontsize=9, loc="center left")
ax.grid(alpha=0.3)
fig.tight_layout()
out = os.path.join(FIGDIR, "proj2050_religious_shares.png")
fig.savefig(out, dpi=150)
plt.close(fig)
print("wrote", out)

# --- 4. per-group panels: group share of total 2025-2050 ----------------------
fig, axes = plt.subplots(2, 4, figsize=(13, 6.4), sharex=True, sharey=True)
for ax, g in zip(axes.flat, GROUPS):
    y = 100 * df[f"group_{g}_share"].to_numpy()
    ax.plot(t, y, color=COLORS[g], lw=2)
    ax.set_title(GD[g], fontsize=10, color=COLORS[g])
    ax.grid(alpha=0.3)
    ax.annotate(f"{y[0]:.1f}%\u2192{y[-1]:.1f}%",
                xy=(0.05, 0.88), xycoords="axes fraction", fontsize=9)
    ax.set_xlim(2025, 2050)
fig.suptitle("Case 2: each group's share of the total population, 2025-2050 "
             "[deepest rule, \u03b1=0.8 (central)]", fontsize=11)
fig.text(0.5, 0.02, "Year", ha="center", fontsize=10)
fig.text(0.01, 0.5, "Share of total population (%)", va="center",
         rotation="vertical", fontsize=10)
fig.tight_layout(rect=[0.03, 0.05, 1, 0.94])
out = os.path.join(FIGDIR, "proj2050_religious_groups.png")
fig.savefig(out, dpi=150)
plt.close(fig)
print("wrote", out)
print("done")
