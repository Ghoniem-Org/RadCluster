#!/usr/bin/env python3
"""Figures for the religious-dimension colonial runs.

Reads outputs/<stamp>_colonial_religious/results_religious_<rule>_alpha<a>.csv.
Writes into <indir>/figures/.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
INDIR = os.path.join(BASE, "outputs", "20260927_161912_colonial_religious")
FIGDIR = os.path.join(INDIR, "figures")
os.makedirs(FIGDIR, exist_ok=True)

GROUPS = ["pro", "cat", "mor", "jew", "mus", "dhr", "oth", "una"]
GD = {"pro": "Protestant", "cat": "Catholic", "mor": "Mormon (LDS)",
      "jew": "Jewish", "mus": "Muslim", "dhr": "Hindu/Buddhist",
      "oth": "Orthodox/other", "una": "Unaffiliated"}
COLORS = {"pro": "#1f77b4", "cat": "#ff7f0e", "mor": "#2ca02c",
          "jew": "#d62728", "mus": "#9467bd", "dhr": "#8c564b",
          "oth": "#e377c2", "una": "#7f7f7f"}
RULES = ["shallowest_inheritance", "deepest_inheritance"]
RD = {"shallowest_inheritance": "Strict",
      "deepest_inheritance": "Lenient"}
ALPHA = 0.8

# Pew 2023-24 Religious Landscape Study anchors (adults; ESTIMATED).
PEW = {"pro": 0.40, "cat": 0.19, "mor": 0.016, "jew": 0.02,
       "mus": 0.01, "dhr": 0.02, "oth": 0.01, "una": 0.29}

dfs = {r: pd.read_csv(os.path.join(
    INDIR, f"results_religious_{r}_alpha{ALPHA}.csv")) for r in RULES}
df = dfs["deepest_inheritance"]
t = df["year"].to_numpy()

# --- 1. religious-share time series + Pew anchors ------------------------------
fig, ax = plt.subplots(figsize=(10, 5.5))
for g in GROUPS:
    ax.plot(t, 100 * df[f"group_{g}_share"], color=COLORS[g], lw=1.8,
            label=GD[g])
    ax.plot(2024, 100 * PEW[g], marker="x", ms=8, mew=2, color=COLORS[g])
ax.set_xlim(1650, 2025)
ax.set_xlabel("Year")
ax.set_ylabel("Share of modeled population (%)")
ax.set_title("Religious composition 1650-2025  [deepest inheritance, α=0.8 (central)]\n"
             "× marks: Pew 2023-24 Religious Landscape Study (adults)")
ax.legend(ncol=4, fontsize=9, loc="center left")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "fig_religious_shares.png"), dpi=150)
plt.close(fig)

# --- 2. deepest-cluster share by group, 2025, 2 rules ---------------------------
x = np.arange(len(GROUPS))
w = 0.30
fig, ax = plt.subplots(figsize=(10.5, 5.2))
for i, r in enumerate(RULES):
    last = dfs[r].iloc[-1]
    vals = [100 * last[f"Cn_share_{g}"] for g in GROUPS]
    ax.bar(x + (i - 1) * w, vals, w, label=RD[r], color=plt.cm.Set2(i))
ax.set_xticks(x)
ax.set_xticklabels([GD[g] for g in GROUPS], rotation=18, ha="right")
ax.set_ylabel("Share of group in deepest cluster C₁₄ in 2025 (%)")
ax.set_title("Deepest-cluster share by religious group, 2025  [α=0.8 (central)]")
ax.legend()
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "fig_religious_cn_by_group.png"), dpi=150)
plt.close(fig)

# --- 3. within-group cluster distribution, 2025 (deepest, 0.8 central) -------------------
fs = np.load(os.path.join(INDIR,
                          "finalstate_religious_deepest_inheritance_alpha0.8.npz"))
Y = fs["Y"][-1].reshape(8, 2, 15)          # [g, sex, d]
W = Y.sum(axis=1)                          # [g, d]
W = W / W.sum(axis=1, keepdims=True)
bins = [(0, 1, "C₀"), (1, 5, "C₁–₄"), (5, 10, "C₅–₉"),
        (10, 14, "C₁₀–₁₃"), (14, 15, "C₁₄")]
fig, ax = plt.subplots(figsize=(10.5, 5.2))
x = np.arange(len(GROUPS))
bot = np.zeros(len(GROUPS))
for lo, hi, lab in bins:
    v = 100 * W[:, lo:hi].sum(axis=1)
    ax.bar(x, v, bottom=bot, label=lab)
    bot += v
ax.set_xticks(x)
ax.set_xticklabels([GD[g] for g in GROUPS], rotation=18, ha="right")
ax.set_ylabel("Share of group (%)")
ax.set_title("2025 cluster distribution within each religious group  "
             "[deepest inheritance, α=0.8 (central)]")
ax.legend(ncol=5, fontsize=9)
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "fig_religious_cluster_dist.png"), dpi=150)
plt.close(fig)

# --- 4. total deepest-cluster share time series, 3 rules ------------------------
fig, ax = plt.subplots(figsize=(10, 5))
for i, r in enumerate(RULES):
    ax.plot(dfs[r]["year"], 100 * dfs[r]["Cn_share"], lw=1.8,
            color=plt.cm.Set2(i), label=RD[r])
ax.set_xlim(1650, 2025)
ax.set_xlabel("Year")
ax.set_ylabel("Deepest-cluster (C₁₄) share of total population (%)")
ax.set_title("Deepest-cluster share over time by bloodline rule  "
             "[religious model, α=0.8 (central)]")
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "fig_cn_total_series.png"), dpi=150)
plt.close(fig)

print("figures written to", FIGDIR)
