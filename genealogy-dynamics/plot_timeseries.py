"""Cluster shares over time (total population = paternal + maternal), with
historical period domains.

For each bloodline rule: x = year (1650-2025), y = share of population (%),
one colored curve per cluster (C0..Cn; total population, summed over the
paternal and maternal populations), vertical shaded bands for historical
periods. The colonial input directory is auto-discovered (newest
outputs/*_colonial) so a hardcoded stamp can never go stale.
"""
import os
import glob
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from genealogy.populations import N_CLUSTERS, MAX_DEPTH, RULE_DISPLAY

BASE = os.path.dirname(os.path.abspath(__file__))


def _newest_colonial_dir():
    cands = sorted(glob.glob(os.path.join(BASE, "outputs", "*_colonial")))
    if not cands:
        raise RuntimeError("no outputs/*_colonial directory found")
    return cands[-1]


COLONIAL_DIR = _newest_colonial_dir()
print("colonial dir:", COLONIAL_DIR)
STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
OUTDIR = os.path.join(BASE, "outputs", f"{STAMP}_timeseries")
os.makedirs(OUTDIR, exist_ok=True)

# (label, start, end, color): warm -> cool progression across US history.
PERIODS = [
    ("Colonial", 1650, 1776, "#ffcccc"),
    ("Early Republic", 1776, 1860, "#ffe0b3"),
    ("Gilded Age", 1861, 1913, "#fff8c9"),
    ("World Wars", 1914, 1945, "#d8ecd8"),
    ("Postwar", 1946, 2025, "#d9e8ff"),
]

CLASS_COLORS = plt.get_cmap("tab20").colors  # 20 distinct colors


def plot(rule, title, fname):
    df = pd.read_csv(os.path.join(COLONIAL_DIR, f"results_{rule}_alpha0.8.csv"))
    years = df["year"].to_numpy()
    fig, ax = plt.subplots(figsize=(13, 6.8))

    for label, s, e, color in PERIODS:
        ax.axvspan(s, e, color=color, alpha=0.55, zorder=0, lw=0)
        ax.text((s + e) / 2, 96.5, label, ha="center", va="top",
                fontsize=8, style="italic", color="#555555", zorder=5)

    for d in range(N_CLUSTERS):
        # total-population share of cluster C_d (paternal + maternal)
        y = (df[f"share_C{d}_pat"].to_numpy()
             + df[f"share_C{d}_mat"].to_numpy()) * 100.0
        ax.plot(years, y, color=CLASS_COLORS[d % 20],
                lw=2.4 if d == 0 else 1.4,
                label="C0 foreign-born" if d == 0
                else (f"C{d} ({MAX_DEPTH}+)" if d == MAX_DEPTH else f"C{d}"),
                zorder=3)

    ax.set_xlim(1650, 2025)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Year")
    ax.set_ylabel("Share of population (%)")
    ax.set_title(title + "\n25-year snapshots 1650–2025 · eps=0.00075/yr, α=0.8 (enforced central endogamy)",
                 fontsize=11)
    ax.grid(True, alpha=0.2)
    ax.legend(title="Cluster (total population)", fontsize=8, title_fontsize=9,
              loc="upper left", bbox_to_anchor=(1.01, 1.0),
              frameon=True)
    fig.tight_layout()
    out = os.path.join(OUTDIR, fname)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


specs = [
    ("shallowest_inheritance",
     f"Cluster shares over time [shallowest inheritance: child = min(1+min(parents), {MAX_DEPTH})]",
     "ts_shallowest_inheritance.png"),
    ("deepest_inheritance",
     f"Cluster shares over time [deepest inheritance: child = min(1+max(parents), {MAX_DEPTH})]",
     "ts_deepest_inheritance.png"),
]
for rule, title, fname in specs:
    print("wrote", plot(rule, title, fname))
print("outdir:", OUTDIR)
