"""Distribution functions: total population vs cluster at 25-year snapshots.

For each bloodline rule, draws one figure: a smooth continuous curve per
snapshot year (1650..2025), colored by year with a colorbar. PCHIP
interpolation is used so curves stay smooth without overshooting the
discrete cluster values. Cluster values are total-population shares
(paternal + maternal).
"""
import os
import glob
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.interpolate import PchipInterpolator

from genealogy.populations import N_CLUSTERS, MAX_DEPTH

BASE = os.path.dirname(os.path.abspath(__file__))


def _newest_colonial_dir():
    cands = sorted(glob.glob(os.path.join(BASE, "outputs", "*_colonial")))
    if not cands:
        raise RuntimeError("no outputs/*_colonial directory found")
    return cands[-1]


COLONIAL_DIR = _newest_colonial_dir()
print("colonial dir:", COLONIAL_DIR)
STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
OUTDIR = os.path.join(BASE, "outputs", f"{STAMP}_distributions")
os.makedirs(OUTDIR, exist_ok=True)

CLASSES = np.arange(N_CLUSTERS)
FINE = np.linspace(0, MAX_DEPTH, 20 * MAX_DEPTH + 1)


def load(rule):
    path = os.path.join(COLONIAL_DIR, f"snapshots_{rule}.csv")
    df = pd.read_csv(path)
    # total-population share per cluster (paternal + maternal)
    pct = np.array([(df[f"share_C{d}_pat"] + df[f"share_C{d}_mat"]).to_numpy()
                    for d in CLASSES]).T * 100.0
    return df["year"].values, pct


def smooth_curve(y):
    f = PchipInterpolator(CLASSES, np.clip(y, 0, None), extrapolate=False)
    return np.clip(f(FINE), 0, None)


def plot(rule, title, fname):
    years, pct = load(rule)
    fig, ax = plt.subplots(figsize=(11, 6.5))
    cmap = plt.get_cmap("viridis")
    norm = Normalize(vmin=years.min(), vmax=years.max())
    for i, yr in enumerate(years):
        ax.plot(FINE, smooth_curve(pct[i]), color=cmap(norm(yr)), lw=1.8)
    # mark the discrete cluster values lightly for the final year
    ax.scatter(CLASSES, pct[-1], color=cmap(norm(years[-1])), s=18, zorder=5)
    ax.set_xlim(0, MAX_DEPTH)
    ax.set_ylim(0, 100)
    ax.set_xticks(CLASSES)
    ax.set_xticklabels(["C0\nforeign"] + [f"C{d}" for d in range(1, N_CLUSTERS)])
    ax.set_xlabel("Cluster (total population)")
    ax.set_ylabel("Share of population (%)")
    ax.set_title(title + f"\n25-year snapshots 1650–2025 · "
                 f"eps=0.00075/yr, alpha=0.6", fontsize=11)
    ax.grid(True, alpha=0.25)
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax)
    cbar.set_label("Snapshot year")
    fig.tight_layout()
    out = os.path.join(OUTDIR, fname)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


specs = [
    ("shallowest_inheritance",
     f"Population share vs cluster [shallowest inheritance: child = min(1+min(parents), {MAX_DEPTH})]",
     "dist_shallowest_inheritance_linear.png"),
    ("deepest_inheritance",
     f"Population share vs cluster [deepest inheritance: child = min(1+max(parents), {MAX_DEPTH})]",
     "dist_deepest_inheritance_linear.png"),
]
for rule, title, fname in specs:
    print("wrote", plot(rule, title, fname))
print("outdir:", OUTDIR)
