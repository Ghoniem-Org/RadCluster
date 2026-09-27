"""Regenerate the colonial per-population heatmaps with historical-era shading.

Reads the 25-year snapshot CSVs (no re-solve needed) and re-plots both
rules x both populations (paternal/maternal) with vertical shaded era
bands. Overwrites heatmap_<rule>_<pop>.png in the newest *_colonial dir.
"""
import os
import glob

import numpy as np
import pandas as pd

from genealogy.visualization import plot_pedigree_heatmap
from genealogy.populations import N_CLUSTERS, MAX_DEPTH, RULES, RULE_DISPLAY

BASE = os.path.dirname(os.path.abspath(__file__))


def _newest_colonial_dir():
    cands = sorted(glob.glob(os.path.join(BASE, "outputs", "*_colonial")))
    if not cands:
        raise RuntimeError("no outputs/*_colonial directory found")
    return cands[-1]


COLONIAL_DIR = _newest_colonial_dir()
print("colonial dir:", COLONIAL_DIR)
VMAX = 0.40
CAP_NOTE = (f"eps=0.00075/yr, \u03b1=0.8 (central), s=0.512 \u00b7 shared scale 0\u2013{VMAX:.0%} "
            "(1650 C0=100% saturates)")
CHILD_FORMULA = {"shallowest_inheritance": "min",
                 "deepest_inheritance": "1+m"}

for rule in RULES:
    df = pd.read_csv(os.path.join(COLONIAL_DIR, f"snapshots_{rule}.csv"))
    snap_years = df["year"].to_numpy()
    for pop in ("pat", "mat"):
        total = df[[f"share_C{d}_pat" for d in range(N_CLUSTERS)]].to_numpy().sum(axis=1) \
              + df[[f"share_C{d}_mat" for d in range(N_CLUSTERS)]].to_numpy().sum(axis=1)
        Z = df[[f"share_C{d}_{pop}" for d in range(N_CLUSTERS)]].to_numpy().T
        out = os.path.join(COLONIAL_DIR, f"heatmap_{rule}_{pop}.png")
        subtitle = (f"child = min(1+{CHILD_FORMULA[rule]}"
                    f"(parents), {MAX_DEPTH}) | " + CAP_NOTE)
        plot_pedigree_heatmap(
            snap_years, Z, out,
            title=(f"Population cluster distribution, 1650-2025 "
                   f"[{RULE_DISPLAY[rule]} rule]"),
            subtitle=subtitle,
            vmin=0.0, vmax=VMAX, population=pop)
        print("wrote", out)
