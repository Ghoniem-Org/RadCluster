#!/usr/bin/env python3
"""Step-6 tail study: 2025 cluster distribution at n = 12/14/18/24 (two-sex
colonial, shallowest/deepest x alpha 0.0/0.8/0.9) and the n=14 family comparison.

Reads the step-6 output dirs (paths passed as args or discovered from
outputs/step6_logs/*.log "outdir:" lines), then:

1. Cap-share vs n table (2025 share of C_n, both rules, all alphas).
2. Sub-cap tail shape: for each n, the last ~6 clusters' shares at alpha=0.8
   (shallowest + deepest) -- does the distribution decay before the cap or pile
   up at it? Quantify: cap share vs sum of the five preceding clusters.
3. Totals-vs-n check (totals must be n-invariant).
4. Writes outputs/step6_figures/step6_tail_shape.png
   (2025 distribution overlays for n=12/14/18/24, shallowest and deepest,
   alpha=0.8) and outputs/step6_figures/step6_n_comparison.png
   (cap share vs n + total-population invariance panel).

Usage: python3 analyze_tail.py
"""
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
LOGDIR = os.path.join(BASE, "outputs", "step6_logs_rev")
FIGDIR = os.path.join(BASE, "outputs", "step6_figures")
os.makedirs(FIGDIR, exist_ok=True)


def discover_dirs():
    """Map n -> outdir from the step-6 run logs' 'outdir:' lines."""
    dirs = {}
    for log, n in (("n12_twosex.log", 12), ("n14_twosex.log", 14),
                   ("n18_twosex.log", 18), ("n24_twosex.log", 24)):
        path = os.path.join(LOGDIR, log)
        if not os.path.exists(path):
            continue
        with open(path) as f:
            for line in f:
                m = re.search(r"outdir:\s*(\S+)", line)
                if m:
                    dirs[n] = m.group(1)
    return dirs


def final_dist(d, n, rule, alpha):
    """2025 cluster-total shares for (rule, alpha) from results CSV."""
    df = pd.read_csv(os.path.join(d, f"results_{rule}_alpha{alpha}.csv"))
    row = df[df["year"] == 2025.0].iloc[0]
    tot = float(row["total_M"])
    return np.array([float(row[f"share_C{dd}"]) for dd in range(n + 1)]), tot


def main():
    dirs = discover_dirs()
    print("dirs:", dirs)
    rules = ["shallowest_inheritance", "deepest_inheritance"]
    alphas = [0.0, 0.8, 0.9]

    # --- 1/2/3: cap-share table + totals --------------------------------------
    print("\n## 2025 cap share (C_n) by n")
    hdr = "rule        | " + " | ".join(f"alpha={a}" for a in alphas)
    print(hdr)
    for rule in rules:
        line = f"{rule:12s}| "
        for a in alphas:
            cells = []
            for n, d in sorted(dirs.items()):
                share, tot = final_dist(d, n, rule, a)
                cells.append(f"n={n}: {share[-1]*100:.3f}%")
            line += " | ".join(cells) + "  || "
        print(line)
    print("\n## 2025 totals (M) by n -- must be n-invariant")
    for rule in rules:
        for a in [0.8]:
            row = [f"n={n}: {final_dist(d, n, rule, a)[1]:.2f}"
                   for n, d in sorted(dirs.items())]
            print(f"{rule} alpha={a}: " + "  ".join(row))
    print("\n## sub-cap tail (alpha=0.8): cap vs sum of preceding 5 clusters")
    for rule in rules:
        for n, d in sorted(dirs.items()):
            share, _ = final_dist(d, n, rule, 0.8)
            cap = share[-1] * 100
            pre5 = share[-6:-1].sum() * 100
            tail = share[-10:-1][::-1] if n >= 10 else share[::-1]
            decays = all(tail[i] >= tail[i + 1] for i in range(len(tail) - 1))
            print(f"{rule} n={n:2d}: cap={cap:6.3f}%  "
                  f"C_(n-5)..C_(n-1)={pre5:6.3f}%  "
                  f"cap/pre5={cap/max(pre5, 1e-12):.2f}  "
                  f"monotone-decay(last10)={decays}")

    # --- figure 1: distribution overlays --------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    colors = {12: "#d95f02", 14: "#1f78b4", 18: "#33a02c", 24: "#6a3d9a"}
    for ax, rule in zip(axes, rules):
        for n in sorted(dirs):
            d = dirs[n]
            share, _ = final_dist(d, n, rule, 0.8)
            xs = np.arange(n + 1)
            ax.semilogy(xs, np.maximum(share, 1e-9), "o-",
                        color=colors[n], ms=3, lw=1.4,
                        label=f"n={n} (cap {share[-1]*100:.2f}%)")
            ax.plot([n], [max(share[-1], 1e-9)], "o", color=colors[n],
                    ms=8, mfc="none", mew=2)
        ax.set_xlabel("cluster depth d")
        ax.set_title(f"{'shallowest' if 'shallow' in rule else 'deepest'} inheritance, "
                     f"2025, alpha=0.8")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("2025 population share (log)")
    fig.suptitle("Cluster-distribution tail, 1650-2025 colonial run "
                 "(two-sex model); ringed points = capped cluster C_n")
    fig.tight_layout()
    p1 = os.path.join(FIGDIR, "step6_tail_shape.png")
    fig.savefig(p1, dpi=150)
    print("\nwrote", p1)

    # --- figure 2: cap share vs n + total invariance --------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ns = sorted(dirs)
    for rule, marker, lab in (("shallowest_inheritance", "o", "shallowest"),
                              ("deepest_inheritance", "s", "deepest")):
        caps = [final_dist(dirs[n], n, rule, 0.8)[0][-1] * 100 for n in ns]
        axes[0].plot(ns, caps, marker + "-", lw=1.8,
                     label=f"{lab} (alpha=0.8)")
    axes[0].set_xlabel("max depth n")
    axes[0].set_ylabel("2025 cap share C_n (%)")
    axes[0].set_title("Cap share vs truncation depth (alpha=0.8)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    for rule, marker, lab in (("shallowest_inheritance", "o", "shallowest"),
                              ("deepest_inheritance", "s", "deepest")):
        tots = np.array([final_dist(dirs[n], n, rule, 0.8)[1] for n in ns])
        ppm = (tots - tots[0]) / tots[0] * 1e6
        axes[1].plot(ns, ppm, marker + "-", lw=1.8, label=lab)
    axes[1].axhline(0, color="k", lw=0.8, ls="--")
    axes[1].set_xlabel("max depth n")
    axes[1].set_ylabel("2025 total vs n=12 (ppm)")
    axes[1].set_title("Total population vs n (solver noise only)")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    fig.suptitle("Truncation-depth study: deeper n redistributes the tail, "
                 "not the total")
    fig.tight_layout()
    p2 = os.path.join(FIGDIR, "step6_n_comparison.png")
    fig.savefig(p2, dpi=150)
    print("wrote", p2)


if __name__ == "__main__":
    main()
