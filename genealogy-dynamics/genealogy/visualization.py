"""Visualization: cluster distributions over time + final snapshot + sex split."""

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import numpy as np

from .populations import N_CLUSTERS, MAX_DEPTH, RULE_DISPLAY

CLASS_COLORS = plt.get_cmap("tab20").colors  # 20 distinct colors


def _cluster_labels_multiline():
    return ["C0\nforeign"] + [f"C{d}" for d in range(1, N_CLUSTERS)]


def plot_stacked_area(t, share, path, title, ylabel="Share of population"):
    """share: (len(t), N_CLUSTERS) fractions of the total population."""
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.stackplot(t, share.T, colors=[CLASS_COLORS[i % 20] for i in range(N_CLUSTERS)],
                 alpha=0.85)
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(0, 1)
    ax.set_xlabel("Year")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    leg_labels = [f"C0 foreign-born"] + [
        f"C{d}" + (f" ({MAX_DEPTH}+)" if d == MAX_DEPTH else "") for d in range(1, N_CLUSTERS)]
    ax.legend(leg_labels, title="Cluster", fontsize=8, title_fontsize=9,
              loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_final_distribution(T, share, path, title):
    """T: final cluster totals (millions); share: fractions."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    a1.bar(range(N_CLUSTERS), T,
           color=[CLASS_COLORS[i % 20] for i in range(N_CLUSTERS)])
    a1.set_xlabel("Cluster")
    a1.set_ylabel("Millions of people")
    a1.set_title("Population by cluster")
    a2.bar(range(N_CLUSTERS), share * 100.0,
           color=[CLASS_COLORS[i % 20] for i in range(N_CLUSTERS)])
    a2.set_xlabel("Cluster")
    a2.set_ylabel("Share of population (%)")
    a2.set_title("Share by cluster")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_birth_rate(t, births, path, title):
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(t, births, color="darkred", lw=2)
    ax.set_xlabel("Year")
    ax.set_ylabel("Births (M/yr)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_mf_split(t, Tpat, Tmat, path, title):
    """Paternal vs maternal total population over time."""
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(t, Tpat, color="#1f77b4", lw=2, label="paternal (P_pat)")
    ax.plot(t, Tmat, color="#d62728", lw=2, label="maternal (P_mat)")
    ax.set_xlabel("Year")
    ax.set_ylabel("Millions of people")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_pn_vs_alpha(alphas, series, path):
    """Deepest-cluster (C_n) share at t=200 yr vs assortativity, all rules.

    series: {rule_key: [shares]}; labels from RULE_DISPLAY.
    """
    from genealogy.populations import RULE_DISPLAY
    markers = {"shallowest_inheritance": "s-",
               "deepest_inheritance": "^-"}
    fig, ax = plt.subplots(figsize=(8, 5))
    for rule, vals in series.items():
        ax.plot(alphas, vals, markers.get(rule, "x-"), lw=2.5, ms=8,
                label=RULE_DISPLAY.get(rule, rule))
    ax.set_xlabel("Assortativity α")
    ax.set_ylabel("Deepest-cluster ($C_n$) share at t=200 yr (%)")
    ax.set_title("Deepest-cluster share vs assortativity")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_pn_time(t, series, path, title):
    """series: {label: C_n share (% over time)}."""
    fig, ax = plt.subplots(figsize=(10, 5))
    for label, y in series.items():
        ax.plot(t, y, lw=2.2, label=label)
    ax.set_xlabel("Year")
    ax.set_ylabel("Deepest-cluster ($C_n$) share (%)")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_scenario_legend(path):
    """Color key mapping tab20 colors to cluster indices."""
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.axis("off")
    for i in range(N_CLUSTERS):
        ax.add_patch(plt.Rectangle((0.05, 0.92 - 0.055 * i), 0.08, 0.04,
                                   color=CLASS_COLORS[i % 20]))
        ax.text(0.16, 0.94 - 0.055 * i,
                f"C{i}" + (" foreign-born" if i == 0
                           else (f" ({MAX_DEPTH}+)" if i == MAX_DEPTH else "")),
                fontsize=9, va="center")
    ax.set_title("Cluster color key")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_pedigree_heatmap(years, Z, path, title, subtitle="", vmin=0.0,
                          vmax=0.40, population=""):
    """Heatmap: cluster share (rows C_0..C_n) vs snapshot year (columns).

    years: 1-D snapshot years; Z: (N_CLUSTERS, n_snaps) shares (fractions).
    population: "" or "pat"/"mat" for per-population labeling.
    """
    fig, ax = plt.subplots(figsize=(13, 6.5))
    im = ax.imshow(Z, aspect="auto", origin="lower", vmin=vmin, vmax=vmax,
                   cmap="viridis",
                   extent=[years[0], years[-1], -0.5, N_CLUSTERS - 0.5])
    ax.set_yticks(range(N_CLUSTERS))
    ax.set_yticklabels(_cluster_labels_multiline(), fontsize=8)
    ax.set_xlabel("Year")
    ax.set_ylabel("Cluster")
    pop = f" [{population} population]" if population else ""
    ax.set_title(title + pop + ("\n" + subtitle if subtitle else ""))
    # historical era shading (matching the timeseries plots)
    for label, s, e, color in [
            ("Colonial", 1650, 1776, "#ffcccc"), ("Early Republic", 1776, 1860, "#ffe0b3"),
            ("Gilded Age", 1861, 1913, "#fff8c9"), ("World Wars", 1914, 1945, "#d8ecd8"),
            ("Postwar", 1946, 2025, "#d9e8ff")]:
        ax.axvspan(s, e, color=color, alpha=0.25, zorder=0, lw=0)
        ax.text((s + e) / 2, N_CLUSTERS - 0.7, label, ha="center", va="top",
                fontsize=7, style="italic", color="#444444", zorder=5)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Share of population")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
