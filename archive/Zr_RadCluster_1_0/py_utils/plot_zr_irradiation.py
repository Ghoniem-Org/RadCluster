"""
plot_zr_irradiation.py
Reads Irradiation_Microstructure_Database_Zr_Zircaloys.xlsx and produces
publication-quality figures in a variety of formats.

Usage:
    python plot_zr_irradiation.py
    python plot_zr_irradiation.py --xlsx /path/to/file.xlsx --outdir figures/

Output figures (PDF + PNG at 300 dpi):
    fig01  <a>-loop density vs temperature (by material)
    fig02  <a>-loop mean size vs temperature (by material)
    fig03  <a>-loop density vs dose (by alloy)
    fig04  <a>-loop mean size vs dose (by alloy)
    fig05  <c>-loop density vs dose (Zircaloy-4 proton)
    fig06  <c>-loop mean size vs dose (Zircaloy-4 proton)
    fig07  <a>-loop size distributions (T comparison)
    fig08  <c>-loop size distribution (Tournadre data)
    fig09  Combined: <a> and <c> loop density on same axes
    fig10  Nb effect: Zy-2 vs ZIRLO density and size
    fig11  Temperature-series overview (loops + SPP state)
    fig12  Dose-series overview (all alloys, fixed T)
"""

import argparse, os, warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

warnings.filterwarnings("ignore")
matplotlib.rcParams.update({
    "font.family":       "serif",
    "font.serif":        ["Times New Roman", "DejaVu Serif"],
    "font.size":         11,
    "axes.labelsize":    12,
    "axes.titlesize":    11,
    "legend.fontsize":   9,
    "xtick.labelsize":   10,
    "ytick.labelsize":   10,
    "axes.linewidth":    1.2,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
})

# ── palette ───────────────────────────────────────────────────────────────────
C = {
    "Zy2":     "#1B4F8A",
    "Zy4":     "#C0392B",
    "ZIRLO":   "#117A65",
    "ZrNb":    "#7D6608",
    "PureZr":  "#6C3483",
    "IonZr":   "#D35400",
    "BWR":     "#1A3A5C",
    "Neutron": "#2C5F8A",
    "Proton":  "#E67E22",
    "Ion":     "#922B21",
}
MARKERS = {"Zy2":"o","Zy4":"s","ZIRLO":"^","ZrNb":"D","PureZr":"v","IonZr":"P"}

def save(fig, outdir, stem):
    for ext in ("pdf","png"):
        fig.savefig(os.path.join(outdir, f"{stem}.{ext}"))
    plt.close(fig)

def bar_kw(color, alpha=0.82):
    return dict(color=color, alpha=alpha, edgecolor="white", linewidth=0.5)

def annotate(ax, text, x=0.97, y=0.96, fontsize=8):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=fontsize,
            ha="right", va="top", color="#555",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1))

# ─────────────────────────────────────────────────────────────────────────────
def load(xlsx, sheet):
    df = pd.read_excel(xlsx, sheet_name=sheet, header=2)
    df.columns = [str(c).strip() for c in df.columns]
    return df.dropna(how="all")

def num(s):
    try:
        return float(str(s).replace("~","").replace("<","").replace(">","")
                    .replace("≤","").replace("≥","").split("–")[0].split("×")[0])
    except:
        return np.nan

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlsx", default="Irradiation_Microstructure_Database_Zr_Zircaloys.xlsx")
    parser.add_argument("--outdir", default="figures_zr")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    xl = args.xlsx
    df_a  = load(xl, "a-Loop Density & Mean Size")
    df_c  = load(xl, "<c>-Loop Density & Mean Size")
    df_b  = load(xl, "Loop Size Dist. Bins")
    df_v  = load(xl, "Void & Cavity Data")
    df_s  = load(xl, "Second Phase Particles")
    df_ts = load(xl, "T-Series (Fixed Dose)")
    df_ds = load(xl, "Dose-Series (Fixed T)")

    # ── helper: extract columns by partial name ───────────────────────────────
    def col(df, partial):
        m = [c for c in df.columns if partial.lower() in c.lower()]
        return m[0] if m else None

    # column refs for a-loop sheet
    a_mat  = col(df_a, "Material")
    a_irr  = col(df_a, "Irrad.")
    a_dose = col(df_a, "Dose")
    a_temp = col(df_a, "Temp.")
    a_dens = col(df_a, "a>-Loop Density") or col(df_a, "Density")
    a_diam = col(df_a, "Mean") or col(df_a, "Diam")
    a_cite = col(df_a, "Cite")

    # column refs for c-loop sheet
    c_mat  = col(df_c, "Material")
    c_dose = col(df_c, "Dose")
    c_temp = col(df_c, "Temp.")
    c_dens = col(df_c, "c>-Loop Density") or col(df_c, "Density")
    c_diam = col(df_c, "Mean") or col(df_c, "Diam")

    # column refs for bin sheet
    b_mat  = col(df_b, "Material")
    b_bw   = col(df_b, "Bin Width")
    b_bc   = col(df_b, "Bin Centre")
    b_freq = col(df_b, "Rel. Freq")
    b_dens = col(df_b, "Abs. Density")
    b_type = col(df_b, "Loop Type")
    b_temp = col(df_b, "Temp")
    b_dose = col(df_b, "Dose")

    # column refs for T-series
    ts_mat  = col(df_ts, "Material")
    ts_irr  = col(df_ts, "Irrad.")
    ts_dose = col(df_ts, "Dose")
    ts_temp = col(df_ts, "Temp.")
    ts_ad   = [c for c in df_ts.columns if "a>-Loop Density" in c or "a>-loop Density" in c]
    ts_as   = [c for c in df_ts.columns if "a>-Loop Mean" in c or "a>-loop Mean" in c]
    ts_cd   = [c for c in df_ts.columns if "c>-Loop Density" in c or "c>-loop Density" in c]
    ts_cs   = [c for c in df_ts.columns if "c>-Loop Mean" in c or "c>-loop Mean" in c]
    ts_sp   = col(df_ts, "SPP")
    ts_grow = col(df_ts, "Growth")

    # column refs for dose series
    ds_mat  = col(df_ds, "Material")
    ds_irr  = col(df_ds, "Irrad.")
    ds_temp = col(df_ds, "Temp")
    ds_dose = col(df_ds, "Dose")

    def asub(mat_contains=None, irr_contains=None):
        m = df_a.copy()
        if mat_contains:
            m = m[m[a_mat].astype(str).str.contains(mat_contains, case=False, na=False)]
        if irr_contains:
            m = m[m[a_irr].astype(str).str.contains(irr_contains, case=False, na=False)]
        m["_T"] = m[a_temp].apply(num)
        m["_D"] = m[a_dose].apply(num)
        m["_N"] = m[a_dens].apply(num)
        m["_d"] = m[a_diam].apply(num)
        return m

    def csub(mat_contains=None):
        m = df_c.copy()
        if mat_contains:
            m = m[m[c_mat].astype(str).str.contains(mat_contains, case=False, na=False)]
        m["_T"] = m[c_temp].apply(num)
        m["_D"] = m[c_dose].apply(num)
        m["_N"] = m[c_dens].apply(num)
        m["_d"] = m[c_diam].apply(num)
        return m

    # ═════════════════════════════════════════════════════════════════════════
    # FIG 01 – <a>-loop density vs temperature
    # ═════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(8, 5))
    specs = [
        ("Zircaloy-2","Proton",C["Zy2"],"o","Zircaloy-2, 2 dpa, proton [Harte 2018]"),
        ("ZIRLO","Proton",     C["ZIRLO"],"^","Low-Sn ZIRLO, 2 dpa, proton [Harte 2018]"),
        ("Crystal bar","Neutron",C["PureZr"],"v","Crystal bar Zr, 1 dpa, neutron [Northwood 1979]"),
    ]
    for mat, irr, col_, mk, lbl in specs:
        d = asub(mat, irr)
        d = d.dropna(subset=["_T","_N"])
        if len(d) < 2: continue
        d = d.sort_values("_T")
        ax.plot(d["_T"], d["_N"], marker=mk, color=col_, linewidth=2, markersize=7, label=lbl)

    ax.set_xlabel("Irradiation temperature (°C)")
    ax.set_ylabel(r"$\langle a\rangle$-loop density ($\times10^{22}$ m$^{-3}$)")
    ax.set_title(r"$\langle a\rangle$-Loop Number Density vs Irradiation Temperature", fontweight="bold")
    ax.axvline(300, color="gray", ls="--", lw=1, alpha=0.6, label="I→V transition ~300°C")
    ax.legend(frameon=False, fontsize=9)
    ax.set_xlim(250, 470)
    ax.set_ylim(0, 8)
    ax.grid(True, ls=":", alpha=0.3)
    annotate(ax, "Nb-free: density ↓ with T\nNb-bearing: nearly constant")
    save(fig, args.outdir, "fig01_aloop_density_vs_T")

    # ═════════════════════════════════════════════════════════════════════════
    # FIG 02 – <a>-loop mean size vs temperature
    # ═════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(8, 5))
    for mat, irr, col_, mk, lbl in specs:
        d = asub(mat, irr)
        d = d.dropna(subset=["_T","_d"])
        if len(d) < 2: continue
        d = d.sort_values("_T")
        ax.plot(d["_T"], d["_d"], marker=mk, color=col_, linewidth=2, markersize=7, label=lbl)

    ax.set_xlabel("Irradiation temperature (°C)")
    ax.set_ylabel(r"Mean $\langle a\rangle$-loop diameter $\bar{d}$ (nm)")
    ax.set_title(r"$\langle a\rangle$-Loop Mean Diameter vs Irradiation Temperature", fontweight="bold")
    ax.annotate("", xy=(450, 30), xytext=(280, 4),
                arrowprops=dict(arrowstyle="->", color=C["Zy2"], lw=1.5))
    ax.text(330, 18, "×7.5 increase\n(Zircaloy-2)", color=C["Zy2"], fontsize=8)
    ax.legend(frameon=False, fontsize=9)
    ax.set_xlim(250, 470); ax.set_ylim(0, 38)
    ax.grid(True, ls=":", alpha=0.3)
    save(fig, args.outdir, "fig02_aloop_size_vs_T")

    # ═════════════════════════════════════════════════════════════════════════
    # FIG 03 – <a>-loop density vs dose
    # ═════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(8, 5))
    dose_specs = [
        ("Zircaloy-2","Proton", C["Zy2"],"o","Zircaloy-2, 350°C, proton [Balogh 2021]"),
        ("Crystal bar","Neutron",C["PureZr"],"v","Crystal bar Zr, 300°C, neutron [Northwood 1979]"),
        ("α-Zr","Ion", C["IonZr"],"P","α-Zr, 300°C, in-situ ion [Idrees 2013]"),
    ]
    for mat, irr, col_, mk, lbl in dose_specs:
        d = asub(mat, irr)
        d = d.dropna(subset=["_D","_N"])
        if len(d) < 1: continue
        d = d.sort_values("_D")
        ax.semilogx(d["_D"], d["_N"], marker=mk, color=col_, linewidth=2, markersize=7, label=lbl)

    ax.set_xlabel("Dose (dpa)")
    ax.set_ylabel(r"$\langle a\rangle$-loop density ($\times10^{22}$ m$^{-3}$)")
    ax.set_title(r"$\langle a\rangle$-Loop Density vs Dose (Fixed Temperature)", fontweight="bold")
    ax.text(3, 7.2, "SXRD total DD is 4–15× higher\n(sub-TEM loops undetected)", fontsize=8,
            color="gray", style="italic")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(True, ls=":", alpha=0.3, which="both")
    save(fig, args.outdir, "fig03_aloop_density_vs_dose")

    # ═════════════════════════════════════════════════════════════════════════
    # FIG 04 – <a>-loop mean size vs dose
    # ═════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(8, 5))
    for mat, irr, col_, mk, lbl in dose_specs:
        d = asub(mat, irr)
        d = d.dropna(subset=["_D","_d"])
        if len(d) < 1: continue
        d = d.sort_values("_D")
        ax.semilogx(d["_D"], d["_d"], marker=mk, color=col_, linewidth=2, markersize=7, label=lbl)

    # fit power law for Zy-2
    phi = np.array([2.3, 4.7, 7.0]); d_fit = np.array([8.0, 9.5, 11.0])
    if len(phi) > 1:
        m = np.polyfit(np.log(phi), np.log(d_fit - 4.5), 1)
        phi_p = np.logspace(np.log10(0.8), np.log10(12), 50)
        d_p = np.exp(m[1]) * phi_p**m[0] + 4.5
        ax.semilogx(phi_p, d_p, "--", color=C["Zy2"], lw=1.2, alpha=0.7,
                    label=f"Power-law fit: m≈{m[0]:.2f}")

    ax.set_xlabel("Dose (dpa)")
    ax.set_ylabel(r"Mean $\langle a\rangle$-loop diameter $\bar{d}$ (nm)")
    ax.set_title(r"$\langle a\rangle$-Loop Mean Size vs Dose (Fixed Temperature)", fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(True, ls=":", alpha=0.3, which="both")
    save(fig, args.outdir, "fig04_aloop_size_vs_dose")

    # ═════════════════════════════════════════════════════════════════════════
    # FIG 05–06 – <c>-loop density and size vs dose
    # ═════════════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    dc = csub("Zircaloy-4").dropna(subset=["_D"]).sort_values("_D")
    dn = csub("Zircaloy-2").dropna(subset=["_D"]).sort_values("_D")

    ax = axes[0]
    dc_n = dc.dropna(subset=["_N"])
    dn_n = dn.dropna(subset=["_N"])
    if len(dc_n): ax.plot(dc_n["_D"], dc_n["_N"], "s-", color=C["Zy4"], lw=2, ms=7,
                          label="Zircaloy-4, 350°C, proton [Tournadre 2012]")
    if len(dn_n): ax.plot(dn_n["_D"], dn_n["_N"], "o--", color=C["Zy2"], lw=2, ms=7,
                          label="Zircaloy-2, ~280°C, neutron [Cockeram 2011]")
    ax.axvline(2, color="gray", ls=":", lw=1.2, alpha=0.7, label="onset ~2 dpa")
    ax.set_xlabel("Dose (dpa)"); ax.set_ylabel(r"$\langle c\rangle$-loop density ($\times10^{21}$ m$^{-3}$)")
    ax.set_title(r"$\langle c\rangle$-Loop Density vs Dose", fontweight="bold")
    ax.legend(frameon=False, fontsize=8); ax.grid(True, ls=":", alpha=0.3)

    ax = axes[1]
    dc_d = dc.dropna(subset=["_d"])
    dn_d = dn.dropna(subset=["_d"])
    if len(dc_d): ax.plot(dc_d["_D"], dc_d["_d"], "s-", color=C["Zy4"], lw=2, ms=7,
                          label="Zircaloy-4, proton [Tournadre 2012]")
    if len(dn_d): ax.plot(dn_d["_D"], dn_d["_d"], "o--", color=C["Zy2"], lw=2, ms=7,
                          label="Zircaloy-2, neutron [Cockeram 2011]")
    ax.set_xlabel("Dose (dpa)"); ax.set_ylabel(r"Mean $\langle c\rangle$-loop diameter $\bar{d}_c$ (nm)")
    ax.set_title(r"$\langle c\rangle$-Loop Mean Size vs Dose", fontweight="bold")
    ax.legend(frameon=False, fontsize=8); ax.grid(True, ls=":", alpha=0.3)

    fig.suptitle(r"$\langle c\rangle$-Component Loop Evolution — Zircaloy-4 & Zircaloy-2",
                 fontweight="bold")
    fig.tight_layout()
    save(fig, args.outdir, "fig05_06_cloop_vs_dose")

    # ═════════════════════════════════════════════════════════════════════════
    # FIG 07 – <a>-loop size distributions (temperature comparison)
    # ═════════════════════════════════════════════════════════════════════════
    zy2_temps = [(280, C["Zy2"], "Zircaloy-2, 280°C"),
                 (350, C["ZIRLO"], "Zircaloy-2, 350°C"),
                 (450, C["Zy4"], "Zircaloy-2, 450°C")]

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), sharey=False)
    fig.suptitle(r"$\langle a\rangle$-Loop Size Distributions: Zircaloy-2 at Fixed Dose (2 dpa), "
                 r"Varying Temperature [Harte 2018]", fontweight="bold")

    for ax, (T, col_, lbl) in zip(axes, zy2_temps):
        sub = df_b[
            df_b[b_mat].astype(str).str.contains("Zircaloy-2", na=False) &
            (df_b[b_temp].apply(num) == T) &
            df_b[b_type].astype(str).str.contains("<a>", na=False)
        ].dropna(subset=[b_bc, b_freq])
        if len(sub) == 0:
            ax.set_title(lbl); continue
        x = sub[b_bc].apply(num).values
        w = sub[b_bw].apply(num).values
        y = sub[b_freq].apply(num).values
        ok = ~(np.isnan(x)|np.isnan(y))
        ax.bar(x[ok], y[ok], width=w[ok]*0.85, color=col_, alpha=0.82,
               edgecolor="white", linewidth=0.5)
        mean_est = np.average(x[ok], weights=y[ok])
        ax.axvline(mean_est, color=col_, lw=1.8, ls="--", alpha=0.9,
                   label=f"Mean≈{mean_est:.0f} nm")
        ax.set_xlabel("Loop diameter (nm)")
        ax.set_ylabel("Relative frequency (%)")
        ax.set_title(lbl)
        ax.legend(frameon=False, fontsize=8)
        annotate(ax, "[CIT-05] Harte 2018")

    fig.tight_layout()
    save(fig, args.outdir, "fig07_aloop_size_distributions_T")

    # ═════════════════════════════════════════════════════════════════════════
    # FIG 08 – <c>-loop size distribution (Tournadre 2012)
    # ═════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(7, 4.5))
    sub_c = df_b[
        df_b[b_mat].astype(str).str.contains("Zircaloy-4", na=False) &
        df_b[b_type].astype(str).str.contains("<c>", na=False)
    ].dropna(subset=[b_bc, b_freq])

    if len(sub_c):
        x = sub_c[b_bc].apply(num).values
        w = sub_c[b_bw].apply(num).values
        y = sub_c[b_freq].apply(num).values
        ok = ~(np.isnan(x)|np.isnan(y))
        ax.bar(x[ok], y[ok], width=w[ok]*0.85, color=C["Zy4"], alpha=0.82,
               edgecolor="white", linewidth=0.5, label="Zircaloy-4, 11.5 dpa, 350°C")
        mean_c = np.average(x[ok], weights=y[ok])
        ax.axvline(mean_c, color=C["Zy4"], lw=2, ls="--",
                   label=f"Mean≈{mean_c:.0f} nm")

    ax.set_xlabel(r"$\langle c\rangle$-loop diameter (nm)")
    ax.set_ylabel("Relative frequency (%)")
    ax.set_title(r"$\langle c\rangle$-Loop Size Distribution — Zircaloy-4 "
                 r"(Proton, 11.5 dpa, 350°C) [Tournadre 2012]", fontweight="bold")
    ax.legend(frameon=False)
    ax.grid(True, ls=":", alpha=0.3)
    annotate(ax, "[CIT-06] Tournadre 2012\nVacancy type; basal plane habit")
    save(fig, args.outdir, "fig08_cloop_size_distribution")

    # ═════════════════════════════════════════════════════════════════════════
    # FIG 09 – Combined <a> and <c> loop density on same axes
    # ═════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(8, 5))

    # a-loops (Zy-2, proton T-series)
    d_azy2 = asub("Zircaloy-2","Proton").dropna(subset=["_T","_N"]).sort_values("_T")
    ax.plot(d_azy2["_T"], d_azy2["_N"], "o-", color=C["Zy2"], lw=2, ms=7,
            label=r"$\langle a\rangle$ Zircaloy-2, proton 2 dpa [Harte 2018]")

    d_azirlo = asub("ZIRLO","Proton").dropna(subset=["_T","_N"]).sort_values("_T")
    ax.plot(d_azirlo["_T"], d_azirlo["_N"], "^-", color=C["ZIRLO"], lw=2, ms=7,
            label=r"$\langle a\rangle$ ZIRLO, proton 2 dpa [Harte 2018]")

    # c-loops (right y-axis)
    ax2 = ax.twinx()
    dc_ts = csub("Zircaloy-4").dropna(subset=["_T","_N"]).sort_values("_T")
    if len(dc_ts):
        ax2.plot(dc_ts["_T"], dc_ts["_N"], "s--", color=C["Zy4"], lw=2, ms=7,
                 label=r"$\langle c\rangle$ Zircaloy-4, proton [Tournadre 2012]")
    ax2.set_ylabel(r"$\langle c\rangle$-loop density ($\times10^{21}$ m$^{-3}$)",
                   color=C["Zy4"])
    ax2.tick_params(axis="y", labelcolor=C["Zy4"])
    ax2.spines["right"].set_color(C["Zy4"])
    ax2.set_ylim(0, 6)

    ax.set_xlabel("Irradiation temperature (°C)")
    ax.set_ylabel(r"$\langle a\rangle$-loop density ($\times10^{22}$ m$^{-3}$)")
    ax.set_title(r"$\langle a\rangle$ and $\langle c\rangle$ Loop Density vs Temperature", fontweight="bold")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1+lines2, labels1+labels2, frameon=False, fontsize=9, loc="center right")
    ax.grid(True, ls=":", alpha=0.3)
    save(fig, args.outdir, "fig09_combined_aloop_cloop_density")

    # ═════════════════════════════════════════════════════════════════════════
    # FIG 10 – Nb effect: Zy-2 vs ZIRLO density and size (2x2 grid)
    # ═════════════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    fig.suptitle("Nb Effect on Dislocation Loop Microstructure: Zircaloy-2 vs Low-Sn ZIRLO\n"
                 "(2 dpa, proton irradiation, [Harte 2018])", fontweight="bold")

    temps = [280, 350, 450]
    labels = ["280°C", "350°C", "450°C"]

    for ax, alloy, col_ in zip([axes[0,0],axes[0,1]],
                                ["Zircaloy-2","ZIRLO"],
                                [C["Zy2"], C["ZIRLO"]]):
        d = asub(alloy, "Proton").dropna(subset=["_T","_N"]).sort_values("_T")
        ax.bar([str(t) for t in d["_T"].values], d["_N"].values,
               color=col_, alpha=0.82, edgecolor="white", linewidth=0.5)
        ax.set_xlabel("Temperature (°C)"); ax.set_ylabel(r"Density ($\times10^{22}$ m$^{-3}$)")
        ax.set_title(f"{alloy} — Loop Density")
        ax.grid(True, ls=":", alpha=0.3, axis="y")

    for ax, alloy, col_ in zip([axes[1,0],axes[1,1]],
                                ["Zircaloy-2","ZIRLO"],
                                [C["Zy2"], C["ZIRLO"]]):
        d = asub(alloy, "Proton").dropna(subset=["_T","_d"]).sort_values("_T")
        ax.bar([str(t) for t in d["_T"].values], d["_d"].values,
               color=col_, alpha=0.82, edgecolor="white", linewidth=0.5)
        ax.set_xlabel("Temperature (°C)"); ax.set_ylabel(r"Mean diameter $\bar{d}$ (nm)")
        ax.set_title(f"{alloy} — Loop Mean Size")
        ax.grid(True, ls=":", alpha=0.3, axis="y")

    fig.tight_layout()
    save(fig, args.outdir, "fig10_Nb_effect_Zy2_vs_ZIRLO")

    # ═════════════════════════════════════════════════════════════════════════
    # FIG 11 – Temperature-series overview (a-loop density + size side by side)
    # ═════════════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Temperature Series Overview — All Materials (Fixed Dose)\n"
                 "Left: Loop Density | Right: Loop Mean Size", fontweight="bold")

    for ax, ycol, ylabel in zip(axes, ["_N","_d"],
        [r"$\langle a\rangle$-loop density ($\times10^{22}$ m$^{-3}$)",
         r"Mean $\langle a\rangle$-loop diameter $\bar{d}$ (nm)"]):
        for mat, irr, col_, mk, lbl in specs:
            d = asub(mat, irr).dropna(subset=["_T", ycol]).sort_values("_T")
            if len(d) < 1: continue
            ax.plot(d["_T"], d[ycol], marker=mk, color=col_, lw=2, ms=7, label=lbl)
        ax.set_xlabel("Irradiation temperature (°C)")
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False, fontsize=8)
        ax.grid(True, ls=":", alpha=0.3)

    fig.tight_layout()
    save(fig, args.outdir, "fig11_T_series_overview")

    # ═════════════════════════════════════════════════════════════════════════
    # FIG 12 – Dose-series overview (density + size on log-dose axis)
    # ═════════════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Dose-Series Overview — All Materials (Fixed Temperature)\n"
                 "Left: Loop Density | Right: Loop Mean Size", fontweight="bold")

    for ax, ycol, ylabel in zip(axes, ["_N","_d"],
        [r"$\langle a\rangle$-loop density ($\times10^{22}$ m$^{-3}$)",
         r"Mean $\langle a\rangle$-loop diameter $\bar{d}$ (nm)"]):
        for mat, irr, col_, mk, lbl in dose_specs:
            d = asub(mat, irr).dropna(subset=["_D", ycol]).sort_values("_D")
            if len(d) < 1: continue
            ax.semilogx(d["_D"], d[ycol], marker=mk, color=col_, lw=2, ms=7, label=lbl)
        ax.set_xlabel("Dose (dpa)")
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False, fontsize=8)
        ax.grid(True, ls=":", alpha=0.3, which="both")

    fig.tight_layout()
    save(fig, args.outdir, "fig12_dose_series_overview")

    print(f"\nAll figures saved to: {os.path.abspath(args.outdir)}/")
    produced = sorted(os.listdir(args.outdir))
    print(f"  {len(produced)} files:")
    for f in produced:
        print(f"    {f}")

if __name__ == "__main__":
    main()
