"""
size_distributions.py
Reads Size_Distribution_Database_FM_Steels.xlsx and produces
publication-quality bar-graph figures of dislocation loop and
void/cavity size distributions — one figure per condition.

Programmatic use (Jupyter):
    from size_distributions import plot_size_distributions
    figs = plot_size_distributions(
        xlsx="../../input/Size_Distribution_Database_FM_Steels.xlsx",
        outdir="../../output/size_distributions",
        font_size=16, save=True, show=True,
    )

CLI use:
    python size_distributions.py --xlsx ... --outdir ...
"""

import argparse
import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore")

# ── Per-condition catalog ────────────────────────────────────────────────────
# Each entry produces ONE single-panel figure. Empty data → skipped.
# fields: kind, material, dose, temp, irr_type (optional filter), color,
#         mean (optional vline), citation, xlim, title_prefix, bimodal (void only)
LOOP_CASES = [
    dict(kind="loop", material="EUROFER97", dose=15,   temp=330, color="#1B4F8A",
         mean=3.4,  citation="[CIT-01] Klimenkov 2012", xlim=(0, 16),
         title="EUROFER97 — 15 dpa / 330 °C (BOR-60, neutron)"),
    dict(kind="loop", material="EUROFER97", dose=32,   temp=335, color="#E67E22",
         mean=4.8,  citation="[CIT-01] Klimenkov 2012", xlim=(0, 16),
         title="EUROFER97 — 32 dpa / 335 °C (BOR-60, neutron)"),
    dict(kind="loop", material="EUROFER97", dose=16.3, temp=250, color="#154360",
         mean=2.5,  citation="[CIT-03] Gaganidze 2011", xlim=(0, 14),
         title="EUROFER97 — 16.3 dpa / 250 °C (BOR-60, neutron)"),
    dict(kind="loop", material="EUROFER97", dose=16.3, temp=300, color="#1A5276",
         mean=4.0,  citation="[CIT-03] Gaganidze 2011", xlim=(0, 14),
         title="EUROFER97 — 16.3 dpa / 300 °C (BOR-60, neutron)"),
    dict(kind="loop", material="EUROFER97", dose=16.3, temp=350, color="#2980B9",
         mean=10.0, citation="[CIT-03] Gaganidze 2011", xlim=(0, 28),
         title="EUROFER97 — 16.3 dpa / 350 °C (BOR-60, neutron)"),
    dict(kind="loop", material="EUROFER97", dose=16.3, temp=415, color="#AED6F1",
         mean=12.5, citation="[CIT-03] Gaganidze 2011", xlim=(0, 28),
         title="EUROFER97 — 16.3 dpa / 415 °C (BOR-60, neutron)"),
    dict(kind="loop", material="T91", dose=15.4, temp=376, color="#117A65",
         mean=8.0,  citation="[CIT-05] Gao 2018", xlim=(0, 50),
         title="T91 — 15.4 dpa / 376 °C (BOR-60, neutron)"),
    dict(kind="loop", material="T91", dose=35.1, temp=376, color="#1E8449",
         mean=14.0, citation="[CIT-05] Gao 2018", xlim=(0, 50),
         title="T91 — 35.1 dpa / 376 °C (BOR-60, neutron)"),
    dict(kind="loop", material="HT9", dose=17.1, temp=377, color="#6E2F5E",
         mean=10.0, citation="[CIT-06] Gao 2019", xlim=(0, 50),
         title="HT9 — 17.1 dpa / 377 °C (BOR-60, neutron)"),
    dict(kind="loop", material="HT9", dose=35.1, temp=377, color="#922B21",
         mean=16.0, citation="[CIT-06] Gao 2019", xlim=(0, 50),
         title="HT9 — 35.1 dpa / 377 °C (BOR-60, neutron)"),
    dict(kind="loop", material="Pure Fe", irr_type="Ion",     color="#E67E22",
         mean=3.5,  citation="[CIT-24] Tanaka 2018", xlim=(0, 10),
         title="Pure Fe — Ion (6.4 MeV Fe³⁺, 30 dpa, 300 °C)"),
    dict(kind="loop", material="Pure Fe", irr_type="Neutron", color="#1B4F8A",
         mean=9.5,  citation="[CIT-08] Lambrecht 2009", xlim=(0, 26),
         title="Pure Fe — Neutron (0.19 dpa, 300 °C)"),
]

VOID_CASES = [
    dict(kind="void", material="EUROFER97", dose=15,   temp=330, color="#1B4F8A",
         citation="[CIT-01] Klimenkov 2012", xlim=(0, 5),
         title="EUROFER97 — 15 dpa / 330 °C (BOR-60, neutron)"),
    dict(kind="void", material="EUROFER97", dose=32,   temp=335, color="#2E86C1",
         citation="[CIT-01] Klimenkov 2012", xlim=(0, 5),
         title="EUROFER97 — 32 dpa / 335 °C (BOR-60, neutron)"),
    dict(kind="void", material="T91", dose=15.4, temp=376, color="#117A65",
         citation="[CIT-05] Gao 2018", xlim=(0, 5),
         title="T91 — 15.4 dpa / 376 °C (BOR-60, neutron)"),
    dict(kind="void", material="HT9", dose=17.1, temp=377, color="#1E8449",
         citation="[CIT-06] Gao 2019", xlim=(0, 5),
         title="HT9 — 17.1 dpa / 377 °C (BOR-60, neutron)"),
    dict(kind="void", material="T91", dose=35.1, temp=376, color="#922B21",
         citation="[CIT-05] Gao 2018", xlim=(0, 18), bimodal=True,
         title="T91 — 35.1 dpa / 376 °C (BOR-60, neutron) — bimodal"),
    dict(kind="void", material="HT9", dose=35.1, temp=377, color="#C0392B",
         citation="[CIT-06] Gao 2019", xlim=(0, 18), bimodal=True,
         title="HT9 — 35.1 dpa / 377 °C (BOR-60, neutron) — bimodal"),
]


def _slug(text):
    s = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return s.lower()


def _apply_style(font_family, font_serif, font_size, axes_linewidth, dpi, save_dpi):
    matplotlib.rcParams.update({
        "font.family":       font_family,
        "font.serif":        list(font_serif),
        "font.size":         font_size,
        "axes.labelsize":    font_size,
        "axes.titlesize":    font_size,
        "legend.fontsize":   font_size,
        "xtick.labelsize":   font_size,
        "ytick.labelsize":   font_size,
        "axes.linewidth":    axes_linewidth,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "figure.dpi":        dpi,
        "savefig.dpi":       save_dpi,
        "savefig.bbox":      "tight",
    })


def _load_bins(xlsx_path, sheet):
    df = pd.read_excel(xlsx_path, sheet_name=sheet, header=2)
    df.columns = [str(c).strip() for c in df.columns]
    return df.dropna(how="all")


def _annotate_cite(ax, text, fontsize):
    ax.text(0.97, 0.96, text, transform=ax.transAxes, fontsize=fontsize,
            ha="right", va="top", color="#555555",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=1))


def _save(fig, outdir, stem, formats):
    for ext in formats:
        fig.savefig(os.path.join(outdir, f"{stem}.{ext}"))


def _filter(df, mat_c, type_c, dose_c, temp_c, bin_c, freq_c,
            material=None, irr_type=None, dose=None, temp=None):
    m = df.copy()
    if material:
        m = m[m[mat_c].astype(str).str.contains(material, case=False)]
    if irr_type:
        m = m[m[type_c].astype(str).str.contains(irr_type, case=False)]
    if dose is not None:
        m = m[pd.to_numeric(m[dose_c], errors="coerce") == dose]
    if temp is not None:
        m = m[pd.to_numeric(m[temp_c], errors="coerce") == temp]
    return m.dropna(subset=[bin_c, freq_c])


def _xywh(df, bin_c, bw_c, freq_c):
    x = pd.to_numeric(df[bin_c], errors="coerce").values
    w = pd.to_numeric(df[bw_c],  errors="coerce").values
    y = pd.to_numeric(df[freq_c], errors="coerce").values
    ok = ~(np.isnan(x) | np.isnan(y))
    return x[ok], w[ok], y[ok]


def plot_size_distributions(
    xlsx="Ferritics_Irradiated_Microstructure_Data.xlsx",
    outdir="figures",
    save=True,
    show=False,
    formats=("pdf", "png"),
    include_overview=True,
    # graphics controls — single font_size drives all text
    font_size=16,
    font_family="serif",
    font_serif=("Times New Roman", "DejaVu Serif"),
    axes_linewidth=1.2,
    dpi=150,
    save_dpi=300,
    cite_fontsize=None,        # defaults to font_size - 4
    figsize=(8, 5),
):
    """
    Build one figure per condition (loop and void size distributions).

    Empty datasets are skipped. Returns dict{stem -> Figure} of figures
    actually produced.

    Set `include_overview=False` to skip the mean-loop-size-vs-T scatter.
    """
    if cite_fontsize is None:
        cite_fontsize = max(7, font_size - 4)

    _apply_style(font_family, font_serif, font_size, axes_linewidth, dpi, save_dpi)

    if save:
        os.makedirs(outdir, exist_ok=True)

    loops = _load_bins(xlsx, "Loop Bins")
    voids = _load_bins(xlsx, "Void Bins")

    def col(df, partial):
        matches = [c for c in df.columns if partial.lower() in c.lower()]
        return matches[0] if matches else None

    l_mat, l_type, l_dose, l_temp = (col(loops, k) for k in
                                     ("Material", "Irrad. Type", "Dose", "Temp"))
    l_bin, l_bw, l_freq = (col(loops, k) for k in
                           ("Bin Centre", "Bin Width", "Rel. Frequency"))

    v_mat, v_type, v_dose, v_temp = (col(voids, k) for k in
                                     ("Material", "Irrad. Type", "Dose", "Temp"))
    v_bin, v_bw, v_freq = (col(voids, k) for k in
                           ("Bin Centre", "Bin Width", "Rel. Frequency"))

    out = {}
    skipped = []

    # ── Loop distributions (one figure per condition) ───────────────────────
    for case in LOOP_CASES:
        d = _filter(loops, l_mat, l_type, l_dose, l_temp, l_bin, l_freq,
                    material=case["material"], irr_type=case.get("irr_type"),
                    dose=case.get("dose"), temp=case.get("temp"))
        x, w, y = _xywh(d, l_bin, l_bw, l_freq)
        if x.size == 0:
            skipped.append(("loop", case["title"]))
            continue

        stem = "loop_" + _slug(case["title"])
        fig, ax = plt.subplots(figsize=figsize)
        ax.bar(x, y, width=w * 0.85, color=case["color"],
               alpha=0.85, edgecolor="white", linewidth=0.6)
        if case.get("mean") is not None:
            ax.axvline(case["mean"], color=case["color"], lw=1.8, ls="--",
                       alpha=0.9, label=f"Mean = {case['mean']:.1f} nm")
            ax.legend(frameon=False)
        ax.set_xlabel("Loop diameter (nm)")
        ax.set_ylabel("Relative frequency (%)")
        ax.set_title(case["title"], fontweight="bold")
        ax.set_xlim(*case["xlim"])
        _annotate_cite(ax, case["citation"], cite_fontsize)
        fig.tight_layout()
        out[stem] = fig

    # ── Void distributions (one figure per condition) ───────────────────────
    for case in VOID_CASES:
        d = _filter(voids, v_mat, v_type, v_dose, v_temp, v_bin, v_freq,
                    material=case["material"], irr_type=case.get("irr_type"),
                    dose=case.get("dose"), temp=case.get("temp"))
        x, w, y = _xywh(d, v_bin, v_bw, v_freq)
        if x.size == 0:
            skipped.append(("void", case["title"]))
            continue

        stem = "void_" + _slug(case["title"])
        fig, ax = plt.subplots(figsize=figsize)
        if case.get("bimodal"):
            small = x < 2.2
            large = ~small
            ax.bar(x[small], y[small], width=w[small] * 0.85,
                   color="#2E86C1", alpha=0.85, edgecolor="white", lw=0.6,
                   label="He-bubbles ($r < r_c$)")
            ax.bar(x[large], y[large], width=w[large] * 0.85,
                   color=case["color"], alpha=0.85, edgecolor="white", lw=0.6,
                   label="Growing voids ($r > r_c$)")
            ax.axvline(2.0, color="gray", lw=1.2, ls=":", alpha=0.8,
                       label="Approx. critical radius")
            ax.legend(frameon=False)
        else:
            ax.bar(x, y, width=w * 0.85, color=case["color"],
                   alpha=0.85, edgecolor="white", linewidth=0.6)
        ax.set_xlabel("Cavity diameter (nm)")
        ax.set_ylabel("Relative frequency (%)")
        ax.set_title(case["title"], fontweight="bold")
        ax.set_xlim(*case["xlim"])
        _annotate_cite(ax, case["citation"], cite_fontsize)
        fig.tight_layout()
        out[stem] = fig

    # ── Overview scatter (kept as one figure) ───────────────────────────────
    if include_overview:
        try:
            summary_df = pd.read_excel(xlsx, sheet_name="Loop Size Distributions", header=3)
            summary_df.columns = [str(c).strip() for c in summary_df.columns]
            summary_df = summary_df.dropna(how="all")

            mean_col  = [c for c in summary_df.columns if "Mean Diam" in c][0]
            temp_col2 = [c for c in summary_df.columns if "Temp" in c][0]
            mat_col2  = [c for c in summary_df.columns if "Material" in c][0]
            type_col2 = [c for c in summary_df.columns if "Irrad. Type" in c][0]

            fig, ax = plt.subplots(figsize=(figsize[0] + 2, figsize[1] + 0.5))
            materials = ["EUROFER97", "T91", "HT9", "Pure Fe", "Fe-9Cr"]
            markers   = ["o", "s", "^", "D", "v"]
            cols_m    = ["#1B4F8A", "#117A65", "#922B21", "#D4AC0D", "#6C3483"]

            plotted_any = False
            for mat, mk, c in zip(materials, markers, cols_m):
                sub = summary_df[summary_df[mat_col2].astype(str).str.contains(mat, case=False)]
                for _, row in sub.iterrows():
                    try:
                        t = float(str(row[temp_col2]).replace("~", ""))
                        m_raw = str(row[mean_col]).replace("~", "").replace("+", "").split()[0]
                        m = float(m_raw)
                        is_ion = "ion" in str(row[type_col2]).lower()
                        fc = "none" if is_ion else c
                        ax.scatter(t, m, marker=mk, color=c, facecolors=fc,
                                   s=80, zorder=5, linewidths=1.5)
                        plotted_any = True
                    except Exception:
                        pass

            if plotted_any:
                legend_patches = [plt.scatter([], [], marker=mk, color=c, s=70, label=mat)
                                  for mat, mk, c in zip(materials, markers, cols_m)]
                legend_patches += [
                    mpatches.Patch(facecolor="gray", edgecolor="gray", alpha=0.8,
                                   label="Filled = neutron"),
                    mpatches.Patch(facecolor="none", edgecolor="gray", linewidth=1.5,
                                   label="Open = ion"),
                ]
                ax.legend(handles=legend_patches, frameon=False, ncol=2, loc="upper left")
                ax.set_xlabel("Irradiation temperature (°C)")
                ax.set_ylabel("Mean loop diameter (nm)")
                ax.set_title("Mean Dislocation Loop Diameter vs Irradiation Temperature",
                             fontweight="bold")
                ax.set_xlim(50, 650); ax.set_ylim(0, 35)
                ax.grid(True, ls=":", alpha=0.4)
                fig.tight_layout()
                out["overview_mean_loop_size_vs_temperature"] = fig
            else:
                plt.close(fig)
                skipped.append(("overview", "mean loop size vs T (no parseable rows)"))
        except Exception as e:
            skipped.append(("overview", f"failed: {e}"))

    if save:
        for stem, fig in out.items():
            _save(fig, outdir, stem, formats=formats)
        print(f"Saved {len(out)} figure(s) to {os.path.abspath(outdir)}/")

    if skipped:
        print(f"Skipped {len(skipped)} figure(s) (no data):")
        for kind, title in skipped:
            print(f"  [{kind}] {title}")

    if not show:
        for fig in out.values():
            plt.close(fig)

    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlsx", default="Ferritics_Irradiated_Microstructure_Data.xlsx")
    parser.add_argument("--outdir", default="figures")
    parser.add_argument("--font-size", type=int, default=16)
    args = parser.parse_args()
    plot_size_distributions(xlsx=args.xlsx, outdir=args.outdir,
                            font_size=args.font_size, save=True, show=False)


if __name__ == "__main__":
    main()
