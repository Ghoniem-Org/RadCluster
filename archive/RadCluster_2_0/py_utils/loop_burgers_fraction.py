"""
loop_burgers_fraction.py
Extracts and plots the fraction of ½⟨111⟩ dislocation loops vs temperature
and irradiation dose from Ferritics_Irradiated_Microstructure_Data.xlsx.

Two data sources are used:
  1. The `Dominant Burgers Vector` column of the
     `N – Dislocation Loops` and `Ion – Dislocation Loops` sheets.
     Free-text descriptions are parsed into numeric f₁₁₁ ∈ [0, 1].
  2. The `111 Fraction` sheet, which contains size-binned histograms.
     For conditions where both ⟨100⟩ and ½⟨111⟩ histograms are present
     (e.g. 350°C / 16 dpa), an integrated bulk fraction is computed.

Programmatic use (Jupyter):
    from loop_burgers_fraction import plot_loop_111_fraction
    figs, df = plot_loop_111_fraction(
        xlsx="../../input/Ferritics_Irradiated_Microstructure_Data.xlsx",
        outdir="../../output/microstructure",
        font_size=16, save=True, show=True,
    )
"""

import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# ── Burgers-vector text parser ──────────────────────────────────────────────
# Returns f₁₁₁ ∈ [0, 1] (or None if unparseable).  May depend on T and dose
# when the source string encodes a transition (e.g. "½⟨111⟩→⟨100⟩").

def _has(s, *needles):
    return any(n in s for n in needles)

def parse_burgers_fraction(bv_str, temp_c=None, dose_dpa=None):
    if bv_str is None or (isinstance(bv_str, float) and np.isnan(bv_str)):
        return None
    s = str(bv_str).strip()
    s_low = s.lower()

    # ── 1. Explicit numeric percentage hints, e.g. "(~40% ⟨100⟩)"
    m = re.search(r"~?(\d{1,3})\s*%\s*(?:⟨|<|a)?\s*1?00", s_low)
    if m:
        f100 = float(m.group(1)) / 100.0
        return max(0.0, min(1.0, 1.0 - f100))
    m = re.search(r"~?(\d{1,3})\s*%\s*(?:⟨|<|a)?\s*111", s_low)
    if m:
        return max(0.0, min(1.0, float(m.group(1)) / 100.0))

    # ── 2. Temperature-conditional transitions
    # e.g. "½⟨111⟩ (≤300°C)→⟨100⟩ (≥350°C)" or "½⟨111⟩→⟨100⟩"
    if "→" in s and "111" in s and "100" in s:
        if temp_c is None:
            return 0.5
        t_lo_match = re.search(r"≤\s*(\d{2,4})\s*°?c", s_low)
        t_hi_match = re.search(r"≥\s*(\d{2,4})\s*°?c", s_low)
        t_lo = float(t_lo_match.group(1)) if t_lo_match else 300.0
        t_hi = float(t_hi_match.group(1)) if t_hi_match else 350.0
        if temp_c <= t_lo: return 1.0
        if temp_c >= t_hi: return 0.0
        return float((t_hi - temp_c) / max(1.0, (t_hi - t_lo)))

    # ── 3. Dose-conditional transitions
    # e.g. "a/2⟨111⟩ dom. (<4 dpa); a⟨100⟩ dom. (>4 dpa)"
    if ("dpa" in s_low) and ("<" in s_low or ">" in s_low) and "111" in s_low and "100" in s_low:
        if dose_dpa is None:
            return 0.5
        d_match = re.search(r"<\s*(\d+(?:\.\d+)?)\s*dpa", s_low)
        d_thr = float(d_match.group(1)) if d_match else 4.0
        return 1.0 if dose_dpa <= d_thr else 0.15

    has_111 = "111" in s_low
    has_100 = "100" in s_low

    # ── 4. Single-population entries
    if has_111 and not has_100:
        return 1.0
    if has_100 and not has_111:
        return 0.0

    # ── 5. Mixed population — qualify by "dominant" wording
    if has_111 and has_100:
        # Find which symbol comes first / which is followed by "dom"
        # split on common separators; check which fragment contains "dom"
        m100_dom = bool(re.search(r"100[^\.;]*?\bdom", s_low))
        m111_dom = bool(re.search(r"111[^\.;]*?\bdom", s_low))
        if m100_dom and not m111_dom: return 0.15
        if m111_dom and not m100_dom: return 0.85
        return 0.5  # mixed, no dominance qualifier

    return None


# ── Build a tidy DataFrame of f_111 vs (T, dose) ─────────────────────────────
def extract_f111_table(xlsx):
    rows = []

    for sheet, irr_kind in (("N – Dislocation Loops",   "Neutron"),
                            ("Ion – Dislocation Loops", "Ion")):
        df = pd.read_excel(xlsx, sheet_name=sheet, header=[1, 2])
        df.columns = [str(c2).replace("\n", " ").strip() for (_, c2) in df.columns]
        for _, r in df.iterrows():
            mat = r.get("Material")
            if pd.isna(mat):
                continue
            T = pd.to_numeric(str(r.get("Temp. (°C)", "")).split("–")[0].strip(), errors="coerce")
            D = pd.to_numeric(str(r.get("Dose (dpa)", "")).split("–")[0].strip(), errors="coerce")
            f = parse_burgers_fraction(r.get("Dominant Burgers Vector"),
                                       temp_c=None if pd.isna(T) else float(T),
                                       dose_dpa=None if pd.isna(D) else float(D))
            if f is None or pd.isna(T):
                continue
            rows.append(dict(material=str(mat).strip(),
                             irradiation=irr_kind,
                             dose_dpa=float(D) if not pd.isna(D) else np.nan,
                             temperature_C=float(T),
                             f_111=float(f),
                             source=sheet,
                             raw_bv=str(r.get("Dominant Burgers Vector"))))

    # ── 111 Fraction sheet: integrated bulk fraction at conditions where both
    # ⟨100⟩ and ½⟨111⟩ histograms are present (e.g. 350°C / 16 dpa).
    try:
        raw = pd.read_excel(xlsx, sheet_name="111 Fraction", header=None)
        # Row 0 has block labels; row 1 has 'Diameter'/'Fraction'; data starts row 2.
        labels = [str(raw.iat[0, j]) if pd.notna(raw.iat[0, j]) else "" for j in range(raw.shape[1])]
        # Group columns by label header (label appears in the first col of each block)
        blocks = {}
        for j, lab in enumerate(labels):
            if lab and lab.lower() != "nan":
                blocks[lab] = (j, j + 1)  # (Diameter col, Fraction col)
        # Sum each block's fractions to get total counts (proxy for relative density)
        block_sum = {}
        for lab, (jd, jf) in blocks.items():
            vals = pd.to_numeric(raw.iloc[2:, jf], errors="coerce").dropna().values
            block_sum[lab] = vals.sum()
        # If both 350C_16dpa_100 and _111 are present, compute integrated f_111
        if "350C_16dpa_100" in block_sum and "350C_16dpa_111" in block_sum:
            s100 = block_sum["350C_16dpa_100"]
            s111 = block_sum["350C_16dpa_111"]
            f = s111 / (s100 + s111) if (s100 + s111) > 0 else np.nan
            rows.append(dict(material="EUROFER97 (111 sheet)",
                             irradiation="Neutron",
                             dose_dpa=16.0,
                             temperature_C=350.0,
                             f_111=float(f),
                             source="111 Fraction (integrated)",
                             raw_bv="histogram-derived"))
    except Exception:
        pass

    out = pd.DataFrame(rows)
    return out


# ── Plotting ─────────────────────────────────────────────────────────────────
def _apply_style(font_size, font_family, font_serif, axes_linewidth, dpi, save_dpi):
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


def _r2(x, y):
    """Coefficient of determination for a least-squares linear fit."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    if x.size < 3 or np.std(x) == 0:
        return 0.0
    a, b = np.polyfit(x, y, 1)
    yhat = a * x + b
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 0.0 if ss_tot == 0 else 1.0 - ss_res / ss_tot


def _spearman_abs(x, y):
    """|Spearman ρ| — captures monotonic (incl. nonlinear) dependence; better
    than R² when the underlying relationship is sigmoidal or step-like."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    if x.size < 3:
        return 0.0
    rx = pd.Series(x).rank().values
    ry = pd.Series(y).rank().values
    if np.std(rx) == 0 or np.std(ry) == 0:
        return 0.0
    return float(abs(np.corrcoef(rx, ry)[0, 1]))


def _scatter(ax, df, xkey, xlabel, color_by="material", marker_by="irradiation"):
    materials = sorted(df["material"].unique())
    cmap = plt.get_cmap("tab10")
    color_map = {m: cmap(i % 10) for i, m in enumerate(materials)}
    marker_map = {"Neutron": "o", "Ion": "^"}

    for mat in materials:
        sub = df[df["material"] == mat]
        for irr, mk in marker_map.items():
            ssub = sub[sub["irradiation"] == irr]
            if ssub.empty:
                continue
            ax.scatter(ssub[xkey], 100.0 * ssub["f_111"],
                       color=color_map[mat], marker=mk,
                       s=110, edgecolors="black", linewidths=0.6,
                       alpha=0.85, label=f"{mat} ({irr})")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Fraction of ½⟨111⟩ loops (%)")
    ax.set_ylim(-5, 105)
    ax.grid(True, ls=":", alpha=0.4)


def plot_loop_111_fraction(
    xlsx="Ferritics_Irradiated_Microstructure_Data.xlsx",
    outdir="figures",
    save=True,
    show=False,
    formats=("pdf", "png"),
    weak_threshold=0.20,  # |Spearman ρ| below this is called "weak"
    force=None,           # "T", "dose", "both", or None for auto
    font_size=16,
    font_family="serif",
    font_serif=("Times New Roman", "DejaVu Serif"),
    axes_linewidth=1.2,
    dpi=150,
    save_dpi=300,
    figsize=(10, 6),
):
    """
    Build figures of f_111 vs (T, dose).

    Strategy:
      - Compute |Spearman ρ| for f_111 vs T and vs dose (rank-based;
        captures monotonic non-linear trends — e.g. the sigmoidal
        ½⟨111⟩→⟨100⟩ transition with T).
      - If both > weak_threshold, plot both.
      - If only one is, plot only that one (the weak axis is omitted,
        per the user's spec).
      - If neither is, plot both anyway as a diagnostic.
      - `force` overrides the auto choice.

    Returns (figs_dict, dataframe).
    """
    _apply_style(font_size, font_family, font_serif, axes_linewidth, dpi, save_dpi)
    if save:
        os.makedirs(outdir, exist_ok=True)

    df = extract_f111_table(xlsx)
    if df.empty:
        print("No parseable f_111 data found.")
        return {}, df

    rho_T    = _spearman_abs(df["temperature_C"], df["f_111"])
    rho_dose = _spearman_abs(df["dose_dpa"],      df["f_111"])
    r2_T     = _r2(df["temperature_C"], df["f_111"])
    r2_dose  = _r2(df["dose_dpa"],      df["f_111"])
    print(f"Parsed {len(df)} data points")
    print(f"   |Spearman ρ|(T)    = {rho_T:.3f}   (linear R² = {r2_T:.3f})")
    print(f"   |Spearman ρ|(dose) = {rho_dose:.3f}   (linear R² = {r2_dose:.3f})")
    print(f"   weak threshold     = {weak_threshold}")

    if force == "T":
        which = ("T",)
    elif force == "dose":
        which = ("dose",)
    elif force == "both":
        which = ("T", "dose")
    else:
        T_strong    = rho_T    >= weak_threshold
        dose_strong = rho_dose >= weak_threshold
        if T_strong and dose_strong:
            which = ("T", "dose")
        elif T_strong:
            which = ("T",)
        elif dose_strong:
            which = ("dose",)
        else:
            print("   → both dependences look weak; plotting both as a diagnostic.")
            which = ("T", "dose")

    out = {}

    if "T" in which:
        fig, ax = plt.subplots(figsize=figsize)
        _scatter(ax, df, "temperature_C", "Irradiation temperature (°C)")
        ax.set_title(f"Fraction of ½⟨111⟩ loops vs Temperature  "
                     f"(|ρ| = {rho_T:.2f})", fontweight="bold")
        ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5),
                  fontsize=max(8, font_size - 4))
        fig.tight_layout()
        out["loop_111_fraction_vs_temperature"] = fig

    if "dose" in which:
        fig, ax = plt.subplots(figsize=figsize)
        _scatter(ax, df, "dose_dpa", "Dose (dpa)")
        ax.set_xscale("log")
        ax.set_title(f"Fraction of ½⟨111⟩ loops vs Dose  "
                     f"(|ρ| = {rho_dose:.2f})", fontweight="bold")
        ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5),
                  fontsize=max(8, font_size - 4))
        fig.tight_layout()
        out["loop_111_fraction_vs_dose"] = fig

    if save:
        for stem, fig in out.items():
            for ext in formats:
                fig.savefig(os.path.join(outdir, f"{stem}.{ext}"))
        print(f"Saved {len(out)} figure(s) to {os.path.abspath(outdir)}/")

    if not show:
        for fig in out.values():
            plt.close(fig)

    return out, df


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--xlsx", default="Ferritics_Irradiated_Microstructure_Data.xlsx")
    p.add_argument("--outdir", default="figures")
    p.add_argument("--font-size", type=int, default=16)
    a = p.parse_args()
    plot_loop_111_fraction(xlsx=a.xlsx, outdir=a.outdir,
                           font_size=a.font_size, save=True, show=False)
