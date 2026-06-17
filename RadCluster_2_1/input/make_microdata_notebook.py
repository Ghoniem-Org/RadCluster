#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate EuroferMicroData.ipynb (density/size maps from EuroferMicrostructure.xlsx)."""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

nb = new_notebook()
cells = []

cells.append(new_markdown_cell(
"# EUROFER / F-M Steel Irradiated Microstructure — Density & Size Maps\n"
"\n"
"Visualises the consolidated experimental database **`EuroferMicrostructure.xlsx`** as\n"
"(dose, temperature) maps. Each panel draws **approximate iso-value contours** (filled\n"
"pale bands of constant defect density or size, obtained by triangulating the scattered\n"
"data) and overlays the **experimental points colour-coded by literature source**.\n"
"\n"
"* x-axis: dose (dpa, log scale)  *  y-axis: irradiation temperature (°C, 250-550)\n"
"* colour = reference,  marker = material;  legends are placed beneath each figure.\n"
"\n"
"Panels: (1) neutron loop density, (2) neutron loop size, (3) neutron void/cavity\n"
"density, (4) radiation-induced precipitate number density, (5) ion-irradiation loop\n"
"density and (6) ion-irradiation void/cavity density."
))

# ---------------- setup cell ----------------
setup = r'''
import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ---- global style: size-16 fonts ----
plt.rcParams.update({
    "font.size": 16, "axes.titlesize": 16, "axes.labelsize": 16,
    "xtick.labelsize": 14, "ytick.labelsize": 14, "figure.dpi": 110,
})

# ---- locate the consolidated workbook ----
CANDIDATES = ["../../input/EuroferMicrostructure.xlsx",
              "../input/EuroferMicrostructure.xlsx",
              "input/EuroferMicrostructure.xlsx",
              "EuroferMicrostructure.xlsx"]
XLSX = next((p for p in CANDIDATES if os.path.exists(p)), CANDIDATES[0])
print("Reading:", os.path.abspath(XLSX))

# ---- numeric parsers (data cells contain ~, <, >, ranges and superscripts) ----
def parse_num(x, geom=False):
    """Plain numbers / ranges (a-b) / '~x' / '<x' / '>x'. Range -> mean (or geom)."""
    if x is None:
        return np.nan
    if isinstance(x, (int, float)):
        return np.nan if (isinstance(x, float) and np.isnan(x)) else float(x)
    s = str(x).strip()
    if not s:
        return np.nan
    s = (s.replace("−", "-").replace("–", "-").replace("—", "-")
           .replace("≤", "").replace("≥", "").replace("~", "")
           .replace("<", "").replace(">", ""))
    s = re.sub(r"[×x]\s*10\s*[⁰-⁹\^\-\d]*", "", s)   # strip sci-notation tails
    nums = [float(n) for n in re.findall(r"\d+\.?\d*", s) if n not in ("", ".")]
    if not nums:
        return np.nan
    if len(nums) == 1:
        return nums[0]
    a, b = nums[0], nums[1]
    if geom and a > 0 and b > 0:
        return (a * b) ** 0.5
    return (a + b) / 2.0

_SUP = {"⁰":"0","¹":"1","²":"2","³":"3","⁴":"4",
        "⁵":"5","⁶":"6","⁷":"7","⁸":"8","⁹":"9","⁻":"-"}
def _desuper(t):
    return "".join(_SUP.get(c, c) for c in t)

def parse_sci(x):
    """Absolute densities in m^-3 written as 7.4x10^23, ~10^23, etc."""
    if x is None:
        return np.nan
    if isinstance(x, (int, float)):
        return np.nan if (isinstance(x, float) and np.isnan(x)) else float(x)
    s = str(x)
    m = re.search(r"(\d+\.?\d*)\s*[×x]\s*10\s*([⁰-⁹⁻\^\-]*\d*)", s)
    if m and m.group(2):
        try:
            return float(m.group(1)) * 10 ** int(_desuper(m.group(2)).replace("^", ""))
        except ValueError:
            pass
    m = re.search(r"10\s*([⁰-⁹⁻]+)", s)
    if m:
        try:
            return 10.0 ** int(_desuper(m.group(1)))
        except ValueError:
            pass
    n = re.findall(r"\d+\.?\d*", s)
    return float(n[0]) if n else np.nan

# value parsers
pdens = lambda v: parse_num(v, geom=True)   # density in units of 1e22 m^-3
psize = lambda v: parse_num(v, geom=False)  # diameter in nm
pprec = parse_sci                           # precipitate density in m^-3

# ---- citation map: Density-workbook [CIT-##] -> "Author Year" ----
_ref = pd.read_excel(XLSX, sheet_name="References", header=None).iloc[2:]
def _short(a):
    a = str(a).replace(" et al.", "").replace(" et al", "")
    a = re.split(r"\s*[&/]\s*", a)[0]   # keep first author (before & or /)
    return a.strip()
def _yr(y):
    try:
        return str(int(float(y)))
    except (ValueError, TypeError):
        return str(y)
DEN_MAP = {}
for _, r in _ref.iterrows():
    dk = r.iloc[7]
    if isinstance(dk, str) and dk.strip().upper().startswith("CIT"):
        DEN_MAP[dk.strip()] = f"{_short(r.iloc[1])} {_yr(r.iloc[2])}"

# ---- generic loader ----
def load(sheet, mat, dose, temp, val, ref, valparser=pdens,
         data_start=3, irrcol=None, irr=None, refmap=None):
    raw = pd.read_excel(XLSX, sheet_name=sheet, header=None).iloc[data_start:]
    df = pd.DataFrame({
        "material": raw.iloc[:, mat].astype(str).str.strip(),
        "dose":     raw.iloc[:, dose].map(parse_num),
        "temp":     raw.iloc[:, temp].map(parse_num),
        "value":    raw.iloc[:, val].map(valparser),
        "ref":      raw.iloc[:, ref].astype(str).str.strip(),
    })
    if irrcol is not None and irr is not None:
        # match by leading word so "Neutron fission" is NOT caught by irr="Ion"
        col = raw.iloc[:, irrcol].astype(str).str.strip().str.lower()
        keep = col.str.startswith(irr.lower())
        df = df[keep.values]
    if refmap is not None:
        df["ref"] = df["ref"].map(
            lambda k: refmap.get(str(k).strip().strip("[]"), str(k).strip()))
    df = df[df["material"].str.lower() != "nan"]
    return df.dropna(subset=["dose", "temp", "value"]).reset_index(drop=True)

# ---- iso-value label formatter ----
def _fmt_pow(v):
    e = int(np.floor(np.log10(v)))
    m = v / 10 ** e
    return (rf"$10^{{{e}}}$" if abs(m - 1) < 0.05
            else rf"${m:.0f}{{\times}}10^{{{e}}}$")

_MARKERS = ["o", "s", "^", "D", "v", "P", "*", "X", "<", ">", "h", "8"]

def make_field(df, levels, title, cbar_label, logz=True, dens_scale=1.0,
               cmap="YlGnBu", xlim=(0.08, 900)):
    """Filled iso-value contours (triangulated) + reference-coloured scatter."""
    df = df.copy()
    df["x"] = np.log10(df["dose"])
    fig, ax = plt.subplots(figsize=(11, 8))

    # contour from points inside the plotting window (avoid distortion)
    cdf = df[df["temp"].between(240, 560)]
    agg = cdf.groupby(["x", "temp"], as_index=False)["value"].mean()
    levels = np.array(sorted(levels), dtype=float)
    L = np.log10(levels) if logz else levels
    if len(agg) >= 3:
        zz = agg["value"].values * dens_scale
        Z = np.log10(zz) if logz else zz
        try:
            tcf = ax.tricontourf(agg["x"], agg["temp"], Z, levels=L,
                                 cmap=cmap, alpha=0.40, extend="both")
            ax.tricontour(agg["x"], agg["temp"], Z, levels=L,
                          colors="0.45", linewidths=0.8)
            cb = fig.colorbar(tcf, ax=ax, pad=0.02, fraction=0.046)
            cb.set_ticks(L)
            cb.set_ticklabels([_fmt_pow(v) for v in levels] if logz
                              else [f"{v:g}" for v in levels])
            cb.set_label(cbar_label)
        except Exception as exc:
            print("  [contour skipped]", exc)

    # scatter: colour = reference, marker = material
    refs = sorted(df["ref"].unique())
    mats = sorted(df["material"].unique())
    cmap_r = plt.get_cmap("tab20", max(len(refs), 1))
    rc = {r: cmap_r(i) for i, r in enumerate(refs)}
    mk = {m: _MARKERS[i % len(_MARKERS)] for i, m in enumerate(mats)}
    for (r, m), s in df.groupby(["ref", "material"]):
        ax.scatter(s["x"], s["temp"], color=rc[r], marker=mk[m], s=150,
                   edgecolor="k", linewidth=0.7, zorder=5)

    ax.set_xlim(np.log10(xlim[0]), np.log10(xlim[1]))
    ax.set_ylim(250, 550)
    ticks = [0.1, 1, 10, 100, 1000]
    ax.set_xticks([np.log10(t) for t in ticks])
    ax.set_xticklabels([str(t) for t in ticks])
    ax.set_xlabel("Dose (dpa)")
    ax.set_ylabel(r"Irradiation temperature (°C)")
    ax.set_title(title)
    ax.grid(alpha=0.25, linestyle=":")

    # legends underneath
    ref_h = [mlines.Line2D([], [], marker="o", ls="", color=rc[r], mec="k",
                           ms=11, label=r) for r in refs]
    mat_h = [mlines.Line2D([], [], marker=mk[m], ls="", color="0.55", mec="k",
                           ms=11, label=m) for m in mats]
    fig.subplots_adjust(bottom=0.40)
    nref_rows = int(np.ceil(len(refs) / 4)) + 1   # +1 for the title row
    leg1 = ax.legend(handles=ref_h, title="Reference (source)", loc="upper left",
                     bbox_to_anchor=(0.0, -0.12), ncol=4, fontsize=12,
                     title_fontsize=13, frameon=True)
    ax.add_artist(leg1)
    ax.legend(handles=mat_h, title="Material", loc="upper left",
              bbox_to_anchor=(0.0, -0.12 - 0.075 * nref_rows), ncol=6,
              fontsize=12, title_fontsize=13, frameon=True)
    print(f"{title}: {len(df)} points, {len(refs)} references")
    return fig, ax
'''
cells.append(new_code_cell(setup.strip("\n")))

# ---- (1) neutron loop density ----
cells.append(new_markdown_cell(
"## (1) Neutron irradiation — dislocation-loop number density\n"
"Source sheet `Loops_DenSize` (rows with `Irrad. Type = Neutron*`). Density is reported "
"in units of $10^{22}\\,$m$^{-3}$; contour bands are drawn at constant absolute density."
))
cells.append(new_code_cell(
'df = load("Loops_DenSize", mat=0, dose=7, temp=8, val=10, ref=19,\n'
'          valparser=pdens, irrcol=5, irr="Neutron", refmap=DEN_MAP)\n'
'levels = [1e20, 3e20, 1e21, 3e21, 1e22, 3e22, 1e23]\n'
'make_field(df, levels,\n'
'           "Neutron — dislocation-loop density",\n'
'           r"Loop density $N_d$ (m$^{-3}$)",\n'
'           logz=True, dens_scale=1e22)\n'
'plt.show()'
))

# ---- (2) neutron loop size ----
cells.append(new_markdown_cell(
"## (2) Neutron irradiation — mean dislocation-loop diameter\n"
"Same conditions as panel (1); contour bands are drawn at constant **mean loop diameter** (nm)."
))
cells.append(new_code_cell(
'df = load("Loops_DenSize", mat=0, dose=7, temp=8, val=15, ref=19,\n'
'          valparser=psize, irrcol=5, irr="Neutron", refmap=DEN_MAP)\n'
'levels = [2, 4, 6, 8, 10, 14, 18, 25]\n'
'make_field(df, levels,\n'
'           "Neutron — mean loop diameter",\n'
'           r"Mean loop diameter $\\bar d$ (nm)",\n'
'           logz=False)\n'
'plt.show()'
))

# ---- (3) neutron voids ----
cells.append(new_markdown_cell(
"## (3) Neutron irradiation — void / cavity number density\n"
"Source sheet `Voids_DenSize` (total cavity density, $10^{22}\\,$m$^{-3}$). "
"Below ~415 °C the population is unimodal (He-stabilised bubbles); a bimodal "
"bubble+void structure develops at higher dose/temperature."
))
cells.append(new_code_cell(
'df = load("Voids_DenSize", mat=0, dose=7, temp=8, val=10, ref=20,\n'
'          valparser=pdens, irrcol=5, irr="Neutron", refmap=DEN_MAP)\n'
'levels = [3e20, 1e21, 3e21, 1e22, 3e22, 1e23]\n'
'make_field(df, levels,\n'
'           "Neutron — void / cavity density",\n'
'           r"Cavity density $N_v$ (m$^{-3}$)",\n'
'           logz=True, dens_scale=1e22)\n'
'plt.show()'
))

# ---- (4) precipitates ----
cells.append(new_markdown_cell(
"## (4) Neutron irradiation — radiation-induced precipitate number density\n"
"Source sheet `N_Precipitates` (number density already in m$^{-3}$). Phases include "
"α′ (Cr-rich bcc), G-phase (Mn/Ni/Si), MX (TaC, VN), M$_{23}$C$_6$ and Si-rich clusters; "
"the reference column already carries author–year labels."
))
cells.append(new_code_cell(
'df = load("N_Precipitates", mat=0, dose=6, temp=7, val=11, ref=15,\n'
'          valparser=pprec, refmap=None)\n'
'levels = [1e21, 3e21, 1e22, 3e22, 1e23, 3e23, 1e24]\n'
'make_field(df, levels,\n'
'           "Neutron — precipitate number density",\n'
'           r"Precipitate density (m$^{-3}$)",\n'
'           logz=True, dens_scale=1.0)\n'
'plt.show()'
))

# ---- (5) ion loop density ----
cells.append(new_markdown_cell(
"## (5) Ion irradiation — dislocation-loop number density\n"
"Source sheet `Loops_DenSize` (rows with `Irrad. Type = Ion*`): self-ion, dual-ion and "
"in-situ Kr data. Note the much higher loop densities reached under self-ion irradiation "
"(e.g. pure Fe, Fe-9Cr) and the ~+60–70 °C temperature shift relative to neutron data."
))
cells.append(new_code_cell(
'df = load("Loops_DenSize", mat=0, dose=7, temp=8, val=10, ref=19,\n'
'          valparser=pdens, irrcol=5, irr="Ion", refmap=DEN_MAP)\n'
'levels = [1e21, 3e21, 1e22, 3e22, 1e23, 3e23]\n'
'make_field(df, levels,\n'
'           "Ion — dislocation-loop density",\n'
'           r"Loop density $N_d$ (m$^{-3}$)",\n'
'           logz=True, dens_scale=1e22)\n'
'plt.show()'
))

# ---- (6) ion voids ----
cells.append(new_markdown_cell(
"## (6) Ion irradiation — void / cavity number density\n"
"Source sheet `Voids_DenSize` (rows with `Irrad. Type = Ion*`): dual-ion T91/F82H and "
"high-dose self-ion HT9. Peak swelling (bimodal cavities) occurs near 450–470 °C."
))
cells.append(new_code_cell(
'df = load("Voids_DenSize", mat=0, dose=7, temp=8, val=10, ref=20,\n'
'          valparser=pdens, irrcol=5, irr="Ion", refmap=DEN_MAP)\n'
'levels = [3e20, 1e21, 3e21, 1e22, 3e22, 1e23]\n'
'make_field(df, levels,\n'
'           "Ion — void / cavity density",\n'
'           r"Cavity density $N_v$ (m$^{-3}$)",\n'
'           logz=True, dens_scale=1e22)\n'
'plt.show()'
))

nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
}
out = "../codes/Notebooks/EuroferMicroData.ipynb"
with open(out, "w", encoding="utf-8") as f:
    nbf.write(nb, f)
print("wrote", out, "with", len(cells), "cells")
