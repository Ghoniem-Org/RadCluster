"""Figures for `experimental_database.tex` -- the group-separated analysis.

The manuscript section treats **neutron** and **ion** irradiation as two separate
populations and never pools them, so the panels produced here are not the ones in
`FerriticSteelsMicroData.ipynb` Sections 5-8 (which fit all irradiation types
together).  Everything upstream of the fitting -- cell parsing, the ion
neutron-equivalent temperature shift, reference numbering, the log-linear
machinery and the loop size-distribution fits -- is *executed out of the
notebook* rather than copied, so the two can never drift apart.  The cells are
selected by content (see `CELL_MARKERS`), not by index, so they survive editing
and re-ordering of the notebook.

Run from anywhere:

    python make_experiment_figures.py

Writes `figures/experiments/*.png` and `figures/experiments/fit_summary.txt`.
"""

import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402  (must follow the backend selection)

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")   # the notebook prints Burgers-vector glyphs

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)                     # the notebook opens the workbook by relative path
NOTEBOOK = "FerriticSteelsMicroData.ipynb"
OUTDIR = os.path.join("figures", "experiments")
os.makedirs(OUTDIR, exist_ok=True)
DPI = 300          # raster resolution of every saved panel


def save(fig, fname):
    """Write one panel and release it.  300 dpi on an 11x8 in canvas is 3300x2400 px."""
    fig.savefig(os.path.join(OUTDIR, fname), dpi=DPI)
    plt.close(fig)

# A code cell is executed if it contains any of these strings.  One marker per
# cell we depend on; each is a definition that would have to be renamed for the
# marker to go stale, at which point the KeyError below is immediate and loud.
CELL_MARKERS = (
    "XLSX = ",                    # workbook path, fonts, T_SHIFT_ION, KB
    "def parse_density",          # free-text cell parsers
    "def build_database",         # sheet map and database assembly -> DB
    "REFERENCES = [",             # reference numbering -> DB.ref_no
    "def loglin_fit",             # log-linear fit, legend placement, markers, select()
    "def load_size_distributions",  # loop size-distribution histograms -> META, HIST, FITS
    "def mfit(",                  # multi-covariate fit, LOO, nested F test
    "def alloy_family",           # family / Cr columns and covariate builders
    "def load_loop_fractions",    # <100> fractions, logit fit, panels -> LOOPFRAC
)


def load_notebook_machinery():
    """Execute the notebook's analysis cells into this module's namespace."""
    nb = json.load(open(NOTEBOOK, encoding="utf-8"))
    n = 0
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        if any(m in source for m in CELL_MARKERS):
            exec(compile(source, f"<{NOTEBOOK}>", "exec"), globals())
            n += 1
    missing = [k for k in ("DB", "mfit", "loglin_fit", "FITS", "LOOPFRAC")
               if k not in globals()]
    if missing:
        raise RuntimeError(f"notebook cells did not define {missing}; check CELL_MARKERS")
    print(f"\nexecuted {n} notebook cells\n")


load_notebook_machinery()

# ---------------------------------------------------------------------------
# typography
# ---------------------------------------------------------------------------
# Manuscript figures are reproduced at a fraction of their native width, so they
# carry larger type than the notebook's inline panels (20 pt / 10 pt).  Everything
# drawn on the axes -- axis labels, tick labels and the reference tags on the
# points -- is set at FS_PLOT; only the legend box is smaller, at FS_LEGEND, since
# it can carry six or more entries and would otherwise dominate the panel.  These
# settings are applied here, after the notebook cells have run, so the notebook's
# own inline figures are unaffected.
FS_PLOT, FS_LEGEND = 24, 14
plt.rcParams.update({
    "font.size": FS_PLOT, "axes.labelsize": FS_PLOT, "axes.titlesize": FS_PLOT,
    "xtick.labelsize": FS_PLOT, "ytick.labelsize": FS_PLOT,
    "legend.fontsize": FS_LEGEND,
})

# The reference tags on the data points are set separately, and smaller.  They are
# drawn one per database row rather than one per panel, so on the dense panels --
# the coverage map, and the 573 K and 650 K columns of the neutron loop panels --
# neighbouring tags collide at the axis size.  16 pt separates the tags that are
# merely close while staying readable against 24 pt axes.  It cannot separate rows
# measured at the same temperature and value, whose markers coincide: those tags
# overlap at any font size and are resolved by the caption, not the panel.
FS_POINT_TAG = 16
PT_LABEL_FS = FS_POINT_TAG   # the notebook's 13 pt, overridden for these figures

# ---------------------------------------------------------------------------
# subsets and groupings
# ---------------------------------------------------------------------------
# Spallation (SINQ, mixed n+p+He) is a neutron-spectrum irradiation, so its rows
# are shown on the neutron panels -- but there are at most two of them per panel
# and their He/dpa ratio is an order of magnitude above any fission spectrum, so
# they are never included in a fit.  Ion panels are ion rows only.
NEUTRON_SERIES = ("Neutron", "Spallation")
GRPCOL = ["#1f4e79", "#c0392b", "#27ae60", "#8e44ad", "#e67e22", "#16a085"]


def subset(defect, quantity, irradiations):
    """Rows of one panel: finite value, finite equivalent temperature, sorted by T."""
    col = QCOL[quantity]                                    # noqa: F821 (from notebook)
    d = DB[(DB.defect == defect)                            # noqa: F821
           & DB.irradiation.isin(irradiations)
           & np.isfinite(DB[col])                           # noqa: F821
           & np.isfinite(DB["T_K_eq"])].copy()              # noqa: F821
    return d.sort_values(["T_K_eq", "material"]).reset_index(drop=True)


def dose_labels(d, thr):
    """Three-way dose label; rows without a reported dose are kept out of the fit."""
    v = d["dose"].to_numpy(float)
    lo, hi = rf"$\leq$ {thr:g} dpa", f"> {thr:g} dpa"
    return np.where(~np.isfinite(v), "dose not reported", np.where(v > thr, hi, lo)), lo


def fit_and_plot(d, quantity, xmode, ylabel, fname, labels=None, ref_level=None,
                 ylog=True, drop_label=None, note="", group_fit=True):
    """One panel: group-intercept / common-slope log-linear fit, drawn and saved.

    `labels` selects the grouping (None = a single group).  With `group_fit=False`
    the labels only colour the markers and a single common line is fitted -- used
    where a partition is worth showing but demonstrably does not reduce the scatter,
    so that the figure does not imply a structure the statistics reject.

    `drop_label` names a group that is plotted but excluded from the fit -- used for
    rows whose dose was not reported, which cannot enter a dose-grouped model without
    inventing a value.  Spallation rows are always excluded from the fit: there are at
    most two per panel and their He/dpa ratio is an order of magnitude above any
    fission spectrum.
    """
    col = QCOL[quantity]                                    # noqa: F821
    y = d[col].to_numpy(float)
    x = (d["T_K_eq"].to_numpy(float) if xmode == "T"
         else 1000.0 / d["T_K_eq"].to_numpy(float))
    labels = (np.array(["all"] * len(d), dtype=object) if labels is None
              else np.asarray(labels, dtype=object))

    use = d["irradiation"].ne("Spallation").to_numpy().copy()
    if drop_label is not None:
        use &= labels != drop_label
    lv_fit = [l for l in sorted(pd.unique(labels[use]).tolist(), key=str)]
    if ref_level in lv_fit:
        lv_fit = [ref_level] + [l for l in lv_fit if l != ref_level]

    dums = [(labels[use] == l).astype(float) for l in lv_fit[1:]] if group_fit else []
    m = mfit(y[use], [x[use]] + dums, "grouped")            # noqa: F821
    base = mfit(y[use], [x[use]], "T only") if dums else m  # noqa: F821

    fig, ax = plt.subplots(figsize=(11, 8))
    lv_all = lv_fit + [l for l in sorted(pd.unique(labels).tolist(), key=str)
                       if l not in lv_fit]
    for i, l in enumerate(lv_all):
        g = labels == l
        c = GRPCOL[i % len(GRPCOL)] if l in lv_fit else "0.45"
        for ser, gg in d[g].groupby("series", sort=False):
            name = (str(ser) if str(l) in ("all", str(ser)) else f"{l} -- {ser}")
            ax.scatter(x[g][(d[g]["series"] == ser).to_numpy()], gg[col], s=200, c=[c],
                       marker=MARK[ser], edgecolors="k", linewidths=1.0,  # noqa: F821
                       zorder=3, label=f"{name} ($n$ = {len(gg)})")
        if m is not None and dums and l in lv_fit:
            k = lv_fit.index(l)
            off = 0.0 if k == 0 else m["beta"][1 + k]
            gf = g & use
            # A group confined to a narrow temperature interval constrains only its
            # intercept.  Clipping its line to that interval draws a stub that reads
            # as a vertical feature, so such a group is drawn dashed across the whole
            # fitted range instead -- the model is defined there, and the dashing says
            # the slope was not measured within this group.
            panel_span = x[use].max() - x[use].min()
            spread = (gf.sum() > 1 and panel_span > 0
                      and (x[gf].max() - x[gf].min()) > 0.15 * panel_span)
            lo, hi = (x[gf].min(), x[gf].max()) if spread else (x[use].min(), x[use].max())
            xs = np.linspace(lo, hi, 100)
            yy = np.exp(m["beta"][0] + off + m["beta"][1] * xs)
            ax.plot(xs, yy, "-" if spread else "--", color=c, lw=3.0, zorder=2)
            ax.fill_between(xs, yy * np.exp(-2 * m["s"]), yy * np.exp(2 * m["s"]),
                            color=c, alpha=0.12, zorder=1)
    if m is not None and not dums:                    # one common line for the panel
        xs = np.linspace(x[use].min(), x[use].max(), 100)
        yy = np.exp(m["beta"][0] + m["beta"][1] * xs)
        ax.plot(xs, yy, "k-", lw=3.0, zorder=2)
        ax.fill_between(xs, yy * np.exp(-2 * m["s"]), yy * np.exp(2 * m["s"]),
                        color="0.5", alpha=0.16, zorder=1)
    annotate_refs(ax, d, xmode, col)                        # noqa: F821
    ax.set_xlabel(r"$T_{\rm eq}$ (K)" if xmode == "T" else r"$1000/T_{\rm eq}$ (K$^{-1}$)")
    ax.set_ylabel(ylabel)
    if ylog:
        ax.set_yscale("log")
    ax.grid(alpha=0.3, ls=":")
    place_legend_inside(ax, fontsize=FS_LEGEND)                                 # noqa: F821
    fig.tight_layout()
    save(fig, fname)

    rep = dict(figure=fname, n_plotted=len(d), n_fitted=int(use.sum()),
               levels=(lv_fit if dums else ["all"]),
               xmode=xmode, note=note, model=m, base=base,
               loo=loo_rms(y[use], [x[use]] + dums),        # noqa: F821
               p_vs_T=(f_pvalue(base, m) if dums else np.nan),  # noqa: F821
               gm=float(np.exp(np.log(y[use]).mean())),
               gsd=float(np.exp(np.log(y[use]).std(ddof=1))) if use.sum() > 1 else np.nan)
    REPORT.append(rep)
    return rep


def fmt(rep):
    """The fit of one panel as text: coefficients, offsets, scatter, diagnostics."""
    m, b = rep["model"], rep["base"]
    L = [f"--- {rep['figure']}   n plotted = {rep['n_plotted']}, n fitted = {rep['n_fitted']}"]
    if rep["note"]:
        L.append(f"    note: {rep['note']}")
    L.append(f"    geometric mean = {rep['gm']:.4g}   GSD (about the mean) = {rep['gsd']:.3f}")
    if m is None:
        return L + ["    model not identifiable"]
    v = "T" if rep["xmode"] == "T" else "1000/T"
    L.append(f"    ln y = a + b*({v})   b = {m['beta'][1]:+.5g} +/- {m['se'][1]:.3g}"
             f"   (t = {m['beta'][1]/m['se'][1]:+.2f})")
    L.append(f"    intercept '{rep['levels'][0]}'  a = {m['beta'][0]:.4f} +/- {m['se'][0]:.3g}")
    for i, l in enumerate(rep["levels"][1:], start=1):
        c, s = m["beta"][1 + i], m["se"][1 + i]
        L.append(f"      offset '{l}' = {c:+.4f} +/- {s:.3g}  -> x{np.exp(c):.3g}"
                 f"   (t = {c/s:+.2f})")
    L.append(f"    s^2 = {m['s2']:.4f}   s = {m['s']:.4f}   GSD about fit = {m['gsd']:.3f}"
             f"   R^2 = {m['R2']:.3f}   adj R^2 = {m['adjR2']:.3f}")
    L.append(f"    AICc = {m['AICc']:.1f}   LOO rms = {rep['loo']:.3f}")
    if len(rep["levels"]) > 1:
        L.append(f"    vs T only:  s^2 {b['s2']:.4f} -> {m['s2']:.4f}"
                 f"  ({100*(1-m['s2']/b['s2']):.0f} % of the variance removed),"
                 f"  F test p = {rep['p_vs_T']:.4f}")
    if rep["xmode"] == "invT":
        E = -m["beta"][1] * 1000 * KB                       # noqa: F821
        sE = abs(m["se"][1]) * 1000 * KB                    # noqa: F821
        L.append(f"    Arrhenius  y = y0 exp(-E/kT):  E_eff = {E:+.4f} +/- {sE:.4f} eV")
    return L


REPORT = []

# ---------------------------------------------------------------------------
# 1. data sources -- coverage of the (dose, temperature) plane
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 8))
STYLE = {("Neutron", "loop"): ("#1f4e79", "o"), ("Neutron", "void"): ("#1f4e79", "s"),
         ("Ion", "loop"): ("#c0392b", "^"), ("Ion", "void"): ("#c0392b", "v"),
         ("Spallation", "loop"): ("#27ae60", "D"), ("Spallation", "void"): ("#27ae60", "P")}
cov = DB[np.isfinite(DB["T_K_eq"]) & np.isfinite(DB["dose"])]          # noqa: F821
for (irr, defect), g in cov.groupby(["irradiation", "defect"], sort=False):
    c, mk = STYLE[(irr, defect)]
    ax.scatter(g["T_K_eq"], g["dose"], s=200, c=[c], marker=mk, edgecolors="k",
               linewidths=1.0, zorder=3,
               label=f"{irr} -- {'loops' if defect == 'loop' else 'cavities'}"
                     f" ($n$ = {len(g)})")
for _, r in cov.iterrows():
    ax.annotate(str(int(r["ref_no"])), (r["T_K_eq"], r["dose"]), textcoords="offset points",
                xytext=(10, 7), fontsize=PT_LABEL_FS, zorder=4)                # noqa: F821
ax.set_xlabel(r"$T_{\rm eq}$ (K)")
ax.set_ylabel(r"dose (dpa)")
ax.set_yscale("log")
ax.grid(alpha=0.3, ls=":")
place_legend_inside(ax, fontsize=FS_LEGEND)                                                        # noqa: F821
fig.tight_layout()
save(fig, "coverage_dose_temperature.png")

# ---------------------------------------------------------------------------
# 2. neutron irradiation
# ---------------------------------------------------------------------------
NL = r"loop number density $N_{\ell}$ (m$^{-3}$)"
NC = r"cavity number density $N_{c}$ (m$^{-3}$)"
DL = r"mean loop diameter $\bar{d}_{\ell}$ (nm)"
DC = r"mean cavity diameter $\bar{d}_{c}$ (nm)"

d = subset("loop", "density", NEUTRON_SERIES)
fit_and_plot(d, "density", "T", NL, "neutron_loop_density_T.png",
             labels=d["family"], ref_level="9Cr F/M", group_fit=False,
             note="alloy family colours the markers only: no grouping tested reduces the "
                  "leave-one-out error below the T-only model, so a single line is fitted")
fit_and_plot(d, "density", "invT", NL, "neutron_loop_density_invT.png",
             labels=d["family"], ref_level="9Cr F/M", group_fit=False)

d = subset("loop", "size", NEUTRON_SERIES)
fit_and_plot(d, "size", "T", DL, "neutron_loop_size_T.png", ylog=False,
             labels=d["family"], ref_level="9Cr F/M",
             note="alloy family is the selected model for this panel")

d = subset("void", "density", NEUTRON_SERIES)
fit_and_plot(d, "density", "T", NC, "neutron_cavity_density_T.png")
fit_and_plot(d, "density", "invT", NC, "neutron_cavity_density_invT.png")

d = subset("void", "size", NEUTRON_SERIES)
fit_and_plot(d, "size", "T", DC, "neutron_cavity_size_T.png", ylog=False)

# ---------------------------------------------------------------------------
# 3. ion irradiation
# ---------------------------------------------------------------------------
d = subset("loop", "density", ("Ion",))
lab, ref = dose_labels(d, 15.0)
fit_and_plot(d, "density", "T", NL, "ion_loop_density_T.png", labels=lab, ref_level=ref,
             drop_label="dose not reported", note="dose split at 15 dpa")
fit_and_plot(d, "density", "invT", NL, "ion_loop_density_invT.png",
             labels=lab, ref_level=ref, drop_label="dose not reported")

d = subset("loop", "size", ("Ion",))
lab, ref = dose_labels(d, 30.0)
fit_and_plot(d, "size", "T", DL, "ion_loop_size_T.png", ylog=False,
             labels=lab, ref_level=ref, drop_label="dose not reported",
             note="dose split at 30 dpa")

d = subset("void", "density", ("Ion",))
lab = np.where(d["cavity_kind"].eq("He bubble"), "He bubble", "void")
fit_and_plot(d, "density", "T", NC, "ion_cavity_density_T.png",
             labels=lab, ref_level="void",
             note="CONFOUNDED: in this panel the six He-bubble rows are exactly the six "
                  "rows at dose <= 30 dpa, so 'He bubble' and 'low dose' are the same split")
fit_and_plot(d, "density", "invT", NC, "ion_cavity_density_invT.png",
             labels=lab, ref_level="void")

d = subset("void", "size", ("Ion",))
lab, ref = dose_labels(d, 30.0)
fit_and_plot(d, "size", "T", DC, "ion_cavity_size_T.png", ylog=False,
             labels=lab, ref_level=ref, drop_label="dose not reported",
             note="dose split at 30 dpa; He-bubble / void is a different and weaker "
                  "split on this panel (see fit_summary)")
d2 = subset("void", "size", ("Ion",))
fit_and_plot(d2, "size", "T", DC, "ion_cavity_size_bubble_T.png", ylog=False,
             labels=np.where(d2["cavity_kind"].eq("He bubble"), "He bubble", "void"),
             ref_level="void", note="the He-bubble / void split of the same panel")

# ---------------------------------------------------------------------------
# 4. loop size distributions
# ---------------------------------------------------------------------------
SDCOL = {"300C_15dpa_area1": "#1f77b4", "300C_15dpa_area2": "#d62728",
         "300C_15dpa_area3": "#2ca02c", "350C_3dpa": "#9467bd",
         "350C_16dpa_100": "#ff7f0e", "350C_16dpa_111": "#8c564b"}


def plot_distributions(T_C, fname):
    sel = META[np.isclose(META["T_C"].astype(float), T_C)]              # noqa: F821
    fig, ax = plt.subplots(figsize=(11, 8))
    for _, m in sel.iterrows():
        k = m["dataset"]
        dd, ff = HIST[k]                                                # noqa: F821
        v = FITS[k]                                                     # noqa: F821
        c = SDCOL.get(k, "k")
        ax.step(dd, ff, where="mid", color=c, lw=2.0, alpha=0.55)
        ax.scatter(dd, ff, s=90, color=c, edgecolors="k", linewidths=0.7, zorder=3,
                   label=f"[{ref_no(m['ref'])}] {m['material']}, {m['dose']:g} dpa, "  # noqa: F821
                         f"{texify(m['population'])}")                  # noqa: F821
        dg = np.linspace(max(dd.min() * 0.5, 0.2), dd.max() * 1.15, 400)
        ax.plot(dg, 100.0 * lognormal_pdf(dg, v["mu"], v["sig"]) * v["bw"],  # noqa: F821
                color=c, lw=3.0, zorder=4)
        ax.axvline(v["mean"], color=c, ls=":", lw=1.8, alpha=0.8)
    ax.set_xlabel(r"loop diameter $d$ (nm)")
    ax.set_ylabel(r"loop number fraction per bin (\%)".replace("\\%", "%"))
    ax.set_xlim(left=0)
    ax.grid(alpha=0.3, ls=":")
    place_legend_inside(ax, fontsize=FS_LEGEND)                                             # noqa: F821
    fig.tight_layout()
    save(fig, fname)


plot_distributions(300.0, "loop_size_distribution_300C.png")
plot_distributions(350.0, "loop_size_distribution_350C.png")

fig, ax = plt.subplots(figsize=(11, 8))
order = sorted(FITS, key=lambda k: (float(META.loc[META.dataset == k, "T_C"].iloc[0]),  # noqa: F821
                                    FITS[k]["mean"]))                   # noqa: F821
for k in order:
    m = META.loc[META.dataset == k].iloc[0]                             # noqa: F821
    v = FITS[k]                                                         # noqa: F821
    dg = np.linspace(0.2, 70, 800)
    ax.plot(dg, lognormal_pdf(dg, v["mu"], v["sig"]),                   # noqa: F821
            ls="-" if float(m["T_C"]) == 300 else "--", lw=3.0, color=SDCOL.get(k, "k"),
            label=f"[{ref_no(m['ref'])}] {m['T_C']:.0f} $^\\circ$C, {m['material']}, "  # noqa: F821
                  f"{texify(m['population'])} -- $\\bar d$ = {v['mean']:.1f} nm")  # noqa: F821
ax.set_xlabel(r"loop diameter $d$ (nm)")
ax.set_ylabel(r"probability density $f(d)$ (nm$^{-1}$)")
ax.set_xlim(0, 70)
ax.grid(alpha=0.3, ls=":")
place_legend_inside(ax, fontsize=FS_LEGEND)                                                 # noqa: F821
fig.tight_layout()
save(fig, "loop_size_distribution_all.png")

# ---------------------------------------------------------------------------
# 5. <100> / 1/2<111> loop character fractions
# ---------------------------------------------------------------------------
# A fraction is bounded on [0, 100], so these panels are fitted on the logit scale
# rather than the log scale used everywhere else (see notebook Section 12).  The six
# measurements reporting exactly 0 % or 100 % are qualitative statements with an
# infinite logit: they are drawn as open symbols and excluded from the quoted fit,
# with the continuity-corrected fit that includes them shown dashed for comparison.
# Each panel is produced twice: once over both alloy classes, and once over
# EUROFER97 alone.  The ODS variants are a different material problem -- the oxide
# dispersion changes the loop population directly -- and pooling them with the
# unstrengthened steel is the single largest source of scatter in these panels, so
# the EUROFER97-only version is the cleaner calibration target.
LF_SUBSETS = (("", LOOPFRAC),                                          # noqa: F821
              ("_eurofer97", LOOPFRAC[LOOPFRAC.alloy == "EUROFER97"]))  # noqa: F821
LF_REPORT = []

for suffix, lf_all in LF_SUBSETS:
    for tag, irr in (("neutron", "Neutron"), ("ion", "Ion")):
        sub_lf = lf_all[lf_all.irradiation == irr]
        fname = f"loop_fraction_{tag}_T{suffix}.png"
        fig, fit, sens = panel_loop_fraction(sub_lf)                   # noqa: F821
        save(fig, fname)
        LF_REPORT.append((fname, len(sub_lf), fit, sens))
    save(panel_loop_fraction_by_source(lf_all),                        # noqa: F821
         f"loop_fraction_by_source{suffix}.png")

# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
lines = ["FIT SUMMARY -- group-separated analysis (neutron and ion never pooled)",
         "=" * 78, ""]
for rep in REPORT:
    lines += fmt(rep) + [""]


# --- which grouping was selected, and against what -------------------------------
def candidates(d, defect):
    """Every grouping tested on a panel.  Degenerate ones are dropped, not fitted."""
    T = cov_T(d)                                                        # noqa: F821
    dose = d["dose"].to_numpy(float)
    c = {"T only": [T],
         "T + ln(dose)": [T, cov_lndose(d)],                            # noqa: F821
         "T + ODS flag": [T, cov_ods(d)],                               # noqa: F821
         "T + alloy family": [T] + group_dummies(                       # noqa: F821
             d["family"].to_numpy(), ref="9Cr F/M")[0],
         "T + Cr wt%": [T, cov_cr(d)]}                                  # noqa: F821
    for thr in (10.0, 15.0, 20.0, 30.0):
        c[f"T + dose > {thr:g} dpa"] = [T, (dose > thr).astype(float)]
    if defect == "void":
        c["T + He-bubble flag"] = [T, cov_bubble(d)]                    # noqa: F821
    out = {}
    for k, cols in c.items():
        if len(d) - len(cols) - 1 < 1:
            continue
        if all(np.isfinite(np.asarray(v, float)).all() and np.ptp(np.asarray(v, float)) > 0
               for v in cols[1:]):
            out[k] = cols
    return out


lines += ["", "MODEL SELECTION -- every grouping tested on each panel",
          "(ranked by leave-one-out prediction error in ln y; `p_vs_T` is the nested",
          " F test against the T-only model. Rows without a reported dose, and",
          " spallation rows, are excluded so that every candidate sees the same data.)",
          "=" * 78]
for defect, dl in (("loop", "loops"), ("void", "cavities")):
    for quantity in ("density", "size"):
        for tag, irr in (("neutron", NEUTRON_SERIES), ("ion", ("Ion",))):
            d = subset(defect, quantity, irr)
            d = d[d.irradiation.ne("Spallation") & np.isfinite(d["dose"])].reset_index(drop=True)
            cand = candidates(d, defect)
            if len(cand) < 2:
                continue
            lines += ["", f"{tag} {dl} -- {quantity}   n = {len(d)}",
                      compare(d, quantity, cand).to_string(index=False)]   # noqa: F821

# --- leverage: the panels that rest on one or two points -------------------------
lines += ["", "", "LOOP CHARACTER FRACTIONS -- logit f = a + b T_eq", "=" * 78,
          "(panels with the `_eurofer97` suffix drop the EUROFER-ODS measurements.)",
          "(f is the <100> percentage of the loop population; s is on the logit",
          " scale. The quoted fit uses only measurements strictly between 0 and 100 %;",
          " `sensitivity` adds the boundary measurements under a 0.5 pp continuity",
          " correction, which is an arbitrary choice and changes the significance.)"]
for fname, n_plot, fit, sens in LF_REPORT:
    lines.append("")
    lines.append(f"--- {fname}   n plotted = {n_plot}, n fitted = {fit.get('n', 0)}")
    if "a" not in fit:
        lines.append("    fewer than 3 interior points or no temperature spread -> no fit")
    else:
        lines.append(f"    a = {fit['a']:+.4f} +/- {fit['se_a']:.4f}   "
                     f"b = {fit['b']:+.5f} +/- {fit['se_b']:.5f}  "
                     f"(t = {fit['b']/fit['se_b']:+.2f})")
        lines.append(f"    s = {fit['s']:.4f} (logit)   R^2 = {fit['R2']:.3f}   "
                     f"T(f = 50 %) = {fit['T50']:.0f} K = {fit['T50']-273.15:.0f} C")
    if "a" in sens:
        lines.append(f"    sensitivity (0 %/100 % included, eps = {LF_EPS} pp): "      # noqa: F821
                     f"n = {sens['n']}, b = {sens['b']:+.5g} +/- {sens['se_b']:.3g} "
                     f"(t = {sens['b']/sens['se_b']:+.2f}), R^2 = {sens['R2']:.3f}, "
                     f"T50 = {sens['T50']-273.15:.0f} C")

lines += ["", "within-source fits (>= 3 interior points at >= 2 temperatures):"]
for key, gg in LOOPFRAC.groupby("ref_key"):                            # noqa: F821
    gi = gg[gg.interior]
    if len(gi) >= 3 and gi["T_K_eq"].nunique() >= 2:
        r = logit_fit(gi["T_K_eq"], gi["frac"])                        # noqa: F821
        if "a" in r:
            lines.append(f"   [{int(gg.ref_no.iloc[0]):2d}] {key:22s} n = {r['n']}, "
                         f"b = {r['b']:+.5f} +/- {r['se_b']:.5f} "
                         f"(t = {r['b']/r['se_b']:+.2f}), R^2 = {r['R2']:.3f}, "
                         f"T50 = {r['T50']-273.15:.0f} C")

lines += ["", "", "LEVERAGE CHECKS", "=" * 78]
for quantity, ylab in (("density", "cavity number density"), ("size", "mean cavity diameter")):
    d = subset("void", quantity, NEUTRON_SERIES)
    d = d[d.irradiation.ne("Spallation")]
    col = QCOL[quantity]                                                # noqa: F821
    full = loglin_fit(d["T_K_eq"], d[col])                              # noqa: F821
    fm = d[d.family.ne("pure Fe")]
    cut = loglin_fit(fm["T_K_eq"], fm[col])                             # noqa: F821
    lines.append(f"neutron {ylab}: all n = {full['n']}, b = {full.get('b', np.nan):+.5g}, "
                 f"R^2 = {full.get('R2', np.nan):.3f}  |  F/M steels only n = {cut['n']}, "
                 f"b = {cut.get('b', np.nan):+.5g}, R^2 = {cut.get('R2', np.nan):.3f}, "
                 f"T range {fm['T_K_eq'].min():.0f}-{fm['T_K_eq'].max():.0f} K")

d = subset("void", "density", ("Ion",))
lines.append("ion cavity density: He-bubble rows = "
             f"{sorted(d.loc[d.cavity_kind.eq('He bubble'), 'dose'].round(1).tolist())} dpa; "
             "void rows = "
             f"{sorted(d.loc[~d.cavity_kind.eq('He bubble'), 'dose'].round(1).tolist())} dpa "
             "-> the two partitions are identical, the effect cannot be attributed")

lines += ["", "LOG-NORMAL LOOP SIZE-DISTRIBUTION FITS", "=" * 78]
for _, m in META.iterrows():                                            # noqa: F821
    v = FITS[m["dataset"]]                                              # noqa: F821
    lines.append(f"{m['dataset']:>18s} | {m['material']:>12s} | {m['T_C']:.0f} C | "
                 f"{m['dose']:g} dpa | {m['population']} | mu = {v['mu']:.4f}  "
                 f"sigma = {v['sig']:.4f} | mean = {v['mean']:.2f} nm  "
                 f"median = {v['median']:.2f}  mode = {v['mode']:.2f}  "
                 f"SD = {v['sd_ln']:.2f} | R^2 = {v['R2']:.3f}")

txt = "\n".join(lines)
open(os.path.join(OUTDIR, "fit_summary.txt"), "w", encoding="utf-8").write(txt + "\n")
print(txt)
print(f"\nwrote {len([f for f in os.listdir(OUTDIR) if f.endswith('.png')])} figures "
      f"and fit_summary.txt to {OUTDIR}")
