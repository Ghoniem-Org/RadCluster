#!/usr/bin/env python3
"""Dose-evolution figures for a RadCluster_2_1 run, with EUROFER97 data overlaid.

Regenerates `plots/density_vs_dose.png` and `plots/size_vs_dose.png` (plus the
`_annotated` variants) for ANY run directory under `output/`.

WHY THIS FILE EXISTS.  The reference figures for run
`20260904_212202_..._I80000V20000_im5vm5` were made by a script that was never
saved next to them -- `output/` is gitignored, so the only surviving relative,
`output/20260904_153005_f281fcf/make_figures.py`, is a DIFFERENT and earlier
variant: it plots the campaign row's `at_dose` trajectory out of
`digital_twin/results/B3_coal.jsonl`, not the run's own solution, and it draws
the marker-on-curve "annotated" style with no experimental bands.  The notebook
never produced these figures at all.  So the headline figures of the reference
case could not be regenerated from anything in version control.  This script
closes that gap and is committed.

The six series are derived from the run's own `plots/plot_data.pkl`, using the
same expressions `digital_twin/run_ensemble.py` uses for the campaign
observables, so a figure and a ledger row for the same vector agree by
construction:

    d_111  = 2*sqrt(mean_n_111 * Omega / (pi*b_111))     [loop, 1/2<111>]
    d_100  = 2*sqrt(mean_n_100 * Omega / (pi*b_100))     [loop, <100>]
    d_cav  = 2*(3*mean_n_v*Omega/(4*pi))^(1/3)           [spherical cavity]

Verified against the reference run: reproduces its provenance addendum's
"Model at key doses" table to every printed digit.

ENCODING (revised 2026-09-06).  Experimental points are coloured by *literature
source*, not by defect population, and every source is named in the legend --
replacing the earlier filled/open "experiment, 330 C" / "experiment, 300/350 C"
pair, which spent two legend entries on temperature and none on provenance.
Defect population is carried by marker SHAPE (o = 1/2<111>, s = <100>,
^ = cavity), and the model legend entries repeat that shape, so a reader can
still tell which curve each measurement belongs against.  The irradiation
temperature of an individual point is no longer encoded; every point in the
window is 300-350 C and the model is at 330 C.

Usage:
    python codes/make_dose_figures.py <run_dir> [--annotated] [--in-place]

    <run_dir>     an output/ run directory (must contain plots/plot_data.pkl)
    --annotated   also write the marker/annotated *_annotated.png variants
    --in-place    write into <run_dir>/plots/ (default: a new timestamped dir,
                  so an existing reference run is never overwritten)
"""
import argparse, pickle, sys, datetime, subprocess, re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerTuple
from matplotlib.colors import to_rgba

MOD = Path(__file__).resolve().parent.parent
DB = MOD.parent / 'docs' / 'Database'

# Categorical identity: one fixed hue per defect population, assigned in a fixed
# order and never cycled.  These are the three slots the reference figures used;
# they are kept verbatim so regenerated figures are comparable to the originals.
C111, C100, CCAV = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#d8d7d2'

# One fixed hue per literature source.  Deliberately drawn from hues the three
# population colours do NOT occupy (no medium blue, orange or green): a source
# colour must never read as a population colour, since the two encodings share
# the same plot.  Order follows docs/Database/SOURCES.md.
SOURCE_COLOR = {
    'Dethloff 2016':  '#111111',
    'Dethloff 2018':  '#8a1c1c',
    'Klimenkov 2011': '#b8860b',
    'Klimenkov 2020': '#7b2d8e',
    'Weiß 2012':      '#d81b60',
    'Chauhan 2021':   '#6b4423',
    'Coppola 2019':   '#4a6572',
}

LW = 4.5          # model curve weight
MS = 260          # experimental marker area (pt^2); ~16 pt across

# Axis labels and tick labels carry the figure at a glance, so they are the
# largest text; the legend stays at 16 pt (the project floor) because it carries
# up to ten entries -- three model series plus one per source.
LABEL_PT, TICK_PT, LEGEND_PT = 28, 24, 16

plt.rcParams.update({
    'font.size': TICK_PT, 'axes.labelsize': LABEL_PT, 'axes.titlesize': LABEL_PT,
    'xtick.labelsize': TICK_PT, 'ytick.labelsize': TICK_PT,
    'legend.fontsize': LEGEND_PT,
    'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': INK2, 'ytick.color': INK2, 'figure.facecolor': 'white',
    'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'axes.linewidth': 1.2, 'lines.linewidth': LW})

XMIN = 1e-3       # both figures start here
NMIN = 1e19       # density floor
DMAX = 18.0       # size figure: drop measurements above this diameter (nm)


def load_model(run_dir: Path):
    """Six observable series + dose axis, from the run's own solution."""
    pkl = run_dir / 'plots' / 'plot_data.pkl'
    if not pkl.exists():
        sys.exit(f'no plots/plot_data.pkl in {run_dir} -- re-run the notebook '
                 f'for this case, or point at a run directory that has one')
    # The pickle holds live InputData / RateEquations objects, so the package
    # must be importable before it can be loaded.
    sys.path[:0] = [str(MOD.parent), str(MOD)]
    d = pickle.load(open(pkl, 'rb'))
    R, E = d['results'], d['input_data'].energetics
    Om = float(R['Omega'])
    b111, b100 = float(E['b_111']) * 1e-9, float(E['b_100']) * 1e-9
    n111 = np.maximum(np.asarray(R['mean_n_111'], float), 0.0)
    n100 = np.maximum(np.asarray(R['mean_n_100'], float), 0.0)
    nv = np.maximum(np.asarray(R['mean_n_v'], float), 0.0)
    return np.asarray(R['dose'], float), {
        'N_loops_111': np.asarray(R['N_loops_111'], float),
        'N_loops_100': np.asarray(R['N_loops_100'], float),
        'N_voids':     np.asarray(R['N_voids'], float),
        'd_111_nm':    2 * np.sqrt(n111 * Om / (np.pi * b111)) * 1e9,
        'd_100_nm':    2 * np.sqrt(n100 * Om / (np.pi * b100)) * 1e9,
        'd_cavity_nm': 2 * (3 * nv * Om / (4 * np.pi)) ** (1 / 3) * 1e9,
    }


def short_source(title):
    """'Dethloff et al. - 2016 - Microstructural ...' -> 'Dethloff 2016'."""
    if not isinstance(title, str) or not title.strip():
        return None
    first = title.split(' - ')[0].replace(' et al.', '').strip()
    year = re.search(r'\b(19|20)\d{2}\b', title)
    return f'{first} {year.group(0)}' if year else first


def attribute(df):
    """Add a short `Source` per row of InterstitialLoop.xlsx.

    The `Paper` column uses MERGED CELLS: one title heads a block of rows and
    pandas leaves NaN on every row but the first, so attribution is a forward
    fill.  The fill is fenced by irradiation type: the sheet has one place --
    Excel row 40 -- where a neutron row sits immediately below an ion block, and
    an unfenced fill would hand it to Boulanger 2009 (an ion study).  A finer
    fence would be wrong: the Dethloff 2018 merge legitimately spans an HFR ->
    BOR-60 facility change.  Fencing leaves row 40 unattributed, exactly what
    SOURCES.md caveat 1 records; it is then assigned to Chauhan 2021 by the
    inference stated there (its 1/2<111> partner, same 6.2 nm, is Chauhan), and
    the assumption is printed at run time.  Fix the merge range to retire this.
    """
    src = df['Paper'].map(short_source)
    fence = df['Type of Irradiation'].astype(str).str.strip().str.lower()
    df = df.assign(Source=src.groupby((fence != fence.shift()).cumsum()).ffill())
    orphan = df['Source'].isna()
    if orphan.any():
        print(f'note: {int(orphan.sum())} loop row(s) sit outside every merged Paper '
              f'block (SOURCES.md caveat 1); assigned to Chauhan 2021 by inference')
        df.loc[orphan, 'Source'] = 'Chauhan 2021'
    return df


def load_experiment():
    """EUROFER97, neutron, 300-350 C.  Provenance: docs/Database/SOURCES.md."""
    def euro(df):
        m = (df['Material'].astype(str).str.contains('urofer97', na=False)
             & df['Type of Irradiation'].astype(str).str.strip().str.lower().eq('neutron')
             & df['Irradiation Temp [C]'].between(300, 350))
        return df[m].copy()
    L = euro(attribute(pd.read_excel(DB / 'InterstitialLoop.xlsx')))
    V = euro(pd.read_excel(DB / 'Void.xlsx'))
    V['Source'] = V['Paper'].map(short_source)
    # A blank Loop Type is a TOTAL loop density; it belongs to neither
    # character-resolved series and would double-count if left in.
    lt = L['Loop Type'].astype(str)
    return ((L[lt.str.contains('111', na=False)], C111, 'o',
             'N_loops_111', 'd_111_nm', r'$\frac{1}{2}\langle111\rangle$ loops'),
            (L[lt.str.contains('100', na=False)], C100, 's',
             'N_loops_100', 'd_100_nm', r'$\langle100\rangle$ loops'),
            (V, CCAV, '^', 'N_voids', 'd_cavity_nm', 'cavities'))


def draw(D, M, SER, which, annotated):
    ycol = 'N' if which == 'N' else 'd'
    ylabel = (r'Number density (m$^{-3}$)' if ycol == 'N'
              else 'Mean diameter (nm)')
    datacol = 'Density [m^-3]' if ycol == 'N' else 'Diameter [nm]'
    fig, ax = plt.subplots(figsize=(11, 8))
    seen = []                       # sources actually drawn in this figure
    shapes = {}                     # source -> marker shapes it contributes

    for df, c, mk, nk, dk, _ in SER:
        key = nk if ycol == 'N' else dk
        # Model: a smooth 2 px line.  Markers only in the annotated variant,
        # where they mark the actual computed doses.
        if annotated:
            ax.plot(D, M[key], '-', color=c, marker=mk, ms=7,
                    mec='white', mew=1.0, zorder=4)
        else:
            ax.plot(D, M[key], '-', color=c, zorder=4)

        s_all = df.dropna(subset=['Dose [dpa]', datacol])
        if ycol == 'd':
            # Diameters above DMAX are a handful of 20-48 nm outliers (the
            # 350 C Klimenkov rows).  On the linear size axis they alone set the
            # scale and press every model curve into the bottom sixth of the
            # panel, so the figure is cut at DMAX.
            s_all = s_all[s_all[datacol] <= DMAX]
        # Experimental band: the rectangle spanning this population's measured
        # dose range x its measured value range.  Experimental coverage is
        # narrow and sits at the right-hand edge, so the band is what stops a
        # reader reading agreement into the unconstrained low-dose decades.
        # Band and curve share one hue -- the band IS that population's data --
        # so it carries a saturated edge in the curve colour rather than a wash
        # that greys out where two bands overlap.  The edge is DASHED: a solid
        # one in the curve's own colour reads as a second model curve.
        if not annotated and len(s_all) > 1:
            x0, x1 = s_all['Dose [dpa]'].min(), s_all['Dose [dpa]'].max()
            y0, y1 = s_all[datacol].min(), s_all[datacol].max()
            if x1 > x0 and y1 > y0:
                ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                                       facecolor=to_rgba(c, 0.10),
                                       edgecolor=to_rgba(c, 0.8), lw=2.0,
                                       ls=(0, (6, 4)), zorder=1))
        # Colour = literature source; shape = defect population.
        for source, sub in s_all.groupby('Source'):
            ax.scatter(sub['Dose [dpa]'], sub[datacol], s=MS, marker=mk,
                       facecolors=SOURCE_COLOR.get(source, INK2),
                       edgecolors='white', linewidths=1.6, zorder=5)
            if source not in seen:
                seen.append(source)
            shapes.setdefault(source, [])
            if mk not in shapes[source]:
                shapes[source].append(mk)

    ax.set_xscale('log')
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel(ylabel)
    ax.set_xlim(left=XMIN)
    if ycol == 'N':
        ax.set_yscale('log')
        # Density floor at 1e19 m^-3.  The <100> population is still nucleating
        # at the left edge and its density there runs orders below the lowest
        # measurement; letting the axis follow it squashes the 1e20-1e22 region
        # where the model meets the data into a sliver.
        ax.set_ylim(bottom=NMIN)
    else:
        ax.set_ylim(0, DMAX)
    ax.grid(True, which='major', color=GRID, lw=0.8, zorder=0)
    ax.grid(True, which='minor', color=GRID, lw=0.4, alpha=0.55, zorder=0)
    ax.set_axisbelow(True)
    # Full box: all four spines, with ticks on all four sides.
    for sp in ax.spines.values():
        sp.set_visible(True)
    ax.tick_params(which='both', top=True, right=True, direction='out')

    # Two legends.  The model series carry both their colour and their marker
    # shape, so the shape key that identifies the experimental points is stated
    # once and the source legend needs only colour.  The source legend sits
    # BELOW the axes: six or seven entries placed inside would land on the
    # 10-32 dpa column where every measurement is, and the only clear interior
    # region is too small to hold them without clipping a curve.
    mh = [plt.Line2D([], [], color=c, lw=LW, marker=mk, ms=13, mec='white')
          for _, c, mk, _, _, _ in SER]
    ax.legend(mh, [s[5] for s in SER], frameon=not annotated,
              loc='best', labelcolor=INK,
              handletextpad=0.8, borderaxespad=0.8)
    order = ([s for s in SOURCE_COLOR if s in seen]
             + [s for s in seen if s not in SOURCE_COLOR])
    # A source's swatch shows the actual marker shapes it contributes -- one
    # marker if it measured a single population, three if it measured all of
    # them -- so no shape appears in the legend that is not on the plot, and
    # none on the plot that is missing from the legend.
    sh = [tuple(plt.Line2D([], [], lw=0, marker=mk, ms=14, mec='white', mew=1.6,
                           mfc=SOURCE_COLOR.get(s, INK2))
                for mk in shapes[s]) for s in order]
    nmax = max(len(h) for h in sh)
    if annotated:
        ax.annotate('line + small markers = model (markers at computed doses);\n'
                    'large markers = experiment',
                    xy=(0.015, 0.985), xycoords='axes fraction', va='top',
                    fontsize=13, color=INK2)
    fig.tight_layout()
    # Placed after tight_layout, in figure coordinates below the axes; the
    # savefig bbox_inches='tight' grows the canvas to take it in.
    fig.legend(sh, order, ncol=4, loc='upper center',
               bbox_to_anchor=(0.5, 0.04), frameon=False, labelcolor=INK,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.35)},
               handlelength=1.0 * nmax, handletextpad=0.6, columnspacing=1.6)
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir', type=Path)
    ap.add_argument('--annotated', action='store_true',
                    help='also write the marker/annotated *_annotated.png variants')
    ap.add_argument('--in-place', action='store_true',
                    help='write into <run_dir>/plots/ instead of a new run dir')
    a = ap.parse_args()

    run = a.run_dir.resolve()
    D, M = load_model(run)
    SER = load_experiment()
    n_pts = {s[5]: len(s[0].dropna(subset=['Dose [dpa]'])) for s in SER}
    print(f'model: {len(D)} dose points, {D[D > 0].min():.2e} - {D.max():.4g} dpa')
    print('experiment: ' + ', '.join(f'{v} {k}' for k, v in n_pts.items()))

    if a.in_place:
        out = run / 'plots'
    else:
        sha = subprocess.run(['git', '-C', str(MOD), 'rev-parse', '--short', 'HEAD'],
                             capture_output=True, text=True).stdout.strip()
        out = MOD / 'output' / f"{datetime.datetime.now():%Y%m%d_%H%M%S}_{sha}_figs" / 'plots'
    out.mkdir(parents=True, exist_ok=True)

    jobs = [('density_vs_dose.png', 'N', False), ('size_vs_dose.png', 'd', False)]
    if a.annotated:
        jobs += [('density_vs_dose_annotated.png', 'N', True),
                 ('size_vs_dose_annotated.png', 'd', True)]
    for name, which, ann in jobs:
        fig = draw(D, M, SER, which, ann)
        fig.savefig(out / name, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print('wrote', out / name)


if __name__ == '__main__':
    main()
