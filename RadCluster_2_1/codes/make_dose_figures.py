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
the titled/marker "annotated" style with no experimental bands.  The notebook
never produced these figures at all.  So the headline figures of the reference
case could not be regenerated from anything in version control.  This script
closes that gap and is committed.

The six series are derived from the run's own `plots/plot_data.pkl`, using the
same expressions `digital_twin/run_ensemble.py` uses for the campaign
observables, so a figure and a ledger row for the same vector agree by
construction:

    d_111  = 2*sqrt(mean_n_111 * Omega / (pi*b_111))     [loop, ½<111>]
    d_100  = 2*sqrt(mean_n_100 * Omega / (pi*b_100))     [loop, <100>]
    d_cav  = 2*(3*mean_n_v*Omega/(4*pi))^(1/3)           [spherical cavity]

Verified against the reference run: reproduces its provenance addendum's
"Model at key doses" table to every printed digit.

Usage:
    python codes/make_dose_figures.py <run_dir> [--annotated] [--in-place]

    <run_dir>     an output/ run directory (must contain plots/plot_data.pkl)
    --annotated   also write the titled/marker *_annotated.png variants
    --in-place    write into <run_dir>/plots/ (default: a new timestamped dir,
                  so an existing reference run is never overwritten)
"""
import argparse, pickle, sys, datetime, subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, LogLocator

MOD = Path(__file__).resolve().parent.parent
DB = MOD.parent / 'docs' / 'Database'

# Categorical identity: one fixed hue per defect population, assigned in a fixed
# order and never cycled.  These are the three slots the reference figures used;
# they are kept verbatim so regenerated figures are comparable to the originals.
C111, C100, CCAV = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#d8d7d2'

# The reference figures were drawn at base 20 pt with a 12 pt legend.  The legend
# is raised to 16 pt here: 16 pt is the project floor for plot text, and 12 pt
# legend entries are the one place the originals fell below it.  Everything else
# matches.  Set LEGEND_PT = 12 to reproduce the originals exactly.
BASE_PT, LEGEND_PT = 20, 16

plt.rcParams.update({
    'font.size': BASE_PT, 'axes.labelsize': BASE_PT, 'axes.titlesize': BASE_PT,
    'xtick.labelsize': BASE_PT, 'ytick.labelsize': BASE_PT,
    'legend.fontsize': LEGEND_PT,
    'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'text.color': INK,
    'xtick.color': INK2, 'ytick.color': INK2, 'figure.facecolor': 'white',
    'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'axes.linewidth': 1.0, 'lines.linewidth': 2.0})


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


def load_experiment():
    """EUROFER97, neutron, 300-350 C.  Provenance: docs/Database/SOURCES.md."""
    def euro(df):
        m = (df['Material'].astype(str).str.contains('urofer97', na=False)
             & df['Type of Irradiation'].astype(str).str.strip().str.lower().eq('neutron')
             & df['Irradiation Temp [C]'].between(300, 350))
        return df[m].copy()
    L = euro(pd.read_excel(DB / 'InterstitialLoop.xlsx'))
    V = euro(pd.read_excel(DB / 'Void.xlsx'))
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
        # Experimental band: the rectangle spanning this population's measured
        # dose range x its measured value range.  Experimental coverage is
        # narrow and sits at the right-hand edge, so the band is what stops a
        # reader reading agreement into the unconstrained low-dose decades.
        if not annotated and len(s_all) > 1:
            x0, x1 = s_all['Dose [dpa]'].min(), s_all['Dose [dpa]'].max()
            y0, y1 = s_all[datacol].min(), s_all[datacol].max()
            if x1 > x0 and y1 > y0:
                ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                                       facecolor=c, alpha=0.15, edgecolor='none',
                                       zorder=1))
        # Filled = measured at 330 C (the modelled condition); open = 300/350 C.
        at = df['Irradiation Temp [C]'] == 330
        for sub, filled in ((df[at], True), (df[~at], False)):
            s = sub.dropna(subset=['Dose [dpa]', datacol])
            if len(s):
                ax.scatter(s['Dose [dpa]'], s[datacol], s=130, marker=mk,
                           facecolors=c if filled else 'none', edgecolors=c,
                           linewidths=2.0, zorder=5)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Dose (dpa)')
    ax.set_ylabel(ylabel)
    ax.set_xlim(left=1e-6)
    if ycol == 'N':
        # Clip the density floor at 1e16 m^-3.  The <100> population is still
        # nucleating below ~3e-5 dpa and its density there runs down to ~1e13;
        # letting the axis follow it costs four decades of vertical space and
        # squashes the entire 1e20-1e22 region where the model meets the data
        # into a sliver.  Nothing is hidden that a reader can act on -- a
        # density three orders below the lowest measurement is not a prediction
        # anyone compares against.
        ax.set_ylim(bottom=1e16)
    if annotated:
        ax.set_title('Defect %s vs dose - EUROFER97, 330 °C'
                     % ('number density' if ycol == 'N' else 'size'), pad=14)
    if ycol == 'd':
        ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 3, 5), numticks=20))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:g}'))
        ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.grid(True, which='major', color=GRID, lw=0.8, zorder=0)
    ax.grid(True, which='minor', color=GRID, lw=0.4, alpha=0.55, zorder=0)
    ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

    # Identity is never colour-alone: every series is named in the legend, and
    # the three populations additionally differ in marker shape.
    h = [plt.Line2D([], [], color=c, lw=2, marker=(mk if annotated else None),
                    ms=8, mec='white') for _, c, mk, _, _, _ in SER]
    h += [plt.Line2D([], [], lw=0, marker='o', ms=11, mfc=INK2, mec=INK2),
          plt.Line2D([], [], lw=0, marker='o', ms=11, mfc='none', mec=INK2, mew=2)]
    lab = [s[5] for s in SER] + ['experiment, 330 °C', 'experiment, 300/350 °C']
    ax.legend(h, lab, frameon=not annotated, loc='best', labelcolor=INK,
              handletextpad=0.8, borderaxespad=0.8)
    if annotated:
        ax.annotate('line + small markers = model (markers at computed doses);\n'
                    'large markers = experiment',
                    xy=(0.015, 0.985), xycoords='axes fraction', va='top',
                    fontsize=13, color=INK2)
    fig.tight_layout()
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir', type=Path)
    ap.add_argument('--annotated', action='store_true',
                    help='also write the titled/marker *_annotated.png variants')
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
