"""
compare_window_modes.py
=======================
Run and compare the two solver modes for the I=V=10000 full_CD system:

  Mode 0 — full_system     : reference (20 006 ODEs, no window)
  Mode 4 — active_window   : two independent sliding windows + OpenMP-parallel
                              RHS (auto-serial when threads=1)

Physical parameters are taken verbatim from the notebook cell:
  T=573 K, G=1e-6 dpa/s, I=V=10000, i_mobile=v_mobile=1,
  i_discrete=I, v_discrete=V, I_bin=V_bin=0 (full per-size),
  he_kinetics=quasi_steady_state, fission cascade.

Window design (derived from I=V=1000 reference run):
  SIA front saturates at n≈62 at t=1e4 s  →  window_w0_i = 100
  VAC front reaches m≈1000 at t=1e4 s     →  window_w0_v = 500
  (expand thresholds / pads tuned to match these observations)

Usage:
  cd EuroferMicrostructure/RadCluster_1_0/codes
  python compare_window_modes.py [--ref-dir <path>] [--skip-mode0]
"""

import sys
import os
import io
import time
import argparse
import importlib
from pathlib import Path

# ── path setup ────────────────────────────────────────────────────────────────
CODES_DIR   = Path(__file__).resolve().parent.parent
MODULE_ROOT = CODES_DIR.parent
REPO_ROOT   = MODULE_ROOT.parent
for p in [str(REPO_ROOT), str(MODULE_ROOT)]:
    if p not in sys.path:
        sys.path.insert(0, p)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Reload py_utils in case stale .pyc files exist
import RadCluster_1_0.py_utils.defect_production  as _dp
import RadCluster_1_0.py_utils.binding_energies   as _be
import RadCluster_1_0.py_utils.bin_moment_rates   as _bmr
import RadCluster_1_0.py_utils.input_data         as _inp
import RadCluster_1_0.py_utils.reaction_rates     as _rr
import RadCluster_1_0.py_utils.rate_equations     as _re
import RadCluster_1_0.py_utils.cpp_bridge         as _cb
import RadCluster_1_0.py_utils.post_process       as _pp
import RadCluster_1_0.py_utils.simulation         as _sim
for _m in [_dp, _be, _bmr, _inp, _rr, _re, _cb, _pp, _sim]:
    importlib.reload(_m)
from RadCluster_1_0.py_utils.simulation import RadClusterSimulation


# ══════════════════════════════════════════════════════════════════════════════
# 1.  Shared physical parameters (mirror notebook exactly)
# ══════════════════════════════════════════════════════════════════════════════

I          = int(1e4)
V          = int(1e4)
i_mobile   = 1
v_mobile   = 1
i_discrete = I
v_discrete = V
I_bin      = 0
V_bin      = 0
C_FLOOR    = 1e-25
HE_KINETICS = 'quasi_steady_state'

PARAM_OVERRIDES = {
    'eta':       0.3,
    'f_cl_i':    0.58,
    'f_cl_v':    0.45,
    'E_m_1D':    0.34,
    'i_mobile':  i_mobile,
    'L_hat':     71.8,
    'c_C':       1.94e-4,
    'E_b_C_SIA': 0.45,
    'rho_d':     1e14,
    'Z_i':       1.1,
    'Z_ii':      1.1,
    'shape_function': 'linear',
    'i_discrete': i_discrete,
    'v_discrete': v_discrete,
    'I_bin':     I_bin,
    'V_bin':     V_bin,
}

BASE_SOLVER_CFG = {
    't_span':   (1e-6, 1e4),
    'n_points': 200,
    'log_time': True,
    'rtol':     1e-6,
    'atol':     1e-25,
}

# ── Window parameters (physics-calibrated) ───────────────────────────────────
# From I=V=1000 analysis:
#   SIA: 99.9% content at n≤24; physical front at n≈62 at t=1e4 s.
#   VAC: void front grows ~linearly; extrapolates to m≈1000 at t=1e4 s for V=10000.
WINDOW_METHOD = {
    'linsol':            'gmres',
    'window_gmres_maxl': 40,
    'window_prec':       1,
    # SIA window
    'window_w0_i':       100,     # initial SIA window: sizes 1..100
    'window_C_expand':   1e-22,   # expand when c_i[x_hi] > 1e-22
    'window_expand_pad': 50,      # grow SIA window by 50 at a time
    # VAC window
    'window_w0_v':       500,     # initial VAC window: sizes 1..500
    'window_C_expand_v': 1e-22,   # expand when c_v[x_hi] > 1e-22
    'window_expand_pad_v': 200,   # grow VAC window by 200 at a time
    # Misc
    'window_check_every': 1,
    'window_N_thresh':    500,
}

MODE4_CFG = {**BASE_SOLVER_CFG, 'solver_method': WINDOW_METHOD}
MODE0_CFG = {**BASE_SOLVER_CFG, 'solver_method': {
    'linsol': 'gmres',
    'window_gmres_maxl': 40, 'window_prec': 1,
}}


# ══════════════════════════════════════════════════════════════════════════════
# 2.  Helper: build a fresh sim with the shared parameters
# ══════════════════════════════════════════════════════════════════════════════

PHYSICS_OPTION = 'full_CD_fission'   # overridden by --physics CLI arg


def build_sim(solver_mode, physics_option=None):
    """Construct an RadClusterSimulation with the reference parameters."""
    if physics_option is None:
        physics_option = PHYSICS_OPTION
    saved = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=I, V=V,
            solver_mode=solver_mode,
            physics_option=physics_option,
            C_floor=C_FLOOR,
            he_kinetics=HE_KINETICS,
            i_mobile=i_mobile,
            v_mobile=v_mobile,
        )
    finally:
        sys.stdout, sys.stderr = saved

    saved2 = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        inp = sim.input_data
        for key, val in PARAM_OVERRIDES.items():
            placed = False
            for d in [inp.production_fission, inp.production_fusion,
                      inp.diffusion, inp.reactions,
                      inp.energetics, inp.dissociation]:
                if key in d:
                    d[key] = val
                    placed = True
            if not placed:
                inp.reactions[key] = val
        if 'i_mobile' in PARAM_OVERRIDES:
            inp.diffusion['i_mobile'] = int(PARAM_OVERRIDES['i_mobile'])
            inp.reactions['i_mobile'] = int(PARAM_OVERRIDES['i_mobile'])
        inp._calculate_derived()
        sim.rebuild_rates()
    finally:
        sys.stdout, sys.stderr = saved2

    return sim


# ══════════════════════════════════════════════════════════════════════════════
# 3.  Run a mode, return (results_dict, wall_clock_s, output_dir)
# ══════════════════════════════════════════════════════════════════════════════

def run_mode(label, solver_mode, solver_cfg):
    print(f'\n{"="*70}')
    print(f'  Running {label}  ({solver_mode} / {PHYSICS_OPTION})')
    print(f'{"="*70}')
    sim = build_sim(solver_mode)
    t0  = time.perf_counter()
    results = sim.run_adaptive(
        solver_config=solver_cfg,
        save_output=True,
        progress_callback=None,
        boundary_threshold=0.1,
        max_doublings=0,
        points_per_segment=10,
    )
    wall = time.perf_counter() - t0
    if results is None:
        print(f'  {label}: FAILED')
        return None, wall, None

    # Locate the just-written output directory (newest)
    out_root = MODULE_ROOT / 'output'
    dirs = sorted(out_root.glob(f'*_{solver_mode}_*'), key=lambda p: p.name)
    out_dir = dirs[-1] if dirs else None

    t_arr = results['t']
    print(f'  {label}: {len(t_arr)} pts, '
          f't=[{t_arr[0]:.2e},{t_arr[-1]:.2e}] s, '
          f'wall={wall:.1f} s')
    print(f'  Swelling(final)={results["swelling"][-1]*100:.5f} %  '
          f'N_loops={results["N_loops"][-1]:.3e} m^-3  '
          f'N_voids={results["N_voids"][-1]:.3e} m^-3')
    print(f'  delta_FP={results["delta_FP"][-1]:.2e}  '
          f'delta_He={results["delta_He"][-1]:.2e}')
    return results, wall, out_dir


# ══════════════════════════════════════════════════════════════════════════════
# 4.  Locate / load existing runs from saved .npy files
# ══════════════════════════════════════════════════════════════════════════════

def _load_run_from_dir(d, solver_mode, label):
    """
    Load a saved run from directory *d*.  Returns (results_dict, dir_name) or
    (None, None) if the directory doesn't contain valid results.
    """
    d = Path(d)
    t_path = d / 'results_t.npy'
    y_path = d / 'results_y.npy'
    if not (t_path.exists() and y_path.exists()):
        return None, None
    t = np.load(t_path)
    y = np.load(y_path, mmap_mode='r')
    print(f'  Found {label}: {d.name}  ({y.shape[1]} pts,  N_eq={y.shape[0]})')
    sim = build_sim(solver_mode)
    re  = sim.rate_equations
    from RadCluster_1_0.py_utils.post_process import calculate_derived_quantities
    saved = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        res = calculate_derived_quantities(t, y, sim.input_data, re)
    finally:
        sys.stdout, sys.stderr = saved
    res['y'] = y
    res['_source_dir'] = d.name
    return res, d.name


def find_reference_run(ref_dir=None):
    """
    Return the results dict for mode-0.  If ref_dir is given, load from there.
    Otherwise scan the output directory for the most recent full_system run
    with N_eq == 20006 (I=V=10000, QSS He) matching PHYSICS_OPTION (or
    bin_moment fallback since they are mathematically identical with I_bin=0).
    Legacy `cpp_full` directory names are also matched.
    """
    out_root = MODULE_ROOT / 'output'

    # Prefer exact physics match, then accept bin_moment (mathematically identical).
    # Both new (`full_system`) and legacy (`cpp_full`) directory prefixes are scanned.
    glob_patterns = [
        f'*_full_system_{PHYSICS_OPTION}_*',
        f'*_cpp_full_{PHYSICS_OPTION}_*',
        '*_full_system_bin_moment_CD_fission_*',
        '*_cpp_full_bin_moment_CD_fission_*',
        '*_full_system_full_CD_fission_*',
        '*_cpp_full_full_CD_fission_*',
    ]
    candidates = []
    for pat in glob_patterns:
        candidates += sorted(out_root.glob(pat), key=lambda p: p.name, reverse=True)
    if ref_dir:
        candidates = [Path(ref_dir)] + candidates

    seen = set()
    for d in candidates:
        if d in seen:
            continue
        seen.add(d)
        y_path = d / 'results_y.npy'
        if not y_path.exists():
            continue
        y = np.load(y_path, mmap_mode='r')
        if y.shape[0] == 20006:      # I=V=10000, QSS-He
            return _load_run_from_dir(d, 'full_system', 'mode-0 reference')
    return None, None


def find_window_run(load_dir=None):
    """
    Return the results dict for active_window.  If load_dir is given, load
    from there; otherwise find the most recent matching output directory for
    PHYSICS_OPTION.  Legacy `sliding_OpenMP` directory names are also matched.
    """
    out_root = MODULE_ROOT / 'output'
    label = 'Mode IV'

    if load_dir:
        return _load_run_from_dir(Path(load_dir), 'active_window', label)

    glob_patterns = [
        f'*_active_window_{PHYSICS_OPTION}_*',
        f'*_sliding_OpenMP_{PHYSICS_OPTION}_*',
    ]
    candidates = []
    for pat in glob_patterns:
        candidates += sorted(out_root.glob(pat), key=lambda p: p.name, reverse=True)
    for d in candidates:
        res, src = _load_run_from_dir(d, 'active_window', label)
        if res is not None:
            return res, src
    return None, None


# ══════════════════════════════════════════════════════════════════════════════
# 5.  Comparison plots
# ══════════════════════════════════════════════════════════════════════════════

COLORS  = {'mode0': '#1f77b4', 'mode4': '#2ca02c'}
LABELS  = {'mode0': 'Mode 0 — full_system (ref)',
           'mode4': 'Mode IV — active_window'}
LS      = {'mode0': '-', 'mode4': '--'}

def _relerr(ref, new):
    """Element-wise relative error |new - ref| / (|ref| + 1e-30)."""
    return np.abs(new - ref) / (np.abs(ref) + 1e-30)


def make_comparison_plots(res_dict, wall_dict, out_path):
    """
    res_dict : {'mode0': results, 'mode4': results}  (any key may be absent)
    wall_dict: {'mode0': float,   'mode4': float}
    out_path : Path to save the figure
    """
    modes = [k for k in ['mode0', 'mode4'] if k in res_dict and res_dict[k] is not None]

    fig = plt.figure(figsize=(18, 22))
    gs  = gridspec.GridSpec(4, 3, figure=fig, hspace=0.42, wspace=0.35)

    # ── (0,0) Swelling ──────────────────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    for m in modes:
        r = res_dict[m]
        ax.loglog(r['dose'], r['swelling']*100,
                  color=COLORS[m], ls=LS[m], lw=1.8, label=LABELS[m])
    ax.set_xlabel('Dose (dpa)');  ax.set_ylabel('Swelling (%)')
    ax.set_title('Void swelling');  ax.legend(fontsize=7)

    # ── (0,1) Number densities ───────────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 1])
    for m in modes:
        r = res_dict[m]
        ax.loglog(r['dose'], r['N_loops'], color=COLORS[m], ls=LS[m], lw=1.8)
        ax.loglog(r['dose'], r['N_voids'], color=COLORS[m], ls=LS[m], lw=1.0, alpha=0.5)
    # phantom lines for legend
    from matplotlib.lines import Line2D
    ax.add_line(Line2D([],[], color='k', lw=1.8, label='Loops (solid)'))
    ax.add_line(Line2D([],[], color='k', lw=1.0, alpha=0.5, label='Voids (faint)'))
    ax.set_xlabel('Dose (dpa)');  ax.set_ylabel('Number density (m⁻³)')
    ax.set_title('N_loops & N_voids');  ax.legend(fontsize=7)

    # ── (0,2) Mean cluster sizes ─────────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 2])
    for m in modes:
        r = res_dict[m]
        ax.loglog(r['dose'], r['mean_n_i'], color=COLORS[m], ls=LS[m], lw=1.8)
        ax.loglog(r['dose'], r['mean_n_v'], color=COLORS[m], ls=LS[m], lw=1.0, alpha=0.5)
    ax.set_xlabel('Dose (dpa)');  ax.set_ylabel('Mean cluster size (defects)')
    ax.set_title('Mean sizes (loops solid, voids faint)')

    # ── (1,0) Point-defect concentrations ───────────────────────────────────
    ax = fig.add_subplot(gs[1, 0])
    for m in modes:
        r = res_dict[m]
        ax.loglog(r['dose'], r['C_i1'], color=COLORS[m], ls=LS[m], lw=1.8)
        ax.loglog(r['dose'], r['C_v1'], color=COLORS[m], ls=LS[m], lw=1.0, alpha=0.5)
    ax.set_xlabel('Dose (dpa)');  ax.set_ylabel('Concentration (m⁻³)')
    ax.set_title('Monomers C_i1 (solid) & C_v1 (faint)')

    # ── (1,1) SIA and VAC total content ─────────────────────────────────────
    ax = fig.add_subplot(gs[1, 1])
    for m in modes:
        r = res_dict[m]
        S_I = r['C_SIA_tot'] * r['Omega']   # SIA content (atom fraction)
        S   = r['swelling']                  # VAC content / void swelling (atom fraction)
        ax.loglog(r['dose'], S_I, color=COLORS[m], ls=LS[m], lw=1.8, label=LABELS[m])
        ax.loglog(r['dose'], S,   color=COLORS[m], ls=LS[m], lw=1.0, alpha=0.5)
    ax.set_xlabel('Dose (dpa)');  ax.set_ylabel('Defect content (atom fraction)')
    ax.set_title('S_I SIA content (solid) & S VAC content (faint)')
    ax.legend(fontsize=7)

    # ── (1,2) Conservation diagnostics ──────────────────────────────────────
    ax = fig.add_subplot(gs[1, 2])
    for m in modes:
        r = res_dict[m]
        ax.loglog(r['dose'], np.maximum(r['delta_FP'], 1e-16),
                  color=COLORS[m], ls=LS[m], lw=1.8, label=r'$\delta_{FP}$ '+LABELS[m])
    ax.axhline(1e-6, color='gray', ls='--', lw=0.8, label='1e-6 threshold')
    ax.set_xlabel('Dose (dpa)');  ax.set_ylabel(r'$\delta_{FP}$ (relative)')
    ax.set_title('Frenkel-pair conservation');  ax.legend(fontsize=6)

    # ── (2,0) Relative errors vs mode 0 ─────────────────────────────────────
    ref = res_dict.get('mode0')
    if ref is not None and 'mode4' in res_dict and res_dict['mode4'] is not None:
        quantities = ['swelling', 'N_loops', 'N_voids', 'mean_n_i', 'mean_n_v', 'delta_FP']
        q_labels   = ['Swelling', 'N_loops', 'N_voids', r'$\bar{n}_i$', r'$\bar{n}_v$', r'$\delta_{FP}$']
        q_colors   = ['C0','C1','C2','C3','C4','C5']

        ax = fig.add_subplot(gs[2, 0])
        wr    = res_dict['mode4']
        t_ref = ref['dose']
        t_win = wr['dose']
        for qn, ql, qc in zip(quantities, q_labels, q_colors):
            try:
                ref_q = np.interp(t_win, t_ref, ref[qn])
                win_q = wr[qn]
                err   = _relerr(ref_q, win_q)
                ax.loglog(t_win, np.maximum(err, 1e-16),
                          color=qc, lw=1.4, label=ql)
            except Exception:
                pass
        ax.axhline(1e-3, color='gray', ls='--', lw=0.8, label='0.1% error')
        ax.set_xlabel('Dose (dpa)')
        ax.set_ylabel('|ε| relative error')
        ax.set_title(f'Relative error: {LABELS["mode4"]}\nvs Mode 0 reference')
        ax.legend(fontsize=6, ncol=2)

    # ── (3,0) Wall-clock time bar chart ─────────────────────────────────────
    ax = fig.add_subplot(gs[3, 0])
    labels_bar = [LABELS[m] for m in modes if m in wall_dict]
    times_bar  = [wall_dict[m]  for m in modes if m in wall_dict]
    colors_bar = [COLORS[m]     for m in modes if m in wall_dict]
    bars = ax.bar(range(len(labels_bar)), times_bar, color=colors_bar, width=0.5)
    ax.set_xticks(range(len(labels_bar)))
    ax.set_xticklabels([l.split('—')[0].strip() for l in labels_bar],
                       rotation=15, ha='right', fontsize=8)
    ax.set_ylabel('Wall-clock time (s)')
    ax.set_title('Solver timing comparison')
    for bar, t in zip(bars, times_bar):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height()*1.02,
                f'{t:.0f}s', ha='center', va='bottom', fontsize=8)
    if times_bar:
        ax.set_ylim(0, max(times_bar)*1.25)

    # ── (3,1) Speedup ratio ──────────────────────────────────────────────────
    ax = fig.add_subplot(gs[3, 1])
    if 'mode0' in wall_dict and wall_dict['mode0'] > 0 and 'mode4' in wall_dict:
        speedup = wall_dict['mode0'] / wall_dict['mode4']
        ax.bar([0], [speedup], color=COLORS['mode4'], width=0.4)
        ax.set_xticks([0])
        ax.set_xticklabels([LABELS['mode4'].split('—')[0].strip()],
                           rotation=15, ha='right', fontsize=8)
        ax.set_ylabel('Speedup vs Mode 0')
        ax.set_title('Speedup factor')
        ax.axhline(1.0, color='gray', ls='--', lw=0.8)
        ax.text(0, speedup*1.02, f'{speedup:.1f}×',
                ha='center', va='bottom', fontsize=9)
    else:
        ax.text(0.5, 0.5, 'Mode 0 timing\nnot available',
                ha='center', va='center', transform=ax.transAxes, fontsize=9)
        ax.set_title('Speedup factor')

    # ── (3,2) Annotation box ─────────────────────────────────────────────────
    ax = fig.add_subplot(gs[3, 2])
    ax.axis('off')
    lines = [
        'Run parameters',
        f'  I = V = {I:,}',
        f'  i_mobile = v_mobile = {i_mobile}',
        f'  G = 1×10⁻⁶ dpa/s',
        f'  T = 573 K',
        f'  t_span = (1e-6, 1e4) s',
        f'  rtol={BASE_SOLVER_CFG["rtol"]:.0e}  atol={BASE_SOLVER_CFG["atol"]:.0e}',
        '',
        'Window parameters (Mode IV)',
        f'  SIA: w0={WINDOW_METHOD["window_w0_i"]}  '
            f'pad={WINDOW_METHOD["window_expand_pad"]}',
        f'  VAC: w0={WINDOW_METHOD["window_w0_v"]}  '
            f'pad={WINDOW_METHOD["window_expand_pad_v"]}',
        f'  C_expand = {WINDOW_METHOD["window_C_expand"]:.0e}',
        f'  linsol = gmres  maxl={WINDOW_METHOD["window_gmres_maxl"]}',
        f'  OMP threads = auto from N_eq (override via OMP_NUM_THREADS)',
        '',
        'From I=V=1000 reference:',
        '  SIA 99.9% content: n≤24',
        '  VAC front at t=1e4s: m≈1000',
    ]
    ax.text(0.05, 0.97, '\n'.join(lines), transform=ax.transAxes,
            fontsize=7.5, va='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    fig.suptitle(
        'Sliding-window solver comparison\n'
        f'Mode 0 (full_system) vs Mode IV (active_window)\n'
        f'I = V = {I:,}   physics: {PHYSICS_OPTION}   fission cascade',
        fontsize=11, y=0.995
    )
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f'\nSaved comparison figure → {out_path}')
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# 6.  SIA & VAC size distribution comparison (final time)
# ══════════════════════════════════════════════════════════════════════════════

def make_distribution_plots(res_dict, out_path):
    modes = [k for k in ['mode0','mode4'] if k in res_dict and res_dict[k] is not None]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for m in modes:
        y = res_dict[m].get('y')
        if y is None:
            continue
        ci = y[:I, -1]
        cv = y[I:I+V, -1]
        n_arr = np.arange(1, I+1)
        m_arr = np.arange(1, V+1)

        # SIA
        mask_i = ci > 1e-30
        axes[0].loglog(n_arr[mask_i], ci[mask_i],
                       color=COLORS[m], ls=LS[m], lw=1.6, label=LABELS[m])
        # VAC
        mask_v = cv > 1e-30
        axes[1].loglog(m_arr[mask_v], cv[mask_v],
                       color=COLORS[m], ls=LS[m], lw=1.6)

    axes[0].set_xlabel('SIA cluster size n')
    axes[0].set_ylabel('Concentration (atom fraction)')
    axes[0].set_title('SIA size distribution at t=10⁴ s')
    axes[0].legend(fontsize=8)

    axes[1].set_xlabel('Vacancy cluster size m')
    axes[1].set_ylabel('Concentration (atom fraction)')
    axes[1].set_title('VAC size distribution at t=10⁴ s')

    fig.suptitle(
        f'Cluster size distributions at final time  —  I=V={I:,}',
        fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f'Saved distribution figure → {out_path}')
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# 7.  Text report
# ══════════════════════════════════════════════════════════════════════════════

def write_report(res_dict, wall_dict, report_path):
    lines = []
    sep   = '=' * 72

    lines += [sep,
              'SLIDING-WINDOW SOLVER COMPARISON REPORT',
              f'RadCluster_1_0  —  {PHYSICS_OPTION}  —  I=V=10,000',
              sep, '']

    # ── Run parameters ───────────────────────────────────────────────────────
    lines += ['PHYSICAL PARAMETERS',
              '-'*40,
              f'  Domain:          I = V = {I:,}',
              f'  Mobility:        i_mobile = v_mobile = {i_mobile}',
              f'  Dose rate:       G = 1×10⁻⁶ dpa/s',
              f'  Temperature:     T = 573 K',
              f'  t_span:          (1×10⁻⁶, 1×10⁴) s',
              f'  He treatment:    quasi_steady_state',
              f'  i_discrete:      {i_discrete}  (full per-size, no bins)',
              f'  v_discrete:      {v_discrete}',
              f'  N_eq (mode 0):   20,006',
              '',
              'SOLVER SETTINGS',
              '-'*40,
              f'  linsol:          gmres   gmres_maxl=40',
              f'  rtol:            {BASE_SOLVER_CFG["rtol"]:.0e}',
              f'  atol:            {BASE_SOLVER_CFG["atol"]:.0e}',
              '',
              'WINDOW PARAMETERS (Mode IV)',
              '-'*40,
              f'  SIA:  w0={WINDOW_METHOD["window_w0_i"]}  '
                  f'C_expand={WINDOW_METHOD["window_C_expand"]:.0e}  '
                  f'pad={WINDOW_METHOD["window_expand_pad"]}',
              f'  VAC:  w0={WINDOW_METHOD["window_w0_v"]}  '
                  f'C_expand_v={WINDOW_METHOD["window_C_expand_v"]:.0e}  '
                  f'pad_v={WINDOW_METHOD["window_expand_pad_v"]}',
              f'  OMP threads:  auto from N_eq (override via OMP_NUM_THREADS)',
              '']

    # ── Timing table ─────────────────────────────────────────────────────────
    lines += ['WALL-CLOCK TIME',
              '-'*40]
    w0 = wall_dict.get('mode0', None)
    for m, lbl in [('mode0','full_system'), ('mode4','active_window')]:
        if m in wall_dict:
            w = wall_dict[m]
            speedup = f'  ({w0/w:.1f}× speedup)' if w0 and m != 'mode0' else ''
            lines.append(f'  {lbl:<22}  {w:>8.1f} s{speedup}')
    lines.append('')

    # ── Accuracy table ───────────────────────────────────────────────────────
    ref = res_dict.get('mode0')
    lines += ['ACCURACY (relative error vs mode 0 at final time)',
              '-'*40]
    if ref is None:
        lines.append('  [mode 0 reference not available]')
    elif 'mode4' not in res_dict or res_dict['mode4'] is None:
        lines.append('  [Mode IV result not available]')
    else:
        keys = ['swelling','N_loops','N_voids','mean_n_i','mean_n_v','delta_FP']
        lines.append(f'  {"Quantity":<14}  {"Mode IV":>12}')
        lines.append('  ' + '-'*30)
        for k in keys:
            try:
                r0   = float(ref[k][-1])
                rw   = float(res_dict['mode4'][k][-1])
                err  = abs(rw - r0) / (abs(r0) + 1e-30)
                lines.append(f'  {k:<14}  {err:>12.3e}')
            except Exception:
                lines.append(f'  {k:<14}  {"N/A":>12}')
    lines.append('')

    # ── Conservation ─────────────────────────────────────────────────────────
    lines += ['CONSERVATION DIAGNOSTICS (at final time)',
              '-'*40]
    for m, lbl in [('mode0','full_system'), ('mode4','active_window')]:
        if m in res_dict and res_dict[m] is not None:
            r = res_dict[m]
            lines.append(f'  {lbl:<22}  '
                         f'delta_FP={r["delta_FP"][-1]:.3e}  '
                         f'delta_He={r["delta_He"][-1]:.3e}')
    lines.append('')

    # ── Summary ──────────────────────────────────────────────────────────────
    lines += ['SUMMARY',
              '-'*40,
              '  Window design is based on I=V=1000 analysis:',
              '    SIA: 99.9% content always in n≤24; active front at n≈62',
              '    VAC: void front grows ~linearly; reaches m≈1000 at t=1e4 s',
              '  For I=V=10000 this gives active window fractions of:',
              f'    SIA: ~{100/I*100:.1f}% of domain at saturation',
              f'    VAC: ~{1200/V*100:.1f}% of domain at saturation',
              '  Expected GMRES cost reduction: ~15× per linear solve.',
              '',
              '  Mode IV should give results within rtol of mode 0.',
              '  Any residual difference is due to boundary truncation, not',
              '  algorithmic error; increase window_expand_pad to reduce it.',
              sep]

    report_path.write_text('\n'.join(lines))
    print(f'Saved text report   → {report_path}')


# ══════════════════════════════════════════════════════════════════════════════
# 8.  Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref-dir',      help='Path to existing mode-0 output directory')
    ap.add_argument('--load-mode4',   help='Load Mode IV from this existing output directory')
    ap.add_argument('--skip-mode0',   action='store_true',
                    help='Skip running mode 0 (use existing reference only)')
    ap.add_argument('--skip-mode4',   action='store_true',
                    help='Skip running Mode IV (load existing if available)')
    ap.add_argument('--wall-mode0',   type=float, default=0.0,
                    help='Known wall-clock time for mode 0 (s), for speedup plots')
    ap.add_argument('--wall-mode4',   type=float, default=0.0,
                    help='Known wall-clock time for Mode IV (s)')
    ap.add_argument('--physics',       default='full_CD_fission',
                    choices=['full_CD_fission','full_CD_fusion',
                             'bin_moment_CD_fission','bin_moment_CD_fusion'],
                    help='Physics option for all modes (default: full_CD_fission)')
    ap.add_argument('--out-dir',      default=str(MODULE_ROOT / 'output' / 'window_comparison'),
                    help='Directory for comparison plots and report')
    args = ap.parse_args()

    # Set global physics option so build_sim / find_* use it
    global PHYSICS_OPTION
    PHYSICS_OPTION = args.physics

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    res_dict  = {}
    wall_dict = {}

    # ── Mode 0 (reference) ───────────────────────────────────────────────────
    if not args.skip_mode0:
        r0, w0, _ = run_mode('Mode 0 — full_system', 'full_system', MODE0_CFG)
        res_dict['mode0']  = r0
        wall_dict['mode0'] = w0
    else:
        print('\nSearching for existing mode-0 reference run …')
        r0, src = find_reference_run(args.ref_dir)
        if r0 is not None:
            res_dict['mode0'] = r0
            wall_dict['mode0'] = args.wall_mode0
            print(f'  Loaded from {src}')
        else:
            print('  No I=V=10000 mode-0 run found.  Relative-error plots will be skipped.')

    # ── Mode IV ──────────────────────────────────────────────────────────────
    if not args.skip_mode4:
        r4, w4, _ = run_mode('Mode IV  — active_window',  'active_window',  MODE4_CFG)
        res_dict['mode4']  = r4
        wall_dict['mode4'] = w4
    else:
        print('\nSearching for existing Mode IV run …')
        r4, src4 = find_window_run(args.load_mode4)
        if r4 is not None:
            res_dict['mode4']  = r4
            wall_dict['mode4'] = args.wall_mode4
            print(f'  Loaded from {src4}')
        else:
            print('  No Mode IV run found.')

    # ── Plots and report ─────────────────────────────────────────────────────
    make_comparison_plots(res_dict, wall_dict, out_dir / 'comparison_overview.png')
    make_distribution_plots(res_dict,          out_dir / 'size_distributions.png')
    write_report(res_dict, wall_dict,           out_dir / 'comparison_report.txt')

    print('\nAll done.  Output in:', out_dir)


if __name__ == '__main__':
    main()
