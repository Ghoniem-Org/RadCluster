"""
test_window_solver.py
=====================
Validates the Phase-I dynamic-window C++ solver against the full (no-window)
C++ solver, using the same physics as the reference run
  output/20260313_085443_26edc1d  (Nv=200, Ni=1000, t=(1e-3,1e6), CVODE BDF BAND)

The full reference run (238 s) is NOT re-run here.  Instead:

  TEST 1 – Transient accuracy (Nv=200 Ni=300 t=1e-3→1e4):
    Full solver vs Window solver (w0_i=60).
    Checks Cv1, Ci1 at all 100 output points.
    Ci20 / mean_i are checked only at EARLY times (t<10 s) when the active
    front lies well within the initial window, so the two solutions are
    analytically identical.

  TEST 2 – Reference configuration (Nv=200 Ni=1000 t=1e-3→1e6, n=1000):
    Window solver only (w0_i=100).
    Correct Cv1 steady-state (3.264e-6) determined by running the full
    solver (238 s) separately; asserted here as (2e-6, 5e-6).
    Physical sanity + speedup check.

Usage:
    cd ClusterDynamics/codes
    python test_window_solver.py
"""

import sys, time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── path setup ────────────────────────────────────────────────────────────────
try:
    _HERE = Path(__file__).resolve().parent
except NameError:
    _HERE = Path().resolve()
BASE_DIR = _HERE.parent
OUTPUT_DIR = BASE_DIR / 'output'
OUTPUT_DIR.mkdir(exist_ok=True)
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from py_utils.simulation import ClusterDynamicsSimulation
from py_utils.cpp_bridge import run_cpp_solver
from py_utils.visualization import create_run_directory

# ── helpers ───────────────────────────────────────────────────────────────────
PASS = '\033[92m✓ PASS\033[0m'
FAIL = '\033[91m✗ FAIL\033[0m'

def check(cond, msg):
    tag = PASS if cond else FAIL
    print(f"  {tag}  {msg}")
    return cond

def rel_err_max(a, b, floor=1e-30):
    """Max relative error where both values are above floor."""
    mask = (np.abs(a) > floor) & (np.abs(b) > floor)
    if not mask.any():
        return 0.0
    return float(np.max(np.abs(a[mask] - b[mask]) / np.abs(b[mask])))

def rel_err_max_early(a, b, t, t_cutoff=10.0, floor=1e-30, sig_frac=None):
    """Max relative error at time points t <= t_cutoff.

    If sig_frac is given, only compare where b > sig_frac * max(b) — avoids
    meaningless large relative errors when b is at the numerical floor.
    """
    idx = np.where(t <= t_cutoff)[0]
    if len(idx) == 0:
        return 0.0
    a_e, b_e = a[idx], b[idx]
    if sig_frac is not None:
        threshold = max(float(np.max(np.abs(b_e))), 1e-300) * sig_frac
        floor = max(floor, threshold)
    return rel_err_max(a_e, b_e, floor)

# ══════════════════════════════════════════════════════════════════════════════
# TEST 1 – Transient accuracy: full vs window (Nv=200, Ni=300)
# ══════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('TEST 1  –  Transient accuracy (Nv=200, Ni=300)')
print('  t=(1e-3, 1e4)  n_pts=100')
print('  Full solver (window_mode=0, band) vs Window solver (w0_v=100, w0_i=60)')
print('  Cv1, Ci1 compared over all time  |  Ci20/mean_i only for t < 10 s')
print('='*70)

NV1, NI1 = 200, 300
T_SPAN1   = (1e-3, 1e4)
N_PTS1    = 100

sim1 = ClusterDynamicsSimulation(Nv=NV1, Ni=NI1)

cfg_full = {
    't_span': T_SPAN1, 'n_points': N_PTS1, 'rtol': 1e-4, 'atol': 1e-20,
    'log_time': True,
    'solver_method': {
        'backend': 'cvode', 'lmm': 'bdf', 'linsol': 'band',
        'mu': NV1+NI1-1, 'ml': NV1+NI1-1,
    },
}
cfg_win = {
    't_span': T_SPAN1, 'n_points': N_PTS1, 'rtol': 1e-4, 'atol': 1e-20,
    'log_time': True,
    'solver_method': {
        'backend': 'cvode', 'lmm': 'bdf',
        'window_mode': 1,
        'window_w0_v': 100, 'window_w0_i': 60,
        'window_C_expand': 1e-18, 'window_expand_pad': 10,
        'window_check_every': 1,
    },
}

t0 = time.perf_counter()
res_full = run_cpp_solver(sim1, cfg_full)
t_full = time.perf_counter() - t0
t0 = time.perf_counter()
res_win = run_cpp_solver(sim1, cfg_win)
t_win = time.perf_counter() - t0

all_pass = True

if res_full is None or res_win is None:
    print(f"  {FAIL}  One or both solvers returned None")
    all_pass = False
else:
    t_f = res_full['time']
    t_w = res_win['time']

    all_pass &= check(len(t_f) == N_PTS1 and len(t_w) == N_PTS1,
                      f"Both produce {N_PTS1} output points")
    all_pass &= check(np.allclose(t_f, t_w, rtol=1e-6), "Time vectors match")

    Cv1_f = res_full['concentrations']['Cv1']
    Cv1_w = res_win['concentrations']['Cv1']
    err_Cv1_all   = rel_err_max(Cv1_w, Cv1_f)
    err_Cv1_early = rel_err_max_early(Cv1_w, Cv1_f, t_f, t_cutoff=10.0)
    all_pass &= check(err_Cv1_all < 0.05,
                      f"Cv1 max rel-error (all t) = {err_Cv1_all:.3e}  (tol 5%)")
    all_pass &= check(err_Cv1_early < 0.005,
                      f"Cv1 max rel-error (t≤10s) = {err_Cv1_early:.3e}  (tol 0.5%)")

    Ci1_f = res_full['concentrations']['Ci1']
    Ci1_w = res_win['concentrations']['Ci1']
    err_Ci1_all   = rel_err_max(Ci1_w, Ci1_f)
    err_Ci1_early = rel_err_max_early(Ci1_w, Ci1_f, t_f, t_cutoff=10.0)
    all_pass &= check(err_Ci1_all < 0.05,
                      f"Ci1 max rel-error (all t) = {err_Ci1_all:.3e}  (tol 5%)")
    all_pass &= check(err_Ci1_early < 0.005,
                      f"Ci1 max rel-error (t≤10s) = {err_Ci1_early:.3e}  (tol 0.5%)")

    # Ci20: at t<=10 s the window (x_hi=60) fully covers Ci20 and the front is
    # at ~x<20; the two solvers are analytically identical → tight tolerance.
    Ci20_f = res_full['concentrations'].get('Ci20')
    Ci20_w = res_win['concentrations'].get('Ci20')
    if Ci20_f is not None and Ci20_w is not None:
        # Compare only where Ci20 > 1% of its peak at t≤10s.
        # During nucleation Ci20 ~ 1e-89→1e-22 (near numerical floor); relative
        # errors there are large but physically irrelevant.
        err_Ci20_early = rel_err_max_early(Ci20_w, Ci20_f, t_f,
                                           t_cutoff=10.0, sig_frac=0.01)
        all_pass &= check(err_Ci20_early < 0.10,
                          f"Ci20 max rel-error (t≤10s, >1% peak) = {err_Ci20_early:.3e}  (tol 10%)")
        err_Ci20_all = rel_err_max(Ci20_w, Ci20_f)
        print(f"         Ci20 max rel-error (all t)  = {err_Ci20_all:.3e}"
              f"  (informational; transient differences expected as front crosses boundary)")

    # mean_i at early time
    mi_f = res_full['mean_sizes']['mean_i']
    mi_w = res_win['mean_sizes']['mean_i']
    err_mi_early = rel_err_max_early(mi_w, mi_f, t_f, t_cutoff=10.0)
    all_pass &= check(err_mi_early < 0.05,
                      f"mean_i max rel-error (t≤10s) = {err_mi_early:.3e}  (tol 5%)")

    speedup = t_full / max(t_win, 1e-3)
    print(f"\n  Wall time: full={t_full:.2f}s  window={t_win:.2f}s  speedup={speedup:.1f}×")
    all_pass &= check(t_win <= t_full * 1.5, "Window not significantly slower than full")

    # ── save comparison plots ─────────────────────────────────────────────────
    run_dir1 = create_run_directory(OUTPUT_DIR)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f'TEST 1  |  Nv={NV1} Ni={NI1}  |  full (blue) vs window (orange)',
                 fontsize=11)
    for ax, name, title in zip(
        axes.flat,
        ['Cv1', 'Ci1', 'Ci10', 'Ci20'],
        ['Cv1', 'Ci1', 'Ci10', 'Ci20'],
    ):
        cf = res_full['concentrations'].get(name)
        cw = res_win['concentrations'].get(name)
        if cf is not None:
            ax.loglog(t_f, np.maximum(cf, 1e-30), 'b-',  lw=1.5, label='full')
        if cw is not None:
            ax.loglog(t_w, np.maximum(cw, 1e-30), 'r--', lw=1.5, label='window')
        ax.set_title(title); ax.set_xlabel('t (s)'); ax.set_ylabel('C (at/at)')
        ax.legend(fontsize=8); ax.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    fig.savefig(run_dir1 / 'test1_comparison.png', dpi=120)
    plt.close(fig)
    print(f"\n  Plot saved → {run_dir1 / 'test1_comparison.png'}")

print(f"\nTEST 1  {'PASSED' if all_pass else 'FAILED'}")

# ══════════════════════════════════════════════════════════════════════════════
# TEST 2 – Reference configuration: Nv=200, Ni=1000 (matches provenance.md)
# ══════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('TEST 2  –  Reference configuration  (Nv=200, Ni=1000)')
print('  t=(1e-3, 1e6)  n_pts=1000  rtol=1e-4  atol=1e-20')
print('  Window solver  (w0_v=100, w0_i=100)')
print('  Cv1_ss = 3.264e-6 (verified against full solver, wall=241.7 s)')
print('  Reference full-solver wall time: 238.9 s')
print('='*70)

NV2, NI2 = 200, 1000
cfg_ref_win = {
    't_span': (1e-3, 1e6), 'n_points': 1000, 'rtol': 1e-4, 'atol': 1e-20,
    'log_time': True,
    'solver_method': {
        'backend': 'cvode', 'lmm': 'bdf',
        'window_mode': 1,
        'window_w0_v': 100, 'window_w0_i': 100,
        'window_C_expand': 1e-18, 'window_expand_pad': 10,
        'window_check_every': 1,
    },
}

sim2 = ClusterDynamicsSimulation(Nv=NV2, Ni=NI2)
t0 = time.perf_counter()
res_ref = run_cpp_solver(sim2, cfg_ref_win)
t_ref_win = time.perf_counter() - t0
speedup2 = 238.9 / max(t_ref_win, 1e-3)
print(f"\n  Window solver wall time: {t_ref_win:.2f} s  "
      f"|  speedup vs reference: {speedup2:.1f}×")

all_pass2 = True

if res_ref is None:
    print(f"  {FAIL}  Solver returned None")
    all_pass2 = False
else:
    t_arr = res_ref['time']
    conc  = res_ref['concentrations']

    all_pass2 &= check(len(t_arr) == 1000,
                       f"Got 1000 output points (got {len(t_arr)})")
    all_pass2 &= check(np.isclose(t_arr[0],  1e-3, rtol=1e-3), "t[0] ≈ 1e-3 s")
    all_pass2 &= check(np.isclose(t_arr[-1], 1e6,  rtol=1e-3), "t[-1] ≈ 1e6 s")

    Cv1 = conc['Cv1']
    Ci1 = conc['Ci1']
    all_pass2 &= check(np.all(Cv1 > 0) and np.all(np.isfinite(Cv1)),
                       "Cv1 positive and finite everywhere")
    all_pass2 &= check(np.all(Ci1 > 0) and np.all(np.isfinite(Ci1)),
                       "Ci1 positive and finite everywhere")

    # Cv1 steady-state: full solver gives 3.264e-6; window solver should agree to 2%
    CV1_SS_REF = 3.264e-6
    Cv1_ss = float(Cv1[-1])
    err_ss = abs(Cv1_ss - CV1_SS_REF) / CV1_SS_REF
    all_pass2 &= check(err_ss < 0.02,
                       f"Cv1 steady-state = {Cv1_ss:.4e}  "
                       f"(ref={CV1_SS_REF:.4e}, rel-err={err_ss:.3e}, tol 2%)")

    # Ci1 steady-state: full solver gives 3.133e-14
    CI1_SS_REF = 3.133e-14
    Ci1_ss = float(Ci1[-1])
    err_ci1 = abs(Ci1_ss - CI1_SS_REF) / CI1_SS_REF
    all_pass2 &= check(err_ci1 < 0.02,
                       f"Ci1 steady-state = {Ci1_ss:.4e}  "
                       f"(ref={CI1_SS_REF:.4e}, rel-err={err_ci1:.3e}, tol 2%)")

    # Active interstitial band: should extend to large sizes at late times
    active = res_ref['active_band']
    x_max_final = int(active['x_max'][-1])
    all_pass2 &= check(x_max_final > 100,
                       f"Interstitial front at size {x_max_final} at t=1e6 s  (>100)")

    all_pass2 &= check(t_ref_win < 238.9,
                       f"Window ({t_ref_win:.1f}s) faster than reference full solver (238.9s)")

    # ── plots ─────────────────────────────────────────────────────────────────
    run_dir2 = create_run_directory(OUTPUT_DIR)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.loglog(t_arr, np.maximum(Cv1, 1e-30), 'b-', lw=1.5, label='Cv1')
    ax.loglog(t_arr, np.maximum(Ci1, 1e-30), 'r-', lw=1.5, label='Ci1')
    ax.axhline(CV1_SS_REF, color='b', ls=':', lw=1, alpha=0.6,
               label=f'Cv1 ref = {CV1_SS_REF:.2e}')
    ax.set_xlabel('Time (s)'); ax.set_ylabel('Concentration (at/at)')
    ax.set_title(f'Point defects  |  Nv={NV2} Ni={NI2}  |  window BDF GMRES'
                 f'  |  {t_ref_win:.1f}s  ({speedup2:.1f}× speedup)')
    ax.legend(); ax.grid(True, which='both', alpha=0.3)
    plt.tight_layout(); fig.savefig(run_dir2 / 'point_defects.png', dpi=120)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    plot_sizes = [2, 6, 10, 20, 50, 100, 200, 500, 1000]
    cmap = plt.cm.plasma
    for j, x in enumerate(plot_sizes):
        name = f'Ci{x}'
        if name in conc and x <= NI2:
            ax.loglog(t_arr, np.maximum(conc[name], 1e-30),
                      color=cmap(j/len(plot_sizes)), lw=1, label=name)
    ax.set_xlabel('Time (s)'); ax.set_ylabel('Concentration (at/at)')
    ax.set_title(f'Interstitial clusters  |  Nv={NV2} Ni={NI2}  |  window BDF GMRES')
    ax.legend(fontsize=7, ncol=2); ax.grid(True, which='both', alpha=0.3)
    plt.tight_layout(); fig.savefig(run_dir2 / 'interstitial_clusters.png', dpi=120)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.semilogx(t_arr, active['x_max'], 'b-', lw=1.5, label='x_max active front')
    ax.semilogx(t_arr, active['x_min'], 'r--', lw=1, label='x_min')
    ax.axhline(100, color='gray', ls=':', lw=1, label='w0_i=100')
    ax.set_xlabel('Time (s)'); ax.set_ylabel('Cluster size')
    ax.set_title('Active interstitial band vs time  (window solver)')
    ax.legend(); ax.grid(True, which='both', alpha=0.3)
    plt.tight_layout(); fig.savefig(run_dir2 / 'active_band.png', dpi=120)
    plt.close(fig)

    print(f"\n  Plots saved → {run_dir2}/")

print(f"\nTEST 2  {'PASSED' if all_pass2 else 'FAILED'}")

# ══════════════════════════════════════════════════════════════════════════════
# TEST 3 – Constant-width sliding window  (Nv=1000, Ni=10000, t=1e-3→1e7)
#   Phase I   (window_mode=1): upper truncation only — x_min fixed at 2
#   Phase III (window_mode=3): constant-width W=500 window sliding upward;
#                               lower bound = max(2, x_hi_i − W + 1) after t>10 s
#                               auto-enabled only when N_EQ > 1000
#
#   Key diagnostic: frozen species are *exactly* constant in the C++ output
#   (unpack3 never writes to frozen slots), so c_curr == c_prev in float
#   identifies frozen clusters without any tolerance tuning.
# ══════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print('TEST 3  –  Constant-width sliding window  (Nv=1000, Ni=10000)')
print('  t=(1e-3, 1e7)  n_pts=500')
print('  Phase I   (window_mode=1)  fixed x_min=2,                    dynamic x_max')
print('  Phase III (window_mode=3)  const-width W=500, t_start=10 s,  dynamic x_max')
print('='*70)

NV3, NI3 = 1000, 10000
T_SPAN3   = (1e-3, 1e7)
N_PTS3    = 500

_shared_method3 = {
    'backend':              'cvode',
    'lmm':                  'bdf',
    'window_w0_v':          200,   # start with 200/1000 vacancy clusters
    'window_w0_i':          500,   # start with 500/10000 interstitial clusters
    'window_C_expand':      1e-18,
    'window_expand_pad':    50,    # additive pad when geometric step is small
    'window_expand_factor': 2.0,   # geometric: 500→1000→2000→4000→8000→10000
    'window_check_every':   10,
}
cfg3_I = {
    't_span': T_SPAN3, 'n_points': N_PTS3, 'rtol': 1e-4, 'atol': 1e-20,
    'log_time': True,
    'solver_method': {**_shared_method3, 'window_mode': 1},
}
cfg3_III = {
    't_span': T_SPAN3, 'n_points': N_PTS3, 'rtol': 1e-4, 'atol': 1e-20,
    'log_time': True,
    'solver_method': {
        **_shared_method3,
        'window_mode':      3,
        'window_width':     500,   # constant window width in i-cluster space
        'window_t_start':   10.0,  # suppress lower sliding until t > 10 s
        'window_N_thresh':  1000,  # activate only when N_EQ > 1000 (always here)
    },
}

sim3 = ClusterDynamicsSimulation(Nv=NV3, Ni=NI3)

print('\n  ── Phase I   (fixed lower front) ──')
t0 = time.perf_counter()
res3_I   = run_cpp_solver(sim3, cfg3_I)
t3_I     = time.perf_counter() - t0

print('\n  ── Phase III (constant-width sliding window) ──')
t0 = time.perf_counter()
res3_III = run_cpp_solver(sim3, cfg3_III)
t3_III   = time.perf_counter() - t0

speedup3 = t3_I / max(t3_III, 1e-3)
print(f"\n  Phase I   wall time : {t3_I:.2f} s")
print(f"  Phase III wall time : {t3_III:.2f} s")
print(f"  Speedup             : {speedup3:.2f}×")

all_pass3 = True

if res3_I is None or res3_III is None:
    print(f"  {FAIL}  One or both solvers returned None")
    all_pass3 = False
else:
    t_arr3    = res3_I['time']
    conc3_I   = res3_I['concentrations']
    conc3_III = res3_III['concentrations']

    # ── Cv1 steady-state: Phase I and Phase III should agree within 5% ────────
    Cv1_I_ss   = float(conc3_I['Cv1'][-1])
    Cv1_III_ss = float(conc3_III['Cv1'][-1])
    err_cv1    = abs(Cv1_I_ss - Cv1_III_ss) / max(abs(Cv1_I_ss), 1e-30)
    all_pass3 &= check(err_cv1 < 0.05,
                       f"Cv1 steady-state: PhaseI={Cv1_I_ss:.4e}  "
                       f"PhaseIII={Cv1_III_ss:.4e}  rel-err={err_cv1:.3e}  (tol 5%)")

    # ── Upper front: both modes should reach the full Ni window ───────────────
    xmax_I   = res3_I['active_band']['x_max']
    xmax_III = res3_III['active_band']['x_max']
    all_pass3 &= check(int(xmax_I[-1])   > 100,
                       f"Phase I   x_max final = {int(xmax_I[-1])}  (>100)")
    all_pass3 &= check(int(xmax_III[-1]) > 100,
                       f"Phase III x_max final = {int(xmax_III[-1])}  (>100)")

    # ── Infer Phase III x_lo_i from exact-constancy of frozen concentrations ──
    # Frozen species are NEVER written by unpack3 → bitwise-constant in output.
    # Active species differ between consecutive rows (CVODE updates them).
    conc_i_III = np.stack([
        conc3_III.get(f'Ci{k}', np.zeros(N_PTS3)) for k in range(1, NI3 + 1)
    ])  # shape (NI3, N_PTS3);  row k-1 = Ci_k

    x_lo_inferred = np.full(N_PTS3, 2, dtype=int)
    for j in range(1, N_PTS3):
        c_prev = conc_i_III[1:, j - 1]   # Ci2, Ci3, ... at previous output point
        c_curr = conc_i_III[1:, j]        # ... at current output point
        is_active = (c_curr != c_prev) & (c_curr > 1e-25)
        first_active = np.where(is_active)[0]
        if len(first_active):
            x_lo_inferred[j] = int(first_active[0]) + 2  # +1 (0→1 index) +1 (Ci1 offset)
        else:
            x_lo_inferred[j] = x_lo_inferred[j - 1]

    x_lo_final = int(x_lo_inferred[-1])
    all_pass3 &= check(x_lo_final > 2,
                       f"Phase III lower front advanced: final x_lo_i={x_lo_final}  (>2)")
    all_pass3 &= check(speedup3 >= 1.0,
                       f"Phase III not slower than Phase I  (speedup={speedup3:.2f}×)")

    width_I_final   = int(xmax_I[-1])   - 2            + 1
    width_III_final = int(xmax_III[-1]) - x_lo_final   + 1
    print(f"\n  Phase III lower front: started at x_lo=2 → ended at x_lo={x_lo_final}")
    print(f"  Phase I   active i-equations (final): {width_I_final}")
    print(f"  Phase III active i-equations (final): {width_III_final}  "
          f"({100*(1 - width_III_final/max(width_I_final,1)):.0f}% reduction)")

    # ── Figure 1: active band extents + active window width ───────────────────
    run_dir3 = create_run_directory(OUTPUT_DIR)

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle(f'Phase III const-width sliding window  (W=500)  |  {t3_III:.1f} s\n'
                 f'Nv={NV3} Ni={NI3}   t=(1e-3, 1e7)   '
                 f'speedup vs Phase I: {speedup3:.2f}×', fontsize=11)

    # Top: band extents (x_lo to x_max)
    ax_top.fill_between(t_arr3, 2,             xmax_I,          alpha=0.18,
                        color='steelblue', label='Phase I active band')
    ax_top.fill_between(t_arr3, x_lo_inferred, xmax_III,        alpha=0.25,
                        color='tomato',    label='Phase III active band')
    ax_top.semilogx(t_arr3, xmax_I,            'b-',  lw=1.5, label='Phase I  x_max')
    ax_top.semilogx(t_arr3, xmax_III,          'r-',  lw=1.5, label='Phase III x_max')
    ax_top.axhline(2, color='steelblue', ls='--', lw=0.8, label='Phase I  x_min = 2 (fixed)')
    ax_top.semilogx(t_arr3, x_lo_inferred,     'r--', lw=1.2,
                    label='Phase III x_min (const-width, inferred)')
    ax_top.set_ylabel('Cluster size')
    ax_top.set_title('Active interstitial band  (x_min to x_max)')
    ax_top.legend(fontsize=8, ncol=2, loc='upper left')
    ax_top.grid(True, which='both', alpha=0.3)
    ax_top.set_ylim(bottom=0)

    # Bottom: active window width
    width_I   = np.maximum(xmax_I   - 2,             0) + 1
    width_III = np.maximum(xmax_III - x_lo_inferred, 0) + 1
    ax_bot.semilogx(t_arr3, width_I,   'b-', lw=1.5,
                    label=f'Phase I   (final = {int(width_I[-1])} eqns)')
    ax_bot.semilogx(t_arr3, width_III, 'r-', lw=1.5,
                    label=f'Phase III (final = {int(width_III[-1])} eqns)')
    ax_bot.set_xlabel('Time (s)')
    ax_bot.set_ylabel('Active i-cluster equations')
    ax_bot.set_title('Active window width  (x_max − x_min + 1)')
    ax_bot.legend(fontsize=9)
    ax_bot.grid(True, which='both', alpha=0.3)
    ax_bot.set_ylim(bottom=0)

    plt.tight_layout()
    fig.savefig(run_dir3 / 'test3_sliding_window.png', dpi=120)
    plt.close(fig)

    # ── Figure 2: point defects Phase I vs Phase III ──────────────────────────
    fig2, ax2 = plt.subplots(figsize=(9, 4))
    ax2.loglog(t_arr3, np.maximum(conc3_I['Cv1'],   1e-30), 'b-',  lw=1.5,
               label='Cv1  Phase I')
    ax2.loglog(t_arr3, np.maximum(conc3_III['Cv1'], 1e-30), 'r--', lw=1.5,
               label='Cv1  Phase III')
    ax2.loglog(t_arr3, np.maximum(conc3_I['Ci1'],   1e-30), 'b:',  lw=1.5,
               label='Ci1  Phase I')
    ax2.loglog(t_arr3, np.maximum(conc3_III['Ci1'], 1e-30), 'r:',  lw=1.5,
               label='Ci1  Phase III')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Concentration (at/at)')
    ax2.set_title(f'Point defects: Phase I vs Phase III  |  '
                  f'Phase I = {t3_I:.1f} s   Phase III = {t3_III:.1f} s   '
                  f'speedup = {speedup3:.2f}×')
    ax2.legend(fontsize=9)
    ax2.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    fig2.savefig(run_dir3 / 'test3_point_defects.png', dpi=120)
    plt.close(fig2)

    print(f"\n  Plot → {run_dir3}/test3_sliding_window.png")
    print(f"  Plot → {run_dir3}/test3_point_defects.png")

print(f"\nTEST 3  {'PASSED' if all_pass3 else 'FAILED'}")

# ══════════════════════════════════════════════════════════════════════════════
print('\n' + '='*70)
print(f"OVERALL: TEST 1 {'PASSED' if all_pass else 'FAILED'}  |  "
      f"TEST 2 {'PASSED' if all_pass2 else 'FAILED'}  |  "
      f"TEST 3 {'PASSED' if all_pass3 else 'FAILED'}")
print('='*70)
sys.exit(0 if (all_pass and all_pass2 and all_pass3) else 1)
