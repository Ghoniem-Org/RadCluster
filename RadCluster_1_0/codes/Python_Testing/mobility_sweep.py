#!/usr/bin/env python3
"""
mobility_sweep.py — diagnostic sweep of SIA / vacancy mobility cutoffs.

Hypothesis under test
---------------------
With the current input (n_max_i=50, m_max_v=5, L_hat=50), mobile SIA
clusters of size 4–50 glide in 1D and are rapidly drained at fixed sinks
(dislocations, GBs, precipitates) and at voids — even when the 3D
coalescence diffusivity D_SIA_eff[n>=2] is zero.  Vacancy clusters lose
mobility at m=6 and become permanent sinks for the residual flux.  The
asymmetry should disappear (and may invert) when i_mobile is reduced
toward 1, when v_mobile is raised, or when L_hat is shortened.

Cases
-----
  A  baseline                   i_mobile=50, v_mobile=5,  L_hat=50
  B  di-SIA only                i_mobile=2,  v_mobile=5,  L_hat=50
  C  monomer-only SIA           i_mobile=1,  v_mobile=5,  L_hat=50
  D  symmetric mobility         i_mobile=50, v_mobile=20, L_hat=50
  E  mono-defect (reference)    i_mobile=1,  v_mobile=1,  L_hat=50
  F  short MFP (1D suppressed)  i_mobile=50, v_mobile=5,  L_hat=5

Output
------
  output/mobility_sweep_<timestamp>/
      summary.csv                — per-case scalars at 0.1 dpa
      mean_sizes_vs_dose.png     — mean_n_i, mean_n_v(t) for all cases
      swelling_vs_dose.png       — swelling(t) for all cases
      conservation_vs_dose.png   — delta_FP vs dose for all cases
      size_dist_at_0.1dpa.png    — c_n, c_m at 0.1 dpa, all cases
"""

import sys, os, io, csv, time
from pathlib import Path
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MODULE_ROOT = Path(__file__).resolve().parent.parent.parent
REPO_ROOT   = MODULE_ROOT.parent
for p in [str(REPO_ROOT), str(MODULE_ROOT)]:
    if p not in sys.path:
        sys.path.insert(0, p)

import importlib
import RadCluster_1_0.py_utils.defect_production as _dp_mod
import RadCluster_1_0.py_utils.binding_energies   as _be_mod
import RadCluster_1_0.py_utils.bin_moment_rates   as _bmr
import RadCluster_1_0.py_utils.input_data         as _inp_mod
import RadCluster_1_0.py_utils.reaction_rates     as _rr_mod
import RadCluster_1_0.py_utils.rate_equations     as _re_mod
import RadCluster_1_0.py_utils.cpp_bridge         as _cb_mod
import RadCluster_1_0.py_utils.post_process       as _pp_mod
import RadCluster_1_0.py_utils.simulation         as _sim_mod
for _m in [_dp_mod, _be_mod, _bmr, _inp_mod, _rr_mod, _re_mod,
           _cb_mod, _pp_mod, _sim_mod]:
    importlib.reload(_m)
from RadCluster_1_0.py_utils.simulation import RadClusterSimulation


# ── Solver / domain configuration ──────────────────────────────────────────
# I,V = 1000 each → N_eq = I + V + 6 ≈ 2006 (full_system + Woodbury).
# This matches the original reference run.  Vacancy clusters in case A may
# bump the V=1000 boundary at 0.01 dpa; delta_FP will quantify any leakage.
I_SIM, V_SIM = 1000, 1000

# Reference dose for the snapshot column in summary.csv / size_dist plot.
# 0.01 dpa is past the nucleation transient (~1e-9 dpa) and well into
# steady-state cluster growth — the SIA/vac size asymmetry is fully developed.
DOSE_SNAPSHOT = 0.01

SOLVER_CONFIG = {
    't_span':   (1e-6, 1e4),    # 1e4 s @ G=1e-6 dpa/s = 0.01 dpa
    'n_points': 120,
    'log_time': True,
    'rtol':     1e-6,
    'atol':     1e-22,
    'solver_method': {
        'linsol':             'gmres',  # required for Woodbury auto-engage
        # window_mode=0 (full_system) so Woodbury preconditioner kicks in
        # whenever i_mobile>=2 or v_mobile>=2  (parameters.h:494-501)
        'window_gmres_maxl':  30,
    },
}

CASES = [
    ('A_baseline',     {'i_mobile': 50, 'v_mobile': 5,  'L_hat': 50.0}),
    ('B_di_SIA',       {'i_mobile': 2,  'v_mobile': 5,  'L_hat': 50.0}),
    ('C_mono_SIA',     {'i_mobile': 1,  'v_mobile': 5,  'L_hat': 50.0}),
    ('D_symmetric',    {'i_mobile': 50, 'v_mobile': 20, 'L_hat': 50.0}),
    ('E_mono_defect',  {'i_mobile': 1,  'v_mobile': 1,  'L_hat': 50.0}),
    ('F_short_MFP',    {'i_mobile': 50, 'v_mobile': 5,  'L_hat': 5.0}),
]


def run_case(name, overrides, quiet=True, heartbeat_every=1.5):
    """Run one case via active_window and return time-resolved results.

    Per-row progress is forwarded to the *real* terminal even while the
    inner solver's stdout is captured.  A line is printed at most once
    every `heartbeat_every` seconds (wall-clock) plus always at the
    final point, to avoid flooding while still confirming progress.
    """
    i_mob = int(overrides['i_mobile'])
    v_mob = int(overrides['v_mobile'])
    L_hat = float(overrides['L_hat'])

    real_stdout = sys.stdout
    real_stderr = sys.stderr
    state = {'last_print': 0.0, 'rows': 0, 'G_dpa': None, 't0': time.time()}

    def _hb(row):
        # called by C++ bridge once per output time step
        state['rows'] += 1
        now = time.time()
        if (now - state['last_print']) < heartbeat_every:
            return
        state['last_print'] = now
        t   = row.get('t', 0.0)
        G   = state['G_dpa'] or 1e-6
        dose = t * G
        ci1 = row.get('c_i1', 0.0)
        cv1 = row.get('c_v1', 0.0)
        elapsed = now - state['t0']
        real_stdout.write(
            f"    [{name}] step {state['rows']:>4d}  "
            f"t={t:.2e}s  dose={dose:.3e} dpa  "
            f"c_i1={ci1:.2e}  c_v1={cv1:.2e}  "
            f"({elapsed:5.1f}s elapsed)\n"
        )
        real_stdout.flush()

    saved_stdout, saved_stderr = sys.stdout, sys.stderr
    if quiet:
        sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=I_SIM, V=V_SIM,
            solver_mode='full_system',
            physics_option='full_CD_fission',
            C_floor=1e-25,
            he_kinetics='quasi_steady_state',
            i_mobile=i_mob,
            v_mobile=v_mob,
        )
        inp = sim.input_data
        inp.diffusion['L_hat']    = L_hat
        inp.diffusion['i_mobile'] = i_mob
        inp.diffusion['v_mobile'] = v_mob
        inp.reactions['i_mobile']   = i_mob
        inp.reactions['v_mobile']   = v_mob
        inp.reactions['i_discrete'] = i_mob
        inp.reactions['v_discrete'] = max(v_mob, 5)
        inp.reactions['I_bin']      = 0
        inp.reactions['V_bin']      = 0
        inp._calculate_derived()
        sim.rebuild_rates()
        state['G_dpa'] = sim.input_data.derived['G']

        results = sim.run_adaptive(
            solver_config=SOLVER_CONFIG, save_output=False,
            progress_callback=_hb, boundary_threshold=0.1,
            max_doublings=1, points_per_segment=10,
        )
    finally:
        if quiet:
            sys.stdout, sys.stderr = saved_stdout, saved_stderr

    if results is None:
        return None

    dose = np.asarray(results['dose'])
    j01  = int(np.argmin(np.abs(dose - DOSE_SNAPSHOT)))

    # Reconstruct per-size distributions at the snapshot dose
    y_all = np.asarray(results['y'])
    s2m   = 1.0 / sim.input_data.derived['Omega']
    I_run = sim.input_data.I
    V_run = sim.input_data.V
    c_n_01 = np.maximum(y_all[0:I_run, j01], 0.0) * s2m            # m^-3
    c_v_01 = np.maximum(y_all[I_run:I_run + V_run, j01], 0.0) * s2m

    return {
        'name':      name,
        'overrides': overrides,
        't':         np.asarray(results['t']),
        'dose':      dose,
        'mean_n_i':  np.asarray(results['mean_n_i']),
        'mean_n_v':  np.asarray(results['mean_n_v']),
        'swelling':  np.asarray(results['swelling']),
        'N_loops':   np.asarray(results['N_loops']),
        'N_voids':   np.asarray(results['N_voids']),
        'delta_FP':  np.asarray(results['delta_FP']),
        'idx_01':    j01,
        'c_n_01':    c_n_01,
        'c_v_01':    c_v_01,
        'I_run':     I_run,
        'V_run':     V_run,
    }


def write_summary_csv(cases, out_dir):
    path = out_dir / 'summary.csv'
    tag  = f"{DOSE_SNAPSHOT:g}"
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['case', 'i_mobile', 'v_mobile', 'L_hat',
                    'dose_dpa@end', f'mean_n_i@{tag}', f'mean_n_v@{tag}',
                    f'swelling@{tag}[%]', f'N_loops@{tag}[m^-3]',
                    f'N_voids@{tag}[m^-3]', f'delta_FP@{tag}'])
        for c in cases:
            if c is None:
                continue
            j  = c['idx_01']
            ov = c['overrides']
            w.writerow([
                c['name'],
                ov['i_mobile'], ov['v_mobile'], ov['L_hat'],
                f"{c['dose'][-1]:.4g}",
                f"{c['mean_n_i'][j]:.3f}",
                f"{c['mean_n_v'][j]:.3f}",
                f"{c['swelling'][j] * 100:.4f}",
                f"{c['N_loops'][j]:.3e}",
                f"{c['N_voids'][j]:.3e}",
                f"{c['delta_FP'][j]:.3e}",
            ])
    return path


def plot_time_series(cases, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
    for c in cases:
        if c is None:
            continue
        axes[0].plot(c['dose'], c['mean_n_i'], label=c['name'])
        axes[1].plot(c['dose'], c['mean_n_v'], label=c['name'])
    for ax, ttl in zip(axes, ['mean_n_i (SIA)', 'mean_n_v (vacancy)']):
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlabel('dose (dpa)'); ax.set_ylabel(ttl)
        ax.grid(True, which='both', alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle('Mean cluster size vs dose — mobility-cutoff sweep')
    fig.tight_layout()
    fig.savefig(out_dir / 'mean_sizes_vs_dose.png', dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for c in cases:
        if c is None:
            continue
        ax.plot(c['dose'], np.asarray(c['swelling']) * 100, label=c['name'])
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('dose (dpa)'); ax.set_ylabel('swelling (%)')
    ax.grid(True, which='both', alpha=0.3); ax.legend(fontsize=8)
    ax.set_title('Swelling vs dose')
    fig.tight_layout()
    fig.savefig(out_dir / 'swelling_vs_dose.png', dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for c in cases:
        if c is None:
            continue
        ax.plot(c['dose'], c['delta_FP'], label=c['name'])
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('dose (dpa)'); ax.set_ylabel('delta_FP')
    ax.axhline(1e-3, ls='--', c='k', alpha=0.5, label='1e-3 (target)')
    ax.grid(True, which='both', alpha=0.3); ax.legend(fontsize=8)
    ax.set_title('Frenkel-pair conservation vs dose')
    fig.tight_layout()
    fig.savefig(out_dir / 'conservation_vs_dose.png', dpi=150)
    plt.close(fig)


def plot_size_distributions(cases, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for c in cases:
        if c is None:
            continue
        ns = np.arange(1, c['I_run'] + 1)
        ms = np.arange(1, c['V_run'] + 1)
        axes[0].plot(ns, c['c_n_01'], label=c['name'], lw=1.0)
        axes[1].plot(ms, c['c_v_01'], label=c['name'], lw=1.0)
    for ax, ttl in zip(axes, [f'SIA c_n at {DOSE_SNAPSHOT} dpa',
                              f'Vacancy c_m at {DOSE_SNAPSHOT} dpa']):
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlabel('cluster size')
        ax.set_ylabel('concentration (m$^{-3}$)')
        ax.grid(True, which='both', alpha=0.3); ax.legend(fontsize=8)
        ax.set_title(ttl)
        ax.set_ylim(bottom=1e10)
    fig.tight_layout()
    fig.savefig(out_dir / f'size_dist_at_{DOSE_SNAPSHOT:g}dpa.png', dpi=150)
    plt.close(fig)


def main():
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    out_dir = MODULE_ROOT / 'output' / f'mobility_sweep_{timestamp}'
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output -> {out_dir}", flush=True)

    results = []
    sweep_t0 = time.time()
    for k, (name, ov) in enumerate(CASES, start=1):
        clock = datetime.now().strftime('%H:%M:%S')
        print(f"\n=== [{k}/{len(CASES)}] {name}  "
              f"(i_mobile={ov['i_mobile']}, v_mobile={ov['v_mobile']}, "
              f"L_hat={ov['L_hat']})  start {clock} ===", flush=True)
        case_t0 = time.time()
        try:
            r = run_case(name, ov, quiet=True)
        except Exception as e:
            print(f"  FAILED: {type(e).__name__}: {e}", flush=True)
            r = None
        case_dt = time.time() - case_t0
        if r is None:
            print(f"  (no result; {case_dt:.1f}s)", flush=True)
        else:
            j = r['idx_01']
            print(f"  -> mean_n_i = {r['mean_n_i'][j]:.2f}   "
                  f"mean_n_v = {r['mean_n_v'][j]:.2f}   "
                  f"swelling = {r['swelling'][j]*100:.4f}%   "
                  f"delta_FP = {r['delta_FP'][j]:.2e}   "
                  f"({case_dt:.1f}s)", flush=True)
        results.append(r)
        # Incremental write so partial results survive a kill / crash
        try:
            write_summary_csv(results, out_dir)
        except Exception:
            pass

    total_dt = time.time() - sweep_t0
    csv_path = write_summary_csv(results, out_dir)
    print(f"\nSweep finished in {total_dt:.1f}s.  Wrote {csv_path}", flush=True)
    plot_time_series(results, out_dir)
    plot_size_distributions(results, out_dir)
    print(f"Plots written to {out_dir}", flush=True)


if __name__ == '__main__':
    main()
