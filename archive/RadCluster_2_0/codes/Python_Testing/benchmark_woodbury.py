#!/usr/bin/env python3
"""
benchmark_woodbury.py — Benchmark Woodbury vs Jacobi preconditioner.

Uses a SHORT time span (t=1e-6 to 0.1) to isolate the preconditioner
effect from the physics stiffness at later times. All 4 cases complete
in under a minute each.
"""

import sys, os, io, time, importlib, subprocess, json
from pathlib import Path

MODULE_ROOT = Path(__file__).resolve().parent.parent.parent
REPO_ROOT   = MODULE_ROOT.parent
for p in [str(REPO_ROOT), str(MODULE_ROOT)]:
    if p not in sys.path:
        sys.path.insert(0, p)

# ── Rebuild C++ solver ───────────────────────────────────────────────────────
build_dir = MODULE_ROOT / 'build'
build_dir.mkdir(exist_ok=True)
cmake_src = MODULE_ROOT / 'cpp_utils'
print("Building C++ solver...", flush=True)
for cmd in [
    ['cmake', '-S', str(cmake_src), '-B', str(build_dir), '-DCMAKE_BUILD_TYPE=Release'],
    ['cmake', '--build', str(build_dir), '--config', 'Release'],
]:
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"Build failed: {res.stderr[-1000:]}", flush=True)
        sys.exit(1)
print("Build OK\n", flush=True)

# ── Import ───────────────────────────────────────────────────────────────────
import RadCluster_2_0.py_utils.defect_production as _dp
import RadCluster_2_0.py_utils.binding_energies as _be
import RadCluster_2_0.py_utils.bin_moment_rates as _bmr
import RadCluster_2_0.py_utils.input_data as _inp
import RadCluster_2_0.py_utils.reaction_rates as _rr
import RadCluster_2_0.py_utils.rate_equations as _re
import RadCluster_2_0.py_utils.cpp_bridge as _cb
import RadCluster_2_0.py_utils.post_process as _pp
import RadCluster_2_0.py_utils.simulation as _sim
for _m in [_dp, _be, _bmr, _inp, _rr, _re, _cb, _pp, _sim]:
    importlib.reload(_m)
from RadCluster_2_0.py_utils.simulation import RadClusterSimulation

# ── Parameters ───────────────────────────────────────────────────────────────
I = 10_000
V = 10_000

BASE_OVERRIDES = {
    'eta': 0.3, 'f_cl_i': 0.0, 'f_cl_v': 0.5,
    'E_m_1D': 0.4, 'L_hat': 71.8,
    'c_C': 1.94e-4, 'E_b_C_SIA': 0.65,
    'rho_d': 1e16, 'Z_i': 1.05, 'Z_ii': 5,
    'i_discrete': I, 'v_discrete': V,
    'I_bin': 0, 'V_bin': 0,
    'shape_function': 'linear',
}


def make_config(prec_type, t_end):
    return {
        't_span': (1e-6, t_end),
        'n_points': 100,
        'log_time': True,
        'rtol': 1e-6,
        'atol': 1e-25,
        'solver_method': {
            'linsol': 'gmres',
            'window_width': 50,
            'concentration_threshold': 1e-20, 'window_pad': 100,
            'prec_type': prec_type,
        }
    }


def create_sim(i_mob, v_mob):
    saved = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=I, V=V, solver_mode='active_window',
            physics_option='full_CD_fission',
            C_floor=1e-25, he_kinetics='quasi_steady_state',
            i_mobile=i_mob, v_mobile=v_mob,
        )
        ov = dict(BASE_OVERRIDES)
        ov['i_mobile'] = i_mob
        inp = sim.input_data
        for key, val in ov.items():
            placed = False
            for d in [inp.production_fission, inp.production_fusion,
                      inp.diffusion, inp.reactions,
                      inp.energetics, inp.dissociation]:
                if key in d:
                    d[key] = val
                    placed = True
            if not placed:
                inp.reactions[key] = val
        inp.diffusion['i_mobile'] = int(ov['i_mobile'])
        inp.reactions['i_mobile'] = int(ov['i_mobile'])
        inp._calculate_derived()
        sim.rebuild_rates()
    finally:
        sys.stdout, sys.stderr = saved
    return sim


def run_case(label, i_mob, v_mob, prec_type, t_end):
    prec_name = 'Woodbury' if prec_type == 1 else 'Jacobi'
    print(f"\n{'='*60}", flush=True)
    print(f"  {label}", flush=True)
    print(f"  i_mobile={i_mob}, v_mobile={v_mob}, prec={prec_name}, t_end={t_end:.0e}", flush=True)
    print(f"{'='*60}", flush=True)

    sim = create_sim(i_mob, v_mob)
    cfg = make_config(prec_type, t_end)

    t0 = time.perf_counter()
    results = sim.run(solver_config=cfg, save_output=False)
    wall = time.perf_counter() - t0

    info = {'label': label, 'i_mobile': i_mob, 'v_mobile': v_mob,
            'prec_type': prec_type, 'wall_time': wall, 't_end': t_end}

    if results is not None:
        info['n_pts'] = len(results['t'])
        info['swelling'] = results['swelling'][-1] * 100
        info['delta_FP'] = results['delta_FP'][-1]
        print(f"  Wall time:  {wall:.2f} s", flush=True)
        print(f"  Swelling:   {info['swelling']:.6f} %", flush=True)
        print(f"  delta_FP:   {info['delta_FP']:.2e}", flush=True)
    else:
        print(f"  FAILED ({wall:.1f} s)", flush=True)
    return info


# ══════════════════════════════════════════════════════════════════════════════
# Run benchmarks at two time horizons
# ══════════════════════════════════════════════════════════════════════════════

all_results = []

for t_end in [1e-1, 1e0]:
    tag = f"t={t_end:.0e}"
    all_results.append(run_case(f"Mono+Woodbury ({tag})",   1, 1, 1, t_end))
    all_results.append(run_case(f"Mono+Jacobi ({tag})",     1, 1, 0, t_end))
    all_results.append(run_case(f"Coal+Woodbury ({tag})",  10, 5, 1, t_end))
    all_results.append(run_case(f"Coal+Jacobi ({tag})",    10, 5, 0, t_end))

# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n\n{'='*60}", flush=True)
print(f"  BENCHMARK SUMMARY  (I={I}, V={V})", flush=True)
print(f"{'='*60}", flush=True)
print(f"  {'Case':<35s} {'Wall(s)':>8s} {'Swell%':>10s} {'dFP':>10s}", flush=True)
print(f"  {'-'*33}  {'-'*8} {'-'*10} {'-'*10}", flush=True)
for r in all_results:
    sw = f"{r.get('swelling',0):.5f}" if 'swelling' in r else 'FAIL'
    dfp = f"{r.get('delta_FP',0):.1e}" if 'delta_FP' in r else 'FAIL'
    print(f"  {r['label']:<35s} {r['wall_time']:>8.2f} {sw:>10s} {dfp:>10s}", flush=True)

# Speedup analysis
print(f"\n  Speedup (Jacobi time / Woodbury time):", flush=True)
for t_end in [1e-1, 1e0]:
    tag = f"t={t_end:.0e}"
    for mob_label, im, vm in [("Monomer", 1, 1), ("Coalescence", 10, 5)]:
        wood = [r for r in all_results if r['i_mobile']==im and r['prec_type']==1 and r['t_end']==t_end]
        jac  = [r for r in all_results if r['i_mobile']==im and r['prec_type']==0 and r['t_end']==t_end]
        if wood and jac and wood[0]['wall_time'] > 0:
            ratio = jac[0]['wall_time'] / wood[0]['wall_time']
            print(f"    {mob_label} ({tag}): {ratio:.2f}x", flush=True)

print(f"\n  Coalescence overhead (Coal time / Mono time):", flush=True)
for t_end in [1e-1, 1e0]:
    tag = f"t={t_end:.0e}"
    for pt, pn in [(1, "Woodbury"), (0, "Jacobi")]:
        mono = [r for r in all_results if r['i_mobile']==1  and r['prec_type']==pt and r['t_end']==t_end]
        coal = [r for r in all_results if r['i_mobile']==10 and r['prec_type']==pt and r['t_end']==t_end]
        if mono and coal and mono[0]['wall_time'] > 0:
            ratio = coal[0]['wall_time'] / mono[0]['wall_time']
            print(f"    {pn} ({tag}): {ratio:.2f}x", flush=True)

print(flush=True)
