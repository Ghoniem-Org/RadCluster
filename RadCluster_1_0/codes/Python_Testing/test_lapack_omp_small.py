#!/usr/bin/env python3
"""
test_lapack_omp_small.py — Quick verification of LAPACK + OMP auto-pick.

Small problem (I=V=1000, N_eq~2002), short time span (t_end=1e-2).
Runs once with Jacobi, once with Woodbury, using the already-built
solver.exe.  Does NOT rebuild — uses whatever's in build/Release/.

Verifies:
  1. Auto-thread picker fires (expect "auto-selected 4 threads ...")
  2. Woodbury works (no "binary built without LAPACK" fallback)
  3. Speedup of Woodbury over Jacobi
"""

import os
# Drop OMP_NUM_THREADS BEFORE importing any solver bridge so the spawned
# solver.exe does not inherit the machine-scope value.
os.environ.pop('OMP_NUM_THREADS', None)

import sys, io, time, importlib
from pathlib import Path

MODULE_ROOT = Path(__file__).resolve().parent.parent.parent
REPO_ROOT   = MODULE_ROOT.parent
for p in [str(REPO_ROOT), str(MODULE_ROOT)]:
    if p not in sys.path:
        sys.path.insert(0, p)

import RadCluster_1_0.py_utils.input_data    as _inp
import RadCluster_1_0.py_utils.reaction_rates as _rr
import RadCluster_1_0.py_utils.rate_equations as _re
import RadCluster_1_0.py_utils.cpp_bridge     as _cb
import RadCluster_1_0.py_utils.post_process   as _pp
import RadCluster_1_0.py_utils.simulation     as _sim
for _m in [_inp, _rr, _re, _cb, _pp, _sim]:
    importlib.reload(_m)
from RadCluster_1_0.py_utils.simulation import RadClusterSimulation

I = 1000
V = 1000

OVERRIDES = {
    'eta': 0.3, 'f_cl_i': 0.2, 'f_cl_v': 0.15,
    'E_m_1D': 0.4, 'i_mobile': 1, 'L_hat': 5,
    'c_C': 0.001, 'E_b_C_SIA': 0.65,
    'rho_d': 1e15, 'Z_i': 1.08, 'Z_ii': 1.01,
    'shape_function': 'linear',
    'i_discrete': I, 'v_discrete': V,
    'I_bin': 0, 'V_bin': 0,
}


def make_sim():
    saved = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=I, V=V, solver_mode='cpp_full',
            physics_option='full_CD_fission',
            C_floor=1e-25, he_options='quasi_steady_state',
            i_mobile=1, v_mobile=1,
        )
        inp = sim.input_data
        for key, val in OVERRIDES.items():
            placed = False
            for d in [inp.production_fission, inp.production_fusion,
                      inp.diffusion, inp.reactions,
                      inp.energetics, inp.dissociation]:
                if key in d:
                    d[key] = val
                    placed = True
            if not placed:
                inp.reactions[key] = val
        inp._calculate_derived()
        sim.rebuild_rates()
    finally:
        sys.stdout, sys.stderr = saved
    return sim


def run_case(prec_type):
    cfg = {
        't_span': (1e-6, 1e-2),
        'n_points': 50,
        'log_time': True,
        'rtol': 1e-6,
        'atol': 1e-25,
        'solver_method': {
            'backend': 'cvode', 'lmm': 'bdf', 'linsol': 'gmres',
            'window_w0_i': 50, 'window_width': 200,
            'window_C_expand': 1e-20, 'window_expand_pad': 100,
            'window_prec': 1, 'window_gmres_maxl': 30,
            'window_N_thresh': 500, 'prec_type': prec_type,
        }
    }
    name = 'Woodbury' if prec_type == 1 else 'Jacobi'
    print(f"\n{'='*60}\n  {name}  (prec_type={prec_type})\n{'='*60}", flush=True)
    sim = make_sim()
    t0 = time.perf_counter()
    results = sim.run(solver_config=cfg, save_output=False)
    wall = time.perf_counter() - t0
    if results is None:
        print(f"  FAILED after {wall:.2f}s", flush=True)
        return None
    md = results.get('metadata', {})
    n_thr = md.get('omp_threads_used', '?')
    print(f"  wall = {wall:.2f}s   omp_threads_used = {n_thr}", flush=True)
    return {'wall': wall, 'n_thr': n_thr, 'prec': name}


print(f"OMP_NUM_THREADS env at script start: "
      f"{os.environ.get('OMP_NUM_THREADS', '<unset>')}", flush=True)
print(f"N_eq for I=V={I}: ~{2 * I + 2}", flush=True)

j = run_case(0)
w = run_case(1)

if j and w:
    print(f"\n{'='*60}\n  SUMMARY\n{'='*60}")
    print(f"  Jacobi   wall = {j['wall']:.2f}s  ({j['n_thr']} threads)")
    print(f"  Woodbury wall = {w['wall']:.2f}s  ({w['n_thr']} threads)")
    if w['wall'] > 0:
        print(f"  Speedup (Jacobi / Woodbury) = {j['wall'] / w['wall']:.2f}x")
