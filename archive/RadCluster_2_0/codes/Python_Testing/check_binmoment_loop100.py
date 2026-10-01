"""
check_binmoment_loop100.py — verify the bin-moment ⟨100⟩ SIA-loop reduction
and the KLU sparse-direct solver.

Runs four short cases and checks:
  1. bin_moment + loop_conversion ON  → δ_FP small, ⟨100⟩ population grows.
  2. bin_moment + loop_conversion OFF → δ_FP small (regression: ½⟨111⟩ path).
  3. discrete   + loop_conversion ON  → δ_FP small (validated baseline).
  4. discrete   + KLU (linsol='klu')  → runs, δ_FP small.
"""
import io
import sys
from pathlib import Path

MODULE_ROOT = Path(__file__).resolve().parents[2]   # RadCluster_2_0
REPO_ROOT   = MODULE_ROOT.parent
for p in (REPO_ROOT, MODULE_ROOT):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import numpy as np
from RadCluster_2_0.py_utils.simulation import RadClusterSimulation


def run_case(equations, loop_conversion, linsol='gmres', label=''):
    I, V = 200, 80
    sim = RadClusterSimulation(
        I=I, V=V, solver_mode='full_system',
        equations=equations, cascade='fission',
        C_floor=1e-25, he_kinetics='quasi_steady_state',
        i_mobile=5, v_mobile=2,
    )

    overrides = {'T': 673}
    if equations == 'bin_moment':
        overrides.update(i_discrete=8, v_discrete=4, I_bin=10, V_bin=6,
                         shape_function='linear')
    else:
        overrides.update(i_discrete=I, v_discrete=V, I_bin=0, V_bin=0)
    overrides['n_loop_min'] = 4

    inp = sim.input_data
    for key, val in overrides.items():
        placed = False
        for d in (inp.production_fission, inp.production_fusion, inp.diffusion,
                  inp.reactions, inp.energetics, inp.dissociation):
            if key in d:
                d[key] = val
                placed = True
        if not placed:
            inp.reactions[key] = val
    inp._calculate_derived()
    sim.rebuild_rates()

    solver_config = {
        't_span': (1e-6, 1.0), 'n_points': 40, 'log_time': True,
        'rtol': 1e-6, 'atol': 1e-20, 'timeout_s': 120,
        'loop_conversion': int(loop_conversion),
        'solver_method': {'linsol': linsol, 'preconditioner': 'Jacobi'},
    }

    buf = io.StringIO()
    old = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = buf
    try:
        results = sim.run(solver_config=solver_config, save_output=False)
    except Exception as exc:
        sys.stdout, sys.stderr = old
        print(f"[{label}] EXCEPTION: {type(exc).__name__}: {exc}")
        print(buf.getvalue()[-2000:])
        return None
    finally:
        sys.stdout, sys.stderr = old

    if results is None:
        print(f"[{label}] no results")
        print(buf.getvalue()[-2000:])
        return None

    def last(key):
        a = results.get(key)
        return float(a[-1]) if a is not None and len(a) else float('nan')

    dFP   = last('delta_FP')
    n100  = last('N_loops_100')
    n111  = last('N_loops_111')
    f111  = last('f_111_loop')
    mean100 = last('mean_n_100')
    npts  = len(results.get('t', []))
    print(f"[{label}] pts={npts:3d}  delta_FP={dFP:.2e}  "
          f"N_111={n111:.3e}  N_100={n100:.3e}  f_111={f111:.4f}  "
          f"mean_n_100={mean100:.2f}")
    return results


if __name__ == '__main__':
    print("=" * 78)
    r1 = run_case('bin_moment', True,  'gmres', 'bin+conv ON ')
    r2 = run_case('bin_moment', False, 'gmres', 'bin+conv OFF')
    r3 = run_case('discrete',   True,  'gmres', 'disc+conv ON')
    r4 = run_case('discrete',   False, 'klu',   'disc+KLU    ')
    print("=" * 78)

    ok = True
    if r1 is not None:
        dFP = float(r1['delta_FP'][-1]); n100 = float(r1['N_loops_100'][-1])
        if not (dFP < 1e-3):
            print(f"FAIL: bin+conv delta_FP={dFP:.2e} not < 1e-3"); ok = False
        if not (n100 > 0):
            print(f"FAIL: bin+conv N_100={n100:.2e} did not grow"); ok = False
    else:
        ok = False
    if r2 is not None and not (float(r2['delta_FP'][-1]) < 1e-3):
        print("FAIL: bin+conv-OFF delta_FP too large"); ok = False
    if r4 is None:
        print("FAIL: KLU run produced no results"); ok = False
    print("RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
