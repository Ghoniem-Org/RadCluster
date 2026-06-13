"""He-conservation audit — runs all four physics options x both He-kinetics
modes at small system size and reports the delta_He / delta_FP diagnostics.

The conservation identity tested (Eq. delta_He):
    dynamic:  c_h(t) + Q(t) + J_He_sink(t) = c_h(0) + Q(0) + int_0^t G_He dt'
    QSS:      Q(t) + J_He_sink(t)          = Q(0)          + int_0^t G_He dt'
where Q = Q_tot (Case 2) or sum_m Q_m (Case 1), and J_He_sink is the CVODE
accounting integral of all He sink losses.  delta_He is the relative defect;
values <~ 1e-6 mean every He channel is balanced, > 1e-3 indicates a leak.

Each configuration runs in its own subprocess (pass "<equations> <cascade>
<he_kinetics>" as argv to run one config) so a failure in one cannot affect
the others.
"""
import sys, io, subprocess, traceback
from pathlib import Path

MODULE_ROOT = Path(__file__).resolve().parent.parent.parent
REPO_ROOT   = MODULE_ROOT.parent

I = V = 200
i_mobile, v_mobile = 5, 3

SOLVER_CONFIG = {
    't_span': (1e-6, 1e3),
    'n_points': 25,
    'log_time': True,
    'rtol': 1e-8, 'atol': 1e-30,
    'solver_method': {'linsol': 'gmres'},
}

BIN_KW = dict(i_discrete=i_mobile, v_discrete=v_mobile, I_bin=15, V_bin=15)


def run_single(eq, cas, hk):
    sys.path.insert(0, str(REPO_ROOT))
    sys.path.insert(0, str(MODULE_ROOT))
    import numpy as np
    from RadCluster_2_0.py_utils.simulation import RadClusterSimulation

    po = ('full_CD_' if eq == 'discrete' else 'bin_moment_CD_') + cas
    _saved = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(I=I, V=V, solver_mode='full_system',
                                   physics_option=po, C_floor=1e-25,
                                   he_kinetics=hk,
                                   i_mobile=i_mobile, v_mobile=v_mobile)
        inp = sim.input_data
        if eq == 'bin_moment':
            for k, v in BIN_KW.items():
                inp.reactions[k] = v
            inp.reactions['shape_function'] = 'linear'
            inp._calculate_derived()
            sim.rebuild_rates()
        results = sim.run_adaptive(
            solver_config=SOLVER_CONFIG, save_output=False,
            boundary_threshold=0.1, max_doublings=0, points_per_segment=10)
    finally:
        sys.stdout, sys.stderr = _saved

    if results is None:
        print('RESULT FAILED solver-returned-None')
        return 1
    dHe = np.asarray(results['delta_He'])
    dFP = np.asarray(results['delta_FP'])
    He_tot = np.asarray(results['C_He_tot'])
    # Ignore the start-up region where the He inventory is still ~0 and the
    # relative diagnostic is dominated by round-off of tiny numbers.
    mask = He_tot > (He_tot[-1] * 1e-6 if He_tot[-1] > 0 else 0.0)
    dHe_max = float(np.max(dHe[mask])) if mask.any() else float(np.max(dHe))
    print(f'RESULT OK {dHe_max:.3e} {dHe[-1]:.3e} {dFP[-1]:.3e}')
    return 0


def main():
    configs = [(eq, cas, hk)
               for eq in ('discrete', 'bin_moment')
               for cas in ('fission', 'fusion')
               for hk in ('quasi_steady_state', 'dynamic')]
    rows = []
    for eq, cas, hk in configs:
        label = f'{eq:>10s}/{cas:<7s}/{hk}'
        proc = subprocess.run([sys.executable, __file__, eq, cas, hk],
                              capture_output=True, text=True, timeout=1800)
        line = next((l for l in proc.stdout.splitlines()
                     if l.startswith('RESULT')), None)
        if line is None or 'FAILED' in line or proc.returncode != 0:
            rows.append((label, None))
            print(f'{label}:  FAILED', flush=True)
            tail = (proc.stderr or proc.stdout).strip().splitlines()
            for l in tail[-4:]:
                print(f'    {l}', flush=True)
            continue
        _, _, dmax, dlast, dfp = line.split()
        rows.append((label, float(dmax)))
        print(f'{label}:  max delta_He = {dmax}   final delta_He = {dlast}'
              f'   final delta_FP = {dfp}', flush=True)

    print('\n' + '=' * 78, flush=True)
    bad = [r for r in rows if r[1] is None or r[1] > 1e-6]
    if not bad:
        print(f'ALL {len(rows)} CONFIGS PASS: max delta_He <= 1e-6.', flush=True)
    else:
        print(f'{len(bad)} of {len(rows)} configs exceed delta_He = 1e-6:', flush=True)
        for r in bad:
            print(f'  {r[0]}: {r[1]}', flush=True)


if __name__ == '__main__':
    if len(sys.argv) == 4:
        sys.exit(run_single(*sys.argv[1:4]))
    main()
