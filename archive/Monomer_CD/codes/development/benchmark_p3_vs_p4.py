import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.stderr.reconfigure(encoding='utf-8', errors='replace')
import sys, time
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from py_utils.simulation    import ClusterDynamicsSimulation
from py_utils.cpp_bridge    import run_cpp_solver

NV = 1000
NI = 50000

COMMON = {
    't_span':   (1e-5, 1e7),
    'n_points': 1000,
    'rtol':     1e-4,
    'atol':     1e-20,
    'log_time': True,
}

PHASE3 = dict(COMMON, solver_method={
    'backend':              'cvode',
    'lmm':                  'bdf',
    'linsol':               'gmres',
    'window_mode':          3,
    'window_w0_v':          500,
    'window_w0_i':          1000,
    'window_C_expand':      1e-18,
    'window_expand_pad':    1000,
    'window_expand_factor': 1.0,
    'window_check_every':   10,
    'window_width':         1000,
    'window_t_start':       10.0,
    'window_N_thresh':      2000,
    'window_prec':          1,
})

PHASE4 = dict(COMMON, solver_method={
    'backend':              'cvode',
    'lmm':                  'bdf',
    'linsol':               'gmres',
    'window_mode':          4,
    'window_omp_threads':   12,
    'window_w0_v':          500,
    'window_w0_i':          1000,
    'window_C_expand':      1e-18,
    'window_expand_pad':    1000,
    'window_expand_factor': 1.0,
    'window_check_every':   10,
    'window_width':         1000,
    'window_t_start':       10.0,
    'window_N_thresh':      2000,
    'window_prec':          1,
})

results = {}
for label, cfg in [("Phase III (sliding_window)", PHASE3),
                   ("Phase IV  (sliding_OpenMP)", PHASE4)]:
    print("=" * 60)
    print(f"Running {label}  NV={NV}  NI={NI}")
    print("=" * 60)
    sim = ClusterDynamicsSimulation(Nv=NV, Ni=NI)
    t0 = time.perf_counter()
    r  = run_cpp_solver(sim, cfg, base_dir=BASE_DIR)
    elapsed = time.perf_counter() - t0
    cv1 = r['concentrations']['Cv1'][-1] if r else float('nan')
    ci1 = r['concentrations']['Ci1'][-1] if r else float('nan')
    band = r.get('active_band', {}) if r else {}
    hi_i = int(band.get('x_max', [0])[-1]) if band else 0
    lo_i = int(band.get('x_min', [0])[-1]) if band else 0
    print(f"  Wall time : {elapsed:.2f} s")
    print(f"  Cv1(end)  : {cv1:.4e}")
    print(f"  Ci1(end)  : {ci1:.4e}")
    if hi_i:
        print(f"  Active Ci : [{lo_i}..{hi_i}]  ({hi_i/NI*100:.1f}% of Ni)")
    results[label] = {'time': elapsed, 'cv1': cv1}

print()
print("=" * 60)
print("BENCHMARK SUMMARY  |  NV=%d  NI=%d  n_points=1000" % (NV, NI))
print("=" * 60)
t3 = results["Phase III (sliding_window)"]['time']
t4 = results["Phase IV  (sliding_OpenMP)"]['time']
cv3 = results["Phase III (sliding_window)"]['cv1']
cv4 = results["Phase IV  (sliding_OpenMP)"]['cv1']
print("  Phase III wall time : %7.2f s" % t3)
print("  Phase IV  wall time : %7.2f s" % t4)
print("  Speedup (III/IV)    :   %.2fx" % (t3/t4))
err = abs(cv4-cv3)/cv3*100 if cv3 else float('nan')
print("  Cv1 rel. error      :   %.4f%%" % err)
