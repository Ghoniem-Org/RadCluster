import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.stderr.reconfigure(encoding='utf-8', errors='replace')
import time
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from py_utils.simulation import ClusterDynamicsSimulation
from py_utils.cpp_bridge  import run_cpp_solver

NV        = 1000
NI        = 10_000   # Ni_initial
NI_MAX    = 150_000  # pre-allocated upper bound

COMMON = {
    't_span':   (1e-5, 1e6),
    'n_points': 20,
    'rtol':     1e-4,
    'atol':     1e-30,
    'log_time': True,
}

WINDOW_COMMON = {
    'backend':              'cvode',
    'lmm':                  'bdf',
    'linsol':               'gmres',
    'window_w0_v':          500,
    'window_w0_i':          1000,
    'window_C_expand':      1e-18,
    'window_expand_pad':    1000,
    'window_expand_factor': 1.0,
    'window_check_every':   1,
    'window_width':         1000,
    'window_t_start':       10.0,
    'window_N_thresh':      2000,
    'window_prec':          1,
    # Dynamic Ni extension
    'Ni_max':               NI_MAX,
    'Ni_extend_margin':     500,   # extend when front is within 500 of P.Ni
    'Ni_extend_tol':        1e-4,  # conservation trigger
}

PHASE3 = dict(COMMON, solver_method=dict(WINDOW_COMMON, window_mode=3))
PHASE4 = dict(COMMON, solver_method=dict(WINDOW_COMMON,
                                          window_mode=4,
                                          window_omp_threads=12))

results = {}
for label, cfg in [("Phase III (sliding_window)", PHASE3),
                   ("Phase IV  (sliding_OpenMP)",  PHASE4)]:
    print("=" * 60)
    print(f"Running {label}  NV={NV}  NI={NI}  NI_MAX={NI_MAX}")
    print("=" * 60)
    sim = ClusterDynamicsSimulation(Nv=NV, Ni=NI)
    t0  = time.perf_counter()
    r   = run_cpp_solver(sim, cfg, base_dir=BASE_DIR)
    elapsed = time.perf_counter() - t0

    if r:
        cv1 = r['concentrations']['Cv1'][-1]
        ci1 = r['concentrations']['Ci1'][-1]
    else:
        cv1 = ci1 = float('nan')

    print(f"  Wall time : {elapsed:.2f} s")
    print(f"  Cv1(end)  : {cv1:.4e}")
    print(f"  Ci1(end)  : {ci1:.4e}")
    results[label] = {'time': elapsed, 'cv1': cv1}

print()
print("=" * 60)
print(f"BENCHMARK SUMMARY  |  NV={NV}  NI_initial={NI}  NI_MAX={NI_MAX}  t_end=1e6  n_points=20")
print("=" * 60)
t3  = results["Phase III (sliding_window)"]['time']
t4  = results["Phase IV  (sliding_OpenMP)"]['time']
cv3 = results["Phase III (sliding_window)"]['cv1']
cv4 = results["Phase IV  (sliding_OpenMP)"]['cv1']
print(f"  Phase III wall time : {t3:7.2f} s")
print(f"  Phase IV  wall time : {t4:7.2f} s")
print(f"  Speedup (III/IV)    :   {t3/t4:.2f}x")
err = abs(cv4-cv3)/abs(cv3)*100 if cv3 and cv3 == cv3 else float('nan')
print(f"  Cv1 rel. error      :   {err:.4f}%")
