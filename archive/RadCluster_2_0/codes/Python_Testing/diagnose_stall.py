#!/usr/bin/env python3
"""Diagnose bin-moment stall: run to t=8000 (just before stall) and inspect state."""

import sys, os, io, importlib, numpy as np
from pathlib import Path

MODULE_ROOT = Path(__file__).resolve().parent.parent.parent
REPO_ROOT   = MODULE_ROOT.parent
for p in [str(REPO_ROOT), str(MODULE_ROOT)]:
    if p not in sys.path:
        sys.path.insert(0, p)

for mod_name in ['defect_production', 'binding_energies', 'bin_moment_rates',
                 'input_data', 'reaction_rates', 'rate_equations',
                 'cpp_bridge', 'post_process', 'simulation']:
    m = importlib.import_module(f'RadCluster_2_0.py_utils.{mod_name}')
    importlib.reload(m)
from RadCluster_2_0.py_utils.simulation import RadClusterSimulation

I, V = 10000, 10000
i_mobile, v_mobile = 10, 5
i_discrete, v_discrete = 200, 50

saved = sys.stdout, sys.stderr
sys.stdout = sys.stderr = io.StringIO()
try:
    sim = RadClusterSimulation(
        I=I, V=V, solver_mode='full_system',
        physics_option='bin_moment_CD_fission',
        C_floor=1e-20, he_kinetics='quasi_steady_state',
        i_mobile=i_mobile, v_mobile=v_mobile,
    )
    overrides = {
        'eta': 0.3, 'f_cl_i': 0.6, 'f_cl_v': 0.15,
        'E_m_1D': 0.4, 'i_mobile': i_mobile, 'L_hat': 71.8,
        'c_C': 1.94e-4, 'E_b_C_SIA': 0.65,
        'rho_d': 1e15, 'Z_i': 1.05, 'Z_ii': 1.2,
        'shape_function': 'linear',
        'i_discrete': i_discrete, 'v_discrete': v_discrete,
        'r_ratio': 1.3,
    }
    inp = sim.input_data
    for key, val in overrides.items():
        placed = False
        for d in [inp.production_fission, inp.production_fusion,
                  inp.diffusion, inp.reactions,
                  inp.energetics, inp.dissociation]:
            if key in d:
                d[key] = val
                placed = True
        if not placed:
            inp.reactions[key] = val
    inp.diffusion['i_mobile'] = i_mobile
    inp.reactions['i_mobile'] = i_mobile
    inp._calculate_derived()
    sim.rebuild_rates()
finally:
    sys.stdout, sys.stderr = saved

re = sim.rate_equations
PM = re.n_mom
print(f"N_eq={re.N_eq} i_d={re.i_discrete} I_bin={re.I_bin} "
      f"v_d={re.v_discrete} V_bin={re.V_bin} PM={PM}", flush=True)

# Run to t=8000 (just before stall at ~9000)
config = {
    't_span': (1e-6, 8e3),
    'n_points': 100,
    'log_time': True,
    'rtol': 1e-5,
    'atol': 1e-20,
    'solver_method': {'linsol': 'dense'},
}

print("Running to t=8000...", flush=True)
results = sim.run(solver_config=config, save_output=False, timeout_s=60)

if results is None:
    print("FAILED — solver returned None", flush=True)
    sys.exit(1)

t = results['t']
y = results['y']
print(f"Completed {len(t)} points, t_final={t[-1]:.4e}", flush=True)

i_d = re.i_discrete
i_VAC = i_d + PM * re.I_bin

# Analyze SIA bin moments
print(f"\n=== SIA bin moments at t={t[-1]:.2e} ===", flush=True)
# Compute bin edges with r_ratio=1.3
edges = []
edge = i_d + 1
for k in range(re.I_bin):
    n_lo = edge
    n_hi = max(int(edge * 1.3), edge + 1)
    n_hi = min(n_hi, I + 1)
    edges.append((n_lo, n_hi))
    edge = n_hi

for k in range(re.I_bin):
    mu0 = y[i_d + PM*k, -1]
    mu1 = y[i_d + PM*k + 1, -1] if PM >= 2 else 0
    n_bar = mu1 / max(abs(mu0), 1e-300)
    n_lo, n_hi = edges[k]
    bw = n_hi - n_lo

    # Check linear reconstruction at bin edges
    S1 = sum(range(n_lo, n_hi))
    S2 = sum(n*n for n in range(n_lo, n_hi))
    det = bw * S2 - S1 * S1
    c_lo = c_hi = 0
    if abs(det) > 1e-30:
        phi0_lo = (S2 - S1 * n_lo) / det
        phi1_lo = (bw * n_lo - S1) / det
        c_lo = phi0_lo * mu0 + phi1_lo * mu1
        phi0_hi = (S2 - S1 * (n_hi-1)) / det
        phi1_hi = (bw * (n_hi-1) - S1) / det
        c_hi = phi0_hi * mu0 + phi1_hi * mu1

    flag = ""
    if mu0 < 0: flag += " MU0_NEG"
    if mu1 < 0 and mu0 > 1e-30: flag += " MU1_NEG"
    if c_lo < -1e-30: flag += f" c_lo={c_lo:.1e}_NEG"
    if c_hi < -1e-30: flag += f" c_hi={c_hi:.1e}_NEG"
    if mu0 > 1e-30 and (n_bar < n_lo * 0.8 or n_bar > n_hi * 1.2):
        flag += f" NBAR_OUT"

    print(f"  bin {k:2d}: [{n_lo:5d},{n_hi:5d}) bw={bw:4d}  "
          f"mu0={mu0:+10.2e}  mu1={mu1:+10.2e}  nbar={n_bar:8.1f}  "
          f"c_lo={c_lo:+.2e}  c_hi={c_hi:+.2e}{flag}", flush=True)

# Same for VAC
v_d = re.v_discrete
vac_start = i_VAC + v_d
print(f"\n=== VAC bin moments at t={t[-1]:.2e} ===", flush=True)
v_edges = []
edge = v_d + 1
for k in range(re.V_bin):
    m_lo = edge
    m_hi = max(int(edge * 1.3), edge + 1)
    m_hi = min(m_hi, V + 1)
    v_edges.append((m_lo, m_hi))
    edge = m_hi

for k in range(re.V_bin):
    mu0 = y[vac_start + PM*k, -1]
    mu1 = y[vac_start + PM*k + 1, -1] if PM >= 2 else 0
    n_bar = mu1 / max(abs(mu0), 1e-300)
    m_lo, m_hi = v_edges[k]
    bw = m_hi - m_lo

    S1 = sum(range(m_lo, m_hi))
    S2 = sum(m*m for m in range(m_lo, m_hi))
    det = bw * S2 - S1 * S1
    c_lo = c_hi = 0
    if abs(det) > 1e-30:
        phi0_lo = (S2 - S1 * m_lo) / det
        phi1_lo = (bw * m_lo - S1) / det
        c_lo = phi0_lo * mu0 + phi1_lo * mu1
        phi0_hi = (S2 - S1 * (m_hi-1)) / det
        phi1_hi = (bw * (m_hi-1) - S1) / det
        c_hi = phi0_hi * mu0 + phi1_hi * mu1

    flag = ""
    if mu0 < 0: flag += " MU0_NEG"
    if mu1 < 0 and mu0 > 1e-30: flag += " MU1_NEG"
    if c_lo < -1e-30: flag += f" c_lo={c_lo:.1e}_NEG"
    if c_hi < -1e-30: flag += f" c_hi={c_hi:.1e}_NEG"
    if mu0 > 1e-30 and (n_bar < m_lo * 0.8 or n_bar > m_hi * 1.2):
        flag += f" NBAR_OUT"

    print(f"  bin {k:2d}: [{m_lo:5d},{m_hi:5d}) bw={bw:4d}  "
          f"mu0={mu0:+10.2e}  mu1={mu1:+10.2e}  nbar={n_bar:8.1f}  "
          f"c_lo={c_lo:+.2e}  c_hi={c_hi:+.2e}{flag}", flush=True)

# Check discrete concentrations
sia_neg = np.sum(y[:i_d, -1] < -1e-30)
vac_neg = np.sum(y[i_VAC:i_VAC+v_d, -1] < -1e-30)
print(f"\nNegative discrete: SIA={sia_neg}/{i_d}, VAC={vac_neg}/{v_d}", flush=True)

# Find the discrete SIA concentration at the boundary (n=i_discrete)
c_boundary = y[i_d - 1, -1]
print(f"c_SIA[{i_d}] = {c_boundary:.3e} (discrete/bin boundary)", flush=True)

print(f"\ndelta_FP = {results['delta_FP'][-1]:.3e}", flush=True)
