#!/usr/bin/env python3
"""
Growth rate diagnostic — examines WHY the SIA peak doesn't move.

At ~0.1 dpa, extracts the SIA distribution and computes:
  - Net growth rate v(n) = K_grow(n)*c_1 - K_shrink(n)*c_v1  [size increments per second]
  - Production contribution at each n
  - Sink loss at each n
  - Monomer SIA and vacancy concentrations

This tells us whether loops are growing (v > 0) or shrinking (v < 0) at each size.
"""

import sys, os, io
from pathlib import Path
import numpy as np

MODULE_ROOT = Path(__file__).resolve().parent.parent.parent
REPO_ROOT   = MODULE_ROOT.parent
for p in [str(REPO_ROOT), str(MODULE_ROOT)]:
    if p not in sys.path:
        sys.path.insert(0, p)

import importlib
import RadCluster_1_0.py_utils.defect_production as _dp_mod
import RadCluster_1_0.py_utils.binding_energies as _be_mod
import RadCluster_1_0.py_utils.bin_moment_rates as _bmr
import RadCluster_1_0.py_utils.input_data as _inp_mod
import RadCluster_1_0.py_utils.reaction_rates as _rr_mod
import RadCluster_1_0.py_utils.rate_equations as _re_mod
import RadCluster_1_0.py_utils.cpp_bridge as _cb_mod
import RadCluster_1_0.py_utils.post_process as _pp_mod
import RadCluster_1_0.py_utils.simulation as _sim_mod
for _m in [_dp_mod, _be_mod, _bmr, _inp_mod, _rr_mod, _re_mod,
           _cb_mod, _pp_mod, _sim_mod]:
    importlib.reload(_m)
from RadCluster_1_0.py_utils.simulation import RadClusterSimulation

Omega  = 1.18e-29
b_111  = 2.482e-10

def n_to_d_nm(n):
    return 2.0 * np.sqrt(n * Omega / (np.pi * b_111)) * 1e9

# ── Run baseline simulation to 0.1 dpa ─────────────────────────────────────
SOLVER_CONFIG = {
    't_span':   (1e-6, 1e5),
    'n_points': 100,
    'log_time': True,
    'rtol':     1e-5,
    'atol':     1e-20,
    'solver_method': {
        'linsol': 'gmres',
        'window_w0_i': 50, 'window_width': 150,
        'window_C_expand': 1e-18, 'window_expand_pad': 10,
        'window_prec': 1,
        'window_gmres_maxl': 20, 'window_N_thresh': 500,
    }
}

I_SIM, V_SIM = 3000, 3000
i_mobile = 10

OVERRIDES = {
    'eta': 0.3, 'f_cl_i': 0.3, 'f_cl_v': 0.25,
    'E_m_1D': 0.4, 'i_mobile': i_mobile, 'L_hat': 71.8,
    'c_C': 1.94e-4, 'E_b_C_SIA': 0.65,
    'rho_d': 1e13, 'Z_i': 1.05, 'Z_ii': 1.2,
    'shape_function': 'linear',
}

print('Running baseline simulation to 0.1 dpa...')
saved = sys.stdout, sys.stderr
sys.stdout = sys.stderr = io.StringIO()
try:
    sim = RadClusterSimulation(
        I=I_SIM, V=V_SIM,
        solver_mode='active_window', physics_option='full_CD_fission',
        C_floor=1e-25, he_kinetics='quasi_steady_state',
        i_mobile=i_mobile, v_mobile=5,
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
    inp.diffusion['i_mobile'] = i_mobile
    inp.reactions['i_mobile'] = i_mobile
    inp.reactions['i_discrete'] = i_mobile
    inp.reactions['v_discrete'] = 5
    inp.reactions['I_bin'] = 0
    inp.reactions['V_bin'] = 0
    inp._calculate_derived()
    sim.rebuild_rates()
    results = sim.run_adaptive(
        solver_config=SOLVER_CONFIG, save_output=False,
        progress_callback=None, boundary_threshold=0.1,
        max_doublings=0, points_per_segment=10,
    )
finally:
    sys.stdout, sys.stderr = saved

if results is None:
    print('Simulation failed!')
    sys.exit(1)

# ── Extract state at 0.1 dpa ───────────────────────────────────────────────
d = inp.derived
G = d['G']
dose = results['dose']
idx_01 = np.argmin(np.abs(dose - 0.1))
actual_dose = dose[idx_01]
t_01 = results['t'][idx_01]

y = results['y'][:, idx_01]
c_i = np.maximum(y[0:I_SIM], 0.0)  # SIA concentrations [at.frac]
c_v = np.maximum(y[I_SIM:2*I_SIM], 0.0)  # VAC concentrations [at.frac]

ci1 = c_i[0]    # SIA monomer
cv1 = c_v[0]    # vacancy monomer

s2m = 1.0 / d['Omega']
print(f'\n{"="*70}')
print(f'  GROWTH RATE DIAGNOSTIC at dose = {actual_dose:.4f} dpa (t = {t_01:.2e} s)')
print(f'{"="*70}')
print(f'  C_i1 (SIA monomer) = {ci1:.4e} at.frac = {ci1*s2m:.4e} m^-3')
print(f'  C_v1 (vac monomer) = {cv1:.4e} at.frac = {cv1*s2m:.4e} m^-3')
print(f'  C_i1/C_v1 ratio    = {ci1/max(cv1, 1e-300):.3f}')

rr = sim.reaction_rates
re = sim.rate_equations

# ── Compute growth rates for each size n ───────────────────────────────────
# At each size n, the net flux is:
#   Growth: K_grow(n) * c_1 * c_n  (SIA monomer absorption makes n → n+1)
#   Shrink: K_shrink(n) * c_v1 * c_n  (vacancy monomer makes n → n-1)
#   Net growth velocity: v(n) = K_grow(n) * c_1 - K_shrink(n) * c_v1
#   This has units of [s^-1] — rate of concentration being transported

print(f'\n{"n":>5} {"d(nm)":>7} {"c_n(m-3)":>11} {"K_grow":>10} {"K_shrink":>10} '
      f'{"v_net":>11} {"G_emit":>10} {"k2_sink":>10} {"dc/dt ratio":>11}')
print('-' * 100)

for n in [2, 3, 4, 5, 8, 10, 11, 12, 15, 17, 20, 30, 50, 100, 200, 500]:
    if n > I_SIM:
        continue
    cn = c_i[n-1]
    K_g = rr.K_SIA_grow[n-1]    # growth rate constant (monomer absorption)
    K_s = rr.K_SIA_shrink[n-1]  # shrinkage rate constant (vacancy absorption)
    G_em = rr.G_SIA[n-1]        # thermal emission rate
    k2 = rr.k2_SIA[n-1]         # fixed sink rate

    v_net = K_g * ci1 - K_s * cv1  # net growth velocity [s^-1]
    # dc_n/dt contributions (just the growth/shrink part)
    growth_flux = K_g * ci1 * cn   # rate of c_n being removed by growth
    shrink_flux = K_s * cv1 * cn   # rate of c_n being removed by shrinkage
    ratio = growth_flux / max(shrink_flux, 1e-300) if shrink_flux > 0 else float('inf')

    print(f'{n:5d} {n_to_d_nm(n):7.2f} {cn*s2m:11.3e} {K_g:10.3e} {K_s:10.3e} '
          f'{v_net:11.3e} {G_em:10.3e} {k2:10.3e} {ratio:11.3f}')

# ── Overall flux balance for SIA monomers ──────────────────────────────────
print(f'\n{"="*70}')
print(f'  SIA MONOMER BALANCE (n=1)')
print(f'{"="*70}')

# Production: P_1 = eta * G * (1 - f_cl_i)
Pr_1 = re.Pr_SIA[0]
print(f'  Production P_1:           {Pr_1:.4e} at.frac/s')

# Sink loss to fixed sinks
sink_loss = rr.k2_SIA[0] * ci1
print(f'  Fixed sink loss (k2*c1):  {sink_loss:.4e} at.frac/s')

# Loss to cluster growth (SIA monomers absorbed by all loops)
growth_loss = ci1 * np.dot(rr.K_SIA_grow, c_i)
print(f'  Loop growth loss:         {growth_loss:.4e} at.frac/s')

# Recombination with vacancies
recom_loss = rr.K_iv * ci1 * cv1
print(f'  V-I recombination:        {recom_loss:.4e} at.frac/s')

# Cavity absorption
K_3D = rr.K_3D_cav_pref
m13 = np.arange(1, V_SIM+1)**(1./3.)
cav_loss = ci1 * K_3D * np.dot(m13, c_v)
print(f'  Cavity absorption:        {cav_loss:.4e} at.frac/s')

total_loss = sink_loss + growth_loss + recom_loss + cav_loss
print(f'  ────────────────────')
print(f'  Total loss:               {total_loss:.4e} at.frac/s')
print(f'  Production/Loss ratio:    {Pr_1/max(total_loss,1e-300):.4f}')
print(f'  Fixed sinks fraction:     {sink_loss/max(total_loss,1e-300)*100:.1f}%')
print(f'  Loop growth fraction:     {growth_loss/max(total_loss,1e-300)*100:.1f}%')
print(f'  Recombination fraction:   {recom_loss/max(total_loss,1e-300)*100:.1f}%')
print(f'  Cavity fraction:          {cav_loss/max(total_loss,1e-300)*100:.1f}%')

# ── Growth front velocity ──────────────────────────────────────────────────
print(f'\n{"="*70}')
print(f'  GROWTH FRONT ANALYSIS')
print(f'{"="*70}')
v_growth = np.array([rr.K_SIA_grow[n-1]*ci1 - rr.K_SIA_shrink[n-1]*cv1 for n in range(1, I_SIM+1)])
print(f'  v_net(n=10) = {v_growth[9]:.4e} s^-1 (growth velocity at mobility cutoff)')
print(f'  v_net(n=50) = {v_growth[49]:.4e} s^-1')
print(f'  v_net(n=100) = {v_growth[99]:.4e} s^-1')
print(f'  v_net(n=500) = {v_growth[499]:.4e} s^-1')
print(f'  At t={t_01:.0e}s, if uniform growth at v_net(100):')
print(f'    distance = v_net * t ≈ {v_growth[99]*t_01:.1f} size increments')
print(f'  → This means the growth front reaches n ≈ {v_growth[99]*t_01:.0f} by 0.1 dpa')
print(f'  → Diameter = {n_to_d_nm(max(1,v_growth[99]*t_01)):.2f} nm')

# ── Distribution shape ─────────────────────────────────────────────────────
print(f'\n{"="*70}')
print(f'  SIA SIZE DISTRIBUTION (top 20 sizes by concentration)')
print(f'{"="*70}')
sorted_idx = np.argsort(c_i)[::-1]
for rank in range(20):
    n = sorted_idx[rank] + 1
    cn = c_i[sorted_idx[rank]]
    print(f'  #{rank+1:2d}  n={n:5d}  d={n_to_d_nm(n):6.2f} nm  c_n={cn*s2m:.3e} m^-3')
