#!/usr/bin/env python3
"""
Bias study — examines how dislocation bias (Z_i) and production bias (f_cl_i, f_cl_v)
affect SIA and vacancy cluster size peaks, and plots the fundamental monomer fluxes.

With the loop geometry correction (A_loop * n^{1/2} for n ≥ 4), the Z_i_loop bias
creates differential SIA vs vacancy absorption that drives net loop growth.
"""

import sys, os, io, time, json
from pathlib import Path
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
import RadCluster_2_0.py_utils.defect_production as _dp_mod
import RadCluster_2_0.py_utils.binding_energies as _be_mod
import RadCluster_2_0.py_utils.bin_moment_rates as _bmr
import RadCluster_2_0.py_utils.input_data as _inp_mod
import RadCluster_2_0.py_utils.reaction_rates as _rr_mod
import RadCluster_2_0.py_utils.rate_equations as _re_mod
import RadCluster_2_0.py_utils.cpp_bridge as _cb_mod
import RadCluster_2_0.py_utils.post_process as _pp_mod
import RadCluster_2_0.py_utils.simulation as _sim_mod
for _m in [_dp_mod, _be_mod, _bmr, _inp_mod, _rr_mod, _re_mod,
           _cb_mod, _pp_mod, _sim_mod]:
    importlib.reload(_m)
from RadCluster_2_0.py_utils.simulation import RadClusterSimulation
from RadCluster_2_0.py_utils.defect_production import FISSION as _FISSION_DICT

# ── Constants ───────────────────────────────────────────────────────────────
Omega  = 1.18e-29
b_111  = 2.482e-10
def n_to_d_nm(n):
    return 2.0 * np.sqrt(np.maximum(n, 1) * Omega / (np.pi * b_111)) * 1e9

# ── Solver config ───────────────────────────────────────────────────────────
I_SIM, V_SIM = 3000, 3000
SOLVER_CONFIG = {
    't_span':   (1e-6, 1e5),   # → 0.1 dpa at G=1e-6
    'n_points': 150,
    'log_time': True,
    'rtol':     1e-5,
    'atol':     1e-20,
    'solver_method': {
        'linsol': 'gmres',
        'window_width': 50,
        'concentration_threshold': 1e-18, 'window_pad': 10,
    }
}

BASE_PARAMS = {
    'E_m_1D': 0.4, 'i_mobile': 10, 'L_hat': 71.8,
    'c_C': 1.94e-4, 'E_b_C_SIA': 0.65,
    'rho_d': 1e13,
    'shape_function': 'linear',
}


def run_case(name, overrides):
    """Run one simulation and return full results + rate info."""
    params = dict(BASE_PARAMS)
    params.update(overrides)

    i_mob = int(params.get('i_mobile', 10))

    # Patch FISSION dict for production parameters
    fission_backup = dict(_FISSION_DICT)
    for key in ['eta', 'f_cl_i', 'f_cl_v', 's_i', 's_v', 'i_cascade', 'v_cascade']:
        if key in params:
            _FISSION_DICT[key] = params[key]

    saved = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=I_SIM, V=V_SIM,
            solver_mode='active_window', physics_option='full_CD_fission',
            C_floor=1e-25, he_kinetics='quasi_steady_state',
            i_mobile=i_mob, v_mobile=5,
        )
        inp = sim.input_data
        for key, val in params.items():
            placed = False
            for d in [inp.production_fission, inp.production_fusion,
                      inp.diffusion, inp.reactions,
                      inp.energetics, inp.dissociation]:
                if key in d:
                    d[key] = val
                    placed = True
            if not placed:
                inp.reactions[key] = val
        inp.diffusion['i_mobile'] = i_mob
        inp.reactions['i_mobile'] = i_mob
        inp.reactions['i_discrete'] = i_mob
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
        _FISSION_DICT.update(fission_backup)

    if results is None:
        return None

    d = inp.derived
    rr = sim.reaction_rates

    # Extract time-resolved monomer concentrations and fluxes
    y_all = results['y']
    n_t = len(results['t'])
    omega_i = d['omega_i_eff']
    omega_v = d['omega_v_eff']
    s2m = 1.0 / d['Omega']

    ci1_vs_t = np.maximum(y_all[0, :], 0.0)       # at.frac
    cv1_vs_t = np.maximum(y_all[I_SIM, :], 0.0)   # at.frac
    flux_i = omega_i * ci1_vs_t * s2m   # SIA monomer flux [m^-3 s^-1]
    flux_v = omega_v * cv1_vs_t * s2m   # VAC monomer flux [m^-3 s^-1]

    # SIA peak at 0.1 dpa
    dose = results['dose']
    idx_01 = np.argmin(np.abs(dose - 0.1))
    c_i_01 = np.maximum(y_all[0:I_SIM, idx_01], 0.0) * s2m
    c_v_01 = np.maximum(y_all[I_SIM:2*I_SIM, idx_01], 0.0) * s2m

    # SIA peak (exclude n=1)
    ipeak = np.argmax(c_i_01[1:]) + 1
    sia_peak_n = ipeak + 1
    sia_peak_d = n_to_d_nm(sia_peak_n)
    sia_peak_c = c_i_01[ipeak]
    sia_mean_d = n_to_d_nm(results['mean_n_i'][idx_01])

    # VAC peak (exclude m=1)
    vpeak = np.argmax(c_v_01[1:]) + 1
    vac_peak_n = vpeak + 1
    vac_peak_d = n_to_d_nm(vac_peak_n)

    # Growth rate diagnostic
    K_g = rr.K_SIA_grow
    K_s = rr.K_SIA_shrink
    ci1_01 = ci1_vs_t[idx_01]
    cv1_01 = cv1_vs_t[idx_01]
    v_net_100 = K_g[99] * ci1_01 - K_s[99] * cv1_01 if I_SIM > 100 else 0
    ratio = (K_g[99] * ci1_01) / max(K_s[99] * cv1_01, 1e-300) if I_SIM > 100 else 0

    return {
        'name': name,
        'params': params,
        'dose': dose,
        't': results['t'],
        'flux_i': flux_i,
        'flux_v': flux_v,
        'ci1_m3': ci1_vs_t * s2m,
        'cv1_m3': cv1_vs_t * s2m,
        'sia_peak_n': sia_peak_n,
        'sia_peak_d': sia_peak_d,
        'sia_peak_c': sia_peak_c,
        'sia_mean_d': sia_mean_d,
        'vac_peak_n': vac_peak_n,
        'vac_peak_d': vac_peak_d,
        'N_loops': results['N_loops'][idx_01],
        'N_voids': results['N_voids'][idx_01],
        'swelling': results['swelling'][idx_01] * 100,
        'v_net_100': v_net_100,
        'growth_ratio': ratio,
        'c_i_01': c_i_01,
        'c_v_01': c_v_01,
        'omega_i': omega_i,
        'omega_v': omega_v,
    }


# ── Case definitions ────────────────────────────────────────────────────────
# Explore: (1) dislocation bias Z_i, (2) production bias f_cl_i/f_cl_v
CASES = {
    # ── Baseline ─────────────────────────────────────────────────────────
    'baseline': {
        'eta': 0.30, 'f_cl_i': 0.30, 'f_cl_v': 0.25,
        'Z_i': 1.05, 'Z_ii': 1.2,
    },

    # ── Dislocation bias scan ────────────────────────────────────────────
    'Zi_1.10': {
        'eta': 0.30, 'f_cl_i': 0.30, 'f_cl_v': 0.25,
        'Z_i': 1.10, 'Z_ii': 1.2,
    },
    'Zi_1.15': {
        'eta': 0.30, 'f_cl_i': 0.30, 'f_cl_v': 0.25,
        'Z_i': 1.15, 'Z_ii': 1.2,
    },
    'Zi_1.20': {
        'eta': 0.30, 'f_cl_i': 0.30, 'f_cl_v': 0.25,
        'Z_i': 1.20, 'Z_ii': 1.2,
    },

    # ── Production bias scan (increase f_cl_i → more SIA in clusters → more free VAC)
    'fcli_0.50': {
        'eta': 0.30, 'f_cl_i': 0.50, 'f_cl_v': 0.25,
        'Z_i': 1.05, 'Z_ii': 1.2,
    },
    'fclv_0.40': {
        'eta': 0.30, 'f_cl_i': 0.30, 'f_cl_v': 0.40,
        'Z_i': 1.05, 'Z_ii': 1.2,
    },

    # ── Combined: strong dislocation bias + production bias ─────────────
    'Zi1.15_fcli0.50': {
        'eta': 0.30, 'f_cl_i': 0.50, 'f_cl_v': 0.25,
        'Z_i': 1.15, 'Z_ii': 1.2,
    },
    'Zi1.15_fclv0.40': {
        'eta': 0.30, 'f_cl_i': 0.30, 'f_cl_v': 0.40,
        'Z_i': 1.15, 'Z_ii': 1.2,
    },
}


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    out_dir = Path(__file__).parent / 'bias_study_output'
    out_dir.mkdir(exist_ok=True)

    print('=' * 70)
    print('  BIAS STUDY — Loop Geometry Corrected')
    print('  (A_loop * n^{1/2} * Z_i for SIA; A_loop * n^{1/2} for vacancy)')
    print('=' * 70)

    all_results = []
    t0 = time.time()

    for i, (name, overrides) in enumerate(CASES.items()):
        print(f'\n[{i+1}/{len(CASES)}] {name}', flush=True)
        for k, v in overrides.items():
            print(f'    {k} = {v}')
        t1 = time.time()
        try:
            res = run_case(name, overrides)
        except Exception as e:
            print(f'    FAILED: {e}')
            import traceback; traceback.print_exc()
            res = None
        elapsed = time.time() - t1

        if res:
            all_results.append(res)
            print(f'    Done ({elapsed:.0f}s) | SIA peak: n={res["sia_peak_n"]}, '
                  f'd={res["sia_peak_d"]:.2f} nm | VAC peak: n={res["vac_peak_n"]} | '
                  f'v_net(100)={res["v_net_100"]:.3e} | ratio={res["growth_ratio"]:.3f}')

    total = time.time() - t0
    print(f'\nTotal time: {total:.0f}s ({total/60:.1f} min)')

    # ═══════════════════════════════════════════════════════════════════════
    # PLOT 1: Fundamental monomer fluxes ω_i·c_i1 and ω_v·c_v1 vs dose
    # ═══════════════════════════════════════════════════════════════════════
    fig1, axes1 = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    for res in all_results:
        dose = res['dose']
        mask = dose > 1e-8
        axes1[0].loglog(dose[mask], res['flux_i'][mask], label=res['name'])
        axes1[1].loglog(dose[mask], res['flux_v'][mask], label=res['name'])

    axes1[0].set_ylabel(r'$\omega_i^{eff} \cdot C_{i1}$  [m$^{-3}$ s$^{-1}$]')
    axes1[0].set_title('SIA monomer flux (drives loop growth)')
    axes1[0].legend(fontsize=7, ncol=2)
    axes1[0].grid(True, alpha=0.3)

    axes1[1].set_ylabel(r'$\omega_v^{eff} \cdot C_{v1}$  [m$^{-3}$ s$^{-1}$]')
    axes1[1].set_xlabel('Dose [dpa]')
    axes1[1].set_title('Vacancy monomer flux (drives loop shrinkage)')
    axes1[1].legend(fontsize=7, ncol=2)
    axes1[1].grid(True, alpha=0.3)

    fig1.tight_layout()
    fig1.savefig(str(out_dir / 'monomer_fluxes.png'), dpi=150)
    print(f'\nSaved: {out_dir / "monomer_fluxes.png"}')

    # ═══════════════════════════════════════════════════════════════════════
    # PLOT 2: Flux ratio ω_i·c_i1 / ω_v·c_v1 vs dose
    # ═══════════════════════════════════════════════════════════════════════
    fig2, ax2 = plt.subplots(figsize=(10, 5))
    for res in all_results:
        dose = res['dose']
        mask = dose > 1e-8
        ratio = res['flux_i'][mask] / np.maximum(res['flux_v'][mask], 1e-30)
        ax2.semilogx(dose[mask], ratio, label=res['name'])
    ax2.axhline(y=1.0, color='red', linestyle='--', linewidth=1.5, label='ratio = 1 (growth/shrink balance)')
    ax2.set_xlabel('Dose [dpa]')
    ax2.set_ylabel(r'$\omega_i C_{i1}\;/\;\omega_v C_{v1}$')
    ax2.set_title('Monomer flux ratio — above 1 means net SIA loop growth')
    ax2.legend(fontsize=7, ncol=2)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, max(2.0, ax2.get_ylim()[1]))
    fig2.tight_layout()
    fig2.savefig(str(out_dir / 'flux_ratio.png'), dpi=150)
    print(f'Saved: {out_dir / "flux_ratio.png"}')

    # ═══════════════════════════════════════════════════════════════════════
    # PLOT 3: SIA size distributions at 0.1 dpa
    # ═══════════════════════════════════════════════════════════════════════
    fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(14, 6))
    ns = np.arange(1, I_SIM + 1)
    d_nm = n_to_d_nm(ns)

    for res in all_results:
        c_i = res['c_i_01']
        mask = c_i > 1e5
        ax3a.semilogy(d_nm[mask], c_i[mask], label=res['name'], linewidth=1.2)
        ax3b.plot(d_nm[mask], c_i[mask], label=res['name'], linewidth=1.2)

    ax3a.set_xlabel('SIA loop diameter [nm]')
    ax3a.set_ylabel('Concentration [m$^{-3}$]')
    ax3a.set_title('SIA size distribution at 0.1 dpa (log scale)')
    ax3a.set_xlim(0, 15)
    ax3a.legend(fontsize=7)
    ax3a.grid(True, alpha=0.3)

    ax3b.set_xlabel('SIA loop diameter [nm]')
    ax3b.set_ylabel('Concentration [m$^{-3}$]')
    ax3b.set_title('SIA size distribution at 0.1 dpa (linear scale)')
    ax3b.set_xlim(0, 10)
    ax3b.legend(fontsize=7)
    ax3b.grid(True, alpha=0.3)

    fig3.tight_layout()
    fig3.savefig(str(out_dir / 'sia_distributions.png'), dpi=150)
    print(f'Saved: {out_dir / "sia_distributions.png"}')

    # ═══════════════════════════════════════════════════════════════════════
    # Summary table
    # ═══════════════════════════════════════════════════════════════════════
    print('\n' + '=' * 100)
    print(f'{"Case":<22} {"Z_i":>5} {"f_cl_i":>7} {"f_cl_v":>7} '
          f'{"SIA pk d":>9} {"SIA pk n":>9} {"VAC pk n":>9} '
          f'{"v_net(100)":>11} {"ratio":>7} {"N_loops":>10} {"Swell%":>8}')
    print('-' * 100)
    for r in all_results:
        p = r['params']
        print(f'{r["name"]:<22} {p["Z_i"]:>5.2f} {p["f_cl_i"]:>7.2f} {p["f_cl_v"]:>7.2f} '
              f'{r["sia_peak_d"]:>9.2f} {r["sia_peak_n"]:>9d} {r["vac_peak_n"]:>9d} '
              f'{r["v_net_100"]:>11.3e} {r["growth_ratio"]:>7.3f} '
              f'{r["N_loops"]:>10.2e} {r["swelling"]:>8.4f}')
    print('=' * 100)

    plt.close('all')
    print('\nDone.')
