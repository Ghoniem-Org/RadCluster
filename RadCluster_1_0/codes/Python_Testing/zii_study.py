#!/usr/bin/env python3
"""
Z_ii sweep — the SIA-SIA coalescence bias factor.

Z_ii multiplies the full K_ii_coal rate (SIA absorption by SIA loops)
but does NOT affect vacancy absorption. This makes it the primary knob
for shifting the growth/shrink balance.

With corrected loop geometry (A_loop * n^{1/2}):
  growth_ratio = Z_ii * Z_i_loop * Di * c_i1 / (Dv * c_v1)

Baseline ratio ≈ 0.742 at Z_ii=1.2.  Need ratio > 1.0 for net growth.
Estimated crossover: Z_ii ≈ 1.6-1.7.
"""

import sys, os, io, time
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
from RadCluster_1_0.py_utils.defect_production import FISSION as _FISSION_DICT

Omega  = 1.18e-29
b_111  = 2.482e-10
def n_to_d_nm(n):
    return 2.0 * np.sqrt(np.maximum(n, 1) * Omega / (np.pi * b_111)) * 1e9

I_SIM, V_SIM = 3000, 3000
SOLVER_CONFIG = {
    't_span':   (1e-6, 1e5),
    'n_points': 150,
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

BASE = {
    'eta': 0.30, 'f_cl_i': 0.30, 'f_cl_v': 0.25,
    'E_m_1D': 0.4, 'i_mobile': 10, 'L_hat': 71.8,
    'c_C': 1.94e-4, 'E_b_C_SIA': 0.65,
    'rho_d': 1e13, 'Z_i': 1.05,
    'shape_function': 'linear',
}


def run_one(name, Z_ii_val):
    params = dict(BASE)
    params['Z_ii'] = Z_ii_val

    fission_bak = dict(_FISSION_DICT)
    for key in ['eta', 'f_cl_i', 'f_cl_v']:
        if key in params:
            _FISSION_DICT[key] = params[key]

    saved = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=I_SIM, V=V_SIM,
            solver_mode='active_window', physics_option='full_CD_fission',
            C_floor=1e-25, he_kinetics='quasi_steady_state',
            i_mobile=10, v_mobile=5,
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
        inp.diffusion['i_mobile'] = 10
        inp.reactions['i_mobile'] = 10
        inp.reactions['i_discrete'] = 10
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
        _FISSION_DICT.update(fission_bak)

    if results is None:
        return None

    d = inp.derived
    rr = sim.reaction_rates
    s2m = 1.0 / d['Omega']
    y = results['y']
    dose = results['dose']
    idx = np.argmin(np.abs(dose - 0.1))

    ci1 = np.maximum(y[0, idx], 0.0)
    cv1 = np.maximum(y[I_SIM, idx], 0.0)
    c_i = np.maximum(y[0:I_SIM, idx], 0.0) * s2m
    c_v = np.maximum(y[I_SIM:2*I_SIM, idx], 0.0) * s2m

    # SIA peak (n >= 2)
    ipk = np.argmax(c_i[1:]) + 1
    sia_pk_n = ipk + 1
    sia_pk_d = n_to_d_nm(sia_pk_n)

    # Mean SIA diameter
    ns = np.arange(1, I_SIM + 1, dtype=float)
    mask2 = c_i[1:] > 1e5
    if mask2.any():
        mean_d = np.sum(n_to_d_nm(ns[1:])[mask2] * c_i[1:][mask2]) / np.sum(c_i[1:][mask2])
    else:
        mean_d = 0.0

    # VAC peak (m >= 2)
    vpk = np.argmax(c_v[1:]) + 1
    vac_pk_n = vpk + 1

    # Growth ratio at n=100
    K_g = rr.K_SIA_grow[99]
    K_s = rr.K_SIA_shrink[99]
    v_net = K_g * ci1 - K_s * cv1
    ratio = (K_g * ci1) / max(K_s * cv1, 1e-300)

    # Flux arrays vs dose
    omega_i = d['omega_i_eff']
    omega_v = d['omega_v_eff']
    flux_i = omega_i * np.maximum(y[0, :], 0.0) * s2m
    flux_v = omega_v * np.maximum(y[I_SIM, :], 0.0) * s2m

    return {
        'name': name, 'Z_ii': Z_ii_val,
        'dose': dose, 'flux_i': flux_i, 'flux_v': flux_v,
        'sia_pk_n': sia_pk_n, 'sia_pk_d': sia_pk_d,
        'sia_mean_d': mean_d,
        'vac_pk_n': vac_pk_n,
        'v_net': v_net, 'ratio': ratio,
        'N_loops': results['N_loops'][idx],
        'swelling': results['swelling'][idx] * 100,
        'c_i': c_i, 'c_v': c_v,
        'ci1_m3': ci1 * s2m, 'cv1_m3': cv1 * s2m,
    }


# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    out_dir = Path(__file__).parent / 'zii_study_output'
    out_dir.mkdir(exist_ok=True)

    Z_ii_values = [1.0, 1.2, 1.5, 1.7, 2.0, 2.5, 3.0]

    print('=' * 70)
    print('  Z_ii SWEEP — SIA-SIA Coalescence Bias')
    print('  Loop geometry corrected (A_loop * n^{1/2})')
    print('=' * 70)

    results = []
    t0 = time.time()
    for i, zii in enumerate(Z_ii_values):
        name = f'Zii_{zii:.1f}'
        print(f'\n[{i+1}/{len(Z_ii_values)}] Z_ii = {zii}', flush=True)
        t1 = time.time()
        try:
            r = run_one(name, zii)
        except Exception as e:
            print(f'    FAILED: {e}')
            import traceback; traceback.print_exc()
            r = None
        dt = time.time() - t1
        if r:
            results.append(r)
            marker = ' *** GROWING ***' if r['ratio'] > 1.0 else ''
            print(f'    Done ({dt:.0f}s) | peak: n={r["sia_pk_n"]}, d={r["sia_pk_d"]:.2f} nm | '
                  f'mean_d={r["sia_mean_d"]:.2f} nm | ratio={r["ratio"]:.3f}{marker}')

    total = time.time() - t0
    print(f'\nTotal: {total:.0f}s ({total/60:.1f} min)')

    # ── Plot 1: Flux ratio vs dose for all Z_ii ─────────────────────────
    fig1, ax1 = plt.subplots(figsize=(10, 5))
    for r in results:
        dose = r['dose']
        mask = dose > 1e-8
        rat = r['flux_i'][mask] / np.maximum(r['flux_v'][mask], 1e-30)
        ax1.semilogx(dose[mask], rat, label=f'Z_ii={r["Z_ii"]:.1f}', linewidth=1.5)
    ax1.axhline(y=1.0, color='red', ls='--', lw=2, label='ratio=1 (balance)')
    ax1.set_xlabel('Dose [dpa]')
    ax1.set_ylabel(r'$\omega_i C_{i1}\;/\;\omega_v C_{v1}$')
    ax1.set_title(r'Monomer flux ratio vs $Z_{ii}$ — above 1 = net loop growth')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, max(3.0, ax1.get_ylim()[1]))
    fig1.tight_layout()
    fig1.savefig(str(out_dir / 'flux_ratio_vs_Zii.png'), dpi=150)
    print(f'\nSaved: {out_dir / "flux_ratio_vs_Zii.png"}')

    # ── Plot 2: SIA size distributions at 0.1 dpa ───────────────────────
    fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(14, 6))
    ns = np.arange(1, I_SIM + 1)
    d_nm = n_to_d_nm(ns)
    for r in results:
        c_i = r['c_i']
        m = c_i > 1e5
        ax2a.semilogy(d_nm[m], c_i[m], label=f'Z_ii={r["Z_ii"]:.1f}', lw=1.3)
        ax2b.plot(d_nm[m], c_i[m], label=f'Z_ii={r["Z_ii"]:.1f}', lw=1.3)
    for ax in [ax2a, ax2b]:
        ax.set_xlabel('SIA loop diameter [nm]')
        ax.set_ylabel('Concentration [m$^{-3}$]')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    ax2a.set_title('SIA distribution at 0.1 dpa (log)')
    ax2a.set_xlim(0, 15)
    ax2b.set_title('SIA distribution at 0.1 dpa (linear)')
    ax2b.set_xlim(0, 10)
    fig2.tight_layout()
    fig2.savefig(str(out_dir / 'sia_dist_vs_Zii.png'), dpi=150)
    print(f'Saved: {out_dir / "sia_dist_vs_Zii.png"}')

    # ── Plot 3: Growth ratio and SIA peak vs Z_ii ────────────────────────
    fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(12, 5))
    zvals = [r['Z_ii'] for r in results]
    ratios = [r['ratio'] for r in results]
    peaks = [r['sia_pk_d'] for r in results]
    means = [r['sia_mean_d'] for r in results]

    ax3a.plot(zvals, ratios, 'bo-', markersize=8, linewidth=2)
    ax3a.axhline(y=1.0, color='red', ls='--', lw=2, label='ratio=1 (balance)')
    ax3a.set_xlabel(r'$Z_{ii}$')
    ax3a.set_ylabel('Growth/Shrink ratio at n=100')
    ax3a.set_title('Growth ratio vs Z_ii')
    ax3a.legend()
    ax3a.grid(True, alpha=0.3)

    ax3b.plot(zvals, peaks, 'rs-', markersize=8, linewidth=2, label='Peak diameter')
    ax3b.plot(zvals, means, 'g^-', markersize=8, linewidth=2, label='Mean diameter')
    ax3b.axhspan(2, 5, alpha=0.15, color='green', label='Target: 2-5 nm')
    ax3b.set_xlabel(r'$Z_{ii}$')
    ax3b.set_ylabel('Diameter [nm]')
    ax3b.set_title('SIA loop size vs Z_ii at 0.1 dpa')
    ax3b.legend()
    ax3b.grid(True, alpha=0.3)

    fig3.tight_layout()
    fig3.savefig(str(out_dir / 'ratio_and_peak_vs_Zii.png'), dpi=150)
    print(f'Saved: {out_dir / "ratio_and_peak_vs_Zii.png"}')

    # ── Summary table ────────────────────────────────────────────────────
    print('\n' + '=' * 95)
    print(f'{"Z_ii":>6} {"SIA pk n":>9} {"SIA pk d":>9} {"Mean d":>8} '
          f'{"VAC pk n":>9} {"ratio":>7} {"v_net(100)":>11} '
          f'{"N_loops":>10} {"Swell%":>8} {"Status":>10}')
    print('-' * 95)
    for r in results:
        status = 'GROWING' if r['ratio'] > 1.0 else 'shrinking'
        print(f'{r["Z_ii"]:>6.1f} {r["sia_pk_n"]:>9d} {r["sia_pk_d"]:>9.2f} '
              f'{r["sia_mean_d"]:>8.2f} {r["vac_pk_n"]:>9d} {r["ratio"]:>7.3f} '
              f'{r["v_net"]:>11.3e} {r["N_loops"]:>10.2e} {r["swelling"]:>8.4f} '
              f'{status:>10}')
    print('=' * 95)

    plt.close('all')
    print('\nDone.')
