#!/usr/bin/env python3
"""
Parameter sweep for RadCluster_1_0
========================================
Objective: find input parameter sets that produce SIA loop size distribution
peak at 2-5 nm diameter by 0.1 dpa.

Diameter formula:  d_nm = 2 * sqrt(n * Omega / (pi * b_111)) * 1e9
  -> d=2 nm => n ~ 66,  d=5 nm => n ~ 413,  d=3 nm => n ~ 149

Strategy: vary the most influential parameters within physical ranges,
run the C++ active_window solver, extract the SIA peak location at
the time step closest to 0.1 dpa.
"""

import sys, os, io, json, time
from pathlib import Path
import numpy as np

# ── Path setup ──────────────────────────────────────────────────────────────
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

# ── Physical constants for diameter conversion ──────────────────────────────
Omega  = 1.18e-29   # m^3  (atomic volume bcc Fe)
b_111  = 2.482e-10  # m    (Burgers vector 1/2<111>)

def n_to_diameter_nm(n):
    """SIA cluster size -> loop diameter in nm."""
    return 2.0 * np.sqrt(n * Omega / (np.pi * b_111)) * 1e9

def diameter_nm_to_n(d_nm):
    """Loop diameter in nm -> SIA cluster size n."""
    d_m = d_nm * 1e-9
    return np.pi * b_111 * (d_m / 2.0)**2 / Omega


# ── Solver config (extend to 0.1 dpa: t_end = 1e5 s at G=1e-6 dpa/s) ─────
# Reduced domain (I=V=3000) and output points for faster sweep
I_SWEEP = 3000   # max SIA cluster size (d=5nm -> n~413, so 3000 is plenty)
V_SWEEP = 3000   # max vacancy cluster size

SOLVER_CONFIG = {
    't_span':   (1e-6, 1e5),   # reach 0.1 dpa
    'n_points': 100,           # fewer points for speed
    'log_time': True,
    'rtol':     1e-5,          # slightly relaxed for speed
    'atol':     1e-20,
    'solver_method': {
        'linsol':            'gmres',
        'window_width':      50,
        'concentration_threshold': 1e-18,
        'window_pad':              10,
    }
}

# ── Base parameters (notebook defaults) ─────────────────────────────────────
BASE_PARAMS = {
    'eta':       0.3,
    'f_cl_i':    0.3,
    'f_cl_v':    0.25,
    'E_m_1D':    0.4,
    'i_mobile':  10,
    'L_hat':     71.8,
    'c_C':       1.94e-4,
    'E_b_C_SIA': 0.65,
    'rho_d':     1e13,
    'Z_i':       1.05,
    'Z_ii':      1.2,
    'shape_function': 'linear',
}

# ── Parameter sets to test ──────────────────────────────────────────────────
# Physics rationale for each run:
#   - More coalescence (Z_ii, i_mobile) → larger loops
#   - Higher clustering (f_cl_i) → more SIA in clusters → faster growth
#   - Lower dislocation density (rho_d) → less sink loss → more growth
#   - Longer mean free path (L_hat) → less 1D sink capture → more coalescence
#   - Lower 1D migration barrier (E_m_1D) → faster 1D → more encounters
#   - Less C trapping (lower c_C or E_b_C_SIA) → faster effective diffusion
#   - Higher temperature → more thermal emission from small, growth of large

RUNS = {
    # ── Run 0: Current baseline (L_hat=71.8, rot=13547) ─────────────────
    'R00_baseline': {},

    # ── Runs 1-5: L_hat scan (key lever for 1D→3D mobility) ─────────────
    # Lower L_hat → more direction changes → higher effective 3D diffusivity
    # for glissile clusters n=4..i_mobile → more coalescence
    'R01_Lhat_30': {'L_hat': 30.0},        # rot=2363
    'R02_Lhat_15': {'L_hat': 15.0},        # rot=592
    'R03_Lhat_8':  {'L_hat': 8.0},         # rot=169
    'R04_Lhat_4':  {'L_hat': 4.0},         # rot=43
    'R05_Lhat_2':  {'L_hat': 2.0},         # rot=11.5

    # ── Runs 6-8: L_hat + coalescence + production ──────────────────────
    'R06_Lhat8_Zii2': {
        'L_hat': 8.0,
        'Z_ii': 2.0,
        'f_cl_i': 0.50,
    },

    'R07_Lhat4_full': {
        'L_hat': 4.0,
        'Z_ii': 1.5,
        'f_cl_i': 0.50,
        'rho_d': 5e12,
    },

    'R08_Lhat2_full': {
        'L_hat': 2.0,
        'Z_ii': 1.5,
        'f_cl_i': 0.50,
        'rho_d': 5e12,
    },
}


def run_single(name, overrides, I=I_SWEEP, V=V_SWEEP):
    """Run one simulation and return the SIA peak info at ~0.1 dpa."""
    params = dict(BASE_PARAMS)
    params.update(overrides)

    i_mob = int(params.get('i_mobile', 10))
    v_mob = 5

    # Patch module-level FISSION dict for production parameters
    # (production_rates() reads from this, not from inp.production_fission)
    _fission_backup = dict(_FISSION_DICT)
    for key in ['eta', 'f_cl_i', 'f_cl_v', 's_i', 's_v', 'i_cascade', 'v_cascade']:
        if key in params:
            _FISSION_DICT[key] = params[key]

    # Suppress output
    saved = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=I, V=V,
            solver_mode='active_window',
            physics_option='full_CD_fission',
            C_floor=1e-25,
            he_kinetics='quasi_steady_state',
            i_mobile=i_mob,
            v_mobile=v_mob,
        )

        # Apply overrides
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
        # Discrete/bin settings (full_CD)
        inp.reactions['i_discrete'] = i_mob
        inp.reactions['v_discrete'] = v_mob
        inp.reactions['I_bin'] = 0
        inp.reactions['V_bin'] = 0
        inp._calculate_derived()
        sim.rebuild_rates()

        results = sim.run_adaptive(
            solver_config=SOLVER_CONFIG,
            save_output=False,
            progress_callback=None,
            boundary_threshold=0.1,
            max_doublings=0,
            points_per_segment=10,
        )
    finally:
        sys.stdout, sys.stderr = saved
        # Restore FISSION dict for next run
        _FISSION_DICT.update(_fission_backup)

    if results is None:
        return None

    # ── Extract SIA peak at closest time to 0.1 dpa ──────────────────────
    dose = results['dose']
    t_arr = results['t']
    G_dpa = inp.derived['G']
    Omega_val = inp.derived['Omega']
    s2m = 1.0 / Omega_val

    # Find index closest to 0.1 dpa
    target_dose = 0.1
    idx_01 = np.argmin(np.abs(dose - target_dose))
    actual_dose = dose[idx_01]

    # Get SIA size distribution at that time
    # results['y'] has shape [N_eq, n_times]; SIA clusters are y[0:I, :]
    y_full = results['y']           # atom fractions [N_eq, n_t]
    c_i_at_01 = np.maximum(y_full[0:I, idx_01], 0.0) * s2m  # m^-3

    N_sia = c_i_at_01.shape[0]
    ns = np.arange(1, N_sia + 1)
    d_nm = n_to_diameter_nm(ns)

    # Find peak: weighted by concentration (exclude n=1 mono-interstitials)
    c_trimmed = c_i_at_01[1:]  # n >= 2
    ns_trimmed = ns[1:]
    d_trimmed = d_nm[1:]

    if np.max(c_trimmed) <= 0:
        return None

    # Peak = argmax of concentration for n >= 2
    ipeak = np.argmax(c_trimmed)
    peak_n = ns_trimmed[ipeak]
    peak_d = d_trimmed[ipeak]
    peak_c = c_trimmed[ipeak]

    # Also compute number-weighted mean diameter (for clusters n >= 2 with c > threshold)
    mask = c_trimmed > 1e10  # above noise floor
    if mask.any():
        mean_d = np.sum(d_trimmed[mask] * c_trimmed[mask]) / np.sum(c_trimmed[mask])
        total_N = np.sum(c_trimmed[mask])
    else:
        mean_d = 0.0
        total_N = 0.0

    # Find the diameter range containing 50% of the concentration (FWHM-like)
    # Sort by concentration to find the dominant sizes
    sorted_idx = np.argsort(c_trimmed)[::-1]
    cum_c = np.cumsum(c_trimmed[sorted_idx])
    half_total = cum_c[-1] * 0.5
    top_half_mask = cum_c <= half_total
    if top_half_mask.any():
        top_half_n = ns_trimmed[sorted_idx[top_half_mask]]
        d_range = (n_to_diameter_nm(top_half_n.min()), n_to_diameter_nm(top_half_n.max()))
    else:
        d_range = (peak_d, peak_d)

    return {
        'name':         name,
        'dose_actual':  actual_dose,
        'peak_n':       int(peak_n),
        'peak_d_nm':    peak_d,
        'peak_c_m3':    peak_c,
        'mean_d_nm':    mean_d,
        'total_N_m3':   total_N,
        'd_range_nm':   d_range,
        'mean_n_i':     results['mean_n_i'][idx_01],
        'N_loops':      results['N_loops'][idx_01],
        'swelling_pct': results['swelling'][idx_01] * 100,
        'delta_FP':     results['delta_FP'][idx_01],
        'params':       params,
    }


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print('=' * 80)
    print('  PARAMETER SWEEP — SIA Loop Size Peak Optimization')
    print('  Target: peak diameter 2-5 nm at 0.1 dpa')
    print(f'  (n_target: {diameter_nm_to_n(2.0):.0f}-{diameter_nm_to_n(5.0):.0f})')
    print('=' * 80)

    all_results = []
    t0_total = time.time()

    for i, (name, overrides) in enumerate(RUNS.items()):
        print(f'\n[{i+1}/{len(RUNS)}] Running {name} ...', flush=True)
        if overrides:
            for k, v in overrides.items():
                print(f'    {k} = {v}')
        else:
            print('    (baseline — no overrides)')

        t0 = time.time()
        try:
            res = run_single(name, overrides)
        except Exception as e:
            print(f'    FAILED: {e}')
            res = None
        elapsed = time.time() - t0

        if res is not None:
            all_results.append(res)
            in_target = 2.0 <= res['peak_d_nm'] <= 5.0
            marker = ' *** IN TARGET ***' if in_target else ''
            print(f'    Done in {elapsed:.1f}s  |  dose={res["dose_actual"]:.4f} dpa')
            print(f'    Peak: n={res["peak_n"]}, d={res["peak_d_nm"]:.2f} nm, '
                  f'C={res["peak_c_m3"]:.2e} m^-3{marker}')
            print(f'    Mean d={res["mean_d_nm"]:.2f} nm, '
                  f'N_loops={res["N_loops"]:.2e} m^-3, '
                  f'swelling={res["swelling_pct"]:.4f}%')
        else:
            print(f'    FAILED ({elapsed:.1f}s)')

    total_time = time.time() - t0_total

    # ── Summary table ─────────────────────────────────────────────────────
    print('\n' + '=' * 80)
    print('  SUMMARY TABLE — sorted by distance to target (3.5 nm)')
    print('=' * 80)
    target_d = 3.5  # nm (center of 2-5 nm range)
    all_results.sort(key=lambda r: abs(r['peak_d_nm'] - target_d))

    print(f'{"Run":<25} {"Peak d (nm)":>11} {"Peak n":>7} {"Mean d (nm)":>11} '
          f'{"N_loops":>10} {"Swell %":>9} {"In target":>10}')
    print('-' * 85)
    for r in all_results:
        in_t = 'YES' if 2.0 <= r['peak_d_nm'] <= 5.0 else 'no'
        print(f'{r["name"]:<25} {r["peak_d_nm"]:>11.2f} {r["peak_n"]:>7d} '
              f'{r["mean_d_nm"]:>11.2f} {r["N_loops"]:>10.2e} '
              f'{r["swelling_pct"]:>9.4f} {in_t:>10}')

    # ── Best result ───────────────────────────────────────────────────────
    print('\n' + '=' * 80)
    if all_results:
        best = all_results[0]
        print(f'  BEST RUN: {best["name"]}')
        print(f'  Peak diameter: {best["peak_d_nm"]:.2f} nm  (target: 2-5 nm)')
        print(f'  Peak cluster size n = {best["peak_n"]}')
        print(f'  Mean diameter: {best["mean_d_nm"]:.2f} nm')
        print(f'  Loop density: {best["N_loops"]:.2e} m^-3')
        print(f'  Swelling: {best["swelling_pct"]:.4f} %')
        print(f'  delta_FP: {best["delta_FP"]:.2e}')
        print(f'\n  Optimal parameter overrides:')
        # Only print non-base values
        for k, v in best['params'].items():
            base_v = BASE_PARAMS.get(k)
            if base_v != v:
                print(f'    {k:>15} = {v}  (base: {base_v})')
            else:
                print(f'    {k:>15} = {v}')
    else:
        print('  No successful runs!')

    print(f'\n  Total sweep time: {total_time:.0f}s ({total_time/60:.1f} min)')
    print('=' * 80)

    # ── Save results to JSON ──────────────────────────────────────────────
    out_path = Path(__file__).parent / 'param_sweep_results.json'
    serializable = []
    for r in all_results:
        sr = dict(r)
        sr['d_range_nm'] = list(sr['d_range_nm'])
        # Convert numpy types
        for k, v in sr.items():
            if isinstance(v, (np.integer,)):
                sr[k] = int(v)
            elif isinstance(v, (np.floating,)):
                sr[k] = float(v)
            elif isinstance(v, dict):
                sr[k] = {kk: (float(vv) if isinstance(vv, (np.floating,)) else
                              int(vv) if isinstance(vv, (np.integer,)) else vv)
                          for kk, vv in v.items()}
        serializable.append(sr)
    with open(out_path, 'w') as f:
        json.dump(serializable, f, indent=2)
    print(f'\nResults saved to {out_path}')
