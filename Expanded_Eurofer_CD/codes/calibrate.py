#!/usr/bin/env python3
"""
calibrate.py — Autonomous parameter calibration for Expanded_Eurofer_CD.

Explores the parameter space using Latin Hypercube Sampling, runs bin-moment
cluster dynamics simulations via the C++ solver, and identifies parameter
combinations that match experimental TEM data at 10 dpa for EUROFER97:

    SIA loop density  ~ 10^22 m^-3       diameter ~ 10 nm
    Void density      ~ 3-5 x 10^21 m^-3 diameter ~ 2-4 nm

Usage
-----
    cd Expanded_Eurofer_CD
    python codes/calibrate.py              # default 20 + 10 runs
    python codes/calibrate.py --phase1 30  # more exploration
    python codes/calibrate.py --timeout 120  # 2 min per run
"""

import sys, os, csv, time, argparse
from pathlib import Path
import numpy as np

# ── Setup path ────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR   = SCRIPT_DIR.parent
sys.path.insert(0, str(BASE_DIR.parent))   # repo root

from Expanded_Eurofer_CD.py_utils.simulation import ExpandedEuroferCDSimulation

# ── Physical constants for size conversion ────────────────────────────────────
OMEGA = 1.18e-29     # atomic volume [m^3]
B_111 = 2.482e-10    # Burgers vector 1/2<111> [m]


# ══════════════════════════════════════════════════════════════════════════════
# Parameter space
# ══════════════════════════════════════════════════════════════════════════════

PARAM_SPEC = [
    # (name,         sheet,                 lo,    hi,   scale)
    # -- Production --
    ('eta',          'production_fission',  0.15,  0.35, 'linear'),
    ('f_cl_i',       'production_fission',  0.10,  0.55, 'linear'),
    ('f_cl_v',       'production_fission',  0.02,  0.15, 'linear'),
    # -- Diffusion: 1D glide & mobility cutoff --
    ('E_m_1D',       'diffusion',           0.01,  0.20, 'linear'),
    ('i_mobile',     'diffusion',           3,     100,  'int'),
    ('L_hat',        'diffusion',           20,    100,  'linear'),
    # -- Diffusion: solute trapping (dominant effect on SIA mobility) --
    ('c_C',          'diffusion',           1e-5,  5e-4, 'log'),
    ('E_b_C_SIA',    'diffusion',           0.20,  0.50, 'linear'),
    # -- Reactions: sinks and coalescence bias --
    ('rho_d',        'reactions',           1e12,  1e14, 'log'),
    ('Z_i',          'reactions',           1.02,  1.20, 'linear'),
    ('Z_ii',         'reactions',           1.0,   10.0, 'log'),
]

PARAM_NAMES = [p[0] for p in PARAM_SPEC]


def latin_hypercube(n_samples, n_dims, rng):
    """Simple LHS: stratified random in [0,1]^d."""
    result = np.zeros((n_samples, n_dims))
    for j in range(n_dims):
        perm = rng.permutation(n_samples)
        for i in range(n_samples):
            result[perm[i], j] = (i + rng.random()) / n_samples
    return result


def unit_to_physical(u, spec):
    """Map u in [0,1] to physical value according to spec."""
    name, sheet, lo, hi, scale = spec
    if scale == 'log':
        return 10 ** (np.log10(lo) + u * (np.log10(hi) - np.log10(lo)))
    elif scale == 'int':
        return int(round(lo + u * (hi - lo)))
    else:
        return lo + u * (hi - lo)


def physical_to_unit(val, spec):
    """Map physical value back to [0,1]."""
    name, sheet, lo, hi, scale = spec
    if scale == 'log':
        return (np.log10(val) - np.log10(lo)) / (np.log10(hi) - np.log10(lo))
    elif scale == 'int':
        return (val - lo) / (hi - lo)
    else:
        return (val - lo) / (hi - lo)


def sample_to_params(unit_vec):
    """Convert a [0,1]^d vector to a dict of physical parameters."""
    return {spec[0]: unit_to_physical(u, spec)
            for u, spec in zip(unit_vec, PARAM_SPEC)}


# ══════════════════════════════════════════════════════════════════════════════
# Objective function
# ══════════════════════════════════════════════════════════════════════════════

# Experimental targets at 10 dpa
TARGET_DOSE      = 10.0      # dpa
TARGET_N_LOOPS   = 1.0e22    # m^-3  (TEM-visible SIA loop density)
TARGET_D_I_NM    = 10.0      # nm    (TEM-visible SIA loop diameter)
TARGET_N_VOIDS   = 4.0e21    # m^-3  (TEM-visible void density)
TARGET_D_V_NM    = 3.0       # nm    (TEM-visible void diameter)

# TEM visibility threshold: clusters below this size are invisible
N_MIN_TEM = 10   # atoms (~1 nm diameter)


def atoms_to_loop_diameter_nm(mean_n_i):
    """Convert mean SIA cluster size (atoms) to loop diameter (nm)."""
    if mean_n_i <= 0:
        return 0.0
    r = np.sqrt(mean_n_i * OMEGA / (np.pi * B_111))
    return 2 * r * 1e9   # m -> nm


def atoms_to_void_diameter_nm(mean_n_v):
    """Convert mean vacancy cluster size (atoms) to void diameter (nm)."""
    if mean_n_v <= 0:
        return 0.0
    r = (3 * mean_n_v * OMEGA / (4 * np.pi)) ** (1.0 / 3.0)
    return 2 * r * 1e9   # m -> nm


def evaluate_cost(results, sim):
    """
    Compute scalar cost measuring deviation from experimental targets at
    TARGET_DOSE dpa.  Uses TEM-visible clusters only (n,m >= N_MIN_TEM).
    Returns (cost, metrics_dict) or (inf, None) on failure.
    """
    if results is None:
        return float('inf'), None

    dose = results['dose']
    if dose[-1] < TARGET_DOSE * 0.9:
        return float('inf'), None

    idx = np.argmin(np.abs(dose - TARGET_DOSE))
    inv_Omega = 1.0 / OMEGA

    # Reconstruct per-size distributions at the target dose
    re_obj = sim.rate_equations
    N = sim.input_data.I
    M = sim.input_data.V
    yj = np.maximum(results['y'][:, idx], 0.0)

    is_bin = hasattr(re_obj, 'bins')
    if is_bin:
        from Expanded_Eurofer_CD.py_utils.bin_moment_rates import \
            distribution_from_moments_hat
        i_d = getattr(re_obj, 'i_discrete', 0)
        v_d = getattr(re_obj, 'v_discrete', 0)
        K_bins = re_obj.I_bin if hasattr(re_obj, 'I_bin') else len(re_obj.bins)
        c_n = np.zeros(N)
        c_n[:i_d] = yj[:i_d]
        if K_bins > 0:
            mom = yj[i_d:]
            mu0 = mom[0::2][:K_bins]
            mu1 = mom[1::2][:K_bins]
            c_binned = distribution_from_moments_hat(mu0, mu1, re_obj.bins, N)
            c_n[i_d:] = c_binned[i_d:]
        i_VAC = getattr(re_obj, 'i_VAC', i_d + 2 * K_bins)
        K_v = getattr(re_obj, 'V_bin', getattr(re_obj, 'K_v', 0))
        c_v = np.zeros(M)
        c_v[:v_d] = yj[i_VAC:i_VAC + v_d]
        if K_v > 0:
            vac_start = i_VAC + v_d
            vmu0 = yj[vac_start::2][:K_v]
            vmu1 = yj[vac_start + 1::2][:K_v]
            c_v_binned = distribution_from_moments_hat(
                vmu0, vmu1, re_obj.vac_bins, M)
            c_v[v_d:] = c_v_binned[v_d:]
        elif v_d == 0:
            c_v = yj[i_VAC:i_VAC + M]
    else:
        i_SIA = getattr(re_obj, 'i_SIA', 0)
        i_VAC = getattr(re_obj, 'i_VAC', N)
        c_n = yj[i_SIA:i_VAC]
        c_v = yj[i_VAC:i_VAC + M]

    # TEM-visible quantities: only clusters with n,m >= N_MIN_TEM
    ns = np.arange(1, N + 1, dtype=float)
    ms = np.arange(1, M + 1, dtype=float)
    tem_i = ns >= N_MIN_TEM
    tem_v = ms >= N_MIN_TEM

    # SIA loops (TEM-visible)
    N_loops_vis = np.sum(c_n[tem_i]) * inv_Omega          # [m^-3]
    content_vis = np.dot(ns[tem_i], c_n[tem_i])
    count_vis   = np.sum(c_n[tem_i])
    mean_n_i_vis = content_vis / max(count_vis, 1e-300)
    d_i_nm = atoms_to_loop_diameter_nm(mean_n_i_vis)

    # Voids (TEM-visible)
    N_voids_vis  = np.sum(c_v[tem_v]) * inv_Omega         # [m^-3]
    vcontent_vis = np.dot(ms[tem_v], c_v[tem_v])
    vcount_vis   = np.sum(c_v[tem_v])
    mean_n_v_vis = vcontent_vis / max(vcount_vis, 1e-300)
    d_v_nm = atoms_to_void_diameter_nm(mean_n_v_vis)

    # Guard against zeros
    if N_loops_vis < 1e10 or N_voids_vis < 1e10:
        return float('inf'), None
    if d_i_nm < 0.01 or d_v_nm < 0.01:
        return float('inf'), None

    # Log-ratio penalty for densities, relative penalty for sizes
    cost = (
        (np.log10(N_loops_vis) - np.log10(TARGET_N_LOOPS)) ** 2
      + ((d_i_nm - TARGET_D_I_NM) / TARGET_D_I_NM) ** 2
      + (np.log10(N_voids_vis) - np.log10(TARGET_N_VOIDS)) ** 2
      + ((d_v_nm - TARGET_D_V_NM) / TARGET_D_V_NM) ** 2
    )

    metrics = {
        'dose_dpa':   dose[idx],
        'N_loops':    N_loops_vis,
        'N_voids':    N_voids_vis,
        'mean_n_i':   mean_n_i_vis,
        'mean_n_v':   mean_n_v_vis,
        'd_i_nm':     d_i_nm,
        'd_v_nm':     d_v_nm,
        'delta_FP':   results['delta_FP'][idx],
        'delta_He':   results['delta_He'][idx],
    }
    return cost, metrics


# ══════════════════════════════════════════════════════════════════════════════
# Simulation runner
# ══════════════════════════════════════════════════════════════════════════════

def run_single(params, timeout_s=180, run_id=0):
    """
    Run one simulation with the given parameter dict.
    Returns (cost, metrics, wall_time) or (inf, None, wall_time).
    """
    t0 = time.time()

    try:
        sim = ExpandedEuroferCDSimulation(
            I=5000, V=5000,
            solver_mode='cpp_full',
            physics_option='bin_moment_CD_fission',
            he_options='quasi_steady_state',
        )

        # ── Apply parameter overrides ────────────────────────────────────────
        inp = sim.input_data
        for pname, sheet, *_ in PARAM_SPEC:
            val = params[pname]
            if sheet == 'production_fission':
                inp.production_fission[pname] = val
            elif sheet == 'diffusion':
                inp.diffusion[pname] = val
            elif sheet == 'reactions':
                inp.reactions[pname] = val

        # Propagate i_mobile and v_mobile to reactions dict too (used by solver)
        if 'i_mobile' in params:
            inp.reactions['i_mobile'] = int(params['i_mobile'])
            inp.diffusion['i_mobile'] = int(params['i_mobile'])

        # Recompute derived quantities and rebuild rate constants
        inp._calculate_derived()
        sim.rebuild_rates()

        # ── Solver config (fast) ─────────────────────────────────────────────
        G = float(inp.reactions.get('G', 1e-6))
        t_end = TARGET_DOSE / G   # simulate exactly to target dose

        solver_config = {
            't_span': (1e-3, t_end),
            'n_points': 50,
            'log_time': True,
            'rtol': 1e-5,
            'atol': 0.0,
        }

        results = sim.run(
            solver_config=solver_config,
            save_output=False,
            timeout_s=timeout_s,
        )

        cost, metrics = evaluate_cost(results, sim)

    except Exception as exc:
        print(f"  [run {run_id}] EXCEPTION: {exc}")
        cost, metrics = float('inf'), None

    wall = time.time() - t0
    return cost, metrics, wall


# ══════════════════════════════════════════════════════════════════════════════
# Logging
# ══════════════════════════════════════════════════════════════════════════════

LOG_FIELDS = (
    ['run_id', 'phase', 'cost', 'wall_s']
    + PARAM_NAMES
    + ['N_loops', 'N_voids', 'mean_n_i', 'mean_n_v',
       'd_i_nm', 'd_v_nm', 'delta_FP', 'delta_He']
)


def append_log(log_path, row_dict):
    """Append one row to the CSV log, creating headers if needed."""
    write_header = not log_path.exists()
    with open(log_path, 'a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=LOG_FIELDS, extrasaction='ignore')
        if write_header:
            w.writeheader()
        w.writerow(row_dict)


def build_row(run_id, phase, params, cost, metrics, wall):
    """Assemble a log row dict."""
    row = {'run_id': run_id, 'phase': phase, 'cost': f'{cost:.6f}',
           'wall_s': f'{wall:.1f}'}
    for p in PARAM_NAMES:
        row[p] = f'{params[p]:.6g}'
    if metrics:
        for k, v in metrics.items():
            row[k] = f'{v:.6g}'
    return row


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='Calibrate Expanded_Eurofer_CD')
    parser.add_argument('--phase1', type=int, default=20,
                        help='Number of Phase 1 (exploration) samples')
    parser.add_argument('--phase2', type=int, default=10,
                        help='Number of Phase 2 (refinement) samples')
    parser.add_argument('--timeout', type=int, default=180,
                        help='Max seconds per run')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    log_path = BASE_DIR / 'output' / 'calibration_log.csv'
    os.makedirs(log_path.parent, exist_ok=True)

    # Remove old log
    if log_path.exists():
        log_path.unlink()

    all_results = []   # (cost, params, metrics)
    run_counter = 0

    # ── Phase 1: Latin Hypercube exploration ──────────────────────────────────
    n1 = args.phase1
    n_dims = len(PARAM_SPEC)
    lhs = latin_hypercube(n1, n_dims, rng)

    print(f"\n{'='*70}")
    print(f"  Phase 1: Latin Hypercube Sampling — {n1} runs")
    print(f"  Timeout: {args.timeout}s per run")
    print(f"  Log: {log_path}")
    print(f"{'='*70}\n")

    for i in range(n1):
        params = sample_to_params(lhs[i])
        param_str = '  '.join(f'{k}={params[k]:.4g}' for k in PARAM_NAMES)
        print(f"[{i+1}/{n1}] {param_str}")

        cost, metrics, wall = run_single(params, timeout_s=args.timeout,
                                         run_id=run_counter)
        all_results.append((cost, params, metrics))
        row = build_row(run_counter, 'phase1', params, cost, metrics, wall)
        append_log(log_path, row)

        status = f"cost={cost:.4f}" if cost < 1e10 else "FAILED"
        summary = ""
        if metrics:
            summary = (f"  N_loops={metrics['N_loops']:.2e}"
                       f"  d_i={metrics['d_i_nm']:.1f}nm"
                       f"  N_voids={metrics['N_voids']:.2e}"
                       f"  d_v={metrics['d_v_nm']:.1f}nm")
        print(f"  -> {status}  ({wall:.0f}s){summary}\n")
        run_counter += 1

    # ── Phase 2: Refinement around best candidates ────────────────────────────
    n2 = args.phase2
    if n2 > 0:
        # Sort by cost, take top 3
        valid = [(c, p, m) for c, p, m in all_results if c < 1e10]
        valid.sort(key=lambda x: x[0])

        if not valid:
            print("\nNo valid Phase 1 runs — skipping Phase 2.")
        else:
            n_best = min(3, len(valid))
            print(f"\n{'='*70}")
            print(f"  Phase 2: Refinement — {n2} runs around top {n_best}")
            print(f"{'='*70}\n")

            # Generate refinement samples: perturb top candidates by +-20%
            refine_samples = []
            for rank in range(n_best):
                _, best_p, _ = valid[rank]
                n_per = max(1, n2 // n_best + (1 if rank < n2 % n_best else 0))
                for _ in range(n_per):
                    perturbed = {}
                    for spec in PARAM_SPEC:
                        pname = spec[0]
                        base_val = best_p[pname]
                        u_base = physical_to_unit(base_val, spec)
                        # Perturb in unit space by +-20%
                        delta = 0.20 * (rng.random() * 2 - 1)
                        u_new = np.clip(u_base + delta, 0.0, 1.0)
                        perturbed[pname] = unit_to_physical(u_new, spec)
                    refine_samples.append(perturbed)

            for i, params in enumerate(refine_samples):
                param_str = '  '.join(f'{k}={params[k]:.4g}'
                                      for k in PARAM_NAMES)
                print(f"[R{i+1}/{len(refine_samples)}] {param_str}")

                cost, metrics, wall = run_single(
                    params, timeout_s=args.timeout, run_id=run_counter)
                all_results.append((cost, params, metrics))
                row = build_row(run_counter, 'phase2', params, cost,
                                metrics, wall)
                append_log(log_path, row)

                status = f"cost={cost:.4f}" if cost < 1e10 else "FAILED"
                summary = ""
                if metrics:
                    summary = (f"  N_loops={metrics['N_loops']:.2e}"
                               f"  d_i={metrics['d_i_nm']:.1f}nm"
                               f"  N_voids={metrics['N_voids']:.2e}"
                               f"  d_v={metrics['d_v_nm']:.1f}nm")
                print(f"  -> {status}  ({wall:.0f}s){summary}\n")
                run_counter += 1

    # ── Final report ──────────────────────────────────────────────────────────
    valid = [(c, p, m) for c, p, m in all_results if c < 1e10]
    valid.sort(key=lambda x: x[0])

    print(f"\n{'='*70}")
    print(f"  CALIBRATION COMPLETE — {run_counter} runs, {len(valid)} valid")
    print(f"{'='*70}")

    if not valid:
        print("\n  No valid runs found. Try relaxing bounds or increasing timeout.")
        return

    print(f"\n  Experimental targets at {TARGET_DOSE} dpa:")
    print(f"    SIA loop density:  {TARGET_N_LOOPS:.1e} m^-3")
    print(f"    SIA loop diameter: {TARGET_D_I_NM} nm")
    print(f"    Void density:      {TARGET_N_VOIDS:.1e} m^-3")
    print(f"    Void diameter:     {TARGET_D_V_NM} nm")

    print(f"\n  Top {min(5, len(valid))} results:\n")
    print(f"  {'Rank':>4}  {'Cost':>8}  {'N_loops':>10}  {'d_i(nm)':>8}"
          f"  {'N_voids':>10}  {'d_v(nm)':>8}")
    print(f"  {'-'*60}")

    for rank, (cost, params, metrics) in enumerate(valid[:5]):
        print(f"  {rank+1:>4}  {cost:8.4f}"
              f"  {metrics['N_loops']:10.2e}  {metrics['d_i_nm']:8.1f}"
              f"  {metrics['N_voids']:10.2e}  {metrics['d_v_nm']:8.1f}")

    # Print best parameter set
    best_cost, best_params, best_metrics = valid[0]
    print(f"\n  Best parameter set (cost={best_cost:.4f}):\n")
    print(f"  {'Parameter':>12}  {'Value':>12}  {'Sheet'}")
    print(f"  {'-'*45}")
    for spec in PARAM_SPEC:
        pname, sheet = spec[0], spec[1]
        val = best_params[pname]
        if isinstance(val, int):
            print(f"  {pname:>12}  {val:>12}  {sheet}")
        else:
            print(f"  {pname:>12}  {val:>12.4g}  {sheet}")

    print(f"\n  Full log: {log_path}\n")


if __name__ == '__main__':
    main()
