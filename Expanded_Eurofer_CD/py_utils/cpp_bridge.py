"""
cpp_bridge.py — Python → C++ solver bridge for Expanded_Eurofer_CD.

Responsibilities
----------------
1. Collect pre-computed rate constants and physics parameters from
   InputData / ReactionRates and write them to a temporary parameter file.
2. Invoke solver.exe with --param_file=<path>.
3. Parse the binary output and reconstruct the standard results dict.

Parameter file format
---------------------
One "key=value" entry per line.  Arrays use indexed entries:
  KVV_0=<val>, KVV_1=<val>, ...

New fields vs. Eurofer_CD
--------------------------
- solver_mode:     integer (0=cpp_full, 3=cpp_sliding_win, 4=sliding_OpenMP)
- physics_option:  integer (0=full_CD_fission, 1=full_CD_fusion,
                             2=bin_moment_fission, 3=bin_moment_fusion)
- A_sph, A_loop, A_1D, B_rot: geometric prefactors
- trap_SIA, trap_VAC, trap_loop: solute trapping sums
- K_1D_pref_k: 1D glide prefactor for SIA cluster k
- Bin-moment parameters: K_bins, r_ratio, n_moments, n1_bin
- n_max_i, m_max_v: mobility cutoffs
"""

import logging
import os
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

_SOLVER_MODE_MAP = {
    'cpp_full':        0,
    'cpp_sliding_win': 3,
    'sliding_OpenMP':  4,
}
_PHYSICS_OPTION_MAP = {
    'full_CD_fission':       0,
    'full_CD_fusion':        1,
    'bin_moment_CD_fission': 2,
    'bin_moment_CD_fusion':  3,
}


def write_param_file(sim, solver_config, path):
    """
    Write all solver parameters to a text file (key=value, one per line).

    Parameters
    ----------
    sim           : ExpandedEuroferCDSimulation — fully initialised
    solver_config : dict  — t_span, rtol, atol, solver_method, etc.
    path          : str or Path
    """
    inp  = sim.input_data
    rr   = sim.reaction_rates
    re_obj = sim.rate_equations

    d      = inp.derived
    method = solver_config.get('solver_method', {})

    N  = inp.N
    M  = inp.M

    lines = []

    # ── Cluster size limits ───────────────────────────────────────────────────
    lines.append(f"N={N}")
    lines.append(f"M={M}")
    lines.append(f"Ni={N}")    # legacy key
    lines.append(f"Nv={M}")    # legacy key
    lines.append(f"Ni_max={N}")

    # ── Solver mode and physics option ────────────────────────────────────────
    sm_int = _SOLVER_MODE_MAP.get(inp.solver_mode, 0)
    po_int = _PHYSICS_OPTION_MAP.get(inp.physics_option, 0)
    lines.append(f"solver_mode_int={sm_int}")
    lines.append(f"physics_option_int={po_int}")
    lines.append(f"window_mode={sm_int}")   # legacy: window_mode mirrors solver_mode

    # ── Geometric rate constant prefactors (Eq. 128) ─────────────────────────
    lines.append(f"A_sph={d['A_sph']:.17e}")
    lines.append(f"A_loop={d['A_loop']:.17e}")
    lines.append(f"B_rot={d['B_rot']:.17e}")
    lines.append(f"L_hat={d['L_hat']:.17e}")

    # ── Solute trapping sums (Eq. 42, 48, 52) ────────────────────────────────
    lines.append(f"trap_SIA={d['trap_SIA']:.17e}")
    lines.append(f"trap_VAC={d['trap_VAC']:.17e}")
    lines.append(f"trap_loop={d['trap_loop']:.17e}")

    # ── Mobility cutoffs ──────────────────────────────────────────────────────
    lines.append(f"n_max_i={d['n_max_i']}")
    lines.append(f"m_max_v={d['m_max_v']}")

    # ── Vacancy cluster rate arrays (0-indexed, size M) ───────────────────────
    for k, v in enumerate(rr.K_VAC_grow):
        lines.append(f"KVV_{k}={v:.17e}")
    for k, v in enumerate(rr.K_VAC_shrink):
        lines.append(f"KVI_{k}={v:.17e}")
    for k, v in enumerate(rr.G_VAC):
        lines.append(f"GVV_{k}={v:.17e}")
    for k, v in enumerate(rr.K_HeV):
        lines.append(f"KHeV_{k}={v:.17e}")
    for k, v in enumerate(re_obj.Pr_VAC):
        lines.append(f"Pr_VAC_{k}={v:.17e}")
    # m^{1/3} factors
    for k in range(M):
        m = k + 1
        lines.append(f"m13_{k}={m**(1.0/3.0):.17e}")

    # ── SIA cluster rate arrays (0-indexed, size N) ───────────────────────────
    for k, v in enumerate(rr.K_SIA_grow):
        lines.append(f"KII_{k}={v:.17e}")
    for k, v in enumerate(rr.K_SIA_shrink):
        lines.append(f"KIV_{k}={v:.17e}")
    for k, v in enumerate(rr.G_SIA):
        lines.append(f"GII_{k}={v:.17e}")
    for k, v in enumerate(rr.k2_SIA):
        lines.append(f"k2_SIA_{k}={v:.17e}")
    for k, v in enumerate(re_obj.Pr_SIA):
        lines.append(f"Pr_SIA_{k}={v:.17e}")

    # ── 1D glide prefactors K_1D_pref[n-1] (Eq. 141) ─────────────────────────
    for k, v in enumerate(rr.K_1D_pref):
        lines.append(f"K_1D_pref_{k}={v:.17e}")

    # Legacy K_IclV separable coefficients (computed from K_1D_pref)
    Omega = float(d['Omega'])
    Di_eff = float(d['Di_eff'])
    r0     = float(d['r0'])
    C0     = 4.0 * np.pi * r0 * Di_eff / Omega
    lines.append(f"K_IclV_ns_0=0.00000000000000000e+00")
    lines.append(f"K_IclV_ni_0=0.00000000000000000e+00")
    for k in range(1, N):
        n = k + 1
        lines.append(f"K_IclV_ns_{k}={C0 * n**(-2.0/3.0):.17e}")
        lines.append(f"K_IclV_ni_{k}={C0 / n:.17e}")

    # ── Scalar physics ────────────────────────────────────────────────────────
    kBT = float(d['kBT'])
    nu_h = float(inp.energetics.get('nu_h', 3.0e12))
    E_m_h = float(inp.energetics.get('E_m_h', 0.06))
    E_b_hV1 = float(inp.energetics.get('E_b_hV_1', 2.30))
    beta_He = nu_h * np.exp(-(E_b_hV1 + E_m_h) / kBT)

    lines.extend([
        f"G_He={re_obj.G_He:.17e}",
        f"k2_disl_v={rr.k2_vac_scalar:.17e}",
        f"k2_disl_i={rr.k2_SIA_scalar:.17e}",
        f"k2_disl_He={rr.k2_He_scalar:.17e}",
        f"Cv_eq={d['Cv_eq']:.17e}",
        f"beta_He={beta_He:.17e}",
        f"delta_He=-0.80000000000000000e+00",    # He pressure coeff [eV]
        f"beta_He_exp=0.70000000000000000e+00",  # He pressure exponent
        f"kBT={kBT:.17e}",
        f"L_He_max=10.00000000000000000e+00",    # fall-back cap
        f"K_iv={rr.K_iv:.17e}",
    ])

    # ── Bin-moment parameters ─────────────────────────────────────────────────
    if hasattr(re_obj, 'K'):
        K_bins = re_obj.K
        r_ratio = re_obj.r
        lines.append(f"K_bins={K_bins}")
        lines.append(f"r_ratio={r_ratio:.17e}")
        lines.append(f"n1_bin={re_obj.n1}")
    else:
        lines.append(f"K_bins=0")
        lines.append(f"r_ratio=2.00000000000000000e+00")
        lines.append(f"n1_bin=1")

    # ── He mode ───────────────────────────────────────────────────────────────
    he_mode_int = 0 if getattr(re_obj, 'he_mode', 'case2') == 'case2' else 1
    lines.append(f"he_mode={he_mode_int}")

    # ── Initial conditions ────────────────────────────────────────────────────
    y0 = re_obj.get_initial_conditions()
    for k, v in enumerate(y0):
        lines.append(f"y0_{k}={v:.17e}")

    C_floor = float(inp.reactions.get('C_floor', 1e-15))
    lines.append(f"C_floor={C_floor:.17e}")

    # Free He mode: 'dynamic'=0 integrates Eq.157; 'quasi_steady_state'=1 uses QSS
    he_options_str = str(inp.reactions.get('he_options', 'dynamic')).lower()
    qss_He_int = 1 if he_options_str == 'quasi_steady_state' else 0
    lines.append(f"qss_He={qss_He_int}")

    # ── Solver settings ────────────────────────────────────────────────────────
    t_span = solver_config.get('t_span', (1e-8, 1e7))
    lines.append(f"t_begin={t_span[0]:.17e}")
    lines.append(f"t_end={t_span[1]:.17e}")
    lines.append(f"n_points={int(solver_config.get('n_points', 200))}")
    lines.append(f"log_time={1.0 if solver_config.get('log_time', True) else 0.0}")
    lines.append(f"rtol={solver_config.get('rtol', 1e-8):.17e}")
    lines.append(f"atol={solver_config.get('atol', 1e-50):.17e}")

    # ── Integration method ────────────────────────────────────────────────────
    _backend_map = {'cvode': 0, 'arkode': 1}
    _lmm_map     = {'bdf': 2, 'adams': 1}
    _linsol_map  = {'dense': 0, 'band': 1, 'banded': 1, 'gmres': 2}

    N_tot = N + M + 1
    lines.append(f"backend={_backend_map.get(str(method.get('backend','cvode')).lower(), 0)}")
    lines.append(f"lmm={_lmm_map.get(str(method.get('lmm','bdf')).lower(), 2)}")
    lines.append(f"linsol={_linsol_map.get(str(method.get('linsol','dense')).lower(), 0)}")
    lines.append(f"mu={int(method.get('mu', N_tot - 1))}")
    lines.append(f"ml={int(method.get('ml', N_tot - 1))}")
    lines.append(f"max_order={int(method.get('max_order', 0))}")
    lines.append(f"ark_table={int(method.get('ark_table', 111))}")

    # ── Dynamic window parameters ──────────────────────────────────────────────
    lines.append(f"window_w0_v={M}")
    lines.append(f"window_w0_i={int(method.get('window_w0_i', N))}")
    lines.append(f"window_C_expand={float(method.get('window_C_expand', 1e-18)):.17e}")
    lines.append(f"window_expand_pad={int(method.get('window_expand_pad', 10))}")
    lines.append(f"window_expand_factor={float(method.get('window_expand_factor', 0.0)):.17e}")
    lines.append(f"window_check_every={int(method.get('window_check_every', 1))}")
    lines.append(f"window_C_contract={float(method.get('window_C_contract', 0.0)):.17e}")
    lines.append(f"window_min_active_i={int(method.get('window_min_active_i', 5))}")
    lines.append(f"window_prec={int(method.get('window_prec', 1))}")
    lines.append(f"window_nuc_guard={float(method.get('window_nuc_guard', 0.0)):.17e}")
    lines.append(f"window_width={int(method.get('window_width', 500))}")
    lines.append(f"window_t_start={float(method.get('window_t_start', 10.0)):.17e}")
    lines.append(f"window_N_thresh={int(method.get('window_N_thresh', 1000))}")
    lines.append(f"window_omp_threads={int(method.get('window_omp_threads', 0))}")
    lines.append(f"window_gmres_maxl={int(method.get('window_gmres_maxl', 20))}")
    lines.append(f"Ni_extend_tol=0.00000000000000000e+00")
    lines.append(f"Ni_extend_margin=0")

    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


# ── Output parsing ─────────────────────────────────────────────────────────────

def _parse_stdout(text, N):
    rows = []
    for line in text.strip().splitlines():
        parts = line.split()
        if len(parts) == 1 + N:
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                pass
    return np.array(rows) if rows else np.empty((0, 1 + N))


# ── Main entry point ──────────────────────────────────────────────────────────

def run_cpp_solver(sim, solver_config, base_dir=None):
    """
    Run the Expanded_Eurofer_CD C++ solver and return the standard results dict.

    Parameters
    ----------
    sim           : ExpandedEuroferCDSimulation
    solver_config : dict
    base_dir      : Path or None

    Returns
    -------
    dict or None
    """
    from .post_process import calculate_derived_quantities

    if base_dir is None:
        base_dir = Path(__file__).resolve().parent.parent

    exe_name  = 'solver.exe' if sys.platform == 'win32' else 'solver'
    build_dir = Path(base_dir) / 'build'
    exe_path  = build_dir / 'Release' / exe_name
    if not exe_path.exists():
        exe_path = build_dir / 'Debug' / exe_name
    if not exe_path.exists():
        exe_path = build_dir / exe_name

    if not exe_path.exists():
        print(f"C++ solver not found at {build_dir}")
        print("  Build with:")
        print(f"    cd {Path(base_dir) / 'cpp_utils'}")
        print( "    cmake -S . -B ../build -DCMAKE_BUILD_TYPE=Release")
        print( "    cmake --build ../build --config Release")
        return None

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt',
                                     delete=False, prefix='expanded_cd_') as tf:
        param_path = tf.name

    bin_path = param_path[:-4] + '.bin'
    re_obj   = sim.rate_equations
    N_tot    = re_obj.N_eq

    try:
        write_param_file(sim, solver_config, param_path)
        print(f"C++ solver: {exe_path.name}  N_eq={N_tot}"
              f"  solver_mode='{sim.input_data.solver_mode}'"
              f"  physics='{sim.input_data.physics_option}'")

        proc = subprocess.Popen(
            [str(exe_path), f'--param_file={param_path}'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        def _fwd_stderr():
            for raw in proc.stderr:
                sys.stderr.write(raw.decode('utf-8', errors='replace'))
                sys.stderr.flush()

        t_fwd = threading.Thread(target=_fwd_stderr, daemon=True)
        t_fwd.start()
        stdout_data = proc.stdout.read()
        proc.wait()
        t_fwd.join()
    finally:
        try:
            os.unlink(param_path)
        except OSError:
            pass

    if proc.returncode != 0:
        print(f"C++ solver failed (exit code {proc.returncode})")
        try:
            os.unlink(bin_path)
        except OSError:
            pass
        return None

    # Parse binary output
    sol_arr = None
    try:
        raw    = np.fromfile(bin_path, dtype=np.float64)
        n_cols = 1 + N_tot
        n_rows = raw.size // n_cols
        if n_rows > 0:
            sol_arr = raw[:n_rows * n_cols].reshape(n_rows, n_cols)
    except Exception:
        pass
    finally:
        try:
            os.unlink(bin_path)
        except OSError:
            pass

    if sol_arr is None or sol_arr.shape[0] == 0:
        text    = stdout_data.decode('utf-8', errors='replace')
        sol_arr = _parse_stdout(text, N_tot)

    if sol_arr.shape[0] == 0:
        print("C++ solver produced no parseable output")
        return None

    n_pts = sol_arr.shape[0]
    print(f"C++ solver completed — {n_pts} time points")

    t = sol_arr[:, 0]
    y = sol_arr[:, 1:].T   # (N_tot, n_pts)

    results = calculate_derived_quantities(t, y, sim.input_data, re_obj)

    sm     = sim.input_data.solver_mode
    po     = sim.input_data.physics_option
    linsol = str(solver_config.get('solver_method', {}).get('linsol', 'dense')).upper()
    results['metadata'] = {
        'solver_stats': {
            'success':       True,
            'message':       f'C++ CVODE BDF {sm}/{po} / {linsol}',
            'n_time_points': n_pts,
        },
    }
    print("Results processing complete.")
    return results
