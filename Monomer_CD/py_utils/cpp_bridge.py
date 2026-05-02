"""
cpp_bridge.py – Python bridge for the ClusterDynamics C++ SUNDIALS solver.

Responsibilities
----------------
1. Collect all pre-computed rate constants and physics parameters from
   InputData / ReactionRates and write them to a temporary parameter file.
2. Invoke solver.exe with --param_file=<path> (avoids Windows CreateProcess
   command-line length limits that arise when Nv/Ni are large).
3. Parse the solver stdout (t + N concentrations per row) and reconstruct
   the standard results dict via post_process.calculate_derived_quantities.

Parameter file format
---------------------
One "key=value" entry per line (blank lines and lines starting with '#'
are ignored).  Arrays are written as individual indexed entries:
  KCV_0=<val>
  KCV_1=<val>
  ...
This mirrors the key names that parameters.h / build_parameters() expect.

The C++ solver output format is one row per time point:
  t  Cv1 Cv2 ... Cv_Nv  Ci1 Ci2 ... Ci_Ni
  (1 + Nv + Ni space-separated values, scientific notation)
"""

import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)


# ── Parameter file writer ─────────────────────────────────────────────────────

def write_param_file(sim, solver_config, path):
    """
    Write all solver parameters to a text file (key=value, one per line).

    Parameters
    ----------
    sim           : ClusterDynamicsSimulation  – fully initialised
    solver_config : dict  – t_span, rtol, atol, log_time, n_points, backend, ...
    path          : str or Path  – file to write
    """
    inp = sim.input_data
    rr  = sim.reaction_rates
    re  = sim.rate_equations

    lines = []

    method_opts = solver_config.get('solver_method', {})

    # ── Cluster size limits ───────────────────────────────────────────────────
    Nv = int(inp.model_params['Nv'])
    Ni = int(inp.model_params['Ni'])
    lines.append(f"Nv={Nv}")
    lines.append(f"Ni={Ni}")

    # ── Rate constant arrays ──────────────────────────────────────────────────
    for k, v in enumerate(rr.KCV):
        lines.append(f"KCV_{k}={v:.17e}")
    for k, v in enumerate(rr.KCI):
        lines.append(f"KCI_{k}={v:.17e}")
    for k, v in enumerate(rr.KLV):
        lines.append(f"KLV_{k}={v:.17e}")
    for k, v in enumerate(rr.KLI):
        lines.append(f"KLI_{k}={v:.17e}")
    for k, v in enumerate(rr.GCV):
        lines.append(f"GCV_{k}={v:.17e}")
    for k, v in enumerate(rr.GLV):
        lines.append(f"GLV_{k}={v:.17e}")

    # ── Extend KLV/KLI/GLV analytically to Ni_max ────────────────────────────
    # KLV[k] = KLV[0]*sqrt(k+1)  (Kl_v scales as x^0.5, exact)
    # KLI[k] = KLI[0]*sqrt(k+1)  (Kl_i scales as x^0.5, exact)
    # GLV[k] = 0 for k >= 1      (diinterstitial dissociation only at k=0)
    _Ni_max = int(method_opts.get('Ni_max', Ni))
    if _Ni_max > Ni:
        _klv0 = float(rr.KLV[0])
        _kli0 = float(rr.KLI[0])
        for k in range(Ni, _Ni_max):
            _x = k + 1
            lines.append(f"KLV_{k}={_klv0 * _x**0.5:.17e}")
            lines.append(f"KLI_{k}={_kli0 * _x**0.5:.17e}")
            lines.append(f"GLV_{k}=0.00000000000000000e+00")
        _c_floor = float(inp.model_params.get('C_floor', 1e-100))
        for k in range(Nv + Ni, Nv + _Ni_max):
            lines.append(f"y0_{k}={_c_floor:.17e}")

    # ── Scalar physics ────────────────────────────────────────────────────────
    p = inp.material_params
    d = inp.derived
    for key, val in [
        ("P_prod",  p['P']),
        ("alpha",   d['alpha']),
        ("Cv_eq",   d['Cv_eq']),
        ("C2v_eq",  d['C2v_eq']),
        ("Z_v",     p['Z_v']),
        ("Z_i",     p['Z_i']),
        ("rho_d",   p['rho_d']),
        ("Dv",      d['Dv']),
        ("D2v",     d['D2v']),
        ("Di",      d['Di']),
        ("K_nuc_i", d['K_nuc_i']),
        ("C_floor", inp.model_params.get('C_floor', 1e-100)),
    ]:
        lines.append(f"{key}={val:.17e}")

    # ── Initial conditions ────────────────────────────────────────────────────
    y0 = re.get_initial_conditions()
    for k, v in enumerate(y0):
        lines.append(f"y0_{k}={v:.17e}")

    # ── Solver settings ────────────────────────────────────────────────────────
    t_span  = solver_config.get('t_span', (1e-6, 1e5))
    lines.append(f"t_begin={t_span[0]:.17e}")
    lines.append(f"t_end={t_span[1]:.17e}")
    lines.append(f"n_points={int(solver_config.get('n_points', 200))}")
    lines.append(f"log_time={1.0 if solver_config.get('log_time', True) else 0.0}")
    lines.append(f"rtol={solver_config.get('rtol', 1e-8):.17e}")
    lines.append(f"atol={solver_config.get('atol', 1e-50):.17e}")

    # ── Integration method ────────────────────────────────────────────────────
    _ark_table_map = {
        'SDIRK_2_1_2':            100, 'SDIRK_5_3_4':           107,
        'KVAERNO_7_4_5':          110, 'ARK548L2SA_DIRK_8_4_5': 111,
        'ESDIRK547L2SA_7_4_5':    121,
    }
    _backend_map = {'cvode': 0, 'arkode': 1}
    _lmm_map     = {'bdf': 2, 'adams': 1}
    _linsol_map  = {'dense': 0, 'band': 1, 'banded': 1, 'gmres': 2}

    lines.append(f"backend={_backend_map.get(str(method_opts.get('backend', 'cvode')).lower(), 0)}")
    lines.append(f"lmm={_lmm_map.get(str(method_opts.get('lmm', 'bdf')).lower(), 2)}")
    lines.append(f"linsol={_linsol_map.get(str(method_opts.get('linsol', 'dense')).lower(), 0)}")
    lines.append(f"mu={int(method_opts.get('mu', re.N - 1))}")
    lines.append(f"ml={int(method_opts.get('ml', re.N - 1))}")
    lines.append(f"max_order={int(method_opts.get('max_order', 0))}")
    _ark_raw = method_opts.get('ark_table', 'ARK548L2SA_DIRK_8_4_5')
    ark_val  = _ark_raw if isinstance(_ark_raw, int) else _ark_table_map.get(str(_ark_raw).upper(), 111)
    lines.append(f"ark_table={ark_val}")

    # ── Dynamic window solver (Phase I + Phase II) ────────────────────────────
    lines.append(f"window_mode={int(method_opts.get('window_mode', 0))}")
    lines.append(f"window_w0_v={int(method_opts.get('window_w0_v', Nv))}")
    lines.append(f"window_w0_i={int(method_opts.get('window_w0_i', Ni))}")
    lines.append(f"window_C_expand={float(method_opts.get('window_C_expand', 1e-18)):.17e}")
    lines.append(f"window_expand_pad={int(method_opts.get('window_expand_pad', 10))}")
    lines.append(f"window_expand_factor={float(method_opts.get('window_expand_factor', 0.0)):.17e}")
    lines.append(f"window_check_every={int(method_opts.get('window_check_every', 1))}")
    # Phase II only
    lines.append(f"window_C_contract={float(method_opts.get('window_C_contract', 0.0)):.17e}")
    lines.append(f"window_min_active_i={int(method_opts.get('window_min_active_i', 5))}")
    lines.append(f"window_prec={int(method_opts.get('window_prec', 0))}")
    lines.append(f"window_nuc_guard={float(method_opts.get('window_nuc_guard', 0.0)):.17e}")
    # Phase III
    lines.append(f"window_width={int(method_opts.get('window_width', 500))}")
    lines.append(f"window_t_start={float(method_opts.get('window_t_start', 10.0)):.17e}")
    lines.append(f"window_N_thresh={int(method_opts.get('window_N_thresh', 1000))}")
    # Phase IV: Multithread-OpenMP
    lines.append(f"window_omp_threads={int(method_opts.get('window_omp_threads', 0))}")
    # Dynamic Ni extension
    lines.append(f"Ni_max={int(method_opts.get('Ni_max', Ni))}")
    lines.append(f"Ni_extend_tol={float(method_opts.get('Ni_extend_tol', 0.0)):.17e}")
    lines.append(f"Ni_extend_margin={int(method_opts.get('Ni_extend_margin', 0))}")

    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


# ── Output parsing ─────────────────────────────────────────────────────────────

def _parse_stdout(text, N):
    """
    Convert solver stdout (space-separated, n_points × (1+N)) to a numpy array.

    Returns
    -------
    numpy.ndarray  shape (n_points, 1+N)
    """
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
    Run the compiled ClusterDynamics C++ solver and return the standard results dict.

    Parameters
    ----------
    sim           : ClusterDynamicsSimulation – fully initialised
    solver_config : dict  – same keys as SOLVER_CONFIG in the notebook
    base_dir      : Path or None  – ClusterDynamics/ root; auto-detected if None

    Returns
    -------
    dict (same format as post_process.calculate_derived_quantities) or None on failure
    """
    from py_utils.post_process import calculate_derived_quantities

    if base_dir is None:
        base_dir = Path(__file__).resolve().parent.parent

    exe_name  = 'solver.exe' if sys.platform == 'win32' else 'solver'
    build_dir = Path(base_dir) / 'build'
    exe_path  = build_dir / 'Debug'   / exe_name
    if not exe_path.exists():
        exe_path = build_dir / 'Release' / exe_name
    if not exe_path.exists():
        exe_path = build_dir / exe_name

    if not exe_path.exists():
        print(f"❌ solver executable not found at {build_dir / exe_name}")
        print("   Build it with:")
        print(f"     cd {Path(base_dir) / 'cpp_utils'}")
        print( "     cmake -S . -B ../build -DCMAKE_BUILD_TYPE=Release")
        print( "     cmake --build ../build --config Release")
        return None

    # Write parameters to a temp file so we don't hit the Windows
    # CreateProcess command-line length limit (WinError 206) with large Nv/Ni.
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt',
                                     delete=False, prefix='cd_params_') as tf:
        param_path = tf.name

    # Binary output file: same name as the param file but with .bin extension.
    # solver.exe derives this path automatically and writes raw float64 rows.
    bin_path = param_path[:-4] + '.bin' if param_path.endswith('.txt') else param_path + '.bin'

    method_opts = solver_config.get('solver_method', {})

    try:
        write_param_file(sim, solver_config, param_path)

        N   = sim.rate_equations.N
        Nv  = int(sim.input_data.model_params['Nv'])
        Ni  = int(sim.input_data.model_params['Ni'])
        _Ni_max = int(method_opts.get('Ni_max', Ni))
        if _Ni_max > Ni:
            N = Nv + _Ni_max   # output has Ni_max columns
        print(f"Running C++ solver ({exe_path.name})  …  "
              f"(Nv={Nv}, Ni={Ni}, N_EQ={N}, param_file)")
        proc = subprocess.run(
            [str(exe_path), f'--param_file={param_path}'],
            capture_output=True,
        )
    finally:
        try:
            os.unlink(param_path)
        except OSError:
            pass

    if proc.returncode != 0:
        print(f"❌ C++ solver failed (exit code {proc.returncode}):")
        stderr_text = proc.stderr.decode('utf-8', errors='replace') if isinstance(proc.stderr, bytes) else proc.stderr
        print(stderr_text[:5000])
        try:
            os.unlink(bin_path)
        except OSError:
            pass
        return None

    stderr_text = proc.stderr.decode('utf-8', errors='replace') if isinstance(proc.stderr, bytes) else proc.stderr
    if stderr_text.strip():
        print("C++ solver info:\n" + stderr_text)

    # Read output: prefer the binary file written by the solver (fast bulk read);
    # fall back to text stdout parsing if the binary file is absent or malformed.
    sol_arr = None
    try:
        raw = np.fromfile(bin_path, dtype=np.float64)
        n_cols = 1 + N
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
        # Binary read failed — fall back to text stdout
        stdout_text = proc.stdout.decode('utf-8', errors='replace') if isinstance(proc.stdout, bytes) else proc.stdout
        sol_arr = _parse_stdout(stdout_text, N)

    if sol_arr.shape[0] == 0:
        print("❌ C++ solver produced no parseable output")
        return None

    n_pts = sol_arr.shape[0]
    print(f"✓ C++ solver completed — {n_pts} time points")

    t  = sol_arr[:, 0]           # (n_pts,)
    y  = sol_arr[:, 1:].T        # (N, n_pts)

    # If Ni_max > Ni, trim the output back to Nv+Ni so post_process works correctly
    # with rate_equations (which has Ni, not Ni_max).  The extended columns are all
    # C_floor values (1e-100) for cluster sizes beyond the initial Ni.
    _N_trim = sim.rate_equations.N   # Nv + Ni (original)
    if y.shape[0] > _N_trim:
        y = y[:_N_trim, :]

    results = calculate_derived_quantities(t, y, sim.input_data, sim.rate_equations)

    backend_name = str(method_opts.get('backend', 'cvode')).lower()
    linsol_label = str(method_opts.get('linsol', 'dense')).upper()
    if backend_name == 'arkode':
        ark_tbl     = str(method_opts.get('ark_table', 'ARK548L2SA_DIRK_8_4_5'))
        backend_str = f'C++ SUNDIALS ARKODE ARKStep DIRK ({ark_tbl}) / {linsol_label}'
    else:
        lmm_label   = str(method_opts.get('lmm', 'bdf')).upper()
        win_mode    = int(method_opts.get('window_mode', 0))
        if win_mode == 4:
            n_thr = int(method_opts.get('window_omp_threads', 0))
            thr_label = f'{n_thr} threads' if n_thr > 0 else 'OMP_NUM_THREADS'
            backend_str = (f'C++ SUNDIALS CVODE {lmm_label} WINDOW-IV (OpenMP, {thr_label}) / GMRES')
        elif win_mode:
            w0_v = int(method_opts.get('window_w0_v', Nv))
            w0_i = int(method_opts.get('window_w0_i', Ni))
            backend_str = (f'C++ SUNDIALS CVODE {lmm_label} WINDOW / GMRES'
                           f' (w0_v={w0_v}, w0_i={w0_i})')
        else:
            backend_str = f'C++ SUNDIALS CVODE {lmm_label} / {linsol_label}'

    results['metadata'] = {
        'solver_stats': {
            'success':       True,
            'message':       backend_str,
            'n_time_points': n_pts,
        },
    }
    print("✓ Results processing complete!")
    return results
