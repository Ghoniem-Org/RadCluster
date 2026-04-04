"""
cpp_bridge.py – Python bridge for the Eurofer_CD C++ SUNDIALS solver.

Responsibilities
----------------
1. Collect all pre-computed rate constants and physics parameters from
   InputData / ReactionRates and write them to a temporary parameter file.
2. Invoke solver.exe with --param_file=<path>.
3. Parse the solver binary output and reconstruct the standard results dict
   via post_process.calculate_derived_quantities.

Parameter file format
---------------------
One "key=value" entry per line.  Arrays are written as individual indexed
entries:  KVV_0=<val>, KVV_1=<val>, ...

This mirrors the key names that parameters.h / build_parameters() expect.

Output state vector order (matches C++ layout):
  Ci1..Ci_Ni  Cv1..Cv_Nv  C_He
  (Ni + Nv + 1 columns after the time column)
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


# ── Parameter file writer ─────────────────────────────────────────────────────

def write_param_file(sim, solver_config, path):
    """
    Write all solver parameters to a text file (key=value, one per line).

    Parameters
    ----------
    sim           : EuroferCDSimulation  – fully initialised
    solver_config : dict  – t_span, rtol, atol, log_time, n_points,
                            solver_method (nested dict), he_mode
    path          : str or Path  – file to write
    """
    inp = sim.input_data
    rr  = sim.reaction_rates
    re  = sim.rate_equations

    lines = []
    method_opts = solver_config.get('solver_method', {})

    Ni = int(inp.model_params['Ni'])
    Nv = int(inp.model_params['Nv'])
    _Ni_max = int(method_opts.get('Ni_max', Ni))

    # ── Cluster size limits ───────────────────────────────────────────────────
    lines.append(f"Ni={Ni}")
    lines.append(f"Nv={Nv}")
    lines.append(f"Ni_max={_Ni_max}")

    # ── Vacancy cluster arrays ────────────────────────────────────────────────
    for k, v in enumerate(rr.KVV):
        lines.append(f"KVV_{k}={v:.17e}")
    for k, v in enumerate(rr.KVI):
        lines.append(f"KVI_{k}={v:.17e}")
    for k, v in enumerate(rr.GVV):
        lines.append(f"GVV_{k}={v:.17e}")
    for k, v in enumerate(rr.KHeV):
        lines.append(f"KHeV_{k}={v:.17e}")
    for k, v in enumerate(re.Pr_VAC):
        lines.append(f"Pr_VAC_{k}={v:.17e}")

    # m^{1/3} factors for the K_IclV separable cross-term
    for k in range(Nv):
        m = k + 1
        lines.append(f"m13_{k}={m**(1.0/3.0):.17e}")

    # ── SIA cluster arrays (for k=0..Ni-1) ───────────────────────────────────
    for k, v in enumerate(rr.KII):
        lines.append(f"KII_{k}={v:.17e}")
    for k, v in enumerate(rr.KIV):
        lines.append(f"KIV_{k}={v:.17e}")
    for k, v in enumerate(rr.GII):
        lines.append(f"GII_{k}={v:.17e}")
    for k, v in enumerate(rr.k2_SIA_cluster):
        lines.append(f"k2_SIA_{k}={v:.17e}")
    for k, v in enumerate(re.Pr_SIA):
        lines.append(f"Pr_SIA_{k}={v:.17e}")

    # ── K_IclV separable coefficients ────────────────────────────────────────
    # K_IclV[n-1,m-1] = K_IclV_ns[n-1] + K_IclV_ni[n-1] * m^{1/3}
    # K_IclV_ns[n-1] = C0 * n^{-2/3}  (C0 = 4π·r0·Di/Ω)
    # K_IclV_ni[n-1] = C0 / n
    # For n=1 both are 0 (mono-SIA excluded from cross-recombination).
    d = inp.derived
    Omega = float(d['Omega'])
    Di    = float(d['Di'])
    r0    = float(d['r0'])
    C0    = 4.0 * np.pi * r0 * Di / Omega   # separable constant

    # k=0 → n=1: excluded (K_IclV[0,:] = 0)
    lines.append(f"K_IclV_ns_0=0.00000000000000000e+00")
    lines.append(f"K_IclV_ni_0=0.00000000000000000e+00")
    for k in range(1, Ni):
        n = k + 1
        lines.append(f"K_IclV_ns_{k}={C0 * n**(-2.0/3.0):.17e}")
        lines.append(f"K_IclV_ni_{k}={C0 / n:.17e}")

    # Extend to Ni_max if needed (analytical continuation: same power-law)
    if _Ni_max > Ni:
        for k in range(Ni, _Ni_max):
            n = k + 1
            # KII[k] = KII[0] * (n/1)^{1/3} * ... (Waite): use linear extrapolation
            # from the last computed value at n=Ni  (r_AB ≈ r0*(1+n^{1/3}), D = Di)
            _kii0 = float(rr.KII[0]) / (2.0 ** (1.0/3.0) + 1.0)  # rough extrap base
            # Simpler: use analytic formula for large n
            r_AB = r0 * (1.0 + n**(1.0/3.0))
            _kii = 4.0 * np.pi * r_AB * Di / Omega
            _kiv0 = float(rr.KIV[0]) / (2.0 * r0) * (r0 * (1.0 + 1.0))
            r_AB_iv = r0 * (1.0**(1.0/3.0) + n**(1.0/3.0))
            _kiv = 4.0 * np.pi * r_AB_iv * float(d['Dv']) / Omega
            lines.append(f"KII_{k}={_kii:.17e}")
            lines.append(f"KIV_{k}={_kiv:.17e}")
            lines.append(f"GII_{k}=0.00000000000000000e+00")  # negligible for large n
            _k2 = float(inp.material_params.get('Z_i', 1.05)) * \
                  float(inp.material_params.get('rho_d', 5e14)) * Di / float(n)
            lines.append(f"k2_SIA_{k}={_k2:.17e}")
            lines.append(f"Pr_SIA_{k}=0.00000000000000000e+00")
            lines.append(f"K_IclV_ns_{k}={C0 * n**(-2.0/3.0):.17e}")
            lines.append(f"K_IclV_ni_{k}={C0 / n:.17e}")

    # ── Scalar physics ────────────────────────────────────────────────────────
    import math
    kB = 8.617333262e-5  # eV/K
    T  = float(inp.material_params['T'])
    kBT = kB * T

    nu_He   = float(inp.material_params.get('nu_He',   6.25e12))
    E_m_He  = float(inp.material_params.get('E_m_He',  0.06))
    E_b_HeV = float(inp.material_params.get('E_b_HeV', 2.60))
    beta_He = nu_He * math.exp(-(E_b_HeV + E_m_He) / kBT)

    from Eurofer_CD.py_utils.binding_energies import _He_fit
    delta_He   = float(_He_fit['delta_He'])
    beta_He_exp= float(_He_fit['beta_He'])
    L_He_max   = float(inp.L_He_max)

    lines.extend([
        f"G_He={re.G_He:.17e}",
        f"k2_disl_v={rr.k2_disl_v:.17e}",
        f"k2_disl_i={rr.k2_disl_i:.17e}",
        f"k2_disl_He={rr.k2_disl_He:.17e}",
        f"Cv_eq={d['Cv_eq']:.17e}",
        f"beta_He={beta_He:.17e}",
        f"delta_He={delta_He:.17e}",
        f"beta_He_exp={beta_He_exp:.17e}",
        f"kBT={kBT:.17e}",
        f"L_He_max={L_He_max:.17e}",
    ])

    # ── Initial conditions ────────────────────────────────────────────────────
    y0 = re.get_initial_conditions()
    # y0 layout (Python): [Ci1..Ci_Ni, Cv1..Cv_Nv, C_He]  — same as C++
    for k, v in enumerate(y0):
        lines.append(f"y0_{k}={v:.17e}")
    # Extend to Ni_max if needed
    C_floor_val = float(inp.model_params.get('C_floor', 1e-100))
    for k in range(len(y0), Nv + _Ni_max + 1):
        lines.append(f"y0_{k}={C_floor_val:.17e}")

    # ── Floor ─────────────────────────────────────────────────────────────────
    lines.append(f"C_floor={C_floor_val:.17e}")

    # ── Solver settings ────────────────────────────────────────────────────────
    t_span = solver_config.get('t_span', (1e-8, 1e7))
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

    N_tot = Nv + _Ni_max + 1
    lines.append(f"backend={_backend_map.get(str(method_opts.get('backend', 'cvode')).lower(), 0)}")
    lines.append(f"lmm={_lmm_map.get(str(method_opts.get('lmm', 'bdf')).lower(), 2)}")
    lines.append(f"linsol={_linsol_map.get(str(method_opts.get('linsol', 'dense')).lower(), 0)}")
    lines.append(f"mu={int(method_opts.get('mu', N_tot - 1))}")
    lines.append(f"ml={int(method_opts.get('ml', N_tot - 1))}")
    lines.append(f"max_order={int(method_opts.get('max_order', 0))}")
    _ark_raw = method_opts.get('ark_table', 'ARK548L2SA_DIRK_8_4_5')
    ark_val  = _ark_raw if isinstance(_ark_raw, int) else _ark_table_map.get(str(_ark_raw).upper(), 111)
    lines.append(f"ark_table={ark_val}")

    # ── Dynamic window parameters ──────────────────────────────────────────────
    lines.append(f"window_mode={int(method_opts.get('window_mode', 0))}")
    lines.append(f"window_w0_v={int(method_opts.get('window_w0_v', Nv))}")
    lines.append(f"window_w0_i={int(method_opts.get('window_w0_i', Ni))}")
    lines.append(f"window_C_expand={float(method_opts.get('window_C_expand', 1e-18)):.17e}")
    lines.append(f"window_expand_pad={int(method_opts.get('window_expand_pad', 10))}")
    lines.append(f"window_expand_factor={float(method_opts.get('window_expand_factor', 0.0)):.17e}")
    lines.append(f"window_check_every={int(method_opts.get('window_check_every', 1))}")
    lines.append(f"window_C_contract={float(method_opts.get('window_C_contract', 0.0)):.17e}")
    lines.append(f"window_min_active_i={int(method_opts.get('window_min_active_i', 5))}")
    lines.append(f"window_prec={int(method_opts.get('window_prec', 1))}")
    lines.append(f"window_nuc_guard={float(method_opts.get('window_nuc_guard', 0.0)):.17e}")
    lines.append(f"window_width={int(method_opts.get('window_width', 500))}")
    lines.append(f"window_t_start={float(method_opts.get('window_t_start', 10.0)):.17e}")
    lines.append(f"window_N_thresh={int(method_opts.get('window_N_thresh', 1000))}")
    lines.append(f"window_omp_threads={int(method_opts.get('window_omp_threads', 0))}")
    lines.append(f"window_gmres_maxl={int(method_opts.get('window_gmres_maxl', 20))}")
    lines.append(f"Ni_extend_tol={float(method_opts.get('Ni_extend_tol', 0.0)):.17e}")
    lines.append(f"Ni_extend_margin={int(method_opts.get('Ni_extend_margin', 0))}")

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
    Run the compiled Eurofer_CD C++ solver and return the standard results dict.

    Parameters
    ----------
    sim           : EuroferCDSimulation – fully initialised
    solver_config : dict  – same keys as SOLVER_CONFIG in the notebook
    base_dir      : Path or None  – Eurofer_CD/ root; auto-detected if None

    Returns
    -------
    dict (same format as post_process.calculate_derived_quantities) or None
    """
    from Eurofer_CD.py_utils.post_process import calculate_derived_quantities

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

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt',
                                     delete=False, prefix='eurofer_cd_params_') as tf:
        param_path = tf.name

    bin_path = (param_path[:-4] + '.bin' if param_path.endswith('.txt')
                else param_path + '.bin')

    method_opts = solver_config.get('solver_method', {})
    Ni  = int(sim.input_data.model_params['Ni'])
    Nv  = int(sim.input_data.model_params['Nv'])
    _Ni_max = int(method_opts.get('Ni_max', Ni))
    N   = _Ni_max + Nv + 1   # total output columns per row

    try:
        write_param_file(sim, solver_config, param_path)
        print(f"Running C++ solver ({exe_path.name})  …  "
              f"(Ni={Ni}, Nv={Nv}, N_EQ={N}, he_mode='{sim.rate_equations.he_mode}')")
        _proc = subprocess.Popen(
            [str(exe_path), f'--param_file={param_path}'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # Forward stderr line-by-line so progress prints appear live in Jupyter
        def _fwd_stderr():
            for raw in _proc.stderr:
                sys.stderr.write(raw.decode('utf-8', errors='replace'))
                sys.stderr.flush()
        _t = threading.Thread(target=_fwd_stderr, daemon=True)
        _t.start()
        _stdout_data = _proc.stdout.read()
        _proc.wait()
        _t.join()

        class _ProcResult:
            returncode = _proc.returncode
            stdout     = _stdout_data
        proc = _ProcResult()
    finally:
        try:
            os.unlink(param_path)
        except OSError:
            pass

    if proc.returncode != 0:
        print(f"❌ C++ solver failed (exit code {proc.returncode})")
        try:
            os.unlink(bin_path)
        except OSError:
            pass
        return None

    # stderr was streamed live to sys.stderr above; nothing to print here

    # Read binary output
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
        stdout_text = (proc.stdout.decode('utf-8', errors='replace')
                       if isinstance(proc.stdout, bytes) else proc.stdout)
        sol_arr = _parse_stdout(stdout_text, N)

    if sol_arr.shape[0] == 0:
        print("❌ C++ solver produced no parseable output")
        return None

    n_pts = sol_arr.shape[0]
    print(f"✓ C++ solver completed — {n_pts} time points")

    t = sol_arr[:, 0]
    y = sol_arr[:, 1:].T   # (N, n_pts)

    # Trim back to Ni+Nv+1 if Ni_max > Ni
    _N_trim = sim.rate_equations.N   # Ni + Nv + 1
    if y.shape[0] > _N_trim:
        y = y[:_N_trim, :]

    results = calculate_derived_quantities(t, y, sim.input_data, sim.rate_equations)

    # Build metadata label
    backend_name = str(method_opts.get('backend', 'cvode')).lower()
    linsol_label = str(method_opts.get('linsol', 'dense')).upper()
    win_mode     = int(method_opts.get('window_mode', 0))

    if backend_name == 'arkode':
        ark_tbl     = str(method_opts.get('ark_table', 'ARK548L2SA_DIRK_8_4_5'))
        backend_str = f'C++ SUNDIALS ARKODE ARKStep DIRK ({ark_tbl}) / {linsol_label}'
    elif win_mode == 0:
        lmm_label   = str(method_opts.get('lmm', 'bdf')).upper()
        backend_str = f'C++ SUNDIALS CVODE {lmm_label} FULL / {linsol_label}'
    elif win_mode == 4:
        n_thr = int(method_opts.get('window_omp_threads', 0))
        thr_label = f'{n_thr} threads' if n_thr > 0 else 'OMP_NUM_THREADS'
        backend_str = f'C++ SUNDIALS CVODE BDF WINDOW-IV (OpenMP, {thr_label}) / GMRES'
    else:
        phase_names = {1: 'WINDOW-I', 2: 'WINDOW-II', 3: 'WINDOW-III'}
        w0_i = int(method_opts.get('window_w0_i', Ni))
        backend_str = (f'C++ SUNDIALS CVODE BDF {phase_names.get(win_mode,"WINDOW")}'
                       f' / GMRES (w0_i={w0_i})')

    results['metadata'] = {
        'solver_stats': {
            'success':       True,
            'message':       backend_str,
            'n_time_points': n_pts,
        },
    }
    print("✓ Results processing complete!")
    return results
