"""
cpp_bridge.py — Python → C++ solver bridge for RadCluster_1_0.

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
- Bin-moment parameters: I_bin, V_bin, i_discrete, v_discrete, r_ratio (computed)
- i_mobile, v_mobile: mobility cutoffs
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


def write_param_file(sim, solver_config, path, y0_override=None):
    """
    Write all solver parameters to a text file (key=value, one per line).

    Parameters
    ----------
    sim           : RadClusterSimulation — fully initialised
    solver_config : dict  — t_span, rtol, atol, solver_method, etc.
    path          : str or Path
    y0_override   : ndarray or None — if provided, use these initial conditions
                    instead of re_obj.get_initial_conditions().  Used by
                    adaptive continuation to resume from a mid-run state.
    """
    inp  = sim.input_data
    rr   = sim.reaction_rates
    re_obj = sim.rate_equations

    d      = inp.derived
    method = solver_config.get('solver_method', {})

    I  = inp.I
    V  = inp.V

    lines = []

    # ── Cluster size limits ───────────────────────────────────────────────────
    lines.append(f"I={I}")
    lines.append(f"V={V}")
    # Legacy keys for backward compat with older C++ builds
    lines.append(f"N={I}")
    lines.append(f"M={V}")
    lines.append(f"Ni={I}")
    lines.append(f"Nv={V}")
    lines.append(f"Ni_max={I}")

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
    lines.append(f"i_mobile={d['i_mobile']}")
    lines.append(f"v_mobile={d['v_mobile']}")
    # Legacy keys
    lines.append(f"n_max_i={d['i_mobile']}")
    lines.append(f"m_max_v={d['v_mobile']}")

    # ── Boundary flux option ─────────────────────────────────────────────────
    # 0 = absorption (open boundary, default), 1 = reflection (closed boundary)
    bf = d.get('boundary_flux', 'absorption')
    lines.append(f"boundary_flux={1 if bf == 'reflection' else 0}")

    # ── Vacancy cluster rate arrays (0-indexed, size V) ───────────────────────
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
    for k in range(V):
        m = k + 1
        lines.append(f"m13_{k}={m**(1.0/3.0):.17e}")

    # ── SIA cluster rate arrays (0-indexed, size I) ───────────────────────────
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
    for k in range(1, I):
        n = k + 1
        lines.append(f"K_IclV_ns_{k}={C0 * n**(-2.0/3.0):.17e}")
        lines.append(f"K_IclV_ni_{k}={C0 / n:.17e}")

    # ── Mobile cluster effective 3D diffusivities (for coalescence) ─────────
    for k, v in enumerate(rr.D_SIA_eff):
        lines.append(f"D_SIA_eff_{k}={v:.17e}")
    for k, v in enumerate(rr.D_VAC_eff):
        lines.append(f"D_VAC_eff_{k}={v:.17e}")
    lines.append(f"A_sph_inv_O23={rr.A_sph_inv_O23:.17e}")
    lines.append(f"A_loop_inv_O23={rr.A_loop_inv_O23:.17e}")
    Z_i_loop = float(inp.reactions.get('Z_i', 1.05))  # loop bias = dislocation bias Z_i
    lines.append(f"Z_i_loop={Z_i_loop:.17e}")
    Z_ii = float(inp.reactions.get('Z_ii', 1.0))
    lines.append(f"Z_ii={Z_ii:.17e}")

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
        f"K_3D_cav_pref={rr.K_3D_cav_pref:.17e}",
    ])

    # ── Bin-moment parameters ─────────────────────────────────────────────────
    I_bin = getattr(re_obj, 'I_bin', getattr(re_obj, 'K', 0))
    V_bin = getattr(re_obj, 'V_bin', getattr(re_obj, 'K_v', 0))
    i_discrete = getattr(re_obj, 'i_discrete', getattr(re_obj, 'n1', 1))
    v_discrete = getattr(re_obj, 'v_discrete', 1)
    r_ratio = getattr(re_obj, 'r', 2.0)

    lines.append(f"I_bin={I_bin}")
    lines.append(f"V_bin={V_bin}")
    lines.append(f"i_discrete={i_discrete}")
    lines.append(f"v_discrete={v_discrete}")
    lines.append(f"r_ratio={r_ratio:.17e}")
    # Legacy keys for backward compat
    lines.append(f"K_bins={I_bin}")
    lines.append(f"K_v_bins={V_bin}")
    lines.append(f"n1_bin={i_discrete}")

    # Shape function: constant=0, linear=1, lognormal=2
    _sf_map = {'constant': 0, 'linear': 1, 'lognormal': 2}
    sf = getattr(re_obj, 'shape_function', 'linear')
    n_mom = getattr(re_obj, 'n_mom', 2)
    lines.append(f"shape_function={_sf_map.get(sf, 1)}")
    lines.append(f"n_mom={n_mom}")

    # ── He mode ───────────────────────────────────────────────────────────────
    he_mode_int = 0 if getattr(re_obj, 'he_mode', 'case2') == 'case2' else 1
    lines.append(f"he_mode={he_mode_int}")

    # ── Initial conditions ────────────────────────────────────────────────────
    y0 = y0_override if y0_override is not None else re_obj.get_initial_conditions()
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

    N_tot = I + V + 1
    lines.append(f"backend={_backend_map.get(str(method.get('backend','cvode')).lower(), 0)}")
    lines.append(f"lmm={_lmm_map.get(str(method.get('lmm','bdf')).lower(), 2)}")
    lines.append(f"linsol={_linsol_map.get(str(method.get('linsol','dense')).lower(), 0)}")
    lines.append(f"mu={int(method.get('mu', N_tot - 1))}")
    lines.append(f"ml={int(method.get('ml', N_tot - 1))}")
    lines.append(f"max_order={int(method.get('max_order', 0))}")
    lines.append(f"hmin={float(method.get('hmin', 0.0)):.17e}")
    lines.append(f"ark_table={int(method.get('ark_table', 111))}")

    # ── Dynamic window parameters ──────────────────────────────────────────────
    lines.append(f"window_w0_v={V}")
    lines.append(f"window_w0_i={int(method.get('window_w0_i', I))}")
    lines.append(f"window_C_expand={float(method.get('window_C_expand', 1e-18)):.17e}")
    lines.append(f"window_expand_pad={int(method.get('window_expand_pad', 10))}")
    # VAC window parameters — default to the SIA values if not explicitly set
    lines.append(f"window_C_expand_v={float(method.get('window_C_expand_v', method.get('window_C_expand', 1e-18))):.17e}")
    lines.append(f"window_expand_pad_v={int(method.get('window_expand_pad_v', method.get('window_expand_pad', 10)))}")
    lines.append(f"window_expand_factor={float(method.get('window_expand_factor', 0.0)):.17e}")
    lines.append(f"window_check_every={int(method.get('window_check_every', 1))}")
    lines.append(f"window_C_contract={float(method.get('window_C_contract', 0.0)):.17e}")
    lines.append(f"window_min_active_i={int(method.get('window_min_active_i', 5))}")
    lines.append(f"window_prec={int(method.get('window_prec', 1))}")
    lines.append(f"window_nuc_guard={float(method.get('window_nuc_guard', 0.0)):.17e}")
    lines.append(f"window_width={int(method.get('window_width', 500))}")
    lines.append(f"window_t_start={float(method.get('window_t_start', 10.0)):.17e}")
    lines.append(f"window_N_thresh={int(method.get('window_N_thresh', 1000))}")
    lines.append(f"window_gmres_maxl={int(method.get('window_gmres_maxl', 20))}")
    lines.append(f"Ni_extend_tol=0.00000000000000000e+00")
    lines.append(f"Ni_extend_margin=0")

    # ── Woodbury preconditioner parameters ─────────────────────────────────────
    # prec_type: 0=Jacobi (legacy), 1=Woodbury (bordered-banded, default for GMRES)
    linsol_int = _linsol_map.get(str(method.get('linsol','dense')).lower(), 0)
    window_mode_int = int(method.get('window_mode', 0))
    # Woodbury only for full solver (window_mode==0) with GMRES — the sliding
    # window already keeps the active system small enough for Jacobi+GMRES.
    prec_type_default = 1 if (linsol_int == 2 and window_mode_int == 0) else 0
    # User-facing name takes priority over the legacy integer.  Accept a few
    # case/spelling variants so the notebook config stays forgiving.
    _prec_name_map = {
        'jacobi':    0,
        'woodbury':  1,
        'woodburry': 1,   # common typo
    }
    if 'preconditioner' in method:
        key = str(method['preconditioner']).strip().lower()
        if key not in _prec_name_map:
            raise ValueError(
                f"Unknown preconditioner='{method['preconditioner']}'. "
                f"Use 'Jacobi' or 'Woodbury'.")
        prec_type_value = _prec_name_map[key]
    else:
        prec_type_value = int(method.get('prec_type', prec_type_default))
    lines.append(f"prec_type={prec_type_value}")
    # prec_bw: half-bandwidth (auto from mobility cutoffs)
    prec_bw_default = max(2 * d['i_mobile'], 2 * d['v_mobile']) + 1
    lines.append(f"prec_bw={int(method.get('prec_bw', prec_bw_default))}")
    # prec_rank: number of mobile species forming the dense border
    prec_rank_default = d['i_mobile'] + d['v_mobile']
    lines.append(f"prec_rank={int(method.get('prec_rank', prec_rank_default))}")

    verbose = 1 if solver_config.get('_verbose', False) else 0
    lines.append(f"verbose={verbose}")

    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


# ── Output parsing ─────────────────────────────────────────────────────────────

def _parse_stdout(text, N_eq):
    rows = []
    for line in text.strip().splitlines():
        parts = line.split()
        if len(parts) == 1 + N_eq:
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                pass
    return np.array(rows) if rows else np.empty((0, 1 + N_eq))


# ── Diagnostic line parser ────────────────────────────────────────────────────

def _parse_kv_line(line):
    """Parse a C++ diagnostic stderr line of the form  key=value  key=value ...
    Returns a dict of {str: float} for every key=value token found."""
    out = {}
    for token in line.split():
        if '=' in token:
            k, _, v = token.partition('=')
            try:
                out[k] = float(v)
            except ValueError:
                pass
    return out


def _make_stderr_handler(progress_callback, info_out=None):
    """
    Return a callable suitable for use as a daemon-thread target that reads
    proc.stderr line by line.

    When progress_callback is None  → forward each line to sys.stderr verbatim.
    When progress_callback is given → also parse [diag] / [ci5_rates] /
    [cv5_rates] lines and call progress_callback(row_dict) once all three
    lines for a given time step have been received.

    When info_out is a dict, the handler stores any solver-emitted metadata
    lines (currently `[OpenMP_threads] N`) into it under a stable key.

    The row_dict passed to the callback contains atom-fraction concentrations
    and atom-fraction/s rates exactly as the C++ solver computed them:
      t, c_i1, c_v1, c_i2, c_v2, c_i5, c_v5, Q_tot, SIA_tot, VAC_tot
      ci5_prod, ci5_emit_in, ci5_emit_out, ci5_grow_in, ci5_grow_out,
      ci5_shrink_in, ci5_shrink_out, ci5_1D_loss, ci5_sink
      cv5_prod, cv5_emit_in, cv5_emit_out, cv5_grow_in, cv5_grow_out,
      cv5_shrink_in, cv5_shrink_out, cv5_1D_loss, cv5_sink
    """
    pending = {}     # accumulates fields for the current time step
    lock    = threading.Lock()

    def _flush():
        if pending:
            try:
                progress_callback(dict(pending))
            except Exception:
                pass
            pending.clear()

    _diag_prefixes = ('[diag]', '[ci5_rates]', '[cv5_rates]')

    def _thread(proc_stderr):
        for raw in proc_stderr:
            line = raw.decode('utf-8', errors='replace')
            stripped = line.strip()

            if info_out is not None and stripped.startswith('[OpenMP_threads]'):
                try:
                    info_out['omp_threads_used'] = int(stripped.split()[-1])
                except (ValueError, IndexError):
                    pass

            if info_out is not None and stripped.startswith('[stats]'):
                try:
                    kv = _parse_kv_line(stripped[len('[stats]'):])
                    info_out['solver_stats_final'] = {
                        k: (int(v) if isinstance(v, float) and v.is_integer() else v)
                        for k, v in kv.items()
                    }
                except Exception:
                    pass

            # Only echo non-diagnostic lines to stderr; diagnostic lines
            # are consumed silently by the progress_callback parser below.
            if not any(stripped.startswith(p) for p in _diag_prefixes):
                sys.stderr.write(line)
                sys.stderr.flush()

            if progress_callback is None:
                continue
            with lock:
                if stripped.startswith('[diag]'):
                    _flush()   # emit previous time step before starting new one
                    kv = _parse_kv_line(stripped[len('[diag]'):])
                    pending.update(kv)
                elif stripped.startswith('[ci5_rates]'):
                    kv = _parse_kv_line(stripped[len('[ci5_rates]'):])
                    pending.update({f'ci5_{k}': v for k, v in kv.items()})
                elif stripped.startswith('[cv5_rates]'):
                    kv = _parse_kv_line(stripped[len('[cv5_rates]'):])
                    pending.update({f'cv5_{k}': v for k, v in kv.items()})
                elif stripped.startswith('Done:'):
                    _flush()   # flush the last time step

    return _thread


# ── Main entry point ──────────────────────────────────────────────────────────

def run_cpp_solver(sim, solver_config, base_dir=None, progress_callback=None,
                   timeout_s=None, y0_override=None):
    """
    Run the RadCluster_1_0 C++ solver and return the standard results dict.

    Parameters
    ----------
    sim               : RadClusterSimulation
    solver_config     : dict
    base_dir          : Path or None
    progress_callback : callable or None
        If provided, called once per output time step with a dict of
        concentrations and rate-breakdown values (all in atom fraction / s).
        The solver's verbose mode is automatically enabled.
    timeout_s         : float or None
        Maximum wall-clock seconds to allow for the C++ solver.  If exceeded,
        the process is killed and None is returned.
    y0_override       : ndarray or None
        Custom initial conditions for adaptive continuation runs.

    Returns
    -------
    dict or None
    """
    from .post_process import calculate_derived_quantities

    if base_dir is None:
        base_dir = Path(__file__).resolve().parent.parent

    # If a callback is requested, enable C++ verbose output automatically
    if progress_callback is not None:
        solver_config = dict(solver_config)
        solver_config['_verbose'] = True

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

    proc = None
    try:
        write_param_file(sim, solver_config, param_path, y0_override=y0_override)
        print(f"C++ solver: {exe_path.name}  N_eq={N_tot}"
              f"  solver_mode='{sim.input_data.solver_mode}'"
              f"  physics='{sim.input_data.physics_option}'")

        proc = subprocess.Popen(
            [str(exe_path), f'--param_file={param_path}'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        solver_info = {}
        stderr_fn = _make_stderr_handler(progress_callback, info_out=solver_info)
        t_fwd = threading.Thread(target=stderr_fn, args=(proc.stderr,), daemon=True)
        t_fwd.start()
        partial = False
        stdout_data = b''
        try:
            stdout_data = proc.stdout.read()
            proc.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
            print(f"C++ solver killed after {timeout_s}s timeout — parsing partial output")
            partial = True
        t_fwd.join(timeout=2)
    except KeyboardInterrupt:
        print("\n*** Ctrl+C — terminating C++ solver, parsing partial output ***")
        partial = True
        if proc is not None:
            try:
                proc.terminate()
                proc.wait(timeout=5)
            except (subprocess.TimeoutExpired, Exception):
                proc.kill()
                proc.wait()
    finally:
        try:
            os.unlink(param_path)
        except OSError:
            pass

    if not partial and proc.returncode != 0:
        print(f"C++ solver failed (exit code {proc.returncode})")

    # Parse binary output (works for both complete and partial/interrupted runs —
    # the C++ solver flushes each row to the .bin file as it's computed)
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
        if stdout_data:
            text    = stdout_data.decode('utf-8', errors='replace')
            sol_arr = _parse_stdout(text, N_tot)

    if sol_arr is None or sol_arr.shape[0] == 0:
        print("C++ solver produced no parseable output")
        return None

    status = "partial" if partial else "completed"
    n_pts = sol_arr.shape[0]
    print(f"C++ solver {status} — {n_pts} time points")

    t = sol_arr[:, 0]
    y = sol_arr[:, 1:].T   # (N_tot, n_pts)

    results = calculate_derived_quantities(t, y, sim.input_data, re_obj)
    results['y'] = y   # raw ODE state [N_eq, n_pts] in atom fraction

    sm     = sim.input_data.solver_mode
    po     = sim.input_data.physics_option
    linsol = str(solver_config.get('solver_method', {}).get('linsol', 'dense')).upper()
    results['metadata'] = {
        'solver_stats': {
            'success':       True,
            'message':       f'C++ CVODE BDF {sm}/{po} / {linsol}',
            'n_time_points': n_pts,
        },
        'omp_threads_used':   solver_info.get('omp_threads_used'),
        'solver_stats_final': solver_info.get('solver_stats_final'),
    }
    print("Results processing complete.")
    return results
