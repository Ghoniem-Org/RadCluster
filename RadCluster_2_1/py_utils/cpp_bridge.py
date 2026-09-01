"""
cpp_bridge.py -- Python -> C++ solver bridge for RadCluster_2_1.

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
- solver_mode:     integer (0=full_system, 4=active_window)
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
import signal as _signal
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

# How long to wait for the C++ solver to finish its current step and emit
# its final [stats] line after we ask it to stop.  The solver checks the
# interrupt flag once per output time-point, so the upper bound is roughly
# one CVode(t_out) call -- typically seconds, but can be tens of seconds
# when a single step has been struggling to converge.
GRACEFUL_SHUTDOWN_TIMEOUT_S = 60.0


def _pin_child_to_cpu_group(pid):
    """Keep a solver.exe on a specific Windows processor group.

    A SINGLE PIN AT SPAWN DOES NOT HOLD, which is the whole reason this is a
    retry loop and not one call.  Measured 2026-08-18: pinning 5 already-running
    solvers moved the load immediately and permanently (group 1: 3.8 % ->
    14.2 %), but the identical call issued straight after Popen left all 18
    solvers of a fresh stage homed to group 0.  The child's own startup -- CRT
    and the OpenMP runtime, which is group-aware -- re-homes the thread after
    CreateProcess returns, so a pin placed before that is simply overwritten.

    So pin, then re-pin over the first seconds of the row, from a daemon thread
    that never blocks the row.  Rows run for hours; a few seconds on the wrong
    socket costs nothing measurable.
    """
    group = os.environ.get("RADCLUSTER_CPU_GROUP", "").strip()
    if not group or sys.platform != "win32":
        return
    def _keep():
        for delay in (0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0):
            if delay:
                time.sleep(delay)
            try:
                if _pin_once(pid, int(group)) == 0:
                    return                    # process gone -- row finished
            except Exception:
                return
    t = threading.Thread(target=_keep, daemon=True)
    t.start()


def _pin_once(pid, group):
    """One pinning pass.  Returns the number of threads moved/confirmed, 0 if
    the process no longer exists.

    NO-OP unless RADCLUSTER_CPU_GROUP is set, so every host that does not opt in
    behaves exactly as before.  Never raises: failing to pin must cost a little
    throughput, never a row.

    WHY IT EXISTS.  A Windows host with more than 64 logical CPUs is split into
    processor GROUPS, and the scheduler homes threads to the group their process
    was assigned at creation.  On MATRIX-PC2 (2x Xeon Gold 6230, two groups of 40
    logical) that is group 0 for every solver, whatever the parent did: measured
    2026-08-18 under a live 20-row stage, each solver.exe had exactly ONE thread,
    all 20 homed to group 0, socket 0 saturated and socket 1 at 3 %.  The run_
    ensemble pool workers themselves DO spread across both groups -- it is only
    the solver, where all the compute is, that piles up.

    The fix is one call per row.  Measured immediately after moving 5 of 20 live
    solver threads to group 1:  group 0  50.7 % -> 37.4 %,  group 1  3.8 % ->
    14.2 %.

    NOT a numerical knob: affinity changes where a thread runs, never what it
    computes.  Rows produced with and without it are directly comparable.
    """
    try:
        import ctypes
        from ctypes import wintypes

        class _THREADENTRY32(ctypes.Structure):
            _fields_ = [("dwSize", wintypes.DWORD), ("cntUsage", wintypes.DWORD),
                        ("th32ThreadID", wintypes.DWORD),
                        ("th32OwnerProcessID", wintypes.DWORD),
                        ("tpBasePri", ctypes.c_long),
                        ("tpDeltaPri", ctypes.c_long),
                        ("dwFlags", wintypes.DWORD)]

        class _GROUP_AFFINITY(ctypes.Structure):
            _fields_ = [("Mask", ctypes.c_size_t), ("Group", ctypes.c_ushort),
                        ("Reserved", ctypes.c_ushort * 3)]

        k = ctypes.WinDLL("kernel32", use_last_error=True)
        k.OpenThread.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        k.OpenThread.restype = wintypes.HANDLE
        k.SetThreadGroupAffinity.argtypes = [wintypes.HANDLE,
                                             ctypes.POINTER(_GROUP_AFFINITY),
                                             ctypes.POINTER(_GROUP_AFFINITY)]
        k.GetActiveProcessorCount.argtypes = [ctypes.c_ushort]
        k.GetActiveProcessorCount.restype = wintypes.DWORD

        g = int(group)
        n_cpu = int(k.GetActiveProcessorCount(ctypes.c_ushort(g)))
        if n_cpu <= 0:
            return 0
        mask = (1 << n_cpu) - 1

        snap = k.CreateToolhelp32Snapshot(0x00000004, 0)   # TH32CS_SNAPTHREAD
        te = _THREADENTRY32()
        te.dwSize = ctypes.sizeof(_THREADENTRY32)
        moved = 0
        ok = k.Thread32First(snap, ctypes.byref(te))
        while ok:
            if te.th32OwnerProcessID == pid:
                # THREAD_SET_INFORMATION | THREAD_QUERY_LIMITED_INFORMATION
                h = k.OpenThread(0x0020 | 0x0800, False, te.th32ThreadID)
                if h:
                    ga = _GROUP_AFFINITY(Mask=mask, Group=g)
                    if k.SetThreadGroupAffinity(h, ctypes.byref(ga), None):
                        moved += 1
                    k.CloseHandle(h)
            ok = k.Thread32Next(snap, ctypes.byref(te))
        k.CloseHandle(snap)
        if moved:
            logging.debug("pinned solver pid %s (%d thread(s)) to CPU group %d",
                          pid, moved, g)
        return moved
    except Exception as exc:                       # never lose a row over this
        logging.debug("CPU-group pin skipped for pid %s: %s", pid, exc)
        return 0


def _send_graceful_interrupt(proc):
    """Ask the C++ subprocess to finish its current integration step and
    exit cleanly via its own signal handler.  On Windows this requires the
    process to have been launched with CREATE_NEW_PROCESS_GROUP."""
    if proc is None or proc.poll() is not None:
        return
    try:
        if sys.platform == 'win32':
            proc.send_signal(_signal.CTRL_BREAK_EVENT)
        else:
            proc.send_signal(_signal.SIGTERM)
    except Exception:
        try:
            proc.terminate()
        except Exception:
            pass

#: Poll interval while waiting on the solver.  Short enough that an abort is
#: acted on promptly, long enough to be free next to a multi-hour integration.
_ABORT_POLL_S = 5.0


def _wait_with_abort_file(proc, timeout_s):
    """proc.wait(timeout=timeout_s), but also abortable from OUTSIDE the process.

    WHY.  A detached campaign worker has no console, so Ctrl+C cannot reach it,
    and on Windows Stop-Process is TerminateProcess -- uncatchable.  Before this,
    the only way to stop such a run was to kill it, which threw away every
    in-flight row: the solver never got its interrupt, so no partial trajectory
    was finalized and nothing was written.  That is how a 12-row batch was lost
    on 2026-08-11 after ~7 h of integration.

    With this, dropping the abort file named by RADCLUSTER_ABORT_FILE asks the
    solver to finish its CURRENT integration step and flush -- the same graceful
    path the wall-clock budget already used, and the reason budget-cut rows keep
    their dose ladder instead of vanishing.  The row comes back marked partial,
    with every rung it reached intact.

    Raises subprocess.TimeoutExpired on the wall-clock budget, exactly as
    proc.wait did, so the caller's existing handling is unchanged.  With no
    abort file configured this is proc.wait plus a 5 s poll, i.e. a no-op.
    """
    abort = os.environ.get('RADCLUSTER_ABORT_FILE')
    if not abort:
        return proc.wait(timeout=timeout_s)
    deadline = None if not timeout_s else (time.monotonic() + float(timeout_s))
    while True:
        remaining = None if deadline is None else deadline - time.monotonic()
        if remaining is not None and remaining <= 0:
            raise subprocess.TimeoutExpired(proc.args, timeout_s)
        try:
            return proc.wait(timeout=min(_ABORT_POLL_S,
                                         remaining if remaining is not None
                                         else _ABORT_POLL_S))
        except subprocess.TimeoutExpired:
            pass
        if os.path.exists(abort):
            print(f"    abort file present ({abort}) -- asking the solver to "
                  f"finalize at the dose reached so far and flush.")
            raise subprocess.TimeoutExpired(proc.args, timeout_s)


_SOLVER_MODE_MAP = {
    'full_system':   0,
    'active_window': 4,
    # Backward-compat aliases (legacy names)
    'cpp_full':       0,
    'sliding_OpenMP': 4,
}
_PHYSICS_OPTION_MAP = {
    'full_CD_fission':       0,
    'full_CD_fusion':        1,
    'bin_moment_CD_fission': 2,
    'bin_moment_CD_fusion':  3,
}


def sia100_block_length(sim, solver_config):
    """Length of the appended ⟨100⟩ SIA block in the C++ state vector.

    The ⟨100⟩ block carries the SAME size reduction as the ½⟨111⟩ population:
    full per-size (length I) in discrete modes, and the discrete-prefix +
    bin-moment vector (i_discrete + n_mom·I_bin) in bin_moment modes.  Zero when
    loop conversion is off.

    Shared by ``write_param_file`` (which must emit a y0 of the full width) and
    ``run_cpp_solver`` (which splits the block back off).  Keeping one rule in
    one place is deliberate: when the two disagreed, the writer emitted a short
    y0 and the reader still found a block, so the mismatch was invisible.
    """
    if not int(solver_config.get(
            'loop_conversion', sim.input_data.reactions.get('loop_conversion', 0))):
        return 0
    re_obj = sim.rate_equations
    if hasattr(re_obj, 'bins'):        # BinMomentRateEquations
        return (int(getattr(re_obj, 'i_discrete', 0))
                + int(getattr(re_obj, 'n_mom', 2))
                * int(getattr(re_obj, 'I_bin', len(getattr(re_obj, 'bins', [])))))
    return int(sim.input_data.I)


def write_param_file(sim, solver_config, path, y0_override=None):
    """
    Write all solver parameters to a text file (key=value, one per line).

    Parameters
    ----------
    sim           : RadClusterSimulation -- fully initialised
    solver_config : dict  -- t_span, rtol, atol, solver_method, etc.
    path          : str or Path
    y0_override   : ndarray or None -- if provided, use these initial conditions
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
    lines.append(f"physics_option_int={po_int}")
    lines.append(f"window_mode={sm_int}")   # solver mode selector (0/4)

    # ── Geometric rate constant prefactors (Eq. 128) ─────────────────────────
    lines.append(f"A_sph={d['A_sph']:.17e}")
    lines.append(f"A_loop={d['A_loop']:.17e}")
    lines.append(f"B_rot={d['B_rot']:.17e}")
    lines.append(f"L_hat={d['L_hat']:.17e}")

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
    # Per-size SIA loop bias Z_i^loop(n).  K_ii_coal recomputes the SIA capture
    # rate inside the C++ from a SCALAR P.Z_i_loop, so a size-dependent bias
    # written only into K_SIA_grow (KII) would be silently ignored -- KII is
    # consumed only for the top-bin term.  Shipping the array is what makes the
    # Wolfer bias reach the growth path.  Absent -> C++ falls back to the scalar.
    _zarr = getattr(rr, "Z_i_loop_arr", None)
    if _zarr is not None:
        for k, v in enumerate(_zarr):
            lines.append(f"Z_i_loop_arr_{k}={v:.17e}")
    for k, v in enumerate(rr.G_SIA):
        lines.append(f"GII_{k}={v:.17e}")
    for k, v in enumerate(rr.k2_SIA):
        lines.append(f"k2_SIA_{k}={v:.17e}")
    for k, v in enumerate(re_obj.Pr_SIA):
        lines.append(f"Pr_SIA_{k}={v:.17e}")

    # ── 1D glide prefactors K_1D_pref[n-1] (Eq. 141) ─────────────────────────
    for k, v in enumerate(rr.K_1D_pref):
        lines.append(f"K_1D_pref_{k}={v:.17e}")

    # ── Mobile cluster effective 3D diffusivities (for coalescence) ─────────
    for k, v in enumerate(rr.D_SIA_eff):
        lines.append(f"D_SIA_eff_{k}={v:.17e}")
    for k, v in enumerate(rr.D_VAC_eff):
        lines.append(f"D_VAC_eff_{k}={v:.17e}")
    # ½⟨111⟩ loop coarsening: glide law continued past i_mobile, used ONLY by
    # the loop–loop coalescence edge (see ReactionRates.D_loop_coal).
    _lc = int(getattr(rr, 'loop_coal', 0))
    lines.append(f"loop_coal={_lc}")
    if _lc:
        for k, v in enumerate(rr.D_loop_coal):
            lines.append(f"D_loop_coal_{k}={v:.17e}")
    lines.append(f"A_sph_inv_O23={rr.A_sph_inv_O23:.17e}")
    lines.append(f"A_loop_inv_O23={rr.A_loop_inv_O23:.17e}")
    # Loop SIA bias — own key, falling back to the network bias Z_i when the
    # workbook has no row (keeps pre-existing workbooks bit-identical).
    _Z_i = float(inp.reactions.get('Z_i', 1.05))
    Z_i_loop = float(inp.reactions.get('Z_i_loop', _Z_i) or _Z_i)
    lines.append(f"Z_i_loop={Z_i_loop:.17e}")
    Z_ii = float(inp.reactions.get('Z_ii', 1.0))
    lines.append(f"Z_ii={Z_ii:.17e}")

    # ── ½⟨111⟩ → ⟨100⟩ loop conversion (optional; appended SIA100 block) ──────
    # Enabled via solver_config['loop_conversion'] or reactions['loop_conversion'].
    # The C++ appends a c_i100 block of length I at the end of the state vector
    # (N_eq += I) and computes the 2-D junction/absorption kernels on the fly.
    loop_conv = int(solver_config.get(
        'loop_conversion', inp.reactions.get('loop_conversion', 0)))
    lines.append(f"loop_conversion={loop_conv}")
    if loop_conv:
        lines.append(f"n_loop_min={int(getattr(rr, 'n_loop_min', 4))}")
        # Effective junction threshold, tied to the mobility cutoff so the
        # channel cannot be silently dead when i_mobile < n_j_min_junc.  Must
        # use the SAME rule as the Python kernel or the two RHS disagree.
        from .reaction_rates import effective_n_j_min
        lines.append("n_j_min_junc=%d" % int(effective_n_j_min(
            inp.reactions.get('n_j_min_junc', 30),
            d['i_mobile'],
            inp.reactions.get('n_j_min_frac', 0.6))))
        lines.append(f"phi_max_junc={float(inp.reactions.get('phi_max_junc', 0.5)):.17e}")
        lines.append(f"sigma_s_junc={float(inp.reactions.get('sigma_s_junc', 0.35)):.17e}")
        # Marian two-step success probability P_success(T) — gates junction +
        # absorption (½⟨111⟩→½⟨110⟩→⟨100⟩, Fig. 3); computed by ReactionRates.
        lines.append(f"loop_conv_psuccess={float(getattr(rr, 'conv_psuccess', 1.0)):.17e}")
        # Absorption-only success gate P_success^abs(T) = A_abs·[1 +
        # exp((ΔH₂^abs − ΔH_rev)/k_BT)]^{-1}.  Decoupled from the junction gate
        # because a mobile ≤ i_mobile cluster joining a several-hundred-atom
        # ⟨100⟩ loop is templated by the host character and should not pay the
        # junction's reversion penalty.  Falls back to conv_psuccess so an older
        # ReactionRates object (no such attribute) is byte-identical.
        lines.append("loop_conv_psuccess_abs=%.17e" % float(
            getattr(rr, 'conv_psuccess_abs',
                    getattr(rr, 'conv_psuccess', 1.0))))
        # Absorption-only scale (1-D-glide enhanced capture of ½⟨111⟩ by
        # sessile ⟨100⟩).  Separate from conv_psuccess, which multiplies the
        # junction yield too and so cannot isolate growth from nucleation.
        lines.append("conv_absorb_boost=%.17e" % float(
            inp.reactions.get('absorb_boost_100', 1.0) or 1.0))
        for k, v in enumerate(rr.Gamma_uni):
            lines.append(f"Gamma_uni_{k}={v:.17e}")
        for k, v in enumerate(rr.K_100_grow):
            lines.append(f"K_100_grow_{k}={v:.17e}")
        for k, v in enumerate(rr.K_100_shrink):
            lines.append(f"K_100_shrink_{k}={v:.17e}")
        for k, v in enumerate(rr.G_100):
            lines.append(f"G_100_{k}={v:.17e}")
        # Loop → network-dislocation loss on the sessile ⟨100⟩ block
        # (loop_network_loss.tex).  Zeros when LOOP_NETWORK_LOSS is off, so the
        # OFF path is byte-identical.  The ½⟨111⟩ side needs no extra key: its
        # Λ_n^net is already folded into k2_SIA (additive-sink-strength rule).
        lam100 = getattr(rr, 'Lambda_net_100', None)
        if lam100 is not None:
            for k, v in enumerate(lam100):
                lines.append(f"lambda_net_100_{k}={v:.17e}")

    # ── Scalar physics ────────────────────────────────────────────────────────
    kBT = float(d['kBT'])
    nu_h = float(inp.energetics.get('nu_h', 3.0e12))
    E_m_h = float(inp.energetics.get('E_m_h', 0.06))
    E_b_hV1 = float(inp.energetics.get('E_b_hV_1', 2.30))
    beta_He = nu_h * np.exp(-(E_b_hV1 + E_m_h) / kBT)

    # He-pressure correction to vacancy emission (GVV_eff in rate_kernels.cpp).
    # The kernel multiplies GVV by exp(−δ·β·(ℓ/m)^β / kBT), i.e. the binding
    # energy gains ΔE_b = δ·β·(ℓ/m)^β.  Fit (δ, β) so this power law matches
    # the virial-EOS pressure work P_He·Ω (the +P·Ω term of E_b_bubble in
    # binding_energies.py) at the run temperature, over the physical loading
    # range ℓ/m ∈ [0.05, 2].  δ > 0 ⟹ He suppresses vacancy emission and
    # stabilizes bubbles, consistent with the Python kernels.
    from .binding_energies import _B2, _B3
    Omega = float(d['Omega'])
    _x = np.geomspace(0.05, 2.0, 40)                       # ℓ/m grid
    _PdV = _x * kBT * (1.0 + _B2 * _x / Omega + _B3 * (_x / Omega) ** 2)  # eV
    _slope, _icept = np.polyfit(np.log(_x), np.log(_PdV), 1)
    beta_He_exp = float(_slope)
    delta_He = float(np.exp(_icept) / beta_He_exp)

    lines.extend([
        f"G_He={re_obj.G_He:.17e}",
        f"k2_disl_v={rr.k2_vac_scalar:.17e}",
        f"k2_disl_i={rr.k2_SIA_scalar:.17e}",
        f"k2_disl_He={rr.k2_He_scalar:.17e}",
        f"Cv_eq={d['Cv_eq']:.17e}",
        f"beta_He={beta_He:.17e}",
        f"delta_He={delta_He:.17e}",        # He pressure coeff [eV], fit to virial P·Ω
        f"beta_He_exp={beta_He_exp:.17e}",  # He pressure power-law exponent
        f"kBT={kBT:.17e}",
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

    # ── Explicit integer bin edges ────────────────────────────────────────────
    # Transmit the per-bin lower/upper edges directly so the C++ consumes the
    # exact partition that Python computed.  Re-deriving edges in C++ from
    # r_ratio via std::floor can diverge from numpy.floor over many bins
    # (different SIA/VAC ratios, FP rounding), silently corrupting the last
    # bins.  Each bin k covers integer sizes [bin_lo, bin_hi-1].
    sia_bins = list(getattr(re_obj, 'bins', []) or [])
    vac_bins = list(getattr(re_obj, 'vac_bins', []) or [])
    for k, (nlo, nhi) in enumerate(sia_bins):
        lines.append(f"sia_bin_lo_{k}={int(nlo)}")
        lines.append(f"sia_bin_hi_{k}={int(nhi)}")
    for k, (mlo, mhi) in enumerate(vac_bins):
        lines.append(f"vac_bin_lo_{k}={int(mlo)}")
        lines.append(f"vac_bin_hi_{k}={int(mhi)}")

    # Shape function: constant=0, linear=1, lognormal=2
    _sf_map = {'constant': 0, 'linear': 1, 'lognormal': 2}
    sf = getattr(re_obj, 'shape_function', 'linear')
    lines.append(f"shape_function={_sf_map.get(sf, 1)}")

    # ── He mode ───────────────────────────────────────────────────────────────
    he_mode_int = 0 if getattr(re_obj, 'he_mode', 'case2') == 'case2' else 1
    lines.append(f"he_mode={he_mode_int}")

    # ── Initial conditions ────────────────────────────────────────────────────
    # The C++ reads y0_k for k < N_eq and falls back to 1e-100 for every key it
    # does not find (parameters.h: `optional_param(p, "y0_"+k, 1e-100)`).  A
    # short y0 therefore does not fail — it silently zeroes the tail of the
    # state.  That is exactly how adaptive continuation erased the whole ⟨100⟩
    # block at every segment restart (resume handed back results['y'][:, -1],
    # from which cpp_bridge had already split the block off).
    n_sia100   = sia100_block_length(sim, solver_config)
    n_expected = re_obj.N_eq + n_sia100
    if y0_override is None:
        # Fresh start: get_initial_conditions() covers N_eq only, so the ⟨100⟩
        # block has to be built here.  Give it the SAME shaped IC that
        # get_initial_conditions builds for the ½⟨111⟩ bins — μ₀ = C_floor and
        # μ₁ = C_floor·midpoint — for the reason stated in its docstring: all
        # moments at C_floor gives μ₁/μ₀ = 1, a mean size outside every bin,
        # which the closure cannot reconstruct sensibly.  Leaving the block flat
        # put all 21 ⟨100⟩ bins in that state for the whole early history (the
        # C++ clamps every entry up to C_floor at t_begin and after each output
        # step, so a flat block does not simply grow out of it).
        C_floor_ic = float(inp.reactions.get('C_floor', 1e-15))
        blk = np.full(n_sia100, C_floor_ic)
        if n_sia100 and hasattr(re_obj, 'bins'):
            i_d, Pm = int(re_obj.i_discrete), int(re_obj.n_mom)
            for k, (nlo, nhi) in enumerate(re_obj.bins):
                mid = 0.5 * (nlo + nhi - 1)
                idx = i_d + Pm * k
                if Pm >= 2:
                    blk[idx + 1] = C_floor_ic * mid
                if Pm >= 3:
                    blk[idx + 2] = C_floor_ic * mid * mid
        y0 = np.concatenate([re_obj.get_initial_conditions(), blk])
    else:
        y0 = np.asarray(y0_override, dtype=float)
        # A resume MUST carry the full width; silence here is what hid the bug.
        if y0.size != n_expected:
            raise ValueError(
                f"write_param_file: y0_override has {y0.size} entries but the "
                f"C++ state vector is {n_expected} wide (N_eq={re_obj.N_eq} + "
                f"{n_sia100} for the appended ⟨100⟩ block).  A short y0 would "
                "be silently filled with 1e-100, erasing the ⟨100⟩ population. "
                "Resume via RadClusterSimulation._resume_state / _expand_state, "
                "which re-attach the block.")
    for k, v in enumerate(y0):
        lines.append(f"y0_{k}={v:.17e}")

    C_floor = float(inp.reactions.get('C_floor', 1e-15))
    lines.append(f"C_floor={C_floor:.17e}")

    # Free He mode: 'dynamic'=0 integrates Eq.157; 'quasi_steady_state'=1 uses QSS
    he_kinetics_str = str(inp.reactions.get('he_kinetics', 'dynamic')).lower()
    qss_He_int = 1 if he_kinetics_str == 'quasi_steady_state' else 0
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
    # Backend is fixed (CVODE BDF); only the linear solver is selectable.
    _linsol_map  = {'dense': 0, 'band': 1, 'banded': 1, 'gmres': 2, 'klu': 3}

    N_tot = I + V + 1
    lines.append(f"linsol={_linsol_map.get(str(method.get('linsol','dense')).lower(), 0)}")
    lines.append(f"mu={int(method.get('mu', N_tot - 1))}")
    lines.append(f"ml={int(method.get('ml', N_tot - 1))}")
    lines.append(f"max_order={int(method.get('max_order', 0))}")
    lines.append(f"hmin={float(method.get('hmin', 0.0)):.17e}")
    lines.append(f"hmax={float(method.get('hmax', 0.0)):.17e}")

    # ── Dynamic window parameters ──────────────────────────────────────────────
    # window_width is shared by SIA and VAC; defaults to max(I, V) (full domain).
    lines.append(f"window_width={int(method.get('window_width', max(I, V)))}")
    lines.append(f"concentration_threshold={float(method.get('concentration_threshold', 1e-18)):.17e}")
    lines.append(f"window_pad={int(method.get('window_pad', 10))}")
    # VAC expansion pad defaults to the SIA pad if not explicitly set.
    lines.append(f"window_pad_v={int(method.get('window_pad_v', method.get('window_pad', 10)))}")
    lines.append(f"window_check_every={int(method.get('window_check_every', 1))}")

    # ── Woodbury preconditioner parameters ─────────────────────────────────────
    # prec_type: 0=Jacobi (legacy), 1=Woodbury (bordered-banded, default for GMRES)
    linsol_int = _linsol_map.get(str(method.get('linsol','dense')).lower(), 0)
    window_mode_int = int(method.get('window_mode', 0))
    # Woodbury only for full solver (window_mode==0) with GMRES -- the sliding
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

    When progress_callback is None  -> forward each line to sys.stderr verbatim.
    When progress_callback is given -> also parse [diag] / [ci5_rates] /
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
    Run the RadCluster_2_1 C++ solver and return the standard results dict.

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
    # Loop conversion appends a ⟨100⟩ SIA block at the end of the C++ state
    # vector, so the solver emits N_eq + len(⟨100⟩) columns.  The ⟨100⟩ block
    # carries the SAME size reduction as the ½⟨111⟩ population: full per-size
    # (length I) in discrete modes, and the discrete-prefix + bin-moment vector
    # (i_discrete + n_mom·I_bin) in bin_moment modes.  Parse the wider rows and
    # split the appended block off before post-processing (which expects the
    # original [SIA | VAC | He | conservation] layout).
    _is_bin = hasattr(re_obj, 'bins')   # BinMomentRateEquations
    _n_sia100 = sia100_block_length(sim, solver_config)
    N_tot   += _n_sia100

    proc = None
    try:
        write_param_file(sim, solver_config, param_path, y0_override=y0_override)
        print(f"C++ solver: {exe_path.name}  N_eq={N_tot}"
              f"  solver_mode='{sim.input_data.solver_mode}'"
              f"  physics='{sim.input_data.physics_option}'")

        popen_kwargs = dict(
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # On Windows, give the child its own process group so we can deliver
        # CTRL_BREAK_EVENT to it without also signaling our own console.
        if sys.platform == 'win32':
            popen_kwargs['creationflags'] = subprocess.CREATE_NEW_PROCESS_GROUP

        proc = subprocess.Popen(
            [str(exe_path), f'--param_file={param_path}'],
            **popen_kwargs,
        )
        _pin_child_to_cpu_group(proc.pid)

        solver_info = {}
        stderr_fn = _make_stderr_handler(progress_callback, info_out=solver_info)
        t_fwd = threading.Thread(target=stderr_fn, args=(proc.stderr,), daemon=True)
        t_fwd.start()
        partial = False
        # Drain stdout in a thread so proc.wait(timeout=...) is actually
        # reachable while the child is still running.  If we instead call
        # proc.stdout.read() inline it blocks until the child closes stdout
        # (i.e. exits), defeating the wall-clock cap entirely.
        stdout_buf = []
        def _drain_stdout(stream, buf):
            try:
                buf.append(stream.read())
            except Exception:
                pass
        t_out = threading.Thread(target=_drain_stdout,
                                 args=(proc.stdout, stdout_buf), daemon=True)
        t_out.start()
        stdout_data = b''
        try:
            _wait_with_abort_file(proc, timeout_s)
        except subprocess.TimeoutExpired:
            print(f"C++ solver hit {timeout_s}s timeout -- asking it to finalize gracefully...")
            partial = True
            _send_graceful_interrupt(proc)
            try:
                proc.wait(timeout=GRACEFUL_SHUTDOWN_TIMEOUT_S)
            except subprocess.TimeoutExpired:
                print(f"    Solver did not finalize within "
                      f"{GRACEFUL_SHUTDOWN_TIMEOUT_S:.0f}s -- forcing kill.")
                proc.kill()
                proc.wait()
        t_fwd.join(timeout=5)
        t_out.join(timeout=5)
        stdout_data = b''.join(stdout_buf) if stdout_buf else b''
    except KeyboardInterrupt:
        print("\n*** Ctrl+C -- asking C++ solver to flush and exit gracefully ***")
        partial = True
        _send_graceful_interrupt(proc)
        if proc is not None:
            try:
                # Wait for the child to finish its current output step and
                # emit the final [stats] line; the stdout drain thread keeps
                # collecting whatever lands.  A second Ctrl+C escalates.
                proc.wait(timeout=GRACEFUL_SHUTDOWN_TIMEOUT_S)
            except (subprocess.TimeoutExpired, KeyboardInterrupt):
                print(f"    Solver did not exit within "
                      f"{GRACEFUL_SHUTDOWN_TIMEOUT_S:.0f}s -- forcing kill.")
                try:
                    proc.kill()
                    proc.wait(timeout=5)
                except Exception:
                    pass
            try:
                t_fwd.join(timeout=5)
                t_out.join(timeout=5)
                stdout_data = b''.join(stdout_buf) if stdout_buf else b''
            except Exception:
                pass
    finally:
        try:
            os.unlink(param_path)
        except OSError:
            pass

    if proc is not None and not partial and proc.returncode != 0:
        print(f"C++ solver failed (exit code {proc.returncode})")

    # ── Parse + post-process under a KeyboardInterrupt shield ────────────────
    # The C++ child flushes each row to the .bin file as it's computed, so
    # there is real data to rescue even after Ctrl+C.  We must NOT lose it to
    # a second interrupt during numpy parsing or post-processing -- stash
    # partial results onto `sim` as soon as they exist so the orchestrator
    # (and the notebook fallback) can recover them.
    sim._partial_results = None
    try:
        # Parse binary output (works for both complete and partial runs)
        sol_arr = None
        try:
            try:
                bin_size = os.path.getsize(bin_path)
            except OSError:
                bin_size = -1
            print(f"Reading bin file: {bin_path} ({bin_size} bytes)")
            raw    = np.fromfile(bin_path, dtype=np.float64)
            n_cols = 1 + N_tot
            n_rows = raw.size // n_cols
            # Reshape FIRST -- losing partial data to a downstream print/encode
            # error would be far worse than a missing log line.
            if n_rows > 0:
                sol_arr = raw[:n_rows * n_cols].reshape(n_rows, n_cols)
            try:
                print(f"  parsed {raw.size} doubles -> {n_rows} rows x {n_cols} cols")
            except Exception:
                pass
        except Exception as exc:
            try:
                print(f"  bin parse failed: {type(exc).__name__}: {exc}")
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

        # The C++ solver also stamps interrupted=1 in its [stats] line when its
        # own signal handler tripped -- pick that up too so we don't miss the
        # case where the timeout/Ctrl+C path wasn't exercised but the OS killed
        # the child via console close.
        final_stats = solver_info.get('solver_stats_final') or {}
        if final_stats.get('interrupted'):
            partial = True

        status = "partial" if partial else "completed"
        n_pts = sol_arr.shape[0]
        print(f"C++ solver {status} -- {n_pts} time points")

        t = sol_arr[:, 0]
        y = sol_arr[:, 1:].T   # (N_tot, n_pts)

        # Split off the appended ⟨100⟩ SIA block so downstream post-processing
        # sees the original [SIA | VAC | He | conservation] layout unchanged.
        # Two views of the block are kept:
        #   y_sia100_raw  — the block AS SOLVED (discrete: per-size; bin_moment:
        #                   discrete-prefix + bin-moment vector). This goes to
        #                   calculate_derived_quantities, which forms the ⟨100⟩
        #                   mean size / density at the MOMENT level (subtracting
        #                   the C_floor IC per bin), exactly as the main ½⟨111⟩
        #                   population is handled.  Reconstructing to per-size
        #                   first over-floors near-empty bins at large n and
        #                   inflates the size-weighted mean during the transient.
        #   y_sia100_full — the per-size [I, n_pts] reconstruction, kept only for
        #                   the size-distribution plots (visualization.py).
        y_sia100     = None   # raw block → moment-level scalars
        y_sia100_full = None  # per-size reconstruction → plots/results
        if _n_sia100:
            y_sia100 = y[-_n_sia100:, :]
            y        = y[:-_n_sia100, :]
            if _is_bin:
                from .bin_moment_rates import reconstruct_distribution
                from .post_process import _floor_bin_moments
                _i_d   = int(getattr(re_obj, 'i_discrete', 0))
                _Pm    = int(getattr(re_obj, 'n_mom', 2))
                _sf    = getattr(re_obj, 'shape_function', 'linear')
                _bins  = re_obj.bins
                _Ibin  = len(_bins)
                _Ifull = int(sim.input_data.I)
                _npts  = y_sia100.shape[1]
                _Cf    = float(sim.input_data.reactions.get('C_floor', 1e-15))
                _c100  = np.zeros((_Ifull, _npts))
                # THE C_floor IC IS REMOVED AT THE MOMENT LEVEL, BEFORE the
                # closure reconstructs per-size values -- the same rule
                # post_process._floor_bin_moments enforces, and for the same
                # reason (its docstring: an empty bin that still carries content
                # is how mean_n_100 reached 2853 on a grid of I = 1000).
                #
                # Reconstructing from RAW moments and subtracting C_floor per
                # size afterwards is NOT equivalent: the linear closure spreads a
                # near-empty bin's floor residual across its sizes, so the clamp
                # leaves a positive size-weighted tail at large n.  The <100>
                # population sits near-empty for most of a run -- it exists only
                # by conversion -- so it is in that regime continuously, and the
                # tail is not a small correction.  Measured on the T3 campaign
                # (2026-08-12): the resulting N_100_vis_1 exceeded the
                # moment-level N_loops_100 in 95.3 % of 1173 rows, median 626x,
                # max 65540x, and produced the d_100 = 22.87 nm ceiling pile-up
                # and f_100 = 0.9999 character inversion.
                #
                # This view was documented as plot-only, but run_ensemble.observe
                # reads it for N_100_vis_* / f_100_tem_*, i.e. for the numbers the
                # calibration is scored on.  Flooring here fixes every consumer at
                # once rather than in each caller.
                for _j in range(_npts):
                    _col = np.maximum(y_sia100[:, _j], 0.0)
                    _c = np.zeros(_Ifull)
                    _c[:_i_d] = np.maximum(_col[:_i_d] - _Cf, 0.0)
                    if _Ibin > 0:
                        _mom = _col[_i_d:_i_d + _Pm * _Ibin]
                        _mu0 = _mom[0::_Pm][:_Ibin].astype(float).copy()
                        _mu1 = (_mom[1::_Pm][:_Ibin].astype(float).copy()
                                if _Pm >= 2 else None)
                        _mu2 = _mom[2::_Pm][:_Ibin] if _Pm >= 3 else None
                        for _kb, (_nlo, _nhi) in enumerate(_bins):
                            _m0, _m1 = _floor_bin_moments(
                                _mu0[_kb], _mu1[_kb] if _Pm >= 2 else 0.0,
                                _nlo, _nhi, _Cf, _Pm)
                            _mu0[_kb] = _m0
                            if _Pm >= 2:
                                _mu1[_kb] = _m1
                        _rec = reconstruct_distribution(_sf, _mu0, _mu1, _mu2,
                                                        _bins, _Ifull)
                        _c[_i_d:] = _rec[_i_d:]
                    _c100[:, _j] = _c
                y_sia100_full = _c100
            else:
                # Discrete mode: the block is already per-size, so the raw and
                # reconstructed views coincide.
                y_sia100_full = y_sia100

        results = calculate_derived_quantities(t, y, sim.input_data, re_obj,
                                               y_sia100=y_sia100)
        results['y'] = y   # raw ODE state [N_eq, n_pts] in atom fraction
        if y_sia100_full is not None:
            results['y_sia100'] = y_sia100_full   # ⟨100⟩ per-size [I, n_pts]
        if y_sia100 is not None:
            # The ⟨100⟩ block AS SOLVED — the only view that can serve as an
            # initial condition for a continuation segment.  y_sia100_full is a
            # closure reconstruction with the C_floor IC already removed, so
            # round-tripping it through the solver is lossy; and results['y']
            # has the block split off entirely, so a caller that resumes from
            # results['y'][:, -1] alone hands the solver a SHORT y0 and every
            # ⟨100⟩ component silently falls back to the parameters.h default
            # of 1e-100 (i.e. the population is erased once per segment).
            results['y_sia100_raw'] = y_sia100    # [n_sia100, n_pts] as solved

        # ── Window-bounds sidecar (active_window mode tracks expansion) ───
        # The C++ solver writes <bin_path>.window.csv with one row per output
        # point: t,x_hi_i,x_hi_v.  We trust this over any heuristic derived
        # from y, because indices outside the window stay at C_floor (the
        # initial value), making y-based detection unreliable.
        win_path = bin_path + '.window.csv'
        try:
            if os.path.exists(win_path):
                wdata = np.genfromtxt(win_path, delimiter=',', skip_header=1)
                if wdata.ndim == 1:
                    wdata = wdata.reshape(1, -1)
                n_w = min(wdata.shape[0], len(t))
                # x_hi_i/x_hi_v are 0-indexed inclusive upper bounds; the
                # active window covers sizes 1..(x_hi+1).
                results['n_active_sia'] = (wdata[:n_w, 1].astype(int) + 1)
                results['n_active_vac'] = (wdata[:n_w, 2].astype(int) + 1)
                results['n_active'] = (results['n_active_sia']
                                       + results['n_active_vac'])
        except Exception as exc:
            print(f"  window sidecar parse failed: {type(exc).__name__}: {exc}")
        finally:
            try:
                os.unlink(win_path)
            except OSError:
                pass

        # Stash immediately so a Ctrl+C during the metadata stamping below
        # still leaves something for the orchestrator/notebook to save.
        sim._partial_results = results
    except KeyboardInterrupt:
        print("\n*** Ctrl+C during result post-processing -- "
              "returning whatever was rescued. ***")
        return sim._partial_results

    try:
        sm     = sim.input_data.solver_mode
        po     = sim.input_data.physics_option
        linsol = str(solver_config.get('solver_method', {}).get('linsol', 'dense')).upper()
        msg = f'C++ CVODE BDF {sm}/{po} / {linsol}'
        if partial:
            msg += ' (partial -- interrupted)'
        results['metadata'] = {
            'solver_stats': {
                'success':       True,
                'message':       msg,
                'n_time_points': n_pts,
                'partial':       bool(partial),
            },
            'omp_threads_used':   solver_info.get('omp_threads_used'),
            'solver_stats_final': final_stats or None,
        }
        print("Results processing complete.")
    except KeyboardInterrupt:
        print("\n*** Ctrl+C during metadata stamping -- "
              "returning rescued results without metadata. ***")
    return results
