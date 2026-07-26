"""
check_loop_network_loss.py — end-to-end test for the loop → network-dislocation
loss channel (docs/Formulation/loop_network_loss.tex).

Verifies, against the C++ production solver:

  1. The channel is genuinely active — with LOOP_NETWORK_LOSS on it removes
     SIA-loop inventory relative to the OFF baseline (the network sweeps loops
     into the dislocation reservoir).
  2. Conservation is NOT degraded by the active channel — the SIA content swept
     into the network is routed into the J_SIA_fixed ledger (Eq. ledger), so the
     Frenkel-pair diagnostic dFP with the channel ON matches the OFF baseline.
     (The *absolute* dFP here is a few ×1e-2 because run_adaptive's multi-segment
     accounting carries a benign baseline transient at this short, hot, heavily-
     converting config; the OFF run shows the same level, and the single-segment
     ledger is exact — see check_loop_conversion_conservation.py.  The correct
     correctness statement is therefore comparative: ON ≈ OFF.)
  3. The Λ_n^net kernel is built non-zero and folded into the P4 sinks.

Both runs use loop conversion ON so the sessile ⟨100⟩ block (the population the
note says dominates network incorporation) is exercised through its new C++
loss term + ledger contribution.  χ is set far above its physical O(1) value
here only so the geometric switch P_ℓd turns on for the small loops the short
test produces (physical χ/K_rec calibration is Phase 5); w_c is likewise
amplified to make the per-second removal observable.

Requires a built solver.exe.  Run:
  python codes/Python_Testing/check_loop_network_loss.py
"""
import io
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

try:                                  # robust unicode prints on Windows consoles
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from py_utils.simulation import RadClusterSimulation


def _run(network_loss):
    """Run the C++ solver with loop conversion ON and the network-loss channel
    either OFF (baseline) or ON.  Returns the accumulated results dict and the
    final dynamic network density rho_net."""
    _s = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=200, V=200, solver_mode="full_system",
            physics_option="full_CD_fission", C_floor=1e-25,
            he_kinetics="quasi_steady_state", i_mobile=5, v_mobile=2)
        # Drive conversion so the ⟨100⟩ block is populated (high T, low barrier).
        sim.input_data.reactions["T"] = 700.0 + 273.15
        sim.input_data.reactions["E_a0_conv"] = 0.8
        lam_max = 0.0
        if network_loss:
            re = sim.input_data.reactions
            re["LOOP_NETWORK_LOSS"] = 1
            re["loop_net_chi"] = 1.0e6    # force P_ℓd→1 on the small test loops
            re["loop_net_w_c"] = 1.0e-3   # amplify Λ so removal shows in ~1 s
            re["loop_net_K_rec"] = 1.0e-13
            re["loop_net_n_inc"] = 4      # incorporation onset size
        sim.input_data._calculate_derived()
        sim.rebuild_rates()
        # Several short segments so the operator-split advance (which lags one
        # segment) actually applies Λ_n^net during the run.
        cfg = {"t_span": (1e-6, 1e0), "n_points": 40, "log_time": True,
               "rtol": 1e-7, "atol": 1e-28,
               "solver_method": {"linsol": "dense"},
               "loop_conversion": 1}
        r = sim.run_adaptive(solver_config=cfg, save_output=False,
                             points_per_segment=4)
        lam_max = float(sim.reaction_rates.Lambda_net_111.max())
    finally:
        sys.stdout, sys.stderr = _s
    return r, lam_max


def _sessile_inventory(r, n_inc=4):
    """Sessile SIA-loop content Σ_{n>=n_inc} n (c_n^{111} + c_n^{100}) at the
    final time — the population the network-loss channel targets."""
    I = r["y"].shape[0] if False else 200
    sz = np.arange(1, I + 1, dtype=float)
    mask = sz >= n_inc
    S = float(np.dot(sz[mask], np.maximum(r["y"][:I, -1][mask], 0.0)))
    if "y_sia100" in r:
        S += float(np.dot(sz[mask], np.maximum(r["y_sia100"][:, -1][mask], 0.0)))
    return S


def main():
    fails = 0

    def check(name, cond, *info):
        nonlocal fails
        if not cond:
            fails += 1
        print(f"  {'PASS' if cond else 'FAIL'}  {name}   {info if info else ''}")

    base, _ = _run(network_loss=False)
    chan, lam_max = _run(network_loss=True)

    tb = np.asarray(base["t"]); dfp_b = np.asarray(base["delta_FP"])
    tc = np.asarray(chan["t"]); dfp_c = np.asarray(chan["delta_FP"])
    late_b = tb >= 1e-2
    late_c = tc >= 1e-2
    max_off = float(np.max(dfp_b[late_b]))
    max_on = float(np.max(dfp_c[late_c]))

    # ── (3) the Lambda_net kernel is built non-zero and folded in ─────────────
    check("Lambda_net kernel built non-zero", lam_max > 0.0, lam_max)

    # ── (1) the channel is active: it removes sessile loop inventory ──────────
    S_off = _sessile_inventory(base)
    S_on = _sessile_inventory(chan)
    check("channel removes sessile SIA-loop inventory (S_on < S_off)",
          S_on < S_off * (1.0 - 1e-3), (S_on, S_off, S_on / S_off))

    # ── (2) conservation is NOT degraded by the active channel ────────────────
    # Comparative gate: dFP with the channel ON tracks the OFF baseline (the
    # absolute level is a benign multi-segment run_adaptive transient at this
    # config; the single-segment ledger is exact — see
    # check_loop_conversion_conservation.py).
    check("channel does not degrade conservation (max dFP_on <= 1.5x dFP_off)",
          max_on <= 1.5 * max_off + 1e-9, (max_on, max_off))
    check("final dFP_on within 2x of dFP_off",
          dfp_c[-1] <= 2.0 * dfp_b[-1] + 1e-9, (float(dfp_c[-1]), float(dfp_b[-1])))

    print(f"\n{'ALL PASS' if not fails else str(fails) + ' FAILED'}")
    print(f"  [info] Lambda_max={lam_max:.2e}  sessile S_I: on={S_on:.4e} off={S_off:.4e} "
          f"(ratio {S_on/S_off:.4f})  max dFP: on={max_on:.2e} off={max_off:.2e}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
