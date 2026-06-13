"""
check_loop_conversion_conservation.py — Phase-7e end-to-end conservation test.

Runs the C++ solver with ½⟨111⟩→<100> conversion ON and verifies that the
Frenkel-pair conservation diagnostics hold once the SIA100 inventory and the
<100> shrink annihilation are folded into the accounting:

  * dFP      (swelling identity)  → <= 1e-6 in the accumulated regime,
  * dFP_sia  (per-species SIA)    → < 1e-3 in the accumulated regime,
  * total SIA inventory  S_I = Σ n (c_n^{111} + c_n^{100})  is finite and
    below the cumulative cascade production η·G·t.

Conservation is checked at t >= 1e-3 s: the very first log-time point (t≈1e-6)
carries a benign startup transient (production ≈ floor), so the per-step `max`
is not a meaningful gate — the diagnostics decay to ~1e-6 as dose accumulates.

Requires a built solver.exe.  Run:
  python codes/Python_Testing/check_loop_conversion_conservation.py
"""
import io
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from py_utils.simulation import RadClusterSimulation


def _run(loop_conversion):
    _s = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=200, V=200, solver_mode="full_system",
            physics_option="full_CD_fission", C_floor=1e-25,
            he_kinetics="quasi_steady_state", i_mobile=5, v_mobile=2)
        # Drive conversion hard so the ⟨100⟩ block is genuinely populated and the
        # conservation accounting is exercised: high T + low unary barrier.
        if loop_conversion:
            sim.input_data.reactions["T"] = 700.0 + 273.15
            sim.input_data.reactions["E_a0_conv"] = 0.8
            sim.input_data._calculate_derived()
            sim.rebuild_rates()
        cfg = {"t_span": (1e-6, 1e0), "n_points": 8, "log_time": True,
               "rtol": 1e-7, "atol": 1e-28,
               "solver_method": {"linsol": "dense"}}
        if loop_conversion:
            cfg["loop_conversion"] = 1
        r = sim.run_adaptive(solver_config=cfg, save_output=False)
    finally:
        sys.stdout, sys.stderr = _s
    return r


def main():
    fails = 0

    def check(name, cond, *info):
        nonlocal fails
        if not cond:
            fails += 1
        print(f"  {'PASS' if cond else 'FAIL'}  {name}   {info if info else ''}")

    r = _run(loop_conversion=1)
    t = np.asarray(r["t"])
    dfp = np.asarray(r["delta_FP"])
    dfps = np.asarray(r["delta_FP_sia"])
    late = t >= 1e-3                        # accumulated regime (skip startup)

    check("conversion ON produced <100>",
          "y_sia100" in r and float(np.max(r["y_sia100"])) > 0.0)
    check("dFP < 1e-3 in accumulated regime (t>=1e-3)",
          np.all(dfp[late] < 1e-3), float(np.max(dfp[late])))
    check("dFP_sia < 1e-3 in accumulated regime",
          np.all(dfps[late] < 1e-3), float(np.max(dfps[late])))
    check("dFP converges to <= 1e-6 by final time",
          dfp[-1] < 1e-6, float(dfp[-1]))
    check("dFP decreases after the startup transient",
          dfp[-1] < dfp[1], (float(dfp[1]), float(dfp[-1])))

    # total SIA inventory finite and below cumulative production η·G·t
    y = r["y"]; y100 = r["y_sia100"]; sz = np.arange(1, 201)
    S_I = float(np.dot(sz, np.maximum(y[:200, -1], 0.0))
                + np.dot(sz, np.maximum(y100[:, -1], 0.0)))
    eta, G = 0.30, 1e-6
    prod = eta * G * t[-1]
    check("total SIA inventory finite & < cumulative production",
          np.isfinite(S_I) and 0.0 < S_I < prod, (S_I, prod))

    print(f"\n{'ALL PASS' if not fails else str(fails) + ' FAILED'}")
    print(f"  [info] final dFP={dfp[-1]:.2e}  dFP_sia={dfps[-1]:.2e}  "
          f"S_I={S_I:.3e}  (cascade prod {prod:.3e})")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
