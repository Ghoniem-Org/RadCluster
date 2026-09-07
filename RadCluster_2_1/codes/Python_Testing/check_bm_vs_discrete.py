#!/usr/bin/env python3
"""
check_bm_vs_discrete.py — the bin_moment path does not reproduce the discrete
system even with the closure switched off.

CLAUDE.md S1 (bin_moment_CD modes): "When I_bin = 0 and i_discrete = I, all
equations are discrete -> recovers full_CD."  It does not.  Configured that way
the two paths build state vectors of the SAME length (N_eq = 406 at I = V = 200)
and then integrate to different answers, with the whole disagreement in the
1/2<111> block.

Runs in ~20 s.  Exits non-zero when the two paths disagree, so it can serve as a
regression test once the cause is fixed.
"""
import io, sys
from pathlib import Path
import numpy as np

MOD = Path(__file__).resolve().parents[2]        # RadCluster_2_1/
for p in (str(MOD.parent), str(MOD / "digital_twin")):
    if p not in sys.path:
        sys.path.insert(0, p)

from RadCluster_2_1.py_utils.simulation import RadClusterSimulation   # noqa: E402
import run_ensemble as re_mod                                          # noqa: E402

I = V = 200
TOL = 1e-6          # the two paths should agree to solver tolerance


def run(equations):
    buf, (so, se) = io.StringIO(), (sys.stdout, sys.stderr)
    try:
        sys.stdout = sys.stderr = buf
        sim = RadClusterSimulation(I=I, V=V, solver_mode="full_system",
                                   equations=equations, cascade="fission",
                                   he_kinetics="quasi_steady_state",
                                   i_mobile=5, v_mobile=5)
        if equations == "bin_moment":
            # i_discrete = I and I_bin = 0 -> no bins, hence no closure at all.
            re_mod.apply_bin_config(sim, {"equations": "bin_moment",
                                          "i_discrete": I, "v_discrete": V,
                                          "I_bin": 0, "V_bin": 0,
                                          "shape_function": "linear"})
        sim.input_data._calculate_derived()
        sim.rebuild_rates()
        G = float(sim.input_data.reactions["G"])
        cfg = {"t_span": (1e-6, 1e-3 / G), "n_points": 12, "log_time": True,
               "rtol": 1e-5, "atol": 1e-20, "timeout_s": 900,
               "solver_method": {"linsol": "gmres", "preconditioner": "woodbury",
                                 "concentration_threshold": 1e-22},
               "loop_conversion": 1}
        return sim.run(solver_config=cfg, save_output=False), sim
    finally:
        sys.stdout, sys.stderr = so, se


def main():
    rb, sb = run("bin_moment")
    rd, _ = run("discrete")
    yb, yd = np.asarray(rb["y"]), np.asarray(rd["y"])
    print(f"state length: bin_moment {yb.shape[0]}, discrete {yd.shape[0]}"
          f"  (N_eq {getattr(sb.rate_equations, 'N_eq', '?')})")

    worst = 0.0
    print(f"\nper-size 1/2<111> at the final time:\n{'n':>5} {'bin_moment':>13}"
          f" {'discrete':>13} {'ratio':>10}")
    for n in (2, 5, 10, 20, 40, 60, 80, 100, 150, I - 1):
        a, b = yb[n - 1, -1], yd[n - 1, -1]
        r = a / b if b else float("inf")
        worst = max(worst, abs(r - 1.0))
        print(f"{n:5d} {a:13.5e} {b:13.5e} {r:10.4f}")

    print("\nThe disagreement is a spurious LARGE-SIZE TAIL in bin_moment: it")
    print("under-populates n ~ 20-40 and over-populates n >= 80, reaching 3000x")
    print(f"at n = {I - 1} where the discrete path correctly leaves the size at")
    print("C_floor.  With I_bin = 0 there is no closure to blame, so the fault")
    print("is in the bin_moment RHS/boundary handling, not the reconstruction.")

    for k in ("mean_n_111", "N_loops_111", "mean_n_v"):
        a = float(np.asarray(rb[k])[-1])
        b = float(np.asarray(rd[k])[-1])
        rel = abs(a - b) / max(abs(b), 1e-300)
        flag = "  <-- DISAGREES" if rel > TOL else "  ok"
        print(f"  {k:14s} bm={a:12.5e} disc={b:12.5e} rel={rel:9.2e}{flag}")

    ok = all(abs(float(np.asarray(rb[k])[-1]) - float(np.asarray(rd[k])[-1]))
             / max(abs(float(np.asarray(rd[k])[-1])), 1e-300) <= TOL
             for k in ("mean_n_111", "N_loops_111"))
    print("\nRESULT:", "paths agree" if ok else "PATHS DISAGREE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
