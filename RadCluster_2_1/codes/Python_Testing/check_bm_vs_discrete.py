#!/usr/bin/env python3
"""
check_bm_vs_discrete.py — the discrete path silently ignores LOOP_COAL.

CLAUDE.md S1 says that at I_bin = 0 and i_discrete = I the bin_moment system IS
the discrete system.  It is not, and the reason is NOT a bin_moment defect:

    P.loop_coal -- 1/2<111> loop-loop coarsening, rate_kernels.cpp:1944 -- is
    implemented ONLY in rhs_bin_moment.  rhs_full_CD contains zero references
    to it.  The workbook ships LOOP_COAL = 1, so the discrete path accepts the
    flag, reports it in provenance, and does nothing with it.

So the discrete arm is not the exact solution of the same equations: it is the
same equations MINUS the loop-coarsening channel.  The bin_moment answer is the
physically intended one.

rate_kernels.cpp:1930 predicts exactly what the discrete arm then shows, and
calls it "the missing edge": without coarsening "the 1/2<111> mean size pins AT
the cutoff" and "the density runs 20-60x over the EUROFER97 data".  Measured
here: discrete d_111 frozen at 1.759 / 1.763 / 1.760 nm across a 30x dose range,
and N_111 3.5x the bin_moment value.

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
# The two paths do identical arithmetic in DIFFERENT ORDERS, so CVODE's adaptive
# stepping cannot make them agree better than its own tolerance.  The runs use
# rtol = 1e-5; 1e-4 leaves an order of magnitude of headroom while still failing
# instantly on a missing channel (which showed as 0.43, not 1e-5).
TOL = 1e-4


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

    print("\nBoth paths must carry loop-loop coarsening (P.loop_coal) and pile a")
    print("product past the grid edge AT the edge.  When rhs_case2 lacked the")
    print("channel the per-size ratio ran 0.49 at n=40 and 906 at n=150; when it")
    print("had the channel but dropped the overflow instead of piling it, the")
    print("last size sat at 2.24x while everything below agreed to 0.04%.")

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
