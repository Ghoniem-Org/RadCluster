"""Is d_111's residual a CLOSURE error, or a PATH difference?

CLAUDE.md S1: with I_bin = 0 and i_discrete = I every equation is discrete, so
bin_moment reduces to full_CD.  Run that against the real discrete arm on an
identical grid at the dose where the problem shows:

    agree    -> the two paths solve the same equations; the -5% on d_111 is
                genuine closure error and must be explained by the closure
    disagree -> a path difference, which would ALSO explain why the error does
                not refine away with bin count

check_bm_vs_discrete.py already asks this, but at 1e-2 dpa.  The M2R deviation
is measured at 20 dpa, three orders of magnitude further out, so a difference
that only opens up as the distribution broadens is invisible to it.
"""
import os, sys, time
sys.path.insert(0, "/Users/ghoni/Documents/GitHub/RadCluster")
sys.path.insert(0, "/Users/ghoni/Documents/GitHub/RadCluster/RadCluster_2_1/digital_twin")
os.environ["OMP_NUM_THREADS"] = "6"
import numpy as np
from RadCluster_2_1.py_utils.simulation import RadClusterSimulation
import run_ensemble as re_mod

I = 1000
DOSE = float(os.environ.get("PDOSE", "20"))


def run(mode):
    sim = RadClusterSimulation(I=I, V=I, solver_mode="full_system",
                               equations=mode, cascade="fission",
                               he_kinetics="quasi_steady_state",
                               i_mobile=5, v_mobile=5)
    cfg = {"equations": mode, "i_discrete": I, "v_discrete": I,
           "I_bin": 0, "V_bin": 0, "shape_function": "linear"}
    re_mod.apply_bin_config(sim, cfg)
    sim.input_data.reactions["LOOP_NETWORK_LOSS"] = 1
    sim.input_data._calculate_derived()
    sim.rebuild_rates()
    G = float(sim.input_data.reactions["G"])
    scfg = {"t_span": (1e-6, DOSE / G), "n_points": 45, "log_time": True,
            "required_times": (20.0 / G,), "rtol": 1e-5, "atol": 1e-20,
            "timeout_s": 7200,
            "solver_method": {"linsol": "gmres", "preconditioner": "woodbury",
                              "concentration_threshold": 1e-22},
            "loop_conversion": 1}
    t0 = time.time()
    res = sim.run_adaptive(solver_config=scfg, save_output=False,
                           timeout_s=7200, max_doublings=0)
    w = time.time() - t0
    cfg2 = dict(cfg); cfg2.update({"I": I, "V": I, "dose": DOSE,
                                   "C_floor": float(sim.input_data.reactions.get("C_floor", 1e-15))})
    o = re_mod.observe(res, sim, cfg2, 1.0)
    return w, o, res


KEYS = ["N_loops_111", "d_111_nm", "N_loops_100", "d_100_nm",
        "N_voids", "d_cavity_nm", "mean_n_111", "delta_FP"]
out = {}
for mode in ("discrete", "bin_moment"):
    w, o, res = run(mode)
    out[mode] = o
    print(f"{mode:11s} wall={w:6.1f}s reached={o['dose_reached']:.4g} dpa")

print(f"\n{'observable':14s} {'discrete':>14s} {'bin_moment':>14s} {'rel diff':>10s}")
for k in KEYS:
    a, b = out["discrete"].get(k), out["bin_moment"].get(k)
    if a is None or b is None:
        continue
    rel = (b - a) / a * 100.0 if a else float("nan")
    print(f"{k:14s} {a:14.6g} {b:14.6g} {rel:+9.3f}%")
