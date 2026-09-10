"""Does the d_111 deficit scale with i_discrete rather than with I_bin?

The M2R ladder held i_discrete = 50 and varied I_bin 2 -> 33.  The deficit was
flat at ~-5 to -7% across all of it, and B8 agrees with B32 far better than
either agrees with the exact arm.  An error that does not respond to the knob
being turned is controlled by a knob that was NOT turned -- and the one held
fixed is the discrete/bin interface at i_discrete.

The I_bin = 0 control passed to 0.000% precisely because i_discrete = I there,
i.e. no interface at all.  So: hold I_bin and sweep i_discrete.  If the
deficit collapses as i_discrete grows, the fault is at the interface (or in the
binned region's treatment of sizes just above it), not in bin resolution.
"""
import os, sys, time
sys.path.insert(0, "/Users/ghoni/Documents/GitHub/RadCluster")
sys.path.insert(0, "/Users/ghoni/Documents/GitHub/RadCluster/RadCluster_2_1/digital_twin")
os.environ["OMP_NUM_THREADS"] = "4"
from RadCluster_2_1.py_utils.simulation import RadClusterSimulation
import run_ensemble as re_mod

I, DOSE = 1000, float(os.environ.get("PDOSE", "0.05"))
i_d, nb = int(sys.argv[1]), int(sys.argv[2])
eq = "discrete" if nb == 0 and i_d >= I else "bin_moment"

sim = RadClusterSimulation(I=I, V=I, solver_mode="full_system", equations=eq,
                           cascade="fission", he_kinetics="quasi_steady_state",
                           i_mobile=5, v_mobile=5)
cfg = {"equations": eq, "i_discrete": i_d, "v_discrete": i_d,
       "I_bin": nb, "V_bin": nb, "shape_function": "linear"}
re_mod.apply_bin_config(sim, cfg)
sim.input_data.reactions["LOOP_NETWORK_LOSS"] = 1
sim.input_data.reactions["LOOP_COAL"] = 1
sim.input_data._calculate_derived(); sim.rebuild_rates()
G = float(sim.input_data.reactions["G"])
scfg = {"t_span": (1e-6, DOSE / G), "n_points": 45, "log_time": True,
        "required_times": (DOSE / G,), "rtol": 1e-5, "atol": 1e-20,
        "timeout_s": 1800,
        "solver_method": {"linsol": "gmres", "preconditioner": "woodbury",
                          "concentration_threshold": 1e-22},
        "loop_conversion": 1}
t0 = time.time()
res = sim.run_adaptive(solver_config=scfg, save_output=False, timeout_s=1800,
                       max_doublings=0)
w = time.time() - t0
cfg2 = dict(cfg); cfg2.update({"I": I, "V": I, "dose": DOSE,
                               "C_floor": float(sim.input_data.reactions.get("C_floor", 1e-15))})
o = re_mod.observe(res, sim, cfg2, 1.0)
neq = int(getattr(sim.rate_equations, "N_eq", -1))
print(f"IDRESULT i_d={i_d:5d} I_bin={nb:3d} N_eq={neq:5d} wall={w:6.1f}s "
      f"N_111={o['N_loops_111']:.4e} d_111={o['d_111_nm']:.4f} "
      f"mean_n_111={o['mean_n_111']:.2f} d_100={o['d_100_nm']:.4f} "
      f"d_cav={o['d_cavity_nm']:.4f}")
