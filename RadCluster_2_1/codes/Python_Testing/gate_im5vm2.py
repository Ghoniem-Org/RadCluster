"""
Stage-2 gate / Stage-4 baseline driver.

Reproduces the RadCluster_1_0 reference run
  20260504_201741_full_system_full_CD_fission_I1000V1000_im5vm2
under RadCluster_2_1, and prints the macroscopic observables so the
post-Stage-2 numbers can be compared against the v1.0 baseline.

v1.0 baseline (for reference):
  C_SIA_tot = 1.3205355691e20   C_VAC_tot = 1.9494579688e20
  C_He_tot  = 8.9886e11         mean_n_i  = 7.6703   mean_n_v = 3.2384
  swelling  = 2.3004e-9         delta_FP  = 9.9995e-5   delta_He = nan
"""
import sys
from pathlib import Path

ROOT = Path(sys.argv[1] if len(sys.argv) > 1
            else r"d:/GitHub/RadCluster/RadCluster_2_1")
sys.path.insert(0, str(ROOT))
print(f"ROOT = {ROOT}")

from py_utils.simulation import RadClusterSimulation

sim = RadClusterSimulation(
    I=1000, V=1000, solver_mode="full_system",
    physics_option="full_CD_fission",
    C_floor=1e-25, he_kinetics="quasi_steady_state",
    i_mobile=5, v_mobile=2,
)

# The v1.0 baseline run was recorded at T=673 K; the workbook has since
# been edited to 573 K.  Pin T to the baseline value for a faithful compare.
sim.input_data.reactions["T"] = 673.0
sim.input_data._calculate_derived()
sim.rebuild_rates()
print(f"T pinned to {sim.input_data.derived['T']} K")

cfg = sim._default_solver_config()
cfg["t_span"]   = (1e-6, 0.01)
cfg["n_points"] = 200
cfg["log_time"] = True
cfg["rtol"]     = 1e-6
cfg["atol"]     = 1e-20
cfg["solver_method"]["linsol"] = "gmres"

results = sim.run_adaptive(
    solver_config=cfg, save_output=False,
    max_doublings=0, points_per_segment=200,
)

print("RESULTS_BEGIN")
if results is None:
    print("RUN FAILED")
else:
    def last(key):
        v = results.get(key)
        if v is None:
            return None
        try:
            return float(v[-1])
        except (TypeError, IndexError):
            try:
                return float(v)
            except (TypeError, ValueError):
                return v
    for k in ["C_SIA_tot", "C_VAC_tot", "C_He_tot", "mean_n_i",
              "mean_n_v", "swelling", "delta_FP", "delta_He"]:
        print(f"{k} = {last(k)!r}")
print("RESULTS_END")
