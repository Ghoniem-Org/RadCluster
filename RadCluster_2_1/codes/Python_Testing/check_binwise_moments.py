"""Bin-by-bin comparison, with NO reconstruction anywhere.

Takes the exact arm's per-size 1/2<111> spectrum, SUMS it onto the binned
run's own bin edges to form the exact mu0 and mu1 per bin, and compares those
against the moments the binned run actually carries.  Reconstruction is what
the closure approximates, so comparing reconstructed spectra would measure the
closure twice; comparing moments measures the dynamics alone.
"""
import glob, sys
import numpy as np
sys.path.insert(0, "/Users/ghoni/Documents/GitHub/RadCluster")
sys.path.insert(0, "/Users/ghoni/Documents/GitHub/RadCluster/RadCluster_2_1/digital_twin")
from RadCluster_2_1.py_utils import visualization as viz
from RadCluster_2_1.py_utils.simulation import RadClusterSimulation
import run_ensemble as _re

OUT = "/Users/ghoni/Documents/GitHub/RadCluster/RadCluster_2_1/output"
I = 1000
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 32


def load(pat):
    f = sorted(glob.glob(f"{OUT}/{pat}/plots/plot_data.pkl"))
    return viz.load_plot_data(f[-1]) if f else None


# Bin edges: rebuild the same layout the run used.
sim = RadClusterSimulation(I=I, V=I, solver_mode="full_system",
                           equations="bin_moment", cascade="fission",
                           he_kinetics="quasi_steady_state",
                           i_mobile=5, v_mobile=5)
_re.apply_bin_config(sim, {"equations": "bin_moment", "i_discrete": 50,
                           "v_discrete": 50, "I_bin": NB, "V_bin": NB,
                           "shape_function": "linear"})
sim.input_data._calculate_derived(); sim.rebuild_rates()
RE = sim.rate_equations
bins = list(RE.bins)
i_d, P = int(RE.i_discrete), int(RE.n_mom)
print(f"layout: i_discrete={i_d} I_bin={len(bins)} P={P} r={RE.r:.4f}")

ref = load("*_M1_D1000_*")
bm = load(f"*_M2R_B{NB}_*")
if ref is None or bm is None:
    print("missing plot_data"); sys.exit(1)

yr = np.asarray(ref["results"]["y"], float)          # [N_eq, n_t] discrete
yb = np.asarray(bm["results"]["y"], float)
C_floor = float(ref["input_data"].reactions.get("C_floor", 1e-15))
jr, jb = -1, -1
print(f"ref dose {float(np.asarray(ref['results']['dose'])[jr]):.4g}  "
      f"bm dose {float(np.asarray(bm['results']['dose'])[jb]):.4g}")

c_ref = np.maximum(yr[:I, jr], 0.0)                  # sizes 1..I
print(f"\n{'bin':>4} {'sizes':>13} {'mu0 exact':>12} {'mu0 bm':>12} {'d%':>8} "
      f"{'mu1 exact':>12} {'mu1 bm':>12} {'d%':>8}")
tot = {"e0": 0.0, "b0": 0.0, "e1": 0.0, "b1": 0.0}
# discrete head, sizes 2..i_d
he = float(np.sum(np.maximum(c_ref[1:i_d] - C_floor, 0.0)))
hb = float(np.sum(np.maximum(yb[1:i_d, jb] - C_floor, 0.0)))
ns = np.arange(2, i_d + 1)
he1 = float(np.dot(ns, np.maximum(c_ref[1:i_d] - C_floor, 0.0)))
hb1 = float(np.dot(ns, np.maximum(yb[1:i_d, jb] - C_floor, 0.0)))
print(f"{'head':>4} {f'2-{i_d}':>13} {he:12.5g} {hb:12.5g} "
      f"{(hb-he)/he*100:+7.2f}% {he1:12.5g} {hb1:12.5g} {(hb1-he1)/he1*100:+7.2f}%")
tot["e0"] += he; tot["b0"] += hb; tot["e1"] += he1; tot["b1"] += hb1

for k, (nlo, nhi) in enumerate(bins):
    # The last bin's nhi can exceed I -- build_bins walks floor(edge*r)+1 and
    # the final edge overshoots the domain.  Clip to what the exact arm
    # actually holds, or the size vector and the spectrum slice disagree.
    hi = min(nhi, I)
    seg = c_ref[nlo - 1:hi]                      # sizes nlo..hi inclusive
    sizes = np.arange(nlo, hi + 1)
    e0 = float(np.sum(seg)); e1 = float(np.dot(sizes, seg))
    off = i_d + P * k
    b0 = float(yb[off, jb]); b1 = float(yb[off + 1, jb]) if P >= 2 else 0.0
    d0 = (b0 - e0) / e0 * 100 if e0 else float("nan")
    d1 = (b1 - e1) / e1 * 100 if e1 else float("nan")
    if e0 > 1e-14 or abs(d0) > 1:
        print(f"{k:4d} {f'{nlo}-{nhi}':>13} {e0:12.5g} {b0:12.5g} {d0:+7.2f}% "
              f"{e1:12.5g} {b1:12.5g} {d1:+7.2f}%")
    tot["e0"] += e0; tot["b0"] += b0; tot["e1"] += e1; tot["b1"] += b1

print(f"\nTOTAL mu0 exact={tot['e0']:.6g} bm={tot['b0']:.6g} "
      f"({(tot['b0']-tot['e0'])/tot['e0']*100:+.2f}%)")
print(f"TOTAL mu1 exact={tot['e1']:.6g} bm={tot['b1']:.6g} "
      f"({(tot['b1']-tot['e1'])/tot['e1']*100:+.2f}%)")
print(f"mean n  exact={tot['e1']/tot['e0']:.4f}  bm={tot['b1']/tot['b0']:.4f} "
      f"({(tot['b1']/tot['b0'])/(tot['e1']/tot['e0'])*100-100:+.2f}%)")
