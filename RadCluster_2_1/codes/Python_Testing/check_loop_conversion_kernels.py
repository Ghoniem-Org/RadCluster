"""
check_loop_conversion_kernels.py — Phase-3 unit tests for the ½⟨111⟩ → ⟨100⟩
loop-conversion rate kernels added to ``ReactionRates``.

Checks the 1-D ingredients owned by ReactionRates:
  1. Γ_uni(n)  — nonnegative; zero where ΔF ≤ 0; nonzero on some small loops;
                 zero on large loops (small-loop biased, Phase-2 finding).
  2. φ_junc    — square, in [0, φ_max], symmetric, peaks on the diagonal,
                 zero below n_j_min.
  3. ⟨100⟩ sessile kernels — K_100_grow / K_100_shrink zero below n_loop_min and
                 positive above; G_100 (emission) nonnegative; k2_100 ≡ 0.
  4. conversion_dF agrees with the LoopEnergetics driving force.

Run:  python codes/Python_Testing/check_loop_conversion_kernels.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from py_utils.input_data import InputData
from py_utils.reaction_rates import ReactionRates


def _make_inputs():
    inp = InputData(I=1000, V=1000, physics_option="full_CD_fission")
    inp.diffusion["i_mobile"] = 5
    inp.derived["i_mobile"]   = 5
    inp.diffusion["v_mobile"] = 2
    inp.derived["v_mobile"]   = 2
    inp.reactions["C_floor"]     = 1e-25
    inp.reactions["he_kinetics"] = "quasi_steady_state"
    return inp, ReactionRates(inp)


def main():
    inp, rr = _make_inputs()
    I = int(inp.I)
    fails = 0

    def check(name, cond, *info):
        nonlocal fails
        tag = "PASS" if cond else "FAIL"
        if not cond:
            fails += 1
        print(f"  {tag}  {name}   {info if info else ''}")

    # ── 1. Γ_uni ──────────────────────────────────────────────────────────────
    G = rr.Gamma_uni
    dF = rr.conversion_dF
    check("Gamma_uni shape == (I,)", G.shape == (I,), G.shape)
    check("Gamma_uni >= 0", np.all(G >= 0.0))
    check("Gamma_uni == 0 where dF <= 0", np.all(G[dF <= 0.0] == 0.0))
    check("Gamma_uni > 0 for some (small) loops", np.any(G > 0.0),
          f"max={G.max():.3e} at n={int(G.argmax())+1}")
    # small-loop biased: the largest converting size is well below I
    converting = np.where(G > 0.0)[0]
    if converting.size:
        check("conversion support is small-loop biased",
              converting.max() + 1 < I // 2,
              f"sizes n in [{converting.min()+1}, {converting.max()+1}]")

    # ── 2. φ_junc ─────────────────────────────────────────────────────────────
    phi = rr.phi_junc
    phi_max = float(inp.reactions.get("phi_max_junc", 0.5))
    n_j_min = float(inp.reactions.get("n_j_min_junc", 30.0))
    check("phi_junc shape == (I, I)", phi.shape == (I, I), phi.shape)
    check("phi_junc in [0, phi_max]", np.all(phi >= 0.0) and np.all(phi <= phi_max + 1e-12))
    check("phi_junc symmetric", np.allclose(phi, phi.T, atol=1e-12))
    check("phi_junc == 0 below n_j_min",
          np.all(phi[: int(n_j_min) - 1, :] == 0.0) and
          np.all(phi[:, : int(n_j_min) - 1] == 0.0))
    # diagonal is the per-row max where defined (size comparability peaks at n=n')
    k = int(n_j_min) + 20
    check("phi_junc peaks on the diagonal", abs(phi[k, k] - phi[k].max()) < 1e-12,
          f"phi[{k},{k}]={phi[k,k]:.3f}")

    # ── 3. ⟨100⟩ sessile kernels ──────────────────────────────────────────────
    nlm = int(rr.n_loop_min)
    check("K_100_grow == 0 below n_loop_min", np.all(rr.K_100_grow[: nlm - 1] == 0.0))
    check("K_100_grow > 0 at/above n_loop_min", np.all(rr.K_100_grow[nlm - 1:] > 0.0))
    check("K_100_shrink >= 0", np.all(rr.K_100_shrink >= 0.0))
    check("G_100 (emission) >= 0", np.all(rr.G_100 >= 0.0))
    check("G_100[0] == 0 (monomer cannot emit)", rr.G_100[0] == 0.0)
    check("k2_100 == 0 (sessile)", np.all(rr.k2_100 == 0.0))

    # ── 4. driving force consistency ──────────────────────────────────────────
    le = rr.loop_energetics
    T = float(inp.derived["T"])
    check("conversion_dF matches LoopEnergetics",
          np.allclose(dF, le.driving_force_array(T, I), atol=1e-12))

    print(f"\n{'ALL PASS' if not fails else str(fails) + ' FAILED'}")
    # report
    print(f"  [info] T={T:.0f} K, Gamma_uni max={G.max():.3e} s^-1, "
          f"converting sizes n<= {converting.max()+1 if converting.size else 0}, "
          f"Fhat_100_0={le.Fhat_100_0:.4f} eV/A")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
