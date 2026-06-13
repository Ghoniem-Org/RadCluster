"""
check_loop_conversion_params.py — Phase-5 test: the ½⟨111⟩→⟨100⟩ conversion
parameters are configurable through the Excel input path.

Verifies that the 9 conversion/junction parameters added to the `Reactions`
sheet:
  1. load into ``InputData.reactions`` with the documented defaults;
  2. are consumed by ``ReactionRates`` (LoopEnergetics T*/n_ref, n_loop_min);
  3. actually drive the physics — overriding ``T_star_conv_C`` shifts the
     calibrated ⟨100⟩ prelog magnitude F̂_100^0 and the conversion support.

Run:  python codes/Python_Testing/check_loop_conversion_params.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from py_utils.input_data import InputData
from py_utils.reaction_rates import ReactionRates

DEFAULTS = {
    "E_a0_conv": 1.8, "gamma_a_conv": 0.03, "nu0_conv": 1.0e13,
    "T_star_conv_C": 450.0, "n_ref_conv": 50, "phi_max_junc": 0.5,
    "sigma_s_junc": 0.35, "n_j_min_junc": 30, "n_loop_min": 4,
}

_fails = 0


def check(name, cond, *info):
    global _fails
    if not cond:
        _fails += 1
    print(f"  {'PASS' if cond else 'FAIL'}  {name}   {info if info else ''}")


def main():
    inp = InputData(I=200, V=200, physics_option="full_CD_fission")

    # 1. loaded with documented defaults
    for k, v in DEFAULTS.items():
        got = inp.reactions.get(k)
        check(f"{k} loaded", got is not None and abs(float(got) - float(v)) < 1e-9,
              got)

    # 2. consumed by ReactionRates
    rr = ReactionRates(inp)
    check("LoopEnergetics T* consumed",
          abs(rr.loop_energetics.T_star_C - DEFAULTS["T_star_conv_C"]) < 1e-9,
          rr.loop_energetics.T_star_C)
    check("n_loop_min consumed", rr.n_loop_min == DEFAULTS["n_loop_min"],
          rr.n_loop_min)

    # 3. params drive the physics — override T* and confirm F̂_100^0 + support move
    f0 = rr.loop_energetics.Fhat_100_0
    supp0 = int((rr.Gamma_uni > 0).sum())
    inp.reactions["T_star_conv_C"] = 520.0
    rr2 = ReactionRates(inp)
    f1 = rr2.loop_energetics.Fhat_100_0
    supp1 = int((rr2.Gamma_uni > 0).sum())
    check("override T* shifts Fhat_100_0", abs(f1 - f0) > 1e-6, (f0, f1))
    check("override T* shifts conversion support", supp1 != supp0, (supp0, supp1))

    print(f"\n{'ALL PASS' if not _fails else str(_fails) + ' FAILED'}")
    return 1 if _fails else 0


if __name__ == "__main__":
    sys.exit(main())
