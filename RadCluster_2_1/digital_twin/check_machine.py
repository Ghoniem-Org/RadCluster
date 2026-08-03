#!/usr/bin/env python
"""Run this on EVERY machine before it joins the campaign.

    python check_machine.py            # compare against the committed reference
    python check_machine.py --write    # (re)generate the reference, once, on one machine

Why this exists
---------------
The C++ solver is NOT in git -- `build/` is ignored, so every machine compiles
its own from cpp_utils/.  Four machines with different compiler or SUNDIALS
versions produce numerically different results while every hash in the
provenance block still looks plausible, and merge_and_sobol would then mix them
into one Sobol estimate.  That is the same class of error as the grid
saturation this project already hit: a believable number with no validity.

So: pin one fixed configuration, run it, and compare against a reference
committed from a known-good machine.  It takes ~30 s and it is the difference
between "the four machines agree" and "the four machines were assumed to
agree".

The check deliberately uses a SMALL, FAST configuration.  It is verifying that
the toolchain agrees, not that the physics is right -- the run is expected to be
inadmissible (the grid is far too coarse), and that is fine.

Tolerance: 1e-6 relative -- the campaign's own integration rtol.  CVODE with
fixed tolerances is bit-deterministic for a FIXED BINARY, so an earlier 1e-9
bar was right as long as every machine ran the same build.  It is unreachable
across ARCHITECTURES: Windows/x86 (MSVC + MKL, /arch:AVX2) and macOS/arm64
(AppleClang + Accelerate, -march=native) reassociate and vectorise the same
reductions differently, which shows up as ~1e-7 drift on the CVODE-integrated
fields while the closed-form fields (Di_eff, Dv_eff, conv_*) still match to
0.0 exactly.  Holding the gate tighter than the integrator's own accuracy
rejects machines whose physics is identical, so the bar is set AT the
integration tolerance: anything above it is a real discrepancy, anything below
it is beneath the resolution the solver claims in the first place.
OMP_NUM_THREADS is pinned to 1 so reduction order cannot vary within a build.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO))
REF = HERE / "machine_reference.json"

# Fixed probe configuration -- do not change without regenerating the reference
# on a known-good machine and saying so in the commit message.
# solver_mode is part of the probe: it must match what the campaign actually
# runs, or the agreement check validates a code path nothing uses.
#
# V=480, not 120.  At V=120 this probe reported mean_n_v = 18.7 against a
# converged 153.3, N_voids 10.6x low, and delta_FP = 4.1e-3 -- above CLAUDE.md
# S8's "coding error" threshold.  The vacancy grid was truncating: mean cavity
# size is ~153 vacancies, so V=120 could not hold a typical cavity at all.
# delta_FP falls to 1.4e-7 at V=480 and is unchanged at V=960 (2026-08-02).
# I=300, not 150, for the same reason on the other axis.  At I=150 delta_FP is
# 4-5e-3 REGARDLESS of V, because mean_n_111 = 35.9 on a 150-point grid is an
# occupancy of 0.24 and the 1/2<111> population itself truncates.  Note the
# asymmetry with the <100> pile-up, where delta_FP stayed at 1e-8: <100> is
# sessile and one-way, so its growth guard halts without losing atoms, whereas
# truncating the mobile 1/2<111> population does lose them.
# I=300, V=480 gives delta_FP = 1.4e-7 in ~35 s.
PROBE = dict(I=300, V=480, i_mobile=10, v_mobile=5, T_K=623.15, dose=0.02,
             rho_d=1.0e14, E_a0_conv=2.2, dH2_abs_conv=0.30, psucc_abs_pref=2.0,
             gamma_a_conv=0.02, phi_max_junc=0.0, f_cl_i=0.25, f_cl_v=0.05,
             rtol=1e-6, solver_mode="active_window")

# Preconditioner follows solver_mode -- never selected independently (CLAUDE.md S9).
PRECOND = {"full_system": "Woodbury", "active_window": "Jacobi"}

# Quantities compared.  Mixed scales on purpose: a state variable, a derived
# rate, an integrated observable and a conservation residual, so a build
# difference cannot hide in a corner of the model.
FIELDS = ["Di_eff", "Dv_eff", "conv_psuccess", "conv_psuccess_abs",
          "N_loops_100", "N_loops_111", "mean_n_100", "mean_n_111",
          "C_SIA_tot", "C_VAC_tot", "swelling", "delta_FP"]
RTOL = 1e-6      # = PROBE["rtol"]; see module docstring -- the gate cannot be
                 # tighter than the integration tolerance it is checking


def sha256_file(p: Path) -> str:
    try:
        return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    except Exception:
        return "unavailable"


def probe() -> dict:
    from RadCluster_2_1.py_utils.simulation import RadClusterSimulation
    import numpy as np
    _s = sys.stdout, sys.stderr
    buf = io.StringIO()
    try:
        sys.stdout = sys.stderr = buf
        sim = RadClusterSimulation(I=PROBE["I"], V=PROBE["V"],
                                   solver_mode=PROBE["solver_mode"],
                                   equations="discrete", cascade="fission",
                                   C_floor=1e-25,
                                   he_kinetics="quasi_steady_state",
                                   i_mobile=PROBE["i_mobile"],
                                   v_mobile=PROBE["v_mobile"])
        d = sim.input_data
        for k in ("T", "rho_d", "E_a0_conv", "dH2_abs_conv", "psucc_abs_pref",
                  "gamma_a_conv", "phi_max_junc"):
            key = "T" if k == "T" else k
            d.reactions[key] = PROBE["T_K"] if k == "T" else PROBE[k]
        d.production_fission["f_cl_i"] = PROBE["f_cl_i"]
        d.production_fission["f_cl_v"] = PROBE["f_cl_v"]
        d._calculate_derived()
        sim.rebuild_rates()
        rr = sim.reaction_rates
        G = float(d.reactions["G"])
        res = sim.run(solver_config={
            "t_span": (1e-6, PROBE["dose"] / G), "n_points": 12,
            "log_time": True, "rtol": PROBE["rtol"], "atol": 1e-20,
            "timeout_s": 900,
            "solver_method": {"linsol": "gmres",
                              "preconditioner": PRECOND[PROBE["solver_mode"]],
                              "concentration_threshold": 1e-22},
            "loop_conversion": 1}, save_output=False)
    finally:
        sys.stdout, sys.stderr = _s

    out = {"Di_eff": float(d.derived["Di_eff"]),
           "Dv_eff": float(d.derived["Dv_eff"]),
           "conv_psuccess": float(rr.conv_psuccess),
           "conv_psuccess_abs": float(rr.conv_psuccess_abs)}
    for k in FIELDS:
        if k in out:
            continue
        v = res.get(k)
        out[k] = float(np.asarray(v, float).ravel()[-1]) if v is not None else None
    return out


def env() -> dict:
    solver = REPO / "RadCluster_2_1" / "build" / "Release" / "solver.exe"
    if not solver.exists():
        solver = REPO / "RadCluster_2_1" / "build" / "solver"
    try:
        gs = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(REPO),
                                     stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        gs = "unknown"
    return {"machine_id": platform.node(), "platform": platform.platform(),
            "python": platform.python_version(), "git_sha": gs,
            "solver_sha256": sha256_file(solver),
            "solver_exists": solver.exists(),
            "workbook_sha256": sha256_file(
                REPO / "RadCluster_2_1" / "input" / "input_parameters.xlsx")}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true",
                    help="generate the reference (run once, on a known-good machine)")
    a = ap.parse_args(argv)

    e = env()
    print(f"machine   {e['machine_id']}  ({e['platform']})")
    print(f"python    {e['python']}")
    print(f"git       {e['git_sha'][:12]}")
    print(f"solver    {e['solver_sha256'][:16]}  exists={e['solver_exists']}")
    print(f"workbook  {e['workbook_sha256'][:16]}")
    if not e["solver_exists"]:
        print("\n  solver binary NOT FOUND. Build it first:")
        print("      cd RadCluster_2_1 && cmake -S cpp_utils -B build && "
              "cmake --build build --config Release")
        return 2

    print("\nrunning probe (I=150, 0.02 dpa, ~30 s) ...")
    got = probe()

    if a.write:
        REF.write_text(json.dumps({"probe": PROBE, "fields": got,
                                   "generated_on": e}, indent=2),
                       encoding="utf-8")
        print(f"\n  wrote {REF.name} from this machine. Commit it, and have every")
        print("  other machine run `python check_machine.py` before joining.")
        for k in FIELDS:
            print(f"    {k:20s} {got[k]!r}")
        return 0

    if not REF.exists():
        print(f"\n  no {REF.name} committed yet. Generate it on the machine you")
        print("  trust with:  python check_machine.py --write")
        return 2

    ref = json.loads(REF.read_text(encoding="utf-8"))
    if ref["probe"] != PROBE:
        print("\n  *** the probe configuration has changed since the reference "
              "was written; regenerate it.")
        return 2

    print(f"\n  reference generated on {ref['generated_on']['machine_id']} "
          f"(git {ref['generated_on']['git_sha'][:12]}, "
          f"solver {ref['generated_on']['solver_sha256'][:16]})")
    if ref["generated_on"]["git_sha"] != e["git_sha"]:
        print("  note: git SHA differs from the reference machine - pull first "
              "if that is not intentional.")

    bad = []
    print(f"\n  {'field':20s} {'this machine':>16s} {'reference':>16s} {'rel diff':>11s}")
    for k in FIELDS:
        g, r = got.get(k), ref["fields"].get(k)
        if g is None or r is None:
            print(f"  {k:20s} {'None' if g is None else g:>16} "
                  f"{'None' if r is None else r:>16}   (skipped)")
            continue
        den = max(abs(r), 1e-300)
        rel = abs(g - r) / den
        flag = "" if rel <= RTOL else "   <-- MISMATCH"
        if rel > RTOL:
            bad.append((k, g, r, rel))
        print(f"  {k:20s} {g:16.9e} {r:16.9e} {rel:11.2e}{flag}")

    if bad:
        print(f"\n  *** {len(bad)} field(s) exceed rtol={RTOL:g}.")
        print("  This machine does NOT reproduce the reference. Its results")
        print("  cannot be pooled with the others in one Sobol estimate.")
        print("  Most likely a different SUNDIALS/compiler version - compare")
        print("  solver_sha256 above, rebuild from the same commit, re-check.")
        return 1
    print(f"\n  all {len(FIELDS)} fields agree within {RTOL:g}.")
    print("  This machine may join the campaign.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
