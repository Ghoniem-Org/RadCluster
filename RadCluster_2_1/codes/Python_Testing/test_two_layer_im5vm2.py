"""
Stage-4 validation — the two-layer RadCluster_2_1 on the im5vm2 case.

This is the capstone regression test for the graph-based two-layer
refactor.  It reproduces the RadCluster_1_0 reference case

    full_system / full_CD_fission, I = V = 1000, i_mobile = 5,
    v_mobile = 2, T = 673 K, G = 1e-6 dpa/s, t in [1e-6, 1e-2] s

and checks two things:

  (A) STRUCTURE — that a RadClusterSimulation now carries the live
      two-layer description: an abstract-core ReactionAdmissibilityGraph
      instantiated by the EUROFER-97 material layer (Layer 1 + Layer 2).

  (B) PHYSICS  — that running through the restructured (core/ + materials/)
      C++ solver still reproduces the locked post-Stage-2 baseline, i.e.
      the refactor is behaviour-preserving.

Reference values are the post-Stage-2 baseline (the G_He single-source
fix is included by design); agreement to 1e-5 relative — comfortably
above the ~5e-7 solver-tolerance noise — counts as a pass.

Run:  python codes/Python_Testing/test_two_layer_im5vm2.py
Exit code 0 = PASS, 1 = FAIL.
"""
import sys
from pathlib import Path

ROOT = Path(r"d:/GitHub/RadCluster/RadCluster_2_1")
sys.path.insert(0, str(ROOT))

from py_utils.core import ReactionAdmissibilityGraph, StateLayout, EdgeClass
from py_utils.simulation import RadClusterSimulation

# ── locked post-Stage-2 baseline (im5vm2 @ T=673) ────────────────────────────
BASELINE = {
    "C_SIA_tot": 1.49422818e20,
    "C_VAC_tot": 2.37571338e20,
    "C_He_tot":  1.58830031e13,
    "mean_n_i":  7.11591710,
    "mean_n_v":  3.09365893,
    "swelling":  2.80334179e-9,
}
RTOL = 1e-5

failures = []


def check(name, ok, detail=""):
    mark = "PASS" if ok else "FAIL"
    print(f"  [{mark}] {name}{(' -- ' + detail) if detail else ''}")
    if not ok:
        failures.append(name)


print("=" * 72)
print("Stage-4 validation: two-layer RadCluster_2_1 on the im5vm2 case")
print("=" * 72)

# ── build the simulation (constructs the two-layer RAG description) ───────────
sim = RadClusterSimulation(
    I=1000, V=1000, solver_mode="full_system",
    physics_option="full_CD_fission",
    C_floor=1e-25, he_kinetics="quasi_steady_state",
    i_mobile=5, v_mobile=2,
)
sim.input_data.reactions["T"] = 673.0
sim.input_data._calculate_derived()
sim.rebuild_rates()

# ── (A) STRUCTURE checks ─────────────────────────────────────────────────────
print("\n(A) Two-layer structure")
rag = getattr(sim, "material_rag", None)
layout = getattr(sim, "state_layout", None)

check("material RAG present", isinstance(rag, ReactionAdmissibilityGraph))
check("state layout present", isinstance(layout, StateLayout))

if isinstance(rag, ReactionAdmissibilityGraph):
    rag_ok = True
    try:
        rag.validate()
    except Exception as exc:                       # noqa: BLE001
        rag_ok = False
        check("RAG validates", False, repr(exc))
    if rag_ok:
        check("RAG validates", True)
    check("2 populations (1 SIA + 1 vacancy bulk)",
          len(rag.populations) == 2,
          f"got {len(rag.populations)}")
    active = {ec for ec in EdgeClass if rag.edges_of_class(ec)}
    expected = {EdgeClass.GROWTH, EdgeClass.SHRINKAGE, EdgeClass.DISSOCIATION,
                EdgeClass.RECOMBINATION, EdgeClass.ANNIHILATION,
                EdgeClass.COALESCENCE, EdgeClass.SOURCE, EdgeClass.SINK}
    check("all P1-P6 edge classes present", expected <= active,
          f"missing {expected - active}")
    check("P7/P8 edges declared (growth>=3 or solute_trapping)",
          len(rag.edges_of_class(EdgeClass.GROWTH)) >= 3 or
          bool(rag.edges_of_class(EdgeClass.SOLUTE_TRAPPING)),
          f"{len(rag.edges)} edges total")
    print(f"      {rag.summary()}")

# ── (B) PHYSICS check — run and compare to the locked baseline ───────────────
print("\n(B) Behaviour-preserving run")
cfg = sim._default_solver_config()
cfg["t_span"]   = (1e-6, 0.01)
cfg["n_points"] = 200
cfg["log_time"] = True
cfg["rtol"]     = 1e-6
cfg["atol"]     = 1e-20
cfg["solver_method"]["linsol"] = "gmres"

results = sim.run_adaptive(solver_config=cfg, save_output=False,
                           max_doublings=0, points_per_segment=200)

if results is None:
    check("solver run completed", False, "run returned None")
else:
    check("solver run completed", True)

    def last(key):
        v = results.get(key)
        try:
            return float(v[-1])
        except (TypeError, IndexError):
            return float(v)

    print(f"  {'quantity':12s} {'baseline':>16s} {'this run':>16s} {'rel.diff':>11s}")
    for key, ref in BASELINE.items():
        val = last(key)
        rel = abs(val - ref) / abs(ref)
        print(f"  {key:12s} {ref:16.6e} {val:16.6e} {rel:11.2e}")
        check(f"{key} within {RTOL:g}", rel < RTOL, f"rel.diff {rel:.2e}")

    # conservation diagnostics: finite and small (delta_He no longer NaN)
    dfp, dhe = last("delta_FP"), last("delta_He")
    import math
    check("delta_FP finite and < 1e-3",
          math.isfinite(dfp) and abs(dfp) < 1e-3, f"delta_FP={dfp:.2e}")
    check("delta_He finite and < 1e-3 (NaN bug fixed)",
          math.isfinite(dhe) and abs(dhe) < 1e-3, f"delta_He={dhe:.2e}")

# ── verdict ──────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
if failures:
    print(f"RESULT: FAIL  ({len(failures)} check(s) failed: {failures})")
    sys.exit(1)
print("RESULT: PASS  -- two-layer RadCluster_2_1 reproduces the im5vm2 "
      "baseline\n          and carries the abstract-core + EUROFER-97 RAG "
      "description.")
sys.exit(0)
