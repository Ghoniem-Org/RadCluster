"""
check_loop_conversion_integration.py — Phase-4 integration test.

Two parts:

  1. Structural — the full EUROFER-97 RAG builds with the ½⟨111⟩/⟨100⟩ split:
     3 populations, the three conversion edges present (unary INTER_POPULATION,
     junction + absorption COALESCENCE), and a SIA100 state-layout block.

  2. Conservation — all SIA-content-conserving edges (GROWTH, self-coalescence,
     junction, unary transformation, ⟨100⟩ absorption) run *together* through the
     GraphWalker on a minimal two-SIA-population RAG and conserve the total SIA
     content  Σ_n n·c_n^(111) + Σ_m m·c_m^(100).  This catches edge-interaction
     bugs that the per-edge Phase-1 tests cannot.

(The full EUROFER walker is not exercised directly: its P8 SOLUTE_TRAPPING edge
is owned by the He mean-loading reduction, not the discrete reference walker.)

Run:  python codes/Python_Testing/check_loop_conversion_integration.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from py_utils.core.cluster_identifier import Polarity, Population
from py_utils.core.rag import ReactionAdmissibilityGraph, Edge
from py_utils.core.edge_classes import EdgeClass
from py_utils.core.state_layout import StateLayout
from py_utils.core.graph_walker import GraphWalker
from py_utils.input_data import InputData
from py_utils.reaction_rates import ReactionRates
from py_utils.materials import build_eurofer_rag

ATOL = 1e-12
_fails = 0


def check(name, cond, *info):
    global _fails
    tag = "PASS" if cond else "FAIL"
    if not cond:
        _fails += 1
    print(f"  {tag}  {name}   {info if info else ''}")


# ── Part 1: structural ────────────────────────────────────────────────────────
def test_structural():
    inp = InputData(I=1000, V=1000, physics_option="full_CD_fission")
    inp.diffusion["i_mobile"] = 5
    inp.derived["i_mobile"]   = 5
    inp.diffusion["v_mobile"] = 2
    inp.derived["v_mobile"]   = 2
    inp.reactions["C_floor"]     = 1e-25
    inp.reactions["he_kinetics"] = "quasi_steady_state"
    rr = ReactionRates(inp)
    rag, layout = build_eurofer_rag(inp, rr, equations="discrete", cascade="fission")

    pop_names = {(p.polarity.label, p.name) for p in rag.populations}
    check("3 populations incl. bulk-111 / bulk-100",
          ("sia", "bulk-111") in pop_names and ("sia", "bulk-100") in pop_names
          and len(rag.populations) == 3, sorted(pop_names))

    labels = {e.label: e for e in rag.edges}
    check("unary INTER_POPULATION edge present",
          labels.get("loop_111to100_unary") is not None and
          labels["loop_111to100_unary"].edge_class is EdgeClass.INTER_POPULATION)
    junc = labels.get("SIA111_junction")
    check("junction COALESCENCE -> bulk-100",
          junc is not None and junc.edge_class is EdgeClass.COALESCENCE and
          junc.product_population is not None and
          junc.product_population.name == "bulk-100")
    absb = labels.get("SIA100_absorb")
    check("absorption COALESCENCE on bulk-100",
          absb is not None and absb.edge_class is EdgeClass.COALESCENCE and
          absb.population.name == "bulk-100" and
          absb.partner_population.name == "bulk-111")
    check("SIA100 state block present", layout.has("SIA100"),
          layout.block("SIA100").length if layout.has("SIA100") else None)
    rag.validate()
    check("rag.validate() passes", True)


# ── Part 2: composed-conservation on a minimal RAG ────────────────────────────
def test_composed_conservation():
    N = 8
    rag = ReactionAdmissibilityGraph("mini")
    p111 = rag.add_population(Population("bulk-111", Polarity.SIA, n_min=1, mobile_max=N))
    p100 = rag.add_population(Population("bulk-100", Polarity.SIA, n_min=1, mobile_max=0))
    rag.set_monomer_population(p111)

    # Synthetic kernels (values irrelevant to conservation, only structure is).
    Kg = np.full(N, 0.05)                       # GROWTH (absorb I_1)
    Kfull = np.full((N, N), 0.02)               # ½⟨111⟩ collision kernel
    phi = np.zeros((N, N)); phi[2:, 2:] = 0.4   # junction branch (comparable, n≥3)
    Gamma = np.zeros(N); Gamma[2:5] = 0.03      # unary, small loops only
    Kabs = np.full((N, N), 0.015)               # ⟨100⟩ absorbs ½⟨111⟩
    rag.register_kernel("Kg", Kg)
    rag.register_kernel("K_self", (1.0 - phi) * Kfull)
    rag.register_kernel("K_junc", phi * Kfull)
    rag.register_kernel("Gamma", Gamma)
    rag.register_kernel("Kabs", Kabs)

    rag.add_edge(Edge(EdgeClass.GROWTH, "grow111", p111, kernel="Kg"))
    rag.add_edge(Edge(EdgeClass.COALESCENCE, "self111", p111, kernel="K_self"))
    rag.add_edge(Edge(EdgeClass.COALESCENCE, "junction", p111,
                      kernel="K_junc", product_population=p100))
    rag.add_edge(Edge(EdgeClass.INTER_POPULATION, "unary", p111,
                      kernel="Gamma", product_population=p100))
    rag.add_edge(Edge(EdgeClass.COALESCENCE, "absorb", p100,
                      kernel="Kabs", partner_population=p111))

    layout = StateLayout()
    layout.add_discrete("SIA111", N, population=p111)
    layout.add_discrete("SIA100", N, population=p100)
    layout.freeze()
    gw = GraphWalker(rag, layout, boundary="reflection")

    # Representative non-negative state on both ladders.
    y = np.zeros(2 * N)
    y[layout.slice("SIA111")] = np.array([1.0, 0.6, 0.4, 0.25, 0.15, 0.08, 0.03, 0.0])
    y[layout.slice("SIA100")] = np.array([0.0, 0.0, 0.2, 0.15, 0.1, 0.05, 0.0, 0.0])
    dydt = gw.assemble(0.0, y)

    sizes = np.arange(1, N + 1)
    d111 = dydt[layout.slice("SIA111")]
    d100 = dydt[layout.slice("SIA100")]
    content_rate = float(np.dot(sizes, d111) + np.dot(sizes, d100))
    check("all SIA-conserving edges compose -> total content flat",
          abs(content_rate) < ATOL, content_rate)
    check("dydt finite", np.all(np.isfinite(dydt)))
    # The conversion edges must move content 111 -> 100 (100 strictly gains).
    check("bulk-100 net content increases",
          float(np.dot(sizes, d100)) > ATOL, float(np.dot(sizes, d100)))


def main():
    print("--- Part 1: structural ---")
    test_structural()
    print("--- Part 2: composed conservation ---")
    test_composed_conservation()
    print(f"\n{'ALL PASS' if not _fails else str(_fails) + ' FAILED'}")
    return 1 if _fails else 0


if __name__ == "__main__":
    sys.exit(main())
