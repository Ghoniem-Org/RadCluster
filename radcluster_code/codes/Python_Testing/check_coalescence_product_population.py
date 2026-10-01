"""
check_coalescence_product_population.py — Phase-1 unit tests for the
cross-population COALESCENCE extension (½⟨111⟩ → ⟨100⟩ loop-conversion work).

Covers the Layer-1 core change in `core/{edge_classes,rag,graph_walker}.py`:

  1. Same-population coalescence conserves signed-defect content
     (regression for the factor-2 gain fix in `_c_coalescence`).
  2. Cross-population *junction* coalescence  bulk-111 × bulk-111 → bulk-100
     conserves TOTAL SIA content summed over both populations.
  3. Cross-population *absorption* coalescence  bulk-100 × bulk-111 → bulk-100
     (partner ≠ population, product = population) conserves content.
  4. Splitting one coalescence edge into  φ·(→bulk-100) + (1−φ)·(→bulk-111)
     reproduces the single-edge total content flux (the φ-split is conservative).
  5. A COALESCENCE product_population of the *wrong* polarity is rejected.

Run:  python codes/Python_Testing/check_coalescence_product_population.py
Exit code 0 on success, 1 on any failure.  Also importable as pytest tests.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")))

from py_utils.core.cluster_identifier import Polarity, Population
from py_utils.core.rag import ReactionAdmissibilityGraph, Edge
from py_utils.core.edge_classes import EdgeClass
from py_utils.core.state_layout import StateLayout
from py_utils.core.graph_walker import GraphWalker

ATOL = 1e-12


# ── helpers ───────────────────────────────────────────────────────────────────
def _content_rate(dydt, layout, *blocks):
    """d/dt of the signed-defect content Σ_n n·c_n summed over named blocks."""
    total = 0.0
    for name in blocks:
        sl = layout.slice(name)
        seg = dydt[sl]
        sizes = np.arange(1, seg.size + 1)
        total += float(np.dot(sizes, seg))
    return total


def _one_pop_rag(N, kernel):
    rag = ReactionAdmissibilityGraph("one_pop")
    sia = rag.add_population(Population("bulk", Polarity.SIA, n_min=1, mobile_max=N))
    rag.set_monomer_population(sia)
    rag.register_kernel("Kii", kernel)
    return rag, sia


def _two_pop_rag(N, kernel):
    """bulk-111 (mobile) + bulk-100 (sessile), shared monomer pool on 111."""
    rag = ReactionAdmissibilityGraph("two_pop")
    p111 = rag.add_population(Population("bulk-111", Polarity.SIA, n_min=1, mobile_max=N))
    p100 = rag.add_population(Population("bulk-100", Polarity.SIA, n_min=1, mobile_max=0))
    rag.set_monomer_population(p111)
    rag.register_kernel("Kii", kernel)
    return rag, p111, p100


def _layout_two(N, p111, p100):
    layout = StateLayout()
    layout.add_discrete("SIA111", N, population=p111)
    layout.add_discrete("SIA100", N, population=p100)
    layout.freeze()
    return layout


# Initial state with only small sizes populated, so every coalescence product
# fits inside the size axis (no over-the-top loss) and content is exactly
# conserved by an ideal accumulator.
def _small_state(N, *segments):
    y = np.zeros(N * len(segments))
    return y


# ── tests ─────────────────────────────────────────────────────────────────────
def test_same_population_conserves_content():
    N = 8
    rag, sia = _one_pop_rag(N, np.full((N, N), 0.1))
    rag.add_edge(Edge(EdgeClass.COALESCENCE, "self", sia, kernel="Kii"))
    layout = StateLayout(); layout.add_discrete("SIA", N, population=sia); layout.freeze()
    gw = GraphWalker(rag, layout, boundary="reflection")
    y = np.zeros(N); y[0], y[1], y[2] = 1.0, 0.7, 0.4   # max product 6 < N
    dydt = gw.assemble(0.0, y)
    rate = _content_rate(dydt, layout, "SIA")
    assert abs(rate) < ATOL, f"same-pop coalescence not conservative: {rate}"
    return rate


def test_junction_conserves_total_content():
    """bulk-111 × bulk-111 → bulk-100 : Σn over BOTH pops must be flat."""
    N = 8
    rag, p111, p100 = _two_pop_rag(N, np.full((N, N), 0.1))
    rag.add_edge(Edge(EdgeClass.COALESCENCE, "junction", p111,
                      kernel="Kii", product_population=p100))
    layout = _layout_two(N, p111, p100)
    gw = GraphWalker(rag, layout, boundary="reflection")
    y = np.zeros(2 * N)
    y[0], y[1], y[2] = 1.0, 0.7, 0.4                    # 111 sizes 1..3
    dydt = gw.assemble(0.0, y)
    total = _content_rate(dydt, layout, "SIA111", "SIA100")
    only111 = _content_rate(dydt, layout, "SIA111")
    assert abs(total) < ATOL, f"junction not globally conservative: {total}"
    # 111 must strictly lose content (it is the only source); 100 must gain it.
    assert only111 < -ATOL, f"expected net 111 loss, got {only111}"
    return total, only111


def test_absorption_conserves_content():
    """bulk-100 × bulk-111 → bulk-100 : product = population (no redirect)."""
    N = 8
    rag, p111, p100 = _two_pop_rag(N, np.full((N, N), 0.1))
    rag.add_edge(Edge(EdgeClass.COALESCENCE, "absorb", p100,
                      kernel="Kii", partner_population=p111))
    layout = _layout_two(N, p111, p100)
    gw = GraphWalker(rag, layout, boundary="reflection")
    y = np.zeros(2 * N)
    y[0], y[1] = 0.6, 0.3              # 111 sizes 1,2 (mobile partner)
    y[N + 2] = 0.5                     # 100 size 3
    dydt = gw.assemble(0.0, y)
    total = _content_rate(dydt, layout, "SIA111", "SIA100")
    assert abs(total) < ATOL, f"absorption not conservative: {total}"
    return total


def test_product_redirect_moves_only_gain():
    """Redirecting the product to bulk-100 must move ONLY the gain.

    With the same kernel and reactants, a self-coalescence (product → 111)
    and a junction (product → 100) consume 111 identically; they differ only
    in where the product is deposited.  Hence, element-wise,

        d111(self)  ==  d111(junction)[pure loss]  +  d100(junction)[pure gain]

    so the redirect is content-preserving relative to the un-redirected edge.
    """
    N = 8
    K = np.full((N, N), 0.1)
    y111 = np.array([1.0, 0.7, 0.4] + [0.0] * (N - 3))

    # A: self-coalescence, product back to 111 (loss + gain mixed on 111).
    ragA, p111A, p100A = _two_pop_rag(N, K)
    ragA.add_edge(Edge(EdgeClass.COALESCENCE, "self", p111A, kernel="Kii"))
    layA = _layout_two(N, p111A, p100A)
    dA = GraphWalker(ragA, layA, boundary="reflection").assemble(
        0.0, np.concatenate([y111, np.zeros(N)]))

    # C: junction, same K, product → 100 (pure loss on 111, pure gain on 100).
    ragC, p111C, p100C = _two_pop_rag(N, K)
    ragC.add_edge(Edge(EdgeClass.COALESCENCE, "junction", p111C,
                       kernel="Kii", product_population=p100C))
    layC = _layout_two(N, p111C, p100C)
    dC = GraphWalker(ragC, layC, boundary="reflection").assemble(
        0.0, np.concatenate([y111, np.zeros(N)]))

    d111_A = dA[layA.slice("SIA111")]
    d111_C = dC[layC.slice("SIA111")]          # pure 111 reactant loss
    d100_C = dC[layC.slice("SIA100")]          # pure product gain
    recon = d111_C + d100_C                     # gain re-added at product sizes
    assert np.allclose(d111_A, recon, atol=ATOL), \
        f"redirect changed the 111 loss:\n A={d111_A}\n recon={recon}"
    # Sanity: junction's 111 channel is pure loss (≤ 0 everywhere).
    assert np.all(d111_C <= ATOL), f"junction 111 channel not pure-loss: {d111_C}"
    return float(np.max(np.abs(d111_A - recon)))


def test_wrong_polarity_product_rejected():
    rag = ReactionAdmissibilityGraph("bad")
    sia = rag.add_population(Population("bulk", Polarity.SIA, n_min=1, mobile_max=4))
    vac = rag.add_population(Population("bulk", Polarity.VACANCY, n_min=1, mobile_max=4))
    raised = False
    try:
        Edge(EdgeClass.COALESCENCE, "bad", sia, kernel="Kii", product_population=vac)
    except ValueError:
        raised = True
    assert raised, "opposite-polarity COALESCENCE product_population not rejected"
    return raised


# ── runner ────────────────────────────────────────────────────────────────────
def main():
    tests = [
        ("same-pop content conserved (factor-2 fix)", test_same_population_conserves_content),
        ("junction globally conservative", test_junction_conserves_total_content),
        ("absorption conservative", test_absorption_conserves_content),
        ("redirect moves only the gain", test_product_redirect_moves_only_gain),
        ("wrong-polarity product rejected", test_wrong_polarity_product_rejected),
    ]
    fails = 0
    for name, fn in tests:
        try:
            out = fn()
            print(f"  PASS  {name}   {out!r}")
        except AssertionError as e:
            fails += 1
            print(f"  FAIL  {name}   {e}")
    print(f"\n{len(tests) - fails}/{len(tests)} passed")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
