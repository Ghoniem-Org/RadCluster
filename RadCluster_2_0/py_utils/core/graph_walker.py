"""
core/graph_walker.py — The master equation as a graph-walking accumulator.

Host-independent.  Implements the "graph walker" of Ghoniem (2026),
*A Generalized Graph-Based Cluster Dynamics Framework*, Section 3.3:

    dc_a/dt = G_a + sum_r S_ar J_r(c) - D_a c_a

The walker visits every admissible edge of a :class:`ReactionAdmissibilityGraph`,
evaluates the reaction flux J_r from the edge's rate kernel, and deposits
the signed stoichiometric contribution S_ar J_r into the affected vertices
of the flat ODE state vector defined by a :class:`StateLayout`.

This is a *reference* discrete (one-ODE-per-size) accumulator: it carries
no host parameters and no EUROFER assumptions — it dispatches purely on
:class:`EdgeClass`.  The size-resolved rate kernels are precomputed
arrays registered on the RAG by a host material (Layer 2).  The
production solver remains the C++ core; this Python walker is the
executable specification the C++ must agree with, and the substrate for
the bin-moment "reconstruct -> walk -> project" reduction.

Sign convention.  Every per-class routine below realises exactly the
stoichiometric contract fixed in ``EDGE_CLASS_SPEC``; because those signs
are structural, summed signed-defect content q = chi*n is conserved by
construction for every intra-host reaction (sources/sinks aside).
"""
from __future__ import annotations

from typing import Dict, Tuple

import numpy as np

from .cluster_identifier import Polarity, Population
from .edge_classes import EdgeClass
from .rag import Edge, ReactionAdmissibilityGraph
from .state_layout import StateLayout


class GraphWalker:
    """Assembles the master-equation RHS by walking a RAG edge-by-edge.

    Parameters
    ----------
    rag : ReactionAdmissibilityGraph
        The host's reaction graph, with rate kernels registered as
        precomputed numpy arrays (indexed 0-based by size-1) or scalars.
    layout : StateLayout
        The ODE-index map.  Every discrete size axis must be a block
        whose ``meta['population']`` is the corresponding
        :class:`Population`.
    boundary : {'reflection', 'absorption'}
        Treatment of growth/coalescence at the top of a size axis.
        ``reflection`` (default) suppresses the over-the-top reaction;
        ``absorption`` lets it fire and drops the product.
    """

    def __init__(self, rag: ReactionAdmissibilityGraph, layout: StateLayout,
                 boundary: str = "reflection"):
        if boundary not in ("reflection", "absorption"):
            raise ValueError(f"unknown boundary mode {boundary!r}")
        self.rag = rag
        self.layout = layout
        self.boundary = boundary

        # (polarity, population name) -> discrete block name
        self._pop_block: Dict[Tuple[Polarity, str], str] = {}
        for b in layout.blocks:
            if b.kind != "discrete":
                continue
            p = b.meta.get("population")
            if isinstance(p, Population):
                self._pop_block[(p.polarity, p.name)] = b.name
        if not self._pop_block:
            raise ValueError(
                "StateLayout has no discrete blocks tagged with a Population")

        self._dispatch = {
            EdgeClass.GROWTH:           self._c_growth,
            EdgeClass.SHRINKAGE:        self._c_ladder_down,
            EdgeClass.RECOMBINATION:    self._c_ladder_down,
            EdgeClass.DISSOCIATION:     self._c_dissociation,
            EdgeClass.SOURCE:           self._c_source,
            EdgeClass.SINK:             self._c_sink,
            EdgeClass.INTER_POPULATION: self._c_inter_population,
            EdgeClass.COALESCENCE:      self._c_coalescence,
            EdgeClass.ANNIHILATION:     self._c_annihilation,
            EdgeClass.SOLUTE_TRAPPING:  self._c_solute_trapping,
        }

    # ── public API ────────────────────────────────────────────────────────
    def assemble(self, t: float, y: np.ndarray) -> np.ndarray:
        """Return dy/dt for state ``y`` at time ``t`` (the graph walk)."""
        y = np.asarray(y, dtype=float)
        if y.shape != (self.layout.N_eq,):
            raise ValueError(
                f"state vector has shape {y.shape}, expected ({self.layout.N_eq},)")
        dydt = np.zeros(self.layout.N_eq)

        # Concentration (read) and rate (accumulate) views per population.
        conc: Dict[Tuple[Polarity, str], np.ndarray] = {}
        dconc: Dict[Tuple[Polarity, str], np.ndarray] = {}
        for key, bname in self._pop_block.items():
            sl = self.layout.slice(bname)
            conc[key] = np.maximum(y[sl], 0.0)   # clamp: concentrations >= 0
            dconc[key] = dydt[sl]                # a view -> accumulates in place

        ctx = _WalkContext(conc, dconc, dydt, self.layout)
        for edge in self.rag.edges:
            self._dispatch[edge.edge_class](edge, ctx)
        return dydt

    # ── helpers ───────────────────────────────────────────────────────────
    def _key(self, p: Population) -> Tuple[Polarity, str]:
        k = (p.polarity, p.name)
        if k not in self._pop_block:
            raise KeyError(
                f"population {p.name!r} ({p.polarity.label}) has no discrete "
                "block in the StateLayout")
        return k

    def _kernel(self, edge: Edge) -> np.ndarray:
        """Return the edge's rate kernel as an array (broadcast if scalar)."""
        return np.asarray(self.rag.kernel(edge.kernel), dtype=float)

    def _monomer(self, ctx: "_WalkContext", polarity: Polarity) -> float:
        """Current concentration of the ``polarity`` point-defect monomer."""
        mp = self.rag.monomer_population(polarity)
        return float(ctx.conc[(polarity, mp.name)][0])

    def _monomer_dydt(self, ctx: "_WalkContext", polarity: Polarity) -> np.ndarray:
        mp = self.rag.monomer_population(polarity)
        return ctx.dconc[(polarity, mp.name)]

    # ── per-class contributions ───────────────────────────────────────────
    def _c_growth(self, edge: Edge, ctx: "_WalkContext") -> None:
        """(chi,n) + chi-monomer -> (chi,n+1).  Flux K[n] c_n c_1."""
        k = self._key(edge.population)
        c, dc = ctx.conc[k], ctx.dconc[k]
        chi = edge.population.polarity
        K = _as_size_array(self._kernel(edge), c.size)
        c1 = self._monomer(ctx, chi)
        flux = K * c * c1
        if self.boundary == "reflection":
            flux[-1] = 0.0                       # no I_{Nmax}+I_1 -> I_{Nmax+1}
        dc -= flux                               # loss at n
        dc[1:] += flux[:-1]                      # gain at n+1
        self._monomer_dydt(ctx, chi)[0] -= flux.sum()   # one monomer per event

    def _c_ladder_down(self, edge: Edge, ctx: "_WalkContext") -> None:
        """(chi,n) + (-chi)-monomer -> (chi,n-1).  SHRINKAGE / RECOMBINATION.

        The opposite-polarity monomer is consumed (annihilated).  At n=1
        the product is the empty cluster (pure recombination).
        """
        k = self._key(edge.population)
        c, dc = ctx.conc[k], ctx.dconc[k]
        chi = edge.population.polarity
        K = _as_size_array(self._kernel(edge), c.size)
        c1_opp = self._monomer(ctx, chi.opposite)
        flux = K * c * c1_opp
        dc -= flux                               # loss at n
        dc[:-1] += flux[1:]                       # gain at n-1 (n>=2)
        self._monomer_dydt(ctx, chi.opposite)[0] -= flux.sum()  # -chi annihilated

    def _c_dissociation(self, edge: Edge, ctx: "_WalkContext") -> None:
        """(chi,n) -> (chi,n-1) + chi-monomer.  First-order thermal emission."""
        k = self._key(edge.population)
        c, dc = ctx.conc[k], ctx.dconc[k]
        chi = edge.population.polarity
        lam = _as_size_array(self._kernel(edge), c.size)
        flux = lam * c
        flux[0] = 0.0                            # a monomer cannot dissociate
        dc -= flux                               # loss at n
        dc[:-1] += flux[1:]                       # residual cluster at n-1
        self._monomer_dydt(ctx, chi)[0] += flux.sum()   # emitted monomer pooled

    def _c_source(self, edge: Edge, ctx: "_WalkContext") -> None:
        """(empty) -> (chi,n).  Cascade / transmutation injection G_n."""
        k = self._key(edge.population)
        ctx.dconc[k] += _as_size_array(self._kernel(edge), ctx.conc[k].size)

    def _c_sink(self, edge: Edge, ctx: "_WalkContext") -> None:
        """(chi,n) -> (sink).  First-order loss D_n c_n at a fixed sink."""
        k = self._key(edge.population)
        D = _as_size_array(self._kernel(edge), ctx.conc[k].size)
        ctx.dconc[k] -= D * ctx.conc[k]

    def _c_inter_population(self, edge: Edge, ctx: "_WalkContext") -> None:
        """(chi,n,p) -> (chi,n,p').  Population transfer, size fixed."""
        ks, kd = self._key(edge.population), self._key(edge.product_population)
        kap = _as_size_array(self._kernel(edge), ctx.conc[ks].size)
        flux = kap * ctx.conc[ks]
        ctx.dconc[ks] -= flux
        ctx.dconc[kd][:flux.size] += flux

    def _c_coalescence(self, edge: Edge, ctx: "_WalkContext") -> None:
        """(chi,n) + (chi,n') -> (chi,n+n').  Same-polarity binary growth.

        Reference double sum; the kernel is a 2-D array K[n-1, n'-1].
        """
        ka = self._key(edge.population)
        kb = self._key(edge.partner_population)
        ca, cb = ctx.conc[ka], ctx.conc[kb]
        da, db = ctx.dconc[ka], ctx.dconc[kb]
        K = np.atleast_2d(self._kernel(edge))
        same = (ka == kb)
        for i in range(ca.size):
            if ca[i] <= 0.0:
                continue
            for j in range(cb.size):
                ij = i + j + 1                   # product size index (n+n'-1)
                if ij >= da.size:
                    if self.boundary == "reflection":
                        continue
                    ij = -1                      # absorption: drop product
                rate = K[i, j] * ca[i] * cb[j]
                if same:
                    rate *= 0.5                  # avoid double counting
                da[i] -= rate
                db[j] -= rate
                if ij >= 0:
                    da[ij] += 2.0 * rate if same else rate

    def _c_annihilation(self, edge: Edge, ctx: "_WalkContext") -> None:
        """(chi,n) + (-chi,n') -> survivor of size |n-n'|.  Reference sum."""
        ka = self._key(edge.population)              # chi side
        kb = self._key(edge.partner_population)      # -chi side
        ca, cb = ctx.conc[ka], ctx.conc[kb]
        da, db = ctx.dconc[ka], ctx.dconc[kb]
        K = np.atleast_2d(self._kernel(edge))
        for i in range(ca.size):
            if ca[i] <= 0.0:
                continue
            for j in range(cb.size):
                rate = K[i, j] * ca[i] * cb[j]
                if rate <= 0.0:
                    continue
                da[i] -= rate
                db[j] -= rate
                d = (i + 1) - (j + 1)            # signed size difference
                if d > 0 and d - 1 < da.size:
                    da[d - 1] += rate            # chi survivor
                elif d < 0 and (-d) - 1 < db.size:
                    db[-d - 1] += rate           # -chi survivor
                # d == 0: mutual annihilation, no survivor

    def _c_solute_trapping(self, edge: Edge, ctx: "_WalkContext") -> None:
        """(chi,n,c) <-> (chi,n,c +/- e_s).

        The composition axis is not materialised by the discrete reference
        walker; in the EUROFER instantiation helium is carried by the
        mean-loading reduction (Case 1/2), which owns this edge class.
        Reached only if a host wires an explicit composition axis.
        """
        raise NotImplementedError(
            "SOLUTE_TRAPPING requires an explicit composition axis; in the "
            "EUROFER instantiation it is handled by the He mean-loading "
            "reduction (core/reductions). Edge: " + edge.label)


class _WalkContext:
    """Mutable scratch passed to the per-edge contribution routines."""

    __slots__ = ("conc", "dconc", "dydt", "layout")

    def __init__(self, conc, dconc, dydt, layout):
        self.conc = conc
        self.dconc = dconc
        self.dydt = dydt
        self.layout = layout


def _as_size_array(kernel, n: int) -> np.ndarray:
    """Broadcast a scalar kernel to length ``n``; pass an array through."""
    arr = np.asarray(kernel, dtype=float)
    if arr.ndim == 0:
        return np.full(n, float(arr))
    if arr.shape[0] != n:
        raise ValueError(
            f"kernel length {arr.shape[0]} does not match size axis {n}")
    return arr
