"""
materials/eurofer97/declaration.py — EUROFER-97 RAG declaration (Layer 2).

Builds the host-specific reaction admissibility graph G_Eur and its
contiguous ODE-index :class:`StateLayout` for the reduced-activation
ferritic-martensitic steel EUROFER-97, following Ghoniem (2026),
*A Generalized Graph-Based Cluster Dynamics Framework for Irradiated
Materials*, Section 4.

This module is a *declaration*, not an integrator.  It

  * declares the EUROFER-97 populations (Eqs. 42-46): one bulk SIA
    population and one bulk vacancy population, plus helium as the single
    resolved gas (Eq. 41);
  * adds one :class:`Edge` family per EUROFER process class P1-P8
    (Section 4.3), each mapped to one of the ten abstract
    :class:`EdgeClass` values;
  * registers, for every edge, the size-resolved rate kernel as the
    *precomputed numpy array* already built by :class:`ReactionRates`
    (the abstract core never re-derives a rate); and
  * builds the :class:`StateLayout` that packs the SIA ladder, the
    vacancy ladder, the helium block, and a small conservation-accounting
    aux block into one flat ODE state vector.

The production solver (C++/CVODE) remains the engine; this layer is the
"subclass" in the paper's abstract-base-class analogy — it says *what*
EUROFER-97 is, never *how* to integrate it.
"""
from __future__ import annotations

import numpy as np

from ...core import (
    Polarity,
    Population,
    EdgeClass,
    Edge,
    ReactionAdmissibilityGraph,
    StateLayout,
    BinMomentReduction,
    HeReductionMode,
)
from ...defect_production import production_rates
from ...binding_energies import Gamma_TM, Gamma_res, ell_max


# ─────────────────────────────────────────────────────────────────────────────
# Rate-kernel construction helpers
#
# ReactionRates precomputes 1-D size-resolved arrays (0-indexed by size-1) and
# a handful of scalars.  The abstract GROWTH / SHRINKAGE / DISSOCIATION / SINK
# / SOURCE edge classes consume a per-size array (or a broadcastable scalar)
# directly.  COALESCENCE and ANNIHILATION need a 2-D kernel K[n-1, n'-1]; we
# assemble those here from the precomputed effective diffusivities so that no
# rate physics is re-derived — only the precomputed pieces are recombined into
# the shape the graph walker expects (paper Eqs. 79-81).
# ─────────────────────────────────────────────────────────────────────────────

def _xi(n: np.ndarray) -> np.ndarray:
    """Capture radius in lattice units, xi_n = (3 n / 8 pi)^(1/3) (paper Eq. 79)."""
    return (3.0 * n / (8.0 * np.pi)) ** (1.0 / 3.0)


def _coalescence_kernel(D_eff: np.ndarray, A_pref: float) -> np.ndarray:
    """Same-polarity binary-growth kernel K[n-1, n'-1] (paper Eqs. 79-80).

    K_{n,n'} = 8 pi (xi_n + xi_n') (D_n + D_n') / Omega^(2/3)

    ``D_eff`` is the precomputed effective 3-D diffusivity array (mobile
    sizes non-zero, sessile sizes zero) from :class:`ReactionRates`;
    ``A_pref`` is ``8 pi / Omega^(2/3)`` assembled from the precomputed
    ``A_sph_inv_O23`` geometric prefactor.  Sessile-sessile pairs vanish
    automatically because their diffusivities are zero.
    """
    n = np.arange(1, D_eff.size + 1, dtype=float)
    xi = _xi(n)
    xi_sum = xi[:, None] + xi[None, :]
    D_sum = D_eff[:, None] + D_eff[None, :]
    return A_pref * xi_sum * D_sum


def _absorb_kernel(D_partner: np.ndarray, A_pref: float) -> np.ndarray:
    """Cross-character absorption kernel K[m-1, n-1] (Marian growth, P3-like).

    For ``<100>_m + <111>_n -> <100>_{m+n}`` (paper Eqs. 79/81 form):

      K_{m,n} = 8 pi (xi_m + xi_n) (D^{100}_m + D^{111}_n) / Omega^(2/3)

    The sessile ⟨100⟩ primary contributes zero diffusivity, so only the mobile
    ½⟨111⟩ partner ``D_partner = D_SIA_eff`` drives the capture.  Primary axis
    m = ⟨100⟩ size, partner axis n = ½⟨111⟩ size (both span 1..I).
    """
    n = np.arange(1, D_partner.size + 1, dtype=float)
    xi = _xi(n)
    xi_sum = xi[:, None] + xi[None, :]
    return A_pref * xi_sum * D_partner[None, :]


def _annihilation_kernel(D_v: np.ndarray, D_i: np.ndarray,
                         A_pref: float) -> np.ndarray:
    """Cross-polarity V-I annihilation kernel K[m-1, n-1] (paper Eq. 81).

    K_{m,n} = 8 pi (xi_m + xi_n) (D_m^v + D_n^i) / Omega^(2/3)

    The primary axis is the vacancy population (size m), the partner axis
    is the SIA population (size n), matching the EUROFER ANNIHILATION edge
    declared below (vacancy population as ``population``, SIA population as
    ``partner_population``).
    """
    m = np.arange(1, D_v.size + 1, dtype=float)
    n = np.arange(1, D_i.size + 1, dtype=float)
    xi_sum = _xi(m)[:, None] + _xi(n)[None, :]
    D_sum = D_v[:, None] + D_i[None, :]
    return A_pref * xi_sum * D_sum


def _source_array(P_indexed: np.ndarray) -> np.ndarray:
    """Convert a 1-indexed production array P[size] to a 0-indexed size array.

    :func:`production_rates` returns arrays of length ``n_max + 1`` whose
    index 0 is unused and index k holds the size-k production rate.  The
    SOURCE edge class consumes a 0-indexed size array (size n -> [n-1]).
    """
    return np.asarray(P_indexed[1:], dtype=float)


def _trap_mutation_kernel(V: int, T: float) -> np.ndarray:
    """Size-resolved trap-mutation rate kernel Gamma_TM[m-1] (paper Eq. 75 / P7).

    Trap mutation V_m(He) -> V_{m+1} + I_1 is gas-pressure-driven; its rate
    ``Gamma_TM(m, ell)`` (binding_energies, paper Eq. 75) depends on both the
    vacancy size ``m`` and the helium loading ``ell``.  This declaration layer
    carries one bulk vacancy population without an explicit (m, ell)
    occupancy axis, so the loading-resolved barrier table cannot be collapsed
    to a 1-D size array without picking a representative ``ell`` per ``m``.

    We register the *spontaneous-mutation envelope*: for each vacancy size m
    we evaluate ``Gamma_TM`` at the maximum admissible helium loading
    ``ell = ell_max(m)`` (the trap-mutation limit, paper Table 27).  Where
    that (m, ell) pair is absent from the barrier table (``E_TM`` table only
    populates the small over-pressurized bubbles He5-7V1 / He4-6), the rate
    is exactly zero — so the returned array is non-zero only on the handful
    of sizes for which an atomistic trap-mutation barrier is known and zero
    elsewhere.  This is a faithful placeholder: the C++ solver does not yet
    evaluate P7 (Stage-2 review), so no run consumes a non-zero value here;
    the array documents *where* trap mutation is admissible.
    """
    m = np.arange(1, V + 1, dtype=int)
    return np.array([Gamma_TM(int(mi), ell_max(int(mi)), T) for mi in m],
                    dtype=float)


def _resolution_kernel(V: int, phi_dot: float, b0: float) -> np.ndarray:
    """Size-resolved radiation re-solution rate kernel Gamma_res[m-1] (P8).

    Radiation re-solution He_ell V_m -> He_{ell-1} V_m + He detraps one
    helium atom per displacement event; its rate ``Gamma_res(ell) = b0 *
    ell * phi_dot`` (binding_energies, paper Eq. 75 / Section 5.3) is
    proportional to the helium loading ``ell``, not the vacancy size ``m``.

    Without an explicit (m, ell) occupancy axis the per-``ell`` factor
    cannot be resolved here; we register the *per-He-atom* re-solution
    coefficient ``b0 * phi_dot`` broadcast over the vacancy size axis (so
    the C++ kernel arithmetic multiplies it by the local loading ell when
    P8 is implemented).  If the displacement rate or b0 is unavailable the
    array degrades to all-zero, which cleanly disables P8.
    """
    coeff = float(b0) * float(phi_dot)
    return np.full(V, coeff, dtype=float)


# ─────────────────────────────────────────────────────────────────────────────
# The EUROFER-97 RAG declaration
# ─────────────────────────────────────────────────────────────────────────────

def build_eurofer_rag(input_data, reaction_rates, *,
                      equations: str = "discrete",
                      cascade: str = "fission"):
    """Build the EUROFER-97 reaction admissibility graph and state layout.

    Parameters
    ----------
    input_data : InputData
        The loaded EUROFER-97 parameter set (sizes ``I``/``V``, mobility
        cutoffs ``i_mobile``/``v_mobile``, derived physics).
    reaction_rates : ReactionRates
        The precomputed rate-constant arrays built from ``input_data``.
        Its attributes are registered verbatim as the RAG kernel library.
    equations : {'discrete', 'bin_moment'}
        Selects the SIA / vacancy state-space representation: one ODE per
        size (``discrete``) or the logarithmic bin-moment reduction
        (``bin_moment``, paper Section 5.4).
    cascade : {'fission', 'fusion'}
        Selects the cascade production spectrum and, through
        :meth:`HeReductionMode.from_cascade`, the helium reduction
        (fission -> Case 2 decoupled; fusion -> Case 1 mean-field).

    Returns
    -------
    (ReactionAdmissibilityGraph, StateLayout)
        The host-specific graph G_Eur with all kernels registered, and the
        contiguous ODE-index layout that packs it.
    """
    equations = str(equations).lower()
    if equations not in ("discrete", "bin_moment"):
        raise ValueError(
            f"equations must be 'discrete' or 'bin_moment', got {equations!r}")

    I = int(input_data.I)                 # max SIA cluster size
    V = int(input_data.V)                 # max vacancy cluster size
    i_mobile = int(input_data.i_mobile)   # max mobile SIA size
    v_mobile = int(input_data.v_mobile)   # max mobile vacancy size

    # ── 1. The graph and its populations (paper Eqs. 41-46) ──────────────────
    rag = ReactionAdmissibilityGraph("EUROFER-97", gas_species=["He"])

    # Two SIA populations + one vacancy population (paper Eqs. 42-43, Eq. 12).
    #
    #   bulk-111 : glissile ½⟨111⟩ loops + 3D-mobile small clusters (n < 4);
    #              holds the SIA monomer pool; receives ALL cascade SIA
    #              production; self-coalesces (split by phi_junc).
    #   bulk-100 : sessile ⟨100⟩ loops (mobile_max = 0); no glide, no
    #              self-coalescence, no cascade source.  n_min is the
    #              loop-onset floor (no ⟨100⟩ loop exists below it).
    #
    # The two are coupled by the ½⟨111⟩ → ⟨100⟩ conversion edges declared in
    # Section 3: a unary INTER_POPULATION transformation (Dudarev), a binary
    # junction COALESCENCE (Marian), and a ½⟨111⟩-absorption COALESCENCE
    # (Marian growth).  See docs/design_notes/loop_111_to_100_conversion.md.
    n_loop_min = int(getattr(reaction_rates, "n_loop_min", 4))
    sia_111 = rag.add_population(
        Population("bulk-111", Polarity.SIA, n_min=1, mobile_max=i_mobile))
    sia_100 = rag.add_population(
        Population("bulk-100", Polarity.SIA, n_min=n_loop_min, mobile_max=0))
    vac_bulk = rag.add_population(
        Population("bulk", Polarity.VACANCY, n_min=1, mobile_max=v_mobile))

    # The size-1 vertex of bulk-111 is the SIA point-defect monomer pool; the
    # sessile bulk-100 loops draw/emit monomers from it automatically (the
    # walker resolves the monomer by polarity).  Vacancy monomer pool is the
    # vacancy bulk.  GROWTH / SHRINKAGE / RECOMBINATION / DISSOCIATION consume.
    rag.set_monomer_population(sia_111)
    rag.set_monomer_population(vac_bulk)

    # ── 2. Rate-kernel library k (precomputed by ReactionRates) ──────────────
    # Every kernel below is an array/scalar already built by ReactionRates;
    # the abstract core never evaluates a rate formula.
    rr = reaction_rates

    # 8 pi / Omega^(2/3): coalescence/annihilation geometric prefactor.
    # ReactionRates exposes A_sph / Omega^(2/3) as A_sph_inv_O23, with
    # A_sph = (48 pi^2)^(1/3); hence 8 pi / Omega^(2/3) = A_sph_inv_O23 *
    # 8 pi / A_sph.
    A_sph = float(input_data.derived["A_sph"])
    A_8pi = rr.A_sph_inv_O23 * (8.0 * np.pi / A_sph)

    # P1 — V-SIA recombination kernel (scalar K_iv, paper Eqs. 52, 78).
    rag.register_kernel("K_iv", float(rr.K_iv))

    # P2/P3 — monomer absorption (growth / shrinkage ladders).
    #   SIA loop grows by absorbing an SIA monomer       -> K_SIA_grow
    #   SIA loop shrinks by absorbing a vacancy monomer  -> K_SIA_shrink
    #   void grows by absorbing a vacancy monomer        -> K_VAC_grow
    #   void shrinks by absorbing an SIA monomer         -> K_VAC_shrink
    rag.register_kernel("K_SIA_grow", np.asarray(rr.K_SIA_grow, dtype=float))
    rag.register_kernel("K_SIA_shrink", np.asarray(rr.K_SIA_shrink, dtype=float))
    rag.register_kernel("K_VAC_grow", np.asarray(rr.K_VAC_grow, dtype=float))
    rag.register_kernel("K_VAC_shrink", np.asarray(rr.K_VAC_shrink, dtype=float))

    # P4 — fixed-sink loss D_alpha (dislocations + GB + precipitates).
    #   SIA: size-resolved array k2_SIA (mobile-only; 0 for large loops).
    #   Vacancy: scalar k2_vac_scalar broadcast over the size axis.
    rag.register_kernel("D_SIA_sink", np.asarray(rr.k2_SIA, dtype=float))
    rag.register_kernel("D_VAC_sink", float(rr.k2_vac_scalar))

    # P5 — thermal monomer emission (dissociation ladders).
    #   SIA loop emits an SIA monomer  -> G_SIA  (alpha_i(n))
    #   void emits a vacancy monomer   -> G_VAC  (alpha_v(m))
    rag.register_kernel("eps_SIA_emit", np.asarray(rr.G_SIA, dtype=float))
    rag.register_kernel("eps_VAC_emit", np.asarray(rr.G_VAC, dtype=float))

    # SIA-SIA and V-V coalescence — 2-D same-polarity binary kernels
    # assembled from the precomputed effective diffusivities (Eqs. 79-80).
    #
    # The ½⟨111⟩-½⟨111⟩ collision kernel K_ii is *split* by the Marian
    # size-comparability branching phi_junc[n,n'] (built by ReactionRates):
    #   K_111_self     = (1 - phi) * K_ii   ->  product stays in bulk-111
    #   K_111_junction =      phi  * K_ii   ->  product is a ⟨100⟩ loop
    # so the total ½⟨111⟩ collision rate is conserved (no double counting).
    # These 2-D same-/cross-polarity kernels are O(I²).  They are consumed
    # ONLY by the Python GraphWalker (reference RHS); the production C++ solver
    # computes K_ii, K_vv, K_vi, and φ_junc on the fly per (n,n') pair and never
    # reads them.  Register them as zero-argument lazy builders so the matrices
    # are materialised only if the GraphWalker actually evaluates the edge —
    # never on a C++ run, where materialising [I,I] at I~1e5 would need ~80 GB.
    D_SIA_eff = np.asarray(rr.D_SIA_eff, dtype=float)
    D_VAC_eff = np.asarray(rr.D_VAC_eff, dtype=float)
    rag.register_kernel(
        "K_111_self",
        lambda: (1.0 - np.asarray(rr.phi_junc, dtype=float))
                * _coalescence_kernel(D_SIA_eff, A_8pi))
    rag.register_kernel(
        "K_111_junction",
        lambda: np.asarray(rr.phi_junc, dtype=float)
                * _coalescence_kernel(D_SIA_eff, A_8pi))
    rag.register_kernel("K_vv_coal",
                        lambda: _coalescence_kernel(D_VAC_eff, A_8pi))

    # ── ½⟨111⟩ → ⟨100⟩ conversion + sessile-⟨100⟩ kernels (loop-conversion) ──
    # Unary (Dudarev) transformation rate Γ_uni(n) (1-D, INTER_POPULATION edge).
    rag.register_kernel("Gamma_uni", np.asarray(rr.Gamma_uni, dtype=float))
    # Marian absorption  ⟨100⟩_m + ½⟨111⟩_n -> ⟨100⟩_{m+n}  (cross-character),
    # gated by the two-step success probability P_success(T) (Marian Fig. 3) —
    # the absorbed ½⟨111⟩ must rotate through the metastable ½⟨110⟩ to ⟨100⟩.
    rag.register_kernel(
        "K_100_absorb",
        lambda: float(rr.conv_psuccess) * _absorb_kernel(D_SIA_eff, A_8pi))
    # Sessile ⟨100⟩ point-defect ladders (loop geometry; monomer-driven).
    rag.register_kernel("K_100_grow",   np.asarray(rr.K_100_grow, dtype=float))
    rag.register_kernel("K_100_shrink", np.asarray(rr.K_100_shrink, dtype=float))
    rag.register_kernel("eps_100_emit", np.asarray(rr.G_100, dtype=float))
    rag.register_kernel("D_100_sink",   np.asarray(rr.k2_100, dtype=float))
    # Vacancy-cluster annihilation of sessile ⟨100⟩ loops (D^{100} = 0).
    I_size = int(input_data.I)
    rag.register_kernel(
        "K_vi_annih_100",
        lambda: _annihilation_kernel(D_VAC_eff, np.zeros(I_size), A_8pi))

    # V-I cluster annihilation — cross-polarity 2-D kernel K[m-1, n-1]
    # (paper Eq. 81).  Vacancy axis primary, SIA axis partner.
    rag.register_kernel(
        "K_vi_annih",
        lambda: _annihilation_kernel(D_VAC_eff, D_SIA_eff, A_8pi))

    # Cascade production sources (paper Eqs. 11-13).  production_rates
    # returns 1-indexed arrays; convert to 0-indexed size arrays.
    G = float(input_data.derived["G"])
    G_He_r = float(input_data.derived["G_He_r"])
    Pr_SIA, Pr_VAC, _G_He = production_rates(G, cascade, I, V, G_He_r)
    rag.register_kernel("G_SIA_cascade", _source_array(Pr_SIA))
    rag.register_kernel("G_VAC_cascade", _source_array(Pr_VAC))

    # P7 — trap-mutation kernel Gamma_TM (paper Eq. 75).  Built here from
    # binding_energies.Gamma_TM because ReactionRates does not expose a
    # size/loading-resolved trap-mutation array.  The kernel is non-zero
    # only on the sizes for which an atomistic E_TM barrier is tabulated
    # (Table 27) and zero elsewhere — see _trap_mutation_kernel.
    T = float(input_data.derived["T"])
    rag.register_kernel("Gamma_TM", _trap_mutation_kernel(V, T))

    # P8 — radiation re-solution kernel Gamma_res / b0 (paper Eq. 75,
    # Section 5.3).  Gamma_res(ell) = b0 * ell * phi_dot is per-He-atom;
    # we register the per-He-atom coefficient b0 * phi_dot broadcast over
    # the vacancy size axis (the loading factor ell is applied by the
    # C++ kernel arithmetic when P8 is implemented).
    spec = str(input_data.derived.get("spectrum", cascade)).lower()
    b0_key = "b0_fission" if "fiss" in spec else "b0_fusion"
    b0_res = float(input_data.reactions.get(
        b0_key, 0.01 if "fiss" in spec else 0.10))
    rag.register_kernel("Gamma_res", _resolution_kernel(V, G, b0_res))

    # ── 3. Edge families (one per EUROFER process class P1-P8) ───────────────
    # Each Edge is a *family*: the graph walker applies its kernel vectorised
    # over the whole size axis, so one Edge stands for O(N) scalar reactions.

    # P1 — V-SIA recombination (cross-polarity, RECOMBINATION).
    # Walked on the SIA ladder: I_n + V_1 -> I_{n-1}, consuming a vacancy
    # monomer.  The n=1 step (I_1 + V_1 -> empty) is pure recombination.
    rag.add_edge(Edge(
        EdgeClass.RECOMBINATION, "P1_recombination", sia_111,
        kernel="K_iv",
        meta={"process": "P1", "note": "V-SIA recombination, paper Eq. 52"}))

    # P2/P3 — monomer absorption: cavity & loop GROWTH and SHRINKAGE.
    # Vacancy side: void absorbs a vacancy monomer (grow, P2v) or an SIA
    # monomer (shrink, P2i/P6).
    rag.add_edge(Edge(
        EdgeClass.GROWTH, "P2v_cavity_growth", vac_bulk,
        kernel="K_VAC_grow",
        meta={"process": "P2", "note": "void + V_1 -> void_{m+1}"}))
    rag.add_edge(Edge(
        EdgeClass.SHRINKAGE, "P2i_cavity_shrink", vac_bulk,
        kernel="K_VAC_shrink",
        meta={"process": "P2/P6", "note": "void + I_1 -> void_{m-1}"}))
    # SIA side: loop absorbs an SIA monomer (grow, P3) or a vacancy
    # monomer (shrink, P3 vacancy side).
    rag.add_edge(Edge(
        EdgeClass.GROWTH, "P3_loop_growth", sia_111,
        kernel="K_SIA_grow",
        meta={"process": "P3", "note": "loop + I_1 -> loop_{n+1}"}))
    rag.add_edge(Edge(
        EdgeClass.SHRINKAGE, "P3_loop_shrink", sia_111,
        kernel="K_SIA_shrink",
        meta={"process": "P3", "note": "loop + V_1 -> loop_{n-1}"}))

    # P4 — absorption at fixed unresolved sinks (SINK): dislocation
    # network, grain boundaries, MX/M23C6 precipitates.
    rag.add_edge(Edge(
        EdgeClass.SINK, "P4_SIA_sink", sia_111,
        kernel="D_SIA_sink",
        meta={"process": "P4", "note": "SIA fixed-sink loss D_i"}))
    rag.add_edge(Edge(
        EdgeClass.SINK, "P4_VAC_sink", vac_bulk,
        kernel="D_VAC_sink",
        meta={"process": "P4", "note": "vacancy fixed-sink loss D_v"}))

    # P5 — thermal monomer emission (DISSOCIATION): loop emits an SIA,
    # void emits a vacancy, into the respective monomer pool.
    rag.add_edge(Edge(
        EdgeClass.DISSOCIATION, "P5i_SIA_emission", sia_111,
        kernel="eps_SIA_emit",
        meta={"process": "P5", "note": "loop_n -> loop_{n-1} + I_1"}))
    rag.add_edge(Edge(
        EdgeClass.DISSOCIATION, "P5v_VAC_emission", vac_bulk,
        kernel="eps_VAC_emit",
        meta={"process": "P5", "note": "void_m -> void_{m-1} + V_1"}))

    # SIA-SIA and V-V coalescence (COALESCENCE, same-polarity binary).
    # ½⟨111⟩ self-coalescence keeps the (1 - phi_junc) branch in bulk-111.
    rag.add_edge(Edge(
        EdgeClass.COALESCENCE, "SIA111_self_coalescence", sia_111,
        kernel="K_111_self",
        meta={"note": "I_n + I_n' -> I_{n+n'} (½<111>), (1-phi) branch, Eq. 79"}))
    rag.add_edge(Edge(
        EdgeClass.COALESCENCE, "VAC_VAC_coalescence", vac_bulk,
        kernel="K_vv_coal",
        meta={"note": "V_m + V_m' -> V_{m+m'}, paper Eq. 80"}))

    # ── ½⟨111⟩ → ⟨100⟩ conversion edges (loop-conversion work) ───────────────
    # Marian junction (binary): ½⟨111⟩_n + ½⟨111⟩_n' -> ⟨100⟩_{n+n'}.  Cross-
    # population COALESCENCE: reactants in bulk-111, product deposited into
    # bulk-100 (the phi_junc branch of the ½⟨111⟩ collision rate).
    rag.add_edge(Edge(
        EdgeClass.COALESCENCE, "SIA111_junction", sia_111,
        kernel="K_111_junction", product_population=sia_100,
        meta={"mechanism": "Marian junction",
              "note": "½<111>_n + ½<111>_n' -> <100>_{n+n'}"}))
    # Dudarev unary transformation: ½⟨111⟩_n -> ⟨100⟩_n (size-fixed, one-way),
    # INTER_POPULATION gated by the thermodynamic driving force ΔF(n,T).
    rag.add_edge(Edge(
        EdgeClass.INTER_POPULATION, "loop_111to100_unary", sia_111,
        kernel="Gamma_uni", product_population=sia_100,
        meta={"mechanism": "Dudarev unary transformation",
              "note": "½<111>_n -> <100>_n, one-way, gated by ΔF(n,T)>0"}))

    # V-I cluster annihilation (ANNIHILATION, cross-polarity binary).
    # Primary = vacancy population, partner = SIA population.
    rag.add_edge(Edge(
        EdgeClass.ANNIHILATION, "VI_annihilation", vac_bulk,
        kernel="K_vi_annih", partner_population=sia_111,
        meta={"note": "V_m + I_n -> survivor of size |m-n|, paper Eq. 81"}))

    # Cascade production (SOURCE, no precursor vertex): displacement
    # cascades inject SIA and vacancy clusters across the size spectrum.
    rag.add_edge(Edge(
        EdgeClass.SOURCE, "cascade_SIA_source", sia_111,
        kernel="G_SIA_cascade",
        meta={"process": "cascade", "note": "SIA cascade injection G_n"}))
    rag.add_edge(Edge(
        EdgeClass.SOURCE, "cascade_VAC_source", vac_bulk,
        kernel="G_VAC_cascade",
        meta={"process": "cascade", "note": "vacancy cascade injection G_m"}))

    # ── Sessile ⟨100⟩ population edge families (loop-conversion work) ─────────
    # ⟨100⟩ loops grow to TEM-visible sizes mainly by ABSORBING mobile ½⟨111⟩
    # clusters (Marian growth): ⟨100⟩_m + ½⟨111⟩_n -> ⟨100⟩_{m+n}.  Cross-
    # population COALESCENCE with the ½⟨111⟩ partner; product = bulk-100 = the
    # primary population, so no product_population redirect is needed.
    rag.add_edge(Edge(
        EdgeClass.COALESCENCE, "SIA100_absorb", sia_100,
        kernel="K_100_absorb", partner_population=sia_111,
        meta={"mechanism": "Marian absorption growth",
              "note": "<100>_m + ½<111>_n -> <100>_{m+n}"}))
    # ⟨100⟩ point-defect ladders (sessile: monomer-driven capture/emission).
    rag.add_edge(Edge(
        EdgeClass.GROWTH, "P3_100_growth", sia_100,
        kernel="K_100_grow",
        meta={"process": "P3", "note": "<100>_n + I_1 -> <100>_{n+1}"}))
    rag.add_edge(Edge(
        EdgeClass.SHRINKAGE, "P3_100_shrink", sia_100,
        kernel="K_100_shrink",
        meta={"process": "P3", "note": "<100>_n + V_1 -> <100>_{n-1}"}))
    rag.add_edge(Edge(
        EdgeClass.DISSOCIATION, "P5_100_emission", sia_100,
        kernel="eps_100_emit",
        meta={"process": "P5", "note": "<100>_n -> <100>_{n-1} + I_1"}))
    rag.add_edge(Edge(
        EdgeClass.SINK, "P4_100_sink", sia_100,
        kernel="D_100_sink",
        meta={"process": "P4", "note": "sessile <100>: D=0 (no fixed-sink loss)"}))
    # Vacancy-cluster annihilation of sessile ⟨100⟩ loops (primary = vacancy).
    rag.add_edge(Edge(
        EdgeClass.ANNIHILATION, "VI_annihilation_100", vac_bulk,
        kernel="K_vi_annih_100", partner_population=sia_100,
        meta={"note": "V_m + <100>_n -> survivor of size |m-n|"}))

    # ── P7 / P8 — gas-pressure-driven and solute-detrapping edges ────────────
    # These two edges complete the EUROFER-97 RAG *structurally* per the
    # paper's master equations 83-86 (Table 1, rows P7 and P8).  They are
    # declared with their kernels registered; the C++ kernel arithmetic that
    # evaluates P7/P8 is the remaining physics-implementation item — the
    # current C++ solver does NOT yet evaluate P7 or P8 (confirmed in the
    # Stage-2 review).  Declaring the edges now makes the RAG the complete,
    # paper-faithful structural description of the EUROFER-97 network; the
    # graph walker / C++ bridge will pick up P7/P8 once their arithmetic
    # lands, with no further change to this declaration.

    # P7 — trap mutation: V_m(He) -> V_{m+1} + I_1  (paper Eq. 75, Table 1).
    # An over-pressurized gas-bearing cavity emits an SIA monomer while
    # GROWING its vacancy content by one (m -> m+1).  The reaction conserves
    # signed lattice-defect content because -m = -(m+1) + 1 (paper Eq. 36).
    #
    # EdgeClass choice — GROWTH.  P7 is a *unary, gas-pressure-driven* size
    # increment of the vacancy cluster (size_shift = +1), so structurally it
    # is the vacancy ladder's GROWTH step: the paper's Table 1 lists P7 as a
    # "TM ladder in m" for c_{m,ell} (gain at m+1, loss at m), exactly the
    # GROWTH stoichiometry.  It is NOT DISSOCIATION (that is size_shift = -1)
    # and NOT SOLUTE_TRAPPING (that fixes size and only shifts the He
    # composition vector).  The dual character of P7 — it ALSO sources an
    # SIA monomer at n = 1 — has no single-class abstraction in the ten-class
    # catalogue; it is recorded in meta['sia_source'] so the graph walker /
    # C++ kernel can add the n = 1 SIA source alongside the vacancy-ladder
    # growth.  (A future refinement could split P7 into a GROWTH edge on the
    # vacancy ladder plus a SOURCE edge on the SIA monomer, both sharing the
    # Gamma_TM kernel; one GROWTH edge with the annotation is kept here so
    # the edge count matches the eight paper process classes P1-P8.)
    rag.add_edge(Edge(
        EdgeClass.GROWTH, "P7_trap_mutation", vac_bulk,
        kernel="Gamma_TM",
        meta={"process": "P7",
              "note": "trap mutation V_m(He) -> V_{m+1} + I_1, paper Eq. 75",
              "sia_source": "emits one I_1 monomer per mutation (n=1 SIA "
                            "source); gas-pressure-driven, kernel Gamma_TM",
              "cpp_status": "kernel registered; C++ arithmetic not yet "
                            "implemented (Stage-2 review)"}))

    # P8 — radiation re-solution: He_ell V_m -> He_{ell-1} V_m + He
    # (paper Eq. 75, Section 5.3, Table 1).  A displacement event athermally
    # ejects one trapped helium atom from a cavity back into the free-He
    # reservoir; the vacancy size m is unchanged and only the helium
    # composition vector decreases by one (ell -> ell-1).
    #
    # EdgeClass choice — SOLUTE_TRAPPING in the detrap direction.  This is
    # exactly the SOLUTE_TRAPPING contract: size fixed, changes_composition
    # True, the unary (detrap) branch of the trap/detrap pair.  gas_species
    # is the helium index (0) in the RAG's single-entry gas list.
    rag.add_edge(Edge(
        EdgeClass.SOLUTE_TRAPPING, "P8_radiation_resolution", vac_bulk,
        kernel="Gamma_res", gas_species=rag.gas_species.index("He"),
        meta={"process": "P8", "direction": "detrap",
              "note": "radiation re-solution He_ell V_m -> He_{ell-1} V_m "
                      "+ He, paper Eq. 75 / Section 5.3",
              "cpp_status": "kernel registered; C++ arithmetic not yet "
                            "implemented (Stage-2 review)"}))

    # ── 4. State layout (contiguous ODE-index map) ───────────────────────────
    layout = StateLayout()

    if equations == "discrete":
        # One ODE per size: SIA ladder 1..I, vacancy ladder 1..V.  Each
        # discrete block carries meta['population'] so the GraphWalker can
        # map (polarity, name) -> block.
        layout.add_discrete("SIA", I, population=sia_111)
        layout.add_discrete("SIA100", I, population=sia_100)
        layout.add_discrete("VAC", V, population=vac_bulk)
        n_vac_classes = V
    else:
        # Bin-moment reduction: discrete prefix + logarithmic bins.  The
        # BinMomentReduction sizes each block; meta['population'] and the
        # reduction object are carried so the reduction wrapper can
        # reconstruct / project around a GraphWalker.
        i_discrete = int(input_data.i_discrete)
        v_discrete = int(input_data.v_discrete)
        I_bin = int(input_data.I_bin)
        V_bin = int(input_data.V_bin)
        shape = str(input_data.shape_function)

        sia_red = BinMomentReduction(I, i_discrete, I_bin, shape)
        sia100_red = BinMomentReduction(I, i_discrete, I_bin, shape)
        vac_red = BinMomentReduction(V, v_discrete, V_bin, shape)
        layout.add_bin_moment(
            "SIA", sia_red.n_discrete, sia_red.n_bins,
            sia_red.moments_per_bin,
            population=sia_111, reduction=sia_red, n_max=I)
        # Sessile ⟨100⟩ loops share the SIA size grid/reduction.  The
        # conversion edges couple the two SIA bin-moment systems; the
        # reduction wrapper / C++ handle the reconstruct→transfer→project
        # (the discrete reference walker uses the discrete layout above).
        layout.add_bin_moment(
            "SIA100", sia100_red.n_discrete, sia100_red.n_bins,
            sia100_red.moments_per_bin,
            population=sia_100, reduction=sia100_red, n_max=I)
        layout.add_bin_moment(
            "VAC", vac_red.n_discrete, vac_red.n_bins,
            vac_red.moments_per_bin,
            population=vac_bulk, reduction=vac_red, n_max=V)
        # The helium block tracks one Q per vacancy class (Case 1) or one
        # scalar (Case 2); the relevant class count is the number of
        # vacancy size classes carrying helium = the bin-moment length.
        n_vac_classes = vac_red.length

    # Helium block — sized by the cascade-selected He reduction
    # (fission -> Case 2 decoupled, length 1; fusion -> Case 1 mean-field,
    # length n_vac_classes).  free_he_tracked=True: free helium c_h is an
    # explicit ODE entry appended by the reduction.
    he_mode = HeReductionMode.from_cascade(cascade)
    he_len = he_mode.he_block_length(n_vac_classes, free_he_tracked=True)
    layout.add_aux("He", he_len,
                   reduction=he_mode.value, free_he_tracked=True,
                   n_vac_classes=n_vac_classes)

    # Conservation-accounting aux block — five cumulative integrals used by
    # the post-processing conservation diagnostics (cumulative bias flux
    # Delta J^d, integrated SIA/vacancy/He production, integrated sink
    # losses; paper Section 4.5 / Eqs. delta_FP, delta_He).
    layout.add_aux("conservation", 5,
                   note="cumulative accounting integrals for delta_FP/delta_He")

    layout.freeze()
    return rag, layout
