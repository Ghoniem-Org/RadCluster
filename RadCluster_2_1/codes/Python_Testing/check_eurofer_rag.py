"""
check_eurofer_rag.py — validation harness for the EUROFER-97 RAG declaration.

Constructs the im5vm2 case (I=V=1000, i_mobile=5, v_mobile=2,
full_CD_fission, quasi-steady-state helium, C_floor=1e-25), builds the
EUROFER-97 reaction admissibility graph and state layout via
``build_eurofer_rag``, and asserts:

  * ``rag.validate()`` passes (every edge's kernel is registered);
  * every edge has a registered kernel of the correct shape;
  * the SIA / vacancy monomer populations are registered;
  * the StateLayout's N_eq is sensible (SIA + VAC ladders + He + aux);
  * the bin-moment variant also builds and validates.

Run:  python check_eurofer_rag.py
"""
import sys
from pathlib import Path

# The EUROFER instantiation imports as RadCluster_2_1.py_utils.materials...
ROOT = Path(r"d:/GitHub/RadCluster/RadCluster_2_1")
sys.path.insert(0, str(ROOT))

from py_utils.input_data import InputData
from py_utils.reaction_rates import ReactionRates
from py_utils.materials import build_eurofer_rag


def make_im5vm2_inputs():
    """Build the im5vm2 InputData + ReactionRates the way RadClusterSimulation does.

    The i_mobile / v_mobile / he_kinetics / C_floor knobs are
    RadClusterSimulation constructor kwargs that it injects into the
    InputData object *after* construction; we replicate that here so the
    standalone harness exercises the exact im5vm2 configuration.
    """
    inp = InputData(I=1000, V=1000, physics_option="full_CD_fission")

    # Injections performed by RadClusterSimulation.__init__ for im5vm2.
    inp.diffusion["i_mobile"] = 5
    inp.derived["i_mobile"]   = 5
    inp.diffusion["v_mobile"] = 2
    inp.derived["v_mobile"]   = 2
    inp.reactions["C_floor"]     = 1e-25
    inp.reactions["he_kinetics"] = "quasi_steady_state"

    rr = ReactionRates(inp)
    return inp, rr


def check(equations):
    print(f"\n{'='*70}\n  EUROFER-97 RAG check — equations={equations!r}\n{'='*70}")
    inp, rr = make_im5vm2_inputs()
    rag, layout = build_eurofer_rag(
        inp, rr, equations=equations, cascade="fission")

    # ── rag.validate() must pass ─────────────────────────────────────────────
    rag.validate()
    print("rag.validate() ......................... OK")

    # ── every edge has a registered, correctly-shaped kernel ─────────────────
    assert len(rag.edges) > 0, "RAG has no edges"
    for e in rag.edges:
        k = rag.kernel(e.kernel)            # raises KeyError if unregistered
        import numpy as np
        arr = np.asarray(k, dtype=float)
        assert np.all(np.isfinite(arr)), f"edge {e.label!r}: non-finite kernel"
        # Binary edges (COALESCENCE / ANNIHILATION) need a 2-D kernel.
        if e.edge_class.value in ("coalescence", "annihilation"):
            assert arr.ndim == 2, \
                f"edge {e.label!r}: {e.edge_class} needs a 2-D kernel"
        # Unary/source size-axis edges need scalar or 1-D kernels.
        else:
            assert arr.ndim in (0, 1), \
                f"edge {e.label!r}: expected scalar/1-D kernel, got {arr.ndim}-D"
    print(f"every edge has a registered kernel ..... OK  ({len(rag.edges)} edges)")

    # ── monomer populations are registered for both polarities ───────────────
    from py_utils.core import Polarity
    sia_mono = rag.monomer_population(Polarity.SIA)
    vac_mono = rag.monomer_population(Polarity.VACANCY)
    assert sia_mono.name == "bulk-111" and sia_mono.polarity is Polarity.SIA
    assert vac_mono.name == "bulk" and vac_mono.polarity is Polarity.VACANCY
    print("SIA + vacancy monomer populations set .. OK")

    # ── gas list is helium-only ──────────────────────────────────────────────
    assert rag.gas_species == ["He"], f"gas list = {rag.gas_species}"
    print("gas list = ['He'] ...................... OK")

    # ── StateLayout sanity ───────────────────────────────────────────────────
    assert layout.has("SIA") and layout.has("SIA100") and layout.has("VAC")
    assert layout.has("He") and layout.has("conservation")
    assert layout.block("conservation").length == 5
    sia_len    = layout.block("SIA").length
    sia100_len = layout.block("SIA100").length
    vac_len    = layout.block("VAC").length
    he_len     = layout.block("He").length
    expected = sia_len + sia100_len + vac_len + he_len + 5
    assert layout.N_eq == expected, \
        f"N_eq {layout.N_eq} != {expected}"
    assert layout.N_eq > 0
    if equations == "discrete":
        # discrete im5vm2: SIA(1000) + SIA100(1000) + VAC(1000)
        #                  + He(2: Q_tot + c_h) + aux(5)
        assert sia_len == 1000 and sia100_len == 1000 and vac_len == 1000, \
            f"discrete ladders should be 1000, got {sia_len}/{sia100_len}/{vac_len}"
        assert he_len == 2, f"fission He block should be 2 (Q_tot + c_h), got {he_len}"
    print(f"StateLayout N_eq sensible .............. OK  (N_eq = {layout.N_eq})")

    # ── reports ──────────────────────────────────────────────────────────────
    print("\n--- rag.summary() ---")
    print(rag.summary())
    print("\n--- state_layout.describe() ---")
    print(layout.describe())
    return rag, layout


def main():
    # Primary requested case: discrete (full_CD_fission).
    check("discrete")
    # Also exercise the bin-moment variant to confirm both layout paths build.
    check("bin_moment")
    print(f"\n{'='*70}\n  ALL CHECKS PASSED\n{'='*70}")


if __name__ == "__main__":
    main()
