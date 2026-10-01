"""
core/cluster_identifier.py — Abstract cluster state-space primitives (Layer 1).

Host-independent.  Implements the indexed cluster state space of
Ghoniem (2026), *A Generalized Graph-Based Cluster Dynamics Framework
for Irradiated Materials*, Section 2.1.

A cluster identifier is the four-tuple

    sigma = (chi, n, p, c)

  chi : polarity index in {-1, +1}  (-1 vacancy-type, +1 SIA-type)
  n   : size >= 1, the count of host Frenkel-pair sites of that polarity
        carried by the cluster
  p   : population label (a :class:`Population`)
  c   : composition vector, the counts of trapped non-host species
        (solutes / gases), e.g. (c_He, c_H)

The polarity index classifies a cluster by the sign of its net
Frenkel-pair content with respect to the host lattice.  It is unaffected
by solute trapping: a helium-vacancy cluster He_j V_n retains chi = -1
regardless of j.

Nothing here is host-specific.  A host material declares *which*
populations and *which* (chi, n, p, c) vertices are admissible; this
module only provides the vocabulary.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from enum import IntEnum
from typing import Tuple


class Polarity(IntEnum):
    """Sign of a cluster's net Frenkel-pair content w.r.t. the host lattice.

    The polarity index splits the state space into two sub-spaces without
    a host-specific taxonomy.  It also serves as the *species-class*
    discriminator of the paper (C_V for vacancy-type, C_I for SIA-type).
    """

    VACANCY = -1
    SIA = +1

    @property
    def label(self) -> str:
        return "vacancy" if self is Polarity.VACANCY else "sia"

    @property
    def signed_defect(self) -> int:
        """Signed lattice-defect content per unit size (-1 vacancy, +1 SIA)."""
        return int(self.value)

    @property
    def opposite(self) -> "Polarity":
        return Polarity.SIA if self is Polarity.VACANCY else Polarity.VACANCY


@dataclass(frozen=True, order=True)
class Population:
    """A defect population within one polarity layer (paper Section 2.1/3.1).

    A population is a subcollection of clusters of the same polarity that
    share a mobility law, orientation state, microstructural environment,
    or sink coupling.  Populations are physically distinct *phases* of the
    cluster (different binding energies, mobility, energetic stability),
    not mere spatial coordinates; transport between populations is itself
    a CD reaction (an inter-population edge).

    The host material declares its population set, e.g.
        vacancy:  {bulk, GB, disl, ppt}
        bcc SIA:  {bulk-111, bulk-100, disl-trapped}
    The abstract core never enumerates populations itself.
    """

    name: str
    polarity: Polarity
    n_min: int = 1          # minimum admissible size for this population
    mobile_max: int = 1     # largest size that is mobile (>= 0)

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("Population.name must be non-empty")
        if self.n_min < 1:
            raise ValueError(
                f"Population {self.name!r}: n_min must be >= 1 (got {self.n_min})")
        if self.mobile_max < 0:
            raise ValueError(
                f"Population {self.name!r}: mobile_max must be >= 0 "
                f"(got {self.mobile_max})")

    def admits_size(self, n: int) -> bool:
        """True if size ``n`` is admissible for this population."""
        return n >= self.n_min

    def is_mobile(self, n: int) -> bool:
        """True if a size-``n`` cluster of this population is mobile."""
        return self.n_min <= n <= self.mobile_max


@dataclass(frozen=True)
class ClusterIdentifier:
    """The four-tuple sigma = (chi, n, p, c) identifying one cluster vertex.

    Instances are immutable and hashable, so they can be used directly as
    keys of the reaction admissibility graph's vertex set.
    """

    polarity: Polarity
    size: int
    population: Population
    composition: Tuple[int, ...] = ()

    def __post_init__(self) -> None:
        if self.population.polarity is not self.polarity:
            raise ValueError(
                f"polarity {self.polarity.label!r} inconsistent with "
                f"population {self.population.name!r} "
                f"(polarity {self.population.polarity.label!r})")
        if self.size < 1:
            raise ValueError(f"cluster size must be >= 1 (got {self.size})")
        if any(int(x) < 0 for x in self.composition):
            raise ValueError(
                f"composition counts must be >= 0 (got {self.composition})")

    # ── derived quantities ────────────────────────────────────────────────
    @property
    def is_monomer(self) -> bool:
        """True for a bare point-defect monomer (size 1, no trapped solute)."""
        return self.size == 1 and not any(self.composition)

    @property
    def signed_defect(self) -> int:
        """Signed lattice-defect content q(sigma) = chi * n (paper Eq. 49).

        Used by the host conservation laws.  A trapped gas atom carries
        zero signed defect, so the composition does not enter here.
        """
        return self.polarity.signed_defect * self.size

    @property
    def total_solute(self) -> int:
        """Total number of trapped non-host atoms across all species."""
        return sum(int(x) for x in self.composition)

    # ── constructors that derive a neighbour vertex ───────────────────────
    def with_size(self, n: int) -> "ClusterIdentifier":
        """Return the same identifier at a different size ``n``."""
        return replace(self, size=n)

    def with_composition(self, c) -> "ClusterIdentifier":
        """Return the same identifier with a different composition vector."""
        return replace(self, composition=tuple(int(x) for x in c))

    def with_population(self, p: Population) -> "ClusterIdentifier":
        """Return the same cluster moved to a different population ``p``."""
        if p.polarity is not self.polarity:
            raise ValueError(
                "cannot move a cluster to a population of opposite polarity")
        return replace(self, population=p)

    def __repr__(self) -> str:  # compact, e.g. <sia n=4 bulk-111 c=()>
        return (f"<{self.polarity.label} n={self.size} "
                f"{self.population.name} c={self.composition}>")
