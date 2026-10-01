"""
core/state_layout.py — Mapping RAG vertices to the flat ODE state vector.

Host-independent.  The reaction admissibility graph defines *which*
clusters exist; integration needs them packed into a single contiguous
state vector  y in R^{N_eq}.  :class:`StateLayout` is the bijection

    (block, local index)  <->  global ODE index

It is built from a list of named :class:`Block` objects laid out
contiguously.  This formalises — in one place — the index bookkeeping
that the legacy ``rate_equations.py`` / ``bin_moment_rates.py`` scattered
across ad-hoc attributes (``i_SIA``, ``i_VAC``, ``i_Q``, ``i_He``,
``i_J_*``); the graph walker and the C++ bridge both consume a
StateLayout instead of re-deriving offsets.

A block is one of:

  * ``discrete``    — one ODE per size for a (polarity, population) axis
  * ``bin_moment``  — discrete prefix + P moments per logarithmic bin
  * ``scalar``      — a single auxiliary ODE (e.g. free He, Q_tot)
  * ``aux``         — a fixed-length auxiliary array (e.g. the five
                      cumulative conservation-accounting integrals)

The layout carries no physics and no host parameters; a host's
state-space-reduction choice (full discrete vs. bin-moment, He Case 1 vs.
Case 2) is expressed purely as the *set of blocks* it adds.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

_KINDS = ("discrete", "bin_moment", "scalar", "aux")


@dataclass(frozen=True)
class Block:
    """One contiguous segment of the ODE state vector.

    Attributes
    ----------
    name : str
        Unique identifier of the block within a layout.
    kind : str
        One of ``discrete``, ``bin_moment``, ``scalar``, ``aux``.
    length : int
        Number of ODE entries the block occupies.
    meta : dict
        Block-specific descriptors, e.g. for ``bin_moment``:
        ``n_discrete``, ``n_bins``, ``moments_per_bin``; for a
        ``discrete`` size axis: ``population``, ``n_max``.
    """

    name: str
    kind: str
    length: int
    meta: Dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.kind not in _KINDS:
            raise ValueError(
                f"Block {self.name!r}: kind must be one of {_KINDS} "
                f"(got {self.kind!r})")
        if self.length < 0:
            raise ValueError(f"Block {self.name!r}: length must be >= 0")


class StateLayout:
    """An ordered set of :class:`Block` s laid out contiguously in y."""

    def __init__(self) -> None:
        self._blocks: List[Block] = []
        self._offset: Dict[str, int] = {}
        self._frozen = False

    # ── construction ──────────────────────────────────────────────────────
    def add_block(self, block: Block) -> Block:
        """Append a block; its global offset is the running total length."""
        if self._frozen:
            raise RuntimeError("cannot add blocks to a frozen StateLayout")
        if block.name in self._offset:
            raise ValueError(f"duplicate block name {block.name!r}")
        self._offset[block.name] = self.N_eq
        self._blocks.append(block)
        return block

    def add_discrete(self, name: str, n_max: int, **meta) -> Block:
        """Add a one-ODE-per-size discrete axis covering sizes 1..n_max."""
        return self.add_block(Block(name, "discrete", n_max,
                                    {"n_max": n_max, **meta}))

    def add_bin_moment(self, name: str, n_discrete: int, n_bins: int,
                       moments_per_bin: int, **meta) -> Block:
        """Add a hybrid discrete-prefix + bin-moment axis."""
        length = n_discrete + moments_per_bin * n_bins
        return self.add_block(Block(name, "bin_moment", length,
                                    {"n_discrete": n_discrete,
                                     "n_bins": n_bins,
                                     "moments_per_bin": moments_per_bin,
                                     **meta}))

    def add_scalar(self, name: str, **meta) -> Block:
        """Add a single auxiliary ODE."""
        return self.add_block(Block(name, "scalar", 1, dict(meta)))

    def add_aux(self, name: str, length: int, **meta) -> Block:
        """Add a fixed-length auxiliary array (e.g. accounting integrals)."""
        return self.add_block(Block(name, "aux", length, dict(meta)))

    def freeze(self) -> "StateLayout":
        """Lock the layout; no further blocks may be added."""
        self._frozen = True
        return self

    # ── queries ───────────────────────────────────────────────────────────
    @property
    def N_eq(self) -> int:
        """Total number of ODEs across all blocks."""
        return sum(b.length for b in self._blocks)

    @property
    def blocks(self) -> Tuple[Block, ...]:
        return tuple(self._blocks)

    def has(self, name: str) -> bool:
        return name in self._offset

    def block(self, name: str) -> Block:
        for b in self._blocks:
            if b.name == name:
                return b
        raise KeyError(f"no block {name!r} in StateLayout")

    def offset(self, name: str) -> int:
        """Global ODE index of the first entry of block ``name``."""
        try:
            return self._offset[name]
        except KeyError:
            raise KeyError(f"no block {name!r} in StateLayout") from None

    def slice(self, name: str) -> slice:
        """Global index ``slice`` spanning block ``name``."""
        off = self.offset(name)
        return slice(off, off + self.block(name).length)

    def index(self, name: str, local: int) -> int:
        """Global ODE index of local entry ``local`` within block ``name``.

        For a ``discrete`` size axis, ``local`` is the 0-based size index
        (size n maps to local n-1).
        """
        blk = self.block(name)
        if not (0 <= local < blk.length):
            raise IndexError(
                f"local index {local} out of range for block {name!r} "
                f"(length {blk.length})")
        return self._offset[name] + local

    def describe(self) -> str:
        lines = [f"StateLayout: N_eq = {self.N_eq}"]
        for b in self._blocks:
            off = self._offset[b.name]
            lines.append(f"  [{off:6d}:{off + b.length:6d}] "
                         f"{b.name:24s} {b.kind:10s} len={b.length}")
        return "\n".join(lines)

    def __len__(self) -> int:
        return self.N_eq

    def __repr__(self) -> str:
        return (f"<StateLayout N_eq={self.N_eq} blocks={len(self._blocks)} "
                f"{'frozen' if self._frozen else 'open'}>")
