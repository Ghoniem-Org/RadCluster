# Genealogy Dynamics

**Orphan branch of the [RadCluster](https://github.com/Ghoniem/RadCluster)
repository.** It shares no history with `main` and contains none of the
RadCluster source; checking it out replaces the working tree entirely.

A two-sex cluster-dynamics model of the U.S. population. Instead of defect
clusters binned by size, people are binned by **pedigree depth** in two
populations — paternal and maternal — and the distribution evolves under
immigration (source), cross-population mating and births (transitions), and
death and emigration (sinks).

All of it lives in **[`genealogy-dynamics/`](genealogy-dynamics/)**. Start with
[its README](genealogy-dynamics/README.md), which carries the full RAG mapping,
the master equation and its exact conservation identity, the two bloodline
rules, and the run modes (colonial replay 1650–2025, hindcast 1980–2024,
projection 2025–2050).

## Why this lives in the RadCluster repository

It is the RadCluster method applied outside radiation damage, not a separate
codebase. It integrates the same master equation,

$$\frac{d\mathbf{c}}{dt} = \mathbf{P}(t) + \mathbf{S}\,\mathbf{J}(\mathbf{c}) - \mathbf{D}(t)\,\mathbf{c}$$

and keeps the same Reaction Admissibility Graph structure — a binned state
space, a source, admissible transitions between bins, and sinks. Only the host
declaration changes:

| RadCluster | Genealogy Dynamics |
|---|---|
| cluster of size *n* (vertex) | pedigree depth *d* × population, paternal or maternal |
| cascade production (source) | immigration, split between the two populations |
| growth *n* → *n*+1 | births: mating flux over paternal × maternal pairs, child bin set by the bloodline rule |
| rate kernels | mating weights and assortativity, fertility, mortality, emigration |
| fixed sinks | death and emigration, per cluster state |
| provenance-stamped `output/` | `outputs/<timestamp>_<scenario>/` with `provenance.md` |

Births conserve people exactly: every column of **S** sums to one, and the
total balance holds to $10^{-15}$.

## Working on this branch

```bash
git checkout Genealogy-Dynamics   # working tree becomes genealogy-dynamics/
git checkout main                 # RadCluster restored
```

Switching here deletes RadCluster's tracked files from the working tree — they
remain safe in `main` and in the object store — and leaves behind everything
`main` did *not* track: build output, run directories, virtual environments,
several hundred megabytes of it. The root `.gitignore` on this branch hides
those leftovers so `git status` stays readable.

Stage explicitly even so — `git add genealogy-dynamics/` rather than
`git add -A` — so an unrelated working tree can never be swept into this
history.

## Related branches

| Branch | Contents |
|---|---|
| `main` | The RadCluster suite: graph-based cluster dynamics for irradiated materials. |
| `Stock-Market-Dynamics` | The same method applied to S&P 500 market regimes. |
