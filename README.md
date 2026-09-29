# Genealogy Dynamics

**Two-sex graph-based cluster dynamics applied to human population genealogy.**

Genealogy Dynamics is a standalone prototype that adapts the Reaction Admissibility
Graph (RAG) cluster-dynamics method of [RadCluster](https://github.com/Ghoniem/RadCluster)
to a population problem. Instead of defect clusters binned by size, people are binned
by **pedigree depth** in two populations — paternal and maternal — and the distribution
evolves under immigration (source), cross-population mating and births (transitions),
and death and emigration (sinks).

The point of the exercise is that nothing about the solver, the conservation
machinery or the graph walker is specific to radiation damage. Changing the host
declaration — what a vertex is, what a source is, which transitions are admissible —
is enough to move the method to an entirely different domain.

- **Branch** — this is an **orphan branch** of the RadCluster repository. It shares no
  history with `main` and contains none of the RadCluster source; checking it out
  replaces the working tree entirely.
- **Languages** — Python (`numpy`, `scipy`, `pandas`, `matplotlib`).
- **License** — MIT.
- **Method reference** — N. M. Ghoniem, *A Generalized Graph-Based Cluster Dynamics
  Framework for Irradiated Materials*, Journal of Nuclear Materials (submitted, 2026);
  a copy is kept at
  [`genealogy-dynamics/docs/ghoniem_cluster_dynamics_framework.pdf`](genealogy-dynamics/docs/ghoniem_cluster_dynamics_framework.pdf).
- **Project documentation** — the authoritative technical README is
  [`genealogy-dynamics/README.md`](genealogy-dynamics/README.md). This file is the
  branch-level entry point.

---

## Contents

1. [Capabilities](#1-capabilities)
2. [Methodology](#2-methodology)
3. [Repository structure](#3-repository-structure)
4. [Installation](#4-installation)
5. [Quick start](#5-quick-start)
6. [Inputs and data sources](#6-inputs-and-data-sources)
7. [Outputs and provenance](#7-outputs-and-provenance)
8. [Verification](#8-verification)
9. [Limitations](#9-limitations)
10. [Relationship to the other branches](#10-relationship-to-the-other-branches)
11. [License and contact](#11-license-and-contact)

---

## 1. Capabilities

| Capability | Where | Notes |
|---|---|---|
| Two-sex cluster state space | `genealogy/populations.py` | pedigree depth $d = 0 \ldots n$ × population $\alpha \in \{\mathrm{pat}, \mathrm{mat}\}$; $n = 14$ by default (`GENEALOGY_MAX_DEPTH`) |
| Swappable inheritance rule | `genealogy/model.py` | kernels, master equation, solver and plots are written against a swappable child-cluster rule $g(q,r)$, so rule comparisons are apples-to-apples |
| Exact population conservation | `genealogy/model.py` | every column of the stoichiometric matrix sums to 1, so births conserve people; the total balance is verified to $10^{-15}$ |
| Assortative mating | `genealogy/assortativity.py` | mating flux interpolates between fully assortative and fully random by a parameter $\alpha$, with a time-dependent schedule |
| Historical replay | `run_colonial.py` | 1650–2025 with anchored historical rates, producing per-population cluster heatmaps |
| Hindcast against real rates | `run_hindcast.py` | 1980–2024 replay on measured immigration, birth and death series, with emigration calibrated rather than assumed |
| Forward projection | `run_case2.py` | 2025–2050 from the 2025 handoff state on projected inputs — explicitly *not* forecasts |
| Regional variant | `genealogy/regional.py`, `run_*_regional.py` | regional disaggregation with its own verification script |
| Religious-affiliation variant | `genealogy/religion.py`, `run_*_religious.py` | affiliation-resolved populations, with an aggregation consistency check |
| Age structure | `genealogy/age.py`, `plot_age_*.py` | age-resolved extension and its comparison panels |
| Sex-specific rates | `genealogy/sex_rates.py` | mortality and emigration may be uniform or sex-specific |
| Stochastic runs | `run_stochastic_colonial.py`, `aggregate_stochastic.py` | ensemble runs with resume support |
| Provenance-stamped output | all run scripts | every run writes a timestamped directory with `provenance.md` alongside results and figures |

---

## 2. Methodology

### 2.1 The master equation

The state $\mathbf{c}(t) \in \mathbb{R}^{2(n+1)}_{\ge 0}$ holds population counts in
millions, paternal block followed by maternal block, and evolves as

$$\frac{d\mathbf{c}}{dt} = \mathbf{P}(t) + \mathbf{S}\,\mathbf{J}(\mathbf{c}) - \mathbf{D}(t)\,\mathbf{c}$$

which is the same equation RadCluster integrates, with the same three terms: a source,
a stoichiometry-weighted sum of reaction fluxes, and a diagonal loss.

- **Source** $\mathbf{P}(t)$ — immigration $I(t)$, split $\sigma / (1-\sigma)$ into the
  depth-0 bin of each population.
- **Loss** $\mathbf{D}(t)$ — diagonal, mortality plus emigration $\mu_\alpha(t) + \varepsilon_\alpha$.
- **Flux** $J_{qr}(\mathbf{c})$ — births over ordered paternal × maternal pairs, with a
  marriage-function total $B(\mathbf{c}) = \min\left(\sum_d \beta^{\rm pat}_d c_{d,{\rm pat}},\ \sum_d \beta^{\rm mat}_d c_{d,{\rm mat}}\right)$
  and an assortativity parameter $\alpha$ mixing the assortative and random limits.
- **Stoichiometry** $\mathbf{S}$ — the child lands in bin $g(q,r)$ set by the
  inheritance rule, split $s / (1-s)$ between the populations with $s = 0.512$, the
  measured NCHS sex ratio at birth.

The exact identity
$\frac{d}{dt}\sum_{d,\alpha} c_{d,\alpha} = I + B - \sum_{d,\alpha}(\mu_\alpha + \varepsilon_\alpha)\,c_{d,\alpha}$
is checked numerically on every run.

### 2.2 The RAG mapping

| Cluster dynamics (RadCluster) | Genealogy dynamics |
|---|---|
| cluster of size *n* (vertex) | cluster state $(d, \alpha)$: pedigree depth × population; $C_0$ is foreign-born |
| cascade production (source) | immigration, split between the two populations |
| monomer absorption, growth $n \to n+1$ | births: mating flux over paternal × maternal pairs, child bin from the rule $g(q,r)$ |
| rate kernels | mating weights and assortativity, fertility, mortality, emigration |
| fixed sinks | death and emigration, per cluster state |
| `input_parameters.xlsx` | projected and measured input series under `data/` |
| provenance-stamped `output/` | `outputs/<timestamp>_<scenario>/` with `provenance.md` |

### 2.3 Inheritance rules

Two child-cluster rules are implemented — *shallowest* and *deepest* inheritance —
and because everything downstream is written against the rule rather than around it,
the two can be compared on identical inputs. See the project README for the rule
definitions and the colonial-run comparison.

---

## 3. Repository structure

```
/
├── README.md                     # this file — branch entry point
├── .gitignore                    # root-anchored; see Section 10
└── genealogy-dynamics/
    ├── README.md                 # authoritative technical documentation
    ├── requirements.txt
    ├── genealogy/                # the package
    │   ├── model.py              # master equation, stoichiometry, inheritance rules
    │   ├── kernels.py            # fertility, mortality, emigration, mating weights
    │   ├── populations.py        # the two-sex cluster state space
    │   ├── assortativity.py      # assortative ↔ random mating mixing
    │   ├── solver.py             # time integration
    │   ├── age.py, regional.py, religion.py, sex_rates.py   # extensions
    │   └── visualization.py
    ├── data/                     # build scripts, measured series, provenance notes
    ├── docs/                     # STEP_LOG.md and the method paper
    ├── run_*.py                  # the run modes (Section 5)
    ├── plot_*.py                 # figure generators
    ├── test_*.py, verify_*.py, check_*.py   # verification (Section 8)
    └── outputs/                  # timestamped run directories (Section 7)
```

---

## 4. Installation

```bash
git clone https://github.com/Ghoniem/RadCluster.git
cd RadCluster
git checkout Genealogy-Dynamics
cd genealogy-dynamics
python -m pip install -r requirements.txt
```

Requirements are Python 3.9 or newer with `numpy`, `scipy`, `pandas` and
`matplotlib`. There is no compiled component: unlike `main`, this branch needs no
C++ toolchain, no SUNDIALS and no Excel reader.

---

## 5. Quick start

```bash
python3 run_case2.py              # 2025–2050 projection, two-sex × both rules
python3 run_case2_religious.py    # same, religious-affiliation model

python3 data/fetch_data.py && python3 data/build_inputs.py   # rebuild measured inputs
python3 data/calibrate_emigration.py                          # grid-search ε on 1980–2024 totals
python3 run_hindcast.py                                       # hindcast 1980–2024 × both rules

python3 run_colonial.py           # 1650–2025 × both rules, α sweep
python3 plot_timeseries.py        # total-share time series from the colonial run
python3 plot_distributions.py     # final cluster distributions
python3 test_symmetry.py          # two-sex vs one-sex reference
```

The three modes answer different questions. `run_colonial.py` replays 1650–2025 on
anchored historical rates. `run_hindcast.py` replays 1980–2024 on measured series and
validates against observed totals. `run_case2.py` projects forward from the 2025 state
on *assumed* inputs — projections, not forecasts.

---

## 6. Inputs and data sources

All series are public and require no credentials. The project README carries the full
thirteen-row source table with years and citations; in outline:

| Class | Series | Source |
|---|---|---|
| Immigration | lawful permanent residents; unauthorized immigrant stock | DHS Office of Homeland Security Statistics |
| Vital rates | crude birth rate, crude death rate, total fertility rate | World Bank WDI |
| Population | total population; international migrant stock | World Bank WDI |
| Historical | colonial estimates 1650–1780; decennial censuses 1790–2020 | U.S. Census Bureau |
| Pre-1960 anchors | white crude birth rates; crude death rate anchors | Haines; Blodget; Death Registration Area |
| Sex ratio at birth | $s = 0.512$ | NCHS/CDC National Vital Statistics Reports |

Rebuild with `python3 data/fetch_data.py && python3 data/build_inputs.py`. Raw
responses are cached under `data/raw/`; cleaned tables carry `quality` and
`quality_note` columns so assumed values are never silently mixed with measured ones.

---

## 7. Outputs and provenance

Each run writes a timestamped directory:

```
outputs/<YYYYMMDD_HHMMSS>_<scenario>/
├── provenance.md      # inputs, parameters, rule, α, git state
├── results.csv        # yearly per-population clusters, totals, paternal/maternal split
└── *.png              # figures
```

Rule comparisons additionally write
`outputs/<YYYYMMDD_HHMMSS>_comparison/` with the deepest-cluster-share overlays and a
`summary.md`. The convention matches `main`'s: a run is not a result until it carries
the provenance record that lets someone else reproduce it.

---

## 8. Verification

Four checks run automatically with the model:

1. **Non-negativity** — the minimum over all cluster states and times is $\ge 0$.
2. **Pure decay** — with immigration and fertility switched off, the numerical
   solution matches $e^{-(\mu+\varepsilon)t}$ to solver tolerance.
3. **Exact population balance** — numerical $d(\text{total})/dt$ equals the analytic
   $I + B - \sum(\mu+\varepsilon)c$ to floating-point precision, not merely to solver
   tolerance.
4. **Symmetry reduction** (`test_symmetry.py`) — with $s = \sigma = 0.5$ the two-sex
   cluster totals reproduce a one-sex reference to better than $10^{-6}$ M.

The extensions carry their own: `verify_regional.py`, `check_religious_aggregation.py`,
`test_age_consistency.py`, `verify_step5_regression.py`, `validate_step8.py`.

---

## 9. Limitations

This is a prototype, and the project README keeps the authoritative list. The
standing caveats are that the base model is well-mixed and national, that the
projection mode runs on assumed rather than forecast inputs, and that assortativity,
the 1980 cluster split and the immigrant sex split remain assumed even in the
hindcast, where emigration alone is calibrated. Read
[`genealogy-dynamics/README.md`](genealogy-dynamics/README.md) before quoting any
number from a run.

---

## 10. Relationship to the other branches

| Branch | Contents |
|---|---|
| `main` | The RadCluster suite — graph-based cluster dynamics for irradiated materials. |
| `Genealogy-Dynamics` | This branch. |
| `Stock-Market-Dynamics` | The same method applied to S&P 500 market regimes. |

```bash
git checkout Genealogy-Dynamics   # working tree becomes genealogy-dynamics/
git checkout main                 # RadCluster restored
```

Switching here deletes RadCluster's tracked files from the working tree — they remain
safe in `main` and in the object store — and leaves behind everything `main` did *not*
track: build output, run directories, virtual environments, several hundred megabytes
of it. `main`'s `.gitignore` leaves with it, so the root `.gitignore` on this branch
exists to hide those leftovers and keep `git status` readable.

Stage explicitly even so — `git add genealogy-dynamics/` rather than `git add -A` — so
an unrelated working tree can never be swept into this history.

---

## 11. License and contact

Released under the [MIT License](https://github.com/Ghoniem/RadCluster/blob/main/LICENSE),
© 2026 Nasr M. Ghoniem.

**Nasr M. Ghoniem** — Mechanical and Aerospace Engineering Department,
University of California, Los Angeles.
