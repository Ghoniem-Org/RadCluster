# RadCluster

**Graph-based cluster dynamics for irradiated materials.**

RadCluster is a physics-based simulation suite for the evolution of radiation-induced
defect populations — self-interstitial atom (SIA) dislocation loops, vacancy and
helium–vacancy cavities, free interstitial helium, and (in development) solute
segregation and precipitation — in crystalline solids under neutron or ion
irradiation. The reference host material is the reduced-activation
ferritic–martensitic (RAFM) steel **EUROFER-97**, the European structural candidate
for fusion first-wall and blanket components.

The distinguishing feature of the code is that the cluster-dynamics master equation is
not hand-written per material. Instead, a material is *declared* as a
**Reaction Admissibility Graph** (RAG) whose edges are elementary reactions drawn from
a fixed catalogue of ten abstract classes; the right-hand side of the ODE system is
then assembled by walking that graph. Adding a new host material is a declaration,
not a rewrite of the solver.

- **Current release** — `v2.2.0`, the active module `radcluster_code/`. This
  release changes what a checkout contains, not what the solver computes: the
  five archived modules, the assistant's local session state and a scratch
  output directory are no longer tracked — 226 files and 27 MB fewer in the
  working tree (240 MB → 212 MB). Note that this does not shrink the *download*:
  git still fetches the full history, which retains every one of those files.
  The physics is unchanged from `v2.1.2`, which carried the
  Frenkel-pair conservation correction (§2.9); `v2.1.1` renamed the module. The
  preceding graph-based baseline is tagged `v2.0.0` (see *Archived modules*
  below — they are no longer carried in a clone).
- **Languages** — Python (model definition, orchestration, post-processing) and
  C++17 (production ODE solver, SUNDIALS/CVODE).
- **License** — MIT.
- **Primary reference** — N. M. Ghoniem, *A Generalized Graph-Based Cluster Dynamics
  Framework for Irradiated Materials*, Journal of Nuclear Materials (submitted, 2026);
  [`docs/Formulation/cluster_dynamics_framework/Generalized_Cluster_Dynamics.pdf`](docs/Formulation/cluster_dynamics_framework/Generalized_Cluster_Dynamics.pdf).
  Section and equation numbers throughout the code and this README refer to that document.

---

## Contents

1. [Capabilities](#1-capabilities)
2. [Methodology](#2-methodology)
3. [Repository structure](#3-repository-structure)
4. [Installation](#4-installation) — [container image](#40-container-image-no-build-required)
5. [Quick start](#5-quick-start)
6. [Inputs](#6-inputs)
7. [Outputs and provenance](#7-outputs-and-provenance)
8. [Examples, tests and verification](#8-examples-tests-and-verification)
9. [Computational aspects](#9-computational-aspects)
10. [Documentation](#10-documentation)
11. [Publications and how to cite](#11-publications-and-how-to-cite)
12. [Repository conventions](#12-repository-conventions)
13. [License and contact](#13-license-and-contact)

---

## 1. Capabilities

| Capability | Where formulated | Notes |
|---|---|---|
| Graph-based RHS assembly | Section 3.3 | the master equation is a host-independent graph-walking accumulator |
| Reaction Admissibility Graph | Section 2.2 | ten elementary edge classes; physically forbidden reactions excluded at graph-construction time |
| Two equation formulations | Sections 4, 5.4 | *discrete* per-size ODEs, or a logarithmic *bin-moment* grouping |
| Two cascade source models | Sections 5.2, 5.3 | *fission* (He Case 2, decoupled scalar inventory) or *fusion* (He Case 1, mean-field loading per void class) |
| Three intra-bin closures | Section 5.5 | piecewise-constant ($P{=}1$), hat-function ($P{=}2$), log-normal ($P{=}3$) |
| Solute trapping in EUROFER-97 | Eq. 125 | effective jump frequencies with Cr, W, Mn, C, N traps |
| Mixed 1D/3D transport | Section 4.2 | glissile SIA clusters with a mean-free-path correction |
| ½⟨111⟩ → ⟨100⟩ loop-character split | `docs/Formulation/dislocation_loops/` | two coupled SIA loop populations; Marian junction/absorption (kinetic) and Dudarev unary (thermodynamic) channels |
| He–vacancy state-space reduction | Section 8 | bubble loading eliminated as a mean field or as a decoupled inventory |
| Helium equation of state | Table 8 | third-order virial EOS for bubble pressure and vacancy/He binding |
| C++ SUNDIALS/CVODE BDF solver | Section 6.2 | dense, banded, GMRES, or KLU sparse-direct linear solver |
| Woodbury preconditioner | Section 6.2 | bordered-banded Sherman–Morrison–Woodbury preconditioner for GMRES |
| Conservation diagnostics | Eqs. 96, 97 | Frenkel-pair ($\delta_{\rm FP}$) and helium ($\delta_{\rm He}$) residuals reported every run |
| Adaptive size domain | — | the tracked cluster-size domain is expanded on demand when the distribution reaches the boundary |
| Provenance-stamped output | — | every run writes a timestamped directory with git SHA, parameters, solution, summary table and figures |
| Experimental microstructure database | `radcluster_code/Eurofer_micro_database/` | digitized loop/cavity densities and sizes for neutron- and ion-irradiated F/M steels, with fitting and plotting notebooks |

**Under development in `radcluster_code`** — a loop → network-dislocation loss channel
with a dynamically evolving network density (so SIA loop number density saturates with
dose), and radiation-induced segregation (RIS) with solute-rich precipitation
(Cr-rich α′ and Mn–Ni–Si nanofeatures). Both are declared as new populations, edges and
kernels in the EUROFER-97 host graph and reuse the existing solver, reductions and
conservation machinery unchanged. See
[`docs/Formulation/dislocation_loops/loop_network_loss.tex`](docs/Formulation/dislocation_loops/loop_network_loss.tex)
and
[`docs/Formulation/segregation_and_precipitation/radcluster_2_1_RIS_plan.tex`](docs/Formulation/segregation_and_precipitation/radcluster_2_1_RIS_plan.tex).

---

## 2. Methodology

### 2.1 The master equation

For every tracked cluster species $a$ the code integrates

$$\frac{dc_a}{dt} = G_a + \sum_r S_{ar}\,J_r(c) - \mathcal{D}_a c_a ,$$

where $G_a$ is the cascade source, $J_r$ the flux of elementary reaction $r$, $S_{ar}$
its stoichiometric coefficient for species $a$, and $\mathcal{D}_a$ the loss rate at
fixed sinks (network dislocations, grain boundaries, precipitates). Concentrations are
atomic fractions; sizes are atom counts; energies are in eV and temperature in K.

### 2.2 The Reaction Admissibility Graph

A cluster is identified by the four-tuple $\sigma = (\chi, n, p, c)$ — polarity
$\chi \in \\{-1,+1\\}$ (vacancy- or interstitial-type), size $n$, population $p$
(e.g. glissile ½⟨111⟩ loops, sessile ⟨100⟩ loops, bulk cavities), and trapped-solute
composition $c$. The model is the labelled multidigraph

$$G = (V, E, \nu, k),$$

with $V$ the admissible cluster vertices, $E$ every admissible elementary reaction,
$\nu$ the stoichiometric action of each edge, and $k$ the rate-kernel library. Each
edge belongs to exactly one of **ten abstract classes**:

| # | Edge class | Action |
|---|---|---|
| 1 | `GROWTH` | same-polarity monomer absorption, $n \to n+1$ |
| 2 | `SHRINKAGE` | opposite-polarity monomer impingement, $n \to n-1$ |
| 3 | `DISSOCIATION` | thermal monomer emission, $n \to n-1$ |
| 4 | `RECOMBINATION` | cross-polarity monomer annihilation |
| 5 | `ANNIHILATION` | cross-polarity binary (cluster–cluster) hyperedge |
| 6 | `INTER_POPULATION` | $p \to p'$ at fixed size and composition (e.g. loop-character conversion) |
| 7 | `SOLUTE_TRAPPING` | trapping / detrapping, $c \to c \pm e_s$ |
| 8 | `COALESCENCE` | same-polarity binary growth |
| 9 | `SOURCE` | cascade injection, no precursor vertex |
| 10 | `SINK` | absorption at a fixed sink |

Each class carries a fixed *stoichiometric contract* — number of precursor vertices,
kinetic order, whether it crosses polarity layers, how it shifts size, population and
composition. Because the sign pattern of $S_{ar}$ is structural, conservation of signed
defect content $q = \chi n$ holds by construction for every intra-host reaction, rather
than being re-derived for each material.

### 2.3 Two-layer architecture

| Layer | Location | Contents |
|---|---|---|
| **1 — abstract core** (host-independent) | `py_utils/core/`, `cpp_utils/core/` | cluster identifier, the ten edge classes, the RAG container and admissibility logic, the graph walker, the vertex → ODE-index `StateLayout`, and the state-space reductions |
| **2 — material instantiation** | `py_utils/materials/eurofer97/`, `cpp_utils/materials/eurofer97/` | declares the EUROFER-97 populations, its admissible vertex and edge sets, and its P1–P8 rate kernels |

The core carries no bcc-Fe or EUROFER assumptions. A host material *declares* its
populations, adds `Edge` families, and registers size-resolved rate kernels as
precomputed arrays; the abstract core never re-derives a rate. Adding a material means
adding a `materials/<name>/` directory in both `py_utils` and `cpp_utils` — the solver,
preconditioner and reductions are inherited unchanged.

The Python graph walker is the *executable specification* that the production C++ kernel
must agree with; it is also the substrate for the bin-moment
"reconstruct → walk → project" reduction. Every `RadClusterSimulation` object carries
the live graph description on `sim.material_rag` and `sim.state_layout`, so inspecting a
run shows exactly which reactions are admissible and which approximations are active.

### 2.4 Cascade source

The survival-corrected production rate is $G = \eta\,G_{\rm NRT}$, and in-cascade
cluster production follows a power law in cluster size,

$$\epsilon_m^{(i)} = C_i\, m^{-s_i}, \qquad
  C_i = \frac{f_i^{\rm cl}}{\sum_{m=2}^{m_1} m^{1-s_i}},$$

with the monomer channel taking the remainder of the defect content. Fission and fusion
spectra differ in survival efficiency $\eta$, clustered fractions $f^{\rm cl}$, slopes
$s$, and maximum in-cascade cluster sizes (Table 2).

### 2.5 Transport and solute trapping

Point-defect diffusivities follow $D_\alpha = a^2 \nu_\alpha \exp(-E_m^\alpha / k_B T)$.
In an alloy, dissolved solutes retard migration through an effective jump frequency

$$\omega_\alpha^{\rm eff} =
  \frac{\nu_\alpha \exp(-E_m^\alpha/k_B T)}
       {1 + \sum_s z_s c_s \exp\!\left(E_b^{\alpha\text{-}s}/k_B T\right)} ,$$

applied to SIAs, vacancies and helium. Glissile SIA clusters migrate
one-dimensionally with $D_n^{\rm 1D} \propto n^{-s}\exp(-E_m^{\rm 1D}/k_B T)$ and a
mean-free-path correction for rotation and trapping; mobility cut-offs
($i_{\rm mobile}$, $v_{\rm mobile}$) set which cluster sizes participate in
coalescence and 1D reactions.

### 2.6 Reaction kernels

The kernel library (processes P1–P8) covers vacancy–SIA recombination, point-defect
absorption by spherical cavities and by dislocation loops (with bias factors),
loss at fixed sinks, thermal emission, cluster–cluster coalescence and annihilation,
and helium trap mutation and re-solution. Geometric prefactors are exact:
$A_{\rm sph} = (48\pi^2)^{1/3}$, $A_{\rm loop} = 8\sqrt{\pi/\sqrt3}$,
$A_{\rm 1D} = 9/8\pi^{2/3}$, $B_{\rm rot} = (4/\pi)(8\pi/3)^{1/3}$.
Detailed balance between growth and emission kernels is documented in
[`docs/Formulation/reaction_kernels/`](docs/Formulation/reaction_kernels).

### 2.7 Energetics

Vacancy binding to cavities uses the capillarity form
$E_b^v(m) = E_f^v - A_{\rm void}[m^{2/3}-(m-1)^{2/3}]$, corrected for gas pressure from
a third-order virial helium EOS, so helium stabilizes bubbles against thermal emission.
Interstitial loop binding is fitted at small size and blended to the continuum limit,
$E_b^{\rm loop}(n) = A\,n^{+B}$, separately for the ½⟨111⟩ and ⟨100⟩ characters.
Helium binding to bubbles combines the solution energy, the gas-pressure term, and an
atomistically fitted correction.

### 2.8 State-space reductions

Two independent reductions keep the system tractable at engineering doses.

**Helium–vacancy loading.** In the *fusion* case (Case 1) the helium content of a
cavity class equilibrates quickly and is replaced by its mean field
$\bar{\ell}(v)$, giving $I + 2V + 1$ equations. In the *fission* case (Case 2) the
loading decouples into a scalar inventory, giving $I + V + 2$ equations. The unreduced
two-dimensional $(v,\ell)$ grid remains available.

**Logarithmic bin moments.** Sizes up to $i_{\rm discrete}$ (resp. $v_{\rm discrete}$)
are tracked one ODE per size; larger sizes are grouped into logarithmic bins carrying
$P$ moments each,

$$\mu_k^{(0)} = \sum_{n\in\mathcal{B}_k} c_n, \quad
  \mu_k^{(1)} = \sum_{n\in\mathcal{B}_k} n\,c_n, \quad
  \mu_k^{(2)} = \sum_{n\in\mathcal{B}_k} n^2 c_n ,$$

closed by a piecewise-constant, hat-function or log-normal intra-bin shape with
truncation error $O((r-1)^{P+1})$ in the bin ratio $r$. The total system size is
$N_{\rm eq} = i_{\rm discrete} + P\,I_{\rm bin} + v_{\rm discrete} + P\,V_{\rm bin} + n_{\rm He}$,
and the discrete formulation is recovered exactly when the bin count is zero.

### 2.9 Post-processing and verification diagnostics

Swelling is reported through the exact identity $S(t) = S_I(t) + \Delta J^d(t)$
(vacancy content equals SIA content plus the cumulative net bias flux to dislocations),
which yields the dimensionless Frenkel-pair residual $\delta_{\rm FP}(t)$; the helium
balance yields $\delta_{\rm He}(t)$. Both should stay below $\sim10^{-6}$ — values
above $10^{-3}$ indicate a coding error, and the test scripts in
`codes/Python_Testing/` gate on them.

**Frenkel-pair correction (v2.1.2).** Before v2.1.2, the regression case `im5vm2` had
$\delta_{\rm FP} = 1.24\times10^{-2}$ at $t = 10^{-2}$ s, above the error level.

*Cause.* In all three C++ right-hand sides (fission, fusion and bin-moment), SIA thermal
emission $I_{n+1}\to I_n + I_1$ moved the emitting cluster down one size, but never added
the emitted monomer to $I_1$. Vacancy emission does this through `emit_mono`. One SIA was
therefore lost per emission event:

- the vacancy arm of the residual closed to $10^{-12}$;
- the time integral of the emitted-monomer flux equals the SIA arm's deficit at every
  output time.

*Fix.* The monomer is now returned. On `im5vm2`, $\delta_{\rm FP}$ falls to
$2.9\times10^{-7}$, $\delta_{\rm He}$ is $2\times10^{-14}$, and $C_{\rm SIA}^{\rm tot}$ rises
by 0.5 %; the vacancy and helium observables are unchanged.

*Option.* The switch is `sia_emission_monomer` in the reaction parameters:

- `"returned"` is the default;
- `"dropped"` restores the earlier arithmetic bit for bit, for reproducing results
  obtained before v2.1.2;
- in the solver parameter file it is the key `sia_emission_monomer` (1 or 0; absent
  means 0).

*Remaining residual.* In long loop-coarsening runs, the loop distribution can outgrow the
size axis. Products beyond $I$ then carry content out of the domain, and $\delta_{\rm FP}$
grows. Check it against the axis length; it falls as $I$ grows. The correction was found
in the migration of RadCluster into
[GSD](https://github.com/Ghoniem/GSD) (v2.0.0, Stage 5), whose declared equations cannot lose
content and match this C++ kernel at round-off.

Size distributions are rendered as TEM-comparable per-bin densities $dc/dn$ or $dc/dD$
(stair segments whose integral is the population in the displayed range), rather than
per-size concentrations, which is the faithful display primitive for a coarse-grained
representation.

---

## 3. Repository structure

```
RadCluster/
├── radcluster_code/         # Active development — recommended entry point
│   ├── CLAUDE.md               # Physics and solver reference (equations, tables, options)
│   ├── py_utils/               # Python package
│   │   ├── core/                   # Layer 1 — abstract, host-independent
│   │   │   ├── cluster_identifier.py, edge_classes.py, rag.py
│   │   │   ├── graph_walker.py, state_layout.py, to_networkx.py
│   │   │   └── reductions/         # bin-moment and He reductions
│   │   ├── materials/eurofer97/    # Layer 2 — EUROFER-97 RAG declaration
│   │   ├── input_data.py           # workbook reader and derived quantities
│   │   ├── reaction_rates.py, rate_equations.py, bin_moment_rates.py
│   │   ├── binding_energies.py, defect_production.py, loop_energetics.py
│   │   ├── simulation.py           # RadClusterSimulation orchestrator
│   │   ├── cpp_bridge.py           # subprocess bridge to the C++ solver
│   │   ├── post_process.py, size_distributions.py, visualization.py
│   ├── cpp_utils/              # C++17 solver
│   │   ├── core/                   # solver.cpp (CVODE driver), rhs_dispatch.cpp,
│   │   │                           # sparse_jacobian.cpp, parameters.h
│   │   ├── materials/eurofer97/    # rate_kernels.cpp (P1–P8 arithmetic)
│   │   └── CMakeLists.txt
│   ├── codes/
│   │   ├── Notebooks/              # RadCluster_2_1.ipynb — the main driver
│   │   ├── Python_Testing/         # conservation gates, convergence and sweep studies
│   │   └── make_*_figures.py       # RAG, dose and size-effect figure generators
│   ├── digital_twin/           # calibration and verification campaign (see §8)
│   ├── Eurofer_micro_database/ # experimental microstructure database and figures
│   ├── input/                  # input_parameters.xlsx and builders
│   ├── build/                  # CMake artifacts (gitignored)
│   └── output/                 # timestamped run directories (gitignored)
├── archive/                # Earlier modules — gitignored, NOT in a clone (recover from history)
├── docs/
│   ├── Formulation/            # derivations, paper sections, study reports (by topic)
│   ├── Database/               # experimental microstructure databases (xlsx)
│   ├── Literature/             # ~60 reference papers (PDF)
│   └── design_notes/
├── scripts/                # repository utilities (h2jupynb: Hoffman2 Jupyter launcher)
├── requirements.txt
└── LICENSE
```

### Module status

| Module | Status | Description |
|---|---|---|
| `radcluster_code/` | **Active (development)** | Graph-based cluster dynamics plus dislocation-network evolution and RIS/precipitation. Recommended starting point for new users. |
| `RadCluster_2_0` | Archived | The published two-layer graph framework — abstract core plus EUROFER-97 declaration. Released as `v2.0.0`; the regression reference for 2_1. |
| `RadCluster_1_0` | Archived | Earlier non-graph generalized cluster dynamics; superseded by 2_0. |
| `Monomer_CD` | Archived | Monomer-mobility cluster-dynamics scaling reference (Ghoniem & Cho, 1979; no helium). |
| `Zr_RadCluster_1_0` | Archived | Zirconium variant of the 1_0 formulation. |
| `Eurofer`, `Eurofer_CD` | Archived | Earlier EUROFER microstructure and cluster-dynamics notebooks. |

#### Archived modules are in the history, not in a clone

The five archived modules receive no further development, and `archive/` is
gitignored, so **a fresh clone does not contain them**. They are 276 MB of
superseded code that nothing in the active module imports, and carrying it in
every clone served no purpose.

Nothing was deleted. They were untracked, not removed, so published results
stay reproducible — recover any of them from the history:

```bash
git checkout e9f92c3 -- archive/                       # all five, as of the last commit carrying them
git checkout e9f92c3 -- archive/RadCluster_2_0/        # or just one
```

`e9f92c3` is the right ref for this because it is the last commit in which the
modules sat under `archive/`. The tags `v2.0.0` and `eurofer_cd-final` predate
the moves into `archive/`, so at those refs the same modules are at the
repository root — `git checkout v2.0.0 -- RadCluster_2_0/` and
`git checkout eurofer_cd-final -- Eurofer_CD/`, which restore to the root path,
not into `archive/`.

### Branches

| Branch | Kind | Contents |
|---|---|---|
| `main` | production | The RadCluster suite — module, documentation, verification campaign. |
| `campaign-verification` | orphan | `claims/` and `results/` — the coordination store the distributed verification campaign writes to; read through `radcluster_code/digital_twin/verification/.sync/`. |
| `CodeDevelopment` | dormant | Fully merged into `main`; retained for reference, not developed. |

`campaign-verification` is an **orphan branch**: it shares no history with `main`
and holds only machine-written campaign state, never source. It is consumed
through a git worktree rather than checked out directly, so that a pull can never
move a file underneath a running solver.

#### Retired side branches: the method applied outside radiation damage

Two further projects were developed on orphan branches of this repository and
have since been retired. They are preserved as tags rather than branches, on the
remote as well as locally:

| Tag | Project |
|---|---|
| `gsd-freeze-genealogy` | Genealogy Dynamics — a two-sex cluster-dynamics model of the U.S. population |
| `gsd-freeze-market` | Stock-Market Dynamics — S&P 500 market-regime cluster dynamics |

Both were ports of this repository's abstract layers, not separate codebases.
Each integrated the same master equation of Section 2.1 and kept the same
Reaction Admissibility Graph structure — a binned state space, a source,
admissible transitions between bins, and sinks. Only the host declaration
changed, which is the point worth recording:

| | Vertices (state) | Source | Transitions | Sinks |
|---|---|---|---|---|
| **RadCluster** | defect clusters binned by size | displacement cascade | growth, shrinkage, coalescence, character conversion | network dislocations, grain boundaries, precipitates |
| **Genealogy** | people binned by pedigree depth, in paternal and maternal populations | immigration | cross-population mating and births | death, emigration |
| **Stock market** | regime bins — momentum × trailing volatility | index entry | measured Markov drift, bilinear herding flux | index exit |

To bring one back as a working branch:

```bash
git checkout -b Genealogy-Dynamics gsd-freeze-genealogy
```

Be aware of what that does to the working tree. Checking out an orphan branch
deletes RadCluster's tracked files — they remain safely in `main` and in the
object store — and leaves behind everything `main` did *not* track: build
output, run directories, virtual environments. `main`'s `.gitignore` leaves with
it, so those leftovers are no longer hidden and can amount to hundreds of
megabytes of apparently untracked files. Both frozen branches carry a
root-anchored `.gitignore` of their own for exactly this reason, and the
practical rule while on one is to **stage explicitly** — `git add
genealogy-dynamics/` — rather than `git add -A`, which would otherwise sweep an
entire unrelated working tree into the wrong history.

---

## 4. Installation

### 4.0 Container image (no build required)

A published image carries the compiled solver, the Python layer and the input
workbook, so nothing has to be built locally:

```bash
docker pull ghcr.io/ghoniem/radcluster:v2.2.0
docker run --rm -it ghcr.io/ghoniem/radcluster:v2.2.0
```

Inside, `solver` is on `PATH` and `/work` holds `py_utils/`, `input/` and
`codes/` with `PYTHONPATH` already set. Output is written inside the container
and is lost when it exits, so mount a host directory for anything you want to
keep:

```bash
docker run --rm -v "$PWD/output:/work/output" ghcr.io/ghoniem/radcluster:v2.2.0     python3 codes/Python_Testing/check_eurofer_rag.py
```

The solver in the image is compiled for `x86-64-v2` rather than the building
machine's own architecture, so it runs on any x86-64 CPU from about 2009
onwards. That costs a little speed against a local `-march=native` build: for
sustained production runs, build from source as below. The image is rebuilt and
published by `.github/workflows/release.yml` whenever a `v*` tag is pushed.

### 4.1 Python environment

```bash
git clone https://github.com/Ghoniem/RadCluster.git
cd RadCluster
python -m pip install -r requirements.txt
python -m ipykernel install --user --name radcluster --display-name "RadCluster"   # optional
```

Requirements: Python 3.9 or newer with `numpy`, `scipy`, `pandas`, `matplotlib`,
`openpyxl`, `jupyter` and `ipykernel`. `networkx` is optional and used only by the
graph-diagnostics export (`py_utils/core/to_networkx.py`); the solver does not need it.

The Python layer alone is sufficient to build a material graph, inspect the model, and
run the reference walker. Production runs use the C++ solver.

### 4.2 C++ solver

Required:

- a C++17 compiler (GCC, Clang, or MSVC) and CMake ≥ 3.15;
- **SUNDIALS 7.1.1** (CVODE). The build expects it in a `Libraries/` directory that is a
  sibling of the repository root — `Libraries/sundials-7.1.1/` with per-platform
  install trees (`mac_install/`, `win_install/`) — or anywhere `find_package(SUNDIALS)`
  can resolve it (a system or vcpkg install works).

Optional, each auto-detected and each degrading gracefully if absent:

| Dependency | Enables | If missing |
|---|---|---|
| LAPACK/BLAS (Accelerate on macOS; OpenBLAS, reference LAPACK, or MKL elsewhere) | Woodbury preconditioner | compiled out; the solver falls back to Jacobi at runtime |
| OpenMP (`brew install libomp` on macOS; bundled with MSVC) | parallel RHS in `active_window` mode | that mode runs single-threaded |
| SuiteSparse/KLU, with SUNDIALS built `-DENABLE_KLU=ON` | `linsol='klu'` sparse-direct solver | that option errors out at runtime with a diagnostic |

Build:

```bash
cd radcluster_code
cmake -S cpp_utils -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

The binary lands in `build/` (or `build/Release/` with multi-config generators) and is
launched by `py_utils/cpp_bridge.py`. The driver notebook builds it automatically on
first use and rebuilds it when the platform changes; set `FORCE_REBUILD = True` in the
build cell after editing C++ sources. On Windows the build also copies the required
SUNDIALS, LAPACK/MKL and KLU runtime DLLs next to the executable, so a fresh checkout
produces a self-contained binary.

Optimization flags are `-O3 -march=native` (or `/arch:AVX2 /GL /fp:precise` with MSVC).
Fast-math is deliberately *not* enabled: it alters IEEE-754 semantics and can make
CVODE's error estimator diverge on these stiff systems.

---

## 5. Quick start

### 5.1 Notebook (recommended)

Open `radcluster_code/codes/Notebooks/RadCluster_2_1.ipynb`. The notebook documents
every control variable inline, builds the C++ solver if needed, runs the case, prints
the end-of-run summary with conservation diagnostics, and writes a timestamped output
directory. A separate cell re-renders all figures from saved data without re-running the
solver, and another inspects the live reaction graph for the configured case.

### 5.2 Python API

```python
from radcluster_code.py_utils.simulation import RadClusterSimulation

sim = RadClusterSimulation(
    I=1000,                     # max SIA cluster size tracked
    V=1000,                     # max vacancy cluster size tracked
    solver_mode='full_system',  # 'full_system' | 'active_window'
    equations='bin_moment',     # 'discrete' | 'bin_moment'
    cascade='fusion',           # 'fusion' (He Case 1) | 'fission' (He Case 2)
    he_kinetics='quasi_steady_state',
    i_mobile=5, v_mobile=2,     # largest mobile SIA / vacancy cluster
)

print(sim.material_rag.summary())     # the EUROFER-97 reaction graph
print(sim.state_layout)               # vertex -> ODE-index map, N_eq

results = sim.run(solver_config={
    't_span': (1e-6, 1e8),            # seconds
    'n_points': 300,
    'log_time': True,
    'rtol': 1e-6, 'atol': 1e-20,
    'solver_method': {'linsol': 'gmres', 'preconditioner': 'Woodbury'},
})
```

`run_adaptive()` is the same call with automatic expansion of the cluster-size domain
when the distribution reaches the top of the grid. Material and environment parameters
come from the workbook (§6); anything can be overridden per run.

The two physics axes are orthogonal, and the four combinations map to the canonical
`physics_option` strings used on disk (`full_CD_fission`, `full_CD_fusion`,
`bin_moment_CD_fission`, `bin_moment_CD_fusion`); either form is accepted by the
constructor.

---

## 6. Inputs

All physics parameters live in one workbook, `radcluster_code/input/input_parameters.xlsx`,
organized into five sheets:

| Sheet | Contents |
|---|---|
| `Production` | survival efficiency, displacement and helium generation rates, clustered fractions, cascade power-law slopes and cut-offs — one column per spectrum (fission / fusion) |
| `Energetics` | lattice and elastic constants, atomic volume, formation energies, surface energy, loop Burgers vectors |
| `Diffusion` | pre-exponentials and migration energies for vacancies, SIAs and helium; 1D-glide parameters; solute concentrations and trap binding energies; mobility cut-offs |
| `Dissociation` | binding-energy models — void capillarity, helium virial EOS coefficients, loop binding fits, trap-mutation barriers |
| `Reactions` | geometric prefactors, dislocation bias factors, sink densities, cluster-size limits, discrete/bin split, solver settings |

Each row carries a human-readable parameter name, a `Symbol` key, a value, units, and a
note citing the equation or table it comes from. The `Symbol` key is what a run-time
override uses:

```python
PARAM_OVERRIDES = {'T': 673.0, 'rho_d': 5e14, 'eta': 0.28}   # temperature in K, m^-2, -
```

Overrides replace workbook values for that run only and are recorded in the run's
provenance file. Permanent changes belong in the workbook.

---

## 7. Outputs and provenance

Each run writes a self-describing directory:

```
output/YYYYMMDD_HHMMSS_<git-hash>/
├── provenance.md     # all input tables with overrides applied, user selections,
│                     # solver configuration, runtime, machine and CVODE statistics
├── results_t.npy     # time points
├── results_y.npy     # full ODE solution
├── summary.csv       # one-line end-of-run summary
├── diagnostics.txt   # conservation residuals and end-state diagnostics
└── plots/            # PNG figures, plus plot_data.pkl for re-plotting
```

`summary.csv` records the solver and physics option, temperature, dose rate, end time
and accumulated dose, total SIA / vacancy / helium concentrations, mean cluster sizes,
swelling, and the two conservation residuals $\delta_{\rm FP}$ and $\delta_{\rm He}$ —
one row per run, so a campaign's summaries concatenate directly into a table. The
solution arrays carry the full time history for any quantity post-processing needs.
Because the git SHA and the complete resolved parameter set are written with every run,
a figure can always be traced back to the exact code and inputs that produced it.

---

## 8. Examples, tests and verification

**Driver notebooks** — `radcluster_code/codes/Notebooks/RadCluster_2_1.ipynb` (full
workflow) and `EuroferExperiments.ipynb` (comparison against measured
microstructures). The 2_0 baseline driver, `RadCluster_2_0.ipynb`, is in the
archived module — restore it from the history as shown under *Module status*.

**Physics and numerics checks** — `radcluster_code/codes/Python_Testing/` holds standalone
scripts that gate the physics rather than merely exercising the code: helium and
Frenkel-pair conservation, bin-moment versus discrete equivalence, discrete-prefix
convergence, loop-conversion kernel and conservation checks, graph-declaration
consistency, linear-solver and window-mode comparisons, Woodbury benchmarks, and
parameter sweeps.

**Experimental database** — `radcluster_code/Eurofer_micro_database/` contains a curated
spreadsheet of loop and cavity number densities and mean sizes for neutron- and
ion-irradiated ferritic–martensitic steels, a notebook that fits and plots them, and the
resulting dose–temperature coverage maps and size distributions used for calibration.

**Digital twin / verification campaign** — `radcluster_code/digital_twin/` drives
systematic calibration and verification: a design generator, an ensemble runner keyed by
a parameter hash, a calibration ledger, a rescoring pass, a verification pass that emits
exactly the runs needed to confirm unverified observables under a change of grid extent,
and report/figure/table generators that write directly into
`docs/Formulation/verification_campaign/`. Campaigns are distributed across several
machines, including a UGE array job for the UCLA Hoffman2 cluster; results append
per-task so a restart skips completed rows and resubmission is always safe.

---

## 9. Computational aspects

**Solver.** Time integration uses SUNDIALS CVODE with variable-order BDF, appropriate
for the extreme stiffness of coupled defect kinetics (rates spanning many decades, and
physical times from microseconds to $10^8$ s). Two modes are available:

| Mode | Description |
|---|---|
| `full_system` | CVODE on the full state vector, with a dense, banded, GMRES, or KLU sparse-direct linear solver. |
| `active_window` | Two independent sliding windows — one over SIA size, one over vacancy size — with an OpenMP-parallel right-hand side. Thread count is auto-selected from the RHS sweep length and can be overridden with `OMP_NUM_THREADS`; the same code path runs serially when OpenMP is unavailable. |

**Preconditioning.** For the GMRES path the Jacobian has the structure $J = T + UV^{T}$,
with $T$ banded (half-bandwidth $\max(2 i_{\rm mobile}, 2 v_{\rm mobile})+1$) and a
rank-$r$ correction, $r = i_{\rm mobile} + v_{\rm mobile}$, from mobile-species coupling.
The Woodbury preconditioner exploits this exactly, using LAPACK banded factorization
(`dgbtrf`/`dgbtrs`) for $T$ and a dense factorization of the $r \times r$ Schur
complement. It is the default for `full_system` + GMRES; for the sliding-window modes the
active system is small enough (50–200 unknowns) that Jacobi + GMRES converges faster than
Woodbury's setup cost can be amortized. A colored finite-difference sparse Jacobian
feeds the KLU path.

**Cost control.** The bin-moment reduction is what makes engineering doses affordable:
the discrete formulation needs one ODE per cluster size, while the bin-moment
formulation integrates $O(10^2)$ unknowns for the same size range at a controlled
truncation error, with the closure order chosen per study. The helium reductions remove
the second $(v,\ell)$ dimension. The adaptive domain expands the tracked size range only
when the population actually reaches the boundary.

**Architecture and portability.** The Python and C++ trees share the same
core/materials split, so a new host material is added in parallel on both sides and the
solver layer is untouched. The Python walker and the C++ kernels are required to agree —
new physics is developed and validated in Python first, then mirrored in C++. The build
is platform-agnostic (macOS, Linux, Windows) with per-platform dependency discovery and,
on Windows, automatic runtime-DLL deployment.

**Reproducibility.** Provenance stamping (git SHA, machine, resolved parameters),
timestamped output directories, resumable campaign runs, and conservation residuals
reported on every run are treated as part of the solver contract rather than as
optional extras.

---

## 10. Documentation

| Location | Contents |
|---|---|
| `docs/Formulation/cluster_dynamics_framework/` | the main framework manuscript and paper sections, the Ghoniem–Cho lineage note, simulation methodology, novelty analyses |
| `docs/Formulation/reaction_admissibility_graph/` | RAG implementation and architecture supplement, EUROFER-97 graph figures, GraphML exports, structural reports |
| `docs/Formulation/reaction_kernels/` | one-dimensional migration reaction frequencies; detailed balance for effective kernels |
| `docs/Formulation/state_space_reduction_and_solvers/` | bin-moment derivation; GPU/CUDA integrator plan |
| `docs/Formulation/dislocation_loops/` | ½⟨111⟩→⟨100⟩ conversion kernels, loop character analysis, loop → network loss |
| `docs/Formulation/segregation_and_precipitation/` | radiation-induced segregation and solute-precipitation plan |
| `docs/Formulation/steel_applications_and_calibration/` | F/M-steel application, digital-twin appendix, calibration roadmap |
| `docs/Formulation/verification_campaign/` | verification plan, campaign results, approximation studies |
| `docs/Database/` | experimental radiation-microstructure databases for F/M steels |
| `docs/Literature/` | the reference library underlying the parameter set |
| `radcluster_code/CLAUDE.md` | the working physics and solver reference: state vector, kernels, reductions, solver options, parameter tables |

LaTeX sources sit next to their compiled PDFs and `.bib` files; each document compiles
from inside its own folder.

---

## 11. Publications and how to cite

If you use RadCluster, please cite the framework manuscript:

> N. M. Ghoniem, *A Generalized Graph-Based Cluster Dynamics Framework for Irradiated
> Materials*, Journal of Nuclear Materials (submitted, 2026).

```bibtex
@article{Ghoniem2026RadCluster,
  author  = {Ghoniem, N. M.},
  title   = {A Generalized Graph-Based Cluster Dynamics Framework for
             Irradiated Materials},
  journal = {Journal of Nuclear Materials},
  year    = {2026},
  note    = {Submitted}
}
```

The formulation continues a long line of work on defect cluster kinetics:

- N. M. Ghoniem, *A Dynamic Rate Theory for the Response of Metals during Steady State
  and Pulsed Irradiation*, Ph.D. thesis, University of Wisconsin–Madison (1977).
- N. M. Ghoniem and D. D. Cho, *The simultaneous clustering of point defects during
  irradiation*, Phys. Status Solidi A **54** (1979) 171–178.
  [doi:10.1002/pssa.2210540122](https://doi.org/10.1002/pssa.2210540122)
- N. M. Ghoniem and S. Sharafat, *A numerical solution to the Fokker–Planck equation
  describing the evolution of the interstitial loop microstructure during irradiation*,
  J. Nucl. Mater. **92** (1980) 121–135.
  [doi:10.1016/0022-3115(80)90148-8](https://doi.org/10.1016/0022-3115(80)90148-8)
- N. M. Ghoniem and D. D. Cho, *The early stages of void and interstitial loop evolution
  in pulsed fusion reactors*, J. Nucl. Mater. **89** (1980) 359–371.
  [doi:10.1016/0022-3115(80)90068-9](https://doi.org/10.1016/0022-3115(80)90068-9)
- N. M. Ghoniem, J. N. Alhajji and D. Kaletta, *The effect of helium clustering on its
  transport to grain boundaries*, J. Nucl. Mater. **136** (1985) 192–206.
  [doi:10.1016/0022-3115(85)90007-8](https://doi.org/10.1016/0022-3115(85)90007-8)
- N. M. Ghoniem, *Stochastic theory of diffusional planar-atomic clustering and its
  application to dislocation loops*, Phys. Rev. B **39** (1989) 11810–11819.
  [doi:10.1103/PhysRevB.39.11810](https://doi.org/10.1103/PhysRevB.39.11810)

Full reference lists, including the atomistic and experimental sources behind the
parameter set, are in the `.bib` files under `docs/Formulation/`.

---

## 12. Repository conventions

- **Module layout is uniform.** Every active module follows the same
  `py_utils/ · cpp_utils/ · codes/ · input/ · output/ · build/` structure, so moving
  between versions requires no re-orientation.
- **Physics lives in the workbook, not in the code.** Parameters carry a `Symbol` key, a
  unit and a source citation; code reads them, it does not hard-code them.
- **Outputs are immutable and stamped.** Runs never overwrite one another;
  `output/` and `build/` are gitignored.
- **New physics is Python-first.** Develop in the graph walker, pass the conservation
  gate, mirror in C++, then calibrate.
- **Archives are read-only, and local.** Superseded modules stay unchanged in
  `archive/`, which is gitignored; recover them from the history as shown under
  *Module status* when a reproducibility question needs them.
- **Local working state is not versioned.** `.claude/` (assistant session state)
  and `Claude outputs/` (a scratch drop for generated drafts) are gitignored.
  `scripts/link-memory.sh` belonged to a retired scheme for sharing that state
  across machines and is kept only because earlier commits reference it.

Issues and pull requests are welcome. Changes that touch reaction kernels or
stoichiometry should come with the relevant conservation check from
`codes/Python_Testing/` and a note of the resulting $\delta_{\rm FP}$ and
$\delta_{\rm He}$.

---

## 13. License and contact

Released under the [MIT License](LICENSE), © 2026 Nasr M. Ghoniem.

**Nasr M. Ghoniem** — Mechanical and Aerospace Engineering Department,
University of California, Los Angeles · ghoniem@ucla.edu
