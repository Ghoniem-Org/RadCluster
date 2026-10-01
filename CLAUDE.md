# RadCluster — Project Overview

Physics-based simulation suite for **EUROFER97 / ferritic-martensitic steel** behaviour under irradiation and thermal loading. Modelled after the structure of the `Fluor_Zr` repository.

## Repository Layout

```
RadCluster/
├── RadCluster_2_1/     # Active development — adds dislocation evolution + RIS (clone of 2_0)
├── archive/            # Archived modules (read-only, kept for reproducibility)
│   ├── RadCluster_2_0/     # Archived 2026-10-01 — released as v2.0.0; baseline for RadCluster_2_1
│   ├── Eurofer/            # Archived earlier microstructure work
│   ├── Eurofer_CD/         # Archived 2026-05-02 — superseded by RadCluster_1_0
│   ├── RadCluster_1_0/     # Archived 2026-06-13 — superseded by RadCluster_2_0
│   ├── Monomer_CD/         # Archived 2026-06-13 — cluster dynamics scaling reference
│   └── Zr_RadCluster_1_0/  # Archived 2026-06-13 — zirconium cluster-dynamics variant
├── docs/               # Shared documentation, literature, databases
│   ├── Database/           # Experimental microstructure databases (xlsx)
│   ├── Formulation/        # Derivations, paper sections, studies — by topic (see Formulation/README.md)
│   ├── Literature/         # Peer-reviewed papers (PDF)
│   └── design_notes/       # Working design notes
├── scripts/            # Repository utilities (link-memory.sh)
├── README.md           # Public-facing introduction: methodology, install, usage, citation
├── requirements.txt
└── LICENSE             # MIT
```

`README.md` is the user-facing entry point; this file (`CLAUDE.md`) is the internal
working reference. When the input workbook, the output artifacts, the solver options or
the module layout change, both need updating — they describe the same contract.

## Module Status

| Module | Status | Description |
|---|---|---|
| `RadCluster_2_1/` | **Active (dev)** | Clone of `RadCluster_2_0` adding two capabilities: **(a)** a dislocation-evolution / loop→network loss edge that saturates loop density, and **(b)** Radiation-Induced Segregation (RIS) + solute precipitation. Plans: [`docs/Formulation/dislocation_loops/loop_network_loss.tex`](docs/Formulation/dislocation_loops/loop_network_loss.tex) and [`docs/Formulation/segregation_and_precipitation/radcluster_2_1_RIS_plan.tex`](docs/Formulation/segregation_and_precipitation/radcluster_2_1_RIS_plan.tex). See `RadCluster_2_1/CLAUDE.md` §0. |
| `archive/RadCluster_2_0/` | Archived 2026-10-01 | Generalized graph-based cluster dynamics (Ghoniem 2026). Two-layer RAG architecture (abstract core + EUROFER-97 host declaration). Released as git tag `v2.0.0`; stable reference for the 2_1 work. Notebooks: `RadCluster_2_0.ipynb` (simulation driver) and `EuroferExperiments.ipynb`. |
| `archive/RadCluster_1_0/` | Archived 2026-06-13 | Superseded by `RadCluster_2_0/`. Earlier (non-graph) generalized cluster dynamics. |
| `archive/Monomer_CD/` | Archived 2026-06-13 | Monomer-mobility cluster dynamics scaling reference (Ghoniem & Cho 1979, no He). |
| `archive/Zr_RadCluster_1_0/` | Archived 2026-06-13 | Zirconium cluster-dynamics variant of `RadCluster_1_0`. |
| `archive/Eurofer/` | Archived | Earlier EUROFER microstructure notebooks. |
| `archive/Eurofer_CD/` | Archived 2026-05-02 | Superseded by `RadCluster_1_0/`. Last active state at git tag `eurofer_cd-final`. |

## Conventions (mirror Fluor_Zr)

### Module structure
Each active module follows the two-layer `core/` + `materials/` split on both the
Python and the C++ side, so adding a host material never touches the solver layer:
```
<Module>/
├── CLAUDE.md           # Physics description, solver notes
├── codes/              # Jupyter notebooks, test scripts, figure generators
│   ├── Notebooks/          # <Module>.ipynb — the run driver
│   ├── Python_Testing/     # conservation gates, convergence and sweep studies
│   └── make_*_figures.py   # RAG / dose / size-effect figure generators
├── py_utils/           # Python utilities package
│   ├── core/               # Layer 1 — abstract, host-independent
│   │   ├── cluster_identifier.py, edge_classes.py, rag.py
│   │   ├── graph_walker.py, state_layout.py, to_networkx.py
│   │   └── reductions/     # bin_moment.py, he_reduction.py
│   ├── materials/          # Layer 2 — host declarations
│   │   └── eurofer97/declaration.py   # build_eurofer_rag()
│   ├── input_data.py       # InputData class (reads the 5-sheet workbook)
│   ├── reaction_rates.py, rate_equations.py, bin_moment_rates.py
│   ├── binding_energies.py, defect_production.py, loop_energetics.py
│   ├── simulation.py       # Orchestrator; writes timestamped output/
│   ├── cpp_bridge.py       # subprocess wrapper for C++ solver
│   └── post_process.py, size_distributions.py, visualization.py
├── cpp_utils/          # C++17 solver (production engine)
│   ├── core/               # solver.cpp (CVODE driver), rhs_dispatch.cpp,
│   │                       # sparse_jacobian.cpp, parameters.h, rate_equations.h
│   ├── materials/eurofer97/rate_kernels.cpp   # P1–P8 kernel arithmetic
│   └── CMakeLists.txt
├── input/              # input_parameters.xlsx (5 sheets)
├── output/             # Timestamped run directories (gitignored)
└── build/              # CMake build artifacts (gitignored)
```

### Excel input format
All modules share the same 5-sheet structure. Every row carries a parameter name, a
`Symbol` key, a value, units, and a note citing the equation or table it comes from;
the `Symbol` key is what `PARAM_OVERRIDES` targets.
- `Production` — survival efficiency `eta`, displacement rate `phi_dot`, He production
  ratio `G_He/G`, clustered fractions `f_cl_i`/`f_cl_v`, power-law exponents `s_i`/`s_v`,
  in-cascade size cut-offs `m1`/`n1`. One column per spectrum (Fission / Fusion).
- `Energetics` — lattice constant `a`, atomic volume `Omega`, elastic constants,
  formation energies, surface energy `gamma_s`, Burgers vectors `b_111`/`b_100`.
- `Diffusion` — `D0`/`E_m` for vacancy, SIA and He; 1D-glide parameters (`E_m_1D`,
  `L_hat`); dissolved solute concentrations and trap binding energies; mobility
  cut-offs `i_mobile`/`v_mobile`.
- `Dissociation` — binding-energy models: void capillarity, He virial EOS coefficients,
  loop binding fits (`A_111`/`B_111`, `A_100`/`B_100`), trap-mutation barriers.
- `Reactions` — geometric prefactors (`A_sph`, `A_loop`, ...), dislocation bias factors
  `Z_i`/`Z_v`, sink densities `rho_d`, cluster-size limits, the discrete/bin split
  (`i_discrete`, `I_bin`, `shape_function`), and solver settings (`linsol`, `rtol`,
  `atol`, `t_begin`, `t_end`, `window_*`, `hmax`, `timeout_s`).

Historical note: the sheet set was once `Material_Environment` / `Physical_Properties` /
`Model_Parameters`. Those names are gone — `InputData` exposes the five sheets above as
`production_fission`, `production_fusion`, `energetics`, `diffusion`, `dissociation`,
`reactions`, plus a computed `derived` dict.

### Output format
Simulation results go in timestamped subdirectories:
```
output/YYYYMMDD_HHMMSS_<git-hash>/
├── provenance.md       # all input tables with overrides applied, user selections,
│                       # solver config, runtime + machine info, CVODE statistics
├── results_t.npy       # time points
├── results_y.npy       # full ODE solution
├── summary.csv         # one-line end-of-run summary (see below)
├── diagnostics.txt     # conservation residuals and end-state diagnostics
└── plots/              # PNG figures (+ plot_data.pkl for the replot cell)
```

`summary.csv` is written by `post_process.summary_csv_row()` and holds a single row:
`solver`, `physics_option`, `T_K`, `G_dpa_s`, `t_end_s`, `dose_dpa`, `C_SIA_tot_m3`,
`C_VAC_tot_m3`, `C_He_tot_m3`, `mean_n_i`, `mean_n_v`, `swelling`, `delta_FP`,
`delta_He` — so a campaign's summaries concatenate directly into one table.

There is no `results.pkl`: the solution is saved as the two `.npy` arrays above.
Each writer is individually guarded, so a failure in one artifact never costs the
others — the raw arrays are written first.

### C++ solvers
- SUNDIALS 7.1.1 (CVODE, variable-order BDF). `CMakeLists.txt` looks first in a
  `Libraries/sundials-7.1.1/` tree that is a sibling of the repository root
  (`mac_install/` / `win_install/` per platform), then falls back to whatever
  `find_package(SUNDIALS)` resolves (system or vcpkg install).
- C++17, CMake ≥ 3.15 (the build uses `list(PREPEND)`); binaries land in `build/`
  or `build/Release/` (gitignored). Build with
  `cmake -S cpp_utils -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --config Release`.
- Optional dependencies, each auto-detected and each degrading gracefully:
  LAPACK/BLAS (Accelerate, OpenBLAS, reference LAPACK or MKL) → Woodbury
  preconditioner, else compiled out and Jacobi is used at runtime; OpenMP → parallel
  RHS in `active_window`, else serial; SuiteSparse/KLU with SUNDIALS built
  `-DENABLE_KLU=ON` → `linsol='klu'`, else that option errors out with a diagnostic.
- `-O3 -march=native` (`/arch:AVX2 /GL /fp:precise` on MSVC). Fast-math is
  deliberately NOT enabled — it changes IEEE-754 semantics and can make CVODE's error
  estimator diverge on these stiff systems.
- On Windows the build copies the required SUNDIALS / LAPACK-MKL / KLU runtime DLLs
  next to the executable, so a fresh checkout yields a self-contained binary.
- Invoked from Python via `py_utils/cpp_bridge.py` subprocess wrapper. The notebook's
  build cell rebuilds automatically when the platform-correct binary is missing; set
  `FORCE_REBUILD = True` there after editing C++ sources.

## Environment

```bash
# Install dependencies
pip install -r requirements.txt

# Register Jupyter kernel (optional)
python -m ipykernel install --user --name radcluster --display-name "RadCluster"

# Link Claude Code's project memory to the copy in this repo (once per machine)
./scripts/link-memory.sh
```

### Project memory

Claude Code's memory for this project is committed at `.claude/memory/` and
symlinked into `~/.claude/projects/<slug>/memory/`, because the harness path is
keyed to the absolute checkout path and would otherwise exist on one machine
only. Run `scripts/link-memory.sh` after cloning; without it a session on this
machine starts with no project memory and writes new memories somewhere git
never sees. `.claude/memory/MEMORY.md` is the index.

## Shared Resources

- `docs/Database/` — experimental radiation microstructure data for ferritic-martensitic steels
- `docs/Formulation/` — rate-equation derivations (canonical)
- `docs/Literature/` — reference papers
