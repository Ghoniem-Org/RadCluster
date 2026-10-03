# Architecture

The specification is the framework manuscript, [Generalized_Cluster_Dynamics.pdf](https://github.com/Ghoniem/RadCluster/blob/main/docs/Formulation/cluster_dynamics_framework/Generalized_Cluster_Dynamics.pdf), with its implementation supplement in `docs/Formulation/reaction_admissibility_graph/supplement_S1/`. This page is a summary.

## The rule

> A material is declared, not coded.

A model is the labelled multidigraph G = (V, E, ν, k):

- **V**: the admissible cluster vertices;
- **E**: every admissible elementary reaction;
- **ν**: the stoichiometric action of each edge, fixed by its class;
- **k**: the library of rate kernels.

The graph is specific to the host in (V, E) and to the calibration in k. The algorithm that turns it into equations is the same for every host.

## Two layers

| Layer | Python | C++ | Contents |
|---|---|---|---|
| 1: core | `py_utils/core/` | `cpp_utils/core/` | the graph machinery in Python; the CVODE driver, preconditioners and sparse Jacobian in C++ |
| 2: material | `py_utils/materials/eurofer97/` | `cpp_utils/materials/eurofer97/` | the EUROFER-97 populations, edges and kernels |

- The Python core carries no assumption about bcc iron or EUROFER-97.
- The C++ core is the solver layer. Its preconditioners are independent of the material. Its parameter structure (`parameters.h`) and its driver (`solver.cpp`) still carry EUROFER-97 fields and name the EUROFER-97 right-hand sides, and the sparse Jacobian assumes the discrete EUROFER-97 layout.
- A host declares its populations, adds edge families and registers rate kernels as precomputed arrays. The core never derives a rate.
- A new host is a new `materials/<name>/` directory on each side. See [Adding a host material](Adding-a-Host-Material).

## The Python package

`radcluster_code/py_utils/`

| Module | Role | Main objects |
|---|---|---|
| `core/cluster_identifier.py` | the cluster four-tuple | `Polarity`, `Population`, `ClusterIdentifier` |
| `core/edge_classes.py` | the ten reaction classes and the stoichiometric contract of each | `EdgeClass`, `EdgeClassSpec`, `EDGE_CLASS_SPEC` |
| `core/rag.py` | the graph: populations, edge families, kernels, admissibility | `Edge`, `ReactionAdmissibilityGraph` |
| `core/graph_walker.py` | the right-hand side, by walking the edges | `GraphWalker.assemble(t, y)` |
| `core/state_layout.py` | which unknown is where in the state vector | `Block`, `StateLayout` |
| `core/reductions/` | the bin-moment and helium reductions | `BinMomentReduction`, `HeReductionMode` |
| `core/to_networkx.py` | the graph as an explicit NetworkX graph, for diagnostics | `to_networkx`, `structural_report`, `write_graphml` |
| `materials/eurofer97/declaration.py` | the EUROFER-97 graph and layout | `build_eurofer_rag` |
| `input_data.py` | the workbook and the derived quantities | `InputData` |
| `reaction_rates.py` | every size-resolved rate constant | `ReactionRates` |
| `binding_energies.py`, `defect_production.py`, `loop_energetics.py` | binding energies, the cascade source, loop free energies | `E_b_void`, `E_b_loop_i`, `production_rates`, `LoopEnergetics` |
| `rate_equations.py`, `bin_moment_rates.py` | the layout of the unknowns, initial conditions, bins and closures | `RateEquations`, `BinMomentRateEquations`, `build_bins`, `reconstruct_distribution` |
| `simulation.py` | the orchestrator | `RadClusterSimulation` |
| `cpp_bridge.py` | writes the parameter file, starts the solver, reads its output | `write_param_file`, `run_cpp_solver` |
| `post_process.py` | observables and conservation residuals | `calculate_derived_quantities`, `summary_csv_row` |
| `visualization.py` | the figures | `save_all_plots`, `set_plot_config` |

The `ode_system` methods of `RateEquations` and `BinMomentRateEquations` are a Python right-hand side kept as a physics reference. A run never calls them.

## The C++ solver

`radcluster_code/cpp_utils/`

| File | Role |
|---|---|
| `core/solver.cpp` | `main`: reads the parameter file, sets up CVODE and its linear solver, integrates, moves the windows, writes the binary output, handles interrupts |
| `core/parameters.h` | `Parameters`: everything the right-hand side needs, unpacked from the parameter file |
| `core/rate_equations.h` | the right-hand-side and preconditioner callbacks, and `UserData` |
| `core/rhs_dispatch.cpp` | the Jacobi and Woodbury preconditioners |
| `core/sparse_jacobian.cpp`, `.h` | the sparsity pattern, its column colouring and the finite-difference sparse Jacobian for KLU |
| `materials/eurofer97/rate_kernels.cpp` | the EUROFER-97 right-hand sides: `rhs_full_CD` (helium Case 1 and Case 2) and `rhs_bin_moment` |
| `CMakeLists.txt` | the build |

The preconditioners reach the right-hand side only through a function pointer, so `rhs_dispatch.cpp` does not depend on the material at compile time.

## How a run flows

1. `InputData` reads the workbook and computes the derived quantities.
2. `ReactionRates` computes the rate-constant arrays.
3. `RateEquations` or `BinMomentRateEquations` fixes the layout of the unknowns and the initial conditions.
4. `build_eurofer_rag` declares the reaction graph and its `StateLayout`. They are attached to the simulation as `sim.material_rag` and `sim.state_layout`.
5. `cpp_bridge.write_param_file` writes one `key=value` line per parameter, with arrays as indexed entries.
6. `solver --param_file=<path>` integrates and writes one row per output time (time, then the state vector) as binary doubles.
7. `cpp_bridge.run_cpp_solver` reads that file; `post_process` computes the observables and residuals.
8. `simulation` writes the run directory; `visualization` writes the figures.

## Python and C++

- **New physics is developed in Python first.** It is written as edges and kernels for the graph walker, it passes the conservation check there, it is mirrored in `rate_kernels.cpp`, and only then is it calibrated.
- **The graph walker is the executable specification.** It is a discrete, one-equation-per-size accumulator that dispatches on the edge class alone. It has no composition axis: helium is carried by the helium reduction, and the walker raises an error on a solute-trapping edge.
- **C++ is the production solver.** Every simulation integrates through `rate_kernels.cpp`, whose right-hand sides are written by hand for EUROFER-97. They are not generated from the graph.
- **The two are not yet identical.** The kernels registered on the graph for cluster–cluster reactions, and the cascade source registered there, are not the ones the C++ solver uses; see [EUROFER-97 reaction graph](EUROFER-97-Reaction-Graph). The checks that use the walker run small graphs built for the purpose; no check compares the walker with the C++ right-hand side on the complete EUROFER-97 graph.
- **The graph is the description of the run.** `sim.material_rag` says which reactions are admissible and which reductions are active. Its construction is advisory: if it fails, a warning is issued and the run proceeds, because the C++ path does not read it.

Two points follow from the last item and are recorded in the code:

- The graph always declares both SIA loop populations. The solver integrates the ⟨100⟩ population only when `loop_conversion` is on.
- Trap mutation (P7) and radiation re-solution (P8) are declared as edges, with their kernels registered. The declaration notes that the C++ kernels do not evaluate them yet.

## Conservation by construction

Each edge class has a fixed stoichiometric contract: how many clusters it consumes, its kinetic order, whether it crosses polarity, and how it shifts size, population and composition. The sign pattern is therefore structural. The signed defect content q = χn, summed over all clusters, is conserved by every reaction inside the host, and changes only through sources and sinks. It is not re-derived for each material.

At run time the same statement is measured by two residuals, δ_FP and δ_He, reported for every run: [Tests and conservation checks](Tests-and-Conservation-Checks).

## Repository conventions

- **Physics lives in the workbook.** Parameters have a symbol, a unit and a source; the code reads them.
- **Outputs are never overwritten**, and `build/` and `output/` are not tracked. The one tracked run is the designated reference run.
- **Archived modules are not in a clone.** `archive/` is ignored; earlier modules are recovered from the history (`git checkout e9f92c3 -- archive/`).
- **Changes to a kernel or to stoichiometry** come with the relevant conservation check and the δ_FP and δ_He it produces.
