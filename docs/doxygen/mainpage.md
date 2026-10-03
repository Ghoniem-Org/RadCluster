# RadCluster code reference {#mainpage}

RadCluster evolves the populations of radiation-induced defects in a crystalline solid: self-interstitial (SIA) dislocation loops, vacancy and helium–vacancy cavities, and free helium. A material is declared as a **Reaction Admissibility Graph** whose edges are elementary reactions of ten abstract classes. A Python reference right-hand side is assembled from that graph; the production equations are integrated by a C++ CVODE solver. The reference host material is the steel EUROFER-97.

This is the reference of the code, generated from its docstrings and header comments. For what the code does, its methodology and its verification, read the [README](https://github.com/Ghoniem/RadCluster/blob/main/README.md); for how to use it, the [wiki](https://github.com/Ghoniem/RadCluster/wiki).

> **Names.** The active module is `radcluster_code/`. Its Python package is `py_utils` and its C++ sources are in `cpp_utils`. Run labels, file headers and some docstrings still say `RadCluster_2_1`, the module's name before version 2.1.1.

## Where to start

| To find | Look in |
|---|---|
| how a cluster and a reaction are described | `py_utils.core`: `ClusterIdentifier`, `Polarity`, `Population`, `EdgeClass`, `EdgeClassSpec`, `Edge`, `ReactionAdmissibilityGraph` |
| how the graph becomes equations | `py_utils.core.graph_walker.GraphWalker`, `py_utils.core.state_layout.StateLayout`, `py_utils.core.reductions` (`BinMomentReduction`, `HeReductionMode`) |
| how EUROFER-97 is declared | `py_utils.materials.eurofer97.declaration.build_eurofer_rag` |
| where the parameters come from | `py_utils.input_data.InputData` (the workbook), `py_utils.reaction_rates.ReactionRates`, `py_utils.binding_energies`, `py_utils.defect_production`, `py_utils.loop_energetics` |
| one call that runs a case | `py_utils.simulation.RadClusterSimulation` (`run`, `run_adaptive`) |
| how Python starts the C++ solver | `py_utils.cpp_bridge` (`write_param_file`, `run_cpp_solver`) |
| the C++ solver | `solver.cpp` (CVODE driver), `parameters.h` (`Parameters`), `rate_equations.h`, `rhs_dispatch.cpp` (preconditioners), `sparse_jacobian.h`, `rate_kernels.cpp` (the EUROFER-97 right-hand sides) |
| what a run reports | `py_utils.post_process`, `py_utils.visualization` |
| the equations, kernels and solver options in one place | the page [Physics and Solver Reference](radcluster_code/CLAUDE.md) |

## The two layers

| Layer | Python | C++ | Contents |
|---|---|---|---|
| 1: core | `py_utils.core` | `cpp_utils/core/` | in Python, host-independent: the cluster identifier, the ten edge classes, the graph, the graph walker, the state layout, the reductions. In C++, the solver layer: the CVODE driver, the preconditioners and the sparse Jacobian |
| 2: material | `py_utils.materials.eurofer97` | `cpp_utils/materials/eurofer97/` | the EUROFER-97 populations, edges and rate kernels (processes P1–P8) |

The Python core carries no assumption about bcc iron or EUROFER-97. A new host material is a new `materials/<name>/` directory on each side.

The Python graph walker is the reference right-hand side, assembled from the graph. The C++ kernels (`rate_kernels.cpp`) are written by hand for EUROFER-97 and are the production solver; `py_utils.simulation.RadClusterSimulation` integrates only through them, and carries the graph as the description of the run. The Python right-hand sides in `py_utils.rate_equations` and `py_utils.bin_moment_rates` (`ode_system`) are kept as a physics reference and are not called by a run.

## Scripts

Python files outside a package appear as top-level namespaces named after the file.

| Group | Folder | Examples |
|---|---|---|
| physics and numerics checks | `radcluster_code/codes/Python_Testing/` | `check_eurofer_rag`, `check_he_conservation`, `check_bm_vs_discrete`, `gate_im5vm2`, `compare_linsol` |
| figure generators | `radcluster_code/codes/` | `make_rag_figure`, `make_rag_hyperedge_figure`, `make_dose_figures`, `make_size_effect_figures` |
| calibration campaign | `radcluster_code/digital_twin/` | `design`, `run_ensemble`, `merge_and_sobol`, `learn`, `plan`, `verify` |
| verification campaign | `radcluster_code/digital_twin/verification/` | `runs`, `campaign`, `gitsync`, `make_tables` |

## What is not here

- The driver notebook, `radcluster_code/codes/Notebooks/RadCluster_2_1.ipynb`: see the wiki page *Running a simulation*.
- The input workbook, `radcluster_code/input/input_parameters.xlsx`, and the experimental database.
- The derivations, the manuscript and the study reports: `docs/Formulation/`.
- Run outputs and campaign results.

## License

MIT; see `LICENSE`. Primary reference: N. M. Ghoniem, *A Generalized Graph-Based Cluster Dynamics Framework for Irradiated Materials*, Journal of Nuclear Materials (submitted, 2026).
