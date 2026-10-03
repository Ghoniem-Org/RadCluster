# Solver

Time integration is by SUNDIALS CVODE with the variable-order BDF method. Defect kinetics are very stiff: rates span many decades, and a run goes from microseconds to 10⁸ s. The solver is a C++17 program, `radcluster_code/cpp_utils/`, started by `py_utils/cpp_bridge.py`.

## Modes

| `solver_mode` | Description |
|---|---|
| `full_system` | CVODE on the full state vector |
| `active_window` | two independent sliding windows, one over SIA size and one over vacancy size. Only the sizes inside a window are advanced; a window widens when the concentration at its edge exceeds a threshold |

The earlier names `cpp_full` and `sliding_OpenMP` are still accepted. When `I` and `V` are both below 500, `active_window` runs as the full solver.

`active_window` is not free in `I`: the right-hand side is restricted to the window, but CVODE's vector arithmetic runs over the whole state vector. For large size ranges the bin-moment equations are the economical choice.

## Linear solvers

| `linsol` | Method | Notes |
|---|---|---|
| `dense` | dense direct | small systems |
| `band` | banded direct | the bandwidths are `solver_method['mu']` and `['ml']` |
| `gmres` | preconditioned GMRES | the choice for large systems. Forced for the bin-moment equations when the default configuration is used |
| `klu` | sparse direct (SuiteSparse KLU) with a coloured finite-difference Jacobian | discrete equations only; needs SUNDIALS built with KLU |

## Preconditioners for GMRES

| `preconditioner` | Description |
|---|---|
| `Jacobi` | intended as diagonal scaling. In the current code the diagonal is never filled, so this option applies no preconditioning: GMRES runs unpreconditioned |
| `Woodbury` | built for the structure of the Jacobian. Needs LAPACK; without it the solver falls back to the option above |

The Jacobian has the structure J = T + U Vᵀ:

- **T** is banded, with half-bandwidth max(2 `i_mobile`, 2 `v_mobile`) + 1.
- **U Vᵀ** has rank r = `i_mobile` + `v_mobile`. It is the coupling of every equation to the concentrations of the mobile species.

The Woodbury preconditioner solves with this structure by the Sherman–Morrison–Woodbury formula: a banded LAPACK factorization of T (`dgbtrf`, `dgbtrs`) and a dense factorization of the r × r Schur complement (`dgetrf`, `dgetrs`).

The band and the coupling columns are built by finite differences of the right-hand side. The ⟨100⟩ loop block is not part of the modelled structure.

**Default.** Give `preconditioner` explicitly. When it is absent, `cpp_bridge` chooses Woodbury for GMRES unless `solver_method['window_mode']` is 4, and only the default configuration built by `RadClusterSimulation` sets that key. The intent is Woodbury for `full_system` and none for `active_window`, where the active system is small (50 to 200 unknowns) and Woodbury's setup cannot be repaid; with your own `solver_config` in `active_window` mode, state `'preconditioner': 'Jacobi'` to get it.

The derivation is in `docs/Formulation/reaction_admissibility_graph/supplement_S1/`.

## Configuration

`solver_config` is a dictionary passed to `run` or `run_adaptive`.

| Key | Meaning | Default in `cpp_bridge` |
|---|---|---|
| `t_span` | start and end times, s | `(1e-8, 1e7)` |
| `n_points` | number of output times | 200 |
| `log_time` | logarithmic spacing of the output times | `True` |
| `rtol`, `atol` | relative and absolute tolerances | `1e-8`, `1e-50` |
| `timeout_s` | wall-clock limit; the solver then stops gracefully and the partial run is kept | none |
| `loop_conversion` | 1: integrate the ⟨100⟩ loop population and the conversion | the workbook value, else 0 |
| `required_times` | times that must be in the output (`run_adaptive`) | none |
| `solver_method` | the dictionary below | |

`solver_method`:

| Key | Meaning | Default |
|---|---|---|
| `linsol` | `dense`, `band`, `gmres`, `klu` | `dense` |
| `preconditioner` | `Jacobi` or `Woodbury` | see above |
| `prec_bw`, `prec_rank` | half-bandwidth and rank for Woodbury | from the mobility cut-offs |
| `hmax`, `hmin` | largest and smallest internal step; 0 means no limit | 0 |
| `max_order` | largest BDF order; 0 means CVODE's default | 0 |
| `window_width` | initial width of a window | max(`I`, `V`) |
| `concentration_threshold` | the edge concentration that widens a window | `1e-18` |
| `window_pad`, `window_pad_v` | sizes added when the SIA or vacancy window widens | 10 |
| `window_check_every` | output steps between window checks | 1 |

When no `solver_config` is given, `RadClusterSimulation` builds one from the `Reactions` sheet where those symbols exist; the shipped workbook has none of them, so its own defaults apply: `linsol` `dense` (`gmres` for the bin-moment equations), `atol` 10⁻²⁰, `window_width` 100, a wall-clock limit of one hour, and an internal step capped at a tenth of the end time, which keeps BDF from stepping across many decades at once near a steady state.

## Threads

The right-hand side has OpenMP-parallel sweeps. The solver chooses the number of threads from the length of the sweep, the larger of the number of equations and `I + V`: one thread below 500, then 2, 4, 8, 12, 16, 20, and all cores from 80 000.

- The number of equations alone is the wrong measure for the bin-moment equations, which integrate about a hundred unknowns but sweep the reconstructed per-size distributions.
- On Apple silicon the choice is capped at the number of performance cores. The parallel regions wait for all threads, so a thread on an efficiency core delays every one.
- `OMP_NUM_THREADS` overrides the choice. Without OpenMP the same code runs on one thread.
- On a Windows machine with more than 64 logical processors, `RADCLUSTER_CPU_GROUP` keeps a solver on one processor group.

## The state vector

The discrete equations, by cascade:

| Segment | `fission` (helium Case 2) | `fusion` (helium Case 1) |
|---|---|---|
| SIA clusters, sizes 1 to I | `y[0 .. I-1]` | `y[0 .. I-1]` |
| cavities, sizes 1 to V | `y[I .. I+V-1]` | `y[I .. I+V-1]` |
| helium in cavities | one scalar, `y[I+V]` | one per size, `y[I+V .. I+2V-1]` |
| free helium | one entry; absent with `he_kinetics='quasi_steady_state'` | the same |
| cumulative fluxes for the residuals | five entries | five entries |
| ⟨100⟩ loops, with the loop split on | appended last; `cpp_bridge.sia100_block_length` gives the length | the same |

- In the bin-moment equations each size block is the discrete sizes followed by P moments per bin.
- `results['y']` is this vector without the ⟨100⟩ block, which is returned separately as `results['y_sia100']`.
- `sim.state_layout.describe()` prints the blocks of the declared graph. That layout always contains a ⟨100⟩ block and a free-helium entry, and orders the blocks differently, so its `N_eq` is not the number of equations the solver integrates. `sim.rate_equations.N_eq` is that number without the ⟨100⟩ block; the solver prints the full count when it starts.

## Numerical safeguards

- **Concentration floor.** `C_floor` is applied to the state after each output step, not inside the right-hand side. A clamp inside the right-hand side puts a kink in the Jacobian at the floor, which breaks the BDF corrector. The right-hand side uses max(y, 0).
- **No fast-math.** The build does not enable it; it changes IEEE-754 arithmetic and can make CVODE's error estimate diverge on these systems.
- **Size boundary.** By default (`boundary_flux = absorption`) a reaction product larger than the grid is lost; with `reflection` the reaction is suppressed. The growth of ⟨100⟩ loops is halted at the top of the grid in either case. `run_adaptive` enlarges the grid when the content near the top exceeds a threshold.
- **Failed steps are not written.** Only output times that the solver reached are in the output.

## The interface between Python and C++

| Step | Detail |
|---|---|
| parameter file | `write_param_file` writes a text file, one `key=value` per line; an array is written as indexed entries (`name_0=`, `name_1=`, …) |
| start | `solver --param_file=<path>` |
| output | a binary file beside the parameter file: one row of doubles per output time, the time followed by the state vector. A text file `<output>.window.csv` records the upper window bounds |
| progress | one `[cvode]` line per output time and a final `[stats]` line on standard error |
| statistics | `steps`, `nfe` (right-hand sides), `nni` (nonlinear iterations), `nli` (linear iterations), `npe` and `nps` (preconditioner setups and solves), `ncfn` and `netf` (convergence and error-test failures), `nlsetup` |
| stop | SIGINT or SIGTERM (a console event on Windows): the solver finishes its step, writes it and exits |

In the parameter file the options are integers: `window_mode` 0 or 4 (the solver mode), `physics_option_int` 0 to 3, `linsol` 0 to 3, `prec_type` 0 or 1.

## Performance notes

- The bin-moment reduction is what makes engineering doses affordable: it integrates of the order of a hundred unknowns where the discrete equations need one per size.
- In bin-moment mode the cost of a right-hand side grows with `I` and `V` even though the number of equations does not.
- `radcluster_code/codes/Python_Testing/compare_linsol.py`, `benchmark_woodbury.py` and `compare_window_modes.py` time the linear solvers, the preconditioners and the two modes.
- A plan for a GPU integrator is in `docs/Formulation/state_space_reduction_and_solvers/Cuda_C++_Plan.md`.
