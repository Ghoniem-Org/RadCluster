# Running a simulation

The driver is one notebook: `radcluster_code/codes/Notebooks/RadCluster_2_1.ipynb`. Open it with a Python kernel that has the requirements installed, and run the cells in order.

| Cell | What it does |
|---|---|
| Introduction (text) | the architecture, the solver and physics options, every control variable, the loop-character split and what is saved |
| Build | finds the solver executable for this platform, or configures and builds it with CMake. Set `FORCE_REBUILD = True` after editing C++ sources |
| Run | plot settings, user selections, solver configuration, initialization, parameter overrides, the run, the save and the end-of-run summary |
| Material options | prints the live reaction graph: populations, edges by class, the conversion edges and the state layout |
| Replot | reads `plots/plot_data.pkl` of a run and redraws every figure with other axis limits and scales, without running the solver |

You normally edit only the Run cell.

## The Run cell

| § | Section | What you set |
|---|---|---|
| 1 | Plot axis controls | `PLOT_CONFIG`: limits and scales for three groups of figures (`concentration`, `scalar`, `size_dist`) |
| 2 | User selections | the switches, the domain, the mobility cut-offs, the discrete/bin split |
| 3 | Solver configuration | `SOLVER_CONFIG` |
| 4 | Initialize | builds `sim = RadClusterSimulation(...)`: reads the workbook and declares the graph |
| 4b | Parameter overrides | `PARAM_OVERRIDES`, applied on top of the workbook |
| 5 | Progress callback | collects one row per output time |
| 6 | Run | `sim.run_adaptive(...)`, then the save and the summary |

### Switches

| Variable | Values | Meaning |
|---|---|---|
| `DEBUG` | `True`, `False` | print every solver, override and diagnostic message |
| `MATERIAL` | `'eurofer97'` | the host material; the only one shipped |
| `SOLVER_MODE` | `'full_system'`, `'active_window'` | see [Solver](Solver) |
| `EQUATIONS` | `'discrete'`, `'bin_moment'` | one equation per size, or the bin-moment grouping |
| `CASCADE` | `'fission'`, `'fusion'` | the helium reduction: Case 2 (a decoupled inventory) for fission, Case 1 (a mean-field loading) for fusion. The cascade source follows `spectrum` in the workbook |
| `SHAPE_FUNCTION` | `'constant'`, `'linear'`, `'lognormal'` | the closure inside a bin: 1, 2 or 3 moments per bin |
| `HE_KINETICS` | `'quasi_steady_state'`, `'dynamic'` | eliminate free helium algebraically (recommended), or keep it as an unknown |
| `SIA_LOOP_SPLIT` | `True`, `False` | integrate ½⟨111⟩ and ⟨100⟩ loops as two populations with conversion between them. See [Loop character and network](Loop-Character-and-Network) |

### Domain

| Variable | Meaning |
|---|---|
| `I`, `V` | the largest SIA and vacancy cluster sizes tracked |
| `i_mobile`, `v_mobile` | the largest mobile SIA and vacancy clusters |
| `i_discrete`, `v_discrete` | sizes tracked one by one; larger sizes are binned |
| `I_bin`, `V_bin` | the target numbers of logarithmic bins above the discrete sizes. The realised number can differ by one |
| `C_FLOOR` | the concentration floor, in atomic fraction |
| `N_LOOP_MIN` | the smallest ⟨100⟩ loop tracked |

For one equation per size, set `EQUATIONS = 'discrete'`, `i_discrete = I`, `v_discrete = V`, `I_bin = V_bin = 0`.

- `V` sets the cost of every right-hand side in bin-moment mode, because the per-size vacancy distribution is reconstructed over the whole range. The number of equations does not change with it.
- The width of a bin matters more than the number of bins. Very wide bins make the linear reconstruction clamp often, and the integrator's step collapses.

### Solver configuration

```python
SOLVER_CONFIG = {
    't_span':   (1e-6, 4.0e8),      # seconds; 40 dpa at 1e-7 dpa/s
    'n_points': 40,
    'log_time': True,
    'rtol':     1e-5,
    'atol':     1e-20,
    'timeout_s': 172800,
    'loop_conversion': int(SIA_LOOP_SPLIT),
    'solver_method': {
        'linsol':                  'gmres',     # 'dense', 'band', 'gmres' or 'klu'
        'preconditioner':          'Woodbury',  # 'Jacobi' or 'Woodbury'
        'window_width':            10,
        'concentration_threshold': 1e-22,
        'window_pad':              20,
    },
}
```

Every key is described in [Solver](Solver).

## The reference case

As shipped, the Run cell reproduces the reference run: bin-moment equations, fission cascade, `I = 80000`, `V = 20000`, `i_mobile = v_mobile = 5`, to 40 dpa. `PARAM_OVERRIDES` is empty apart from the grid and bin controls, so the workbook alone sets the physics. The notebook records about 50 minutes for it on 12 threads.

The finished run is tracked in the repository, so a checkout shows what a run writes:

```
radcluster_code/output/20260906_063055_full_system_bin_moment_CD_fission_I80000V20000_im5vm5/
```

For a quick trend-only test of about three minutes, the Run cell lists the settings to change: `I = 4e4`, `v_mobile = 1`, `i_discrete = 40`, `I_bin = 20`, `V_bin = 24`, `rtol = 1e-6`, `t_span = (1e-6, 2.0e7)` and `timeout_s = 3600`. Trends on that grid are indicative; absolute densities and sizes are not.

## `run` and `run_adaptive`

| Method | What it does |
|---|---|
| `sim.run(solver_config)` | integrates the fixed domain once |
| `sim.run_adaptive(solver_config, boundary_threshold, max_doublings, points_per_segment)` | integrates in segments of `points_per_segment` output points. After each segment it checks the fraction of the content at the top of the size grid. Above `boundary_threshold` it keeps the results up to the first exceedance, doubles `I`, `V` or both, maps the state onto the larger domain and resumes from there. The run never restarts from the beginning |

- `max_doublings = 0` keeps the domain fixed.
- With the loop → network loss on, the network dislocation density is advanced between segments. `run_adaptive` is therefore the method to use when it is on, even with a fixed domain. `n_points` and `points_per_segment` then set the number of updates.
- `solver_config['required_times']` lists times that must be in the output. Each becomes a segment boundary, so the state there is written exactly. The output times are logarithmic, so a dose to be reported is otherwise not on the grid.

Both return the `results` dictionary: [Outputs and provenance](Outputs-and-Provenance).

## Stopping a run

A run is stopped by asking the solver to stop. It checks for the request between output times, and every output time it has reached is already written, so the trajectory up to that point is kept. If the solver does not exit within 60 seconds it is killed; the rows already written are still read.

| How | What happens |
|---|---|
| Ctrl+C, or *Interrupt Kernel* | the solver is asked to stop; completed segments are kept and saved. A second Ctrl+C forces the stop |
| `timeout_s` is reached | the same graceful stop |
| the file named by the environment variable `RADCLUSTER_ABORT_FILE` appears | the same graceful stop, for a detached run that has no console |

- A run whose solver was interrupted is marked `_PARTIAL` in the name of its directory and in `provenance.md`.
- Under `run_adaptive`, a `timeout_s` given in `solver_config` applies to each segment; given as the argument `timeout_s=`, it is the budget of the whole run. A budget that runs out between two segments ends the run without the `_PARTIAL` mark, as does a CVODE failure: compare the last dose with the one requested.

## Threads

The solver chooses its number of OpenMP threads from the length of the right-hand-side sweep, the larger of the number of equations and `I + V`. `OMP_NUM_THREADS` overrides the choice; the Run cell removes that variable so that the solver chooses. The number used is written to `provenance.md`.

## Without the notebook

The notebook only calls the Python package. The same run from a script is in [Quick start](Quick-Start). The scripts in `radcluster_code/codes/Python_Testing/` are complete examples, for instance `gate_im5vm2.py`.

## Replotting

The Replot cell reads `plots/plot_data.pkl` of the latest run, or of the run named by `PLOT_DATA_PATH`. Each figure has its own limits and scales there. The file is a pickle of live Python objects: it may not load after the classes change, and it is not tracked by git.

`radcluster_code/codes/make_dose_figures.py` redraws the density and size figures against dose for any run directory, with the EUROFER-97 measurements overlaid.
