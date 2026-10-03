# Outputs and provenance

## The run directory

A run with `save_output=True` writes one directory under `radcluster_code/output/`. Its name is built from the run itself:

```
<YYYYMMDD_HHMMSS>[_<tag>]_<solver mode>_<physics option>_I<I>V<V>_im<i_mobile>vm<v_mobile>[_PARTIAL]
```

for example

```
20260906_063055_full_system_bin_moment_CD_fission_I80000V20000_im5vm5
20261002_151035_smoke_full_system_bin_moment_CD_fission_I4000V2000_im5vm1
```

- The time is local time, so the names sort in the order of the runs.
- `<tag>` is `sim.run_tag`, when it is set before the run. Use it to tell apart runs that differ only in their bin layout.
- `_PARTIAL` marks a run that was interrupted or reached its time limit.

Runs are named by their start time to the second, so one run does not overwrite another. `output/` is not tracked by git, except the designated reference run.

```
output/<run>/
├── provenance.md       every input table with the overrides applied, the user selections,
│                       the solver configuration, the machine and the CVODE statistics
├── results_t.npy       the output times, in seconds
├── results_y.npy       the solution: one row per unknown, one column per output time
│                       (without the ⟨100⟩ loop block; see below)
├── summary.csv         one row: the end-of-run summary
├── diagnostics.txt     rate constants, the two conservation residuals, the end state
└── plots/              the figures (PNG), and plot_data.pkl for the Replot cell
```

Each file is written separately, so a failure in one does not cost the others. A Ctrl+C during the save lets the whole save finish and is raised afterwards; a second one aborts.

## provenance.md

| Section | Contents |
|---|---|
| (1) Material Data | the six input tables (`Production` for fission and for fusion, `Energetics`, `Diffusion`, `Dissociation`, `Reactions`) with the overrides applied, and the derived quantities |
| (2) User Selections | solver mode, physics option, `I`, `V`, the mobility cut-offs, the discrete/bin split, the shape function, helium kinetics, the concentration floor, temperature, dose rate, and the path of the workbook |
| (3) Solver Configuration | `solver_config`, with `solver_method` flattened |
| (4) Run Statistics | wall-clock time, number of output times, host name, platform, processor, Python version, memory, OpenMP threads used, whether the run completed, and the CVODE counters (`steps`, `nfe`, `nni`, `nli`, `npe`, `nps`, `ncfn`, `netf`, `nlsetup`) |

For a `run_adaptive` run, the statistics of section (4) (wall-clock time, number of output times, CVODE counters) are those of the last segment, not of the whole run.

`provenance.md` does not record the git commit of the code. The campaign records it for each of its rows; see [Calibration and verification campaign](Calibration-and-Verification-Campaign).

## summary.csv

One row, written by `post_process.summary_csv_row`, so the summaries of many runs concatenate into one table.

| Column | Meaning |
|---|---|
| `solver` | solver mode and physics option |
| `physics_option` | for example `bin_moment_CD_fission` |
| `T_K`, `G_dpa_s` | temperature and displacement rate |
| `t_end_s`, `dose_dpa` | the last output time and the dose there |
| `C_SIA_tot_m3`, `C_VAC_tot_m3`, `C_He_tot_m3` | total SIA, vacancy and helium contents, per m³ |
| `mean_n_i`, `mean_n_v` | mean SIA and vacancy cluster sizes |
| `swelling` | void swelling, as a fraction |
| `delta_FP`, `delta_He` | the Frenkel-pair and helium conservation residuals |

An array quantity is reported at its last finite value; a field is empty only when no finite value exists.

## The results dictionary

`run` and `run_adaptive` return a dictionary. Arrays have one entry per output time.

| Key | Contents |
|---|---|
| `t`, `dose` | time (s) and dose (dpa) |
| `y` | the solution, shape (number of unknowns, number of times). With the loop split on it does not contain the ⟨100⟩ block |
| `y_sia100`, `y_sia100_raw` | the ⟨100⟩ loop population: per size, shape (`I`, number of times), and as solved. Not written to `results_y.npy` |
| `C_i1`, `C_v1`, `C_He_free` | free SIA, vacancy and helium concentrations |
| `C_SIA_tot`, `C_VAC_tot`, `C_He_tot` | total contents |
| `mean_n_i`, `mean_n_v`, `mean_n_111`, `mean_n_100` | mean cluster sizes, in total and for each loop character |
| `N_loops`, `N_loops_111`, `N_loops_100`, `N_voids` | number densities of loops and of voids |
| `f_111_loop` | the fraction of the SIA content in ½⟨111⟩ loops |
| `rho_net` | the network dislocation density |
| `swelling` | void swelling |
| `delta_FP`, `delta_He` | the conservation residuals |
| `delta_FP_sia`, `delta_FP_vac` | the residuals of the SIA and vacancy balances taken separately, for diagnosis |
| `J_SIA_fixed`, `J_VAC_fixed`, `J_SIA_mutual`, `J_VAC_mutual`, `J_He_sink` | the cumulative fluxes the residuals are built from |
| `metadata` | solver statistics, threads used, whether the run is partial |

Units. The solution `y` and the cumulative fluxes are in atomic fractions. The concentrations and number densities (`C_*`, `N_*`) are per m³: atomic fractions divided by the atomic volume, `results['Omega']`.

## Figures

`visualization.save_all_plots` writes, among others:

| File | Shows |
|---|---|
| `point_defects.png`, `totals.png`, `he_content.png` | monomers, total contents and helium against dose |
| `swelling.png`, `network_density.png` | swelling and the network dislocation density against dose |
| `number_densities.png`, `mean_sizes.png` | loop and void densities and mean sizes, all sizes |
| `number_densities_tem.png`, `mean_sizes_tem.png` | the same for clusters large enough to be seen by TEM |
| `loop_fraction.png`, `loop100_dist_evolution_*.png` | the ½⟨111⟩ fraction and the ⟨100⟩ loop distribution, with the loop split on |
| `sia_dist_tem_size.png`, `sia_dist_tem_diameter.png`, `vac_dist_tem_size.png`, `vac_dist_tem_diameter.png` | size distributions as densities per bin |
| `sia_dist_diameter_smooth.png`, `vac_dist_diameter_smooth.png` | diameter distributions as smooth curves at fixed doses |
| `sia_bin_conc.png`, `vac_bin_conc.png`, `vac_bin_content.png` | concentration and content of each discrete size and bin against dose |
| `sia_balance.png`, `vac_balance.png` | the fractions of the defects produced that are in clusters, lost at fixed sinks, and annihilated. They should add up to one |

**Size distributions are densities per bin.** The distribution figures draw `dc/dn` against size, or `dc/dD` against diameter, as stair segments whose integral is the population in the displayed range. A bin-moment run knows the content of a bin, not the concentration at each integer size, so a density per bin is the faithful display and is what a TEM histogram measures.

## Reading a run back

```python
import numpy as np
t = np.load('output/<run>/results_t.npy')
y = np.load('output/<run>/results_y.npy')      # y[k, j]: unknown k at time t[j]
```

The order of the unknowns is described in [Solver](Solver). For derived quantities without rerunning, use the notebook's Replot cell or `visualization.load_plot_data`.

## The reference run

One finished run is tracked: `output/20260906_063055_full_system_bin_moment_CD_fission_I80000V20000_im5vm5/`, without its `plot_data.pkl`. The verification campaign reads it by name (`digital_twin/verification/runs.py` and `make_tables.py`). When it is superseded, change the exemption in `.gitignore` and those two files together.
