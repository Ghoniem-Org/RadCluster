# Tests and conservation checks

The checks are standalone scripts in `radcluster_code/codes/Python_Testing/`. They gate the physics, not only the code: most of them run the solver and test a conservation residual or an equivalence between two formulations. There is no pytest suite.

## Run a check

```bash
cd radcluster_code
PYTHONPATH=$PWD:$PWD/.. python codes/Python_Testing/check_he_conservation.py
```

- Most scripts find the module from their own location. A few insert a path of the machine they were written on; setting `PYTHONPATH` to `radcluster_code` and to the repository root, as above, makes them run anywhere.
- `check_binwise_moments.py`, `check_idiscrete_convergence.py` and `check_path_equivalence_20dpa.py` also import the campaign worker, so add `radcluster_code/digital_twin` to `PYTHONPATH` for them.
- `gate_im5vm2.py` takes the module folder as its argument: `python codes/Python_Testing/gate_im5vm2.py $PWD`.
- Checks that run the solver need it built: [Installation](Installation).

## The two residuals

Every run reports them, and the checks gate on them. See [Physics and equations](Physics-and-Equations).

| Residual | Tests | Expected |
|---|---|---|
| δ_FP | the swelling identity: vacancy content = SIA content + net flux to sinks | below about 10⁻⁶ |
| δ_He | helium in the system = helium produced − helium lost at sinks | below about 10⁻⁶ |

A value above 10⁻³ calls for an explanation. A change to a rate kernel or to stoichiometry comes with the relevant check and the two residuals it produces.

## Graph and declaration

| Script | Needs the solver | Checks |
|---|---|---|
| `check_eurofer_rag.py` | no | the EUROFER-97 graph builds and validates for both equation forms: every edge has a kernel of the right shape, the monomer populations are set, the layout is sensible. The container image runs it at build time |
| `check_coalescence_product_population.py` | no | on small graphs built for the test: coalescence conserves signed defect content, including when the product goes to another population; a product population of the wrong polarity is rejected |
| `check_loop_conversion_integration.py` | no | the EUROFER-97 graph has three populations and the three conversion edges; on a small graph of two SIA populations, all SIA edges together conserve the total SIA content through the graph walker |
| `test_two_layer_im5vm2.py` | yes | a simulation carries the live graph; it also compares the regression case with a baseline recorded for an earlier workbook |

## Conservation

| Script | Checks |
|---|---|
| `check_he_conservation.py` | δ_He and δ_FP for the four physics options with both helium kinetics, eight cases, each in its own process |
| `gate_im5vm2.py` | the regression case `im5vm2` (`I = V = 1000`, discrete, fission, 673 K, to 10⁻² s); prints the observables and both residuals for comparison with the baseline |
| `check_loop_conversion_conservation.py` | δ_FP with the loop conversion on |
| `check_loop_network_loss.py` | the loop → network channel removes loops and does not degrade δ_FP |

## Loop conversion

| Script | Checks |
|---|---|
| `check_loop_conversion_kernels.py` | the unary rate, the junction yield and the ⟨100⟩ kernels |
| `check_loop_conversion_params.py` | the conversion parameters load from the workbook, reach the rates and change the physics |
| `calibrate_loop_conversion.py` | fits the crossover temperature to the measured ½⟨111⟩ fraction |
| `sweep_loop100_and_network.py` | the ½⟨111⟩ fraction against temperature, and the network density against dose |

## Bin moments against the discrete equations

| Script | Question |
|---|---|
| `check_bm_vs_discrete.py` | with no bins and `i_discrete = I`, do the bin-moment equations equal the discrete ones? |
| `check_path_equivalence_20dpa.py` | the same at 20 dpa, where the distribution is broad |
| `check_binwise_moments.py` | bin by bin, are the moments carried by a binned run those of the exact distribution? No reconstruction is involved, so the comparison measures the dynamics and not the closure |
| `check_idiscrete_convergence.py` | does the residual error scale with `i_discrete` or with the number of bins? |
| `check_binmoment_loop100.py` | the bin-moment reduction of the ⟨100⟩ population, and the KLU solver |

The study these belong to is reported in `docs/Formulation/verification_campaign/`.

## Solver

| Script | Compares |
|---|---|
| `compare_linsol.py` | six combinations of equations, linear solver and preconditioner, with their timings (`compare_linsol_report.md`) |
| `benchmark_woodbury.py` | the Woodbury and Jacobi preconditioners on a short time span |
| `compare_window_modes.py` | `full_system` against `active_window` at `I = V = 10000` |
| `test_lapack_omp_small.py` | that LAPACK is linked and the thread choice works |
| `test_solver.py`, `diagnose_stall.py` | CVODE diagnostics for short runs and for a stalled bin-moment run |

## Parameter studies

`bias_study.py`, `zii_study.py`, `mobility_sweep.py`, `param_sweep.py`, `growth_diagnostic.py` and `calibrate.py` are studies, not gates: they vary bias factors, mobility cut-offs and other parameters and report how the size distributions respond.

## Results with the workbook of version 2.2.0

Run on Linux with SUNDIALS 7.1.1, without KLU, on 2 October 2026:

| Script | Result |
|---|---|
| `check_eurofer_rag.py` | all checks passed |
| `check_coalescence_product_population.py` | 5 of 5 passed |
| `check_loop_conversion_integration.py` | all passed |
| `check_he_conservation.py` | all 8 cases passed; the largest δ_He is 1.1 × 10⁻¹¹ |
| `check_loop_conversion_conservation.py` | all passed; final δ_FP 1.6 × 10⁻¹⁵ |
| `gate_im5vm2.py` | δ_FP = 2.9 × 10⁻⁷, δ_He = 2.0 × 10⁻¹⁴ |
| `check_loop_conversion_params.py` | 3 assertions fail: they expect `E_a0_conv`, `gamma_a_conv` and `phi_max_junc` at their values before the calibration |
| `check_loop_conversion_kernels.py` | 1 assertion fails: the junction yield below `n_j_min_junc`, whose effective value is now tied to `i_mobile` |
| `check_loop_network_loss.py` | 1 assertion fails: on its short test case the channel changes the loop inventory by less than the test requires |

The three failures are assertions written for earlier parameter values, not conservation failures. The first two scripts do not run the solver. In the third, δ_FP stays below 3 × 10⁻¹¹ with the channel on and off, and the script passes with the workbooks that preceded the calibration.

## Machine check for campaigns

`radcluster_code/digital_twin/check_machine.py` runs a fixed probe and compares twelve quantities with a committed reference at a relative tolerance of 10⁻⁶. A machine that fails it was built differently and must not contribute rows. See [Calibration and verification campaign](Calibration-and-Verification-Campaign).
