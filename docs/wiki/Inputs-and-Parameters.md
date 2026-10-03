# Inputs and parameters

The physics parameters are in one workbook, `radcluster_code/input/input_parameters.xlsx`, which the code reads at the start of every run.

## The workbook

Each row has a parameter name, a `Symbol`, a value and its units; most rows also have a note with the equation, table or calibration the value comes from. The `Symbol` is the key the code and the overrides use. The `Production` sheet has two value columns, one per spectrum.

| Sheet | Contents | Examples of symbols |
|---|---|---|
| `Production` | the cascade source, for fission and for fusion: survival efficiency, clustered fractions, power-law slopes, largest in-cascade cluster sizes | `eta`, `f_cl_i`, `f_cl_v`, `s_i`, `s_v`, `m1`, `n1` (read as `i_cascade`, `v_cascade`) |
| `Energetics` | lattice and elastic constants, atomic volume, Burgers vectors, formation and migration energies, surface energy, attempt frequencies | `a`, `Omega`, `b_111`, `E_f_v`, `E_m_v`, `E_m_i`, `E_m_h`, `gamma_s`, `nu_v`, `nu_i` |
| `Diffusion` | one-dimensional glide; solute concentrations and trap binding energies; the SIA mobility cut-off | `E_m_1D`, `s_1D`, `L_hat`, `B_rot`, `c_Cr`, `c_C`, `E_b_C_SIA`, `i_mobile` |
| `Dissociation` | the binding-energy models: void capillarity, loop binding fits, helium binding, trap-mutation barriers | `lambda`, `A_111`, `B_111`, `E_b_i2`, `A_100`, `B_100`, `n_tr`, `alpha_He` |
| `Reactions` | bias factors, sink densities, the irradiation condition, the loop-conversion and network parameters, and switches | `Z_i`, `Z_v`, `rho_d`, `d_g`, `rho_p`, `r_p`, `T`, `G`, `spectrum`, `LOOP_COAL`, `LOOP_NETWORK_LOSS` |

The irradiation condition is in the `Reactions` sheet: temperature `T` (K), displacement rate `G` (dpa/s), helium production `G_He_r` (appm/dpa) and `spectrum`. `spectrum` (`fission` or `fusion`) selects the column of the `Production` sheet that the solver uses.

Units are given in each row: lengths in m (the lattice constant and Burgers vectors in nm), energies in eV, temperature in K (the conversion crossover `T_star_conv_C` in °C), time in s. Concentrations are atomic fractions; sizes are numbers of atoms.

**Not every row is read.** Some rows document a value that the code computes or fixes itself, and overriding them changes nothing:

| Row | What the code uses |
|---|---|
| `A_sph`, `A_loop`, `A_1D` | computed from their formulas in `input_data.py` |
| `D0_v`, `D0_i`, `D0_h` | the diffusivities are a² ν exp(−E_m/k_BT), from the attempt frequencies and migration energies |
| `B2`, `B3`, `mu_He` | constants in `binding_energies.py` |
| `b_100` | set to the lattice constant `a` |
| `E_b_v2` | not read |

To see whether a symbol is used, search `py_utils/` for it.

## How it is read

`py_utils.input_data.InputData` reads the five sheets into dictionaries and computes the derived quantities.

| Attribute | Contents |
|---|---|
| `production_fission`, `production_fusion` | the two columns of the `Production` sheet |
| `energetics`, `diffusion`, `dissociation`, `reactions` | one dictionary per sheet, keyed by `Symbol` |
| `derived` | computed quantities: `kBT`, effective diffusivities and jump frequencies (`Di_eff`, `omega_i_eff`, …), the solute trapping sums (`trap_SIA`, `trap_VAC`, `trap_loop`), the equilibrium vacancy concentration `Cv_eq`, `G`, `G_He`, the geometric prefactors |
| `I`, `V`, `i_mobile`, `v_mobile`, `i_discrete`, `v_discrete`, `I_bin`, `V_bin`, `shape_function` | the domain and the discrete/bin split |
| `solver_mode`, `physics_option`, `equations`, `cascade` | the selected mode and physics |

`py_utils.reaction_rates.ReactionRates` then computes the size-resolved rate constants from an `InputData`. These arrays are what the C++ solver receives through its parameter file, and most of what the reaction graph registers as kernels.

Symbols are renamed on reading: `N` → `I`, `M` → `V`, `n_max_i` → `i_mobile`, `m_max_v` → `v_mobile`, `m1` → `i_cascade`, `n1` → `v_cascade`.

## Overrides for one run

An override replaces a workbook value for one run. In the notebook it is an entry of `PARAM_OVERRIDES`, keyed by `Symbol`:

```python
PARAM_OVERRIDES = {'T': 673.0, 'rho_d': 5e14, 'eta': 0.28}    # K, m^-2, -
```

The notebook applies it as follows. A script does the same.

```python
inp = sim.input_data
for key, value in PARAM_OVERRIDES.items():
    placed = False
    for sheet in (inp.production_fission, inp.production_fusion, inp.diffusion,
                  inp.reactions, inp.energetics, inp.dissociation):
        if key in sheet:
            sheet[key] = value
            placed = True
    if not placed:
        inp.reactions[key] = value
inp._calculate_derived()
sim.rebuild_rates()
```

- A symbol present in several sheets is replaced in each. A symbol found nowhere is added to `reactions`.
- `_calculate_derived()` and `rebuild_rates()` are both needed. Without them the solver runs with the old rate constants.
- The values written to `provenance.md` are those after the overrides.
- Permanent changes belong in the workbook.

> An active override silently departs from the calibrated workbook. The notebook's `PARAM_OVERRIDES` is empty for that reason; its comments record what each lever does.

## What the constructor takes

`RadClusterSimulation(...)` takes only these:

| Argument | Meaning |
|---|---|
| `I`, `V` | the largest SIA and vacancy cluster sizes |
| `solver_mode` | `'full_system'` or `'active_window'` |
| `equations`, `cascade` | the two physics choices, given together |
| `physics_option` | the same choice as one string: `full_CD_fission`, `full_CD_fusion`, `bin_moment_CD_fission`, `bin_moment_CD_fusion`. Not together with `equations` and `cascade` |
| `excel_file` | another workbook |
| `C_floor` | the concentration floor |
| `he_kinetics` | `'dynamic'` or `'quasi_steady_state'` |
| `i_mobile`, `v_mobile` | the mobility cut-offs |

The earlier names `N`, `M`, `n_max_i` and `m_max_v` are also accepted. Anything else raises `TypeError`. In particular the bin-moment layout (`i_discrete`, `v_discrete`, `I_bin`, `V_bin`, `shape_function`) is set through `input_data.reactions`, as an override.

Read the layout back from `sim.rate_equations`, not from `sim.reaction_rates`: the first owns `i_discrete`, `I_bin`, `shape_function`, `n_mom` and `N_eq`.

`make_physics_option(equations, cascade)` and `split_physics_option(physics_option)` in `py_utils` convert between the two forms.

## Switches

These are keys of `input_data.reactions`. All but the last are rows of the `Reactions` sheet.

| Symbol | Meaning |
|---|---|
| `spectrum` | `fission` or `fusion`: the cascade source |
| `LOOP_COAL`, `loop_coal_pref` | loop–loop coalescence of ½⟨111⟩ loops by glide, and its prefactor |
| `LOOP_NETWORK_LOSS` | the loop → network-dislocation loss channel. It needs `run_adaptive` |
| `VOID_NETWORK_LOSS` | sweeping of cavities by the network |
| `Z_loop_model` | `0`: a constant loop bias; `1`: a size-dependent one |
| `sia_emission_monomer` | `returned` (the default when the key is absent, as it is in the workbook) or `dropped`. `dropped` reproduces results obtained before version 2.1.2, in which the monomer emitted by a loop was not returned to the monomer population |
| `boundary_flux` | `absorption` (the default when absent): a reaction product larger than the size grid is lost. `reflection`: the reaction is suppressed |

The loop-conversion parameters are described in [Loop character and network](Loop-Character-and-Network).

## Other files in `input/`

| File | Contents |
|---|---|
| `EuroferMicrostructure.xlsx` | the consolidated experimental microstructure workbook, built by `build_eurofer_microstructure.py` |
| `input_parameters.BACKUP-*.xlsx` | the workbook before each recalibration |
| `make_microdata_notebook.py` | generates a notebook of density and size maps from the microstructure workbook |

`py_utils/create_excel.py` writes a workbook with the manuscript's nominal parameter values to `py_utils/input/input_parameters.xlsx`. Its docstring says it overwrites `input/input_parameters.xlsx`; do not copy its output over the calibrated workbook.
