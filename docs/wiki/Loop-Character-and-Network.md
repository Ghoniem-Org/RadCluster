# Loop character and network

Body-centred cubic iron and ferritic–martensitic steels grow interstitial loops of two Burgers vectors: glissile **½⟨111⟩** loops, which glide one-dimensionally, and sessile **⟨100⟩** loops. Cascades produce almost only ½⟨111⟩ clusters, so the ⟨100⟩ loops must form during the evolution of the microstructure. They grow in with dose and dominate at high temperature.

This page describes three channels that act on the loops:

- the conversion of ½⟨111⟩ loops to ⟨100⟩ loops;
- the coalescence of ½⟨111⟩ loops by glide;
- the loss of loops to the dislocation network, with a network density that evolves.

References: [docs/Formulation/dislocation_loops](https://github.com/Ghoniem/RadCluster/tree/main/docs/Formulation/dislocation_loops) (`loop_111_to_100_conversion.tex`, `derivation_loop_conversion_kernels.tex`, `loop_network_loss.tex`), and the notes in `docs/design_notes/`.

## Two populations

With `SIA_LOOP_SPLIT = True` in the notebook, or `solver_config['loop_conversion'] = 1`, the solver integrates two SIA populations:

| Population | Loops | Mobile | Source |
|---|---|---|---|
| `bulk-111` | ½⟨111⟩, hexagonal on {110} | up to `i_mobile` | cascades |
| `bulk-100` | ⟨100⟩, square on {100}, from size `n_loop_min` | no | conversion only |

With the split off, one ½⟨111⟩ population is integrated, and the post-processing reports a ½⟨111⟩ fraction of one and no ⟨100⟩ loops. The graph declares both populations in either case.

## Three conversion channels

| Channel | Reaction | Edge | Role |
|---|---|---|---|
| unary (thermodynamic) | ½⟨111⟩ₙ → ⟨100⟩ₙ | `loop_111to100_unary`, inter-population | carries the temperature trend |
| junction (kinetic) | ½⟨111⟩ₙ + ½⟨111⟩ₙ′ → ⟨100⟩ₙ₊ₙ′ | `SIA111_junction`, coalescence | nucleates ⟨100⟩ loops from collisions of loops of comparable size |
| absorption (kinetic) | ⟨100⟩ₘ + ½⟨111⟩ₙ → ⟨100⟩ₘ₊ₙ | `SIA100_absorb`, coalescence | the main growth route of ⟨100⟩ loops |

Every channel conserves the total SIA content of the two populations: the junction and absorption products carry the summed sizes, and the unary transfer keeps the size. The Frenkel-pair residual is therefore unaffected.

### Unary transformation

A single loop reorients at the rate

```math
\Gamma_{\rm uni}(n,T) = \nu_0\,\exp\Bigl(-\frac{E_a(n)}{k_BT}\Bigr)\,\max\Bigl[0,\;1 - \exp\Bigl(-\frac{\Delta F(n,T)}{k_BT}\Bigr)\Bigr], \qquad E_a(n) = E_a^0 + \gamma_a\,\frac{P_{111}(n)}{b_{111}}
```

- ΔF(n, T) is the free-energy difference between a ½⟨111⟩ and a ⟨100⟩ loop of n interstitials, computed by `py_utils.loop_energetics.LoopEnergetics`. The rate is zero where ΔF ≤ 0.
- The ⟨100⟩ loop energy softens with temperature as (1 − T/T_c)^¼, with T_c = 1185 K, so the conversion turns on above a crossover temperature.
- The magnitude of the ⟨100⟩ energy is not taken from elastic constants. It is derived from one target: the temperature `T_star_conv_C` at which ΔF = 0 for a loop of size `n_ref_conv`.
- P₁₁₁(n) is the loop perimeter, so the barrier grows with size.

### Junction

Two mobile ½⟨111⟩ loops that collide form a ⟨100⟩ loop only when their sizes are comparable. The yield is

```math
\varphi(n,n') = P_{\rm success}(T)\,\varphi_{\max}\,\exp\Bigl[-\frac{(\ln(n/n'))^2}{2\sigma_s^2}\Bigr] \quad \text{for } \min(n,n') \ge n_{j,\min}
```

and zero otherwise. It splits the ½⟨111⟩ collision kernel: the fraction φ forms a ⟨100⟩ loop, and the rest is ordinary coalescence, so the total collision rate is unchanged. The smallest junction size is tied to `i_mobile`, because a junction partner must be mobile.

P_success is the probability that the intermediate configuration goes on to ⟨100⟩ rather than back to ½⟨111⟩; it is set by two barriers, `dH2_conv` and `dH_rev_conv`.

### Absorption

A sessile ⟨100⟩ loop captures mobile ½⟨111⟩ clusters. The kernel is the collision kernel with the diffusivity of the ½⟨111⟩ partner alone, times its own success probability (barrier `dH2_abs_conv`, prefactor `psucc_abs_pref`), which is separate from the junction's.

### Parameters

These are symbols of the `Reactions` sheet (the last row is in `Dissociation`) and can be overridden. Their values are the calibrated ones in the workbook.

| Symbol | Meaning |
|---|---|
| `n_loop_min` | the smallest ⟨100⟩ loop tracked |
| `E_a0_conv`, `gamma_a_conv`, `nu0_conv` | the unary barrier, its slope with perimeter, and the attempt frequency |
| `T_star_conv_C`, `n_ref_conv` | the crossover temperature (°C) and the reference size that calibrate ΔF |
| `phi_max_junc`, `sigma_s_junc`, `n_j_min_junc` | the peak junction yield, its tolerance in log size, and the smallest junction size |
| `dH2_conv`, `dH_rev_conv` | the two barriers of P_success for the junction |
| `dH2_abs_conv`, `psucc_abs_pref` | the barrier and prefactor of P_success for absorption |
| `A_100`, `B_100` | the binding-energy fit of ⟨100⟩ loops, for their thermal emission |

### What the results contain

| Key | Meaning |
|---|---|
| `f_111_loop` | the fraction of the SIA content in ½⟨111⟩ loops |
| `N_loops_111`, `N_loops_100` | the number densities of each character |
| `mean_n_111`, `mean_n_100` | the mean sizes |

Figures: `loop_fraction.png` and `loop100_dist_evolution_*.png`.

## Coalescence of ½⟨111⟩ loops

Clusters up to `i_mobile` always coalesce. Without more, a loop that grows past `i_mobile` never coalesces again, and the mean loop size stays near the cut-off. `LOOP_COAL = 1` continues the one-dimensional glide law beyond `i_mobile` for the coalescence of ½⟨111⟩ loops with one another, and for that channel only: larger loops still do not reach fixed sinks or cavities. `loop_coal_pref` multiplies the diffusivity of this channel.

The calibration notes in the notebook record it as the lever on the mean ½⟨111⟩ loop size, and not a free one: the SIA content is conserved into loops, so smaller loops are more numerous.

## Loss of loops to the network

Without this channel the network dislocation density is a fixed sink, and large sessile loops accumulate without bound. With `LOOP_NETWORK_LOSS = 1`:

- **Loss.** A loop of size n is swept into the network at the frequency Λₙ = v_net ρ_net w_c P(n). The network climb velocity v_net follows from the net flux of SIAs over vacancies to dislocations; w_c is a capture width; P(n) is a geometric factor that rises from 0 to 1 as the loop reaches the spacing of the network. The frequency does not depend on the loop's own diffusivity, so it acts on sessile loops.
- **Network density.** ρ_net grows by the line length of the loops it incorporates and recovers:

```math
\frac{d\rho_{\rm net}}{dt} = \frac{\pi}{\Omega}\sum_n d_n\,\Lambda_n\,c_n - K_{\rm rec}\,\rho_{\rm net}^{3/2}
```

- **Operator splitting.** ρ_net is held constant within an integration segment and advanced between segments by one explicit step, after which the rates are rebuilt. The Jacobian keeps its banded plus low-rank structure, and the Woodbury preconditioner needs no change. For this reason the channel needs `run_adaptive`; `n_points` and `points_per_segment` set the number of updates.
- **Conservation.** The n interstitials of an absorbed loop are added to the cumulative fixed-sink flux, so δ_FP is preserved.
- **Bounds.** ρ_net stays between `rho_d` and `loop_net_rho_max`.

| Symbol | Meaning |
|---|---|
| `LOOP_NETWORK_LOSS` | 0: off, the behaviour of version 2.0; 1: on |
| `loop_net_chi` | the geometric range of the interaction |
| `loop_net_w_c` | the capture width; blank for the default |
| `loop_net_K_rec` | the recovery coefficient |
| `loop_net_rho_max` | the ceiling of ρ_net |
| `loop_net_n_inc`, `loop_net_xi` | the smallest loop size swept, and a floor at small sizes |

The result `rho_net` and the figure `network_density.png` show the network density against dose.

`VOID_NETWORK_LOSS` enables the corresponding sweeping of cavities by the network (`void_net_chi`, `void_net_w_c`, `void_net_m_inc`).

## Checks

| Script | Checks |
|---|---|
| `check_loop_conversion_kernels.py` | the unary rate, the junction yield and the ⟨100⟩ kernels |
| `check_loop_conversion_params.py` | that the workbook parameters reach the rates |
| `check_coalescence_product_population.py` | that coalescence into another population conserves content |
| `check_loop_conversion_integration.py` | all SIA edges together, through the graph walker |
| `check_loop_conversion_conservation.py` | δ_FP with the conversion on, through the C++ solver |
| `check_loop_network_loss.py` | that the network channel removes loops and leaves δ_FP unchanged |
| `sweep_loop100_and_network.py` | the ½⟨111⟩ fraction against temperature, and ρ_net against dose |

All are in `radcluster_code/codes/Python_Testing/`.

## In development

Radiation-induced segregation and the precipitation of solute-rich phases (Cr-rich α′, Mn–Ni–Si) are planned as new populations, edges and kernels of the EUROFER-97 graph, with the solver unchanged. The plan is `docs/Formulation/segregation_and_precipitation/radcluster_2_1_RIS_plan.tex`.
