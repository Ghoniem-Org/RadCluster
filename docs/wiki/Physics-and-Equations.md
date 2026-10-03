# Physics and equations

This page summarizes what the code integrates, in the forms the code uses. The derivations are in the framework manuscript, and [radcluster_code/CLAUDE.md](https://github.com/Ghoniem/RadCluster/blob/main/radcluster_code/CLAUDE.md) lists the equations with their labels in that document. Where that file and the code differ, this page follows the code; the differences are listed at the end.

Units: concentrations are atomic fractions, sizes are numbers of atoms, energies are in eV, temperature in K and time in s. Rate constants carry s⁻¹.

## The master equation

For every tracked cluster species $a$,

```math
\frac{dc_a}{dt} = G_a + \sum_r S_{ar}\,J_r(c) - \mathcal{D}_a\,c_a
```

| Symbol | Meaning |
|---|---|
| $G_a$ | the cascade source |
| $J_r$ | the flux of elementary reaction $r$ |
| $S_{ar}$ | the stoichiometric coefficient of species $a$ in reaction $r$ |
| $\mathcal{D}_a$ | the loss rate at fixed sinks: network dislocations, grain boundaries, precipitates |

## Clusters

A cluster is identified by the four-tuple σ = (χ, n, p, c):

| Part | Meaning |
|---|---|
| χ | polarity: −1 for vacancy type, +1 for interstitial type. Trapped solute does not change it |
| n | size: the number of point defects of that polarity |
| p | population: clusters of one polarity that share a mobility law and an energetics, for example glissile ½⟨111⟩ loops, sessile ⟨100⟩ loops, bulk cavities |
| c | composition: the numbers of trapped solute or gas atoms |

The signed defect content of a cluster is q = χn.

## The ten reaction classes

Every edge of the graph belongs to exactly one class.

| # | Class | Action | Order |
|---|---|---|---|
| 1 | `GROWTH` | absorption of a monomer of the same polarity, n → n+1 | 2 |
| 2 | `SHRINKAGE` | impingement of a monomer of the opposite polarity, n → n−1 | 2 |
| 3 | `DISSOCIATION` | thermal emission of a monomer, n → n−1; the kernel contains a binding energy | 1 |
| 4 | `RECOMBINATION` | annihilation with a monomer of the opposite polarity | 2 |
| 5 | `ANNIHILATION` | two clusters of opposite polarity; the survivor has size \|n − n′\| | 2 |
| 6 | `INTER_POPULATION` | p → p′ at fixed size and composition | 1 |
| 7 | `SOLUTE_TRAPPING` | trapping or detrapping of one solute atom, c → c ± e_s | 2 or 1 |
| 8 | `COALESCENCE` | two clusters of the same polarity, n + n′; the product may go to another population of that polarity | 2 |
| 9 | `SOURCE` | cascade injection; the only class with no precursor | 0 |
| 10 | `SINK` | absorption at a fixed sink | 1 |

## Cascade source

The production rate after in-cascade recombination is $G = \eta\,G_{\rm NRT}$. Clusters are produced directly in cascades with a power law in size:

```math
\epsilon_m = C\,m^{-s}, \qquad C = \frac{f^{\rm cl}}{\sum_{m=2}^{m_1} m^{1-s}}, \qquad P_m = \epsilon_m\,G \quad (m \ge 2)
```

The monomers take the rest of the defects, so that the source conserves atoms exactly, including when the size grid is shorter than the largest cascade cluster:

```math
P_1 = G\Bigl(1 - \sum_{m\ge 2} m\,\epsilon_m\Bigr)
```

There is one set of parameters for SIAs and one for vacancies, and one column of the `Production` sheet for each spectrum: survival efficiency η, clustered fractions, slopes and largest in-cascade sizes. The workbook's `spectrum` selects the column.

## Transport

Point-defect diffusivities are $D_\alpha = a^2\,\nu_\alpha \exp(-E_m^\alpha/k_BT)$.

**Solute trapping.** In the alloy, dissolved solutes (Cr, W, Mn, C, N) slow the migration of SIAs and vacancies:

```math
D_\alpha^{\rm eff} = \frac{a^2\,\nu_\alpha \exp(-E_m^\alpha/k_BT)}{1 + \sum_s z_s\,c_s \exp(E_b^{\alpha s}/k_BT)}, \qquad \alpha = i, v
```

Helium is not trapped by solutes in the code.

**One-dimensional glide.** Glissile SIA clusters migrate along their Burgers vector, slowed by the solutes that trap loops:

```math
D_n^{\rm 1D} = \frac{3a^2\nu_0}{2\,n^{s}} \exp\Bigl(-\frac{E_m^{\rm 1D}}{k_BT}\Bigr)\,\frac{1}{1 + \sum_s z_s\,c_s \exp(E_b^{{\rm loop},s}/k_BT)}
```

with a mean free path $\hat L = L/a$ between changes of direction.

**Mobility cut-offs.** `i_mobile` and `v_mobile` are the largest mobile SIA and vacancy clusters. With both equal to 1, only monomers move. With `LOOP_COAL` on, ½⟨111⟩ loops larger than `i_mobile` still coalesce with one another by glide; see [Loop character and network](Loop-Character-and-Network).

## Reaction kernels

The EUROFER-97 kernel library has eight processes. A capture or emission constant is a volumetric rate divided by the atomic volume, so it has the form of a geometric factor times $D/\Omega^{2/3}$, in s⁻¹.

The geometric prefactors $A_{\rm sph} = (48\pi^2)^{1/3}$ and $A_{\rm loop} = 8\sqrt{\pi/\sqrt{3}}$ are computed from these formulas. $B_{\rm rot}$ is read from the workbook.

| Process | Kernel |
|---|---|
| P1: vacancy–SIA recombination | $K_{iv} = 4\sqrt{3}\,\pi\,(D_i^{\rm eff} + D_v^{\rm eff})/\Omega^{2/3}$ |
| P2: point-defect absorption by a spherical cavity of size m | $A_{\rm sph}\,m^{1/3}\,D_\alpha^{\rm eff}/\Omega^{2/3}$ |
| P3: SIA absorption by a loop of size n | $A_{\rm loop}\,n^{1/2}\,Z_i^{\rm loop}\,D_i^{\rm eff}/\Omega^{2/3}$ |
| P4: fixed sinks | dislocations $Z_\alpha\,\rho_d\,D_\alpha$; grain or lath boundaries $Z^{gb}_\alpha\,\pi^2 D_\alpha/d_g^2$; precipitates $4\pi\,Z^p_\alpha\,\rho_p\,r_p\,D_\alpha$ |
| P5: thermal emission of a vacancy from a cavity | $A_{\rm sph}\,(m-1)^{1/3}\,D_v^{\rm eff}\,e^{-E_b^v/k_BT}/\Omega^{2/3}$ |
| P5: thermal emission of an SIA from a ½⟨111⟩ cluster | $A_{\rm sph}\,(n-1)^{1/3}\,D_i^{\rm eff}\,e^{-E_b^i/k_BT}/\Omega^{2/3}$ |
| P5: thermal emission of an SIA from a ⟨100⟩ loop | $A_{\rm loop}\,(n-1)^{1/2}\,D_i^{\rm eff}\,e^{-E_b^{100}/k_BT}/\Omega^{2/3}$ |
| P6: glissile SIA cluster with a cavity | $A_{\rm sph}\,m^{1/3}\,D_n^{\rm 1D} / [\Omega^{2/3}\,(1 + B_{\rm rot}\,\hat L^2\,m^{-1/3})]$ |
| P7: trap mutation | $\nu_0 \exp(-E_{\rm TM}(m,\ell)/k_BT)$ |
| P8: radiation re-solution | $b_0\,\ell\,\dot\phi$ |

Clusters of sizes 1 to 3 diffuse three-dimensionally; larger mobile clusters glide, and react with cavities through P6; larger clusters still are sessile.

**Helium release.** The C++ kernels release trapped helium at one rate for all cavities, $\beta_{\rm He} = \nu_h \exp[-(E_b^{hV} + E_m^h)/k_BT]$ times the helium content, with $E_b^{hV}$ the binding of a helium atom to a vacancy (`E_b_hV_1`). P7 and P8 are declared in the graph and are not evaluated by the C++ kernels.

**Cluster–cluster reactions.** The kernels registered on the graph, which the Python graph walker uses, have one form for coalescence and for annihilation:

```math
K_{n,n'} = \frac{8\pi\,(\xi_n + \xi_{n'})\,(D_n + D_{n'})}{\Omega^{2/3}}, \qquad \xi_n = \Bigl(\frac{3n}{8\pi}\Bigr)^{1/3}
```

The C++ solver computes its own forms (`K_ii_coal`, `K_vv_coal`, `K_vi_coal` in `rate_kernels.cpp`): the geometry of the target (a loop from size 4, a sphere below), a bias factor for loops, and the diffusivity of the mobile partner. The two sets are not the same; the C++ forms are what a simulation integrates.

Helium–helium clustering in the lattice is not included: the binding of two interstitial helium atoms in bcc iron is negligible.

## Energetics

| Quantity | Model |
|---|---|
| vacancy binding to a void | capillarity, $E_b^v(m) = E_f^v - A_{\rm void}\,[m^{2/3} - (m-1)^{2/3}]$ with $A_{\rm void} = 4\pi\gamma_s r_0^2$, plus an atomistic correction that decays with size |
| vacancy binding to a bubble | the helium in a cavity raises the binding of its vacancies, so helium stabilizes bubbles against emission. The C++ kernels apply this as a correction to the emission rate that depends on the helium-to-vacancy ratio |
| helium pressure | a third-order virial equation of state |
| SIA binding to a loop | $E_b^{\rm loop}(n) = A\,n^{B}$ at small sizes, with a positive exponent, blended to the continuum limit; separate fits for ½⟨111⟩ and ⟨100⟩ loops |

## State-space reductions

Two independent reductions keep the system small enough for engineering doses.

### Helium loading of cavities

The unreduced state of a cavity is two-dimensional: m vacancies and ℓ helium atoms.

| `cascade` | Reduction | What is tracked | Equations |
|---|---|---|---|
| `fusion` | Case 1, mean field: helium equilibrates fast, so each cavity size has a mean loading | the cavity concentration and the helium content of each size | I + 2V + 1 |
| `fission` | Case 2, decoupled: helium is a weak perturbation | the cavity concentrations and one scalar inventory. For the pressure correction the kernels give a cavity of size m the loading $\ell(m) = \bar\ell\,m^{2/3}$, with $\bar\ell$ the inventory divided by the number of cavities | I + V + 2 |

With `he_kinetics='quasi_steady_state'` the free helium concentration is computed from $dc_h/dt = 0$ instead of being integrated, which removes one unknown. Helium migrates with a very small energy, so it equilibrates quickly.

### Logarithmic bin moments

Sizes up to `i_discrete` (and `v_discrete`) keep one equation each. Larger sizes are grouped into logarithmic bins $\mathcal{B}_k$, each carrying P moments:

```math
\mu_k^{(0)} = \sum_{n \in \mathcal{B}_k} c_n, \qquad \mu_k^{(1)} = \sum_{n \in \mathcal{B}_k} n\,c_n, \qquad \mu_k^{(2)} = \sum_{n \in \mathcal{B}_k} n^2 c_n
```

| `shape_function` | P | Closure inside a bin | Truncation error in the bin ratio r, per the manuscript |
|---|---|---|---|
| `constant` | 1 | piecewise constant | O((r−1)²) |
| `linear` | 2 | hat function | O((r−1)³) |
| `lognormal` | 3 | log-normal; falls back to linear for a near-monodisperse bin | O((r−1)⁴) |

The reduction is defined as three steps at every evaluation of the right-hand side (`py_utils.core.reductions.BinMomentReduction`):

1. **Reconstruct** the per-size distribution from the moments with the closure.
2. **Evaluate** the per-size rates.
3. **Project** the per-size rates onto the moments: $d\mu_k^{(p)}/dt = \sum_{n\in\mathcal{B}_k} n^p\,dc_n/dt$.

The C++ right-hand side `rhs_bin_moment` follows this scheme for the point-defect ladders. It evaluates two channels on the moments themselves: the coarsening of loops beyond `i_mobile`, where each bin is taken at its mean size, and the products of mobile-cluster coalescence that land in the bins.

Both the SIA and the vacancy populations are binned. The number of equations is

```math
N_{\rm eq} = i_{\rm discrete} + P\,I_{\rm bin} + v_{\rm discrete} + P\,V_{\rm bin} + n_{\rm He}
```

plus the cumulative fluxes used by the diagnostics, and a second SIA block for the ⟨100⟩ population when the loop split is on. With no bins and `i_discrete = I`, the discrete equations are recovered.

## Conservation diagnostics

Let $S_I = \sum_n n\,c_n$ be the SIA content and $S = \sum_m m\,c_m$ the vacancy content, which is the swelling. The solver integrates the cumulative losses of SIAs and of vacancies at the fixed sinks, $J_I$ and $J_V$. Mutual annihilation removes equal numbers of both, so

```math
S(t) - S_I(t) = J_I(t) - J_V(t)
```

and the Frenkel-pair residual is the violation of this identity:

```math
\delta_{\rm FP}(t) = \frac{\lvert S - S_I - J_I + J_V \rvert}{S + S_I + J_I + J_V}
```

The helium residual δ_He compares the helium in the system with what was produced minus what was lost at sinks.

| Value | Reading |
|---|---|
| below about 10⁻⁶ | every channel is balanced |
| above 10⁻³ | look for a cause |

A large δ_FP is not always a coding error. What happens at the top of the size grid matters:

- **Products that leave the grid.** With the default boundary treatment (`boundary_flux = absorption`), a reaction whose product is larger than `I` or `V` loses that product, and δ_FP grows. This is the usual cause in long loop-coarsening runs; it falls as the grid grows. With `reflection` such reactions are suppressed, and nothing is lost.
- **The ⟨100⟩ loops.** Their growth is halted at the top of the grid without losing atoms. A ⟨100⟩ distribution stacked against the top of the grid therefore leaves δ_FP small. Judge the adequacy of the grid for them by the content near the top, not by δ_FP.

Before version 2.1.2 the C++ kernels did not return the monomer emitted by a loop to the monomer population, which was the whole of the residual of the regression case. The key `sia_emission_monomer` keeps the earlier arithmetic available.

## Where the code and CLAUDE.md differ

`radcluster_code/CLAUDE.md` writes some equations differently from the code. The code is what runs.

| Item | CLAUDE.md | Code |
|---|---|---|
| solute trapping of helium | applies to helium | helium is not trapped (`input_data.py`) |
| emission of an SIA from a ½⟨111⟩ cluster | loop geometry, $A_{\rm loop}(n-1)^{1/2}$ | spherical geometry, $A_{\rm sph}(n-1)^{1/3}$ (`reaction_rates.py`) |
| kernels | written with the jump frequency ω = D/a² | computed as $D/\Omega^{2/3}$, a constant factor $a^2/\Omega^{2/3}$ larger |
| grain-boundary sink | $4\pi^2 a^2\omega/d_g^2$, no bias | $Z^{gb}\pi^2 D/d_g^2$ |
| helium shared among cavities (Case 2) | in proportion to $m^{1/3} c_m$ | loading proportional to $m^{2/3}$ (`rate_kernels.cpp`) |
| which populations are binned | SIAs only | SIAs and vacancies |
| $B_{\rm rot}$ | $(4/\pi)(8\pi/3)^{1/3} \approx 2.627$ | 2.627, read from the workbook; the expression evaluates to 2.586 |
| $A_{\rm sph}$, $A_{\rm 1D}$ | ≈ 7.818 and ≈ 2.632 | $(48\pi^2)^{1/3} = 7.796$; $9/(8\pi^{2/3}) = 0.524$, which is not used |
