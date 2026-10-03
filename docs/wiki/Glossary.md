# Glossary

## Abbreviations

| Abbreviation | Meaning |
|---|---|
| CD | cluster dynamics |
| RAG | Reaction Admissibility Graph: the declared graph of a host material |
| SIA | self-interstitial atom |
| FP | Frenkel pair: one vacancy and one self-interstitial |
| RAFM | reduced-activation ferritic–martensitic (steel) |
| dpa | displacements per atom, the unit of dose |
| NRT | Norgett–Robinson–Torrens displacement rate, before in-cascade recombination |
| BDF | backward differentiation formula, the stiff integrator of CVODE |
| GMRES | generalized minimal residual method, the iterative linear solver |
| KLU | the sparse direct linear solver of SuiteSparse |
| SMW | Sherman–Morrison–Woodbury formula |
| QSS | quasi-steady state |
| TEM | transmission electron microscopy |
| RIS | radiation-induced segregation |
| LF, HF | low and high fidelity, in the calibration campaign |

## Terms

| Term | Meaning |
|---|---|
| polarity | χ = −1 for a vacancy-type cluster, +1 for an interstitial-type cluster |
| population | clusters of one polarity that share a mobility law and an energetics: `bulk-111`, `bulk-100`, `bulk` |
| cluster identifier | the four-tuple (polarity, size, population, composition) |
| signed defect content | q = χn; conserved by every reaction inside the host |
| vertex | an admissible cluster |
| edge | an admissible elementary reaction. In the code an `Edge` is a family: one reaction class applied along the whole size axis of a population |
| edge class | one of the ten abstract kinds of reaction |
| stoichiometric contract | what an edge class fixes: number of precursors, kinetic order, whether it crosses polarity, and its shifts of size, population and composition |
| kernel | the size-resolved rate constant of an edge, registered on the graph as a precomputed array |
| monomer population | the population whose size-1 vertex is the point-defect pool of its polarity |
| graph walker | the routine that assembles the right-hand side by visiting every edge |
| state layout | the map from the blocks of unknowns to positions in the state vector |
| host material | the material a graph is declared for; EUROFER-97 here |
| layer 1, layer 2 | the abstract core, and a material's declaration |
| process P1–P8 | the eight EUROFER-97 kernel families: recombination, cavity absorption, loop absorption, fixed sinks, thermal emission, glissile cluster–cavity reaction, trap mutation, radiation re-solution |
| discrete equations | one equation per cluster size (`full_CD` in option names) |
| bin moment | a moment of the size distribution over a logarithmic bin: μ⁽⁰⁾ number, μ⁽¹⁾ content, μ⁽²⁾ |
| closure, shape function | the assumed distribution inside a bin: `constant`, `linear`, `lognormal` |
| reconstruct, walk, project | the three steps of a bin-moment right-hand side |
| helium Case 1 | the mean-field reduction: one helium content per cavity size; selected by `cascade='fusion'` |
| helium Case 2 | the decoupled reduction: one scalar helium inventory; selected by `cascade='fission'` |
| spectrum | the workbook switch that selects the cascade source: the fission or the fusion column of the `Production` sheet |
| mobility cut-off | `i_mobile`, `v_mobile`: the largest mobile SIA and vacancy clusters |
| glissile, sessile | able, or unable, to glide. ½⟨111⟩ loops are glissile; ⟨100⟩ loops are sessile |
| fixed sink | a sink that is not tracked as a cluster: network dislocations, boundaries, precipitates |
| bias | the preference of a sink for SIAs over vacancies, `Z_i`, `Z_v` |
| network density | ρ_net, the density of network dislocations; evolved between segments when the loop → network loss is on |
| swelling | the vacancy content of the cavities, as a volume fraction |
| δ_FP, δ_He | the Frenkel-pair and helium conservation residuals |
| trap mutation | an over-pressurized bubble gains a vacancy by emitting an SIA |
| re-solution | a displacement event returns a trapped helium atom to the lattice |
| concentration floor | `C_floor`, the smallest concentration written |
| active window | the range of sizes advanced in `active_window` mode |
| segment | one solver call of `run_adaptive`, of `points_per_segment` output times |
| override | a workbook value replaced for one run |
| provenance | what a run records about its inputs, configuration and machine |
| reference run | the tracked run that published figures are checked against |
| digital twin | the calibration and verification campaign in `radcluster_code/digital_twin/` |
| design | the table of parameter vectors a campaign runs, generated once and committed |
| row | one run of a campaign: one parameter vector at one condition |
| usable row | a row that ran and reached its dose; the rule `merge_and_sobol.py` applies |
| extent-verified | an observable that keeps its value when the size grid is enlarged |
| lever | a parameter, in the calibration guide |

## Names in the code

| Name | Meaning |
|---|---|
| `I`, `V` | the largest SIA and vacancy cluster sizes (`N`, `M` in older code) |
| `i_cascade`, `v_cascade` | the largest clusters produced in a cascade (`m1`, `n1` in older code) |
| `im5vm2` | a case label: `i_mobile = 5`, `v_mobile = 2`. The regression case is `im5vm2` with `I = V = 1000` |
| `RadCluster_2_1` | the name of `radcluster_code/` before version 2.1.1 |
| `full_CD` | the discrete equations, in `physics_option` strings |
