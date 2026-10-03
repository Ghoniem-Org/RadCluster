# EUROFER-97 reaction graph

EUROFER-97 is declared in `radcluster_code/py_utils/materials/eurofer97/declaration.py` by one function:

```python
rag, layout = build_eurofer_rag(input_data, reaction_rates, equations='discrete', cascade='fission')
```

It returns the `ReactionAdmissibilityGraph` and its `StateLayout`. It is a declaration, not an integrator: it says what EUROFER-97 is, and nothing about how to integrate it. Helium is the one resolved gas.

```
RAG 'EUROFER-97': 3 populations, 24 edges, 24 kernels
```

## Populations

| Polarity | Name | Smallest size | Largest mobile size | Contents |
|---|---|---|---|---|
| SIA | `bulk-111` | 1 | `i_mobile` | glissile ½⟨111⟩ loops and the small three-dimensionally mobile clusters. It holds the SIA monomers and receives all SIAs produced in cascades |
| SIA | `bulk-100` | `n_loop_min` (4) | 0 | sessile ⟨100⟩ loops. No glide, no coalescence among themselves, no cascade source |
| vacancy | `bulk` | 1 | `v_mobile` | voids and helium bubbles. It holds the vacancy monomers |

The size-1 vertex of `bulk-111` is the SIA monomer pool, and that of `bulk` the vacancy monomer pool. Growth, shrinkage, recombination and emission edges of any population draw from, or return to, the pool of the right polarity.

## Edges

Each edge is a family: one class of reaction along the whole size axis of a population.

| Label | Class | Population | Kernel | Reaction |
|---|---|---|---|---|
| `P1_recombination` | recombination | `bulk-111` | `K_iv` | Iₙ + V₁ → Iₙ₋₁ |
| `P2v_cavity_growth` | growth | `bulk` | `K_VAC_grow` | Vₘ + V₁ → Vₘ₊₁ |
| `P2i_cavity_shrink` | shrinkage | `bulk` | `K_VAC_shrink` | Vₘ + I₁ → Vₘ₋₁ |
| `P3_loop_growth` | growth | `bulk-111` | `K_SIA_grow` | Iₙ + I₁ → Iₙ₊₁ |
| `P3_loop_shrink` | shrinkage | `bulk-111` | `K_SIA_shrink` | Iₙ + V₁ → Iₙ₋₁ |
| `P4_SIA_sink` | sink | `bulk-111` | `D_SIA_sink` | loss of mobile SIA clusters at fixed sinks |
| `P4_VAC_sink` | sink | `bulk` | `D_VAC_sink` | loss of vacancy clusters at fixed sinks. The kernel is one scalar; the C++ applies it to mobile sizes only |
| `P5i_SIA_emission` | dissociation | `bulk-111` | `eps_SIA_emit` | Iₙ → Iₙ₋₁ + I₁ |
| `P5v_VAC_emission` | dissociation | `bulk` | `eps_VAC_emit` | Vₘ → Vₘ₋₁ + V₁ |
| `SIA111_self_coalescence` | coalescence | `bulk-111` | `K_111_self` | Iₙ + Iₙ′ → Iₙ₊ₙ′, the part that stays ½⟨111⟩ |
| `VAC_VAC_coalescence` | coalescence | `bulk` | `K_vv_coal` | Vₘ + Vₘ′ → Vₘ₊ₘ′ |
| `SIA111_junction` | coalescence, product in `bulk-100` | `bulk-111` | `K_111_junction` | ½⟨111⟩ₙ + ½⟨111⟩ₙ′ → ⟨100⟩ₙ₊ₙ′ |
| `loop_111to100_unary` | inter-population | `bulk-111` → `bulk-100` | `Gamma_uni` | ½⟨111⟩ₙ → ⟨100⟩ₙ |
| `VI_annihilation` | annihilation | `bulk` with `bulk-111` | `K_vi_annih` | Vₘ + Iₙ → the survivor, of size \|m − n\| |
| `cascade_SIA_source` | source | `bulk-111` | `G_SIA_cascade` | cascade injection of SIA clusters |
| `cascade_VAC_source` | source | `bulk` | `G_VAC_cascade` | cascade injection of vacancy clusters |
| `SIA100_absorb` | coalescence | `bulk-100` with `bulk-111` | `K_100_absorb` | ⟨100⟩ₘ + ½⟨111⟩ₙ → ⟨100⟩ₘ₊ₙ |
| `P3_100_growth` | growth | `bulk-100` | `K_100_grow` | ⟨100⟩ₙ + I₁ → ⟨100⟩ₙ₊₁ |
| `P3_100_shrink` | shrinkage | `bulk-100` | `K_100_shrink` | ⟨100⟩ₙ + V₁ → ⟨100⟩ₙ₋₁ |
| `P5_100_emission` | dissociation | `bulk-100` | `eps_100_emit` | ⟨100⟩ₙ → ⟨100⟩ₙ₋₁ + I₁ |
| `P4_100_sink` | sink | `bulk-100` | `D_100_sink` | sessile loops do not diffuse to fixed sinks; the kernel holds the loss to the network when `LOOP_NETWORK_LOSS` is on, and is zero otherwise |
| `VI_annihilation_100` | annihilation | `bulk` with `bulk-100` | `K_vi_annih_100` | Vₘ + ⟨100⟩ₙ → the survivor |
| `P7_trap_mutation` | growth | `bulk` | `Gamma_TM` | Vₘ(He) → Vₘ₊₁ + I₁ |
| `P8_radiation_resolution` | solute trapping (detrapping) | `bulk` | `Gamma_res` | Heₗ Vₘ → Heₗ₋₁ Vₘ + He |

By class: 4 growth, 3 shrinkage, 3 dissociation, 1 recombination, 2 annihilation, 1 inter-population, 1 solute trapping, 4 coalescence, 2 source, 3 sink.

Three edges couple the two loop populations: the junction, the unary transformation and the absorption. See [Loop character and network](Loop-Character-and-Network).

**P7 and P8.** Trap mutation has two effects, a larger cavity and one new SIA monomer; no single class describes both, so it is declared as a growth edge on the vacancy ladder with the SIA source recorded in the edge's `meta`. Both edges have their kernels registered. The declaration records that the C++ kernels do not evaluate them yet. The Python graph walker has no composition axis and raises an error on the re-solution edge, so the complete graph cannot be walked as declared. The checks that use the walker run small graphs of SIA populations built for the purpose.

## Kernels

The abstract core never computes a rate. Most kernels are arrays or scalars already computed by `ReactionRates`, registered under a name; the rest are assembled in the declaration from precomputed pieces.

- One-dimensional kernels are indexed by size − 1.
- The two-dimensional kernels of coalescence and annihilation, K[n−1, n′−1], are assembled in the declaration from the effective diffusivities, and registered as functions that build the array on first use. Only the Python graph walker asks for them; the C++ solver computes its own cluster–cluster rates, so a production run never allocates them.
- The junction yield splits the ½⟨111⟩ collision kernel in two: `K_111_junction` is the fraction φ that forms a ⟨100⟩ loop, and `K_111_self` the remainder.
- `Gamma_TM` and `Gamma_res` are built in the declaration from `binding_energies`.
- The cascade sources `G_SIA_cascade` and `G_VAC_cascade` are computed in the declaration by `production_rates` for the `cascade` argument, with that module's default spectrum parameters. The solver's source uses the workbook's `Production` sheet for the workbook's `spectrum`. The two differ when the workbook departs from the module defaults, as the calibrated workbook does.

## State layout

These are the layouts printed by `check_eurofer_rag.py` (`I = V = 1000`, `i_mobile = 5`, `v_mobile = 2`, fission). For the discrete equations:

```
StateLayout: N_eq = 3007
  [     0:  1000] SIA                      discrete   len=1000
  [  1000:  2000] SIA100                   discrete   len=1000
  [  2000:  3000] VAC                      discrete   len=1000
  [  3000:  3002] He                       aux        len=2
  [  3002:  3007] conservation             aux        len=5
```

For the bin-moment equations each size block is a discrete prefix followed by the moments of each bin, and the three populations share the same reduction:

```
StateLayout: N_eq = 69
  [     0:    21] SIA                      bin_moment len=21
  [    21:    42] SIA100                   bin_moment len=21
  [    42:    62] VAC                      bin_moment len=20
  [    62:    64] He                       aux        len=2
  [    64:    69] conservation             aux        len=5
```

- `He` holds the helium unknowns of the reduction selected by the cascade: one scalar inventory for fission, one content per cavity class for fusion, plus free helium.
- `conservation` holds the five cumulative integrals used by the conservation residuals.

This layout describes the declared graph. It always contains the `SIA100` block and a free-helium entry, whatever `loop_conversion` and `he_kinetics` are, so its `N_eq` is not the number of equations the solver integrates. The solver's vector is given in [Solver](Solver).

## Inspecting the graph

Every simulation carries the graph:

```python
rag, layout = sim.material_rag, sim.state_layout
print(rag.summary())
print(layout.describe())
for p in rag.populations:
    print(p.polarity.label, p.name, p.n_min, p.mobile_max)
for e in rag.edges:
    print(e.label, e.edge_class.value, e.population.name, e.kernel, e.meta.get('note', ''))
rag.validate()                      # every edge has a registered kernel
k = rag.kernel('K_SIA_grow')        # a numpy array, indexed by size - 1
```

Without a simulation object:

```python
from radcluster_code.py_utils.input_data import InputData
from radcluster_code.py_utils.reaction_rates import ReactionRates
from radcluster_code.py_utils.materials.eurofer97 import build_eurofer_rag

inp = InputData(I=200, V=200, physics_option='full_CD_fission')
rag, layout = build_eurofer_rag(inp, ReactionRates(inp), equations='discrete', cascade='fission')
```

## The graph as an explicit network

In a run the graph stays implicit: an edge is a family, and no vertex is ever created. For structural questions, `py_utils/core/to_networkx.py` writes it out as a NetworkX `MultiDiGraph`, truncated at a largest size. NetworkX is needed only for this.

```python
from radcluster_code.py_utils.core import to_networkx, structural_report, write_graphml

G = to_networkx(rag, n_max=20)           # 59 vertices, 3118 arcs
print(structural_report(rag, n_max=20))
write_graphml(G, 'eurofer_rag_n20.graphml')
```

The report gives the numbers of vertices and arcs, the arcs of each class, whether the graph is connected, the vertices that cannot be reached from the monomers (a population with no production channel is a declaration error) and the vertices with no outgoing arc. For EUROFER-97 both sets are empty.

Almost all arcs are coalescence and annihilation, which are hyperedges with two reactants: 2784 of the 3118 at a largest size of 20.

## Figures

| Script | Draws |
|---|---|
| `radcluster_code/codes/make_rag_figure.py` | the graph as three lanes, one per population, with size along the horizontal axis and one arc per edge family; `--n-max` sets the largest size |
| `radcluster_code/codes/make_rag_hyperedge_figure.py` | the binary reactions as tables in the (n, m) plane |

Both are generated from the materialized graph, not from a transcribed list of edges, so they cannot drift from the declaration. Their output, GraphML exports and a structural report are in `docs/Formulation/reaction_admissibility_graph/figures/`.
