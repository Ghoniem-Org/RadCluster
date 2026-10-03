# Adding a host material

A host material is a declaration: its populations, its admissible reactions and its rate kernels. The abstract Python core, the reductions and the preconditioners are reused unchanged; the C++ side needs the material's kernel file and its parameters. EUROFER-97 is the worked example: `py_utils/materials/eurofer97/declaration.py`.

## Where it goes

| Side | Location | Contents |
|---|---|---|
| Python | `radcluster_code/py_utils/materials/<name>/` | `__init__.py` and `declaration.py`, with one function `build_<name>_rag(input_data, reaction_rates, *, equations, cascade)` that returns the graph and its state layout |
| C++ | `radcluster_code/cpp_utils/materials/<name>/` | `rate_kernels.cpp`, with the right-hand sides declared in `core/rate_equations.h` |
| parameters | a workbook with the same five sheets | every parameter with its symbol, unit and source |

## What it must have

| Item | Location |
|---|---|
| the populations and the monomer population of each polarity | `declaration.py` |
| one edge family per admissible reaction | `declaration.py` |
| a registered kernel for every edge | `declaration.py`, from the arrays of the rates object |
| the state layout | `declaration.py` |
| the C++ kernels | `cpp_utils/materials/<name>/rate_kernels.cpp`, added to `CMakeLists.txt`. The right-hand sides are written by hand for the material; they are not generated from the graph |
| a check that the graph builds and validates | `codes/Python_Testing/`, after `check_eurofer_rag.py` |
| the conservation residuals of a regression case | recorded with the change |

## Declare the graph

A minimal host with two populations and seven reactions. It runs as written.

```python
import numpy as np
from radcluster_code.py_utils.core import (Polarity, Population, EdgeClass, Edge,
                                           ReactionAdmissibilityGraph, StateLayout)

def build_myhost_rag(I=50, V=50, i_mobile=1, v_mobile=1):
    rag = ReactionAdmissibilityGraph('MyHost')

    # 1. populations, and the monomer pool of each polarity
    sia = rag.add_population(Population('bulk', Polarity.SIA, n_min=1, mobile_max=i_mobile))
    vac = rag.add_population(Population('bulk', Polarity.VACANCY, n_min=1, mobile_max=v_mobile))
    rag.set_monomer_population(sia)
    rag.set_monomer_population(vac)

    # 2. kernels: arrays indexed by size - 1, or scalars
    n = np.arange(1, I + 1, dtype=float)
    m = np.arange(1, V + 1, dtype=float)
    source = np.zeros(I); source[0] = 1e-7
    rag.register_kernel('K_SIA_grow', 1e3 * n ** 0.5)
    rag.register_kernel('K_SIA_shrink', 1e1 * n ** 0.5)
    rag.register_kernel('K_VAC_grow', 1e1 * m ** (1 / 3))
    rag.register_kernel('K_iv', 5e2)
    rag.register_kernel('G_SIA', source)
    rag.register_kernel('G_VAC', source.copy())
    rag.register_kernel('D_SIA', np.where(n <= i_mobile, 1e-2, 0.0))

    # 3. edges: one family per reaction class and population
    rag.add_edge(Edge(EdgeClass.GROWTH, 'loop_growth', sia, kernel='K_SIA_grow'))
    rag.add_edge(Edge(EdgeClass.SHRINKAGE, 'loop_shrink', sia, kernel='K_SIA_shrink'))
    rag.add_edge(Edge(EdgeClass.GROWTH, 'void_growth', vac, kernel='K_VAC_grow'))
    rag.add_edge(Edge(EdgeClass.RECOMBINATION, 'recombination', vac, kernel='K_iv'))
    rag.add_edge(Edge(EdgeClass.SOURCE, 'cascade_SIA', sia, kernel='G_SIA'))
    rag.add_edge(Edge(EdgeClass.SOURCE, 'cascade_VAC', vac, kernel='G_VAC'))
    rag.add_edge(Edge(EdgeClass.SINK, 'sink_SIA', sia, kernel='D_SIA'))

    # 4. layout: one block per population
    layout = StateLayout()
    layout.add_discrete('SIA', I, population=sia)
    layout.add_discrete('VAC', V, population=vac)
    layout.freeze()

    rag.validate()
    return rag, layout
```

```
RAG 'MyHost': 2 populations, 7 edges, 7 kernels; edge classes active: {'growth': 2,
'shrinkage': 1, 'recombination': 1, 'source': 2, 'sink': 1}
```

In a real declaration the kernels are the arrays of a rates object computed from the workbook, as in `build_eurofer_rag(input_data, reaction_rates, *, equations, cascade)`.

- Population names are unique within a polarity; a vacancy `bulk` and an SIA `bulk` can coexist.
- A two-dimensional kernel for coalescence or annihilation is K[n−1, n′−1]. Register it as a function without arguments, so that it is built only if the graph walker asks for it.
- `Edge` checks the contract of its class when it is created: an inter-population edge needs a product population of the same polarity, an annihilation partner must be of the opposite polarity, a solute-trapping edge needs a gas species.
- `add_edge` rejects a duplicate label and a population that is not registered. `validate` checks that every edge has its kernel.
- For the bin-moment equations, use `layout.add_bin_moment(name, n_discrete, n_bins, moments_per_bin, population=..., reduction=BinMomentReduction(...))`.

## Choosing the edge class

| The reaction | Class |
|---|---|
| a cluster absorbs a monomer of its own polarity | `GROWTH` |
| a cluster is hit by a monomer of the opposite polarity | `SHRINKAGE` |
| a cluster emits a monomer thermally | `DISSOCIATION` |
| a monomer annihilates with one of the opposite polarity | `RECOMBINATION` |
| two clusters of opposite polarity react | `ANNIHILATION`, with `partner_population` |
| two clusters of the same polarity merge | `COALESCENCE`; `product_population` sends the product to another population |
| a cluster changes population at fixed size | `INTER_POPULATION`, with `product_population` |
| a cluster traps or releases a solute atom | `SOLUTE_TRAPPING`, with `gas_species` |
| a cascade creates a cluster | `SOURCE` |
| a cluster is lost at a fixed sink | `SINK` |

A reaction that needs two classes, such as trap mutation, is declared as the class of its main effect, with the second effect recorded in the edge's `meta`.

## Check it in Python

```python
from radcluster_code.py_utils.core import GraphWalker, structural_report

rag, layout = build_myhost_rag()
walker = GraphWalker(rag, layout, boundary='reflection')
y = np.random.default_rng(0).random(layout.N_eq) * 1e-9
dydt = walker.assemble(0.0, y)            # the right-hand side, by walking the edges

print(structural_report(rag, n_max=20))   # connectivity; vertices no reaction reaches
```

- `boundary='reflection'` suppresses a reaction whose product would exceed the size axis; `'absorption'` lets it occur and drops the product.
- **Conservation.** Walk the graph without its source and sink edges, and sum q = χn over the rates: Σ n·dc/dt over the SIA block minus the same over the vacancy block. It must vanish to round-off. For the host above it is 10⁻²⁹ against terms of 10⁻¹². `check_coalescence_product_population.py` and `check_loop_conversion_integration.py` do this for loop edges, on small graphs.
- **Structure.** The report must list no vertex that is unreachable from the monomers and none without an outgoing arc.
- **Solute trapping.** The walker has no composition axis, and raises `NotImplementedError` on a `SOLUTE_TRAPPING` edge. In EUROFER-97 helium is carried by the helium reduction, so the complete EUROFER-97 graph, which declares the re-solution edge, cannot be walked as declared. The checks walk small graphs of SIA populations built for the purpose.

## Mirror it in C++

- The right-hand side has the CVODE signature `int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)`; `user_data` is a `UserData` holding the `Parameters`.
- The Python side passes every precomputed rate through the parameter file; the right-hand side does arithmetic only.
- The preconditioners call the right-hand side through `UserData.rhs_fn` and need no change. `solver.cpp` selects the right-hand side by name, `parameters.h` holds the fields the kernels read, and the sparse Jacobian for KLU assumes the discrete EUROFER-97 layout; a second material extends these.
- The Woodbury preconditioner assumes a banded Jacobian plus a low-rank coupling to the mobile species. A model that keeps that structure reuses it as it is.

## Rules

- **Python first.** Develop the physics as edges and kernels for the graph walker; pass the conservation check; mirror it in C++; then calibrate.
- **Do not touch the core for a material.** If a reaction fits none of the ten classes, that is a change to the framework, with its own contract in `EDGE_CLASS_SPEC`.
- **Parameters go in the workbook**, with a symbol, a unit and a source.
- **A new channel has a switch**, so that earlier results stay reproducible: `LOOP_NETWORK_LOSS = 0` restores the fixed network, and `sia_emission_monomer = 'dropped'` the arithmetic before version 2.1.2.
- **A change to a kernel or to stoichiometry** comes with the δ_FP and δ_He it produces.

## Before committing

```bash
cd radcluster_code
PYTHONPATH=$PWD python codes/Python_Testing/check_eurofer_rag.py
python codes/Python_Testing/check_he_conservation.py
python codes/Python_Testing/gate_im5vm2.py $PWD
(cd .. && doxygen docs/doxygen/Doxyfile)      # the code reference still builds without warnings
```
