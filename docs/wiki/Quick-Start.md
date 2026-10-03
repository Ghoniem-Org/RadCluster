# Quick start

This page runs one small case from Python. It takes a few seconds with the solver built ([Installation](Installation)). For the full workflow, use the notebook: [Running a simulation](Running-a-Simulation).

## Where to run from

The package is imported as `radcluster_code.py_utils`, so the repository root must be on Python's path. Either start Python in the repository root, or set `PYTHONPATH`:

```bash
cd RadCluster
PYTHONPATH=$PWD python my_case.py
```

In the container image the package is imported as `py_utils`, because `/work` is the module itself.

## A first run

```python
from radcluster_code.py_utils.simulation import RadClusterSimulation

sim = RadClusterSimulation(
    I=1000, V=1000,                 # largest SIA / vacancy cluster size tracked
    solver_mode='full_system',      # 'full_system' | 'active_window'
    equations='discrete',           # 'discrete' | 'bin_moment'
    cascade='fission',              # helium reduction: 'fission' (Case 2) | 'fusion' (Case 1)
    he_kinetics='quasi_steady_state',
    i_mobile=5, v_mobile=2,         # largest mobile SIA / vacancy cluster
)

print(sim.material_rag.summary())   # the EUROFER-97 reaction graph
print(sim.state_layout.describe())  # the blocks of the declared graph

results = sim.run(solver_config={
    't_span': (1e-6, 1e-2),         # seconds
    'n_points': 40,
    'log_time': True,
    'rtol': 1e-6, 'atol': 1e-20,
    'solver_method': {'linsol': 'gmres', 'preconditioner': 'Woodbury'},
})

print(results['dose'][-1], results['swelling'][-1])
print(results['delta_FP'][-1], results['delta_He'][-1])
```

The constructor reads `radcluster_code/input/input_parameters.xlsx`, computes the rate constants and declares the reaction graph. It prints:

```
RAG 'EUROFER-97': 3 populations, 24 edges, 24 kernels; edge classes active: {'growth': 4,
'shrinkage': 3, 'dissociation': 3, 'recombination': 1, 'annihilation': 2,
'inter_population': 1, 'solute_trapping': 1, 'coalescence': 4, 'source': 2, 'sink': 3}
```

`cascade` selects the helium reduction. The cascade source (its survival efficiency, clustered fractions and sizes) is the column of the `Production` sheet named by the workbook's `spectrum`; for a fusion case set both.

`run` writes the parameter file, starts the C++ solver, reads its binary output and computes the derived quantities. With `save_output=True` (the default) it writes a run directory:

```
radcluster_code/output/20261002_150956_full_system_full_CD_fission_I1000V1000_im5vm2/
```

See [Outputs and provenance](Outputs-and-Provenance) for what is in it and in `results`.

## Change a parameter for one run

The physics comes from the workbook. A value is changed for one run by writing it into the dictionary of its sheet and rebuilding the rates. The key is the workbook's *Symbol* column.

```python
inp = sim.input_data
inp.reactions['T'] = 673.0          # temperature, K
inp.reactions['rho_d'] = 5e14       # network dislocation density, m^-2
inp._calculate_derived()
sim.rebuild_rates()
```

The notebook does this for a whole dictionary, `PARAM_OVERRIDES`: [Inputs and parameters](Inputs-and-Parameters).

## The bin-moment equations

`equations='bin_moment'` tracks sizes up to `i_discrete` one by one and groups the larger sizes into logarithmic bins. The layout is set through the same dictionaries, not through the constructor:

```python
sim = RadClusterSimulation(I=4000, V=2000, solver_mode='full_system',
                           equations='bin_moment', cascade='fission',
                           he_kinetics='quasi_steady_state', i_mobile=5, v_mobile=1)

inp = sim.input_data
inp.reactions.update({'i_discrete': 40, 'v_discrete': 5,
                      'I_bin': 12, 'V_bin': 12, 'shape_function': 'linear'})
inp._calculate_derived()
sim.rebuild_rates()

re = sim.rate_equations             # read the layout back from here
print(re.i_discrete, re.I_bin, re.v_discrete, re.V_bin, re.n_mom, re.N_eq)
```

- `I_bin` and `V_bin` are targets. The realised count can differ by one; read it back from `sim.rate_equations`.
- Passing `i_discrete`, `I_bin` or `shape_function` to the constructor raises `TypeError`, with this procedure in the message.

## Long runs

`run_adaptive` integrates in short segments. Between segments it can double the size domain when the distribution reaches the top of the grid, and it advances the network dislocation density. The notebook uses it for every run.

```python
sim.run_tag = 'smoke'               # optional: appears in the run directory's name
results = sim.run_adaptive(
    solver_config={'t_span': (1e-6, 1e3), 'n_points': 30, 'log_time': True,
                   'rtol': 1e-5, 'atol': 1e-20,
                   'loop_conversion': 1,          # the ½⟨111⟩ / ⟨100⟩ loop split
                   'solver_method': {'linsol': 'gmres', 'preconditioner': 'Woodbury'}},
    boundary_threshold=0.05,        # adapt when more than 5 % of the content is at the boundary
    max_doublings=0,                # 0: keep the domain fixed
    points_per_segment=10,
)
```

## Inspect the model without running

```python
from radcluster_code.py_utils.core import EdgeClass

rag = sim.material_rag
for p in rag.populations:
    print(p.polarity.label, p.name, p.n_min, p.mobile_max)
for e in rag.edges:
    print(e.edge_class.value, e.label, e.kernel)
print(rag.edges_of_class(EdgeClass.COALESCENCE))
```

See [EUROFER-97 reaction graph](EUROFER-97-Reaction-Graph).
