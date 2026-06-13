# RadCluster — Project Overview

Physics-based simulation suite for **EUROFER97 / ferritic-martensitic steel** behaviour under irradiation and thermal loading. Modelled after the structure of the `Fluor_Zr` repository.

## Repository Layout

```
RadCluster/
├── RadCluster_2_1/     # Active development — adds dislocation evolution + RIS (clone of 2_0)
├── RadCluster_2_0/     # Graph-based cluster dynamics — stable baseline / reference
├── archive/            # Archived modules (read-only, kept for reproducibility)
│   ├── Eurofer/            # Archived earlier microstructure work
│   ├── Eurofer_CD/         # Archived 2026-05-02 — superseded by RadCluster_1_0
│   ├── RadCluster_1_0/     # Archived 2026-06-13 — superseded by RadCluster_2_0
│   ├── Monomer_CD/         # Archived 2026-06-13 — cluster dynamics scaling reference
│   └── Zr_RadCluster_1_0/  # Archived 2026-06-13 — zirconium cluster-dynamics variant
├── docs/               # Shared documentation, literature, databases
│   ├── Database/           # Experimental microstructure databases (xlsx)
│   ├── Formulation/        # Rate-equation derivations and code notes (PDF)
│   └── Literature/         # Peer-reviewed papers (PDF)
└── requirements.txt
```

## Module Status

| Module | Status | Description |
|---|---|---|
| `RadCluster_2_1/` | **Active (dev)** | Clone of `RadCluster_2_0` adding two capabilities: **(a)** a dislocation-evolution / loop→network loss edge that saturates loop density, and **(b)** Radiation-Induced Segregation (RIS) + solute precipitation. Plans: [`docs/Formulation/loop_network_loss.tex`](docs/Formulation/loop_network_loss.tex) and [`docs/Formulation/radcluster_2_1_RIS_plan.tex`](docs/Formulation/radcluster_2_1_RIS_plan.tex). See `RadCluster_2_1/CLAUDE.md` §0. |
| `RadCluster_2_0/` | **Active (baseline)** | Generalized graph-based cluster dynamics (Ghoniem 2026). Two-layer RAG architecture (abstract core + EUROFER-97 host declaration). Stable reference for the 2_1 work. Notebooks: `RadCluster_2_0.ipynb` (simulation driver) and `EuroferExperiments.ipynb`. |
| `archive/RadCluster_1_0/` | Archived 2026-06-13 | Superseded by `RadCluster_2_0/`. Earlier (non-graph) generalized cluster dynamics. |
| `archive/Monomer_CD/` | Archived 2026-06-13 | Monomer-mobility cluster dynamics scaling reference (Ghoniem & Cho 1979, no He). |
| `archive/Zr_RadCluster_1_0/` | Archived 2026-06-13 | Zirconium cluster-dynamics variant of `RadCluster_1_0`. |
| `archive/Eurofer/` | Archived | Earlier EUROFER microstructure notebooks. |
| `archive/Eurofer_CD/` | Archived 2026-05-02 | Superseded by `RadCluster_1_0/`. Last active state at git tag `eurofer_cd-final`. |

## Conventions (mirror Fluor_Zr)

### Module structure
Each active module follows:
```
<Module>/
├── CLAUDE.md           # Physics description, solver notes
├── codes/              # Jupyter notebooks and Python test scripts
│   ├── Notebooks/
│   └── Python_Testing/
├── py_utils/           # Python utilities package
│   ├── __init__.py
│   ├── input_data.py   # InputData class (reads 3-sheet Excel)
│   ├── rate_equations.py
│   ├── simulation.py   # Orchestrator; writes timestamped output/
│   ├── cpp_bridge.py   # subprocess wrapper for C++ solver
│   └── visualization.py
├── cpp_utils/          # C++ solver (optional, for performance)
│   ├── CMakeLists.txt
│   ├── parameters.h
│   └── rate_equations.cpp / .h
├── input/              # input_parameters.xlsx (3 sheets)
├── output/             # Timestamped run directories (gitignored)
└── build/              # CMake build artifacts (gitignored)
```

### Excel input format
All modules share the same 3-sheet structure:
- `Material_Environment` — T (K), G (dpa/s), rho (m⁻²), stress, ...
- `Physical_Properties` — Omega, a, b_111, b_100, migration energies, ...
- `Model_Parameters` — Z-factors, cluster cutoffs, capture radii, ...

### Output format
Simulation results go in timestamped subdirectories:
```
output/YYYYMMDD_HHMMSS_<git-hash>/
├── provenance.md       # timestamp, git SHA, key parameters
├── results.pkl         # full ODE solution (binary)
├── summary.csv         # tabulated macroscopic quantities
└── plots/              # PNG figures
```

### C++ solvers
- Use SUNDIALS 7.1.1 (CVODE)
- Built with CMake; binaries land in `build/` (gitignored)
- Invoked from Python via `py_utils/cpp_bridge.py` subprocess wrapper

## Environment

```bash
# Install dependencies
pip install -r requirements.txt

# Register Jupyter kernel (optional)
python -m ipykernel install --user --name radcluster --display-name "RadCluster"
```

## Shared Resources

- `docs/Database/` — experimental radiation microstructure data for ferritic-martensitic steels
- `docs/Formulation/` — rate-equation derivations (canonical)
- `docs/Literature/` — reference papers
