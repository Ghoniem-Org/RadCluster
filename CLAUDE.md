# EuroferMicrostructure — Project Overview

Physics-based simulation suite for **EUROFER97 / ferritic-martensitic steel** behaviour under irradiation and thermal loading. Modelled after the structure of the `Fluor_Zr` repository.

## Repository Layout

```
EuroferMicrostructure/
├── EuroferExperiments/     # Radiation microstructure evolution (active)
├── Expanded_Eurofer_CD/    # Cluster dynamics — current research code (active)
├── Full_CD/                # Cluster dynamics scaling reference (active)
├── ClusterDynamics/    # Cluster-size resolved rate equations (placeholder)
├── Creep/              # Thermal + irradiation creep (placeholder)
├── EuroferProps/       # Material property fitting utilities (placeholder)
├── archive/            # Archived modules (read-only, kept for reproducibility)
│   └── Eurofer_CD/         # Archived 2026-05-02 — superseded by Expanded_Eurofer_CD
└── Docs/               # Shared documentation, literature, databases
    ├── Database/       # Experimental microstructure databases (xlsx)
    ├── Formulation/    # Rate-equation derivations and code notes (PDF)
    └── Literature/     # Peer-reviewed papers (PDF)
```

## Module Status

| Module | Status | Description |
|---|---|---|
| `EuroferExperiments/` | **Active** | Mean-field rate equations for loop + void evolution in EUROFER97 |
| `Expanded_Eurofer_CD/` | **Active** | Generalized cluster dynamics (Ghoniem 2026) — current research code |
| `Full_CD/` | **Active** | Cluster dynamics scaling reference (Ghoniem & Cho 1979, no He) |
| `ClusterDynamics/` | Placeholder | Cluster-size-resolved model (adapt from Fluor_Zr/ClusterDynamics) |
| `Creep/` | Placeholder | Dislocation-mechanics creep model (adapt from Fluor_Zr/Creep) |
| `EuroferProps/` | Placeholder | Property fitting utilities (adapt from Fluor_Zr/ZrProps) |
| `archive/Eurofer_CD/` | Archived 2026-05-02 | Superseded by `Expanded_Eurofer_CD/`. Last active state at git tag `eurofer_cd-final`. |

## Conventions (mirror Fluor_Zr)

### Module structure
Each module follows:
```
<Module>/
├── CLAUDE.md           # Physics description, solver notes
├── code/               # Jupyter notebooks (main + development/)
├── py_utils/           # Python utilities package
│   ├── __init__.py
│   ├── input_data.py   # InputData class (reads 3-sheet Excel)
│   ├── rate_equations.py
│   ├── reaction_rates.py
│   ├── simulation.py   # Orchestrator; writes timestamped output/
│   ├── post_process.py # Macroscopic quantity derivation
│   └── visualization.py
├── cpp_utils/          # C++ solver (optional, for performance)
│   ├── CMakeLists.txt
│   ├── parameters.h
│   ├── rate_equations.cpp / .h
│   └── ode_solver.cpp  (Creep only)
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
- Use SUNDIALS 7.1.1 (CVODE for EuroferExperiments/ClusterDynamics, ARKODE for Creep)
- Built with CMake; binaries land in `build/` (gitignored)
- Invoked from Python via `py_utils/cpp_bridge.py` subprocess wrapper

## Environment

```bash
# Install dependencies
pip install -r requirements.txt

# Register Jupyter kernel (optional)
python -m ipykernel install --user --name eurofer_micro --display-name "EuroferMicrostructure"
```

## Shared Resources

- `Docs/Database/` — experimental radiation microstructure data for ferritic-martensitic steels
- `Docs/Formulation/` — rate-equation derivations (canonical; supersedes `docs/Rate Equations/`)
- `Docs/Literature/` — reference papers
