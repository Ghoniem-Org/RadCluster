# Docs — Shared Documentation

Shared reference materials for the EuroferMicrostructure project.

## Structure

```
Docs/
├── Database/       # Experimental radiation microstructure databases
├── Formulation/    # Derivations, paper sections, studies — organized by topic
│                   # (canonical — supersedes legacy docs/Rate Equations/)
└── Literature/     # Peer-reviewed papers cited across modules
```

## Database/

Excel databases of experimental TEM/APT measurements in ferritic-martensitic steels:

| File | Contents |
|---|---|
| `FerriticSteels_RadiationDatabase_Combined.xlsx` | Combined database: loop density, size, void swelling vs dose/temperature |
| `InterstitialLoop.xlsx` | Interstitial loop data subset |
| `MicroData.xlsx` | General microstructure data |
| `Void.xlsx` | Void/bubble data subset |

## Formulation/

Organized by topic; see [`Formulation/README.md`](Formulation/README.md) for the full index.

| Folder | Topic |
|---|---|
| `cluster_dynamics_framework/` | Generalized graph-based CD framework, paper drafts, novelty analyses |
| `reaction_admissibility_graph/` | RAG figures/exports and Supplementary S1 (architecture, solver coupling) |
| `reaction_kernels/` | 1D reaction kernels, detailed balance |
| `state_space_reduction_and_solvers/` | Bin moments, GPU integrator plan |
| `dislocation_loops/` | Loop conversion, loop->network loss |
| `segregation_and_precipitation/` | RIS and precipitation plan |
| `steel_applications_and_calibration/` | F/M-steel application, digital twin, calibration roadmap |
| `verification_campaign/` | Verification study: plan, results, figures, tables |

## Literature/

Published papers on radiation effects in ferritic-martensitic steels (EUROFER97, F82H, T91):
- TEM loop and void characterization studies
- Ion vs neutron irradiation comparisons
- Atom probe tomography (APT) results
