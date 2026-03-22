# Docs — Shared Documentation

Shared reference materials for the EuroferMicrostructure project.

## Structure

```
Docs/
├── Database/       # Experimental radiation microstructure databases
├── Formulation/    # Rate-equation derivations and code documentation (PDF)
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

Key reference PDFs for rate-equation models:

| File | Contents |
|---|---|
| `Rate_Equations.pdf` | General rate-equation theory |
| `Rate_Equations_Code.pdf` | Code implementation notes |
| `Rate_equation_Eurofer.pdf` | EUROFER97-specific formulation |
| `Ghoniem1985--*.pdf` | Ghoniem 1985 helium clustering paper |
| `Ghoniem-Thesis.pdf` | Ghoniem thesis: cluster dynamics framework |
| `Gao-Ghoniem.pdf` | Gao-Ghoniem sink strength expressions |

## Literature/

Published papers on radiation effects in ferritic-martensitic steels (EUROFER97, F82H, T91):
- TEM loop and void characterization studies
- Ion vs neutron irradiation comparisons
- Atom probe tomography (APT) results
