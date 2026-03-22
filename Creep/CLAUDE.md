# Creep — Thermal and Irradiation Creep Model

**Placeholder module** — to be implemented by adapting `Fluor_Zr/Creep/` for ferritic-martensitic steel.

## Physics

Dislocation-mechanics-based model for:

- **Thermal creep**: dislocation glide + climb driven by applied stress
- **Irradiation creep**: stress-biased absorption of radiation-induced point defects at dislocations (SIPA mechanism)
- **Irradiation growth**: stress-free dimensional change driven by anisotropic defect absorption

State variables:
| Variable | Description |
|---|---|
| ε_cr | Creep strain |
| ρ_m | Mobile dislocation density |
| ρ_f | Forest dislocation density |

## Coupling

Creep module can consume microstructure output from EuroferExperiments (defect concentrations vs dose) as input to irradiation creep terms.

## Status

All files in `py_utils/` and `cpp_utils/` are **placeholders**.

## To implement

1. Port model from `Fluor_Zr/Creep/py_utils/` and `cpp_utils/`
2. Replace Zr creep parameters with EUROFER97 values
3. Build C++ ARKODE solver: configure `cpp_utils/CMakeLists.txt`
4. Validate against creep data in literature

## References

- See `../Docs/Formulation/` for thermal creep model derivation
- See `../Docs/Literature/` for EUROFER97 creep data
