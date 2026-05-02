# [ARCHIVED 2026-05-02] Eurofer_CD — Cluster Dynamics for bcc Fe / EUROFER97

> **This module is archived.** It has been superseded by `RadCluster_1_0/`
> (originally archived as superseded by `Expanded_Eurofer_CD/`, which was
> renamed 2026-05-02),
> which generalizes the master equations (Ghoniem 2026 formulation), adds
> bin-moment reduction, the Woodbury preconditioner, adaptive domain doubling,
> and full conservation diagnostics. This directory is preserved for
> reproducibility of prior results only — no further development. The last
> active commit is reachable via the `eurofer_cd-final` git tag.

---

# Eurofer_CD — Cluster Dynamics for bcc Fe / EUROFER97

Mean-field cluster-size–resolved rate equations for radiation-induced defect
cluster evolution in **EUROFER97 / ferritic-martensitic (F/M) steel**, including
explicit helium–vacancy cluster kinetics.

---

## Physics

### Cluster species and state vector

Each cluster is labelled by a pair (n, ℓ) where:
- n > 0: SIA cluster of size n (interstitial loops)
- n < 0: vacancy cluster of size |n| (voids / vacancy loops)
- n = 0, ℓ = 1: free He atom
- (m, ℓ): He–vacancy cluster with m vacancies and ℓ He atoms

**State vector** (ODE unknowns):

| Index       | Variable        | Description                     |
|-------------|-----------------|----------------------------------|
| 0 … Ni−1   | C_{n=1..Ni}     | SIA clusters                    |
| Ni … Ni+Nv−1| C_{-m=1..Nv}   | Vacancy clusters (pure)         |
| Ni+Nv       | C_{0,1}         | Free He                         |

---

### Rate equation structure (PDF Section 5.7)

The master equation for each species includes:

- **Cascade production**: power-law spectrum ε_m^(i), ε_m^(v) (Section 1)
- **3D Waite capture**: K_AB = 4π·r_AB·(D_A + D_B) (Section 5.3, eq. 87)
- **1D glide**: Glissile SIA clusters (n ≥ n_1D = 4) with effective 1D→3D rate (Section 3)
- **Thermal emission**: α_void(m) = K_vv(m)·C_v^eq_surf(m) (capillary model, Section 4)
- **He capture and emission**: He–vacancy cluster kinetics
- **Dislocation sinks**: Z_i·ρ_d·D_i·C_i and Z_v·ρ_d·D_v·(C_v − C_v^eq)

---

### Binding energies (Section 4)

**Void / vacancy cluster binding** (capillary approximation, eq. 62–63):
```
E_b_void(m) = E_f_v − A_void · [m^(2/3) − (m−1)^(2/3)]
A_void = 4π · γ_s · r_0^2  [eV]
r_0 = (3Ω/4π)^(1/3)
```
For m = 1: E_b = E_f_v (monomer formation energy).

**Interstitial loop binding** (blended power-law / continuum, eq. 83–85):
```
E_b_loop(n) = E_b_inf − (E_b_inf − E_b_2i) · exp(−(n−2)/n_trans)
E_b_2i ≈ 0.80 eV  (bcc Fe, atomistic)
E_b_inf ≈ 1.80 eV  (large ½⟨111⟩ loop limit)
```

**He–vacancy cluster binding** (Caturla et al. atomistic data + fit):
```
E_b_He(m, ℓ) = E_b_HeV + δ_He · (ℓ/m)^β_He
E_b_HeV ≈ 2.60 eV  (He in mono-vacancy)
```

---

### He–vacancy state-space reduction (Section 5.6.5)

Two physics-based options controlled by `he_mode` in Model_Parameters:

**`decoupled`** (fission, low He/dpa ≈ 0.5–1 appm/dpa):
- He is not tracked as a separate cluster species.
- Mean He loading ⟨ℓ⟩_m from quasi-steady-state He mass balance.
- He shifts effective void binding energy:
  ```
  E_b_eff(m) = E_b_void(m) + ⟨ℓ⟩_m · ∂E_b_He/∂ℓ
  ```
- Reduces He-vacancy state space to **zero extra equations**.

**`fast_eq`** (general, including fusion with He/dpa ≈ 10 appm/dpa):
- He distribution within each void class equilibrates rapidly.
- Track marginal C_{-m}^tot = Σ_ℓ C_{-m,ℓ} (one equation per void class).
- ⟨ℓ⟩_m from QSS He mass balance:
  ```
  ⟨ℓ⟩_m = K_HeV(m)·C_He / β_He_emit
  ```
- He-pressure correction applied to effective void emission rate.

---

## Material parameters (bcc Fe / EUROFER97)

| Symbol    | Value      | Units  | Source                    |
|-----------|-----------|--------|---------------------------|
| a         | 2.87 Å    | m      | bcc Fe lattice            |
| Ω         | 1.18×10⁻²⁹| m³    | a³/2 (bcc)               |
| E_f_v     | 1.73      | eV     | DFT (Malerba 2021)        |
| E_m_v     | 0.55      | eV     | DFT / experiment          |
| E_m_i     | 0.013     | eV     | DFT (fast 1D glide)       |
| E_b_2i    | 0.80      | eV     | atomistic                 |
| E_m_He    | 0.06      | eV     | DFT (Borodin 2014)        |
| E_b_HeV   | 2.60      | eV     | Caturla (2005)            |
| γ_s       | 1.50      | J/m²   | Capillary model           |

---

## ODE system

Solved with `scipy.integrate.solve_ivp` (LSODA) in segmented mode.
Default: N = Ni + Nv + 1 equations (100 + 100 + 1 = 201).

---

## File structure

```
Eurofer_CD/
├── CLAUDE.md                    (this file)
├── create_excel.py              (generates input/input_parameters.xlsx)
├── py_utils/
│   ├── __init__.py
│   ├── binding_energies.py      E_b_void, E_b_loop, E_b_He, capture_radius
│   ├── defect_production.py     Cascade production spectra (fission / fusion)
│   ├── input_data.py            InputData class — reads 3-sheet Excel
│   ├── reaction_rates.py        ReactionRates — Waite capture + emission rates
│   ├── rate_equations.py        RateEquations — ODE RHS with he_mode dispatch
│   ├── simulation.py            EuroferCDSimulation — segmented LSODA solver
│   ├── post_process.py          Derived quantities: totals, swelling, He content
│   └── visualization.py         EuroferCDVisualizer — standard figure set
├── codes/
│   └── eurofer_cd.ipynb         Main simulation notebook
├── input/
│   └── input_parameters.xlsx    3-sheet parameter file (Eurofer97 defaults)
└── output/                      Timestamped run directories (gitignored)
```

---

## Running

```python
import sys
sys.path.insert(0, '../')   # from Eurofer_CD/codes/

from Eurofer_CD.py_utils.simulation import EuroferCDSimulation

sim = EuroferCDSimulation(Ni=100, Nv=100, he_mode='decoupled')
results = sim.run_simulation(t_span=(1e-8, 1e6))
```

## Key references

- Ghoniem, N.M. (2024), *Formulation of Cluster Dynamics Equations for
  Irradiated Ferritic-Martensitic Steels* (internal report; `docs/Formulation/Rate_Equations.pdf`)
- Malerba, L. et al. (2021), J. Nucl. Mater. 543, 152463
- Stoller, R.E. (2000), J. Nucl. Mater. 276, 22–32
- Nordlund, K. et al. (2018), Nature Commun. 9, 1084
- Caturla, M.J. et al. (2005), J. Nucl. Mater. 336, 73–82
- Waite, T.R. (1957), Phys. Rev. 107, 463
