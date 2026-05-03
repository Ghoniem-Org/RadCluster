# EuroferExperiments — Radiation Microstructure Model

Mean-field rate-equation model for radiation-induced microstructure evolution in **EUROFER97 / ferritic-martensitic (FM) steel**.

## Physics

Tracks the evolution of point defects and extended defects under irradiation:

| Species | Symbol | Description |
|---|---|---|
| Vacancies | Cv | Mobile mono-vacancies |
| Interstitials | Ci | Mobile mono-interstitials |
| Di-interstitials | C2i | Immobile clusters |
| Tri-interstitials | C3i | Immobile clusters |
| 1/2⟨111⟩ loops | CiL_111 (density), ril_111 (radius) | Glissile interstitial loops |
| ⟨100⟩ loops | CiL_100 (density), ril_100 (radius) | Sessile interstitial loops |
| Voids | C_void (density), r_void (radius) | Vacancy clusters / voids |
| Interstitial traps | Ctrap_i | Point defects trapped at solutes/precipitates |
| Vacancy traps | Ctrap_v | Point defects trapped at solutes/precipitates |

## Sink strengths

- Network dislocations (density ρ): vacancy and interstitial sinks with bias Z_N_i
- Dislocation loops: self-consistent sink terms
- Voids: spherical sinks
- Traps: fixed concentration CT0 with capture radius r_trap

## ODE system

Solved with `scipy.integrate.solve_ivp` (LSODA method by default).
State vector `y` has one entry per species above.

## File structure

```
EuroferExperiments/
├── CLAUDE.md               (this file)
├── code/
│   ├── EuroferMicro.ipynb          # Main simulation notebook
│   ├── EuroferMicroData.ipynb      # Experimental data analysis
│   └── EuroferRadiationMicro.ipynb # Radiation effects comparison
├── py_utils/
│   ├── __init__.py
│   ├── input_data.py       # InputData class — reads 3-sheet Excel
│   ├── rate_equations.py   # RateEquations class — ODE RHS
│   ├── reaction_rates.py   # ReactionRates class — pre-computed rate constants
│   ├── simulation.py       # EuroferMicroSimulation — orchestrator + plotting
│   ├── post_process.py     # (placeholder) separate post-processing
│   └── visualization.py    # (placeholder) separate plotting utilities
├── input/
│   ├── input_parameters.xlsx        # Ion irradiation conditions
│   └── input_parameters_Neutron.xlsx
└── output/                 # Simulation results (gitignored)
```

## Running

```python
import sys
sys.path.insert(0, '../')   # from EuroferExperiments/code/
from EuroferExperiments.py_utils.simulation import run_eurofermicro_simulation

SIMULATION_CONFIG = {
    't_begin': 1e-6, 't_end': 1e8, 'n_points': 1000,
    'method': 'LSODA', 'rtol': 1e-6, 'atol': 1e-20, 'log_time': True
}
results, sim = run_eurofermicro_simulation(SIMULATION_CONFIG)
```
# Sink Strength of Tempered Martensite under Fission Neutron Irradiation

## Overview

There is no single closed-form equation specifically for "tempered martensite" as a
microstructural entity — but the sink strength of the tempered martensite microstructure
can be built up systematically from the contributions of its individual microstructural
components using standard rate theory.

In rate theory the **total sink strength** $k^2$ (units m⁻²) governs how fast point
defects are removed from the matrix. The vacancy and SIA concentrations in steady state
are inversely proportional to $k^2$, so a high sink strength means low vacancy
supersaturation and hence suppressed void nucleation. Contributions add:

$$k^2_\text{total} = k^2_\text{disl} + k^2_\text{GB} + k^2_\text{voids} + k^2_\text{prec} + k^2_\text{loops}$$

---

## 1. Dislocations (dominant term in fresh tempered martensite)

For straight dislocations:

$$k^2_\text{disl} = Z_d \, \rho_d$$

where:
- $\rho_d$ = dislocation line density (m⁻²)
- $Z_d$ = capture efficiency (bias factor)
  - Vacancies: $Z_d^V \approx 1.0$–$1.1$
  - SIAs: $Z_d^I \approx 1.1$–$1.2$
  - Dislocation bias: $B = Z_d^I / Z_d^V - 1 \approx 0.02$–$0.10$

In freshly tempered martensite $\rho_d \sim 10^{14}$–$10^{15}$ m⁻², giving
$k^2_\text{disl} \sim 10^{14}$–$10^{15}$ m⁻². This is the **largest single term**
and the primary reason F/M steels resist swelling.

---

## 2. Grain and Lath Boundaries (planar sinks)

For a polycrystal with planar grain boundary sinks (spherical grain approximation,
perfect sink, Brailsford & Bullough 1972):

$$k^2_\text{GB} \approx \frac{6 \, S_\text{GB}}{d_g}$$

where:
- $S_\text{GB}$ = boundary capture efficiency ($\approx 1$ for perfect sinks)
- $d_g$ = relevant grain / lath diameter (m)

In tempered martensite the relevant length scales are:

| Structural unit | Typical size |
|---|---|
| Prior austenite grain | 10–50 µm |
| Martensite packet | 2–10 µm |
| Martensite lath width | 0.2–0.5 µm |

**Lath boundaries dominate** because they are the finest scale. Using
$d_\text{lath} \sim 0.3$ µm:

$$k^2_\text{lath} \approx \frac{6}{0.3 \times 10^{-6}} \sim 2 \times 10^{7} \text{ m}^{-2}$$

This is much smaller than the dislocation term at beginning of life, but becomes
relatively more important as irradiation anneals out the dislocation network at
higher doses.

---

## 3. Precipitates and Oxide Particles (spherical sinks)

For a population of spherical sinks of radius $r_p$ and number density $N_p$:

$$k^2_\text{prec} = 4\pi \, r_p \, N_p \, Z_p$$

where:
- $Z_p$ = interface capture efficiency
  - Incoherent M₂₃C₆, MX: $Z_p \sim 0.1$–$0.5$
  - ODS oxide particles: $Z_p \to 1$

Typical values in tempered martensite: $r_p \sim 10$–$50$ nm,
$N_p \sim 10^{21}$–$10^{22}$ m⁻³, giving
$k^2_\text{prec} \sim 10^{13}$–$10^{14}$ m⁻².

Significant but smaller than the dislocation term unless the precipitate density is
very high (as in ODS steels, where $k^2_\text{prec}$ can rival $k^2_\text{disl}$).

---

## 4. Dislocation Loops (radiation-induced, growing with dose)

Loops are treated identically to straight dislocations:

$$k^2_\text{loops} = Z_L \, \rho_L$$

where the loop line density equivalent is:

$$\rho_L = \frac{\pi}{4} \, d_L^2 \, N_L$$

with $d_L$ = mean loop diameter (m) and $N_L$ = loop number density (m⁻³).

Under irradiation, as loops nucleate and grow, $k^2_\text{loops}$ increases
and eventually rivals or exceeds the initial dislocation contribution, partially
compensating for the loss of network dislocations by climb and recovery.

---

## 5. Voids (once nucleated, they feed back as sinks)

$$k^2_\text{voids} = 4\pi \, r_v \, N_v$$

with $Z_v \approx 1$ (voids are **unbiased** sinks — they capture vacancies and SIAs
with equal probability). Voids grow because the dislocation bias channels excess SIAs
to dislocations, leaving a net vacancy flux to voids. The void sink strength grows
with dose once the incubation period ends.

---

## Full Expression

Combining all terms:

$$\boxed{k^2_\text{total} = Z_d^{V,I}\,\rho_d + Z_L^{V,I}\,\rho_L + 4\pi r_p N_p Z_p + \frac{6}{d_\text{lath}} + 4\pi r_v N_v}$$

### Relative magnitudes at start of irradiation (T91 / HT9)

| Sink | $k^2$ (m⁻²) | Notes |
|---|---|---|
| Network dislocations | $10^{14}$–$10^{15}$ | Dominant; anneals with dose |
| Lath boundaries | $10^{7}$–$10^{8}$ | Secondary; microstructurally stable |
| M₂₃C₆ / MX precipitates | $10^{13}$–$10^{14}$ | Depends on $Z_p$ and $N_p$ |
| Irradiation loops (growing) | $10^{12}$ → $10^{14}$ | Increases with dose |
| Voids (post-incubation) | $10^{11}$–$10^{13}$ | Small until steady-state swelling |

---

## Swelling Suppression Condition

The condition for void swelling suppression is that the **recombination rate**
dominates over the bias-driven net flux to voids:

$$\frac{k^2_\text{total}}{k^2_\text{voids}} \gg \frac{B}{1 - B}$$

Since $B \sim 0.02$–$0.05$ in bcc Fe, even a modest total sink strength relative
to the void sink strength is sufficient to suppress growth. The dense lath dislocation
network at $\rho_d \sim 10^{14}$–$10^{15}$ m⁻² is what drives this during the
incubation period.

---

## Evolution with Irradiation Dose

The sink strength of tempered martensite is **not static**. With increasing dose:

1. **Network dislocation recovery** — dislocations climb to sinks and annihilate;
   $\rho_d$ decreases, reducing $k^2_\text{disl}$. This is the primary driver of
   the eventual onset of steady-state swelling after the incubation period.

2. **Loop accumulation** — irradiation-induced a⟨100⟩ and ½a⟨111⟩ loops nucleate
   and grow, partially replacing the lost network dislocation sink strength.

3. **Precipitate evolution** — M₂₃C₆ and MX coarsen (radiation-enhanced diffusion),
   reducing $N_p$ and hence $k^2_\text{prec}$. At very high dose (>250 dpa in T91/T92
   under ion irradiation) M₂X precipitates nucleate on dislocation lines and may add
   new sink sites.

4. **Void growth** — as $k^2_\text{voids}$ grows, voids increasingly act as
   recombination sites, partially compensating for lost dislocations — but also
   consuming more vacancies and growing themselves.

The net result is that $k^2_\text{total}$ decreases with dose in the incubation
regime, and swelling eventually accelerates when $k^2_\text{disl}$ has degraded
sufficiently.

---

## Void Nucleation and Stabilisation Summary

### Nucleation mechanism
1. **Cascade collapse** — vacancy-rich cascade core collapses into embryo clusters
   (~2–10 vacancies) on ~10 ps timescale
2. **Thermal migration** — surviving vacancies migrate (migration energy in Fe ~0.55 eV)
   and accumulate at heterogeneous nucleation sites
3. **Heterogeneous nucleation** — preferred sites are dislocation lines, precipitate
   interfaces (M₂₃C₆, MX, G-phase), and lath boundaries

### Stabilisation by helium
A pure vacancy embryo in bcc Fe is **not thermodynamically stable** at reactor
temperatures — thermal emission dissolves it. Stabilisation requires **helium**:

- He is produced by (n,α) transmutation (primarily on ⁵⁸Ni, ⁶⁰Ni, B impurities)
- He is insoluble in Fe (binding energy to vacancy ~1 eV) and is trapped by vacancies
- Internal He pressure suppresses thermal vacancy emission from the embryo
- Even 1–3 He atoms per cluster dramatically reduce the nucleation barrier
- **Bubble-to-void transition**: once a He bubble reaches a critical radius where the
  vacancy flux in exceeds thermal emission, it converts to a growing void

Typical He production rates in fast reactors: **0.1–2 appm He/dpa**. Mixed-spectrum
and thermal reactors produce 10–100× more He/dpa, dramatically increasing void
nucleation density.

---

## Typical Void Data from Fission Neutron Irradiation

| Material | Reactor | Dose (dpa) | T (°C) | Cavity diam. (nm) | Cavity density (m⁻³) | Swelling (%) |
|---|---|---|---|---|---|---|
| EUROFER97 | BOR-60 | 15 | 330 | 2.6 | 1.4×10²⁰ | <0.01 |
| EUROFER97 | BOR-60 | 32 | 335 | 1.6 | 1.1×10²⁰ | <0.01 |
| T91 | BOR-60 | 15.4 | 376 | <2 (bimodal) | ~10²⁰–10²¹ | ~0.01 |
| T91 | BOR-60 | 35.1 | 376 | bimodal: <2 + >2 | ~5×10²⁰ | ~0.02 |
| HT9 | BOR-60 | 17.1 | 377 | <2 | ~3×10²⁰ | 0.02 |
| HT9 | BOR-60 | 35.1 | 377 | bimodal | ~2×10²⁰ | 0.07 |
| HT9 | EBR-II/FFTF | ~200 | 420 | — | — | 0.09–1.02 |
| Pure Fe | Multiple | 1 | 400 | ~5 | ~10²³ | 0.07 |

**Key observations:**
- Cavities in F/M steels are very small (1–3 nm) compared to austenitic SS at
  comparable conditions
- Bimodal distributions (He bubbles + sparse larger voids) emerge at higher doses,
  reflecting two-stage bubble→void transition
- Cavity density (10²⁰–10²¹ m⁻³) is 1–3 orders of magnitude below loop density
- Peak swelling temperature: 400–450°C
- Incubation period: ~50–100 dpa before measurable steady-state swelling
- Steady-state swelling rates: HT9 ~0.033%/dpa; T91 ~0.007%/dpa; T92 ~0.002%/dpa
  (from ion irradiation data; neutron data consistent at lower doses)

---

## Key References

| Reference | Topic |
|---|---|
| Mansur, L.K. (1994) *J. Nucl. Mater.* 216, 97–123 | Rate theory framework; sink strength formalism |
| Brailsford & Bullough (1972) *J. Nucl. Mater.* 44, 121 | Original sink strength derivation |
| Garner & Toloczko (1994) *J. Nucl. Mater.* 212–215, 289 | Swelling suppression in F/M steels |
| Klueh & Harries (2001) *High-Chromium Ferritic and Martensitic Steels for Nuclear Applications*, ASTM | Comprehensive microstructure–swelling treatment |
| Gao et al. (2018) *Acta Mater.* | T91 BOR-60 void and loop data |
| Gao et al. (2019) *J. Nucl. Mater.* | HT9 BOR-60 void and loop data |
| Getto et al. (2015) *J. Nucl. Mater.* 462, 458 | HT9/T91/T92 ion irradiation steady-state swelling rates |
| Was, G.S. (2017) *Fundamentals of Radiation Materials Science*, Springer | Textbook treatment of sink strength and void swelling |
- Ghoniem (1985): rate equations for interstitial loop evolution
- Gao & Ghoniem: void and loop sink strengths
- See `../Docs/Formulation/` for full derivations
- See `../Docs/Literature/` for EUROFER97 experimental data
