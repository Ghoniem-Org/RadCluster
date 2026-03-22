# ClusterDynamics — Cluster-Size Resolved Rate Equations

Cluster-size resolved (CD) rate equation model for radiation-induced defect clustering in **EUROFER97 / ferritic-martensitic (FM) steel**.

Extends the mean-field `EuroferExperiments` model by tracking the full number-density distribution of vacancy and interstitial clusters up to a user-specified maximum size `N_max`.

## Physics

### Species tracked

| Array | Symbol | Description |
|---|---|---|
| `fv[n]` | $f_v(n)$ | Vacancy clusters of size $n = 1 \ldots N_v$ |
| `fi[n]` | $f_i(n)$ | Interstitial clusters of size $n = 1 \ldots N_i$ |
| `CiL_111` | $C_{iL}^{111}$ | Number density of 1/2⟨111⟩ interstitial loops |
| `CiL_100` | $C_{iL}^{100}$ | Number density of ⟨100⟩ interstitial loops |
| `CiL_i_111` | $C_{iL,i}^{111}$ | Interstitial atoms in 1/2⟨111⟩ loops |
| `CiL_i_100` | $C_{iL,i}^{100}$ | Interstitial atoms in ⟨100⟩ loops |
| `C_void` | $C_\text{void}$ | Void number density |
| `r_void` | $r_\text{void}$ | Mean void radius |

Total state vector length: `N_v + N_i + 6`.

### Mobility assumptions

- $f_v(1)$: mobile (mono-vacancy)
- $f_v(n \geq 2)$: immobile
- $f_i(1)$: mobile (mono-interstitial)
- $f_i(2)$: glissile (di-interstitial, migration energy $E_{m,2i}$)
- $f_i(n \geq 3)$: immobile; once $n > N_{loop}$ they nucleate loops

### Cluster flux equations (size-space master equations)

For vacancy clusters of size $n$ ($1 \leq n \leq N_v$):

$$
\frac{df_v(n)}{dt} = J_v(n-1 \to n) - J_v(n \to n+1) + S_v(n)
$$

where the net flux in size space is:

$$
J_v(n \to n+1) = \beta_v(n)\,f_v(n) - \alpha_v(n+1)\,f_v(n+1)
$$

- **Absorption rate** $\beta_v(n) = 4\pi r_n D_v f_v(1) / \Omega$ — vacancy monomer capture by a cluster of size $n$, radius $r_n = r_0 n^{1/3}$.
- **Emission rate** $\alpha_v(n) = \beta_v(n-1)\exp(-E_b(n)/k_BT)$ — thermal emission of a vacancy monomer from size $n$, binding energy $E_b(n)$.

For interstitial clusters the same structure applies with $D_i$ and binding energies $E_{b,i}(n)$.

### Recombination between mobile species

$$
S_v(1) = G_v - \sum_n \beta_{iv}(n)\,f_i(1)\,f_v(1)
$$

where $\beta_{iv} = \omega_{iv}(D_i + D_v)/\Omega$ is the i–v recombination rate constant.

### Sink terms

Mobile species ($v_1$, $i_1$) are absorbed by:
- Network dislocations (density $\rho$, bias $Z_N$)
- Interstitial loops (density $C_{iL}$, mean radius $r_{iL}$)
- Voids (density $C_v$, radius $r_\text{void}$)

Equivalent sink concentrations carried over from `EuroferExperiments`:

$$
C_v^s = \frac{a^2}{z_c}\left[\rho_N + 2\pi\left(\frac{r_{111} C_{iL}^{111} + r_{100} C_{iL}^{100}}{\Omega}\right)\right]
$$

### Loop evolution

Loops nucleate when $f_i(N_{loop})$ grows, and evolve by absorption of $i_1$ and emission of $v_1$:

$$
\frac{dC_{iL,i}^{hkl}}{dt}
= \left(\frac{\ell'_{hkl}}{\ell}\right)
\sqrt{C_{iL,i}^{hkl}\,C_{iL}^{hkl}}
\left(Z_i^{hkl}\,\phi_i - Z_v^{hkl}\,\phi_v\right)
$$

### Void evolution

Void nucleation from vacancy cluster cascade term and growth:

$$
\frac{dr_\text{void}}{dt} = \frac{a^2}{z_c\,r_\text{void}}
\left[\omega_v f_v(1) - \omega_i f_i(1)
- \omega_v e_v\!\left(e^{2\gamma\Omega/r_\text{void}k_BT}-1\right)\right]
$$

## ODE system

The full ODE system is built in `py_utils/rate_equations.py`:

```
state = [fv[1..Nv], fi[1..Ni], CiL_111, CiL_100, CiL_i_111, CiL_i_100, C_void, r_void]
len   = Nv + Ni + 6
```

Solved with `scipy.integrate.solve_ivp` (LSODA). Cluster size cutoffs `N_v`, `N_i`, and `N_loop` are read from `Model_Parameters` sheet.

## File structure

```
ClusterDynamics/
├── CLAUDE.md               (this file)
├── code/
│   └── ClusterDynamics.ipynb   # Main simulation notebook
├── py_utils/
│   ├── __init__.py
│   ├── input_data.py       # InputData — reads 3-sheet Excel, adds CD cutoffs
│   ├── reaction_rates.py   # ClusterRates — absorption/emission coefficients
│   ├── rate_equations.py   # ClusterDynamicsODE — builds and evaluates ODE RHS
│   ├── simulation.py       # ClusterDynamicsSimulation — orchestrator
│   ├── post_process.py     # Derive mean sizes, densities, size distributions
│   └── visualization.py   # Standardised plot functions
├── cpp_utils/              # C++ CVODE solver (placeholder — adapt after Python validated)
├── input/
│   └── input_parameters.xlsx   # 3-sheet Excel (same layout as EuroferExperiments)
└── output/                 # Timestamped run directories (gitignored)
```

## Input Excel sheets

Same 3-sheet structure as EuroferExperiments, with additional entries in `Model_Parameters`:

| Parameter | Notation | Description |
|---|---|---|
| Max vacancy cluster size | `N_v` | Tracks $f_v(1 \ldots N_v)$ |
| Max interstitial cluster size | `N_i` | Tracks $f_i(1 \ldots N_i)$ |
| Loop nucleation cutoff | `N_loop` | $f_i(n > N_{loop})$ feeds loop density |
| Capture radius (cluster) | `r_cap_cd` | Capture radius for cluster–cluster reactions |

## Running

```python
import sys
sys.path.insert(0, '../')   # from ClusterDynamics/code/
from ClusterDynamics.py_utils.simulation import run_cluster_dynamics_simulation

CONFIG = {
    't_begin': 1e-1, 't_end': 1e6, 'n_points': 500,
    'method': 'LSODA', 'rtol': 1e-6, 'atol': 1e-20, 'log_time': True
}
results, sim = run_cluster_dynamics_simulation(CONFIG)
```

## Key references

- Ghoniem & Cho (1974): cluster dynamics framework
- Bullough & Hayns (1975): point defect rate equations
- Stoller (1996): void swelling in ferritic steels
- Souidi et al. (2006): cascade-produced cluster distributions in Fe
- See `../Docs/Formulation/` for full rate-equation derivations
- See `../Docs/Database/` for EUROFER97 experimental data
