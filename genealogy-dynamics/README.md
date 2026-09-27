# Genealogy Dynamics — two-sex cluster-dynamics model of the U.S. population

A standalone prototype that adapts the graph-based (RAG) cluster-dynamics
method to human population genealogy: instead of defect clusters binned by
size, people are binned by **pedigree depth (cluster)** in two populations,
**paternal** ($P_{\rm pat}$) and **maternal** ($P_{\rm mat}$), and the
distribution evolves under immigration (source), cross-population
mating/births (transitions), and death/emigration (sinks).

Two bloodline rules are implemented (see "Bloodline rules" below);
everything downstream (kernels, master equation, solver, plots) is written
against a swappable child-cluster rule, so the comparison is apples-to-apples.

> **Modes.** `run_case2.py` / `run_case2_religious.py` project 2025–2050
> from the Case-1 2025 state with projected inputs (NOT forecasts — inputs
> are assumed/projected, see `data/projection_inputs.py`).
> `run_hindcast.py` replays
> 1980–2024 with real historical rates (see "Hindcast" below); only
> assortativity, the 1980 cluster split, and the immigrant sex split remain
> assumed (emigration is calibrated). `run_colonial.py` replays 1650–2025
> with anchored historical rates and produces per-population cluster
> heatmaps (see "Colonial run" below). Case 1 = 1650–2025 historical;
> Case 2 = 2025–2050 projection.

## The RAG mapping

| Cluster dynamics (RadCluster) | Genealogy dynamics (this app) |
|---|---|
| Cluster of size *n* (vertex) | Cluster state $(d,\alpha)$: pedigree depth $d=0..n$ × population $\alpha\in\{{\rm pat},{\rm mat}\}$ (vertex); $C_0$ = foreign-born; meaning of $C_1..C_n$ depends on the bloodline rule |
| Cascade production (source) | Immigration $I$ splits $\sigma/(1-\sigma)$ into $(C_0,{\rm pat})$, $(C_0,{\rm mat})$ |
| Monomer absorption (growth n → n+1) | Births: mating flux $J_{qr}$ over paternal×maternal pairs; child cluster from the rule $g(q,r)$; births split $s/(1-s)$ between $(C_g,{\rm pat})$, $(C_g,{\rm mat})$ |
| Rate kernels (`reaction_rates.py`) | Mating weights $m^{\rm pat}_q, m^{\rm mat}_r$ (assortativity $\alpha$), fertility $\beta^\alpha_d$, mortality $\mu_\alpha$, emigration $\varepsilon_\alpha$ |
| Fixed sinks | Death + emigration, per cluster state |
| `input_parameters.xlsx` | Assumed/projected input methods in `data/projection_inputs.py` / `genealogy/kernels.py` |
| Provenance-stamped `output/` | `outputs/<timestamp>_<scenario>/` with `provenance.md`, plus `outputs/<timestamp>_comparison/` for rule overlays |

## State and master equation

State $\mathbf{c}(t)\in\mathbb{R}^{2(n+1)}_{\ge 0}$ (millions), paternal block
then maternal block; $n$ defaults to 14 (`GENEALOGY_MAX_DEPTH` env var):

$$\frac{d\mathbf{c}}{dt} = \mathbf{P}(t) + \mathbf{S}\,\mathbf{J}(\mathbf{c}) - \mathbf{D}(t)\,\mathbf{c}$$

- **Source** $\mathbf{P}(t)$: immigration $I(t)$ splits $\sigma/(1-\sigma)$
  ($\sigma=0.5$ assumed) into $C_0$ of each population.
- **Loss** $\mathbf{D}(t)$: diagonal, $\mu_\alpha(t)+\varepsilon_\alpha$
  (uniform across sexes by default; sex-specific rates supported).
- **Flux** $J_{qr}(\mathbf{c}) = B(\mathbf{c})\,[\alpha\,\delta_{qr}
  m^{\rm pat}_q + (1-\alpha)\,m^{\rm pat}_q m^{\rm mat}_r]$ over ordered
  paternal×maternal pairs, with
  $B(\mathbf{c}) = \min(\sum_d\beta^{\rm pat}_d c_{d,{\rm pat}},
  \sum_d\beta^{\rm mat}_d c_{d,{\rm mat}})$ (marriage-function style).
- **Stoichiometry** $\mathbf{S}$: $S_{(g,{\rm pat}),(q,r)} = s$,
  $S_{(g,{\rm mat}),(q,r)} = 1-s$ where $g = g(q,r)$ is the child-cluster
  rule; $s = 0.512$ (measured NCHS sex ratio at birth, ≈1.05 male/female).
  Every column of $\mathbf{S}$ sums to 1, so births conserve people exactly.

Exact total balance:
$\frac{d}{dt}\sum_{d,\alpha} c_{d,\alpha} = I + B - \sum_{d,\alpha}
(\mu_\alpha+\varepsilon_\alpha)\,c_{d,\alpha}$ (verified to $10^{-15}$).

**Symmetry reduction.** With $s=\sigma=0.5$, uniform rates across sexes,
and a symmetric initial condition, the two populations stay identical and
the cluster totals $T_d = c_{d,{\rm pat}}+c_{d,{\rm mat}}$ satisfy exactly
the earlier one-population equations. A 200-year side-by-side integration
agrees to $4.3\times10^{-9}$ M (see `test_symmetry.py`).

## Inheritance rules

**1. `shallowest_inheritance`.** A child inherits the *shallower*
parental bloodline: $g(q,r) = \min(n, 1+\min(q,r))$. Every ancestor on
both sides must be U.S.-born; one immigrant anywhere breaks the line,
and $C_n$ replenishes only through deep×deep intermarriage.

**2. `deepest_inheritance`.** A child inherits the *deeper* parental
bloodline: $g(q,r) = \min(n, 1+\max(q,r))$. One deep parent suffices to
keep the lineage deep; pedigree depth ratchets upward and accumulates
at the $C_n$ cap.

The two rules agree when $q=r$ and when both parents are foreign-born
($g(0,0)=1$); they differ exactly when the parents sit at different
depths. They are the extreme interpretations of pedigree depth, and
bracketing them brackets the plausible answers.

## Rule comparison (colonial 1650–2025, verified 2026-09-27)

Two-sex model, $n=14$, identical demography; deepest-cluster ($C_{14}$)
share of the 2025 population:

| α | $C_{14}$ share (shallowest) | $C_{14}$ share (deepest) |
|---|---|---|
| 0.0 (random-mating null) | 0.00% | 71.88% |
| 0.8 (central) | 0.45% | 18.35% |
| 0.9 (strong-endogamy bound) | 1.46% | 9.44% |

The totals are identical (259.0 M): the inheritance rule redistributes
mass without changing it. Mechanism: shallowest inheritance recycles
deep lineages downward through deep–shallow pairings; deepest
inheritance never pulls a lineage down, so depth ratchets upward to the
cap. Higher assortativity slows the ratchet (fewer mixed-depth
matings), lowering the cap share under the deepest rule.

## Hindcast 1980–2024 (real historical rates)

`run_hindcast.py` replays 1980→2024 with time-dependent $I(t)$, $\beta(t)$,
$\mu(t)$ from public sources (both rules; totals are rule-independent to
$9\times10^{-4}$ M). Rates are linearly interpolated between Jan-1 anchors.

### Data sources (all public, no credentials; accessed 2026-09-25)

| # | Series | Years | Source |
|---|---|---|---|
| 1 | Persons obtaining lawful permanent resident status (fiscal year → calendar year t) | 1980–2024 | DHS Office of Homeland Security Statistics, Yearbook Table 1 |
| 2 | Crude birth rate (per 1,000) → β(t) = 2·CBR/1000 | 1980–2024 | World Bank WDI SP.DYN.CBRT.IN |
| 3 | Crude death rate (per 1,000) → μ(t) = CDR/1000 | 1980–2024 | World Bank WDI SP.DYN.CDRT.IN |
| 4 | Total fertility rate (cross-check only) | 1980–2024 | World Bank WDI SP.DYN.TFRT.IN |
| 5 | Total population (validation target) | 1980–2024 | World Bank WDI SP.POP.TOTL |
| 6 | International migrant stock, 5-yr anchors (C0 validation) | 1990–2024 | World Bank WDI SM.POP.TOTL |
| 7 | Decennial foreign-born (1980 IC) | 1980/1990/2000 | U.S. Census Bureau, WP-81 |
| 8 | Unauthorized immigrant stock, Jan-1 (stock→flow for I(t)) | 1980–2024 | DHS OHSS, Estimates of the Unauthorized Immigrant Population |
| 9 | Colonial population estimates (1650–1780) | 1650–1780, decennial | U.S. Census Bureau, Historical Statistics, Ser. Z 1-19 |
| 10 | Decennial census totals (validation anchors) | 1790–2020 | U.S. Census Bureau |
| 11 | White crude birth rates (pre-1960 CBR anchors) | 1800–1950, decennial | Haines (EH.net encyclopedia, Table 1) |
| 12 | Crude death rate anchors (pre-1960 CDR) | 1800/1850/1900/1930 | Blodget/Haines; Death Registration Area |
| 13 | Sex ratio at birth (→ s = 0.512) | annual, 1940–2024 | NCHS/CDC, National Vital Statistics Reports, Births: Final Data (Vol. 75, No. 2 for 2024); ratio 1,046–1,059 males per 1,000 females |

Reproduce: `python3 data/fetch_data.py && python3 data/build_inputs.py`
(raw cache in `data/raw/`; clean tables in `data/inputs.csv` (1650–2024,
with `quality`/`quality_note` columns), `data/validation.csv`,
`data/validation_long.csv`, `data/initial_condition_1980.csv`,
`data/emigration_calibration.txt`). The Census API now requires a key,
so Census-published aggregates are used instead.

### Rate construction

- **Immigration:** I(t) = DHS LPR count / 1e6 M/yr (fiscal year mapped to
  calendar year) **plus** an unauthorized inflow estimated from DHS OHSS
  stock reports: I_unauth(t) = max(0, S(t)−S(t−1)) + μ(t)·S(t−1).
  Total unauthorized gross inflow 1980–2024 ≈ 15.2M (approximate).
  Two documented biases: (i) stock declines are floored, so the proxy misses
  gross entries during decline years; (ii) DHS LPR counts include adjustments
  of status, so some legalized unauthorized immigrants are counted in *both*
  series (modest double count, 2000–2024).
- **Fertility:** the model's birth flux is $B = \min(\sum\beta^{\rm pat}
  c_{\rm pat}, \sum\beta^{\rm mat} c_{\rm mat})$, so a *uniform*
  β(t) = 2·CBR(t)/1000 gives $B \approx$ CBR/1000 × total — total births
  track the crude birth rate by construction (up to the min with the
  maternal side; see the marriage-function note below).
- **Mortality:** μ(t) = CDR(t)/1000, uniform across sexes and clusters
  (no age structure). Pre-1960 CDR anchored at 1800 ~24, 1850 ~23,
  1900 17.2, 1930 11.3, linear between; 1650–1799 assumed 28/1000.
- **Pre-1820 immigration:** no DHS series exists, so a *demographic residual*
  is used: for each decade between population anchors, the constant annual
  inflow that carries the anchor forward under the model's own vital rates
  (this captures all net migration, including enslaved imports, which DHS
  misses by construction; net migration is floored at zero per decade).
- **Still stylized:** assortativity α = 0.6, immigrant male fraction
  σ = 0.5, and the 1980 cluster split (C0 = 14.08M is Census data; total
  227.225M from WDI; split evenly M/F). No annual US emigration series
  exists, so ε is a **calibrated constant** (see below); no dataset measures
  pedigree depth.
- **Quality flags:** every row of `data/inputs.csv` carries
  `quality` ∈ {measured, interpolated, assumed} = weakest link across the
  immigration and vital-rate sources for that year, plus a `quality_note`.

### Emigration calibration

`data/calibrate_emigration.py` grid-searches ε ∈ [0, 0.004] (model totals
are rule-independent, so the shallowest rule is used) to minimize the sum
of squared errors of total population vs World Bank SP.POP.TOTL over
1980–2024, after adding the unauthorized inflows.

- **Calibrated ε = 0.00075/yr** (≈210k/yr at current population).
- The calibrated value is hardcoded in `run_hindcast.py` (EPS) and reused
  for the colonial residual in `data/build_inputs.py`; the calibration
  record lives in `data/emigration_calibration.txt`.

### Validation results

| Check | Result |
|---|---|
| Total population 2024 | model vs actual **−2.87M, −0.8%**; max abs err over 1980–2024 = 2.87M (2024) |
| Rule-independence of totals | max |shallowest − deepest| total = $7\times10^{-4}$ M |
| Exact total balance (time-dependent) | max deviation $1.3\times10^{-15}$ |

**Marriage-function note.** The $\min(\cdot,\cdot)$ marriage function
interacts with the measured $s=0.512 \neq 0.5$: because 51.2% of births are
male, the maternal side is always the binding constraint, so realized births
run ≈2% below the symmetric rate. Over the 44-year hindcast this barely
registers (−0.9M vs the symmetric limit), but compounded over the 375-year
colonial run it leaves the 2025 total at 300M vs 336M in the symmetric
($s=\sigma=0.5$) limit, and the anchor fit degrades to −13.4% worst error
since 1960 (year 2000). A less harsh two-sex marriage function (e.g.
harmonic mean) would close most of this gap; the $\min$ form is kept because
it is the specified kernel, and the gap is reported rather than hidden.

## Colonial run 1650–2025 (cluster distributions)

`run_colonial.py` replays 1650→2025 with the anchored historical rates from
`data/inputs.csv` (2025 holds 2024 values), both bloodline rules, calibrated
ε = 0.00075/yr, α ∈ {0.0, 0.6, 0.9}, and writes
`outputs/<timestamp>_colonial/` containing:

- `heatmap_<rule>_pat.png` / `heatmap_<rule>_mat.png` — **the four
  graphs**: per-population cluster shares at 25-year snapshots (1650, 1675,
  …, 2000, 2025), x = snapshot year, y = $C_0$–$C_n$, color = population share.
- `results_<rule>.csv` — annual totals + per-population cluster shares
  (α=0.6); `results_<rule>_alpha{a}.csv` — α variants.
- `summary_all.csv` — 2025 totals, paternal/maternal split, $C_n$ shares,
  exact-balance deviations for every rule×α.
- `colonial.md` — key numbers + plain-language reading guide.

### Initial condition and spin-up

1650: the entire colonial population (50,368 people, Census HSUS Z 1-19)
in $C_0$, split evenly M/F; $C_1..C_n$ = 0. This is a **spin-up
approximation**, not a claim about 1650 pedigrees — treat roughly the first
100 years as transient while pedigree depth builds from the cold start.

### 2025 results

| rule | α=0.0 | **0.8 (central)** | 0.9 |
|---|---|---|---|
| shallowest inheritance $C_{14}$ share | 0.00% | **0.45%** | 1.46% |
| deepest inheritance $C_{14}$ share | 71.88% | **18.35%** | 9.44% |

Paternal share of the 2025 total: **51.0%** in every run — near 50/50, as
expected from the measured $s=0.512$ sex ratio at birth and the assumed
$\sigma=0.5$ immigrant split. The 2025 total is 300.2M in every run
(rule-independent); the symmetric ($s=\sigma=0.5$) limit gives 336.2M —
see the marriage-function note above.

**Timescale note.** The model has no age structure, so lineages climb the
cluster ladder faster than real ~27-year generations (newborns reproduce
immediately in the ODE). Read clusters as an *ordinal depth index*; the
robust outputs are the qualitative patterns and the rule contrast, not the
exact timing of when a cluster fills.

## How to run

```bash
pip install -r requirements.txt
python3 run_case2.py              # Case 2: 2025-2050 projection, two-sex x both rules
python3 run_case2_religious.py    # Case 2: 2025-2050 projection, religious model
python3 data/fetch_data.py && python3 data/build_inputs.py   # real data
python3 data/calibrate_emigration.py  # grid-search eps on 1980-2024 totals
python3 run_hindcast.py         # hindcast 1980-2024 x both rules + validation
python3 run_colonial.py         # 1650-2025 x both rules, alpha sweep
python3 plot_timeseries.py      # total-share time series from the colonial run
python3 plot_distributions.py   # final cluster distributions
python3 test_symmetry.py        # two-sex (s=sigma=0.5) vs one-sex reference
```

`run_case2.py` projects 2025–2050 × two rules (`shallowest_inheritance`,
`deepest_inheritance`) × α = 0.0/0.8/0.9 from the Case-1 2025 handoff state,
each writing
`outputs/<YYYYMMDD_HHMMSS>_<scenario>/` with `results.csv` (yearly
per-population clusters, totals, paternal/maternal split), plots, and
`provenance.md`; plus `outputs/<YYYYMMDD_HHMMSS>_comparison/` with the
deepest-cluster-share overlays and `summary.md`.

## Sanity checks (run automatically)

1. **Non-negativity:** min over all cluster states/times ≥ 0.
2. **Pure decay:** with I = 0 and β = 0, the numerical solution must match
   exp(−(μ+ε)t) to solver tolerance.
3. **Exact population balance:** numerical d(total)/dt must equal analytic
   I + B − Σ(μ+ε)·c to floating-point precision (not just solver tolerance).
4. **Symmetry reduction** (`test_symmetry.py`): with s=σ=0.5 the two-sex
   cluster totals reproduce a one-sex reference to < 1e-6 M.

## Limitations (prototype)

- Well-mixed national model: no regions, no age structure. (No age
  structure means lineages climb cluster depth faster than real ~27-year
  generations — clusters are an ordinal depth index.)
- Two-sex structure but uniform vital rates across sexes by default;
  sex-specific β/μ/ε are supported but uncalibrated, as is the immigrant
  sex split σ = 0.5 (assumed).
- The $\min(\cdot,\cdot)$ marriage function is harsh when $s \neq 0.5$
  (see note above); a harmonic-mean form is a candidate refinement.
- Case-2 inputs are assumed/projected (`data/projection_inputs.py`: β/μ OLS
  trends, immigration frozen at the 2024 measured 1.705 M/yr, disaffiliation
  held at the calibrated rate); the hindcast/colonial runs use the DHS
  LPR series plus a net-stock-based unauthorized inflow proxy (floored
  declines; adjustments-of-status double count — see above).
- Emigration is a calibrated constant (no annual US series exists).
- Pre-1960 vital rates are anchored interpolations (Haines CBR; CDR anchors
  at 1800/1850/1900/1930) and 1650–1799 values are assumed; pre-1820
  immigration is a demographic residual against population anchors. Every
  input row carries a quality flag (`data/inputs.csv`).
- The 1650 all-in-$C_0$ initial condition is a spin-up device; the first
  ~100 years are transient.
- Child-cluster rule is a modeling choice (two definitions implemented), not empirical.
- Next steps: regional/cultural compartments (the RAG extends naturally to
  (region, cluster) vertices), age structure, CPS parental-nativity
  calibration of α, better emigration data if a defensible series appears.
