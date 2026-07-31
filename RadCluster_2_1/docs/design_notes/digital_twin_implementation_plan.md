# RadCluster_2_1 Digital Twin — Implementation Plan

**Source specification:** `docs/design_notes/Ferritic_Martensitic_Steel_Microstructure.pdf`
(Ghoniem 2026) — §2.4 (UQ workflow), Appendix E (input parameters),
Appendix F (parameter-identification algorithm), Appendix G (experimental
microstructure database), Appendix H (digital-twin construction), and in
particular Tables H.37 (target measurements) and H.38 (prior ranges).

**Target code:** `RadCluster_2_1/` at commit `b17fe85`.

**Status:** plan, with the revision-3 input-path fixes **implemented** (§4/T0.5
tasks b, c, d). Everything else below is still plan only.

**Revision 3 (2026-07-30)** — author directives, implemented rather than
specified:

| Directive | Outcome |
|---|---|
| Add `lambda`, `A_void_0`, `phi_max_junc`, `loop_net_*` to `input_parameters.xlsx` at nominal values | `lambda`, `A_void_0`, `phi_max_junc` were **already present** — the defect was that the first two were never *read*; now they are. The `loop_net_*` family + `LOOP_NETWORK_LOSS` **added** to the workbook and to `create_excel.py`, at nominals equal to the code defaults. §2.4 |
| Read an input di-interstitial binding energy, overriding the `A_111 = 3` default | `dissociation!E_b_i2 = 0.80 eV` added; overrides `A_111` (and rescales `A_100` to preserve the character ratio). §2.2(c) |
| The `T*` issue may be resolved in LF runs | §2.3 item 3 rewritten — `T*` is no longer a pre-Tier-1 blocker; Tier 2 screening decides it. |

Files touched: [binding_energies.py](../../py_utils/binding_energies.py),
[reaction_rates.py](../../py_utils/reaction_rates.py),
[simulation.py](../../py_utils/simulation.py),
[create_excel.py](../../py_utils/create_excel.py),
[input/input_parameters.xlsx](../../input/input_parameters.xlsx).

**Revision 2 (2026-07-30)** — author directives incorporated. Changes from
revision 1, each carried through §0–§9:

| # | Directive | Where it lands |
|---|---|---|
| 1 | HF tier is bounded by the **existing neutron data envelope**: `T ≤ 400 °C`, `d ≤ 30 dpa`, and a hard **7200 s wall-clock cap per run** — whichever binds first | §0-1, §4/Tier 0.4, §4/Tier 3, §5, §8-1 |
| 2 | `Z_i^loop` gets its own rate-constant entry **and its own workbook key** | §2.1 #16, §2.2, §4/T0.5 |
| 3–4 | `A_V`, `B_V` are **withdrawn** and replaced by the **blended-capillary** parameters of the reference document, Eq. (B.14)–(B.15): `γ_s`, `E_b^v(2)`, `λ` | §2.1 #10–11, §2.2 |
| 5 | Fusion prediction **drops `η_H`** — decision, no longer an assumption | §1, §4/Tier 7, §8-7 |
| 6 | Prior ranges revised: `f_V^cl ∈ 0.2–0.7`, `f_I^cl ∈ 0.05–0.25`, `E_m^I ∈ 0.3–0.8 eV`, `E_m^V ∈ 0.5–1.0 eV`, `T* ∈ 550–610 K` | §2.1, §2.3 |
| 7 | **Di-interstitial binding energy** `E_b^i(2) ∈ 0.6–1.2 eV` added to `θ` | §2.1 #14, §2.2 |
| 8 | **Loop-conversion and network block** added to `θ`: `w_c`, `E_a0_conv`, `dH2_conv`, `φ_max`, `χ` — ranges taken from §3 of the reference document (*Comparison Between Modeling and Experiment*) | §2.4 |

Note that directives 3 and 4 are sequential edits to the same slot: `A_V`/`B_V`
are **not** added to the code or the workbook. The net instruction is the
substitution recorded in §2.2.

---

## 0. Executive summary

The appendix specifies a complete Bayesian digital twin: Sobol design → LF/HF
RadCluster ensembles → autoregressive multi-fidelity GP → posterior over an
18-parameter vector → sequential assimilation → fusion-condition prediction with
credible intervals. RadCluster_2_1 already supplies most of the *forward* side of
that pipeline — a driveable simulator, conservation diagnostics that are exactly
the `Φ_phys` penalties of Eq. (H.20), a TEM-visible observation operator buried
inside `visualization.py`, and a fitted target database with prediction bands.
What is missing is the *inverse* side: there is no parameter I/O layer, no
ensemble driver, no emulator, no calibration, and no cost model.

Three findings from reading the code against the appendix change the shape of the
plan and should be settled before any ensemble is launched:

1. **The high-fidelity tier is now defined by the data envelope, not by
   Appendix F's 100 dpa.** Appendix F asks for `d_HF ≈ 100 dpa` at full bin
   resolution, `N_HF ≃ 3p–8p` = 54–144 runs. The measured cost of the current
   `bin_moment` + `full_system` + Woodbury path (`output/sweep_summary.log`) is
   117 s to 3 dpa at 550 °C but **4892 s to reach only 0.05 dpa at 450 °C**, and
   6047 s to 0.5 dpa at 400 °C. The stiff mid-temperature window — which is
   exactly where the EUROFER97 database lives (250–350 °C) — is two to three
   orders of magnitude short of 100 dpa. **The HF tier is therefore scoped to
   the envelope the experiments actually occupy — `T ≤ 400 °C` and
   `d ≤ 30 dpa` — under a fixed 7200 s per-run wall-clock cap, and each run
   terminates at whichever of the two limits it reaches first.** There is no
   scientific loss: nothing in the neutron database lies outside that box
   (the extremes are Materna-Morris 400 °C / 17.2 dpa and the 335 °C / 32 dpa
   point, which sits marginally over and is handled as the one extrapolation
   case). Tier 0 still measures cost and reach, but its output is now *which
   conditions complete inside the box*, not a redefinition of what "HF" means.

2. **Three parameter slots in Eq. (H.8) required an author decision; all three
   are now settled.** `Z_i^loop` had no workbook key and was aliased to `Z_i` in
   [reaction_rates.py:132](../../py_utils/reaction_rates.py#L132) behind an
   explicit `TODO(Stage3)` — it now gets its own entry in the interstitial
   absorption rate `K_loop` and its own `reactions`-sheet key, so loop growth is
   independently controllable. `A_V` and `B_V` never existed in the code, whose
   void binding is the **blended capillary model** of the reference document,
   Eqs. (B.14)–(B.15) — they are withdrawn from `θ` and replaced by that model's
   own three parameters `(γ_s, E_b^v(2), λ)`. A di-interstitial binding energy
   `E_b^i(2) ∈ 0.6–1.2 eV` is added as the small-`n` anchor of the SIA loop
   binding law. See §2.2 for the code-side work each implies.

3. **The nominal workbook still sits outside its own priors, and the revised
   ranges move which parameters are affected.** Under the revision-2 ranges of
   §2.3 the conflicts are `f_cl_v` (workbook 0.05, prior 0.2–0.7), `T*`
   (workbook 723 K, prior 550–610 K), and the SIA loop binding amplitude
   `A_111 = 3.0` in the workbook, which implies `E_b^i(2) = 2.29 eV` — nearly
   twice the top of the new 0.6–1.2 eV band. A prior that excludes the code's
   own default is a defect in one of the two, and calibration will silently
   inherit whichever is wrong. §2.3.

Everything else in Appendix H maps cleanly onto existing machinery. The plan
below is organised as seven tiers of runs (§4), each with a named deliverable,
and closes with a traceability table (§7) showing where every concept in
Appendix H — the posterior, the discrepancy term, the Wasserstein distance, the
Sobol indices, the acquisition function, the acceptance criteria — becomes a
concrete artefact on disk.

---

## 1. What already exists

| Appendix H requirement | Existing implementation | Gap |
|---|---|---|
| Forward map `M(θ; ξ)`, Eq. (H.4) | `RadClusterSimulation.run()` / `run_adaptive()` ([simulation.py:792](../../py_utils/simulation.py#L792)) driving `build/Release/solver.exe` via `cpp_bridge.py` | None — but see cost, §4.1 |
| State `x(d,T;θ)`, Eq. (H.1) | `results['y']` + `rho_net` + `J_*_fixed` cumulative sink integrals | Complete; `ΔJ_g^fix` is `J_He_sink` |
| Observation operator `H[x]`, Eq. (H.4) | `plot_number_densities_tem()` / `plot_mean_sizes_tem()` ([visualization.py:733,795](../../py_utils/visualization.py#L733)) — TEM cutoff `_N_MIN_TEM = 10` atoms, loop `d = 2√(nΩ/πb)`, cavity `d = 2(3mΩ/4π)^⅓` | **Buried in the plotting layer.** Must be lifted into a testable `extract_observables.py`. §3 |
| Observables `N_ℓ, d̄_ℓ, N_c, d̄_c, f₁₁₁, f₁₀₀, S`, Eq. (H.5) | All present in the `results` dict from `calculate_derived_quantities()`: `N_loops`, `N_loops_111`, `N_loops_100`, `mean_n_i`, `mean_n_111`, `mean_n_100`, `N_voids`, `mean_n_v`, `swelling`, `f_111_loop` | Only the TEM filter and nm conversion are missing at this level |
| Size distributions `p_ℓ(R)`, `p_c(R)` | `reconstruct_distribution()` in `bin_moment_rates.py`; stair-rendered per-bin density in `visualization.py` §8 | Needs to be emitted as a normalised histogram on the *experimental* bin edges |
| Conservation penalties `δ_FP`, `δ_g`, Eq. (H.20) | `results['delta_FP']`, `delta_FP_sia`, `delta_FP_vac`, `delta_He` | None — read directly |
| Experimental targets `D_exp`, Eq. (H.9) | `Eurofer_micro_database/FerriticSteelsMicroData.ipynb` → `DB` (68 rows), `LOOPFRAC` (26 rows), `META`/`HIST`/`FITS` (6 histograms); bands in Table G.36 | Needs a stable machine-readable export, §3.2 |
| Log-space residual, Eq. (H.12) | Fit coefficients `(a, c_g, b, s)` per panel in `figures/experiments/fit_summary.txt` | Hard-code Table G.36 into `targets.yaml`; do not re-fit at calibration time |
| Prior families, §H.2 | — | New (`parameters.yaml`) |
| Sobol design, GP, calibration, acquisition | Partial precedent only: `codes/Python_Testing/calibrate.py` (11-param LHS + threshold screening), `param_sweep.py`, `sweep_loop100_and_network.py` | New (`surrogate.py`, `calibrate.py`, `active_refinement.py`) |
| Cost model `C_CPU(θ, ξ)` for Eq. (H.28) | — | New; Tier 0.2 |

**Hydrogen is out of scope — decided.** Eq. (H.32) asks for predictions at
`η_H = 40–50 appm/dpa`, but RadCluster_2_1 has **no hydrogen population** —
`grep` finds only helium (`G_He_r`, `c_h`). Revision 2 settles this: **`η_H` is
dropped from the fusion condition vector `ξ`, and no hydrogen species is added.**
The fusion prediction of Tier 7 is helium-only and is labelled as such in every
artefact it produces. This is a scoping decision, not a deferral — `predict_fusion.py`
must not accept an `eta_H` key at all, so the omission cannot be silently
reintroduced later. Table H.37/H.38 need no change (neither contains `η_H`); the
erratum is against Eq. (H.32) alone.

---

## 2. Parameter vector: audit and mapping

### 2.1 The mapping table

`θ` of Eq. (H.8) **as revised**, against the actual workbook keys. Sheet names are
the `InputData` attribute names, which are the keys `run_ensemble.py` will write
to. Ranges in **bold** are revision-2 overrides of Table H.38; rows marked ★ are
new parameters not in Eq. (H.8) at all.

**Production**

| # | Symbol | Workbook key | Sheet | Range | Prior | Notes |
|---|---|---|---|---|---|---|
| 1 | `η_FP` | `eta` | `production_fission` | 0.20–0.35 | TN | nominal 0.30 — OK. **Was unreadable until 2026-07-30**, see below |
| 2 | `f_I^cl` | `f_cl_i` | `production_fission` | **0.05–0.25** | B/U | workbook 0.25, but the *code* ran 0.58 — three-way conflict, §2.3 |
| 3 | `f_V^cl` | `f_cl_v` | `production_fission` | **0.2–0.7** | B/U | workbook 0.05, code ran 0.15 — **both outside prior**, §2.3 |

> **The whole Production sheet was a dead code path until 2026-07-30.**
> `production_rates()` read `eta`, `f_cl_i`, `f_cl_v`, `s_i`, `s_v`,
> `i_cascade`, `v_cascade` from hard-coded module dicts in
> `defect_production.py`; `InputData.production_fission` was loaded, stored and
> written to `provenance.md` but never consulted. Parameters #1–#3 were
> therefore **unsamplable**, and Tier 2 would have returned `S_i = 0` for all
> three. Found and fixed during the anchor scan —
> [`anchor_scan_350C_1dpa.md`](anchor_scan_350C_1dpa.md) §5. Backward
> compatibility verified bit-for-bit.

**Mobility**

| # | Symbol | Workbook key | Sheet | Range | Prior | Notes |
|---|---|---|---|---|---|---|
| 4 | `E_m^I` | `E_m_i` | `energetics` + `diffusion` | **0.3–0.8 eV** | TN | nominal 0.34; **written in two sheets** — set both |
| 5 | `E_m^V` | `E_m_v` | `energetics` + `diffusion` | **0.5–1.0 eV** | TN | nominal 0.67; same |
| 6 | `E_m^He` | `E_m_h` | `energetics` + `diffusion` | 0.06–0.09 eV | TN | nominal 0.06; same |
| 7 | `L_1D` | `L_hat` | `diffusion` | 1–10³ nm | LN | **unit change**: `L_hat = L/a`, so sample `L_hat ∈ [3.5, 3.5×10³]`; nominal 50 |
| 8 | `i_mob` | `i_mobile` | ctor kwarg → `diffusion`+`derived` | 1–100 | cat | constructor argument, not workbook-only |
| 9 | `v_mob` | `v_mobile` | ctor kwarg → `diffusion`+`derived` | 1–5 | cat | same; drives Woodbury rank, so it moves cost too |

**Dissociation energetics** — the `A_V`/`B_V` slot is replaced by the three
parameters of the blended capillary model, Eqs. (B.14)–(B.15) of the reference
document. See §2.2.

| # | Symbol | Workbook key | Sheet | Range | Prior | Notes |
|---|---|---|---|---|---|---|
| 10 | `γ_s` | `gamma_s` | `energetics` (read) / `dissociation` (duplicate) | **1.7–2.3 J/m²** | TN | nominal 2.0; sets `E_b^cont`, Eq. (B.15) |
| 11 | `E_b^v(2)` | `E_b_v2` | `dissociation` | **0.10–0.35 eV** | TN | nominal 0.22 (DFT divacancy); sets the Eq. (B.14) atomistic amplitude |
| 11a ★ | `λ` | `lambda` | `dissociation` | **0.4–0.8 vac⁻¹** | TN | nominal 0.5756 = ln100/8; blend decay. **Currently hard-coded**, §2.2 |
| 12 | `A_I` | `A_111` | `dissociation` | 1.0–3.5 eV | TN | nominal 3.0 — **but see #14**; also set `A_100` |
| 13 | `B_I` | `B_111` | `dissociation` | 0.3–0.7 | TN | nominal 0.3873; also `B_100` |
| 14 ★ | `E_b^i(2)` | `E_b_i2` | `dissociation` (**new key**) | **0.6–1.2 eV** | TN | di-interstitial binding; **reparameterises `A_111`**, §2.2 |
| 15 | `A_HeV` | `E_b_hV_1` | `energetics` | 1.5–3.0 eV | TN | nominal 2.30 |

**Sinks and bias**

| # | Symbol | Workbook key | Sheet | Range | Prior | Notes |
|---|---|---|---|---|---|---|
| 16 | `Z_i^d` | `Z_i` | `reactions` | 1.02–1.15 | U | nominal 1.05 |
| 17 | `Z_i^loop` | `Z_i_loop` | `reactions` (**new key**) | 1.05–1.30 | U | nominal 1.10; **de-aliased from `Z_i`**, §2.2 |

**Loop conversion and network evolution** — ranges from §3 of the reference
document; see §2.4 for the evidence behind each.

| # | Symbol | Workbook key | Sheet | Range | Prior | Notes |
|---|---|---|---|---|---|---|
| 18 | `T*` | `T_star_conv_C` | `reactions` | **550–610 K** | TN | **key is in °C, prior is in K** — `param_io` converts; nominal 723 K, **outside prior**, §2.3 |
| 19 | `K_rec` | `loop_net_K_rec` | `reactions` (**new key**) | 10⁻³–10¹ | LN | live only when `LOOP_NETWORK_LOSS = 1`; default 0.0 |
| 20 ★ | `w_c` | `loop_net_w_c` | `reactions` (**new key**) | **1–200 × `b_111`** | LN | capture width; §2.4 |
| 21 ★ | `E_a0_conv` | `E_a0_conv` | `reactions` | **1.4–2.0 eV** | TN | unary direct-rotation barrier — *the* `f₁₀₀` lever; §2.4 |
| 22 ★ | `ΔH₂` | `dH2_conv` | `reactions` | **0.35–0.70 eV** | TN | Marian gate; nominal 1.0 leaves Mechanism B off; §2.4 |
| 23 ★ | `φ_max` | `phi_max_junc` | `reactions` (**new key**) | **0.1–1.0** | U | junction peak yield; §2.4 |
| 24 ★ | `χ` | `loop_net_chi` | `reactions` (**new key**) | **1–60** | LN | elastic capture range; the channel is identically zero below ≈30, §2.4 |

**Held fixed, not sampled** (recorded in `parameters.yaml` with `active: false`
and a reason, so the choice is auditable):

| Key | Value | Why fixed |
|---|---|---|
| `gamma_a_conv` | 0.02 | Table 11 of the reference: physical value; lowering it walks into the stiffness wall without moving `f₁₀₀` (it scales with `P(n)`, and no population lives at large `n`) |
| `dH_rev_conv` | 0.30 eV | Enters `P_succ` only through `Δ = ΔH₂ − ΔH_rev`; non-identifiable against #22 |
| `n_j_min_junc` | 30 | Marian's junction onset (34–37); a mechanism constant, not a fit knob |
| `nu0_conv` | 10¹³ s⁻¹ | Debye frequency |
| `loop_net_rho_max` | 5×10¹⁴ m⁻² | EUROFER network plateau, measured not fitted |
| `n_ref_conv` | 50 | `ΔF` calibration anchor size; redundant with `T*` |

**Count.** `p = 24` active parameters. This is above the "p ≃ 16–20" of
Eq. (H.2), which is expected — Eq. (H.8) predates both the loop-conversion and
loop→network channels. The Tier-2 screening is what brings it back down; §4/Tier 2
states the expectation that 6–10 survive.

Because `A_100`/`B_100` are separate keys, sampling `A_I`/`B_I` must decide
whether the two loop characters share a binding law. Recommendation: sample
`E_b^i(2)` and `B_111` (see §2.2), derive `A_111`, and carry the `⟨100⟩` values as
fixed *offsets* from the Appendix E values (`A_100/A_111 = 0.9545`,
`B_100/B_111 = 0.9246`), so the relative stability of the two characters — which
is what sets `f₁₀₀` — is not an independent free parameter absorbing the
loop-fraction data. With `E_a0_conv` and `ΔH₂` now in `θ` (#21, #22), that
protection matters more, not less: those two are the intended `f₁₀₀` levers and
must not compete with a free loop-binding ratio.

### 2.2 Parameter slots that required a decision — and the code work each implies

All four items below are **settled** in revision 2. Each carries a concrete
code-side task; together they are Tier 0.5.

#### (a) `Z_i^loop` gets its own entry in the interstitial absorption rate — *directive 2*

[reaction_rates.py:132](../../py_utils/reaction_rates.py#L132) currently reads:

```python
Z_i_loop = float(re.get('Z_i', 1.10))   # loop bias factor ≈ same as Z_i
# TODO(Stage3): Z_i_loop currently aliases Z_i; give it its own Excel key.
```

`Z_i_loop` multiplies the loop capture rate
`K_loop(n) = A_loop · n^{1/2} · Z_i^loop · D_i / Ω^{2/3}`
([reaction_rates.py:184–185](../../py_utils/reaction_rates.py#L184)) and is used
at five sites (185, 261, 264, 275, 480). It is the direct gain on SIA capture by
loops — the dominant loop-growth term — so a dedicated key is exactly the
"control of loop growth" the directive asks for.

**Change:**

1. `reaction_rates.py:132` → `Z_i_loop = float(re.get('Z_i_loop', 1.10))`, drop
   the `TODO`.
2. `create_excel.py` REACTIONS, in the *Dislocation Network Sink* block:
   `('Loop interstitial bias factor', 'Z_i_loop', 1.10, '−', 'Eq. P3_i; independent of Z_i (network)')`.
3. Mirror in `cpp_utils/` — the C++ solver reads the same workbook via
   `cpp_bridge.py`, so the parameter must be threaded through
   `parameters.h` / the EUROFER-97 `rate_kernels` alongside `Z_i`.
4. Regression: with `Z_i_loop = Z_i = 1.10` the run must be **bit-identical** to
   the pre-change baseline. That is the acceptance test.

**This must be done before Tier 2.** `Z_i^d` and `Z_i^loop` are the pair whose
*ratio* sets the loop-versus-network partition of the SIA flux; sampling them as
one parameter collapses that degree of freedom and will bias the posterior on
`η_FP` and `f_I^cl` to compensate.

#### (b) `A_V`, `B_V` are withdrawn; the blended capillary model takes the slot — *directives 3 + 4*

Table H.38 named a "vacancy binding-law amplitude and exponent" `A_V m^{-B_V}`,
by analogy with `A_111 n^{-B_111}` for loops. **No such form exists in
RadCluster, and none will be added.** The reference document's own void model is
the atomistic–continuum blend, Eqs. (B.14)–(B.15):

```
E_b^v(α) = E_b^cont(α) + [E_b^v(2) − E_b^cont(2)] · exp(−λ(α−2))
E_b^cont(α) = E_f^v − Ω·2γ_s / R(α),        R(α) = a (3α/8π)^{1/3}
λ = ln(100)/8 ≈ 0.5756 vac⁻¹   (1 % correction at α = 10)
```

Its free parameters are `(γ_s, E_b^v(2), λ)` — rows #10, #11, #11a of §2.1.
These are the three quantities that actually control vacancy emission at the
small sizes where cavity nucleation is decided, they each have literature bounds,
and two of the three are already workbook keys.

**Erratum against Table H.38:** the rows "Vacancy binding-law amplitude `A_V`
0.5–2.0 eV" and "Vacancy binding-law exponent `B_V` 0.3–0.7" are struck and
replaced by the three rows above.

**Code work — this was not purely a `parameters.yaml` edit. Now IMPLEMENTED.**
[binding_energies.py:38–39](../../py_utils/binding_energies.py#L38) hard-coded
both the decay constant and the atomistic amplitudes as module-level constants:

```python
_A_void = {0: 1.2353, 1: 2.9064, 2: 3.4147, 3: 2.1504, 4: -0.1590}  # eV
_lambda_void = 0.5756   # decay constant [vac^-1], = ln(100)/8
```

and `E_b_void(m, E_f_v, gamma_s, Omega)` took no `λ` argument. So `lambda` and
`A_void_0` **sat in the workbook but were never read** — varying them would have
changed the input file and nothing else, which is the worst possible failure mode
for a calibration (it screens as `S_i = 0`, indistinguishable from a genuinely
inert parameter).

Fixed: `E_b_void` now takes `lambda_void=` and `A_void_0=` keyword arguments,
defaulting to the module constants so every existing caller is unchanged, and
`reaction_rates._precompute` supplies them from `dissociation!lambda` and
`dissociation!A_void_0`. Verified live — perturbing `λ` 0.5756 → 0.40 moves the
vacancy emission array `G_VAC(m=2..6)` from `[0.183, 0.269, 0.343, 0.278, 0.169]`
to `[0.012, 0.010, 0.016, 0.022, 0.023]`, and `A_void_0` moves it independently.
No C++ change is needed: the solver consumes **precomputed emission arrays**
written by `cpp_bridge.py`, so any binding-energy parameter that reaches the
Python kernels reaches the C++ run automatically.

**Two items deliberately left open**, both pre-existing and both out of the scope
of the directives:

- **Form discrepancy.** The code uses `A_void[0] · exp(−λ(m−1))` (amplitude free,
  decay from `m = 1`) while Eq. (B.14) uses
  `[E_b^v(2) − E_b^cont(2)] · exp(−λ(α−2))` (amplitude *derived* from the DFT
  divacancy binding, decay from `α = 2`). Adopting Eq. (B.14) would make
  `E_b_v2` the sampled parameter and `A_void_0` a derived diagnostic, removing a
  redundant degree of freedom from `θ`. Until that is done, sample
  `(γ_s, λ, A_void_0)` and treat `E_b_v2` as inactive — otherwise `E_b_v2` and
  `A_void_0` are two names for one quantity and the posterior will be
  unidentifiable along their ridge.
- **`E_b_bubble` carries no atomistic correction at all** — it is pure capillary
  plus He pressure, so `E_b_bubble(m, ℓ=0) ≠ E_b_void(m)`. The reference
  Eq. (B.18) *does* include `+A(m)e^{−λ(α−1)}` for bubbles. This discontinuity at
  `ℓ → 0` predates the twin and should be resolved before Tier 3, because the
  fusion (Case 1) path passes through low-`ℓ` bubbles continuously.

#### (c) Di-interstitial binding energy enters `θ` — *directive 7* — **IMPLEMENTED**

`E_b^i(2) ∈ 0.6–1.2 eV` (row #14). The small-`n` SIA loop binding is the power
law in [binding_energies.py:256](../../py_utils/binding_energies.py#L256),
blended to the continuum at `n_tr = 25`. **The exponent is POSITIVE** — binding
*increases* with loop size:

```
E_b^fit(n) = A_111 · n^{+B_111}      ⟹      E_b^i(2) = A_111 · 2^{+B_111}
```

> **Erratum against `CLAUDE.md` §10**, which writes `E_b^loop(n) = A_111 n^{−B_111}`.
> The code uses `+B_111`, and the positive sign is the physically correct one
> (larger loops bind SIAs more strongly). Revision 1 of this plan inherited the
> `CLAUDE.md` sign and consequently mis-stated the current di-interstitial
> binding as 2.29 eV; the correct legacy value is **3.92 eV**.

**Reparameterise rather than add a fourth free constant.** Set/sample
`(E_b^i(2), B_111)` and derive

```
A_111 ← E_b^i(2) · 2^{−B_111}
```

`A_111` has no independent measurement while `E_b^i(2)` does (DFT, bcc Fe,
≈0.8 eV), so the prior constrains the quantity the atomistics actually report,
at the size where loop nucleation is decided.

**Implemented as follows** (this is now live in the code, not a plan item):

- `binding_energies.A_111_from_E_b_i2(E_b_i2, B_111)` — the inversion above.
- `reaction_rates._precompute` reads `dissociation!E_b_i2`. When present and
  positive it **overrides** the workbook `A_111`, and **rescales `A_100` by the
  same factor** so the `⟨100⟩/⟨111⟩` amplitude ratio — which sets the relative
  stability of the two loop characters, and hence `f₁₀₀` — is preserved rather
  than drifting as a free rider. A blank, absent or non-positive cell falls back
  to the legacy path with `A_111`/`A_100` used verbatim.
- New `dissociation` row `E_b_i2 = 0.80 eV` in both `create_excel.py` and the
  live `input/input_parameters.xlsx`.
- The override is exact: `E_b_loop_i(2)` returns 0.8001 eV, the 1e-4 offset
  being the continuum admixture `(1−w)` at `n = 2`.

**Measured consequence — a real physics change, and it FAILS the 10 dpa
production gate.** At `T = 573 K`, with `B_111 = 0.3873`:

| | legacy `A_111 = 3.0` | override `E_b_i2 = 0.80` |
|---|---|---|
| `A_111` | 3.0 eV | **0.6116 eV** |
| `E_b^i(2)` | **3.92 eV** | 0.80 eV |
| SIA emission `G_SIA(2)` | 7.3×10⁻²⁷ s⁻¹ | **21.5 s⁻¹** |

The legacy value makes di-interstitials effectively indissociable; 0.80 eV is
what DFT reports for bcc Fe. But a 10 dpa run at the known-good production
configuration shows the substitution **destroys the microstructure** — see
[`validation_10dpa_revision3.md`](validation_10dpa_revision3.md):

| | legacy | `E_b_i2 = 0.80` |
|---|---|---|
| `δ_FP` @ 10 dpa | **9.29×10⁻⁶** | **3.16×10⁻¹** |
| `C_SIA_tot` | 1.27×10²⁶ m⁻³ | 1.42×10²² m⁻³ |
| `mean_n_i` | 205.5 SIAs | 13.9 SIAs (0.92 nm) |
| trajectory | evolves to 10 dpa | frozen after ~2.6 dpa |

**Root cause:** the exponent is positive, so anchoring at `n = 2` rescales the
*whole* small-`n` branch by the same factor (0.204 at `n = 2..10`). Fixing the
di-interstitial costs a factor ~5 of binding at every nucleation-relevant size;
emission outruns capture and loops never grow past ~14 SIAs, contradicting the
3–12 nm EUROFER97 anchors.

**Status: mechanism shipped, value disabled.** The `E_b_i2` row exists in the
workbook but is **blank**, so the default is the validated legacy path. Enabling
it requires re-fitting **both** `A_111` and `B_111` to DFT across small `n` — a
two-point solve targeting `E_b(2) = 0.80`, `E_b(25) ≈ E_f^i` gives
`A_111 ≈ 0.5566`, `B_111 ≈ 0.5233` as a starting point — plus re-validation at
this configuration. **Until then `E_b^i(2)` must not be treated as an
independently samplable parameter**: row #14 of §2.1 is inactive, and `θ` reverts
to sampling `(A_111, B_111)` jointly with the constraint that their implied
`E_b^i(2)` lands in 0.6–1.2 eV.

#### (d) Consequence for `p`

`p = 24` — see the count note at the end of §2.1. Appendix F's `N_HF = 3p–8p`
sizing is applied to `p_act` after Tier-2 screening, not to `p`, so no run count
in Appendix F/H changes.

### 2.3 Nominal values that fall outside their own priors

Re-audited against the revision-2 ranges. Workbook values read from
`input/input_parameters.xlsx` at commit `b17fe85`.

| Parameter | `input_parameters.xlsx` | `CLAUDE.md` §2 | Revision-2 prior | Verdict |
|---|---|---|---|---|
| `f_cl_i` (fission) | 0.25 | 0.58 | **0.05–0.25** | **Three-way conflict, now consequential.** Until 2026-07-30 the code used 0.58 and the workbook's 0.25 was inert, so *every existing result in this repository was produced with 0.58*. With the sheet wired, the default silently becomes 0.25. |
| `f_cl_v` (fission) | **0.05** | 0.15 | **0.2–0.7** | Same — code ran 0.15, workbook says 0.05, prior demands ≥0.2. **All three disagree.** Highest-priority reconciliation. |
| `T_star_conv_C` | **450 °C = 723 K** | — | **550–610 K** (277–337 °C) | Workbook outside by ~120 K. **Deferred to the LF ensemble** — see below. |
| `A_111` → `E_b^i(2)` | 3.0 eV → **3.92 eV** | — | **0.6–1.2 eV** | **RESOLVED** — `E_b_i2 = 0.80 eV` is now a workbook key and overrides `A_111` (§2.2(c)). The legacy figure is 3.92 eV, not the 2.29 eV of revision 2: the loop-binding exponent is *positive*. `create_excel.py`'s `A_111 = 0.7501` gives 0.98 eV, inside the band — so the shipped generator was right and the live workbook's 3.0 was the outlier. |
| `E_m_i` | 0.34 | 0.34 | **0.3–0.8 eV** | Inside — the widened range resolves the revision-1 tightness at the lower bound. |
| `E_m_v` | 0.67 | 0.67 | **0.5–1.0 eV** | Inside. |
| `dH2_conv` | **1.0 eV** | — | **0.35–0.70 eV** | Outside by design: §3 of the reference shows `ΔH₂ = 1.0` leaves Mechanism B **switched off** (`P_succ ∼ 10⁻⁵`). The shipped default is a placeholder, and the prior encodes the finding that it must come down. Set the nominal to 0.55. |
| `gamma_a_conv` | 0.03 | — | fixed at **0.02** | Table 11 of the reference adopts 0.02. Workbook ships 0.03; align it. |

Four distinct actions fall out:

1. **`f_cl_v`.** Reconcile 0.05 / 0.15 / [0.2, 0.7] against Appendix A. This
   matters because `f_V^cl` sets the cavity nucleation rate and therefore `N_c`,
   whose band is `e^{2s} = 7.2` — wide enough that a mis-set prior will not be
   caught by the data. The anchor scan measured its influence: `f_cl_v`
   0.05→0.30 is the **single largest lever on `d̄_ℓ`** (+55 %) and moves `N_c`
   by ×6.
2. **`f_cl_i`.** Reconcile 0.25 / 0.58 / [0.05, 0.25] and **write the chosen
   value to the workbook** — this is now load-bearing, not documentation. The
   scan shows `f_cl_i` 0.25→0.10 gives `d̄_ℓ` +30 % and `N_ℓ` ÷1.29, and is the
   **stiffest** direction in the whole prior (`nfe` 6102 → 11063).
3. **`T*` — deferred to the LF ensemble, by author direction.** No pre-Tier-1
   reconciliation is required. The `⟨100⟩` fraction is the only observable in the
   database that pins `T*`, and the §3 campaign evidence says it is a *weak*
   lever over most of the space anyway (Table 7: `T*` 250 → 180 °C left `f₁₀₀`
   bit-identical, because loops live at `n̄ ≈ 52` and never reach the `ΔF > 0`
   cutoff). Rather than pick a value on thin grounds, **the Tier-2 LF screening
   resolves it empirically**: `T*` is sampled over the full 550–610 K prior and
   the ensemble reports whether `f₁₀₀` responds at all. Two admissible outcomes,
   both informative:
   - `S_i^T(T*) ≈ 0` on the `f₁₀₀` channel → `T*` is inert at the loop sizes the
     model populates. Fix it at the prior mean, record it as screened-out, and
     let `E_a0_conv` (#21) carry the loop-fraction signal. This is the expected
     outcome and is consistent with both §3 and Table 11's remark that `T*` need
     only "lie below the desired crossover".
   - `S_i^T(T*)` large → the `ΔF` gate *is* biting, the workbook's 723 K is
     genuinely misplaced, and the posterior locates it.

   The unit trap remains and is handled in `param_io.py`, not by editing the
   workbook: **the prior is in K, the workbook key `T_star_conv_C` is in °C.**
   `run_ensemble.py` converts on write and asserts the written cell is in
   °C. The workbook default stays at 450 °C until the ensemble says otherwise.
4. **`E_b^i(2)` / `A_111` — resolved, no action.** See §2.2(c). The 3.0-vs-0.7501
   discrepancy is now moot: `E_b_i2` overrides `A_111` in both.

**Action, before Tier 1:** reconcile items 1 and 2 (`f_cl_v`, and the `f_cl_i`
documentation) against Appendix A. Item 3 is deferred to Tier 2 by design and
item 4 is closed. `run_ensemble.py` still gets the startup assertion that every
nominal lies inside its own prior (§8-6) — with `T_star_conv_C` explicitly
whitelisted as a known, deliberate exception until Tier 2 reports.

### 2.4 The loop-conversion and network block — *directive 8*

Eq. (H.8) predates both the `½⟨111⟩→⟨100⟩` conversion channel and the
loop→network loss channel, so it carries only `T*` and `K_rec` from a
five-parameter mechanism set. Section 3 of the reference document
(*Comparison Between Modeling and Experiment*) is a measured calibration campaign
over exactly these knobs, and it is the source for every range below. The
evidence is worth stating, because §3's central finding is that the *obvious*
levers do not work and the non-obvious one does.

**What §3 measured.** All runs: C++ solver, fission, `G = 10⁻⁶ dpa/s`,
`full_CD_fission`, `I = 200`, `i_mobile = 10`, `T* = 250 °C`, `γ_a = 0.02`.

| Lever | Effect on `f₁₀₀(500 °C)` | Reading |
|---|---|---|
| dose 3 → 30 dpa | 0.3249 → 0.3249, bit-identical | `f₁₀₀` is **efficiency-limited, not dose-limited** |
| `T*` 250 → 180 °C | unchanged | the `ΔF > 0` cutoff rises 255 → 493 sizes, but loops live at `n̄ ≈ 52` and never reach it |
| `γ_a` 0.02 → 0.015 | stiffness wall, run dose-starves at 0.27 dpa | `γ_a` scales with the perimeter `P(n)`, so it only touches large `n` where no population exists |
| `E_a0_conv` 1.8 → 1.6 eV | **0.325 → 0.948** | **the effective lever** — size-independent, acts at `n ≈ 10–50` where the loops are; `n̄₁₁₁` collapses 52 → 13.5 |
| `ΔH₂` 1.0 → 0.55 eV | 0.325 → 0.342 | Mechanism B engaged but weak (`P_succ ≈ 2×10⁻²`); +5 % |

**#21 `E_a0_conv` — TN, 1.4–2.0 eV, nominal 1.6.** The unary direct-rotation
barrier. Two independent constraints bracket it: measured `f₁₀₀(T)` at
`E_a0 = 1.6` is 0.158 (350 °C) → 0.948 (500 °C), which is the observed steep
grading; and Marian's MD gives the direct rotation as "in excess of 2 eV", so the
band should not extend above 2.0. Below ≈1.4 the crossover region becomes
numerically intractable (`E_a0 = 1.4` at 500 °C dose-starved at 0.01 dpa), which
is a *sampling* hazard, not just a physics bound — Tier 2 must log the
dose-starvation rate against `E_a0_conv` and treat a systematically failing
region as excluded prior support, per Appendix F step 2.

**#22 `ΔH₂` (`dH2_conv`) — TN, 0.35–0.70 eV, nominal 0.55.** The Marian two-step
gate `P_succ(T) = [1 + exp((ΔH₂ − ΔH_rev)/k_BT)]⁻¹`. Its range is set by an
intrinsic trade-off that §3 tabulates directly — the same combination
`Δ = ΔH₂ − ΔH_rev` sets *both* the magnitude and the temperature sensitivity:

| `ΔH₂` | `P_succ(400 °C)` | `P_succ(500 °C)` | ratio |
|---|---|---|---|
| 1.00 (shipped) | 5.7×10⁻⁶ | 2.7×10⁻⁵ | 4.76× |
| 0.70 | 1.0×10⁻³ | 2.5×10⁻³ | 2.44× |
| 0.55 | 1.3×10⁻² | 2.3×10⁻² | 1.73× |
| 0.40 | 1.5×10⁻¹ | 1.8×10⁻¹ | 1.20× |
| 0.35 | 2.97×10⁻¹ | 3.21×10⁻¹ | 1.08× |

Above 0.70 the channel is effectively off; below 0.35 it is strong but
temperature-blind and merely lifts `f₁₀₀` uniformly. **Mechanism B cannot, by
construction, be the origin of the observed temperature dependence** — which is
precisely why both it and `E_a0_conv` must be in `θ`: §3's conclusion is that
*Mechanism A sets where in temperature the conversion turns on; Mechanism B sets
whether it survives to technologically relevant dose.* A twin with only the unary
channel reproduces the crossover but predicts a spurious `f₁₀₀` decay with dose
(measured: 450 °C peak 0.705 → 0.256 as `n̄₁₁₁` grows 18 → 148); a twin with only
the Marian channels is dose-stable but flat in `T`. The calibration needs both,
and `ΔH₂` is the only free knob in Mechanism B once `ΔH_rev` is fixed.

**#23 `φ_max` (`phi_max_junc`) — U, 0.1–1.0, nominal 0.5.** The junction peak
yield. It scales Mechanism B's magnitude *without* touching its temperature slope,
so it is the one parameter that can break the `ΔH₂` magnitude/sensitivity
trade-off. §3's outstanding-work item (i) is exactly the strong-Marian test
`(ΔH₂ = 0.40, φ_max = 1)`; putting `φ_max` in `θ` lets the calibration run that
test continuously instead of at two points. It needs a workbook key —
`create_excel.py` defines it but the shipped workbook has no such row.

**#20 `w_c` (`loop_net_w_c`) — LN, 1–200 × `b_111`, nominal 50.** The
loop→network capture width. §3's companion result tabulates the measured
response at 450 °C, 10 dpa, `I = 1000`, `χ = 50`, `K_rec = 0`,
`ρ_max = 5×10¹⁴ m⁻²` (`ρ̇ ∝ w_c` verified linear):

| `w_c` | `ρ_net` initial | `ρ_net` final | note |
|---|---|---|---|
| `b_111` (physical) | 5×10¹³ | — | ~0.1 %/dpa; would saturate at ~160 dpa |
| `50 b_111` | 5.00×10¹³ | 1.85×10¹⁴ | ×3.7 |
| `150 b_111` | 5.00×10¹³ | 4.56×10¹⁴ | ×9.1, reaches the plateau |

The physical width gives a real but negligible drift, because the net climb
velocity is small (`v_net ≈ 5.5×10⁻¹³ m/s` — the SIA and vacancy fluxes to the
climbing network very nearly cancel). At `150 b_111` the network density rises an
order of magnitude and saturates at the experimental EUROFER plateau, which is
the intended "loop density saturates with dose" behaviour. §3 states plainly that
"the precise physical `(w_c, K_rec)` pair remains an offline calibration" — that
offline calibration *is* this twin, and it is the reason `w_c` must be sampled
rather than pinned at the amplified demo value. With `K_rec` in `θ` but `w_c`
fixed, the calibration cannot recover the loop-density saturation that objective
(a) of RadCluster_2_1 exists to produce.

**#24 `χ` (`loop_net_chi`) — LN, 1–60, nominal 50.** The elastic capture range.
This is in `θ` for a defensive reason rather than a physical one: §3 records that
the geometric capture switch `P_ℓd` **only opens when `χ · d_loop ≳ L_ℓd =
ρ_net^{−1/2}`**, so with few-nm loops against a ~100 nm network spacing it needs
`χ ≳ 30`, and *at the physical `χ ∼ 1` the gain is identically zero*. Sampling
`χ` and `w_c` jointly is the only way the posterior can distinguish "the channel
is weak" from "the channel is switched off by the geometry of the test grid".
Reduced-`I` grids silently disable it, which makes this a live hazard for the
Tier-2 LF ensemble specifically — **the LF grid chosen in T0.2 must be checked
for a nonzero `ρ̇_net` before Tier 2 is authorised**, or the whole screening will
report `w_c`, `χ` and `K_rec` as inert.

**Operational prerequisite (easy to re-trip).** §3 records two conditions that
must hold before `ρ_net` moves at all, and both are driver-side:

1. The driver must call `sim.run_adaptive()`, **not** `sim.run()` — the
   operator-split `ρ_net` update executes only *between* integration segments, so
   a single-shot run leaves it exactly constant even with `LOOP_NETWORK_LOSS = 1`.
   `run_ensemble.py` must therefore use `run_adaptive` unconditionally and assert
   `LOOP_NETWORK_LOSS == 1` implies more than one segment.
2. `χ · d_loop ≳ ρ_net^{−1/2}` must be satisfiable on the chosen grid (above).

**Workbook keys — audited and now complete.** A full re-audit of
`input/input_parameters.xlsx` corrects revision 2's claim that four conversion
keys were missing:

| Key | Revision-2 claim | Actual | Action taken |
|---|---|---|---|
| `E_a0_conv`, `dH2_conv`, `dH_rev_conv`, `gamma_a_conv`, `nu0_conv`, `T_star_conv_C`, `n_ref_conv` | present | present | none |
| `phi_max_junc`, `sigma_s_junc`, `n_j_min_junc`, `n_loop_min` | **"absent"** | **present** (Reactions rows 41–44) | claim withdrawn |
| `lambda`, `A_void_0` | present but unread | confirmed — present, unread | **now read**, §2.2(b) |
| `LOOP_NETWORK_LOSS`, `loop_net_chi`, `loop_net_w_c`, `loop_net_K_rec`, `loop_net_rho_max`, `loop_net_xi`, `loop_net_n_inc` | absent | confirmed absent from workbook *and* from `create_excel.py` | **added to both** |

The `loop_net_*` block was added on the principle **"the workbook nominal is the
code default"**, so adding the rows changes no result:

| Key | Value | Code default it mirrors |
|---|---|---|
| `LOOP_NETWORK_LOSS` | 0 | off — legacy 2_0 behaviour |
| `loop_net_chi` | 1.0 | physical; *the channel is identically zero here* |
| `loop_net_K_rec` | 0.0 | no recovery |
| `loop_net_rho_max` | 1.0×10¹⁶ m⁻² | runaway guard (EUROFER plateau is ~5×10¹⁴) |
| `loop_net_xi` | 0.0 | off |
| `loop_net_w_c` | **blank** | per-character Burgers vector (dynamic default) |
| `loop_net_n_inc` | **blank** | `i_mobile` (dynamic default) |

The two keys with *dynamic* defaults are shipped blank, and the readers were
hardened so a blank cell falls back to the code default rather than producing
`NaN` — reading them with a bare `float()`/`int()` would have turned an empty
cell into `NaN` and propagated it silently into `Λ_n^net`.

**The §3 trap reproduced, as a regression check.** With the shipped nominals the
channel is verifiably inert, and it opens exactly where §3 says it does:

| `I` | `χ` | `w_c` | `max Λ_net` | sizes with `Λ > 0` |
|---:|---:|---|---:|---:|
| 1000 | 1 | `b_111` | 0 | 0 / 1000 |
| 1000 | 50 | `b_111` | 1.63×10⁻⁹ | 372 / 1000 |
| 1000 | 50 | `150 b_111` | 2.46×10⁻⁷ | 372 / 1000 |
| 2000 | 50 | `150 b_111` | 2.46×10⁻⁷ | 1372 / 2000 |

`Λ ∝ w_c` is linear to 3 digits (150× the width → 150× the rate), matching §3's
"`ρ̇ ∝ w_c` verified linear". Note the first row: at the *physical* `χ = 1` the
channel is identically zero on any grid — which is precisely why T0.4b exists.

---

## 3. Modules to build

Directory `RadCluster_2_1/digital_twin/`, mirroring the Appendix H.11 checklist
one-for-one so the correspondence is auditable:

```
digital_twin/
├── parameters.yaml           # H.11-1: 20 × {symbol, workbook key, sheet, units,
│                             #          nominal, bounds, prior family, active}
├── targets.yaml              # H.11-2: Table G.36 bands + Table H.37 rows
├── experiments.yaml          # H.11-2: per-row ξ_j, z_j, σ_j, g_j from the workbook
├── conditions.yaml           # the run grid of §4 (which ξ at which fidelity)
├── param_io.py               #         θ (unit cube) → InputData mutations
├── run_ensemble.py           # H.11-3: LF/HF driver, conservation gate, cost log
├── extract_observables.py    # H.11-4: the observation operator H[x]
├── surrogate.py              # H.11-5: AR multi-fidelity GP, Appendix F steps 5–10
├── calibrate.py              # H.11-6: posterior (SMC or EnKI)
├── active_refinement.py      # H.11-7: acquisition A = I / C_CPU
├── predict_fusion.py         # H.11-8: posterior predictive at DEMO conditions
└── report/                   # every table and figure of §7
```

### 3.1 `extract_observables.py` — the observation operator

This is the highest-value single module, because every downstream stage consumes
its output and because the appendix's comparison is only meaningful if the
simulated observable is defined the same way the microscopist defined it.

```
H(results, input_data, rate_eq, dose_query) -> dict
```

Decisions that must be made explicitly and recorded in `provenance`:

1. **TEM visibility cutoff.** `visualization.py` uses `_N_MIN_TEM = 10` atoms
   (≈1 nm loop). Real WBDF detection limits are 1–1.5 nm and vary by study;
   Dethloff's black-dot convention moves one EUROFER97 measurement from 72 % to
   45 % (Appendix G.5). Treat `d_min` as a **nuisance parameter** sampled with the
   rest (U[0.8, 1.5] nm) rather than a hard-wired constant. It is cheap: it is
   applied in post-processing, so one simulation yields the whole `d_min` curve.
2. **Diameter conversion.** Loop `d = 2√(nΩ/πb)` with `b = b_111` or `b_100` per
   character; cavity `d = 2(3mΩ/4π)^{1/3}`. Already correct in
   `visualization.py:831–848`; move verbatim.
3. **Which swelling.** `results['swelling']` is the vacancy-inventory `S(t)` of the
   conservation identity (CLAUDE.md §8) — total vacancy content, including
   sub-visible clusters. Table H.37 wants `S = (π/6) N_c d̄_c³` from the
   *visible* cavity population. **These are different numbers.** Emit both, as
   `S_inventory` and `S_visible`, and compare only `S_visible` to the database.
4. **Loop fraction convention.** The code's `f_111_loop` is content-weighted
   (Eq. 56, `Σ n c_n`); most database entries are number fractions from loop
   counting. Emit both `f_100_content` and `f_100_number`; compare `f_100_number`
   to the database and carry the difference between them as a lower bound on
   `σ_model` for that observable.
5. **Size distributions on experimental bin edges.** Reconstruct `c_n`, convert to
   diameter, rebin onto the edges in `META`/`HIST` for the matching condition,
   normalise, return. This is what `W₁` consumes.

Unit test: replay a saved run from `output/20260626_162101_.../results_y.npy`
through `H` and assert the numbers reproduce the corresponding panel in
`plots/`. That guarantees no drift between the twin and the notebook.

### 3.2 `experiments.yaml` / `targets.yaml`

`FerriticSteels_parsed_database.csv` has been removed from the working tree (it
was a stale derived artefact). The twin must therefore **regenerate its target
file from the workbook**, using the same content-marker exec pattern already
proven in `Eurofer_micro_database/make_experiment_figures.py` — that keeps the
twin, the notebook and the manuscript on one parser and makes drift impossible:

```python
CELL_MARKERS = ("XLSX = ", "def parse_density", "def build_database",
                "REFERENCES = [", "def load_size_distributions",
                "def load_loop_fractions")
```

`targets.yaml` carries **Table G.36 verbatim** — the eight `(a, b, c_g, s)`
tuples — as frozen constants with a provenance line. Do not re-fit at calibration
time; the fits are a published result and refitting inside the loop would let the
targets move as the code changes.

**The −60 °C ion shift must not reach the simulator.** `T_eq = T − 60 °C` is a
data-side device for putting ion and neutron points on one axis (Appendix G.1),
and Appendix G already notes it is *exactly collinear* with an ion dummy. The
twin simulates the **actual** `T` and the **actual** `φ̇` of each experiment —
the dose-rate difference is physics the model is supposed to produce, not a
correction to apply. `T_eq` is retained in `experiments.yaml` for plotting only,
and `run_ensemble.py` must assert it is never written to `reactions['T']`.

---

## 4. Sequence of runs

Seven tiers. Each is a gate: do not start tier *n* until tier *n−1*'s deliverable
exists and its acceptance test passes.

### Tier 0 — Numerics, cost, and code fixes (no calibration)

**T0.1 Verified-tolerance ladder.** §2.4 of the paper insists physical inference
must never be contaminated by numerical approximation, and fixes `(rtol, atol)`
for the whole campaign. Establish them:

| Axis | Values | Conditions |
|---|---|---|
| `rtol` | 1e-4, 1e-5, 1e-6 | 300 °C / 400 °C / 500 °C, `1e-6` dpa/s, to 1 dpa |
| `atol` | 1e-20, 1e-25 | same |
| `C_floor` | 1e-15, 1e-25 | same |

**The `rtol` rung is already measured** — see
[`validation_10dpa_revision3.md`](validation_10dpa_revision3.md) §8 and
`output/20260730_160946_deltaFP_study/`. A controlled sweep at fixed temperature,
grid and physics gives the **converged** `δ_FP`:

| `rtol` | `δ_FP` @1 dpa, 623 K | @10 dpa |
|---|---|---|
| 1e-4 | 9.7×10⁻⁶ | 9.0×10⁻⁷ |
| 1e-5 | 1.5×10⁻⁷ | 3.9×10⁻⁶ |
| **1e-6** | **1.2×10⁻⁸** | 3.6×10⁻⁶ |
| 1e-7 | 1.4×10⁻⁷ | 3.9×10⁻⁶ |

**Adopt `rtol = 1e-6` for HF.** It is the optimum, not merely the tightest —
1e-7 is *worse* (roundoff-limited).

**LF may run at `rtol = 1e-4`.** At the new `d_LF = 0.1 dpa` the measured
`δ_FP` is **1.5×10⁻⁵ at 623 K and 1.5×10⁻⁶ at 723 K**, both comfortably inside
the 1e-4 health check of §6. Since conservation is no longer an acceptance
criterion, there is nothing to gain from a tighter LF tolerance. Measured cost
at 623 K to 10 dpa: **309 s at `rtol = 1e-4` versus 516 s at 1e-6**, so the loose
setting is ~1.7× cheaper — on top of the 10× saved by the shorter dose. This is
what makes the reduced Tier-2 design (§4/Tier 2) comfortably a workstation job.

Caveat, from the same figure: the tolerance curves separate by roughly one
decade of `δ_FP` per decade of `rtol`, and `δ_FP` climbs with dose regardless.
Neither matters while the criterion is 1e-4, but if conservation is ever
reinstated as a scored quantity the LF and HF tiers would no longer be
comparable and the AR-GP `ρ` would absorb the difference.

Two corrections to earlier revisions of this plan:
- The claim that "a 100× tolerance relaxation cost ~4.7 orders of magnitude" was
  **wrong** — those runs differed in temperature and mobility cutoffs as well as
  tolerance. `rtol = 1e-4` is admissible.
- **`δ_FP` must be read as the converged (final) value, never the max.** The
  maximum always lands on the first output point and is a fixed startup artifact
  (`1.764e-01`) independent of every parameter.

**`δ_FP` grows roughly linearly with dose** (1.2×10⁻⁸ at 1 dpa → 3.6×10⁻⁶ at
10 dpa), so **acceptance criterion 3 (§6) should be dose-scaled**, not a flat
1e-6 — a 30 dpa HF run cannot be held to the same absolute figure as a 1 dpa LF
run. Note also that `C_floor = 1e-15` pinned the entire state at the floor on a
1 dpa discrete run — far too high to be a candidate; use 1e-25 and 1e-40.

Accept the loosest setting for which every observable in §3.1 changes by < 1 %
against the tightest. **Deliverable:** `report/T0_tolerances.csv` + a frozen
`solver_config` block reused verbatim at both fidelities.

**T0.2 Grid-convergence and the LF/HF definition.** Sweep
`(I, V) ∈ {(300,300), (500,500), (1000,1000), (2000,2000)}` and
`(i_discrete, I_bin) ∈ {(20,10), (50,20), (100,30)}` at the same three
temperatures, recording observables *and wall-clock*. LF is then the coarsest
grid whose observables at `d_LF` are within a factor 1.3 of the finest; HF is
`(1000,1000)/(50,20)` or better. **Deliverable:** `report/T0_grid.csv`, and the
two entries `fidelity.LF` / `fidelity.HF` in `conditions.yaml`.

**T0.3 Cost model.** Fit `log C_CPU = f(T, d_target, N_eq, i_mobile)` to the
T0.1/T0.2 timings plus the existing `output/*.log` history. Required by
Eq. (H.28); also the input to the honest answer on how big `N_LF` can be.
**Deliverable:** `report/T0_cost_model.json`, and a go/no-go statement of the
form *"at the verified tolerance, N_LF = X LF runs cost Y core-hours"*.

**T0.4 Reachability inside the data envelope.** The critical unknown, now posed
against a bounded box rather than an open-ended dose target. The **campaign
envelope** is fixed by the existing neutron data:

```
T ≤ 400 °C          d ≤ 30 dpa          wall-clock ≤ 7200 s per run
```

and every run terminates at whichever limit it reaches first. Nothing in the
neutron database lies outside it — the hottest neutron row is 400 °C
(Materna-Morris, 17.2 dpa) and the highest dose is the 335 °C / 32 dpa point,
2 dpa over the cap and handled as the single extrapolation case of §4/Tier 3.
The two 450–500 °C conditions that dominated the §3 conversion campaign are
therefore **out of scope for the twin**, which removes the most expensive and
most dose-starvation-prone part of the parameter space at no cost to the data
comparison.

At the LF grid, integrate to the largest dose reachable within 7200 s at 250,
300, 330, 350, 375, 400 °C. Existing evidence says the upper end of that window
is the worst case (0.05 dpa in 4892 s at 450 °C and 0.5 dpa in 6047 s at 400 °C,
both at the HF grid). **Deliverable:** `report/T0_reach.csv` — dose reached, and
dose reached per CPU-minute, versus temperature, at both LF and HF grids, with an
explicit `envelope_met: {dose|clock}` column recording which limit bound each
run.

If the mid-temperature window cannot reach `d_LF = 0.1 dpa` cheaply, the campaign
stops here and the fix is a solver problem, not a UQ problem: candidates are the
`active_window` solver mode, `he_kinetics='quasi_steady_state'`, a looser
`C_floor`, or Woodbury `prec_bw`/`prec_rank` tuning. Do not proceed by relaxing
tolerances — that is precisely what §2.4 of the paper forbids.

**T0.4b Network-channel liveness.** New gate, from §2.4. On the LF grid selected
in T0.2, run one condition with `LOOP_NETWORK_LOSS = 1` through
`run_adaptive()` and assert `dρ_net/dt ≠ 0`. If the LF loops stay too small for
`χ·d_loop ≳ ρ_net^{−1/2}`, the channel is silently off and parameters #20, #24
and `K_rec` will screen as inert for a purely numerical reason. **Deliverable:**
one line in `report/T0_grid.csv` per candidate LF grid: `rho_net_live: yes/no`.
An LF grid that fails this is not usable, regardless of how well its observables
converge.

**T0.5 Code fixes.** In dependency order:

| | Task | From | Status | Acceptance |
|---|---|---|---|---|
| (a) | `Z_i_loop` its own key in `K_loop` + `reactions` row + C++ mirror | §2.2(a), directive 2 | **open** | bit-identical run at `Z_i_loop = Z_i` |
| (b) | `lambda` / `A_void_0` threaded from the workbook into `E_b_void` | §2.2(b), directives 3–4 | **done** | `G_VAC` verified to move with each |
| (c) | `E_b_i2` key + `A_111`/`A_100` override | §2.2(c), directive 7 | **mechanism done, value disabled** | `E_b_loop_i(2) = 0.8001 eV`; blank cell ⇒ legacy. **0.80 eV fails the 10 dpa gate** — needs an `(A_111, B_111)` refit, §2.2(c) |
| (c2) | Re-fit `(A_111, B_111)` (and `A_100`/`B_100`) to DFT small-`n` binding, then re-validate at 10 dpa | validation | **open** | `δ_FP < 1e-6` at the production config with `E_b^i(2)` in 0.6–1.2 eV |
| (c3) | Fix `run_adaptive` truncation of loop-conversion output arrays | validation §5 | **open** | `f_111_loop`/`N_loops_100` length == `len(t)`; the 4 broken plots render |
| (d) | `loop_net_*` + `LOOP_NETWORK_LOSS` rows added to the workbook and `create_excel.py`; blank-cell readers hardened | §2.4, directive 8 | **done** | nominals = code defaults ⇒ no result change; `Λ_net` finite on the blank path |
| (e) | `extract_observables.py` lifted out of `visualization.py` and unit-tested against a saved run | §3.1 | open | reproduces the saved run's plotted panel |
| (f) | `param_io.py` round-trip test (`θ → workbook → InputData → θ`) | §3 | open | exact for continuous keys, exact for `cat` |
| (g) | Resolve the remaining §2.3 conflicts (`f_cl_v`, `f_cl_i` docs) | §2.3 | open | startup assertion: every nominal inside its own prior, `T_star_conv_C` whitelisted |
| (h) | Eq. (B.14) form for `E_b_void`; atomistic correction added to `E_b_bubble` | §2.2(b), open items | open | `E_b_bubble(m, 0) == E_b_void(m)` |

**No C++ mirror was needed for (b)–(d).** `cpp_bridge.py` exports *precomputed*
emission and rate arrays (`GII_k`, `Pr_SIA_k`, `k2_SIA_k`, `lambda_net_100_k`),
so any binding-energy or sink parameter that reaches the Python kernels reaches
the C++ solver automatically. Task (a) is the exception — `Z_i_loop` is exported
as a scalar and is aliased to `Z_i` at a **second** site,
[cpp_bridge.py:179](../../py_utils/cpp_bridge.py#L179), which must be fixed
together with `reaction_rates.py:132`.

**Verification standard used for (b)–(d), and required for (a).** A parameter is
accepted as wired only when perturbing it is shown to move a downstream array —
not merely when the key parses. Both silent-failure modes were found this way and
neither would have been caught by a schema check: `lambda`/`A_void_0` parsed
correctly for the whole life of the project while being ignored, and a blank
`loop_net_w_c` parsed to `NaN` rather than to its default.

### Tier 1 — Nominal anchor and prior-predictive check

One HF run per anchor condition at the **nominal** `θ`. Anchor set — the
EUROFER97/EUROFER-ODS rows of the database, which are the only ones where alloy,
loop character, cavity, and (for two) a size distribution coexist:

| # | ξ | Source | Targets available |
|---|---|---|---|
| N1 | EUROFER97, n, 250 °C, 16.3 dpa | ref 2 | `N_ℓ = 1.0e22`, `d̄_ℓ = 5.0 nm`, `f₁₀₀ = 10 %` |
| N2 | EUROFER97, n, 300 °C, 15 dpa | refs 3, 4 | `N_ℓ`, `d̄_ℓ`, `d̄_c = 2.5 nm`, `f₁₀₀ = 27–78 %`, **3 histograms** |
| N3 | EUROFER97, n, 330 °C, 15 dpa | ref 1 | `N_ℓ = 1.4e22`, `d̄_ℓ = 3.4`, `N_c = 1.4e21`, `d̄_c = 2.6`, `f₁₀₀ = 72 %` |
| N4 | EUROFER97, n, 335 °C, 32 dpa | ref 1 | `N_ℓ = 1.7e22`, `d̄_ℓ = 4.8`, `N_c = 1.1e21`, `d̄_c = 1.6`, `f₁₀₀ = 27 %` |
| N5 | EUROFER97, n, 350 °C, 16.3 dpa | ref 2 | `N_ℓ = 8.0e21`, `d̄_ℓ = 12.5`, `f₁₀₀ = 79 %` |
| N6 | EUROFER97, n, 400 °C, 17.2 dpa | ref 41 | `f₁₀₀ = 87 %` |
| N7 | EUROFER-ODS, n, 350 °C, 16.3 dpa | ref 5 | both characters resolved, **2 histograms** |
| I1 | EUROFER97, i, 350 °C, 3 dpa | ref 25 | `d̄_ℓ = 9.9 nm`, **1 histogram** |
| I2 | EUROFER97, i, 400 °C, 26 dpa | ref 17 | `d̄_c = 2.0 nm` (dual-beam) |
| I3 | EUROFER97, i, 330 °C, 16 dpa | ref 51 | `f₁₀₀ = 29 %` |
| I4 | EUROFER97, i, 400 °C, 16 dpa | ref 51 | `f₁₀₀ = 22 %` |
| I5 | EUROFER97, i, 400 °C, 26 dpa | ref 49 | `f₁₀₀ = 87 %` |

N3/N4 are the campaign's only genuine EUROFER97 **dose trajectory** at fixed
temperature (15 → 32 dpa, one laboratory, one material) and therefore the only
place where the sequential assimilation of §H.7 can be exercised on real data.
That is thin, and Tier 5 says so.

Ion conditions carry their **actual** dose rate — Boulanger/Serruys Kr²⁺ and the
dual-beam experiments run 2–4 orders above neutron — read from the workbook, not
inferred.

**Envelope check.** Every anchor condition above lies inside the
`T ≤ 400 °C, d ≤ 30 dpa` box of T0.4 except **N4 (335 °C, 32 dpa)**, which is
2 dpa over. N4 is retained — it is one half of the campaign's only dose
trajectory (below) — and is run to the 30 dpa cap, with the residual 2 dpa
carried by the discrepancy GP exactly as the general dose gap is in Tier 3. It is
the only anchor to which that applies, and the deliverable must say so per-row.

**Deliverables:** `report/T1_anchor.csv` (one row per condition: simulated
observables, experimental value, band, `|ln y_sim − ln y_fit| / s`);
`report/T1_prior_predictive.png` — Tier-1 predictions overlaid on the Appendix G
group-separated panels. **Acceptance:** not that the model fits — it very likely
will not — but that (i) every run passes the `< 1e-4` conservation *health
check* on the converged value (§6),
(ii) every observable is finite and physical, (iii) the residual table is
complete. This is the measurement of `δ^model` that sets `σ_model` in Eq. (H.13).

### Tier 2 — LF Sobol screening ensemble → `S_i`, `S_i^T`

Purpose: Eq. (H.24)–(H.25), demote the sloppy parameters before they can absorb
scatter.

- **Design:** Saltelli extension of a Sobol sequence, `N(p+2)` evaluations. With
  `p = 24` (§2.1) and a base **`N = 16`**, that is 416 model evaluations *per
  condition*.
- **Conditions:** three only — N2 (300 °C), N5 (350 °C), I1 (ion, 350 °C, high
  dose rate) — chosen to span the neutron temperature window and to include one
  ion point so the dose-rate-sensitive parameters (`L_1D`, `i_mob`) are exercised.
  All three are inside the T0.4 envelope. **≈ 1250 LF runs.**
- **`d_LF` = 0.1 dpa** (author decision, 2026-07-30; was 1 dpa).

**Why the cut — the cost is now measured, not assumed.** The anchor scan
([`anchor_scan_350C_1dpa.md`](anchor_scan_350C_1dpa.md)) timed 1 dpa at
`I=V=1000`, `im20/vm5`, `rtol=1e-6` at **250–450 s single-threaded**. At the old
design (`N = 64`, `d_LF = 1 dpa`) that is ~420 core-hours — the expensive branch
of the feasibility gate below, requiring a cluster. Dropping to `N = 16` and
`d_LF = 0.1 dpa` brings it to roughly **1250 runs × ~40 s ≈ 14 core-hours**,
trivially parallel on one workstation.

**What 0.1 dpa costs scientifically, stated plainly.** Sobol indices at 0.1 dpa
measure sensitivity in the *nucleation-dominated transient*, not at the
15–30 dpa of the database. Parameters that only matter once loops coarsen — the
loop→network channel most of all — will screen artificially low. Mitigations:
(i) carry the loop-evolution parameters (#19–#20, #24) into Tier 3 regardless of
their Tier-2 index, and (ii) note from the anchor scan that `f₁₀₀` is
*efficiency-limited, not dose-limited* (3 → 30 dpa bit-identical), so the
loop-character channel is largely unaffected by the truncation.
- **Prerequisite from §2.4:** the LF grid must have passed the T0.4b
  network-liveness check. Screening `w_c`, `χ` and `K_rec` on a grid where
  `P_ℓd ≡ 0` returns `S_i = 0` for all three and looks exactly like a genuine
  screening result.
- **Fidelity:** LF grid from T0.2, `d_LF` from T0.4 (nominally 0.1 dpa).
- **Feasibility gate:** at 10 s/run this is ~12 core-hours and trivially
  parallel; at 300 s/run it is 350 core-hours and needs either a cluster or a
  reduced design. T0.3 decides. If the budget forces a cut, cut in this order:
  (i) `N = 64 → 32`; (ii) three conditions → two; (iii) two-wave screening —
  a cheap 8-parameter wave on the production/mobility group first, then a second
  wave on the retained set plus energetics.
- **Failure handling:** Appendix F step 2 — *discard and flag* runs that fail the
  tolerance check; never relax the tolerance. Log the failure rate per parameter
  region; a region that fails systematically is itself a result (it means those
  `θ` are not integrable at the fixed tolerance and must be excluded from the
  prior support, with that exclusion recorded).

**Deliverables:** `report/T2_sobol_indices.csv` (`S_i`, `S_i^T` per parameter per
observable per condition); `report/T2_sobol_heatmap.png` (24 × 7 grid, the
figure that makes identifiability visible at a glance);
`report/T2_active_set.yaml` — the retained `θ_act` with `S_tol` recorded, and the
inactive parameters split into *fixed at nominal* versus *marginalised with prior*.

**Expectation to state up front:** Appendix G.6 already warns that only three of
eight panels resolve a temperature slope at 2σ and that the microstructure over
this database is controlled far more by dose, helium and alloy class than by
temperature. Screening will almost certainly retain **6–10** parameters, not 24.
That is the correct outcome, not a failure.

Two specific predictions worth recording now, so that confirming them counts as
evidence the screening is working rather than as a surprise:

- **`T*` (#18) should screen out.** §3 measured `T*` 250 → 180 °C as
  bit-identical in `f₁₀₀`, because the `ΔF > 0` cutoff sits far above the sizes
  the loops actually occupy. If `T*` comes back with a large `S_i^T` on the
  loop-fraction channel, something is wrong with the `ΔF` implementation, not
  with the screening.
- **`E_a0_conv` (#21) should dominate the `f₁₀₀` channel.** §3's measured
  0.325 → 0.948 response to a 0.2 eV change is the largest single-parameter
  effect anywhere in the campaign record.

### Tier 3 — HF design and the multi-fidelity anchor

Appendix F steps 3–4.

- `N_HF = 3p_act–8p_act`. With `p_act = 8`: **24–64 HF runs**, as a maximin
  nested subset of the Tier-2 design.
- Conditions: the full Tier-1 anchor set (12 conditions) for the 24 core points,
  the three screening conditions for the remainder.

**The HF envelope — revision 2.** Appendix F's `d_HF ≈ 100 dpa` is withdrawn.
The HF tier is bounded by the **existing neutron data**, and every HF run is
stopped by whichever of these binds first:

| Limit | Value | Rationale |
|---|---|---|
| Temperature | `T ≤ 400 °C` | hottest neutron row in the database (Materna-Morris, 17.2 dpa); above it there is nothing to compare against, and it is where the solver is most expensive |
| Dose | `d ≤ 30 dpa` | brackets the 15–32 dpa band the neutron rows occupy |
| Wall clock | **7200 s per run** | hard cap, enforced by the driver, not by the solver's own step control |

`d_HF` is therefore *per condition*: the dose actually reached when the run hit
one of the three limits. Three consequences, all of which the deliverable must
report explicitly:

1. **Every HF run records which limit stopped it** (`envelope_met` column). A run
   stopped by the clock at 0.5 dpa and a run stopped by the dose cap at 30 dpa
   are not the same object and must not be pooled without that label.
2. **Where `d_HF` falls short of the measured dose, the gap is carried by the
   discrepancy GP `δ(θ)`** — which is exactly what Appendix F step 11 is for. Do
   not compare a 1-dpa simulation to a 15-dpa measurement without that term. The
   §3 finding that `f₁₀₀` is *efficiency-limited, not dose-limited* (3 → 30 dpa
   bit-identical) is a useful partial reprieve: for the loop-fraction observable
   specifically, a dose-truncated run may still be informative. It is **not** a
   reprieve for `N_ℓ`, `d̄_ℓ`, `N_c`, `d̄_c` or `S`, all of which are strongly
   dose-dependent.
3. **A clock-limited run is not a failed run.** It is discarded only if it also
   fails the conservation gate (`δ_FP`, `δ_He < 1e-6`). Appendix F step 2's
   discard-and-flag applies to tolerance failures, not to the wall-clock cap.

The 7200 s cap also makes the Tier-3 cost bounded *a priori*: 64 HF runs × 12
conditions × 7200 s = 154 core-hours worst case, which is the first hard cost
number in the plan.

**Deliverables:** `report/T3_hf_runs.csv` (with `d_HF`, `envelope_met`, wall clock
and conservation diagnostics per run); `report/T3_reach.md` — the honest
statement of `d_HF` achieved per condition versus the measured dose, and the
fraction of runs that hit the clock rather than the dose cap. If that fraction is
high in the 350–400 °C band, it is the T0.4 solver problem resurfacing and Tier 4
should not start.

### Tier 4 — Emulator and calibration

Appendix F steps 5–10, then Eq. (H.14).

1. `log` transform + standardisation of every output channel (F-5).
2. LF GP with anisotropic Matérn-5/2, marginal-likelihood hyperparameters (F-6).
3. Discrepancy GP on `r = y_HF − ρ μ_LF`, `ρ` estimated jointly (F-7).
4. Fused predictor `M̂_HF = ρ M̂_LF + δ`, `σ²_HF = ρ²σ²_LF + σ²_δ` (F-8).
5. **Validation gate (F-9):** hold out 20 % of HF runs; require standardised
   errors of unit variance and nominal coverage of the predictive intervals. If
   coverage fails, increase `N_HF` and return to Tier 3. This gate is not
   optional — an over-confident emulator produces an over-confident posterior,
   and the whole point of the twin is the interval, not the mean.
6. Posterior by sequential Monte Carlo (preferred — it gives the evidence and
   handles multimodality) or ensemble Kalman inversion (cheaper, and the
   appendix's own suggestion for the sequential mode). Likelihood exactly
   Eq. (H.19): log-space `χ²` over the scalar set `S`, `w_W Σ W₁` over the
   distributional set `D`, plus `Φ_phys` from Eq. (H.20) using the run's own
   `δ_FP`, `δ_He`.
7. `σ_j²` assembled per Eq. (H.13) from all four terms: `σ_exp` (per-row, from the
   workbook where reported), `σ_db = s_j` (Table G.36), `σ_em` (emulator, from
   step 4), `σ_model` (from Tier 1).

**Deliverables:** `report/T4_posterior.nc` (samples);
`report/T4_corner.png` (pairwise marginals — the figure that shows sloppiness and
parameter correlation directly); `report/T4_posterior_predictive.png` (the
Appendix G panels with the posterior predictive band drawn over the `±2s` data
band — the single figure that answers "is the twin consistent?");
`report/T4_wasserstein.csv` (per-histogram `W₁` before and after calibration);
`report/T4_calibration_table.md` (prior vs posterior median and 90 % CI per
parameter).

### Tier 5 — Sequential assimilation

Eq. (H.26)–(H.27). Exercise the assimilation on the one available trajectory:

- Prior = Tier-4 posterior. Assimilate N3 (15 dpa) → propagate → predict N4
  (32 dpa) → assimilate → report the posterior shift.
- Report `p(θ, x_k | z_{1:k})` at each step and the **predictive** check: did the
  15-dpa-only posterior cover the 32-dpa measurement at its nominal level?

**This is a two-point trajectory from one laboratory.** It demonstrates the
machinery; it does not validate it. State that plainly in the deliverable and
list what would: a 5-point dose series in one material at one temperature, which
does not exist in the current database and is the top recommendation Tier 6's
experiment-selection should produce.

**Deliverable:** `report/T5_assimilation.md` + `report/T5_dose_trajectory.png`.

### Tier 6 — Active refinement and experiment selection

Eq. (H.28)–(H.29). `A(θ, ξ) = I(θ, ξ) / C_CPU(θ, ξ)` with `I` = expected
reduction in posterior predictive variance at DEMO conditions, `C_CPU` from the
T0.3 model. Two uses:

- **Next simulation:** append, refit from F-7, re-run Tier 4. Two or three
  iterations.
- **Next experiment:** maximise over `ξ` (T, dose, dose rate, He/dpa) with `C`
  replaced by an experimental-cost proxy. This is the output most likely to be
  useful outside the code, and the Appendix G confounds give it obvious
  candidates to test — a matched ODS/non-ODS pair from one campaign, and a
  helium-decoupled ion series that breaks the He-bubble/low-dose collinearity.

**Deliverables:** `report/T6_acquisition.png`; `report/T6_recommended_experiments.md`.

### Tier 7 — Fusion-condition prediction

Eq. (H.32) **as amended by directive 5**: `T = 300–550 °C`, `d = 0.1–150 dpa`,
`η_He = 10–15 appm/dpa`, `cascade='fusion'` (Case 1 mean-field He),
`physics_option='bin_moment_CD_fusion'`. **`η_H` is dropped from `ξ` — there is
no hydrogen in this model and none will be added.**

Posterior predictive by propagating ~200 posterior samples through the emulator,
with a small number of full HF confirmation runs at the corners. Report **credible
intervals, never curves** (acceptance criterion 5), for: `N_ℓ`, `d̄_ℓ`, `N_c`,
`d̄_c`, `S`, `f₁₀₀`, and `P(bimodal cavity regime)` — the last as the posterior
probability that the cavity distribution is bimodal, which is what
`reconstruct_distribution` plus a mode test gives directly.

**Hydrogen — scoping decision, stated in the deliverable.** `η_H = 40–50 appm/dpa`
appears in Eq. (H.32) and has no forward model in RadCluster_2_1. Revision 2
settles it: the Tier-7 prediction is **helium-only**, `predict_fusion.py` does not
accept an `eta_H` key, and every Tier-7 artefact carries a `he_only: true` field
in its provenance so a downstream reader cannot mistake the omission for a null
result. What is being given up should be named rather than hidden: H–He synergy
in bubble nucleation, and H trapping at loops and at the α′/MNS features of
objective (b), are outside the twin's predictive scope. If a hydrogen population
is ever added, this tier is the one that must be re-run — nothing upstream of it
depends on `η_H`.

**Extrapolation beyond the calibration envelope.** Tier 7 predicts to 550 °C and
150 dpa, while Tier 3 calibrates only inside `T ≤ 400 °C, d ≤ 30 dpa`. That is a
deliberate extrapolation of roughly 150 °C and a factor 5 in dose, and it is
carried entirely by the emulator's own uncertainty plus the discrepancy GP. The
Tier-7 bands must therefore be plotted with the calibration box drawn on the same
axes, and the deliverable must state that outside that box the intervals are
prior-and-discrepancy-dominated. This is the single largest interpretive caveat
in the campaign and belongs in the figure, not the caption.

**Deliverables:** `report/T7_fusion_prediction.csv`;
`report/T7_fusion_bands.png`; `report/T7_swelling_risk.png`.

---

## 5. Run-count and cost summary

| Tier | Runs | Fidelity | Purpose | Blocking on |
|---|---:|---|---|---|
| T0.1 tolerances | 18 | HF grid, 1 dpa | fix `(rtol, atol, C_floor)` | — |
| T0.2 grid | 36 | mixed | define LF/HF | T0.1 |
| T0.4 reachability | 12 | LF + HF | dose-per-CPU-minute vs T inside the envelope | T0.2 |
| T0.4b network liveness | 2 | LF | `ρ̇_net ≠ 0` on the LF grid | T0.2 |
| T1 anchor | 12 | HF | prior predictive, `σ_model` | T0.5 |
| T2 screening | ~1250 | LF, 0.1 dpa | `S_i`, `S_i^T` → `θ_act` | T0.3 gate, T0.4b |
| T3 HF design | 24–64 | HF | AR-GP anchor | T2 |
| T4 calibration | 0 (emulator) | — | posterior | T3 validation gate |
| T5 assimilation | 0 | — | sequential update | T4 |
| T6 refinement | 3 × 8 | HF | acquisition loop | T4 |
| T7 prediction | ~8 confirm | HF | fusion bands (He-only) | T4, T6 |
| **Total simulator runs** | **≈ 1430** | | | |

The count rose from ≈4380 because `p` went 20 → 24 (§2.1), which costs
`3 × 64 × 4 = 768` extra LF evaluations in the Saltelli design. The 5000 LF
screening runs are 96 % of the count and essentially all of the schedule risk.
**T0.3 must produce a per-run LF cost before Tier 2 is authorised**, and the
fallback ladder in §4/Tier 2 must be exercised rather than the tolerance relaxed.

**The 7200 s cap is what makes the HF side bounded.** Every HF tier now has a
computable worst case rather than an open-ended one:

| Tier | HF runs | Worst-case core-hours at the cap |
|---|---:|---:|
| T1 anchor | 12 | 24 |
| T3 design | 64 × 12 conditions | 1536 |
| T6 refinement | 24 | 48 |
| T7 confirmation | 8 | 16 |
| **HF total** | | **≈ 1620 core-hours** |

T3 dominates and is the place to cut if the budget binds: 24 core points instead
of 64 brings it to 576 core-hours. Note that the worst case assumes *every* run
hits the clock; T0.4 will say what fraction actually does, and at 300 °C and
below the historical evidence (117 s to 3 dpa at 550 °C, but the cold end is the
stiff one) suggests the true figure is closer to the cap than to the floor.

---

## 6. Acceptance gates (Appendix H.10)

Encoded as executable checks in `report/acceptance.py`, run after Tier 4 and
again after Tier 7:

| # | Criterion | Test | Data |
|---|---|---|---|
| 1 | Statistical consistency | `\|ln y_sim − ln y_fit\| ≤ 2s` for every row | Table G.36 bands |
| 2 | Distributional consistency | `W₁ ≤ W_tol` per histogram; `W_tol` set from the *between-analysis* spread of the three Dethloff 300 °C histograms — i.e. the tolerance is the experimental reproducibility, not a guess | `META`/`HIST` |
| 3 | ~~Conservation~~ **(removed from acceptance — see below)** | — | — |
| 4 | Predictive validation | held-out coverage ≈ nominal | Tier-4 holdout, Tier-5 32-dpa point |
| 5 | Extrapolation discipline | Tier-7 output contains no deterministic curve | artefact inspection |

### Criterion 3 withdrawn — conservation is a health check, not an acceptance gate

**Author decision, 2026-07-30.** `δ_FP` and `δ_He` are **removed from the
acceptance function and from the likelihood**. Any value **below 1e-4** indicates
correct integration and a correctly assembled reaction set; there is no
information about `θ` in *how far* below.

Rationale, from the measurements of
[`validation_10dpa_revision3.md`](validation_10dpa_revision3.md) §8:

- The converged `δ_FP` sits at 1e-8–4e-6 across every tolerance, temperature and
  dose tested — always passing, so as a gate it never discriminates.
- It varies with **dose** (~linearly) and with **integrator tolerance**, i.e. it
  measures numerics, not physics. Feeding it into a likelihood over `θ` would
  reward parameter regions that happen to integrate cleanly, which is a
  numerical artefact masquerading as evidence.
- Its maximum is a fixed startup artefact (`1.764e-01`), so a naive
  implementation of the old criterion would have failed every run.

**Replacement:** a per-run **health check**, logged but not scored —
`δ_FP < 1e-4` and `δ_He < 1e-4` on the **converged (final)** value, never the
max. A run exceeding 1e-4 is flagged as a *numerical* failure and discarded per
Appendix F step 2; it is not evidence against its `θ`.

**Consequence for Eq. (H.20):** the `Φ_phys` conservation penalty term is
dropped from the negative log-likelihood of Eq. (H.19). The remaining physical
constraints on `θ` come from the log-space `χ²` over `S` and the `W₁` terms over
`D`.

Criterion 2's tolerance definition is worth stating explicitly because it is the
one place the database gives a free calibration of "how close is close enough":
Dethloff 2018 measured the same specimen three ways and got mean loop diameters
of 4.68, 3.43 and 4.76 nm. A twin that lands anywhere in that spread has matched
the experiment as well as the experiment matches itself.

---

## 7. Traceability — where each Appendix H concept surfaces

| Appendix H | Concept | Artefact |
|---|---|---|
| H.1, Eq. (H.1) | State `x(d,T;θ)` | `results_y.npy` + `rho_net` + `J_*` in every run dir |
| H.1, Eq. (H.5) | Observable vector | `extract_observables.py`; columns of `T1_anchor.csv` |
| H.1, Eq. (H.6) | Posterior `p(θ\|D)` | `T4_posterior.nc`, `T4_corner.png` |
| H.1, Eq. (H.7) | Posterior predictive | `T4_posterior_predictive.png`, `T7_fusion_bands.png` |
| H.2, Eq. (H.8) | 24-parameter vector + priors | `parameters.yaml`, `T4_calibration_table.md` |
| §3 (ref. doc) | Loop-conversion / network ranges | `parameters.yaml` rows #18–#24; §2.4 |
| H.3, Eq. (H.12) | Log-space residual | `T1_anchor.csv` residual column |
| H.3, Eq. (H.11) | Band acceptance `e^{±2s}` | acceptance criterion 1 |
| H.3, Eq. (H.14) | Wasserstein `W₁` | `T4_wasserstein.csv` |
| H.4, Eq. (H.13) | Four-term `σ_j²` | `targets.yaml` + emulator + `T1` `σ_model` |
| H.4, Eq. (H.19) | Negative log-likelihood | `calibrate.py` |
| H.4, Eq. (H.20) | ~~`Φ_phys` penalty~~ **dropped from the likelihood** | `δ_FP`, `δ_He` logged as a health check only, §6 |
| H.5, Eq. (H.23) | AR multi-fidelity GP | `surrogate.py`; `T3_hf_runs.csv` |
| H.6, Eq. (H.24–25) | Sobol `S_i`, `S_i^T` | `T2_sobol_indices.csv`, `T2_sobol_heatmap.png` |
| H.6 | Active set `θ_act` | `T2_active_set.yaml` |
| H.7, Eq. (H.26–27) | Sequential update | `T5_assimilation.md`, `T5_dose_trajectory.png` |
| H.8, Eq. (H.28–29) | Acquisition `A = I/C` | `T6_acquisition.png`, `T0_cost_model.json` |
| H.8 | Experiment selection | `T6_recommended_experiments.md` |
| H.9, Fig. H.25 | Offline workflow | tiers T0–T4 |
| H.9, Fig. H.26 | Online assimilation | tier T5 |
| H.10 | Acceptance criteria | `report/acceptance.py` |
| H.11 | Implementation checklist | `digital_twin/` §3 |
| H.12, Table H.37 | Target measurements | `targets.yaml` |
| H.12, Table H.38 | Prior ranges | `parameters.yaml` (with the §2.2/§2.3 errata and the revision-2 overrides) |
| §2.4, Fig. 2 | Four-stage pipeline | T2 → T3/T4 → T4 → T6 |

---

## 8. Risks, ranked

1. **Dose reach (high, but no longer campaign-ending).** 100 dpa at HF is not
   attainable in the 250–400 °C window. Revision 2 removes the risk of an
   unbounded chase by *defining* the HF tier as the data envelope
   (`T ≤ 400 °C`, `d ≤ 30 dpa`, 7200 s) — the campaign can now always terminate,
   and the open question becomes how much dose it buys rather than whether it
   finishes. Residual risk: if the 350–400 °C band reaches only ~0.1 dpa in
   7200 s, the discrepancy GP is being asked to bridge two orders of magnitude in
   dose, which it cannot do credibly. Mitigation: T0.4 measures this before
   anything is committed; solver-side fixes (`active_window`, QSS helium,
   preconditioner tuning) are tried before any tolerance is touched; if the gap
   is still two orders, report the twin as calibrated on the low-dose transient
   only and say so.
2. **Screening cost (high, schedule).** ~5000 LF runs is the whole budget, and
   `p = 24` made it 18 % worse than revision 1. Mitigation: the three-step
   fallback ladder, decided by the T0.3 cost model. Note that the ladder's first
   rung (`N = 64 → 32`) recovers more than the parameter-count increase cost.
2b. **Silently inert parameters (was high — now largely closed).** A parameter
   the ensemble varies and the solver ignores produces `S_i = 0`, which is
   indistinguishable from a genuine screening result and will be believed. Two
   instances were found and fixed in revision 3: `lambda`/`A_void_0` parsed
   correctly but were hard-coded in `binding_energies.py` and ignored; a blank
   `loop_net_w_c` parsed to `NaN` rather than to its default. `E_b_i2` and the
   `loop_net_*` family now exist as real rows. **One instance remains open:**
   `Z_i_loop` is still aliased to `Z_i` at *two* sites
   (`reaction_rates.py:132`, `cpp_bridge.py:179`). Mitigation: the perturbation
   test of T0.5 — vary each parameter and confirm a downstream array moves —
   applied to `Z_i_loop` when task (a) lands, and to every `active: true`
   parameter before Tier 2.
3. **Weak data (high, scientific).** 68 rows, 8 fitted panels, 3 resolved slopes,
   6 histograms, one 2-point dose trajectory. The posterior will be
   prior-dominated for most of `θ`. Mitigation: report it — a wide posterior
   honestly obtained is the correct product, and Tier 6's experiment
   recommendations are the constructive response.
4. **Confounded groups (medium).** Appendix G identifies two that cannot be
   resolved on this data: ion cavity density (He-bubble ≡ ≤30 dpa, exactly the
   same six rows) and neutron cavity density/size (whole trend on two pure-Fe
   rows, dropping which reverses both slopes). Mitigation: down-weight those
   rows, and never let the twin claim to have "learned" the helium effect from
   the ion cavity panel.
5. **Observation-operator mismatch (medium).** Content- vs number-weighted `f₁₀₀`,
   black-dot conventions, TEM cutoff. Mitigation: emit both conventions, sample
   `d_min`, fold the residual into `σ_model`.
6. **Parameter/workbook conflicts (medium, silent).** §2.3 — now four, including
   the new `A_111` / `E_b^i(2)` disagreement, and one where the workbook and
   `create_excel.py` disagree with *each other*. Mitigation: resolve before
   Tier 1; add a startup assertion that every nominal lies inside its own prior.
7. **No hydrogen (low, scoping — closed).** Eq. (H.32) asks for `η_H`; the model
   has none, and revision 2 decides none will be added. Mitigation: Tier-7 output
   is labelled helium-only in provenance, and `predict_fusion.py` rejects an
   `eta_H` key outright so the gap cannot be papered over later. Residual
   exposure: H–He synergy in bubble nucleation is simply absent from the
   prediction.
8. **Fusion extrapolation beyond the calibration box (medium, new).** Tier 7
   predicts to 550 °C / 150 dpa from a posterior calibrated at `≤ 400 °C` /
   `≤ 30 dpa`. The bands there are prior- and discrepancy-dominated. Mitigation:
   draw the calibration box on every Tier-7 figure; state it in the deliverable
   rather than the caption.

---

## 9. Recommended immediate next steps

Updated for revision 3. Steps 2–4 of revision 2 are **done**; what remains:

1. **Give `Z_i_loop` its own entry and key** — `reaction_rates.py:132` **and**
   `cpp_bridge.py:179` (two alias sites), plus the `reactions` workbook row.
   §2.2(a). Acceptance: bit-identical at `Z_i_loop = Z_i`, then a perturbation
   test showing `K_SIA_grow`/`K_SIA_loop` move when it is varied alone.
2. **Re-run the affected calibrations.** `E_b_i2 = 0.80 eV` changes the SIA
   emission rates by ~27 orders of magnitude at `n = 2` (§2.2(c)), so any result
   tuned against `A_111 = 3.0` — in particular the §3 loop-conversion campaign
   and the `E_a0_conv ≈ 1.6 eV` recommendation — is no longer anchored. Re-run
   at least one §3 condition (500 °C, 3 dpa, `E_a0 = 1.6`) and confirm the
   `f₁₀₀` response is unchanged in character before the §2.4 priors are trusted.
   **This is the highest-value check in the list.**
3. **Resolve the remaining §2.3 conflicts** — `f_cl_v` against Appendix A, and
   the `f_cl_i` = 0.58 / 0.25 documentation disagreement. `T*` is now deferred to
   Tier 2 by design; `E_b^i(2)` is closed.
4. **Close the two `E_b_void` items of §2.2(b)** — adopt the Eq. (B.14) form, and
   give `E_b_bubble` the atomistic correction so `E_b_bubble(m, 0) = E_b_void(m)`.
   Needed before Tier 3, since the fusion path passes continuously through
   low-`ℓ` bubbles.
5. **Lift the observation operator** out of `visualization.py` into
   `digital_twin/extract_observables.py`, with the saved-run regression test.
6. **Run Tier 0.4 and 0.4b** — the reachability table inside the
   `T ≤ 400 °C, d ≤ 30 dpa, 7200 s` envelope, plus the network-liveness check.
   Fourteen runs, and together they decide whether the campaign is executable and
   whether the loop→network parameters can be screened at all.
