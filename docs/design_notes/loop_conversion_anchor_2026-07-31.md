# Loop-conversion anchor runs — re-anchoring the conversion parameter ranges

**Date:** 2026-07-31
**Purpose:** finalize the loop-conversion parameter ranges for the digital-twin
plan, using the current experimental database and the twin's own temperature
envelope.
**Supersedes (in part):** the calibration premise of
[`loop_conversion_calibration.md`](loop_conversion_calibration.md) — see §1.
**Companion:** [`tier2_tier3_parameter_selection.md`](tier2_tier3_parameter_selection.md),
[`loop_growth_grid_limit.md`](loop_growth_grid_limit.md).

---

## 0. Why the existing calibration had to be redone

Three things invalidated the §3 record as a source of twin priors:

1. **It was run outside the twin's envelope.** The §3 campaign
   (`digital_twin_implementation_plan.md` §2.4) was measured at 450–500 °C. The
   twin's HF envelope is now `T ≤ 400 °C`, fixed by the neutron database. The
   two conditions that dominated the campaign are out of scope.
2. **It was calibrated against a different reading of the database.**
   `loop_conversion_calibration.md` set `E_a0_conv = 1.8 eV` to reproduce
   "f₁₁₁ ≈ 1.0 at 250–300 °C → ≈ 0 by 335–400 °C". The database as it now
   stands does not say that (§1).
3. **Its grid was too small for the interesting region.** At `I = 200` with
   `i_mobile = 10`, the runs where conversion is weak — the ones that set the
   upper end of the `E_a0_conv` range — reach `n̄₁₁₁ ≈ 50–81`, i.e. occupancy
   0.27–0.41 against the 0.1 admissibility bound, with `δ_FP` up to 1.8e-01.

---

## 1. The experimental target, restated

`Loop_fractions` sheet of `FerriticSteels_RadiationDatabase.xlsx`, EUROFER97
only (ODS variants excluded), ⟨100⟩ **number** fraction:

| T [°C] | dose [dpa] | f₁₀₀ | irradiation |
|---:|---:|---:|---|
| 250 | 13.4 | 0.10 | neutron |
| 300 | 14.6 | 0.27 | neutron |
| 300 | 15.0 | 0.78, 0.77 | neutron (Dethloff 2016 / 2018) |
| 330 | 15.0 | 0.72, 0.45, 0.73 | neutron (Chauhan 2021) |
| 330 | 32.0 | 0.27, 0.31 | neutron |
| 350 | 17.4 | 0.79 | neutron |
| 400 | 17.2 | 0.87 | neutron |
| 415 | 18.1 | 0.80 | neutron |
| 330 | 16.0 | 0.29 | ion (Kaiser 2018) |
| 400 | 16.0 | 0.22 | ion |
| 400 | 26.0 | 0.87 | ion (Brimbal 2015) |

In-envelope (`T ≤ 400 °C`) neutron: **n = 11, range 0.10–0.87, median 0.72**.

**⟨100⟩ is the majority loop character at 300–400 °C in most measurements** —
the opposite of the premise `E_a0_conv = 1.8` was fitted to.

Three further EUROFER97 rows (300 °C/15.0, 300 °C/16.3, 350 °C/16.3) carry
`f₁₀₀ = 0`, but each is a `Loop Type = 1/2<111>` row with **no paired ⟨100⟩
row** — "⟨100⟩ not reported" is indistinguishable from "⟨100⟩ not present" in
this encoding. They are excluded from the band above and must be flagged in
`targets.yaml` rather than silently read as zeros. Note that the 300 °C/15.0
condition has *both* a 0.0 row and Dethloff's 0.78, so at least one of the two
readings is a convention artefact.

**Scatter is the dominant feature, not the temperature trend.** At 300 °C the
reported values span 0.27–0.78; at 330 °C, 0.45–0.73 at 15 dpa. Any acceptance
function that treats a single f₁₀₀ number as a tight target will be fitting
inter-laboratory convention differences. The 330 °C pair also shows f₁₀₀
*falling* with dose (0.72 → ~0.29 from 15 → 32 dpa), which is worth noting
because the twin plan (§2.4) treats a dose-decaying f₁₀₀ as a *spurious*
signature of a unary-only model.

---

## 2. THE BLOCKING FINDING — the model's ⟨100⟩ loops are invisible

This outranks every parameter range below and should be read first.

| | model (all runs, `E_a0` 1.5–1.8, 300–400 °C) | experiment (EUROFER97 + variants) |
|---|---|---|
| mean ⟨100⟩ diameter | **0.74 – 0.81 nm** | **2.8 – 180 nm, median 7.9** (n = 27) |

The model's ⟨100⟩ population sits at `n₁₀₀ ≈ 11` atoms and **does not grow**.
Applying the same TEM detection cutoff to both loop characters:

| `E_a0` | `T` | `f₁₀₀` number | `f₁₀₀` content | `f₁₀₀` TEM > 1.0 nm | `f₁₀₀` TEM > 1.5 nm |
|---:|---:|---:|---:|---:|---:|
| 1.5 | 400 °C | 0.973 | 0.771 | 0.803 | **0.019** |
| 1.6 | 400 °C | 0.827 | 0.180 | 0.301 | **0.002** |
| 1.6 | 350 °C | 0.579 | 0.099 | 0.127 | **0.000** |
| 1.7 | 300 °C | 0.060 | 0.020 | 0.008 | **0.000** |

**At a realistic 1.5 nm TEM cutoff the model predicts `f₁₀₀ ≈ 0` at every
parameter value tested, against a measured 0.10–0.87.**

Consequence for the twin: **`f₁₀₀` is not currently a usable calibration
target.** Calibrating `E_a0_conv` against the raw number fraction would produce
a tight, confident posterior on a parameter that cannot reproduce the observable
once the measurement convention is applied — the failure is in the ⟨100⟩ *size*
distribution, not in the *fraction*, and no conversion parameter addresses it.

### 2.1 RESOLVED — the cause is a missing mechanism, not a parameter

A further ~25 runs (2026-07-31, later the same day) eliminated every candidate
parameter and identified the structural cause. **The ⟨100⟩ *content* is right;
the *number* is 125× too high; the size follows arithmetically.**

| quantity | model / experiment |
|---|---|
| ⟨100⟩ content | ≈ correct |
| `N₁₀₀` | **125×** too high |
| `d₁₀₀` | 0.10× — and `√125 = 11.2` vs the measured ratio **10.0** |

Because `d ∝ √n` and `n ∝ content/N`, a number density 125× too high forces a
diameter ~11× too small *given the correct content*. That single fact explains
every negative result below.

**Four mechanisms tested, all eliminated:**

| mechanism | knob, and how far it was pushed | `d₁₀₀` |
|---|---|---|
| junction nucleation | `i_mobile` 10→50 (threshold live, `n_j_min_eff` = 30) | nothing populates `n ≥ 60`; `φ_max` 0.5 vs 1.0 bit-identical |
| unary rate | `E_a0` 1.6→2.0 — unary rate down **1700×** | 0.79 → 0.83 nm |
| ΔF gate | `T*` 450→150 — gate `n ≤ 35` → `n ≤ 117` | 0.79 → 0.80 nm |
| **absorption** | `absorb_boost_100` 1→6568 (= `rot_factor`, i.e. pure 1-D glide with the isotropic penalty fully removed) | **0.79 → 0.80 nm** |
| monomer bias | `Z_i_loop` 1.05→1.30, full dose | 0.79 → 0.79 nm |

Growth boosts cannot work: the absorbable ½⟨111⟩ inventory is fixed, and spread
over 4.7×10²³ loops m⁻³ it gives each loop a few atoms however fast the transfer.
Suppression cannot work either: raising `E_a0` cuts number **and** content
together (`f₁₀₀` 0.637 → 0.007), leaving the size unchanged. Fewer-but-larger
requires **redistribution**, and the model has no redistribution channel.

**The structural gap.** The `sia_100` population has *two sources* (junction,
unary) and no removal mechanism at all:

- **no ⟨100⟩×⟨100⟩ coalescence** — deliberate (§4.2 of
  [`loop_111_to_100_conversion.md`](loop_111_to_100_conversion.md): "`bulk-100`
  gets no SOURCE and no self-COALESCENCE"), on the grounds that two sessile
  loops cannot collide diffusively;
- **no dissolution** — measured emission/growth ratio at 350 °C is
  **1e-29 (n=4)** to **1e-54 (n=20)**. With `E_b^100(n) = A_100·n^{B_100}`,
  `A_100 = 3.0`, a 4-atom ⟨100⟩ loop binds at ≈4.9 eV and is effectively
  immortal;
- the only optional sink is the loop→network channel, off by default and
  already measured weak ([`loop_growth_grid_limit.md`](loop_growth_grid_limit.md) §5).

So `N₁₀₀` is a monotonically accumulating counter of conversion events. Nothing
in the current mechanism set can reduce it at fixed content.

**Consequence for the twin: `f₁₀₀`, `N₁₀₀` and `d₁₀₀` are not calibratable**
with the present physics. They should be excluded from the acceptance function —
not down-weighted — until a ⟨100⟩ number-reducing mechanism exists. Calibrating
against them now would drive `E_a0_conv` to whatever value minimises a residual
the model cannot physically attain.

**Candidate mechanisms to add** (in rough order of physical defensibility):

1. **⟨100⟩ loop impingement/coalescence.** Sessile loops cannot collide
   *diffusively*, but they can grow into one another. At `N₁₀₀ ~ 10²³ m⁻³` the
   mean spacing is ~2 nm, comparable to the loop size — impingement is not a
   small correction at these densities. This is the most likely missing term.
2. **⟨100⟩ absorption by the network** — `lambda_net_100_k` already exists in
   `cpp_bridge`; enabling `LOOP_NETWORK_LOSS` gives ⟨100⟩ a removal path.
   Previously measured weak, but it was never tested as a *number* control.
3. **Revisit `A_100`.** `A_100 = 3.0` (matching `A_111`) makes small ⟨100⟩ loops
   unconditionally stable. Table 18's value is **0.7160**, which would make
   `E_b^100(4) ≈ 1.2 eV` instead of 4.9 eV and open an Ostwald-ripening path.
   The workbook value may simply be wrong.

**Caveat on the `i_mobile` runs.** The `i_mobile ≥ 30` cases dose-starved at
0.007–0.035 dpa (against 0.3) *and* were grid-limited (`n̄₁₁₁ = 350` at
`I = 800`, occupancy 0.44). Their qualitative reading — nothing populates ⟨100⟩
above `n = 60` — stands, but no number from them is quotable. Separately this
establishes that **`i_mobile = 50` with conversion is computationally out of
reach** in discrete mode at an admissible grid: 5400 s bought 0.007 dpa. That is
a harder constraint on the twin than the 7d defect and is not fixed by fixing it.

---

## 3. `E_a0_conv` — the only effective lever

Dose-matched at 0.1 dpa (the LF-tier dose), discrete, `i_mobile = 10`,
number-weighted `f₁₀₀`, `I = 600`:

| `E_a0` | 300 °C | 350 °C | 400 °C |
|---:|---:|---:|---:|
| 1.5 | 0.533 | **0.797** | 0.973 |
| 1.6 | 0.261 | 0.579 | **0.827** |
| 1.7 | 0.060 | 0.351 | 0.676 |
| 1.8 | 0.007 | 0.122 | 0.486 |
| **expt** | **0.27–0.78** | **0.79** (330 °C: 0.45–0.73) | **0.87** |

At 3 dpa (`I = 200`) the same sweep runs 1.2 → 2.0 and shows the saturation at
both ends: `f₁₀₀` ≈ 1.00 for `E_a0 ≤ 1.4` at every temperature, and ≈ 0.02 at
2.0. **The informative support is ≈1.45–1.85**; the plan's 1.4–2.0 spent both
ends in flat regions.

**Best fit is dose-dependent**, which matters because LF and HF run at different
doses: at 0.1 dpa the data prefer `E_a0 ≈ 1.5`; at 3 dpa, ≈1.6. The 350/400 °C
columns are dose-converged (identical at 0.1 and 3 dpa), but **300 °C is not**
(0.260 at 0.1 dpa vs 0.356 at 3 dpa) — the slowest condition has not
equilibrated by the LF dose. The multi-fidelity model must absorb that as a
genuine LF/HF discrepancy at low temperature, not as noise.

Adopted: **1.45–1.85, nominal 1.6.**

## 4. Grid robustness of `f₁₀₀`

Eight `I = 200` vs `I = 600` pairs at matched dose (four discarded — the `I=600`
partner starved):

| verdict | count | condition |
|---|---|---|
| stable (`Δf₁₀₀ ≤ 0.023`) | 5 | wherever `n̄₁₁₁ ≲ 60` |
| **grid-sensitive** (`Δ` 0.066–0.156) | 3 | 350–400 °C at low `E_a0`, where `n̄₁₁₁` grows 85 → 281 |

So `f₁₀₀` inherits the occupancy rule: it is trustworthy exactly where the loop
distribution fits the grid. Worst case `Δf₁₀₀ = 0.156` — still below the
experimental scatter at one temperature (0.27–0.78), so grid resolution is not
the limiting uncertainty on this observable, but the grid-sensitive corner must
be flagged rather than averaged in.

## 5. The other three knobs

| Parameter | Measurement | Verdict |
|---|---|---|
| `φ_max` (#23) | 0.1 vs 1.0 gave **bit-identical** `f₁₀₀` | **Code defect, not physics.** `n_j_min_junc = 30` exceeds `i_mobile`, so `Θ(min(n,n′) ≥ n_j_min)` is never satisfied and the junction channel is identically zero. Since the twin samples `i_mobile` over 5–50, `φ_max` would have screened inert over most of its prior. **Fixed:** `n_j_min_eff = min(n_j_min_junc, ⌈0.6·i_mobile⌉)` — exactly 30 at `i_mobile = 50`, so production is unchanged; post-fix `φ_junc` scales 10× with `φ_max` |
| `ΔH₂` (#22) | 1.00 → 0.4813, 0.70 → 0.4819, 0.55 → 0.4912, 0.40 → 0.6585 | Flat above ≈0.45. Narrow to **0.35–0.50**, nominal 0.40 |
| `T*` (#18) | 250 → 450 °C moves `f₁₀₀` by 2.7 % (0.846 → 0.824) | **Drop from `θ`.** Not bit-identical as §3 reported, but far too weak to spend a dimension on |

**Dose independence.** `f₁₀₀` at 350 °C is identical to four decimals at 0.3 and
3 dpa, confirming §3's "efficiency-limited, not dose-limited" in-envelope. The
LF tier at 0.1 dpa can therefore see the conversion signal — unlike the
loop→network block, which is dose-driven and must wait for Tier 3.

## 6. Two defects found, one fixed

**(a) `run_adaptive` truncated seven output series — FIXED.**
`simulation._TS_KEYS` omitted `N_loops_111/100`, `mean_n_111/100`, `f_111_loop`
**and** `delta_FP_sia`/`delta_FP_vac`, so `_merge_results` kept only the last
segment. The conservation *arms* were affected too, which means earlier
diagnostics that read them under `run_adaptive` were reading a fragment. Also
added 2-D handling so `y_sia100` survives a domain doubling. This closes plan
item T0.5(c3).

**(b) `bin_moment` + `loop_conversion` is broken — OPEN, blocking.**
At the runaway base it returned `mean_n₁₀₀ = 2853` on a grid of `I = 1000`, and
`N₁₀₀ = 5.0×10²⁴ m⁻³` — 17 % of every atom in the lattice sitting in ⟨100⟩ loops
at 1 dpa. Both `I = 1000` and `I = 2000` produced the same pathological state,
which is why a naive grid check reports it as "grid-independent (×1.002)".

Isolated at `I = 300`, 350 °C, `E_a0 = 1.6`:

| run | `mean_n₁₀₀` | invariant `n₁₀₀ ≤ I` | `δ_FP` |
|---|---:|---|---:|
| strong conversion, discrete | 299.8 | ok (pinned *at* the ceiling) | 1.2e-01 |
| strong conversion, bin_moment | **433.0** | **VIOLATED** | 1.6e-01 |
| weak conversion, bin_moment | 11.8 | ok | 8.0e-06 |

So the defect manifests only under strong conversion; at normal production
values `bin_moment` is numerically clean. Root cause is on record: build step
**7d** (bin-moment ⟨100⟩ reconstruct→transfer→project) of
[`loop_111_to_100_conversion.md`](loop_111_to_100_conversion.md) §7 was never
implemented — "optional; discrete conversion is complete & validated" — and it
fails **silently** instead of refusing to run. An earlier symptom was already
visible: discrete and bin_moment disagreed by 27 % on `f₁₀₀` (0.113 vs 0.144) at
450 °C.

**Until 7d is fixed, `loop_conversion=1` requires `equations='discrete'`.** This
collides with T0.2: the grid-converged production configuration is `bin_moment`.

## 7. Does conversion bound loop growth?

No — on the evidence available, it relabels the runaway rather than stopping it.

- Raising `I` 200 → 600 did **not** reduce occupancy (0.25–0.40 → 0.06–0.48);
  the ½⟨111⟩ loops grew to the new ceiling.
- Under strong conversion the ⟨100⟩ population reaches the ceiling instead
  (`mean_n₁₀₀ = 299.8` at `I = 300`), growing by Marian absorption of mobile
  ½⟨111⟩, which has no size limit of its own.

This makes conversion the **second** candidate mechanism to fail this test,
after the loop→network channel
([`loop_growth_grid_limit.md`](loop_growth_grid_limit.md) §5). The occupancy
admissibility rule remains the operative control.

**Not fully established.** The clean version of this test — discrete mode,
grid-doubling, at `i_mobile = 50` — has not been run, for the same reason as §2:
it needs either T0.5(c5) fixed or a discrete run at a grid large enough to be
affordable. Stated as the current best reading, not a settled result.

## 8. What this changes in the twin plan

Applied to `digital_twin_implementation_plan.md` revision 4 (§2.4bis):

1. `E_a0_conv` 1.4–2.0 → **1.45–1.85**, nominal 1.6
2. `ΔH₂` 0.35–0.70 → **0.35–0.50**, nominal 0.40
3. `T*` **dropped** from `θ`
4. `φ_max` retained, now identifiable (fix c4); flagged as a screening canary
5. `K_rec` 10⁻³–10¹ → **10⁻⁶–10⁻²**; `K_rec`/`w_c`/`χ` deferred to **Tier 3**
6. Observation operator must emit **`f_100_tem(d_min)`** alongside number and
   content forms, and the calibration must consume the TEM-filtered one
7. `targets.yaml` must flag the three `f₁₀₀ = 0` rows as "not reported"
8. New blocking code item T0.5(c5); T0.5(c3) and (c4) closed

**Open, in priority order:**

1. Re-run the `f₁₀₀` anchor at `i_mobile = 50` to test whether the ⟨100⟩ size
   deficit (§2) is physics or a configuration artefact. **This gates whether
   `f₁₀₀` can be a calibration target at all.**
2. Fix 7d (T0.5(c5)) — required for (1) at a useful grid.
3. Re-run the discrete grid-doubling test of §7 at `i_mobile = 50`.
4. `WEAK_discrete` vs `WEAK_binmoment` agreement check was dropped when the job
   was killed for CPU oversubscription; redo to bound the bin_moment error at
   production values.
