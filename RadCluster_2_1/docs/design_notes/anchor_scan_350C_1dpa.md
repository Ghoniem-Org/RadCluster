# Anchor Parameter Scan — 350 °C, 1 dpa

**Date:** 2026-07-30
**Question:** which parameters move the simulation toward the EUROFER97 data, and
which cause stiffness?
**Condition:** EUROFER97 neutron anchor, 350 °C, `G = 1e-6` dpa/s, fission.
**Budget:** 1 dpa (`t_end = 1e6` s) or 1000 s wall clock, whichever first.
**Config:** `I=V=1000`, `i_mobile=20`, `v_mobile=5`, `i_discrete=50`,
`I_bin=V_bin=20`, linear shape, `full_system`+`bin_moment`, GMRES/Woodbury,
`rtol=1e-6`, `atol=1e-20`, `C_floor=1e-25`, He QSS, conversion OFF, network OFF.

15 runs total (10 + 5 follow-ups). Every run reached 1.0 dpa; none hit the cap.
Logs: `output/log_files/20260730_*anchor*`. Raw: `anchor_scan{,2}.json`.

---

## 0. Headline

1. **No admissible parameter closes the loop-density gap.** Every run except one
   sits at **N_ℓ ≈ 50–190× above** the fitted band and **d̄_ℓ ≈ 0.3–0.5× below**
   it. The best single-parameter move is ~1.5×; the gap is ~120×. **This is
   structural, not parametric** — see §4.
2. **`E_b(2)` is a threshold, not a dial.** 3.92 → 2.0 → 1.2 eV changes the
   answer by <0.1 %. The entire response is a bifurcation between 0.8 and 0.6 eV.
3. **The Production sheet was dead.** `eta`, `f_cl_i`, `f_cl_v`, `s_i`, `s_v`,
   `i_cascade`, `v_cascade` were read from hard-coded module dicts and the
   workbook was ignored. Found by this scan (two runs came back byte-identical
   to the reference, solver counters included). **Now fixed** — §5.
4. **Stiffness tracks the monomer fraction**, not the binding law: the stiffest
   runs are the ones that put more of the cascade into monomers (low `f_cl_i`)
   or more into vacancy clusters (high `f_cl_v`).

---

## 1. The harmonised binding law

`E_b(2)` is selected and `(A_111, B_111)` are pinned at **both** ends of the DFT
branch `E_b^fit(n) = A·n^{+B}` — at `n = 2` (di-interstitial) and at `n = 25`
(blend centre, → the continuum asymptote `E_f^i = 3.64 eV`):

```
B = ln(E_f^i / E_b2) / ln(25/2)        A = E_b2 · 2^(−B)
```

This is the correction to the revision-3 failure, where `E_b(2)` was changed by
sliding `A_111` alone — which rescaled the whole branch by a constant factor and
suppressed nucleation ([`validation_10dpa_revision3.md`](validation_10dpa_revision3.md)).

| `E_b(2)` | `A_111` | `B_111` | E_b at n = 2 / 5 / 10 / 20 / 25 / 40 |
|---:|---:|---:|---|
| 3.92 (legacy) | 3.0000 | 0.3873 | 3.92 / 5.59 / 7.31 / **8.75** / 6.60 / 2.91 |
| 2.00 | 1.6969 | 0.2371 | 2.00 / 2.49 / 2.93 / 3.36 / 3.20 / 2.89 |
| 1.20 | 0.8850 | 0.4393 | 1.20 / 1.79 / 2.43 / 3.23 / 3.20 / 2.89 |
| 0.80 | 0.5278 | 0.5999 | 0.80 / 1.39 / 2.10 / 3.13 / 3.20 / 2.90 |
| 0.60 | 0.3658 | 0.7138 | 0.60 / 1.15 / 1.89 / 3.06 / 3.20 / 2.90 |

The harmonised curves are identical at large `n` and differ only where the
di-interstitial matters — clean isolation. The legacy curve is the outlier, with
an unphysical 8.75 eV hump at `n = 20` that the blend has to fall off a cliff to
correct.

---

## 2. Results

Targets are the Table G.36 neutron fits at 350 °C: `N_ℓ = 5.03e21 m⁻³`
(band 11.6×), `d̄_ℓ = 6.47 nm` (2.3×), `N_c = 2.81e21 m⁻³` (7.2×),
`d̄_c = 2.43 nm` (2.1×).

### 2.1 Binding law (scan 1)

| Run | `E_b(2)` | N_ℓ ratio | d̄_ℓ ratio | N_c ratio | d̄_c ratio | `δ_FP` | nfe | ncfn |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| R01 legacy | 3.92 | ×119 | ×0.33 | ×0.072 | ×0.23 | 1.8e-1 | 6155 | 22 |
| R02 | 2.00 | ×119 | ×0.33 | ×0.072 | ×0.23 | 1.8e-1 | 6030 | 19 |
| R03 | 1.20 | ×119 | ×0.33 | ×0.072 | ×0.23 | 1.8e-1 | 6102 | 18 |
| R04 | 0.80 | ×61 | ×0.30 | ×0.078 | ×0.24 | 1.9e-1 | 5698 | 36 |
| **R05** | **0.60** | **×0.30** | ×0.17 | ×16 | **×1.02** | **5.9e-1** | 5996 | 40 |

### 2.2 Sinks, mobility, production

| Run | change | N_ℓ ratio | d̄_ℓ ratio | N_c ratio | nfe | netf |
|---|---|---:|---:|---:|---:|---:|
| R03 / F2 | reference (`f_cl_i`=0.58, `f_cl_v`=0.15) | ×119 | ×0.33 | ×0.072 | 6102 | 431 |
| R06 | `Z_i` 1.05→1.15 | ×124 | ×0.32 | ×0.072 | 6062 | 450 |
| R07 | `E_m_i` 0.34→0.50 eV | ×119 | ×0.33 | ×0.072 | 6185 | 454 |
| R09 | `ρ_d` 1e13→1e14 m⁻² | ×184 | **×0.41** | ×0.070 | 6819 | 532 |
| F1 | `f_cl_i`=0.25, `f_cl_v`=0.05 (workbook) | ×63 | ×0.35 | ×0.025 | 8418 | 558 |
| F3 | `f_cl_i` 0.25→0.10 | ×49 | **×0.46** | ×0.025 | **11063** | **792** |
| F4 | `η` 0.30→0.20 | ×51 | ×0.35 | ×0.016 | 7805 | 515 |
| **F5** | `f_cl_v` 0.05→0.30 | ×94 | **×0.54** | **×0.150** | 10670 | 784 |

---

## 3. Which parameters matter

**For `d̄_ℓ`** (the tightest band, `s = 0.416`, and wrong in every run):

| rank | parameter | effect on `d̄_ℓ` |
|---|---|---|
| 1 | `f_cl_v` 0.05→0.30 | 2.27 → 3.53 nm (**+55 %**) |
| 2 | `f_cl_i` 0.25→0.10 | 2.27 → 2.96 nm (+30 %) |
| 3 | `ρ_d` ×10 | 2.13 → 2.66 nm (+25 %) |
| — | `η`, `Z_i`, `E_m_i`, `E_b(2)`≥0.8 | ≤3 % |

**For `N_ℓ`** (needs ÷120): `f_cl_i` ÷1.29, `η` ÷1.24, `ρ_d` ×1.55,
`f_cl_v` ×1.49, `Z_i` ×1.04. Stacking the favourable moves buys ~÷1.6.
**Nothing in the admissible set approaches ÷120.**

**Two null results worth keeping.** `E_m_i` 0.34→0.50 eV changes nothing despite
a verified 25× drop in `D_i` (`Di_eff` 1.56e-12 → 6.11e-14 m²/s): loop number
density here is **cascade-production-limited, not diffusion-limited** — clusters
are born from cascades and never dissolve. And `E_b(2)` above ~1 eV is
irrelevant: once the di-interstitial is bound, *how* bound does not matter.

---

## 4. Why the loop density is 120× high — a structural gap

Loops are born from the cascade spectrum and, in this configuration, essentially
nothing removes them: conversion is off, coalescence is weak at these sizes, and
**the loop→network loss channel is off by default and identically zero at the
physical `χ = 1`**. `N_ℓ` therefore accumulates ∝ dose while `d̄_ℓ` stays pinned
near the cascade-born size.

That is exactly the failure mode RadCluster_2_1 objective (a) exists to remove
(`CLAUDE.md` §0a: "SIA loop number density saturates with dose (currently it
grows unbounded — there is no loop↔network coupling)"). **The scan is direct
evidence that objective (a) is required before the loop observables can be
calibrated at all** — no amount of tuning `η`, `f_cl`, `Z_i`, `ρ_d` or the
binding law substitutes for it.

R05 (`E_b(2) = 0.6`) is the only run that reaches the data on any observable
(`d̄_c` ×1.02), and it does so by destroying the loop population and failing
conservation (`δ_FP = 0.59`). It is not a solution; it marks the bifurcation.

---

## 5. Defect found and fixed: the Production sheet was never read

`production_rates()` took `eta`, `f_cl_i`, `f_cl_v`, `s_i`, `s_v`, `i_cascade`,
`v_cascade` from the module-level `FISSION`/`FUSION` dicts in
`defect_production.py`. `InputData.production_fission` was loaded from Excel,
stored, and written into `provenance.md` — but never consulted. The workbook's
Production sheet was decorative.

Detected because R08 (`f_cl_i`) and R10 (`η`) returned **byte-identical solver
counters** to the unperturbed reference.

**Consequences.**
- θ parameters #1–#3 of the twin (`η_FP`, `f_I^cl`, `f_V^cl`) were **unsamplable**.
  Tier-2 screening would have reported `S_i = 0` for all three — indistinguishable
  from a real result.
- It settles the §2.3 `f_cl_i` conflict: the code has always run with the
  Table 2 values **0.58 / 0.15**; the workbook's 0.25 / 0.05 never took effect.
- Every prior result in this repository was produced with 0.58 / 0.15.

**Fix.** `production_rates(..., spec_over=)` merges the workbook over the module
defaults; only present, finite, numeric keys override, so text cells (`'user'`,
`'0.5−1'`) fall back safely. Both call sites
(`rate_equations.py`, `bin_moment_rates.py`) pass the matching sheet.

**Backward compatibility verified:** run F2 (workbook set to 0.58 / 0.15)
reproduces reference R03 bit-for-bit — `N_ℓ = 6.003e23`, `steps = 4310`,
`nfe = 6102`, `ncfn = 18`, `netf = 431`.

> **Default change — needs the author's decision.** With the sheet now live, the
> shipped workbook's `f_cl_i = 0.25`, `f_cl_v = 0.05` take effect and change the
> default result: `N_ℓ` ÷1.9, `N_c` ÷2.9 (F1 vs F2). The workbook, `CLAUDE.md`
> Table 2 (0.58 / 0.15) and the mandated twin priors (`f_cl_i` 0.05–0.25,
> `f_cl_v` 0.2–0.7) are three-way inconsistent. **The workbook has not been
> edited** — whichever set is correct should be written to it deliberately.

---

## 6. Stiffness

Ranked by `nfe` (machine-independent; wall clock is confounded by concurrency):

| nfe | run | driver |
|---:|---|---|
| 11063 | F3 `f_cl_i`=0.10 | most cascade mass into **monomers** |
| 10670 | F5 `f_cl_v`=0.30 | more vacancy clusters |
| 8418 | F1 workbook production | ditto, milder |
| 7805 | F4 `η`=0.20 | — |
| 6819 | R09 `ρ_d`=1e14 | more sink work |
| 6102 | R03/F2 reference | — |
| 5698 | R04 `E_b(2)`=0.8 | — |

**Stiffness is driven by the monomer/small-cluster population, not the binding
law.** Lowering `f_cl_i` (more monomers) nearly doubles the RHS evaluations.
Low `E_b(2)` roughly doubles nonlinear convergence failures (`ncfn` 18 → 36 → 40)
without raising total work — a different, milder mechanism.

For the twin: the Tier-2 LF ensemble will be **most expensive at the low-`f_cl_i`
end of the mandated 0.05–0.25 prior**, which is the opposite of where the current
workbook sits. Budget accordingly.

---

## 7. Recommendations

1. **Commission the loop→network channel before calibrating loop observables.**
   §4. Without it `N_ℓ` cannot be brought within 2 orders of the data by any
   admissible parameter set.
2. **Resolve the three-way `f_cl` inconsistency** (§5) and write the chosen
   values to the workbook.
3. **Treat `E_b(2)` as a threshold parameter**, not a continuous one. Anywhere in
   1.0–3.9 eV is observationally equivalent here; below ~0.8 eV the model
   bifurcates and fails conservation. For the twin, either fix it at ~1.2 eV or
   sample it with the explicit expectation that `S_i ≈ 0` above the threshold.
4. **Re-run this scan with the network channel on** (`χ ≥ 30`, `w_c` swept) —
   that is the parameter most likely to move `N_ℓ` by the required factor, and it
   was not testable in this batch.
5. ~~`δ_FP ≈ 0.18` in every run — far above the gate, needs investigation.~~
   **RETRACTED 2026-07-30.** The `δ_FP` column in §2 is
   `np.nanmax(delta_FP)`, and the maximum always falls on the **first output
   point** (`t = 1e-6` s), where nothing has nucleated and the conservation ratio
   is evaluated on a near-zero denominator. Its value is `1.764e-01` in every
   run *and in every run of two later studies* — invariant to dose, grid,
   binding law and production parameters, because it is a property of the
   initial condition alone.
   The **converged** `δ_FP` for these runs is ~1.5e-07 (measured at the R01
   configuration), comfortably inside the gate. **The conservation of these runs
   is fine; the reporting was not.** See
   [`validation_10dpa_revision3.md`](validation_10dpa_revision3.md) §8. Read the
   §2 `δ_FP` column as "startup transient, constant" and ignore it; `summary.csv`
   in each run directory carries the correct converged value.
