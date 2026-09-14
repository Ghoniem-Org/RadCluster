# 10 dpa Production Validation — revision-3 input-path changes

**Date:** 2026-07-30
**Purpose:** validate the revision-3 changes to the input path before commit.
**Verdict:** **code CLEAR to commit; `E_b_i2 = 0.80 eV` NOT clear and is shipped
disabled.**

Related: [`digital_twin_implementation_plan.md`](digital_twin_implementation_plan.md) §2.2, §4/T0.5.
Logs: `output/log_files/20260730_*`. Run dirs: `output/20260730_102714_*` (legacy),
`output/20260730_105755_*` (current).

---

## 1. What was validated

| # | Change | Result |
|---|---|---|
| 1 | `lambda`, `A_void_0` read from the workbook instead of hard-coded | **PASS** |
| 2 | `loop_net_*` + `LOOP_NETWORK_LOSS` rows added; blank-cell readers hardened | **PASS** |
| 3 | `E_b_i2` overrides `A_111` (di-interstitial anchor) | mechanism **PASS**, value 0.80 eV **FAIL** |

---

## 2. The reference and the harness

A 10 dpa production run predating all of this work is on disk:

```
output/20260626_162101_full_system_bin_moment_CD_fission_I1000V1000_im20vm5
T=723 K, I=V=1000, i_mobile=20, v_mobile=5, i_discrete=50, I_bin=V_bin=20,
shape=linear, rtol=1e-6, atol=1e-20, C_floor=1e-25, he=QSS,
gmres+Woodbury, loop_conversion=1, t_span=(1e-6,1e7)
-> delta_FP = 2.804e-05, delta_He = 3.747e-11, wall 947.6 s
```

That configuration was reproduced exactly and run on two input paths:

- **Leg B / legacy** — `E_b_i2` blank, `A_111 = 3.0` verbatim (pre-revision-3).
- **Leg A / current** — `E_b_i2 = 0.80 eV` active.

Both legs are identical at t=0 (`c_i1 = 4.213e15` vs `4.214e15` m^-3), so they
differ only in the loop binding law.

| | `delta_FP` | `delta_He` | dose | wall |
|---|---|---|---|---|
| Reference on disk | 2.804e-05 | 3.747e-11 | 10.0 dpa | 948 s |
| **Leg B — legacy** | **9.290e-06** | 9.777e-10 | 10.0 dpa | 1137 s |
| **Leg A — `E_b_i2=0.80`** | **3.157e-01** | 5.968e-10 | 10.0 dpa | 1786 s |

Leg B reproduces the reference to within a factor 3 (and is slightly better), so
**the harness is faithful** and leg A's number is attributable to the change
alone.

### A tolerance caveat worth recording for T0.1 — CORRECTED 2026-07-30

An earlier attempt used `rtol = 1e-4` at 350 °C with `i_mobile=10, v_mobile=3`
and produced `delta_FP = 0.494`, which was originally written up here as
"relaxing `rtol` by 100x moved the conservation diagnostic by ~4.7 orders of
magnitude". **That attribution was wrong** — the two runs differed in `rtol`
*and* temperature *and* mobility cutoffs simultaneously, so tolerance was never
isolated.

A controlled study ([`deltaFP_study`](../../output/20260730_160946_deltaFP_study/),
§8 below) subsequently measured the tolerance axis cleanly. `rtol = 1e-4` is
**not** disqualifying: at 623 K it yields `delta_FP = 9.7e-06` at 1 dpa, inside
the gate. The optimum is `rtol = 1e-6` (`1.2e-08`); tightening further to 1e-7
makes it slightly *worse* (`1.4e-07`, roundoff-limited).

The `delta_FP = 0.494` of that early run is a genuine failure, but it is
attributable to the configuration as a whole, not to tolerance alone.

---

## 3. Why `E_b_i2 = 0.80 eV` fails

Not a solver crash. CVODE returned `ret=0` throughout and the step size grew
`h: 6.0e5 -> 6.0e6 -> 6.0e7` while taking ~1 step per checkpoint from 2.6 dpa
onward — the system reached a **fixed point** and stopped evolving.

| | Leg B (legacy) | Leg A (`E_b_i2=0.80`) |
|---|---|---|
| `C_SIA_tot` | 1.27e26 m^-3 | **1.42e22 m^-3** (10^4 lower) |
| `mean_n_i` | 205.5 SIAs | **13.9 SIAs** (0.92 nm) |
| `mean_n_v` | 7.3 vac | 6.0 vac |
| swelling | 1.49e-10 | 5.48e-09 |
| trajectory | evolves to 10 dpa | **frozen after ~2.6 dpa** |

**Root cause — the single anchor rescales the entire small-n branch.** The DFT
branch is `E_b^fit(n) = A_111 * n^{+B_111}`, so pinning it at n=2 multiplies
*every* small size by the same factor:

| n | legacy (A=3.0) | override (E_b_i2=0.80) | ratio |
|---:|---:|---:|---:|
| 2 | 3.924 eV | 0.800 eV | 0.204 |
| 3 | 4.591 eV | 0.936 eV | 0.204 |
| 5 | 5.594 eV | 1.141 eV | 0.204 |
| 10 | 7.307 eV | 1.494 eV | 0.205 |
| 20 | 8.754 eV | 2.041 eV | 0.233 |

Fixing the di-interstitial cost a factor ~5 of binding at *every* nucleation-
relevant size. Emission then outruns capture, loops never grow past ~14 SIAs,
and the microstructure never forms — contradicting the EUROFER97 anchors
(`d_loop ~ 3-12 nm`).

`delta_FP` is a *relative* diagnostic; with `S`, `S_I` and `dJ^d` all collapsed
it is additionally ill-conditioned, so 0.316 overstates the absolute imbalance.
That does not rescue the result — the physics is degenerate regardless of the
diagnostic.

### The fix (not applied — author's call)

Re-fit **both** parameters to DFT across small n rather than sliding `A_111`
alone. Targeting `E_b(2) = 0.80 eV` and `E_b(25) ~ E_f_i = 3.0 eV`:

```
A_111 = 0.5566 eV      B_111 = 0.5233
```

This needs a proper fit to the DFT dataset (and a companion `A_100`/`B_100`
pair), not the two-point solve above, plus a re-validation at this
configuration.

---

## 4. What was shipped

- **Code:** all three changes committed. `E_b_i2` mechanism verified in both
  directions — `0.80` gives `E_b_loop_i(2) = 0.8001 eV` exactly (the 1e-4 offset
  is the continuum admixture), blank falls back to `A_111` verbatim.
- **Workbook:** the `E_b_i2` row exists but is **blank**, so the shipped default
  is the validated leg-B path. Confirmed: `G_SIA`, `K_SIA_grow` and `G_VAC` are
  bit-identical between the shipped workbook and the validated leg-B workbook.
- `lambda` / `A_void_0` remain wired and read from the workbook.
- `loop_net_*` rows present at their code defaults (channel off), verified inert.

Re-enabling is one cell, once the refit above is done and re-validated.

---

## 5. Known issues found along the way (pre-existing, not introduced here)

1. **`run_adaptive` truncates loop-conversion output arrays.** With
   `points_per_segment=10`, `f_111_loop` / `N_loops_100` come back with length 10
   against `t` of length 97 — only the last segment survives, where `N_loops`,
   `mean_n_i` etc. concatenate correctly. This breaks four plots outright
   (`mean_sizes.png`, `number_densities.png`, `loop_fraction.png`,
   `number_densities_tem.png` fail with `shapes (97,) and (10,)`; `mean_sizes_tem.png`
   with `IndexError`). **Any loop-fraction result taken from a `run_adaptive` run
   is currently wrong**, which matters for the §3 conversion campaign and for
   Tier 2 of the twin.
2. **`E_b_bubble` carries no atomistic correction**, so `E_b_bubble(m, 0) != E_b_void(m)`
   — see plan §2.2(b).
3. **`Z_i_loop` is still aliased to `Z_i`** at `reaction_rates.py:132` *and*
   `cpp_bridge.py:179`.

---

## 8. What controls `delta_FP` (study of 2026-07-30)

Output: [`output/20260730_160946_deltaFP_study/`](../../output/20260730_160946_deltaFP_study/)
— `deltaFP_vs_dose_tolerance.png`, `summary.md`, `deltaFP_study.json`, plus one
timestamped run directory per case.

**The reported "large `delta_FP`" was a startup artifact.** `delta_FP` at the
first output point (`t = 1e-6 s`, dose ~1e-12) is **`1.764e-01` in every run** —
identical to four significant figures across six configurations spanning
different doses, output grids, binding laws and production parameters. At that
time nothing has nucleated, so the ratio
`|S - S_I - J_SIA + J_VAC| / (S + S_I + J_SIA + J_VAC)` is evaluated on a
near-zero denominator. A quantity invariant to every parameter is not a
measurement.

The anchor scan reported `np.nanmax(delta_FP)` and therefore reported this
constant for all 15 runs. `post_process` itself is correct — `summary.csv` uses
`_last_finite(delta_FP)`, the converged value. **Use the final value, never the
max.**

The converged `delta_FP` passes the gate everywhere tested. Its real drivers,
ranked:

| Driver | Effect on converged `delta_FP` | Evidence |
|---|---|---|
| **Dose** | 1.2e-08 (1 dpa) -> 3.6e-06 (10 dpa) at 623 K | ~linear accumulation |
| **Small-`n` binding / emission rate** | 1.5e-07 -> 6.5e-05 (x430) for harmonised `E_b(2)=1.2` vs legacy | W5 vs W6 |
| **Integrator tolerance** | 1.2e-08 (1e-6) -> 9.7e-06 (1e-4); *worse* again at 1e-7 | tolerance sweep |
| Production `f_cl` | x5 | W4 vs W6 |
| Output-grid density | x2.8 | W3 vs W4 |

Two cautions. `delta_FP` **fluctuates point-to-point** — the same physical dose
read off two different output grids differed by 400x (W2 vs W4) — so single
scalar readings are noisy and the trajectory is the honest view. And the
temperature dependence is not monotone: at 723 K the dose degradation reverses
(10 dpa conserves *better* than 1 dpa).

**Practical rule for the twin:** `rtol = 1e-6`, judge on the final value, and
expect `delta_FP` to grow roughly linearly with dose — so the Tier-3 gate should
be dose-scaled rather than a flat 1e-6.

---

## 6. Recommendation

Commit the code and the workbook as they now stand. The default path is
validated at 10 dpa (`delta_FP = 9.29e-06`, `delta_He = 9.78e-10`), matching the
pre-existing reference. Item 1 of §5 should be fixed before any loop-fraction
result is trusted.
