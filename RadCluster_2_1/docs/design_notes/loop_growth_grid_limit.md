# Loop growth is unbounded — `d̄_ℓ` is set by the SIA grid, not by physics

**Date:** 2026-07-30
**Status:** blocking finding for the anchoring exercise and for Tier 1 of the twin.
**Data:** `output/20260730_170756_anchor3/` (17-run anchoring scan),
`output/20260730_200514_dfp_diag/` (8-run conservation diagnostic).

---

## 1. Finding

At the anchoring conditions (350 °C, 1 dpa, `i_mobile = 50`, conversion OFF,
loop→network loss OFF) the SIA loop distribution **has no intrinsic upper
bound**. It grows until it reaches the top of the tracked size grid `I`, and
then loses content off the boundary.

Mean loop size scales **exactly linearly with `I`**, and `d̄_ℓ` as `sqrt(I)`:

| `I` | mean loop size | ratio | `d̄_ℓ` [nm] | ratio | `N_ℓ` [m⁻³] | `δ_FP` (SIA arm) |
|---:|---:|---:|---:|---:|---:|---:|
| 1000 | 299 | — | 4.25 | — | 3.44e21 | 1.86e-01 |
| 2000 | 550 | 1.84 | 5.77 | 1.36 | 2.71e21 | 1.75e-01 |
| 4000 | 1018 | 1.85 | 7.85 | 1.36 | 2.10e21 | 1.65e-01 |

`1.36 ≈ sqrt(1.85)`, exactly as expected for `d ∝ sqrt(n)`. Conservation does not
improve with grid size — content flows off the top at whatever ceiling is set.

**Consequence: `d̄_ℓ` and `N_ℓ` from this configuration are readouts of `I`, not
predictions.** Any agreement with the experimental band is coincidental — at
`I = 4000`, `d̄_ℓ = 7.85 nm` *overshoots* the 6.47 nm target, so some intermediate
`I` would match it exactly and mean nothing.

---

## 2. What it is NOT

The 8-run diagnostic (config `C_fi0.05_fv0.7`: `i_mobile=50`, `rho_d=1e14`,
`f_cl_i=0.05`, `f_cl_v=0.70`) separated the SIA and vacancy conservation arms and
varied each grid independently:

| Candidate | Test | Verdict |
|---|---|---|
| Integrator tolerance | `rtol` 1e-6 → 1e-7 | **ruled out** — bit-identical (`δ_FP` 2.34e-01 both) |
| Vacancy grid | `V` 1000 → 4000 at `I`=1000 | **ruled out** — bit-identical |
| SIA grid | `I` 1000 → 4000 at `V`=1000 | **CONFIRMED** — mean loop ∝ `I` |
| `f_cl_v` | control at 0.05 vs 0.70 | **ruled out as a cause** — a confound |

The `f_cl_v` correlation seen in the anchoring scan was **spurious**. High
`f_cl_v` grows loops (`d̄_ℓ` 2.31 → 4.25 nm), and bigger loops hit the ceiling
sooner — that is the whole of the effect. The failing arm simply follows whichever
population is pressing its own ceiling:

| Run | mean loop / `I` | mean cavity / `V` | failing arm |
|---|---|---|---|
| `f_cl_v` = 0.70 | **299** / 1000 | 60 / 1000 | SIA, 1.9e-01 |
| `f_cl_v` = 0.05 control | 12 / 1000 | **654** / 1000 | vacancy, 7.4e-03 |

So grid adequacy is **parameter-dependent in both directions**: raising the
clustered fractions drives loops toward the `I` ceiling; lowering them starves
nucleation so the few cavities that form grow toward the `V` ceiling. **No single
fixed grid serves the whole prior.**

---

## 3. Two runs discarded

`G2` (`I=V=4000`, `rtol=1e-6`) and `T2` (`I=V=4000`, `rtol=1e-7`) both hit the
2400 s wall cap, reporting `steps = 0, nfe = 0` — the solver was killed before
emitting its statistics. Their `δ_FP` of 9.3e-09 and 7.0e-13 reflect *stopping
early at low dose*, not good conservation, and their `d̄_ℓ` of 11.26 and 6.36 nm
are snapshots at different unknown doses. **Both are excluded.** The `I`-ladder
above uses `G3` (`I=4000, V=1000`), which completed in 1767 s.

---

## 4. Why this blocks the anchoring exercise

The 17-run anchoring scan found two configurations inside the experimental band
(`C_fi0.05_fv0.7`, RMS 1.49; `E2_combo`, RMS 1.61). **Both are grid artefacts**
— they were run at `I = 1000`, where the ceiling happened to truncate the loop
distribution near the measured size.

More fundamentally: the model as configured contains **no mechanism that
terminates loop growth**. Coalescence (`i_mobile`) merges loops without limit;
nothing removes large ones. `i_mobile` looked like the cleanest lever in the scan
(`N_ℓ` ÷240, `d̄_ℓ` ×3.75, monotone over 5→50) precisely because it accelerates an
unbounded process whose end point the grid was silently supplying.

**Until a loop-terminating mechanism exists, `d̄_ℓ` and `N_ℓ` cannot be calibrated
at any grid size.**

---

## 5. Recommendation

1. **Commission the loop→network loss channel (objective (a)) before any further
   anchoring.** This is the mechanism whose absence the finding identifies. It is
   already implemented and parameterised (`loop_net_w_c`, `loop_net_chi`,
   `loop_net_K_rec`, `loop_net_rho_max`) but is off by default and identically
   zero at the physical `χ = 1` — see the digital-twin plan §2.4. Its calibration
   is no longer optional tuning; it is what makes the loop observables
   well-posed.
2. **Add a grid-adequacy assertion to every run**: flag when the tracked
   population's upper tail carries non-negligible weight at the top bin. A run
   whose distribution touches the ceiling is not a physical result. This is
   cheap and would have caught the whole problem immediately.
3. **Re-run the anchoring scan after (1)**, with the grid check from (2) active.
   Levers that remain meaningful and well-conserved in the present data —
   `Z_i` (cavity/loop partition, monotone, `δ_FP` ~1e-8), `Z_i_loop` (raises loop
   count and size together), `ρ_d` (shortens the transient, improves conservation
   ~300×) — are unaffected by this finding and their measured responses stand.
4. For the twin, T0.2 must define LF/HF grids **per parameter region**, not
   globally, or use `run_adaptive`'s domain doubling — which requires first
   fixing its loop-conversion array-truncation bug
   ([`validation_10dpa_revision3.md`](validation_10dpa_revision3.md) §5).
