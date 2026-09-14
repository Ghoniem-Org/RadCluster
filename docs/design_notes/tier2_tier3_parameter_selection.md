# Parameter selection for Tier 2 (LF screening) and Tier 3 (HF design)

**Date:** 2026-07-31
**Basis:** ~60 runs at 350 °C — anchoring scans, the `δ_FP` tolerance/dose study,
the SIA-grid diagnostic, and two loop→network commissioning scans.
**Output dirs:** `output/20260730_170756_anchor3/`,
`output/20260730_200514_dfp_diag/`, `output/20260731_012347_loopnet/`,
`output/20260731_055301_loopnet/`.

Companion: [`loop_growth_grid_limit.md`](loop_growth_grid_limit.md),
[`anchor_scan_350C_1dpa.md`](anchor_scan_350C_1dpa.md).

---

## 0. The admissibility rule that has to come first

**Grid adequacy is a precondition, not a diagnostic.** Define

```
occupancy = mean_n_i / I          (and the vacancy analogue mean_n_v / V)
```

Runs with `occupancy > ~0.1` are **inadmissible** — their `d̄_ℓ` is a readout of
`I`, not a prediction. Measured, at the runaway base and with the loop→network
channel at four capture widths:

| | `d̄_ℓ` at `I`=1000 → 2000 | ratio |
|---|---|---|
| channel OFF | 4.25 → 5.77 nm | ×1.358 |
| `w_c` = 150 `b₁₁₁` | 4.29 → 5.83 nm | ×1.359 |
| `w_c` = 500 `b₁₁₁` | 4.36 → 5.94 nm | ×1.362 |
| `w_c` = 2000 `b₁₁₁` | 4.01 → 5.50 nm | ×1.372 |

Every ratio is `sqrt(2)`. **The loop→network channel does not bound loop growth**
— not even at 2000× the physical capture width with `ρ_net` driven to its
ceiling. An earlier recommendation in `loop_growth_grid_limit.md` §5 to
commission this channel as the fix is **withdrawn**: it is not sufficient.

**But the runaway is confined to one corner of the prior.** At the normal
production values (`f_cl_i` = 0.25, `f_cl_v` = 0.05) the distribution is
grid-converged: `d̄_ℓ` 2.31 → 2.30 nm across the same grid doubling, occupancy
0.088, `δ_FP` ~1e-8. The runaway appears at low `f_cl_i` / high `f_cl_v`
(occupancy 0.30), which is exactly where the two apparently band-matching
anchoring runs sat — hence their agreement was an artefact.

So the twin does **not** need new physics before proceeding. It needs the
occupancy check enforced, and runs that fail it discarded as numerically
inadmissible (Appendix F step 2) rather than scored. Implement in
`run_ensemble.py`; it is a two-line check that would have prevented every false
result in this campaign.

---

## 1. Tier 2 (LF screening, 0.1 dpa) — recommended active set

Ranked by measured effect on the observables, restricted to parameters that are
identifiable, well-conserved, and affordable at LF.

### Include — strong, clean, cheap

| Parameter | Evidence | Notes |
|---|---|---|
| `i_mobile` | `N_ℓ` ÷240 and `d̄_ℓ` ×3.75 monotone over 5→50 | **Strongest single lever.** But it accelerates the runaway — must be paired with the occupancy check |
| `f_cl_i` | `d̄_ℓ` +30 %, `N_ℓ` ÷1.29 (0.25→0.10) | Also the **stiffest** direction (`nfe` 6102 → 11063). Budget for it |
| `f_cl_v` | `d̄_ℓ` +55 %, `N_c` ×6 (0.05→0.30) | Largest lever on `d̄_ℓ`; drives the runaway at the top of its range |
| `ρ_d` | `N_ℓ` ×1.5, `d̄_ℓ` ×1.24 per decade | Bonus: shortens the transient, `δ_FP` improves ~300×, `nfe` 5000 → 2283 |
| `Z_i` | `N_c` ×3.1, `N_ℓ` ÷2 over 1.02–1.15 | Monotone, `δ_FP` ~1e-8, cheap. Cavity/loop partition |
| `Z_i_loop` | `N_ℓ` ×2.7 and `d̄_ℓ` ×1.8 over 1.05–1.30 | **Only lever that raises loop count and size together** — newly separable (de-aliased 2026-07-30) |

### Exclude from Tier 2 — measured inert or ill-posed at 0.1 dpa

| Parameter | Why |
|---|---|
| `E_m^I` | **No effect** on `N_ℓ` despite a verified 25× drop in `D_i` — loop density is cascade-production-limited, not diffusion-limited |
| `E_b^i(2)` / `A_111`,`B_111` | **Threshold, not a dial.** 3.92 → 2.0 → 1.2 eV changes nothing (<0.1 %); everything happens in a bifurcation below 0.8 eV. Fix at ~1.2 eV or expect `S_i` ≈ 0 |
| `T*` | `ΔF` cutoff sits far above the sizes loops occupy — §3 of the reference measured 250 → 180 °C as bit-identical |
| `w_c`, `χ`, `K_rec` | The whole loop→network block. See §2 — they need dose, and Tier 2 at 0.1 dpa cannot see them |

Excluding four of the plan's 24 leaves **`p_act` ≈ 20** entering Tier 2, of which
6–10 are expected to survive screening.

---

## 2. The loop→network block belongs in Tier 3, not Tier 2

First quantitative characterisation of these parameters (350 °C, 1 dpa,
`δ_FP` 1e-8–1e-9 throughout):

| Parameter | Character | `ρ_net` response |
|---|---|---|
| `w_c` | **linear gain**, no saturation over 1–500 `b₁₁₁` | 1.00 → 1.07 → 1.21 → 1.69 (×10¹⁴) |
| `χ` | **saturating gate**, flat above ~50 | 1.03 (30) → 1.21 (50) → 1.23 (100) |
| `K_rec` | **potent, saturates below 1e-3** | 1.21 (0) → 1.00 (1e-3) → 1.00 (1e-1) |

Three consequences:

1. **`K_rec`'s prior is in the wrong place.** The plan has it log-uniform over
   `1e-3 – 1e1`; `1e-3` already cancels network growth completely, so **the
   entire prior lies in the saturated regime** and the posterior would be flat —
   a non-result that reads as a finding. Move the range to roughly `1e-6 – 1e-2`.
   Note also that `ρ_net` is clamped below at `rho_d`
   (`simulation.py`), so while the gain is weak `K_rec` is *structurally*
   unidentifiable — it can only push down, into the clamp.
2. **`w_c` is inert at its physical value** (`ρ_net` unchanged at `b₁₁₁`). It only
   acts as an amplified effective parameter; that must be stated in
   `parameters.yaml`, not buried.
3. **`χ` and `w_c` interact** — `χ` gates, `w_c` scales. Main-effect `S_i` will
   misrepresent both; this block requires total-effect indices `S_i^T`.

These are dose-driven (`ρ_net` compounds), so screen them at Tier 3 conditions
(≤30 dpa), never at 0.1 dpa.

---

## 3. Tier 3 (HF design) priorities

1. **`Z_i` × `Z_i_loop` jointly.** Their *ratio* sets the loop/network flux
   partition, they move loops and cavities in opposite directions, they are
   cheap and exceptionally well-conserved. Newly independent — this is the first
   design that can exercise them.
2. **`ρ_d`.** The only parameter that improves numerical behaviour while moving
   the observables; drives the system toward quasi-steady state, which shortens
   HF runs inside the 7200 s cap.
3. **The loop→network block** (`w_c`, `χ`, `K_rec`) at full dose, with the
   corrected `K_rec` prior.
4. **`f_cl_i`, `f_cl_v`** — strongest on the observables but the stiffest and the
   runaway drivers. Admissible only with the occupancy check enforced.

**Grid sizing must be per-region, not global.** Occupancy is parameter-dependent
in *both* directions: raising the clustered fractions drives loops toward the `I`
ceiling; lowering them starves nucleation so the few cavities that form grow
toward the `V` ceiling (measured: mean cavity 654 against `V` = 1000). A single
fixed LF/HF grid cannot serve the whole prior. Either size `I`, `V` per sampled
point from a cheap pre-pass, or use `run_adaptive` domain doubling — which first
requires fixing its loop-conversion array-truncation bug.

---

## 4. Open items

- **Nothing bounds loop growth in the runaway corner.** Not new physics for the
  twin to add, but the corner must be excluded by the occupancy rule, and the
  reason recorded: with `i_mobile` large, coalescence merges loops without limit
  and no implemented mechanism removes them fast enough.
- **`run_adaptive` truncates loop-conversion output arrays** (`f_111_loop`,
  `N_loops_100` come back length 10 against `t` of 97). Any loop-fraction result
  from an adaptive run is currently wrong.
- **`E_b_bubble` carries no atomistic correction**, so `E_b_bubble(m,0) ≠ E_b_void(m)`.
- The two anchoring runs that appeared to land in the experimental band
  (`C_fi0.05_fv0.7`, `E2_combo`) are **artefacts** and must not be used as
  starting points.
