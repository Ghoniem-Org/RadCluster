# Why d_111 does not converge with bin count

**2026-09-09.** Campaign M2R showed five of six observables converging to under
0.3% against the exact arm while `d_111` sat at −5 to −6.6% and did **not**
improve with bin count — the finest rung (B32, N_eq 238) was the worst of nine.
This is the investigation.

## Result

**The ½⟨111⟩ mean-size error is controlled by `i_discrete`, not by `I_bin`.**
The M2R ladder held `i_discrete = 50` fixed and varied only the bin count, so it
turned the one knob that does not control this error.

At fixed `I_bin = 16`, I = V = 1000, 0.05 dpa, exact `d_111` = 2.9787 nm:

| `i_discrete` | `d_111` (nm) | deviation | % of ½⟨111⟩ **content** above `i_discrete` |
|---|---|---|---|
| 25  | 2.7232 | −8.6% | 96.6% |
| 50  | 2.7744 | −6.9% | 91.2% |
| 100 | 2.8128 | −5.6% | 83.4% |
| 200 | 2.8692 | −3.7% | 70.6% |
| 400 | 2.9403 | −1.3% | 47.7% |
| 800 | 2.9813 | **+0.1%** | 8.0% |

The error tracks the fraction of the population the closure has to represent.
It is ordinary closure error — it simply converges along the axis that was held
fixed.

Conversely, at fixed `i_discrete = 400`, varying `I_bin`:

| `I_bin` | 4 | 8 | 33 |
|---|---|---|---|
| `d_111` (nm) | 2.9615 | 2.9503 | 2.9245 |

More bins make it slightly **worse**. Refining the bins does not refine this.

## Why ½⟨111⟩ and not the others

The exact ½⟨111⟩ mean size is ~147 SIAs — roughly **3× `i_discrete = 50`** — so
the distribution's peak sits inside the binned region and the closure carries
it. `d_100` and `d_cavity` show the same monotone improvement with `i_discrete`
(−0.67% → +0.05%, and +0.52% → −0.01% across the same sweep) but with ~10× less
amplitude. Nothing is special about ½⟨111⟩ except how far its characteristic
size sits above the discrete core.

## What was ruled out

**Not a code-path difference.** With `I_bin = 0, i_discrete = I` — the
bin-moment path carrying no closure — it reproduces the discrete arm at 20 dpa
to within rounding on all six observables (`d_111` −0.000%, `N_111` −0.002%,
`delta_FP` +0.001%). `check_bm_vs_discrete.py` asks this at 1e-2 dpa; this
extends it three decades, to the dose where the deviation is measured.

**Not `LOOP_COAL`.** It is the one channel implemented twice (per-size pair sum
in `rhs_case2`, class-pair moment sum in `rhs_bin_moment`) and the only one that
moves content up in size space by large jumps, so it was the leading suspect.
With the channel **off in both arms** the deficit is unchanged: `d_111` −7.3%
(B8) and −7.0% (B32) at 0.05 dpa, against −7.0% and −7.9% with it on.

**Not a domain-edge artefact.** The exact spectrum decays smoothly to the top
of the grid (size 1000 holds 0.01% of the content, no pile-up), so `loop_coal`'s
edge-piling convention is not implicated and I = 1000 is adequate for ½⟨111⟩.

**Not slow accumulation.** The deficit is already −8.4% at 0.005 dpa, the first
meaningful output point, and drifts only to −6.6% by 20 dpa.

## Bin-by-bin evidence

Summing the exact per-size spectrum onto the binned run's own bin edges — moments
against moments, no reconstruction anywhere, so the closure is not measured
twice (B32, 20 dpa):

| sizes | μ₀ deficit | μ₁ deficit |
|---|---|---|
| 94–449 | −1.3 … −3.6% | tracks μ₀ |
| 493–541 | −4.8% | −4.9% |
| 594–652 | −13.5% | −13.6% |
| 715–785 | −31.7% | −31.8% |
| 862–946 | −72.4% | −72.7% |

Sizes ≥ 500 carry **44% of the ½⟨111⟩ content**, so a deficit there lands
directly on the mean.

## A separate bug found on the way

**`C_SIA_tot` is inflated 17–25× in fine-bin runs.** When the top bin empties,
its μ₀ decays to the C_floor level while its μ₁ does not, leaving an
inconsistent moment pair in the *solver state*:

| run | I_bin | last μ₀ | last μ₁ | μ₁/μ₀ | `C_SIA_tot` (exact 3.39e24) |
|---|---|---|---|---|---|
| B8  | 8  | 1.08e-09 | 8.65e-07 | 799 | 3.25e24 |
| B16 | 17 | 2.11e-12 | 9.54e-04 | 4.5e8 | **8.42e25** |
| B32 | 33 | 5.46e-12 | 6.47e-04 | 1.2e8 | **5.82e25** |

A mean size of 1.2e8 in a bin spanning 946–1001 is impossible. `post_process`
already guards the *mean-size* path — `_floor_bin_moments` clamps μ₁/μ₀ back
inside the bin, a guard added after `mean_n_100` once reached 2853 on a grid of
1000 — but `SIA_content_from_mu1` sums the **raw** μ₁, so the corruption flows
straight into `C_SIA_tot`. The guard hides the symptom in one consumer while
another drinks it neat.

Scope: the six observables and `delta_FP` are **not** affected (they take the
clamped path, and `delta_FP` is not formed from `C_SIA_tot`); `swelling` is the
vacancy inventory and is also unaffected. What is wrong is `C_SIA_tot` itself,
which is exported and is `S_I` in the swelling identity S = S_I + ΔJ^d
(CLAUDE.md §8). **Not fixed here** — reported only.

## Consequences for the study

1. **The M2R convergence table is valid but mislabelled.** It is a bin-count
   convergence study, and five observables converge in bin count. `d_111` does
   not, because bin count is not its convergence axis.
2. **A binning ladder must sweep `i_discrete` as well.** Table 4's rungs varied
   both together, which confounded them; M2R varied only `I_bin`, which hid one.
   Neither design can attribute the error. Sweep them independently.
3. **`i_discrete` should be chosen against the population's mean size**, not as
   a fixed number. `i_discrete = 50` with a mean of 147 puts the peak of the
   distribution inside the closure.

## Reproducing

* `codes/Python_Testing/check_idiscrete_convergence.py <i_discrete> <I_bin>`
* `codes/Python_Testing/check_binwise_moments.py <I_bin>`
* `codes/Python_Testing/check_path_equivalence_20dpa.py`
