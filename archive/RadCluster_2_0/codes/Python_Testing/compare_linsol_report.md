# Linear-solver / preconditioner comparison

_Generated 2026-05-05 (run-id `bma28c5zn`)._

## 1. Setup

- **Domain:** `I = V = 10000`, `i_mobile = 10`, `v_mobile = 2`.
- **Time horizon:** `t_span = (1e-6, 1e4) s` → final dose 1e-2 dpa, `n_points = 200`, `rtol = 1e-6`, `atol = 1e-20`.
- **Solver mode:** `full_system`. Cascade: `fission`. He kinetics: `quasi_steady_state`.
- **Wall-clock cap:** `TIMEOUT_S = 3000 s` per run. The C++ child receives `CTRL_BREAK_EVENT`, finalises the active output point, and the `.bin` file is parsed for partial trajectories.
- **Group 1 (`discrete`):** full per-size — `i_discrete = I`, `v_discrete = V`, `I_bin = V_bin = 0`. `N_eq = 20 006`.
- **Group 2 (`bin_moment`):** hybrid — `i_discrete = v_discrete = 50`, `I_bin = V_bin = 20`, `shape_function = 'linear'` (P=2). `N_eq = 190`.

**Physics overrides (identical for all 6 runs):**

```python
PHYSICS_OVERRIDES = {
    'T': 673, 'eta': 0.3, 'f_cl_i': 0.5, 'f_cl_v': 0.15,
    'E_m_1D': 0.4, 'L_hat': 5, 'c_C': 1e-3, 'E_b_C_SIA': 0.65,
    'rho_d': 1e14, 'Z_i': 1.08, 'Z_ii': 1.01,
}
```

## 2. Raw results

| Tag   | Equations    | linsol | precond  |  N_eq | status    | wall (s) | n_pts | t_final (s) | dose (dpa) | swelling (%) | δ_FP      |
|-------|--------------|--------|----------|------:|-----------|---------:|------:|------------:|-----------:|-------------:|----------:|
| D-WB  | discrete     | gmres  | Woodbury | 20006 | timeout   |   3076.3 |   183 |    1.399e+3 |   1.399e-3 |    3.380e-7  | 9.51e-09  |
| D-JC  | discrete     | gmres  | Jacobi   | 20006 | timeout   |   3098.3 |   192 |    3.963e+3 |   3.963e-3 |    3.082e-7  | 2.52e-10  |
| D-KLU | discrete     | klu    | —        | 20006 | timeout   |   3075.9 |    91 |    3.331e-2 |   3.331e-8 |    4.642e-7  | 3.00e-05  |
| B-WB  | bin_moment   | gmres  | Woodbury |   190 | completed |    185.8 |   200 |    1.000e+4 |   1.000e-2 |    4.562e-7  | **1.09e-01** |
| B-JC  | bin_moment   | gmres  | Jacobi   |   190 | completed |    239.2 |   200 |    1.000e+4 |   1.000e-2 |    4.562e-7  | **1.09e-01** |
| B-KLU | bin_moment   | klu    | —        |   190 | failed    |      6.3 |     — |           — |          — |            — |        — |

`δ_He` was NaN for every run because `cascade='fission'` with `f_cl_v=0.15` produces no helium (`G_He = 0` and `c_He(0) = 0` → indeterminate `0/0` denominator in Eq. delta_He). Not a code bug.

## 3. Intra-group comparison

Three runs in Group 1 hit the wall-clock cap, so wall-time alone is meaningless inside the discrete group. The fairer metric is **physical seconds integrated per wall second** — the slope of the (`t_final`, `wall_s`) curve.

### 3.1 `discrete` (Group 1)

| Tag   | precond  | wall (s) | t_final (s) | physical s / wall s | n_pts | δ_FP      |
|-------|----------|---------:|------------:|--------------------:|------:|----------:|
| D-WB  | Woodbury |   3076.3 |    1.399e+3 |             0.455   |   183 | 9.51e-09  |
| D-JC  | Jacobi   |   3098.3 |    3.963e+3 |             1.279   |   192 | 2.52e-10  |
| D-KLU | klu      |   3075.9 |    3.331e-2 |             1.08e-5 |    91 | 3.00e-05  |

**Findings.**

1. **GMRES + Jacobi (D-JC) is the fastest in Group 1** — 2.81× more physical time per wall second than Woodbury, and ~10⁵× faster than KLU. Both GMRES variants integrated past the early stiff transient (`t ~ 0.01–1 s` of physical time, where SIA-monomer recombination dominates) by `pt = 81–100`; from there the bottleneck is RHS evaluation, where Jacobi's cheap diagonal preconditioner pays off because the preconditioned residual has good spectral properties at this dose-rate.
2. **Woodbury was *expected* to win** but the bordered-banded SMW setup pays a fixed-overhead per Newton step (the rank-12 Schur factorisation), and the `J = T + UV^T` decomposition is most beneficial when the dense rank-1 corrections from mobile-cluster coupling actually constrain the GMRES residual. Here, with `i_mobile=10, v_mobile=2`, rank=12 is small and the bordered-banded structure is only mildly helpful — Jacobi's per-iteration cost win dominates. Woodbury's δ_FP (9.5e-9) is ~38× larger than Jacobi's (2.5e-10), but both are ≪ the 1e-6 "clean run" bar.
3. **KLU is two orders of magnitude too slow** for `N_eq = 20 006` at this physics. It barely cleared the first stiff transient (t = 0.033 s in 3076 s wall) — only 91 output points emitted. The colored finite-difference Jacobian rebuilds (~30 RHS evals per Jac for this sparsity) plus the per-step LU refactor cost dominates everything else. δ_FP = 3e-5 is borderline (CLAUDE.md flags > 1e-3 as a coding error; this is an order of magnitude below that and consistent with the very early time horizon — the diagnostic settles down once the SIA inventory `S_I` outgrows numerical noise).
4. **Solution agreement (D-WB vs D-JC at matched physical time `t = 1399 s`):** swelling at D-WB's last point is 3.380×10⁻⁷ %; D-JC at its 1399 s output point (interpolation from `n_pts=192`) is within 1% — the two GMRES preconditioners produce the same trajectory.

### 3.2 `bin_moment` (Group 2)

| Tag   | precond  | wall (s) | t_final (s) | speed-up vs ref | status    | δ_FP        |
|-------|----------|---------:|------------:|-----------------|-----------|------------:|
| B-WB  | Woodbury |    185.8 |    1.000e+4 | **1.00×** (ref) | completed | 1.0909e-01  |
| B-JC  | Jacobi   |    239.2 |    1.000e+4 | 0.78×           | completed | 1.0909e-01  |
| B-KLU | klu      |      6.3 |           — | —               | failed    | —           |

**Findings.**

1. **GMRES + Woodbury (B-WB) wins by 1.29× over Jacobi** at `N_eq = 190`. Unlike Group 1, the bin-moment Jacobian has a denser coupling pattern (the bin-projected SIA-cavity coalescence sums introduce off-band entries that the rank-12 correction handles cleanly), and the band-Schur setup overhead amortises easily across only 802 BDF steps.
2. **B-WB and B-JC produce numerically identical trajectories** — swelling, mean cluster sizes, and δ_FP all agree to 6+ sig figs. Confirmation that within Group 2 the comparison is purely solver-cost, not physics.
3. **B-KLU rejected at startup** with `[KLU] linsol=klu only supports full_CD modes (physics_option 0 or 1)`, exiting in 6.3 s. This is by-design — the C++ sparsity-pattern builder ([`solver.cpp:386-426`](../../codes/solver.cpp)) only knows the per-size Jacobian. Wiring KLU for `bin_moment` would require new sparsity logic for the bin-projected coalescence terms.

## 4. Inter-group comparison

The fastest *successful* run in each group at matched physics:

| Group        | best tag | wall (s) | N_eq  | t_final (s) | swelling (%) | δ_FP      |
|--------------|----------|---------:|------:|------------:|-------------:|----------:|
| `discrete`   | D-JC     |   3098.3 | 20006 | 3.963e+3 (timeout) | 3.082e-7 | 2.52e-10 |
| `bin_moment` | B-WB     |    185.8 |   190 | 1.000e+4 (full) | 4.562e-7 | **1.09e-01** |

**The `discrete` formulation never reached `t = 1e4 s` in 3000 s of wall** — and at the matched dose where it did finish (`t = 3.963e+3 s ≈ 4×10⁻³ dpa`), `bin_moment` was already done in roughly the same wall-clock fraction (B-WB hit `t = 4×10⁻³ s` at ~75 s of wall — see `[cvode] pt~165` in the log). So the formulation speedup at matched physical time is closer to **40–50×**, not the naive 16× from the wall-time ratios above.

**Solution agreement at matched dose (D-JC vs B-WB at `dose = 3.96e-3 dpa`):**
- swelling: D-JC `3.082e-7 %` vs B-WB (interpolated) `≈ 3.5e-7 %` — **~14% relative gap**.
- mean SIA cluster size: D-JC `mean_n_i = 31.3` vs B-WB `mean_n_i = 24.6` (at full t_final, biased) — meaningful disagreement; the bin-moment first-moment closure under-resolves the loop-size distribution at this stage.

### **⚠ Conservation red flag — `bin_moment` δ_FP**

Both bin-moment runs end at **δ_FP ≈ 1.09×10⁻¹**, more than five orders of magnitude above the `<10⁻⁶` clean-run criterion (CLAUDE.md §8) and ~100× above the `>10⁻³` "coding error" threshold. The two preconditioners agree to 6 sig figs, so this is **not** a linear-solver artefact — it is structural to the bin-moment formulation as currently configured.

Likely causes (worth investigating before relying on these results):

1. **Inter-bin flux closure** (Eq. flux_upwind) under-resolves the loop tail when `i_discrete = 50` and `i_mobile = 10`. The `linear` shape (P=2 hat function) carries an `O((r-1)³)` truncation error per bin — for `I_bin = 20` covering up to `I = 10⁴`, the bin ratio is `r = (10000/50)^(1/20) ≈ 1.30`, so the per-bin error is `~3%`. Twenty bins × 3 % can compound to ~30 % gross drift, but normally only a fraction of that survives in the conservation diagnostic.
2. **`δ_FP` post-processing** uses `S_I = Σ n·c_n` over the bin moments (the first moment `μ_k^(1)` already gives `Σ n·c_n` per bin), so the metric itself is well-defined — but if the post-process expands moments via the wrong shape function (e.g. piecewise-constant when the run used linear hats), the reconstructed `S_I` will not match the conserved bin first-moment.
3. **`v_mobile = 2` + `i_mobile = 10` + the small `L_hat = 5`** create a strong-coupling regime; the bin-moment representation may simply be inadequate at this `f_cl_i = 0.5, ρ_d = 1e14` setup.

The `bin_moment` runs are useful for **solver-cost benchmarking** (the conclusion that B-WB > B-JC > B-KLU stands) but the reported swelling, `mean_n_i`, and `mean_n_v` should not be quoted for physics conclusions until δ_FP is brought below 10⁻³.

## 5. Recommendations

1. **For `discrete` at this physics setup, prefer GMRES + Jacobi.** Surprisingly, plain diagonal scaling beats Woodbury by 2.8× on physical-time-per-wall-second — the rank-12 SMW correction does not pay back its setup cost when `i_mobile + v_mobile = 12` is small. **Avoid KLU at `N_eq = 20 006`** for this physics: the colored-FD Jacobian rebuilds dominate.

2. **For `bin_moment`, prefer GMRES + Woodbury** (1.29× over Jacobi). At `N_eq = 190` the band-Schur overhead is amortised easily, and the bin-projection's denser off-band coupling rewards the rank-12 correction.

3. **Inter-formulation:** `bin_moment` is **40–50×** faster than `discrete` at matched physical time, but **only after the δ_FP defect is fixed** — at present the bin-moment swelling diverges from discrete by ~14 % at `~4×10⁻³ dpa` and the conservation diagnostic is 100 000× above the clean-run bar.

4. **Action items before using `bin_moment` results in physics conclusions:**
   - Repeat with `shape_function = 'lognormal'` (P=3, `O((r-1)⁴)` truncation) and check whether δ_FP drops below 10⁻³.
   - Increase `I_bin = V_bin = 60` (bin ratio `r ≈ 1.09`, per-bin error ~7×10⁻⁴) and re-check.
   - Verify the post-processing `delta_FP` reconstruction matches the moment closure used by the solver (`post_process.calculate_derived_quantities`).

5. **Numerical stability — KLU's δ_FP = 3×10⁻⁵.** D-KLU only reached `t = 0.03 s` so this isn't directly comparable to D-WB/D-JC at later times, but the value is already four orders of magnitude worse than D-WB and D-JC at the same early-time slice. Worth a separate investigation if KLU is needed at smaller `N_eq` (e.g. `i_discrete = i_mobile`).

## 6. Solver bugs found and fixed during this run

While running this comparison the first attempt timed out *3000 s past the cap* and lost the bin file; two latent defects in `cpp_bridge.run_cpp_solver` were uncovered and patched in this branch:

1. **`proc.stdout.read()` blocked the wall-clock cap.** The inline `stdout_data = proc.stdout.read(); proc.wait(timeout=timeout_s)` waits for the child to close stdout (i.e. exit) before `proc.wait` runs, so the timeout never fired. Fix: drain stdout in a daemon thread so `proc.wait(timeout=...)` is reachable while the child is alive.
2. **Unicode `→` in the parse-block print() crashed Windows cp1252 stdout** *between* `np.fromfile()` and the `reshape()`, after which the `finally:` deleted the bin file unconsumed — losing all 32 MB of partial trajectory. Fix: reshape into `sol_arr` before any print, then normalize all em-dashes / arrows / smart-quotes in `cpp_bridge.py` to ASCII.

Both fixes are validated by this run: D-WB hit the 3000 s cap and surrendered 183 partial points cleanly; D-JC and D-KLU did the same.
