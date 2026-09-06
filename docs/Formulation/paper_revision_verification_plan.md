# Plan — Verification & approximation study for *Generalized Cluster Dynamics*

**Purpose:** answer the reviewer's §6.2 / Figure 8 critique with a systematic
study of (a) physics approximations and (b) numerical approximations, measured
against an unapproximated reference, scored on the six experimentally
observable parameters.

**Status:** plan only. No runs launched, no code changed.
**Anchor:** the tracked reference run
`output/20260906_063055_full_system_bin_moment_CD_fission_I80000V20000_im5vm5`.

---

## 0. What the reviewer is right about — and the answer already in the repo

The reviewer makes three claims. Two are correct, and the third has a clean
measured answer that the paper simply never reported.

### 0.1 "Why would a preconditioner change results by an order of magnitude?"

It does not, and we can prove it to nine significant figures.

`prec_type` is read at `cpp_utils/core/parameters.h:706` and used only by
`prec_setup_*` / `prec_solve_*` in `cpp_utils/core/rhs_dispatch.cpp`. It selects
the right-preconditioner for SPGMR (`solver.cpp:405`) and enters no rate kernel,
no Jacobian entry and no state-vector layout. It cannot change a converged
answer, and the existing study confirms it.

From `codes/Python_Testing/compare_linsol_results.json`, the one pair where
**both arms actually ran to completion** (bin_moment, N_eq = 190, t_final =
10000 s, dose 0.01 dpa):

| tag | preconditioner | status | mean_n_i | mean_n_v | swelling % | wall (s) |
|---|---|---|---|---|---|---|
| B-WB | Woodbury | completed | 24.56392587969788 | 5.341081396954862 | 4.562396098e-7 | 185.8 |
| B-JC | Jacobi | completed | 24.56392589863494 | 5.341081573868007 | 4.562396216e-7 | 239.2 |

Agreement to 9 significant figures in `mean_n_i`, 7 in `mean_n_v`, 8 in
swelling. **Only the wall-clock differs (1.29x).** This is the null result the
paper should have shown, and it is exactly what the reviewer expects.

### 0.2 So where did the order of magnitude come from?

From plotting runs that never reached the same dose. The generating script says
so in its own docstring (`codes/Python_Testing/compare_linsol.py:7-9`):

> "Each run is capped at TIMEOUT_S wall-clock seconds; **if the cap is hit, the
> C++ solver finalises the current output step gracefully and the partial
> trajectory is kept.**"

The three discrete arms of that study all hit the cap, at wildly different
places:

| tag | linsol / prec | status | t_final (s) | dose reached (dpa) | mean_n_i |
|---|---|---|---|---|---|
| D-WB | GMRES + Woodbury | timeout | 1398.7 | 1.40e-3 | 31.18 |
| D-JC | GMRES + Jacobi | timeout | 3962.7 | 3.96e-3 | 31.32 |
| D-KLU | KLU direct | timeout | 0.033 | 3.33e-8 | 8.19 |

D-JC reached **2.83x the dose** of D-WB; D-KLU stopped **five decades** short.
Curves drawn from these and laid on one axis differ by orders of magnitude
because they are at different doses, not because the linear algebra differs.
The reviewer's instinct — "this suggests something is fundamentally wrong or
not properly presented" — is correct, and the answer is *not properly
presented*.

Note also that D-WB and D-JC agree on `mean_n_i` to 0.4% **despite** the 2.8x
dose difference, which is further evidence the preconditioner is inert.

### 0.3 Secondary defects in §6.2 to fix while we are there

- **Caption/text dose mismatch.** The Fig. 8 caption says "at 0.1 dpa" for both
  panels; the text (p. 50) describes a "10 dpa horizon" and quotes a 9198 s
  discrete reference reaching 10 dpa. Neither dose matches
  `compare_linsol_results.json`, whose `t_span` is (1e-6, 1e4) s = **0.01 dpa**.
  At least two different studies are being described as one figure.
- **Conflated axes.** One legend mixes discretization (discrete vs bin-moment),
  closure order (P = 1/2/3), linear solver (GMRES/KLU), preconditioner
  (Woodbury/Jacobi) and solver mode (sliding window). These live at different
  levels; the reviewer is right that no conclusion can be drawn.
- **Unreported non-convergence.** `delta_FP` = 0.109 for B-WB/B-JC — a 10.9%
  Frenkel-pair conservation error against the model's own 1e-2 gate. Those two
  arms agree with each other but are not converged physics. Agreement between
  preconditioners proves the linear-algebra null; it does not license the run
  as accurate.
- **Metric mismatch.** §6.2 scores `N_loops`, `<n_i>`, `<m_v>`, `C_He`, while
  §6.1 (Figs. 5, 7) scores densities and sizes against TEM. The study should be
  scored on the same six observables as the rest of the paper.

**Deliverable from §0 (no new runs):** a short subsection stating the
preconditioner null result with the table above, and a retraction/replacement of
the old Fig. 8b.

---

## 1. The reference, and the affordability wall

### 1.1 What "without approximations" means here

The unapproximated calculation is: **fully discrete** (one ODE per cluster
size, no bin-moment closure), **per-size helium** (Case 1, `Q_m` resolved per
cavity class), **dynamic free helium** (not QSS), full mobility, no sliding
window, on a domain large enough that neither population piles at the ceiling.

### 1.2 The production reference we are anchoring to

From the tracked run's `provenance.md`:

| | |
|---|---|
| equations | `bin_moment`, shape `linear` (P = 2) |
| grid | `i_discrete` 100, `I_bin` 18, `v_discrete` 5, `V_bin` 20 |
| domain | I = 80000, V = 20000 |
| mobility | `i_mobile` 5, `v_mobile` 5 |
| solver | `full_system`, GMRES + Woodbury, rtol 1e-5 |
| physics | fission, G = 1e-7 dpa/s, 330 °C, `loop_conversion` 1, `lnl` 1, `VOID_NETWORK_LOSS` 1 |
| cost | **1317 s** at 12 OpenMP threads, to 40 dpa |

N_eq ~ 200. The equivalent **discrete** system at this domain is
N_eq ~ 100,000.

### 1.3 Why a discrete reference at the production domain is impossible

Measured anchor (`results/G5a_disc_V1k.manifest.json`), discrete at
**I = 4000, V = 1000**:

- **24,471 s** of wall-clock, and it reached only **6.06 of 15 dpa** before the
  21,600 s budget expired.
- The matched bin-moment run (`G5c_bm_V1k`) did the full 15 dpa in **59.9 s** —
  a **409x speed-up**.

Scaling to the production domain is worse than linear. `prec_setup_woodbury`
(`rhs_dispatch.cpp:168-204`) rebuilds the banded block by finite differences,
one RHS evaluation per column — **O(N) RHS evals per Jacobian setup, hence
O(N^2) overall**. Going from N ~ 5,000 to N ~ 100,000 is 20x in N and ~400x in
setup cost. A discrete run at I = 80000 / V = 20000 is on the order of **1e7 s**
and would still not finish. It is out of reach and will stay out of reach.

**Consequence — state this in the paper.** The verification study is performed
on a *reduced verification domain* where the discrete reference is affordable,
and the transfer of its conclusions to the production domain is established
separately (§2.3). Pretending otherwise is what got Fig. 8 into trouble.

### 1.4 The exact-convergence trick that makes the ladder rigorous

The bin-moment closure applies **only above `i_discrete`**. Sizes
1...`i_discrete` are integrated per-size and are exact. Therefore:

> **At `i_discrete` = I and `v_discrete` = V, the bin-moment system *is* the
> discrete system.**

So `i_discrete` is not a free choice — it is the **convergence parameter of the
closure**, and refining it walks continuously from the production
configuration to the unapproximated reference *on the same domain, with the
same physics, in the same code path*. This is the textbook grid-refinement
argument the reviewer is asking for, and it removes any "different code,
different answer" objection.

### 1.5 Existing data is a pilot, not a result

The G5 pair above was run at `i_mobile` = 40, `I_bin` = 20, `V_bin` = 29. The
reference run is at `i_mobile` = 5, `I_bin` = 18, `V_bin` = 20. `i_mobile` is
**physics, not numerics** (`run_ensemble.py:1355-1359`: 30 -> 50 moves N_100
+109%, N_111 -18%, d_111 -9%), so **G5 cannot be quoted against the current
reference.** The whole ladder must be re-run at the reference's physics.

What the G5 pilot does establish is that the method works. Comparing the two
arms at matched `at_dose` rungs (identical actual dose, not just identical rung
label):

| dose (dpa) | N_111_vis | d_111 | N_100_vis | d_100 | N_voids | d_cavity |
|---|---|---|---|---|---|---|
| 0.0043 | +0.02% | -0.02% | +0.03% | +0.01% | +0.01% | +0.07% |
| 0.0106 | +0.02% | -0.04% | +0.05% | +0.01% | +0.03% | +0.17% |
| 0.0263 | -0.17% | -0.08% | +3.33% | -0.14% | -0.13% | -0.42% |
| 0.0651 | -5.99% | -1.03% | **+17.37%** | -0.07% | +0.05% | +0.06% |
| 0.1611 | -2.22% | -0.82% | **+20.27%** | -0.49% | -0.72% | -2.81% |
| 0.3990 | **-15.90%** | -2.88% | +14.25% | -0.73% | -0.27% | -1.84% |
| 0.9880 | +3.70% | +0.02% | +9.46% | +1.61% | -0.50% | -3.96% |
| 2.4464 | +1.09% | +0.10% | -0.14% | +0.46% | -0.19% | -1.93% |

Read honestly: the closure is essentially exact below ~0.03 dpa, develops a
**transient excursion of +17-20% in N_100 and -16% in N_111 through the loop
nucleation burst (0.05-0.4 dpa)**, and recovers to a few percent by 2.4 dpa.
That is a real, interesting, publishable result — and it is *stronger* than the
paper's present unqualified claim of "within <~3%", which the data does not
support at all doses.

---

## 2. The study

Three ladders. Every rung is one row of `run_ensemble.py`, scored by the
existing `at_dose` machinery. No new solver code is needed for the numerical
ladder; one small change is needed for the physics ladder (§4).

### 2.1 Ladder N — numerical approximations

Fixed: verification domain, reference physics (`i_mobile` 5, `v_mobile` 5,
fission, G = 1e-7, 330 °C, `loop_conversion` 1, `lnl` 1, `VOID_NETWORK_LOSS` 1).

| rung | knob | values | note |
|---|---|---|---|
| N0 | **reference** | `--equations discrete` | `i_discrete` = I by construction |
| N1 | `i_discrete` | I/2, I/4, I/8, I/16, 100 | the closure convergence parameter (§1.4) |
| N2 | `v_discrete` | V, V/4, 20, 5 | same, vacancy side |
| N3 | `I_bin` / `V_bin` | 10, 18, 20, 29, 40 | bin count at fixed `i_discrete` |
| N4 | `--shape-function` | `constant` (P=1), `linear` (P=2), `lognormal` (P=3) | intra-bin closure |
| N5 | `--solver-mode` | `full_system`, `active_window` | phase-space windowing |
| N6 | `--rtol` | 1e-4, 1e-5, 1e-6, 1e-7 | integrator tolerance |
| N7 | domain | (I, V) doubling pairs | extent verification, `learn.EXTENT_TOL` |

Presented two ways, per the reviewer's request:
- **Cumulative** — start at N0 and switch on one approximation at a time until
  the production configuration is reached. This is literally "approximations
  introduced stepwise".
- **One-at-a-time** — each knob moved alone from the reference. Isolates each
  effect and exposes interactions when the two tables disagree.

### 2.2 Ladder P — physics approximations

| rung | knob | values | note |
|---|---|---|---|
| P0 | **reference** | Case 1 He, dynamic free He | per-size `Q_m`, unapproximated |
| P1 | `he_model` | Case 1 -> Case 2 | per-size `Q_m` -> scalar `Q_tot`. **Requires §4.1** |
| P2 | `he_kinetics` | dynamic -> `quasi_steady_state` | free-He QSS |
| P3 | `loop_conversion` | 1 -> 0 | <111>/<100> partitioning off |
| P4 | `E_a0_conv`, `phi_max_junc` | prior box | partitioning *within* the channel |
| P5 | `i_mobile` | 1, 5, 10, 40 | dominant phase-space knob per §6.2 |
| P6 | `v_mobile` | 1, 5 | |
| P7 | `LOOP_COAL` | on/off | coalescence |
| P8 | `LOOP_NETWORK_LOSS`, `VOID_NETWORK_LOSS` | on/off | loss channels |

### 2.3 Ladder C — conditions (the reviewer's explicit ask)

Ladders N and P are repeated at each condition, so the paper can state whether
the recommended configuration is condition-independent:

| condition | cascade | G_He/G | G (dpa/s) | T |
|---|---|---|---|---|
| C-ref | fission | 0.5-1 | 1e-7 | 330 °C |
| C-rate | fission | 0.5-1 | 1e-6, 1e-5 | 330 °C |
| C-fus | fusion | ~10 | 1e-7, 1e-6 | 330 °C |
| C-temp | fission | 0.5-1 | 1e-7 | 400 °C |

This is what turns a single-point comparison into "careful comparisons across
various test cases under different conditions".

**Transfer to the production domain (§1.3).** For each ladder we additionally
run the *chosen* configuration and its nearest neighbour at the production
domain (I = 80000, V = 20000) and check that the *differences between rungs*
are preserved even though the absolute reference is unavailable there. That is
the honest bridge: we verify the ranking transfers, not the absolute error.

---

## 3. Metrics and the comparison protocol

### 3.1 The six observables

| symbol | field in results row | band (330 °C, 15 dpa) |
|---|---|---|
| N_111 | `N_111_vis` (TEM-visible >= 1 nm) | >= 1.73e21 m^-3 |
| d_111 | `d_111_nm` | 3.4 - 7 nm |
| N_100 | `N_100_vis` | >= 4.67e21 m^-3 |
| d_100 | `d_100_nm` | 3.4 - 7 nm |
| N_void | `N_voids` | >= 3.6e20 m^-3 |
| d_void | `d_cavity_nm` | <= 2.9 nm |

Report **signed % deviation from the ladder's own reference rung**, not band
membership — band membership hides the size of the error, which is the whole
point of a verification study.

### 3.2 Guardrail columns — every table carries these

A verification table that reports only the metric is how Fig. 8 went wrong.
Each row must also carry:

| column | why |
|---|---|
| `dose_reached` / `at_dose` rung | **the anti-Fig-8 column** |
| `delta_FP` | conservation; > 1e-2 means the row is not converged physics |
| `occ_v`, `pile_111` | grid saturation — is the answer the ceiling, not the model? |
| N_eq | the cost/accuracy trade |
| wall_s, threads | reproducible cost |
| `run_cfg_sha` | provenance |

### 3.3 The four rules that answer the reviewer

1. **No run that failed to reach the comparison dose may appear as a curve or a
   metric.** It appears in a "did not reach" column, with the dose it did
   reach. This alone fixes Fig. 8b.
2. **Compare only at matched actual dose**, using `at_dose` and checking the
   rung's `dose` field, not its label. The ladder assigns the last output point
   at or below the target, so rungs 0.4-0.9 can all collapse onto dose 0.399
   (visible in §1.5). A table row labelled "0.9 dpa" that reports 0.399 dpa is
   the same class of error the reviewer caught. Either raise the output-point
   count or print the actual dose in the table.
3. **One axis per figure.** Discretization, closure order, linear solver,
   preconditioner and solver mode never share a legend.
4. **Performance-only knobs get a null test, not a curve.** Preconditioner,
   `--omp-threads` and `linsol` must be shown to agree to solver tolerance; if
   they ever do not, that is a bug report, not a result.

---

## 4. Code changes required

### 4.1 Unweld helium model from cascade spectrum — *required for P1, and the only real blocker*

Today `he_mode` is derived from the spectrum, in both back-ends:

- `py_utils/rate_equations.py:116-128` — `if 'fusion' in po: he_mode='case1' else 'case2'`
- `py_utils/bin_moment_rates.py:676-687` — same
- `cpp_utils/core/parameters.h:423` — **asserts** `he_mode == physics_option % 2`
  and aborts otherwise

So Case 1 is reachable only under fusion and Case 2 only under fission. **The
comparison the reviewer asks for cannot currently be run**: changing the He
model necessarily changes the He/dpa ratio too, and the two effects are
confounded. Any table produced without this fix would be answering a different
question.

The fix is small, because the machinery already exists on both sides —
`rhs_case1` and `rhs_case2` are both implemented in `rate_kernels.cpp` (lines
345, 883), and C++ already carries `he_mode` as its own parameter
(`parameters.h:414`). Required:

1. Add an explicit `he_model` input key, defaulting to the current
   spectrum-derived value so **every existing result is bit-identical**.
2. Relax the `parameters.h:423` assertion to a warning when the override is set
   explicitly; keep it as a hard error otherwise.
3. Thread it through `cpp_bridge.py`'s parameter file and add `--he-model` to
   `run_ensemble.py`.
4. Regression-test: fission + Case 2 and fusion + Case 1 must reproduce the
   current reference to machine precision.

Estimated: ~1 day including the regression test. Everything else in the plan
runs on today's code.

### 4.2 Reporting helpers (small)

- `verification_table.py` — read result `.jsonl`, pick a common `at_dose` rung,
  emit the LaTeX table with guardrail columns. Reuses `learn.py`'s loaders.
- Extend `verify.py`'s pairing logic to the ladder (it already enforces the
  `theta_hash` rule that V12/V13/V15 violated — see its docstring).
- A `--strict-dose` flag on the table generator that **refuses** to emit a row
  whose `dose_reached` is below target, so rule 3.3.1 is mechanical rather than
  a matter of discipline.

---

## 5. Run matrix and cost

Verification domain to be fixed by a calibration run (§6, step 2); the working
assumption is **I = 4000, V = 1000**, where discrete is known to cost 24,471 s
for 6.06 dpa and should reach 15 dpa in roughly 30 h.

| ladder | rows | back-end | est. cost |
|---|---|---|---|
| N0 reference, 4 conditions | 4 | discrete | ~120 h (parallel across machines) |
| N1-N7 | ~35 | bin_moment | ~3 h |
| P0-P8 | ~30 | bin_moment (+ 2 discrete spot-checks) | ~15 h |
| C ladder (N+P at 4 conditions) | ~130 | bin_moment | ~12 h |
| Production-domain transfer | ~12 | bin_moment | ~5 h |
| Null tests (prec / omp / linsol) | ~12 | both | ~2 h |
| **total** | **~225 rows** | | **~160 core-hours** |

Affordable. The measured cost model in `calibration_ledger.json` (`cost_model`)
sizes this across the four registered machines; the discrete references are the
only long poles and they are embarrassingly parallel (one per condition).

---

## 6. Execution order

1. **Write §0 now.** The preconditioner null result and the Fig. 8 diagnosis
   need no new runs — the data is in `compare_linsol_results.json`. Draft the
   replacement subsection and the reviewer response letter first, so the rest
   of the work is confirming a stated position rather than searching for one.
2. **Fix the verification domain.** One discrete run at I = 4000 / V = 1000 at
   reference physics with a 36 h budget. Confirm it reaches 15 dpa and that
   `pile_111` shows the loops are not truncated. If it truncates, drop the
   comparison dose to where I = 4000 is still clean rather than growing the
   domain — cost is O(N^2).
3. **Ladder N.** Pure re-runs on existing code. Produces Table 1 and Table 2.
4. **§4.1 code change**, with regression test.
5. **Ladder P.** Produces Table 3.
6. **Ladder C** and the production-domain transfer. Produces Table 4.
7. **Null tests.** Produces Table 5.
8. **Rewrite §6.2** and replace Fig. 8.

Steps 3 and 4 are independent and can run concurrently.

---

## 7. Tables for the paper

| # | title | rows | columns |
|---|---|---|---|
| 1 | Stepwise introduction of numerical approximations | cumulative N0 -> production | 6 observables (% dev), N_eq, wall, delta_FP, dose |
| 2 | Isolated effect of each numerical knob | one-at-a-time | same |
| 3 | Stepwise introduction of physics approximations | cumulative P0 -> production | same |
| 4 | Condition transferability | chosen config at each (He/dpa, G, T) | 6 observables (% dev), rank preserved? |
| 5 | Null knobs | preconditioner, threads, linsol | agreement in significant figures, wall-clock ratio |

Table 5 is the direct reply to the reviewer and should be small, early, and
unambiguous.

Figures: replace Fig. 8 with **(8a)** cumulative numerical error vs dose for the
six observables, and **(8b)** accuracy-vs-cost — % error against wall-clock,
one point per configuration, the Pareto front labelled. Every curve terminates
where its run actually stopped.

---

## 8. Risks

| risk | mitigation |
|---|---|
| Discrete reference will not reach 15 dpa even at I = 4000 | Compare at the highest dose it does reach and say so; the ladder is still valid at that dose. Do **not** grow the domain — cost is O(N^2). |
| delta_FP > 1e-2 on the reference itself (0.060 on the current production run) | Report it. A verification study measures closure error against the discrete solution of *the same equations*; conservation error is a separate, disclosed defect. Do not conflate them. |
| Transient N_100 excursion (§1.5) is domain- or `i_mobile`-dependent | It is the most interesting result in the study — characterise it across the C ladder rather than averaging it away. |
| §4.1 changes existing numbers | Default preserves current behaviour; regression test asserts bit-identity before any ladder row is run. |
| Reviewer asks for the discrete reference at production scale | Answer with the measured O(N^2) cost and the extent-transfer argument (§2.3). This is a physics-code reality, not an evasion — state the number. |
