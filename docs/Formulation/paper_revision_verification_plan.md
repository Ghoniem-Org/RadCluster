# Plan — Verification & approximation study for *Generalized Cluster Dynamics*

**Purpose:** answer the reviewer's §6.2 / Figure 8 critique with a systematic
study of physics and numerical approximations, scored on the six experimentally
observable parameters.

**Status:** plan. No runs launched, no solver code changed.
**Anchor:** the tracked reference run
`output/20260906_063055_full_system_bin_moment_CD_fission_I80000V20000_im5vm5`.
**Target machine:** MacBook Pro (16 cores, 128 GB) — the machine that produced
the reference run.

> **Revision note (supersedes the first version of this file).** Three things in
> the original were wrong and are corrected here.
> 1. It claimed the Woodbury setup is O(N^2). It is not —
>    `rhs_dispatch.cpp:162-171` uses **Curtis-Powell-Reid coloring**, so a
>    Jacobian rebuild costs `(2*prec_bw+1) + prec_rank` = **33 RHS calls
>    independent of N**. Cost is roughly linear in N.
> 2. It proposed a reduced verification domain of I = 4000 / V = 1000 at
>    15.7 dpa. That domain is **physically inadequate**: the reference's mean
>    cavity is 1071 vacancies at 15.72 dpa, so V = 1000 sits *below the mean*
>    and `d_void` would read the ceiling, not the model.
> 3. It gave N_eq ~ 320 for the reference. The actual state vector is **189**.
>
> The study is consequently restructured: **no fully discrete run at the
> production domain**. Validation (Tables 1-3) and verification (Table 4) are
> separated.

---

## 0. What the reviewer is right about — and the answer already in the repo

Two of the reviewer's claims are correct; the third has a measured answer the
paper never reported.

### 0.1 "Why would a preconditioner change results by an order of magnitude?"

It does not, and we can show it to nine significant figures.

`prec_type` is read at `parameters.h:706` and consumed only by
`prec_setup_*` / `prec_solve_*` in `rhs_dispatch.cpp`. It selects the SPGMR
right-preconditioner (`solver.cpp:405`) and enters no rate kernel, no Jacobian
entry, no state layout. From `codes/Python_Testing/compare_linsol_results.json`,
the one pair where **both arms ran to completion** (bin_moment, N_eq = 190,
t_final = 10000 s, dose 0.01 dpa):

| tag | preconditioner | status | mean_n_i | mean_n_v | swelling % | wall (s) |
|---|---|---|---|---|---|---|
| B-WB | Woodbury | completed | 24.56392587969788 | 5.341081396954862 | 4.562396098e-7 | 185.8 |
| B-JC | Jacobi | completed | 24.56392589863494 | 5.341081573868007 | 4.562396216e-7 | 239.2 |

Nine significant figures on `mean_n_i`, seven on `mean_n_v`, eight on swelling.
**Only the wall-clock differs (1.29x).**

### 0.2 So where did the order of magnitude come from?

From plotting runs that never reached the same dose. `compare_linsol.py:7-9`
says so in its own docstring: on hitting the wall-clock cap "the C++ solver
finalises the current output step gracefully and **the partial trajectory is
kept**." All three discrete arms hit the cap, in different places:

| tag | linsol / prec | status | t_final (s) | dose reached (dpa) | mean_n_i |
|---|---|---|---|---|---|
| D-WB | GMRES + Woodbury | timeout | 1398.7 | 1.40e-3 | 31.18 |
| D-JC | GMRES + Jacobi | timeout | 3962.7 | 3.96e-3 | 31.32 |
| D-KLU | KLU direct | timeout | 0.033 | 3.33e-8 | 8.19 |

D-JC reached 2.83x the dose of D-WB; D-KLU stopped five decades short. Curves
built from these and drawn on one axis differ by orders of magnitude because
they are at different doses. Presentation defect, not a solver defect. (D-WB and
D-JC agree on `mean_n_i` to 0.4% *despite* the 2.8x dose gap — further evidence
the preconditioner is inert.)

### 0.3 Secondary defects in §6.2 to fix while we are there

- **Caption/text dose mismatch.** The Fig. 8 caption says "at 0.1 dpa" for both
  panels; the text (p. 50) describes a 10 dpa horizon and a 9198 s discrete
  reference. Neither matches `compare_linsol_results.json`, whose `t_span` is
  (1e-6, 1e4) s = **0.01 dpa**. Two studies are being presented as one figure.
- **Conflated axes.** One legend mixes discretization, closure order,
  linear solver, preconditioner and solver mode. These live at different levels.
- **Unreported non-convergence.** `delta_FP` = 0.109 for B-WB/B-JC, against the
  model's own 1e-2 gate. Those arms agree with each other but are not converged
  physics. Preconditioner agreement proves the linear-algebra null; it does not
  license the run as accurate.
- **Metric mismatch.** §6.2 scores `N_loops`, `<n_i>`, `<m_v>`, `C_He`; §6.1
  scores densities and sizes against TEM. Use the six observables throughout.

**Deliverable from §0 — no new runs required.** Table 5 below plus a rewritten
subsection. Draft this first, so the rest of the work confirms a stated position
rather than searching for one.

---

## 1. Structure of the study, and why it is split

### 1.1 The constraint that forces the split

The reference domain (I = 80000, V = 20000) is not a luxury:

- The reference's mean cavity at 15.72 dpa is **1071 vacancies**
  (`occ_v` = 1071/20000 = 0.054, matching the campaign report).
- The ledger already records `pile_111` = 1.000 at **I = 20000** — "I=20000
  truncates the loops at 15 dpa."

So any domain small enough to make a fully discrete run cheap is too small to
be compared with experiment at 15.72 dpa. A discrete run **at** the production
domain is N_eq ~ 180 000 and, extrapolating the measured G5 point, costs
somewhere between one week and several months — an unattended risk with no
proportionate return.

### 1.2 The split

| | domain | dose | reference | claim |
|---|---|---|---|---|
| **Tables 1-3** — validation | I = 80000, V = 20000 | 15.72 dpa | most-resolved bin-moment rung | model reproduces experiment; closure has converged in its own refinement parameter |
| **Table 4** — verification | I = 8000, V = 20000 | 0.2 / 1 / 2 dpa | **fully discrete** | closure reproduces the exact solution of the same equations |

State the residual gap honestly in the paper: Tables 1-3 demonstrate
*self-convergence*, which in principle could converge to a wrong limit; Table 4
is what excludes that, by measuring the closure against the exact answer in a
regime where the exact answer is computable. Suggested wording:

> *The closure is verified against the exact solution of the same equations
> where that solution is computable (Table 4). At production conditions, where
> it is not, convergence is demonstrated in the closure's own refinement
> parameter (Table 1).*

### 1.3 The design rule that ties the two together

The bin-moment grid contains two independent things: **how much of the
distribution is resolved exactly** (`i_discrete`) and **how coarse the bins are**
(`r = (I/i_discrete)^(1/I_bin)`, derived — see `run_ensemble.bin_layout`). Vary
both at once and the ladder is uninterpretable.

**Hold `r` at the production value and vary only `i_discrete`.** Production is
`r_i` = 1.4497 and `r_v` = 1.5138. This fixes every bin count in the study.

It also makes Table 4 transfer. The production closure resolution is
`i_discrete/I` = 100/80000 = **1/800** and `v_discrete/V` = 5/20000 = **1/4000**.
At I = 8000 the matching SIA rung is `i_discrete` = 10 — and holding `r_i` fixed
returns `I_bin` = 18, *exactly the production value*. Keeping V = 20000 in
Table 4 makes the vacancy side match production identically (`v_discrete` = 5,
`V_bin` = 20). So Table 4's last column is an exact analogue of the production
closure on both axes, and the closure error it measures is the closure error the
production run carries.

### 1.4 Reference run facts (from its `provenance.md`)

| | |
|---|---|
| equations | `bin_moment`, `shape_function` linear (P = 2) |
| grid | `i_discrete` 100, `I_bin` 18, `v_discrete` 5, `V_bin` 20 |
| domain | I = 80000, V = 20000 |
| N_eq | **189** (⟨100⟩ block carried separately) |
| mobility | `i_mobile` 5, `v_mobile` 5 |
| solver | `full_system`, GMRES + Woodbury, rtol 1e-5, atol 1e-20 |
| physics | fission, G = 1e-7 dpa/s, 330 °C, `loop_conversion` 1, `lnl` 1, `VOID_NETWORK_LOSS` 1, `void_net_chi` 1165, `f_cl_v` 0.65 |
| output grid | 37 log-spaced points, `t_span` (1e-6, 4e8) s; **15.72 dpa is point 35** |
| cost | 1317 s at 12 OpenMP threads, to 40 dpa |
| observables at 15.72 dpa | N_100 4.715e21, d_100 5.939, N_111 3.498e21, d_111 4.637, N_void 3.609e20, d_void 2.890 |
| conservation | `delta_FP` 0.0604 at 15.72 dpa, 0.0696 at 40 dpa |

---

## 2. The tables

### 2.1 Table 1 — numerical approximations (closure convergence)

All columns: I = 80000, V = 20000, `full_system`, linear closure, fission,
G = 1e-7, 330 °C, scored at 15.72 dpa. Deviations quoted **relative to B1**, the
most-resolved rung.

| rung | `i_discrete` | `I_bin` | `v_discrete` | `V_bin` | N_eq (main block) |
|---|---|---|---|---|---|
| B1 | 6400 | 7 | 1600 | 6 | ~8 030 |
| B2 | 1600 | 11 | 400 | 9 | ~2 050 |
| B3 | 400 | 14 | 100 | 13 | ~560 |
| B4 | 100 | 18 | 5 | 20 | **189** (= production) |

| | Experiment | B1 | B2 | B3 | B4 = production |
|---|---|---|---|---|---|
| N_111 (10^21 m^-3) | >= 1.73 | *(ref)* | | | 3.498 |
| d_111 (nm) | 3.4 - 7 | *(ref)* | | | 4.637 |
| N_100 (10^21 m^-3) | >= 4.67 | *(ref)* | | | 4.715 |
| d_100 (nm) | 3.4 - 7 | *(ref)* | | | 5.939 |
| N_void (10^20 m^-3) | >= 3.6 | *(ref)* | | | 3.609 |
| d_void (nm) | <= 2.9 | *(ref)* | | | 2.890 |
| *deviation from B1* | | *(ref)* | *%* | *%* | *%* |
| N_eq | — | ~8 030 | ~2 050 | ~560 | 189 |
| wall (s) | — | | | | 1317 |
| dose reached (dpa) | — | | | | 40.0 |
| delta_FP | — | | | | 0.060 |
| occ_v, pile_111 | — | | | | 0.054, — |

The lower block is not decoration. `dose reached` is the row that makes
Figure 8's failure mode impossible. `occ_v` and `pile_111` are the quantitative
demonstration that the domain is adequate — the claim that killed the original
reduced-domain design.

### 2.2 Table 2 — intra-bin closure and tolerance, one at a time

Each column changes exactly one knob from B4; deviations relative to B4, since
these sit *inside* the closure rather than on the discrete/binned axis.

| | B4 | P=1 constant | P=3 lognormal | `I_bin` 10 | `I_bin` 40 | rtol 1e-4 | rtol 1e-6 |
|---|---|---|---|---|---|---|---|
| six observables | as above | | | | | | |
| wall (s) | 1317 | | | | | | |
| dose reached (dpa) | 40.0 | | | | | | |
| delta_FP | 0.060 | | | | | | |

This is where the existing §6.2 claims about P = 1 and P = 3 are either
supported or withdrawn. They are kept out of Table 1 deliberately: shape
function is not on the discrete-vs-binned axis, and mixing the two is what made
the old figure unreadable.

### 2.3 Table 3 — helium

No experimental column — the fusion arms have no matching EUROFER dataset, so
this is model-to-model. All at the B4 grid, 15.72 dpa.

| | fission QSS (= B4) | fission dynamic | fusion QSS | fusion dynamic | fusion Case 1 | fusion Case 2 |
|---|---|---|---|---|---|---|
| six observables | as above | | | | | |
| C_He,tot (m^-3) | 1.18e13 | | | | | |
| wall (s) | 1317 | | | | | |
| delta_He | 9.4e-11 | | | | | |

Three requirements for this table.

1. **`C_He,tot` must be a row.** It is the quantity the He approximations act
   on; without it a flat table reads as a null result rather than a
   demonstration.
2. **The fusion columns are not optional.** Under fission the reference traps
   1.18e13 m^-3 and §6.2 already reports a ~2000x fission/fusion separation, so
   a He block evaluated only at fission will be flat across all six observables
   and will read as evasion. The flatness is a *result* — the QSS reduction is
   safe at fission precisely because there is no helium — but only when the
   fusion columns sit beside it.
3. **The last two columns need the §3.1 code change.** Columns 1-4 are runnable
   on today's code.

### 2.4 Table 4 — binning verification against the exact solution

**I = 8000, V = 20000**, fission, G = 1e-7, 330 °C, `full_system`, reference
physics. SIA domain reduced (that is where the cost is — ⟨111⟩ and ⟨100⟩ are two
blocks); vacancy domain kept at production so the vacancy closure is tested at
the production fraction. Deviations relative to the discrete column.

| rung | `i_discrete` | `I_bin` | `v_discrete` | `V_bin` | N_eq |
|---|---|---|---|---|---|
| D (exact) | 8000 | — | 20000 | — | ~36 000 |
| C1 | 800 | 6 | 1600 | 6 | ~2 430 |
| C2 | 200 | 10 | 400 | 9 | ~640 |
| C3 | 50 | 14 | 100 | 13 | ~210 |
| **C4** | **10** | **18** | **5** | **20** | **~100** |

C4 carries `i_discrete/I` = 1/800 and `v_discrete/V` = 1/4000 — identical to the
production run on both axes (§1.3). It is the column the transfer argument
rests on.

| | D (exact) | C1 | C2 | C3 | C4 (= production resolution) |
|---|---|---|---|---|---|
| N_111, d_111, N_100, d_100, N_void, d_void | *(ref)* | | | | |
| *deviation from D* | *(ref)* | *%* | *%* | *%* | *%* |
| N_eq | ~36 000 | ~2 430 | ~640 | ~210 | ~100 |
| wall (s) | | | | | |
| dose reached (dpa) | | | | | |
| delta_FP, occ_v, pile_111 | | | | | |

Two additions specific to this table.

- **Report it at three doses (0.2, 1, 2 dpa), not one.** The `at_dose` ladder
  gives them from the same trajectories at no extra cost, and the closure error
  is strongly dose-dependent: the G5 pilot showed **+17 to +20% in N_100 and
  -16% in N_111 through the loop nucleation burst at 0.05-0.4 dpa**, recovering
  to a few percent by 2.4 dpa. Collapsing that to a single dose reports the
  recovered value and hides the excursion — the one thing a hostile reader could
  later say was concealed. A small companion figure of deviation vs dose costs
  nothing.
- **Add the P = 1 and P = 3 columns here too.** Here they can be scored against
  the exact answer rather than against B4.

**Known limitation to state:** Table 4 runs at 2 dpa, so it verifies the closure
in the nucleation-and-early-growth regime, not at 15.72 dpa where coalescence
and the sweeping channels dominate. That gap is the price of eliminating the
production-domain discrete run, and it should be declared rather than glossed.

### 2.5 Table 5 — null knobs (the direct reply to the reviewer)

| | Woodbury | Jacobi | agreement | 1 thread | 12 threads | agreement |
|---|---|---|---|---|---|---|
| six observables | | | *sig. figs* | | | *sig. figs* |
| wall (s) | | | *ratio* | | | *ratio* |

Small, early, unambiguous. Seed it with the B-WB/B-JC pair from §0.1 and extend
to the B4 configuration.

---

## 3. Protocol

### 3.1 Required code change — unweld helium model from cascade spectrum

Needed only for Table 3's last two columns. Today `he_mode` is derived from the
spectrum in all three places:

- `py_utils/rate_equations.py:116-128` — `if 'fusion' in po: case1 else case2`
- `py_utils/bin_moment_rates.py:676-687` — same
- `cpp_utils/core/parameters.h:423` — **asserts** `he_mode == physics_option % 2`
  and aborts otherwise

So Case 1 is reachable only under fusion, Case 2 only under fission: changing
the He grouping necessarily changes the He/dpa ratio, and the two effects are
confounded. Both RHS paths already exist (`rate_kernels.cpp:345` and `:883`) and
C++ already carries `he_mode` as its own parameter (`parameters.h:414`), so:

1. Add an explicit `he_model` input key defaulting to the spectrum-derived value
   so **every existing result is bit-identical**.
2. Relax the `parameters.h:423` assertion to a warning when the override is set
   explicitly; keep it a hard error otherwise.
3. Thread it through `cpp_bridge.py` and add `--he-model` to `run_ensemble.py`.
4. Regression test: fission+Case 2 and fusion+Case 1 reproduce the current
   reference to machine precision.

~1 day including the test. Nothing else in the plan needs solver changes.

### 3.2 The output grid must be identical across every column

`15.72 dpa` is output point 35 of the reference's 37-point log grid over
`t_span` (1e-6, 4e8) s. **Every Table 1-3 run must use the same `t_span` and
`n_points`**, or "15.72 dpa" is not the same point in different columns and the
comparison silently drifts — the same class of error as Figure 8.

Practical consequence: run to **40 dpa**, not 15, and read point 35.
`conditions_S8.json` carries `dose_dpa: 15.0`, which *overrides* `--dose`
(`cond.get('dose_dpa', cfg['dose'])`), so a new conditions file with
`dose_dpa: 40.0` is required.

### 3.3 Report the realised bin layout, not the requested one

`run_ensemble.bin_layout()` warns that the realised bin count is **not** the
requested one — `r` is derived and `build_bins()` walks `floor(edge*r)`. Every
table must print the realised `I_bin` / `V_bin` / `r` read back off
`rate_equations`, not the CLI values.

### 3.4 The four rules that answer the reviewer

1. **No run that failed to reach the comparison dose may appear as a curve or a
   metric.** It appears in a "did not reach" column with the dose it did reach.
2. **Compare only at matched actual dose.** Use `at_dose` and check the rung's
   `dose` field, not its label — the ladder assigns the last output point at or
   below the target, so several rung labels can collapse onto one output point.
3. **One axis per figure.** Discretization, closure order, linear solver,
   preconditioner and solver mode never share a legend.
4. **Performance-only knobs get a null test, not a curve.** If preconditioner,
   thread count or `linsol` ever move a converged answer, that is a bug report.

### 3.5 Reporting helper

`verification_table.py`: read result `.jsonl`, select a common `at_dose` rung,
emit the LaTeX table with guardrail rows, reusing `learn.py`'s loaders. Give it
a `--strict-dose` flag that **refuses** to emit a row whose `dose_reached` is
below target, so rule 3.4.1 is mechanical rather than a matter of discipline.

---

## 4. Runs — MacBook Pro

16 cores, 128 GB. Single row at a time with `--workers 1 --omp-threads 12`,
matching the reference run's configuration. The reference cost 1317 s that way.

Template (Table 1, rung B3):

```bash
cd RadCluster_2_1/digital_twin
python run_ensemble.py \
    --design design/T1_B3.csv \
    --conditions conditions_S8_40dpa.json \
    --spec parameters_S4.json \
    --out results/T1_B3_machine0.jsonl \
    --machine 0 --of 1 \
    --equations bin_moment \
    --i-discrete 400 --i-bin 14 --v-discrete 100 --v-bin 13 \
    --shape-function linear \
    --I 80000 --V 20000 --dose 40.0 --lnl 1 --rtol 1e-5 \
    --solver-mode full_system \
    --timeout-s 86400 --workers 1 --omp-threads 12
```

Table 4's discrete rung differs only in `--equations discrete --I 8000
--V 20000 --dose 2.0` (with the bin flags dropped).

| block | runs | estimate | notes |
|---|---|---|---|
| Table 1 (B1-B3; B4 exists) | 3 | ~10 h | B1 at N_eq ~8 000 dominates |
| Table 2 | 6 | ~2 h | P=1 may be unstable, P=3 may time out — both are results |
| Table 3 | 5 | ~3 h | last 2 columns blocked on §3.1 |
| Table 4 | 6-7 | ~4-8 h | discrete arm ~36 000 eq at 2 dpa dominates |
| Table 5 | 4 | ~2 h | doubles as the §0 evidence |
| **total** | **~25** | **~1 day of compute** | fully parallel across machines if wanted |

The Table 4 discrete arm is the only run whose cost is genuinely uncertain: it
is ~4x G5's N_eq at a third of G5's dose, with 12 threads where G5 had 1. G5
measured 24 471 s single-threaded to 6.06 dpa at N_eq ~9 000. **Run it first**
— if it lands far from the estimate, the rest of the schedule is wrong too.

---

## 5. Execution order

1. **Write §0 and Table 5.** No new runs — the data is in
   `compare_linsol_results.json`. Draft the replacement subsection and the
   reviewer response letter first.
2. **Table 4 discrete arm.** The cost probe and the study's only real unknown.
3. **Table 4 remaining rungs** (cheap) — gives the verification result early.
4. **Table 1** B1-B3.
5. **§3.1 code change** + regression test (independent of 3-4; can overlap).
6. **Tables 2 and 3.**
7. **Rewrite §6.2**, replace Fig. 8 with: (8a) closure deviation vs dose for the
   six observables from Table 4; (8b) accuracy vs wall-clock, one point per
   configuration, Pareto front labelled. Every curve terminates where its run
   actually stopped.

About four days end to end, of which roughly one is compute.

---

## 6. Risks

| risk | mitigation |
|---|---|
| Reviewer objects that Tables 1-3 lack a true reference | This is the declared V&V split (§1.2). Table 4 measures the closure error at *exactly* the production resolution (§1.3), and the domain constraint that forces the split is quantitative: mean cavity 1071 vacancies, `pile_111` = 1.000 at I = 20000. |
| Table 4's 2 dpa does not cover the 15.72 dpa regime | Declared as a limitation (§2.4). Report Table 4 at three doses so the dose trend is visible and extrapolatable. |
| Table 4 discrete arm costs far more than estimated | It runs first (§5.2). If it overruns, cut I to 4000 for Table 4 only — at 2 dpa that is defensible on occupancy grounds even though it is not at 15.7 dpa. |
| `delta_FP` ~ 0.06 on the reference itself | Report it. A verification study measures closure error against the same equations; conservation error is a separate, disclosed defect. Do not conflate them. |
| P = 1 / P = 3 fail to run | That is the result §6.2 already claims. Record the failure mode and the dose reached; do not plot a partial trajectory. |
| §3.1 changes existing numbers | Default preserves current behaviour; regression test asserts bit-identity before any ladder row runs. |
