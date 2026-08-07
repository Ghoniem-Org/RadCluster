# Digital-twin campaign harness — 4-machine protocol

> **Start here: `campaign_control_rev6.ipynb` — the ONLY driver.** Open it on
> each machine. It
> configures the run, auto-builds the solver, verifies this machine agrees with
> the reference, launches the worker, gives a live monitor against the
> experimental data, and provides graceful stop/restart. The command lines
> below are what it drives — use them directly if you prefer a shell.
>
> **Keep it the only driver.** A second campaign notebook has twice appeared
> beside it under a different name, and both times the copies drifted: the first
> (2026-08-06, commit `11dc308`) still had `MACHINE = 2` hard-coded and pointed
> at the superseded T2 design, so opening the obvious file would have launched
> the wrong campaign with two machines claiming index 2; the second
> (2026-08-07) lacked the per-machine `timeout_s` override and the Windows
> detachment fix, so it would have relaunched Matrix-PC straight back into
> 100 % starvation. Neither copy announced it was stale — they carried the same
> title and both auto-detected. If you need a variant, branch; do not add a file.

Three modules plus a validation test. Everything else in the plan's §3 listing
(`surrogate.py`, `calibrate.py`, …) is downstream of these and not yet built.

```
digital_twin/
├── parameters.json      # H.11-1 spec: ranges, priors, tier, fixed values
├── conditions.json      # the three Tier-2 conditions (N2, N5, I1)
├── design.py            # generate the ONE canonical Saltelli design
├── run_ensemble.py      # the per-machine worker
├── merge_and_sobol.py   # concatenate + S_i, S_i^T with pairwise deletion
├── test_sobol.py        # Ishigami validation of the estimator
├── design/              # committed design files (+ .meta.json sidecars)
└── results/             # one .jsonl per machine
```

## Bringing a new machine online

Git carries the design, the workbook and the code — but **not the solver**:
`build/` is ignored, so each machine compiles its own from `cpp_utils/`. Four
machines with different SUNDIALS or compiler versions produce numerically
different results while every provenance hash still looks plausible, and the
merge step would pool them into one Sobol estimate.

So on each machine, in order:

```bash
git pull
cd RadCluster_2_1 && cmake -S . -B build && cmake --build build --config Release
cd digital_twin && python check_machine.py
```

`check_machine.py` runs a fixed 30 s probe and compares 12 quantities — two
diffusivities, both success gates, four population readouts, two inventories,
swelling and `delta_FP` — against `machine_reference.json` at `rtol = 1e-9`.
CVODE is deterministic for a fixed binary (verified: 0.00e+00 on a repeat run),
so a nonzero diff means a different build, not floating-point noise.

**Do not start a machine that fails this check** — its rows cannot be pooled
with the others. Rebuild from the same commit and re-check.

## Running it

**Once, on one machine — then commit `design/`:**

```bash
python design.py --tier 2 --N 16 --out design/T2_design_v1.csv
```

`1104 rows = 16 × (21+2) × 3 conditions`, ~276 per machine.

**On each machine k = 0,1,2,3:**

```bash
python run_ensemble.py --design design/T2_design_v1.csv \
                       --machine k --of 4 --workers 10 \
                       --I 800 --V 600 --dose 0.1
```

Machine `k` runs every row with `row_id % 4 == k`. No coordination, no shared
filesystem, restartable — re-running skips rows already in its `.jsonl`.

## Tallying across machines

Each machine writes **its own** file — `T2_design_v1_machine0.jsonl` …
`machine3.jsonl` — so the four never touch the same bytes and **can never
merge-conflict**. That makes git the transport:

```bash
# on each machine, when it finishes (or at any checkpoint)
git add digital_twin/results/T2_design_v1_machine$K.jsonl
git commit -m "T2 results, machine $K"
git push
```

```bash
# anywhere, to tally
git pull
python merge_and_sobol.py --design design/T2_design_v1.csv --results results/
```

Merging is keyed on `row_id`, so it is order-independent and idempotent — pull
again later and re-run to fold in whichever machines have since finished. Rows
are appended as they complete, so you can tally a partial campaign at any time;
`merge_and_sobol.py` reports coverage and refuses to estimate an index from too
few usable base points.

A shared network drive works equally well — point every machine's `--out` at it
and skip the commits. Manual copy is the fallback; nothing depends on the files
arriving together.

**Stopping a distributed campaign — the stop flag does NOT travel.**
`campaign_ops.request_stop()` writes a local `CAMPAIGN_STOP` sentinel that only
the workers on *that* machine poll. It is not in git, and a running worker never
pulls. To halt all four machines you must, on **each** of them:

```python
import campaign_ops as ops; ops.request_stop('reason here')
```

or kill the `run_ensemble.py` PID directly. Either is safe: rows are appended as
they complete, so nothing already written is lost and a restart resumes by
skipping `row_id`s already present. This had to be done by hand on 2026-08-03
and is worth automating if the campaign ever exceeds four machines.

**What the tally checks for you.** It prints a `PROVENANCE SPLIT` warning if the
machines disagree on `git_sha`, `solver_sha256`, `workbook_sha256` or
`design_sha256`, and it reports missing `row_id`s bucketed by `row_id % 4` — a
whole residue class missing means one machine never reported, which is
immediately visible rather than silently shrinking the sample.

## Why the design is generated once and committed

Two independent reasons, both fatal if ignored:

1. **Saltelli indices are computed from paired rows.** `S_j` pairs `f(A_i)`
   against `f(AB_i^(j))`. If two machines sample independently, the pairing is
   destroyed and the indices are silently *wrong* — not noisy, wrong.
   `run_ensemble.py` refuses to start if the design file's SHA-256 does not
   match its sidecar.

2. **Tier 3 is a nested subset of Tier 2.** The shared `θ` points between LF and
   HF *are* the multi-fidelity discrepancy model, so they must coincide exactly.
   A Sobol sequence is nested by construction; Latin hypercube is not and cannot
   be extended without rebuilding. `θ` is referenced by stable `row_id`.

## The admissibility block — read this before trusting any index

`δ_FP` is **blind to grid truncation**. On 2026-08-01 it held at `1e-8` while
99.96 % of the ⟨100⟩ content was stacked against the top of the grid, because
the `if (n < I)` growth guard in `rate_kernels.cpp` halts growth without losing
atoms. Conservation is a health check, not a grid-adequacy check.

Across ~1100 rows with `θ` spanning its full prior box, some cells **will**
saturate the grid. Every row therefore carries:

| field | meaning |
|---|---|
| `occ_111`, `occ_100` | mean size / `I` — **recorded only, no longer a reject rule** (see below) |
| `pile_111`, `pile_100` | content fraction in the top 2 % of the grid |
| `d_over_ceiling_100` | `d_100 / d(n=I)` |
| `dose_reached`, `starved` | did it reach the requested dose |
| `delta_FP`, `delta_He` | conservation |
| `topbin_111`, `topbin_100`, `topbin_v` | `bin_moment`: content fraction in the **top bin**, from raw μ₁ |
| `grid_limited`, `grid_converged` | did any truncation measure trip |
| `admissible` | **since 2026-08-05: the row ran and reached the dose. Nothing else.** |

> ## ⚠ Truncation no longer gates anything — author's decision, 2026-08-05, plan §12(q)
>
> **The campaign estimates the relative RANKING of parameters, not the absolute
> value of any observable, and a ranking survives a truncated tail.** Every field
> in the table above is still measured and written on every row; none of them may
> reject a row. `admissible` now means `not starved`. The grid is fixed at
> **`I = 30000`, `V = 5000`** for all T2/T3 runs and no further convergence study
> will be run — plan §11(i)-1/2/3 are **withdrawn, not deferred**.
>
> **`δ_FP` was withdrawn with them.** Kept as a gate at `1e-2` it would have gone
> on rejecting truncated rows under a different name — **341 of the 459 pooled
> rows (74 %) sit at `δ_FP ≥ 1e-2`**.
>
> **`δ_FP` is not a proxy for truncation on either axis** — measured 2026-08-05,
> and this corrects an earlier claim in this file that it tracked ⟨100⟩
> truncation:
>
> | run | `topbin_100` | `topbin_v` | `δ_FP` |
> |---|---|---|---|
> | row 24, I=10000, V=10000 | 0.998 | 3.6e-24 | **0.907** |
> | row 825, I=30000, V=5000 | 0.975 | 4.3e-24 | **9.85e-04** |
> | row 229, I=30000, V=5000 | 0.995 | **0.565** | **3.95e-01** |
>
> Rows 24 and 825 are both near-totally ⟨100⟩-truncated with no vacancy
> truncation, and their `δ_FP` differs by three orders of magnitude. What moves
> it here is the **vacancy** axis — consistent with plan §11(j), "`δ_FP` is blind
> to ⟨100⟩ truncation". Its response is θ- and condition-dependent, which is a
> reason to record it, not to gate on it.
>
> **Re-scoring the 516 pooled v2 rows, at zero CPU:** 12 admissible as shipped →
> 74 under the corrected rules → **459** with truncation ungated.
>
> **Before publishing any ranking, run the report both ways:**
> ```bash
> python merge_and_sobol.py --design design/T3_rev6.csv --results results/
> python merge_and_sobol.py --design design/T3_rev6.csv --results results/ --require-converged
> ```
> `--require-converged` restores the pre-2026-08-05 rule exactly. If the ranking
> is unchanged, truncation did not buy it anything — that is the empirical form
> of the claim this decision rests on. The default report prints what fraction of
> the rows it used are truncated, so a ranking is never read without that number.
>
> ## ⚠ A timed-out row is KEPT — author's decision, 2026-08-05, plan §12(s)
>
> `--timeout-s` does not kill the solver: `cpp_bridge` sends a graceful
> interrupt and the solver **flushes its trajectory**, so a row that runs out of
> time returns everything up to the dose it reached. Discarding those was
> throwing away finished work — 57 of the 516 pooled v2 rows died for this alone.
>
> **But such a row must not be read at whatever dose it stopped at.** Dose moves
> the observables far harder than any numerical choice we have measured — for
> row 24, `d_100` goes 5.36 nm → 26.505 nm from 0.1 to 1.0 dpa (**5×**), against
> the **0.03 %** the whole `I=50000 → 200000` grid refinement moved it. And a row
> is slow *because of its physics*, so `dose_reached` is correlated with θ:
> pooling end-of-run values across doses would attribute "how far this run got"
> to the parameters. That is a worse failure than the truncation above, because
> it is aligned with the estimand.
>
> **So every row records a dose ladder** (`at_dose`, rungs 0.1 … 1.0 dpa) —
> free, because the `n_points = 40` trajectory is already computed. Screen at a
> common dose:
> ```bash
> python merge_and_sobol.py --design design/T3_rev6.csv --results results/ --at-dose 0.8
> ```
> A timed-out row then contributes to every rung it reached and is absent above
> that — **pairwise, exactly like a missing `AB` row**, which the estimator
> already handles and reports through `n_eff`. The report prints the ladder
> coverage before estimating anything, so you can see which rung is worth using.
>
> Without `--at-dose` a starved row is still excluded rather than silently
> pooled at the wrong dose, and the report names the rung to use instead.
>
> **The timeout is now a budget knob, not a cliff.** Below it a row used to be
> worth 100 % and above it 0 %; now it degrades gracefully. Choose it from the
> cost distribution so the common rung lands high, rather than chasing a p99.
>
> **The residual risk, stated plainly:** truncation depth is strongly θ-dependent
> (at `I = 3200`, row 825 sits at `pile_100 = 1.3e-09` and row 229 at `0.616` —
> same grid, same dose), so it is not a uniform distortion. Where it saturates an
> observable it compresses variance for large-loop θ. The ⟨100⟩ observables are
> **already withheld** under `bin_moment` for an unrelated closure bias
> (`BIN_MOMENT_BLOCKED`), leaving `N_loops_111`, `d_111_nm`, `N_voids`,
> `d_cavity_nm`.
>
> **Only half of those are truncation-safe** — corrected 2026-08-05, plan §12(t2).
> The ½⟨111⟩ pair is: `topbin_111 ≤ 1.2e-02` on every row measured, because that
> block is populated across its range rather than peaked and advecting. **The
> void pair is not.** Row 229 reaches `mean_n_v = 2847` against `V = 5000` —
> 1.8× headroom, `topbin_v = 0.565`, `d_cavity = 4.00 nm` against a 4.83 nm
> ceiling — and `mean_n_v` has never converged on the vacancy axis
> (276 → 709 → 2847 at `V` = 600 → 2400 → 5000). So `N_voids` and `d_cavity_nm`
> will be ranked partly from rows sitting against the vacancy ceiling.
> `--require-converged` matters most for exactly those two.

A row that is **inadmissible is not a failed row** — it ran fine and conserves,
but its observables measure the grid rather than the physics. The smoke test at
`I = 150` produced exactly this: `pile = 0.996`, correctly rejected.

> **Superseded by §12(q) above — retained because it explains the fields.**
> **Revised 2026-08-03 after the T2-v2 stop — plan §11(c).** Two of these gates
> were rejecting good rows. Of 275 completed v2 rows, 100 were genuinely
> grid-converged and only 12 were scored admissible.
>
> - **`occ > 0.10` is withdrawn as a reject criterion.** Row 24 (`occ_100 = 0.171`,
>   `pile_100 = 2.3e-10`) returns `d100 = 5.3573` and `N100 = 3.99268e21` at both
>   `I = 3200` and `I = 12800` — identical. Occupancy is not a truncation test.
>   The same argument was already accepted for `occ_v` in plan §10(j)-1.
> - **`δ_FP < 1e-6` equals the solver's own `rtol`** and was rejecting for solver
>   noise; 88 % of grid-converged rows failed it. Moving to `1e-2` for the
>   bin-moment campaign — plan §11(h).
> - **`pile_v ≤ 0.05` is too loose to guard the vacancy axis.** At `I = 12800`,
>   `V = 600`, row 229 *passed* `pile_v` at 0.026 while `mean_n_v` was 61 % low
>   and `N_100` 20 % off; `δ_FP` caught it at 0.530, dropping to 7.1e-6 at
>   `V = 2400`. **`δ_FP` is the working vacancy guard, not `pile_v`.**
>
> Apply admissibility **post-hoc in `merge_and_sobol.py`**, not in the worker, so
> a threshold can be revised without re-running anything.

## Failure handling — pairwise, not global

The plan says "discard and flag". For Saltelli that is incomplete: `A_i` and
`B_i` are shared across all `p` parameters, so deleting a row globally biases
*every* index. `merge_and_sobol.py` therefore uses **pairwise deletion** —
if `AB_i^(j)` is unusable, base index `i` is dropped from parameter `j`'s
estimator only. Each index reports its own `n_eff`.

**Read `n_eff` before reading any index.** An `S_i^T` computed from 3 of 16 base
points is not a screening result. Validated in `test_sobol.py`: dropping 30 % of
one parameter's AB rows leaves the other parameters' indices bit-identical.

## Provenance

Every row records `git_sha`, `solver_sha256`, `workbook_sha256`,
`design_sha256`, `machine_id`. `merge_and_sobol.py` prints a **PROVENANCE
SPLIT** warning if the four machines disagree on any of them — four machines
with different solver builds produce silently incomparable results, which is
the same class of error as the grid saturation above.

Set `OMP_NUM_THREADS=1` (the runner does this) and use many single-threaded
workers rather than few threaded ones: better throughput for an ensemble, and it
removes OpenMP reduction-order nondeterminism between machines.

## Gotchas that have already cost real work

**`RadClusterSimulation.__init__` silently discards unknown kwargs.** Its
signature ends in `**legacy_kw`, and anything not named there is swallowed and
dropped — no error, no warning. Passing `i_discrete=…`, `I_bin=…` or
`shape_function=…` to the constructor does **nothing**. On 2026-08-02 this
produced nine "refinement" runs that came back bit-identical and were nearly
read as evidence of a structural closure failure. Set them through
`input_data.reactions[...]` then `_calculate_derived()` + `rebuild_rates()`.

**The bin configuration lives on `sim.rate_equations`, not
`sim.reaction_rates`.** `reaction_rates` is a `ReactionRates` in both modes;
`rate_equations` is the `BinMomentRateEquations` that owns `i_discrete`,
`I_bin`, `r`, `shape_function`, `n_mom`, `N_eq`. Querying the wrong object
returns `?`/`nan` rather than raising.

**The `bin_moment` default binning is very coarse:** `i_discrete=10, I_bin=6,
r=1.849` — six bins spanning sizes 11…I. Any accuracy comparison that does not
state its binning is comparing against that.

**`δ_FP` cannot see grid truncation** (see above) — judge grid adequacy by
`pile` and `d/d_ceiling`.

## Before a production launch — open items

1. **Two parameters are `REVISION_PENDING`** (`design.py` warns): `E_a0_conv`
   and `dH2_abs_conv`. Their ranges predate the 2026-08-01 absorption-gate work.
   `E_a0_conv`'s rev-4 range 1.45–1.85 was set with ⟨100⟩ absorption throttled;
   with it open, on-target `N_100` needs ≈2.2–2.25.
2. **`E_a0_conv` may not need sampling at all.** `N_100 ∝ exp(−E_a0/k_BT)`
   exactly (41.3× per 0.2 eV, matching `1/k_BT` to three digits), so it is
   analytically invertible from a single run — a candidate for demotion to a
   derived quantity, which removes the strongest degenerate direction from the
   design.
3. **`ΔH₂^abs` and `A_abs` are degenerate at fixed `T`** — only their product
   (the gate) is identifiable from one temperature. `psucc_abs_pref` is
   therefore held fixed at 2.0 and only `dH2_abs_conv` is sampled; the two
   separate solely through the `f_100(T)` spread across the three conditions.
4. **T0.2 (LF grid) must be re-run with the absorption channel open.** The
   plan's ~40 s/run cost model comes from an anchor scan taken with it
   throttled. Accept a grid only if `pile ≈ 0` at the largest-`d_100` corner of
   the prior box.
5. **T0.4b network liveness.** At `χ = 1`, `P_ℓd ≡ 0` for the loop sizes that
   exist (`ρ_net = 1e14` → `L_ld = 100 nm` vs `χd = 6.5 nm`, i.e. `tanh(−315)`).
   Screening `w_c`, `χ`, `K_rec` on such a grid returns `S_i = 0` for all three
   and looks exactly like a genuine screening result. `χ` nominal has been moved
   1.0 → 5.0 accordingly.
6. ~~**`bin_moment` + `loop_conversion` is broken** (T0.5(c5)) and fails
   silently.~~ **Stale — corrected in plan §10(c).** The recorded diagnosis was
   wrong. The bin-moment ⟨100⟩ reconstruct→transfer→project block *is*
   implemented; the real defect was `post_process.py` flooring `μ₀` and `μ₁`
   independently, which corrupted **every** `bin_moment` run through all three
   populations, conversion on or off. Fixed 2026-08-02 by a shared
   `_floor_bin_moments()`. ⟨100⟩ under `bin_moment` at genuine refinement is
   therefore **untested, not broken** — it is gate §11(i)-1 of the plan and the
   assumption the rev-6 campaign rests on. `bin_moment` results predating
   2026-08-02 remain suspect.

## Status — T2-v2 stopped 2026-08-03, superseded by plan §11

The v2 campaign (`discrete`, `I=3200`, `V=600`, 0.1 dpa) was halted at 275 of
1104 rows having produced **12 admissible rows from 182 core-hours**. Root cause:
the grid was certified against two values of a *single* parameter
(`dH2_abs_conv`) and then run over a 21-parameter box in which `L_hat` — spanning
three decades and never varied in that study — is the strongest driver of grid
demand.

The replacement campaign is specified in **plan §11**: `bin_moment`,
`I = V = 10000`, `i_mobile = 50`, `v_mobile = 5`, 25 + 25 bins, 1 dpa, `p = 19`.
**Do not launch it before the six pre-flight gates in §11(i) pass** — in
particular the ⟨100⟩ closure validation and a regenerated
`machine_reference.json`, which cannot certify a bin-moment build as it stands.

Also note, from plan §11(c)-4: **`active_window` is not free in `I`.** A
bit-identical trajectory costs 14× more going from `I = 3200` to `I = 25600`,
because CVODE's vector arithmetic runs over the full `N_eq ≈ 2I + V` even though
the RHS correctly masks to the sliding window. This is why the campaign moves to
`bin_moment` rather than simply enlarging the discrete grid.
