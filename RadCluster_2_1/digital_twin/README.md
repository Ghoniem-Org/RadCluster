# Digital-twin campaign harness — 4-machine protocol

> **Start here: `campaign_control.ipynb`.** Open it on each machine. It
> configures the run, auto-builds the solver, verifies this machine agrees with
> the reference, launches the worker, gives a live monitor against the
> experimental data, and provides graceful stop/restart. The command lines
> below are what it drives — use them directly if you prefer a shell.

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
| `occ_111`, `occ_100` | mean size / `I` — project rule: `> 0.1` is suspect |
| `pile_111`, `pile_100` | content fraction in the top 2 % of the grid |
| `d_over_ceiling_100` | `d_100 / d(n=I)` |
| `dose_reached`, `starved` | did it reach the requested dose |
| `delta_FP`, `delta_He` | conservation |
| `admissible` | all of the above passed |

A row that is **inadmissible is not a failed row** — it ran fine and conserves,
but its observables measure the grid rather than the physics. The smoke test at
`I = 150` produced exactly this: `pile = 0.996`, correctly rejected.

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
6. **`bin_moment` + `loop_conversion` is broken** (T0.5(c5)) and fails silently.
   Discrete grids large enough for the calibrated ⟨100⟩ sizes cost ~3 h/run,
   which does not scale to this campaign. This is the gating item.
