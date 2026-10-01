# Cross-machine score analysis

Generated 2026-08-19 from all `results/*.jsonl` at the time of writing.
Regenerate by re-running the analysis in the commit that added this file.

## Per-machine contribution

| machine_id | rows | valid | stages | median h/row | best (log-dist) | best row |
|---|---|---|---|---|---|---|
| Nasr-Workstation | 579 | 9 | 41 | 3.35 | 2/6 (1.312) | E02_V0_Ea2.30_Eb0.75 |
| MacBook-Pro.local | 516 | 19 | 4 | 1.02 | 3/6 (1.240) | V109_Zv_lo0.9 |
| Mac.san.rr.com | 264 | 3 | 3 | 1.58 | 2/6 (1.312) | S1511_eta0.249_fcli0.124_Ebi20.75 |
| MATRIX-PC2 | 121 | 31 | 4 | 8.67 | 3/6 (1.259) | S1700_imobile5_si1_phimaxjunc0.05 |

`Mac.san.rr.com` and `MacBook-Pro.local` are the SAME physical machine under two
hostnames (registry index 0).  Its real totals are 780 rows / 22 valid.  This is
also why `plan.py estimate_cost` never matches either against the registry name
"MacBook Pro" and silently falls back to the scaled worst case.

MATRIX-PC2 has the highest valid-row yield (31/121 = 26%) despite being 8.5x the
slowest per row: its stages are better targeted, not luckier.

## Reproducibility

A replicate requires the SAME theta_hash AND condition AND run_cfg_sha.
theta_hash alone covers only the parameter vector, so matching on it yields 284
spurious "replicates" -- the same theta under different irradiation conditions
(N2/N5/I1) legitimately gives different answers.  Requiring all three keys
reduces 284 to 2.  Anyone analysing on theta_hash alone would conclude the
solver is wildly non-deterministic.  It is not.

| replicate | machines | max rel deviation |
|---|---|---|
| theta 17913ff5 / S330 / cfg 7f4cd1c5 | MATRIX-PC2 vs MacBook-Pro.local | 9.85e-07 |
| theta e57aeab4 / N5 / cfg 482d5955 | Nasr-Workstation vs MATRIX-PC2 | 0.792 |

The first is the S1700 base reproduced across OS and architecture to ~1e-6 on
all six observables.

The second is NOT a divergence.  Both rows are `starved`: they hit the
wall-clock budget and were cut at DIFFERENT doses -- 0.242 dpa in 28862 s on
Nasr-Workstation, 0.492 dpa in 50008 s on MATRIX-PC2.  Same trajectory sampled
at two points, not the same calculation twice.  Both are also grid_limited with
delta_FP ~0.55-0.61, so `classify()` rejects them on GRID/NOCONV/STARV and
neither has ever entered a verdict.

CONCLUSION: of the two true replicates in the campaign, the one that ran to
completion agrees to 1e-6 across machines.  There is no evidence of
cross-machine numerical divergence.
