# Calibration Guide - EUROFER97 digital twin

*Derived by `learn.py`; content last changed 2026-08-18 22:42:23Z (re-running with no new results leaves this file untouched). Do not hand-edit: edits are overwritten. Durable notes belong in `calibration_ledger.json` under `notes`, which is preserved across regenerations.*

## Goal

Find one parameter vector that puts all six observables inside their experimental ranges at EUROFER97 330 C / 15 dpa, neutron.

Condition: **EUROFER97, 330 C, 15.0 dpa, neutron** - targets in `targets_330C_15dpa.json`.

## Where the campaign stands

Best valid row: **E00_V0_Ea2.10_Eb0.75** from `S14_calib` - **2/6** observables in range, log-distance 1.312.

| observable | model | target | range | ratio | in range |
|---|---|---|---|---|---|
| N_100 | 6.16e+22 | 4.97e+21 | 4.67e+21 - 9e+21 | 12.4 x | **no** |
| d_100 | 4.86 | 6.2 | 3.4 - 7 | 0.785 x | yes |
| N_111 | 8.92e+21 | 1.93e+21 | 1.73e+21 - 1.1e+22 | 4.62 x | yes |
| d_111 | 1.07 | 6.2 | 3.4 - 7 | 0.172 x | **no** |
| N_void | 2.57e+18 | 1.5e+21 | 3.6e+20 - 3.01e+21 | 0.0017 x | **no** |
| d_void | 0.561 | 2.6 | 2.12 - 2.9 | 0.216 x | **no** |

### Defect inventory (SIA content locked in loops)

Computed from the row own n<->d relation, so no lattice constant is assumed.

| character | model N x nbar | experiment N x nbar | ratio |
|---|---|---|---|
| <100> | 2.78e+25 | 3.65e+24 | 7.63 x |
| <111> | 1.68e+23 | 1.23e+24 | 0.137 x |
| **total** | 2.8e+25 | 4.87e+24 | **5.74 x** |

## What each lever does

A lever is **dead** when its full tested span moves every observable by less than 5%. Only pairs of rows differing in exactly one lever are used, and only rows that are real measurements (full dose, grid-clean, converged).

`inconclusive` means the span was attempted but no pair of valid rows bracketed it - usually because the rows piled at the grid ceiling or never reached dose. An inconclusive lever is an OPEN question, not a closed one.

| lever | tested span | verdict | valid pairs | peak response | moves | stages |
|---|---|---|---|---|---|---|
| `f_cl_v+E_b_v2+s_v` | f_cl_v 0.317->0.7; E_b_v2 0.213->0.35; s_v 1.9->2.5 | **live** (drives rows off-grid) | 1 | 674.9% | d_100, N_void, N_111 | S14_calib |
| `E_a0_conv` | E_a0_conv 1.8->2.3 | **dead** (drives rows off-grid) | 1 | 1.3% | - | S10_calib, S12_calib, S13_calib, S14_calib, S8_15dpa |
| `E_b_i2` | E_b_i2 0.55->0.875 | **inconclusive** (drives rows off-grid) | 0 | - | - | S10_calib, S12_calib, S13_calib, S14_calib, S8_15dpa |
| `Z_i+Z_p_i+Z_p_v+Z_gb_i+Z_gb_v` | Z_i 1.2->1.35; Z_p_i 1.2->1.45; Z_p_v 0.93->0.98; Z_gb_i 1.2->1.4; Z_gb_v 0.93->0.98 | **inconclusive** (drives rows off-grid) | 0 | - | - | S12_calib |
| `d_g+rho_p` | d_g 1e-06->4e-06; rho_p 1e+19->1e+20 | **inconclusive** (drives rows off-grid) | 0 | - | - | S13_calib |
| `f_cl_i` | f_cl_i 0.05->0.124 | **inconclusive** (drives rows off-grid) | 0 | - | - | S8_15dpa |
| `i_mobile` | i_mobile 40->100 | **inconclusive** (drives rows off-grid) | 0 | - | - | S8_15dpa |

### Never varied

These columns exist in the design but have never taken more than one value, so the campaign has no evidence about them:

`B_111`, `E_b_hV_1`, `E_m_h`, `E_m_i`, `E_m_v`, `L_hat`, `Z_i_loop`, `Z_v`, `dH2_abs_conv`, `dH2_conv`, `eta`, `gamma_s`, `lambda`, `loop_net_K_rec`, `loop_net_chi`, `loop_net_w_c`, `phi_max_junc`, `r_p`, `rho_d`, `s_i`

## Stage history

Only stages that ran to 15.0 dpa are in scope; the rest are listed for provenance but cast no vote on any lever.

| stage | dose | scope | rows | valid | swept | best in-range |
|---|---|---|---|---|---|---|
| `S10_calib` | 15 | in | 12 | 0 | `E_b_i2`, `E_a0_conv` | - |
| `S11_active_window` | 15 | in | 1 | 0 | - | - |
| `S11_full_system` | 15 | in | 1 | 0 | - | - |
| `S11eq_active_window` | 0.002 | out | 1 | 0 | - | - |
| `S11eq_full_system` | 0.002 | out | 1 | 1 | - | 1 |
| `S12_calib` | 15 | in | 12 | 0 | `E_b_i2`, `Z_i+Z_p_i+Z_p_v+Z_gb_i+Z_gb_v`, `E_a0_conv` | - |
| `S13_calib` | 15 | in | 12 | 0 | `E_b_i2`, `E_a0_conv`, `d_g+rho_p` | - |
| `S14_calib` | 15 | in | 5 | 3 | `f_cl_v+E_b_v2+s_v`, `E_b_i2`, `E_a0_conv` | 2 |
| `S1_calib` | 1 | out | 2 | 0 | - | - |
| `S1_lnl0` | 1 | out | 12 | 1 | - | 0 |
| `S1_lnl1` | 1 | out | 4 | 0 | - | - |
| `S1_smoke` | 15 | in | 1 | 0 | - | - |
| `S2_cap` | 1 | out | 6 | 0 | - | - |
| `S2_nocap` | 1 | out | 10 | 0 | - | - |
| `S3_nucleation` | 0.02 | out | 8 | 6 | `f_cl_i`, `d_g`, `rho_p+r_p`, `s_i` | 1 |
| `S4_bias` | 0.02 | out | 9 | 8 | `f_cl_i`, `Z_i`, `Z_i_loop`, `Z_p_i+Z_p_v`, `Z_gb_i+Z_gb_v`, `rho_d` | 0 |
| `S4_smoke` | 0.002 | out | 2 | 0 | - | - |
| `S5_binmoment` | 0.005 | out | 2 | 2 | - | 1 |
| `S5_discrete` | 0.005 | out | 2 | 2 | - | 1 |
| `S6_binding` | 0.02 | out | 9 | 7 | `B_111`, `E_b_i2`, `Z_i+Z_i_loop+Z_p_i+Z_p_v+Z_gb_i+Z_gb_v` | 0 |
| `S7_binding` | 0.02 | out | 6 | 0 | - | - |
| `S8_15dpa` | 15 | in | 6 | 0 | `f_cl_i`, `E_b_i2`, `E_a0_conv`, `i_mobile` | - |
| `S9_grid_A` | 0.02 | out | 1 | 1 | - | 0 |
| `S9_grid_B` | 0.02 | out | 1 | 0 | - | - |
| `S9_grid_C` | 0.02 | out | 1 | 0 | - | - |
| `S9_grid_D` | 0.02 | out | 1 | 0 | - | - |
| `T10_ea0_scan` | 2 | out | 24 | 0 | `E_a0_conv` | - |
| `T11_row38_scan` | 2 | out | 12 | 0 | `E_a0_conv` | - |
| `T13_nref` | 16.3 | out | 12 | 0 | `E_a0_conv`, `n_ref_conv` | - |
| `T3_rev6_machine0` | 15 | in | 480 | 0 | - | - |
| `T3_rev6_machine1` | 15 | in | 55 | 0 | - | - |
| `T3_rev6_machine2` | 15 | in | 165 | 0 | - | - |
| `T3_rev6_machine2_t0` | 15 | in | 80 | 0 | - | - |
| `T3_rev6_machine2_t1` | 15 | in | 80 | 0 | - | - |
| `T3_rev6_machine2_t2` | 15 | in | 80 | 0 | - | - |
| `T3_rev6_machine3` | 15 | in | 233 | 0 | - | - |
| `T4_lnl_check_n` | 1 | out | 1 | 0 | `loop_net_K_rec+loop_net_w_c+loop_net_chi` | - |
| `T5_cost218` | 1 | out | 2 | 0 | - | - |
| `T6_fixcheck` | 1 | out | 2 | 0 | `L_hat+Z_i_loop` | - |
| `T7_cost` | 1 | out | 4 | 2 | - | 1 |
| `T9_crossover_batch1` | 16.3 | out | 12 | 0 | - | - |
| `T9_ladder_validation` | 0.05 | out | 4 | 2 | - | 1 |
| `T9_screen_partial` | 2 | out | 13 | 1 | - | 0 |

## Cost model (measured)

| machine | rows timed | median row | min | max |
|---|---|---|---|---|
| MATRIX-PC2 | 81 | 6.5 h | 2.94 h | 14.02 h |
| Mac.san.rr.com | 258 | 1.58 h | 0.33 h | 4.73 h |
| MacBook-Pro.local | 480 | 1.02 h | 0.23 h | 9.2 h |
| Nasr-Workstation | 568 | 3.35 h | 0.04 h | 16.07 h |

`plan.py` sizes a stage from this table and the machine slot count, and refuses to propose a design whose estimated cost exceeds the machine row budget.

## Deferred observables

The planner will not propose levers for these until the entry is removed from `calibration_ledger.json` under `policy.deprioritized_observables`. They still appear in every score above - they are deferred, not ignored.

- **N_void** - Missed by ~300x AND nearly unresponsive to its own governing parameters: the S14 vacancy triple moved N_void only 2.57e18 -> 9.08e18 (3.5x) while moving loop content 160x. A residual that large with a response that small is evidence of a structural defect in the cavity channel, not a parameter that needs tuning. Re-enable once cavity nucleation is shown to respond at all.
- **d_void** - Pinned at 0.56-0.57 nm across every row of every stage, including a 160x swing in loop content. Same reasoning as N_void.

## Next stage

**S16** - worst residuals are N_100 (12.4x), d_111 (0.172x), N_111 (4.62x); sweeping B_111, L_hat, E_m_i, which the ledger has not retired

Sweeping: `B_111`, `L_hat`, `E_m_i`

Design: `design/S16_calib.csv` (36 rows)

Run it with:

```bash
python run_ensemble.py \
    --design design/S16_calib.csv \
    --conditions conditions_S8.json \
    --spec parameters_S4.json \
    --out results/S16_calib_machine1.jsonl \
    --machine 0 --of 1 \
    --equations bin_moment --i-discrete 100 --i-bin 36 \
    --v-discrete 5 --v-bin 20 --allow-mixed \
    --I 80000 --V 2000 --dose 15.0 --lnl 1 --rtol 1e-5 \
    --solver-mode full_system \
    --timeout-s 72000 --workers 40 --omp-threads 1
```

## Multi-machine protocol

Results travel through git; `campaign_ops.sync_results()` is safe to call while a run is in flight (rows are appended, so the newest lands in the next sync).

```bash
git pull --rebase origin main            # 1. take everyone results
python learn.py                          # 2. re-derive ledger + this guide
python plan.py                           # 3. propose this machine next stage
#    ... review the printed command, then run it ...
python -c "import campaign_ops as c; c.sync_results()"   # 4. publish rows
```

Two files must never be committed - both are in `.gitignore` and both have caused an outage before: `CAMPAIGN_STOP` (every machine that pulled it halted at startup, commit 8a600c7) and `results/*.pid`.

