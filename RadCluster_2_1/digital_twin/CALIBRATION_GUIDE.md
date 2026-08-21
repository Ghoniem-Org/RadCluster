# Calibration Guide - EUROFER97 digital twin

*Derived by `learn.py`; content last changed 2026-08-21 14:24:38Z (re-running with no new results leaves this file untouched). Do not hand-edit: edits are overwritten. Durable notes belong in `calibration_ledger.json` under `notes`, which is preserved across regenerations.*

## Goal

Find one parameter vector that puts all six observables inside their experimental ranges at EUROFER97 330 C / 15 dpa, neutron.

Condition: **EUROFER97, 330 C, 15.0 dpa, neutron** - targets in `targets_330C_15dpa.json`.

## Where the campaign stands

Best valid row: **7403** from `W3a_V5k_machine0` - **4/6** observables in range, log-distance 0.811.

| observable | model | target | range | ratio | in range |
|---|---|---|---|---|---|
| N_100 | 8.77e+21 | 4.97e+21 | 4.67e+21 - 9e+21 | 1.76 x | yes |
| d_100 | 6.72 | 6.2 | 3.4 - 7 | 1.08 x | yes |
| N_111 | 1.43e+22 | 1.93e+21 | 1.73e+21 - 1.5e+22 | 7.42 x | yes |
| d_111 | 1.05 | 6.2 | 3.4 - 7 | 0.169 x | **no** |
| N_void | 3.84e+19 | 1.5e+21 | 3.6e+20 - 3.01e+21 | 0.0256 x | **no** |
| d_void | 2.3 | 2.6 | 2.12 - 2.9 | 0.883 x | yes |

### Defect inventory (SIA content locked in loops)

Computed from the row own n<->d relation, so no lattice constant is assumed.

| character | model N x nbar | experiment N x nbar | ratio |
|---|---|---|---|
| <100> | 7.55e+24 | 3.65e+24 | 2.07 x |
| <111> | 2.61e+23 | 1.23e+24 | 0.213 x |
| **total** | 7.81e+24 | 4.87e+24 | **1.6 x** |

## Can any lever still close this?

How far each observable must move from the best row, against the largest change any SINGLE lever has actually produced across its full tested span. Levers can compound, so a `no` is a strong hint that the residual is structural rather than a proof that it is.

| observable | ratio | needs | best single lever | that lever | single lever enough? |
|---|---|---|---|---|---|
| N_100 | 1.76 x | +76% | +16301094% | `B_111` | yes |
| d_100 | 1.08 x | +8% | +693% | `f_cl_v+E_b_v2+s_v` | yes |
| N_111 | 7.42 x | +642% | +113054% | `B_111` | yes |
| d_111 | 0.169 x | +491% | +89% | `s_i` | **no** |
| N_void | 0.0256 x | +3806% | +4111087% | `E_m_v` | yes |
| d_void | 0.883 x | +13% | +875% | `f_cl_v+E_m_v+gamma_s` | yes |

## What each lever does

A lever is **dead** when its full tested span moves every observable by less than 5%. Only pairs of rows differing in exactly one lever are used, and only rows that are real measurements (full dose, grid-clean, converged).

`inconclusive` means the span was attempted but no pair of valid rows bracketed it - usually because the rows piled at the grid ceiling or never reached dose. An inconclusive lever is an OPEN question, not a closed one.

| lever | tested span | verdict | valid pairs | peak response | moves | stages |
|---|---|---|---|---|---|---|
| `B_111` | B_111 0.2->0.6 | **live** (drives rows off-grid) | 23 | 16301094.4% | N_100, N_111, d_100, d_111, d_void | S16_calib_machine1, S20_calib_machine0, V13_calib_machine0 |
| `E_m_v` | E_m_v 0.5->1 | **live** (drives rows off-grid) | 18 | 4111086.8% | N_void, N_100, d_void, N_111, d_100 | V12_calib_machine0, V1_calib_machine0, V2_calib_machine0, V5_calib_machine0, V6_calib_machine0, V7_calib_machine0, V8_calib_machine0 |
| `E_m_i` | E_m_i 0.25->0.6 | **live** (drives rows off-grid) | 22 | 41468.3% | N_100, N_111, d_100, N_void, d_111, d_void | S16_calib_machine1, S20_calib_machine0, V13_calib_machine0, V15_calib_machine0 |
| `L_hat` | L_hat 10->3e+03 | **live** (drives rows off-grid) | 24 | 30324.9% | N_100, N_111, d_100, d_111, d_void | S16_calib_machine1, S20_calib_machine0 |
| `f_cl_v+E_m_v+gamma_s` | f_cl_v 0.317->0.55; E_m_v 0.594->0.8; gamma_s 2.22->2.33 | **live** | 1 | 29161.7% | N_void, d_void, d_100, N_111 | V9_calib_machine0 |
| `f_cl_v+E_m_v+gamma_s+dH2_abs_conv+rho_d` | f_cl_v 0.317->0.55; E_m_v 0.594->0.8; gamma_s 2.22->2.33; dH2_abs_conv 0.26->0.32; rho_d 1e+14->5e+14 | **live** | 1 | 3156.0% | N_void, N_111, d_void, N_100, d_111, d_100 | V11_calib_machine0 |
| `rho_d` | rho_d 1e+13->1e+16 | **live** (drives rows off-grid) | 15 | 1059.9% | N_100, N_111, d_100, N_void, d_void, d_111 | V12_calib_machine0, V1_calib_machine0, V5_calib_machine0, V6_calib_machine0, V7_calib_machine0, V8_calib_machine0, V9_calib_machine0 |
| `E_b_i2` | E_b_i2 0.55->0.875 | **live** (drives rows off-grid) | 6 | 991.2% | N_100, N_void, N_111, d_void, d_111, d_100 | S10_calib, S12_calib, S13_calib, S14_calib, S15_calib_machine0, S8_15dpa, V13_calib_machine0, V15_calib_machine0 |
| `f_cl_v+E_b_v2+s_v` | f_cl_v 0.317->0.7; E_b_v2 0.213->0.35; s_v 1.9->2.5 | **live** (drives rows off-grid) | 2 | 692.8% | d_100, N_void, N_111 | S14_calib |
| `f_cl_v` | f_cl_v 0.05->0.8 | **live** (drives rows off-grid) | 7 | 449.6% | N_void, N_100, N_111, d_100, d_void | S18_calib_machine2, V12_calib_machine0, V1_calib_machine0, V8_calib_machine0 |
| `s_i` | s_i 1->3 | **live** | 12 | 436.5% | d_100, N_111, d_111, N_100 | S17_calib_machine1 |
| `f_cl_i` | f_cl_i 0.02->0.124 | **live** (drives rows off-grid) | 1 | 426.8% | d_100, N_100, d_void, N_111, N_void, d_111 | S15_calib_machine0, S8_15dpa, V10_calib_machine0 |
| `eta` | eta 0.05->0.249 | **live** (drives rows off-grid) | 3 | 397.2% | N_void, N_111, N_100, d_100, d_111 | S15_calib_machine0 |
| `gamma_s` | gamma_s 1.5->3 | **live** (drives rows off-grid) | 20 | 395.8% | d_100, N_100, N_void, N_111, d_void, d_111 | V10_calib_machine0, V1_calib_machine0, V2_calib_machine0, V5_calib_machine0, V6_calib_machine0, V7_calib_machine0, V8_calib_machine0 |
| `dH2_abs_conv` | dH2_abs_conv 0.26->0.55 | **live** (drives rows off-grid) | 19 | 173.8% | N_void, N_111, N_100, d_100, d_void, d_111 | V10_calib_machine0, V13_calib_machine0, V15_calib_machine0, V9_calib_machine0 |
| `phi_max_junc` | phi_max_junc 0.05->0.6 | **live** | 6 | 163.9% | N_100, d_100, N_111, d_111 | S17_calib_machine1 |
| `i_mobile` | i_mobile 5->100 | **live** (drives rows off-grid) | 8 | 132.7% | N_100, N_111, d_100, d_111 | S17_calib_machine1, S8_15dpa, V13_calib_machine0 |
| `loop_net_w_c` | loop_net_w_c 5->2.7e+03 | **live** | 32 | 96.0% | N_100, d_100, N_void, N_111, d_void, d_111 | V11_calib_machine0, V13_calib_machine0, V15_calib_machine0 |
| `Z_i` | Z_i 1->1.35 | **live** | 2 | 86.6% | N_111, N_100, d_100, N_void, d_void, d_111 | V5_calib_machine0, V6_calib_machine0 |
| `Z_v` | Z_v 0.9->1.35 | **live** | 3 | 65.2% | N_111, N_100, d_100 | V1_calib_machine0, V5_calib_machine0, V6_calib_machine0 |
| `E_b_hV_1` | E_b_hV_1 2.21->3 | **live** | 5 | 28.7% | N_100, d_100 | S18_calib_machine2, V11_calib_machine0 |
| `s_v` | s_v 2.5->3 | **live** (drives rows off-grid) | 2 | 23.2% | N_100, N_void, d_100 | V12_calib_machine0, V8_calib_machine0 |
| `loop_net_chi` | loop_net_chi 50->1e+03 | **live** | 7 | 22.3% | N_100, d_100 | V11_calib_machine0 |
| `E_a0_conv` | E_a0_conv 1.8->2.3 | **dead** (drives rows off-grid) | 2 | 2.5% | - | S10_calib, S12_calib, S13_calib, S14_calib, S8_15dpa |
| `dH2_conv` | dH2_conv 0.439->0.75 | **dead** | 3 | 0.0% | - | V9_calib_machine0 |
| `Z_i+Z_p_i+Z_p_v+Z_gb_i+Z_gb_v` | Z_i 1.2->1.35; Z_p_i 1.2->1.45; Z_p_v 0.93->0.98; Z_gb_i 1.2->1.4; Z_gb_v 0.93->0.98 | **inconclusive** (drives rows off-grid) | 0 | - | - | S12_calib |
| `d_g+rho_p` | d_g 1e-06->4e-06; rho_p 1e+19->1e+20 | **inconclusive** (drives rows off-grid) | 0 | - | - | S13_calib |
| `E_b_v2` | E_b_v2 0.1->0.45 | **dead** | 5 | 0.0% | - | S18_calib_machine2, V1_calib_machine0 |
| `f_cl_v+E_m_v+rho_d` | f_cl_v 0.317->0.55; E_m_v 0.594->0.8; rho_d 1e+14->5e+14 | **inconclusive** | 0 | - | - | V10_calib_machine0 |

### Never varied

These columns exist in the design but have never taken more than one value, so the campaign has no evidence about them:

`E_m_h`, `Z_i_loop`, `lambda`, `loop_net_K_rec`, `r_p`

## Stage history

Only stages that ran to 15.0 dpa are in scope; the rest are listed for provenance but cast no vote on any lever.

| stage | dose | scope | rows | valid | swept | best in-range |
|---|---|---|---|---|---|---|
| `AUTOVERIFY_7100` | 15 | in | 1 | 1 | - | 4 |
| `G1a_vb45_machine0` | 15 | in | 1 | 1 | - | 1 |
| `G1b_V80k_machine0` | 15 | in | 1 | 1 | - | 1 |
| `G1c_V320k_machine0` | 15 | in | 1 | 1 | - | 1 |
| `G2a_vb45_machine0` | 15 | in | 1 | 0 | - | - |
| `G2b_V80k_machine0` | 15 | in | 1 | 1 | - | 3 |
| `G3a_V20k_machine0` | 15 | in | 5 | 5 | - | 4 |
| `G3b_V80k_machine0` | 15 | in | 5 | 5 | - | 4 |
| `G4_I160k_machine0` | 15 | in | 1 | 1 | - | 4 |
| `G5a_disc_V1k` | 15 | in | 1 | 0 | - | - |
| `G5b_disc_V4k` | 15 | in | 1 | 0 | - | - |
| `G5c_bm_V1k` | 15 | in | 1 | 0 | - | - |
| `G5d_bm_V4k` | 15 | in | 1 | 0 | - | - |
| `R1_repro_machine0` | 15 | in | 3 | 3 | - | 2 |
| `S10_calib` | 15 | in | 12 | 0 | `E_b_i2`, `E_a0_conv` | - |
| `S11_active_window` | 15 | in | 1 | 0 | - | - |
| `S11_full_system` | 15 | in | 1 | 0 | - | - |
| `S11eq_active_window` | 0.002 | out | 1 | 0 | - | - |
| `S11eq_full_system` | 0.002 | out | 1 | 1 | - | 1 |
| `S12_calib` | 15 | in | 12 | 0 | `E_b_i2`, `Z_i+Z_p_i+Z_p_v+Z_gb_i+Z_gb_v`, `E_a0_conv` | - |
| `S13_calib` | 15 | in | 12 | 0 | `E_b_i2`, `E_a0_conv`, `d_g+rho_p` | - |
| `S14_calib` | 15 | in | 12 | 4 | `f_cl_v+E_b_v2+s_v`, `E_b_i2`, `E_a0_conv` | 2 |
| `S15_calib_machine0` | 15 | in | 6 | 3 | `eta`, `f_cl_i`, `E_b_i2` | 2 |
| `S16_calib_machine1` | 15 | in | 28 | 21 | `E_m_i`, `L_hat`, `B_111` | 2 |
| `S17_calib_machine1` | 15 | in | 13 | 13 | `phi_max_junc`, `i_mobile`, `s_i` | 3 |
| `S18_calib_machine2` | 15 | in | 4 | 4 | `f_cl_v`, `E_b_v2`, `E_b_hV_1` | 2 |
| `S1_calib` | 1 | out | 2 | 0 | - | - |
| `S1_lnl0` | 1 | out | 12 | 1 | - | 0 |
| `S1_lnl1` | 1 | out | 4 | 0 | - | - |
| `S1_smoke` | 15 | in | 1 | 0 | - | - |
| `S20_calib_machine0` | 15 | in | 13 | 12 | `E_m_i`, `L_hat`, `B_111` | 3 |
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
| `T10_ea0_scan` | 2 | out | 24 | 2 | `E_a0_conv` | 1 |
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
| `T7_cost` | 1 | out | 4 | 3 | - | 2 |
| `T9_crossover_batch1` | 16.3 | out | 12 | 0 | - | - |
| `T9_ladder_validation` | 0.05 | out | 4 | 2 | - | 1 |
| `T9_screen_partial` | 2 | out | 13 | 1 | - | 0 |
| `TEST_modes` | 15 | in | 1 | 0 | - | - |
| `V10_calib_machine0` | 15 | in | 12 | 6 | `f_cl_i`, `f_cl_v+E_m_v+rho_d`, `gamma_s`, `dH2_abs_conv` | 3 |
| `V11_calib_machine0` | 15 | in | 13 | 13 | `f_cl_v+E_m_v+gamma_s+dH2_abs_conv+rho_d`, `E_b_hV_1`, `loop_net_chi`, `loop_net_w_c` | 3 |
| `V12_calib_machine0` | 15 | in | 5 | 3 | `f_cl_v`, `E_m_v`, `rho_d`, `s_v` | 4 |
| `V13_calib_machine0` | 15 | in | 9 | 9 | `E_m_i`, `B_111`, `E_b_i2`, `dH2_abs_conv`, `i_mobile`, `loop_net_w_c` | 2 |
| `V14a_V20k_machine0` | 15 | in | 3 | 3 | - | 4 |
| `V14b_V80k_machine0` | 15 | in | 3 | 2 | - | 4 |
| `V15_calib_machine0` | 15 | in | 13 | 13 | `E_m_i`, `E_b_i2`, `dH2_abs_conv`, `loop_net_w_c` | 4 |
| `V1_calib_machine0` | 15 | in | 13 | 10 | `f_cl_v`, `E_m_v`, `gamma_s`, `E_b_v2`, `Z_v`, `rho_d` | 3 |
| `V2_calib_machine0` | 15 | in | 11 | 7 | `E_m_v`, `gamma_s` | 4 |
| `V5_calib_machine0` | 15 | in | 4 | 4 | `E_m_v`, `gamma_s`, `Z_i`, `Z_v`, `rho_d` | 3 |
| `V6_calib_machine0` | 15 | in | 9 | 9 | `E_m_v`, `gamma_s`, `Z_i`, `Z_v`, `rho_d` | 3 |
| `V7_calib_machine0` | 15 | in | 13 | 10 | `E_m_v`, `gamma_s`, `rho_d` | 3 |
| `V8_calib_machine0` | 15 | in | 13 | 12 | `f_cl_v`, `E_m_v`, `gamma_s`, `rho_d`, `s_v` | 4 |
| `V9_calib_machine0` | 15 | in | 13 | 13 | `f_cl_v+E_m_v+gamma_s`, `dH2_conv`, `dH2_abs_conv`, `rho_d` | 4 |
| `W1a_V20k_machine0` | 15 | in | 13 | 7 | - | 3 |
| `W2a_V5k_machine0` | 15 | in | 2 | 2 | - | 0 |
| `W3a_V5k_machine0` | 15 | in | 6 | 6 | - | 4 |
| `W3b_V20k_machine0` | 15 | in | 3 | 3 | - | 3 |
| `W5a_V5k_machine0` | 15 | in | 6 | 6 | - | 4 |
| `W5b_V20k_machine0` | 15 | in | 2 | 2 | - | 3 |
| `W6a_V5k_machine0` | 15 | in | 5 | 4 | - | 2 |
| `W6b_V20k_machine0` | 15 | in | 3 | 3 | - | 2 |

## Cost model (measured)

Completed rows only. A row cut at the budget measures the timeout, not the cost, so timeouts are counted in their own column.

| machine | completed | median row | min | max | timed out |
|---|---|---|---|---|---|
| MATRIX-PC2 | 117 | 8.26 h | 2.94 h | 18.87 h | 5 |
| Mac.san.rr.com | 366 | 1.43 h | 0.04 h | 5.51 h | 1 |
| MacBook-Pro.local | 587 | 1.02 h | 0.02 h | 9.2 h | 2 |
| Nasr-Workstation | 415 | 3.35 h | 0.15 h | 19.79 h | 164 |

`plan.py` sizes a stage from this table and the machine slot count, and refuses to propose a design whose estimated cost exceeds the machine row budget.

## Unaffordable lever values

Measured, not assumed: a level counts here only if it produced ZERO full-dose rows while another level of the same lever - everything else held equal - produced at least one. `plan.py` will not place a design point on these values.

| lever | unaffordable at | evidence |
|---|---|---|
| `E_b_i2` | 0.6 | S14_calib: 0 of 6 rows reached dose at E_b_i2=0.6 |

## Parameters the physics never reads

A design can set these and the solver ignores them, so sweeping one produces bit-identical rows. `plan.py` refuses to propose them. These are BUGS to fix, not physics conclusions - remove the entry once the code reads the parameter.

- **`E_b_v2`** - S18 rows 4801 and 4803 differ only in E_b_v2 (0.10 vs 0.45 eV, the full prior box) and agree in 63 of 67 numeric fields BIT-FOR-BIT -- N_voids, mean_n_v and S_inventory to the last digit; only row_id bookkeeping and wall_s differ. At 603 K a 0.35 eV binding change should move di-vacancy emission by exp(0.35/0.052) ~ 840x. Confirmed in source: `E_b_v2` appears only in py_utils/create_excel.py, which BUILDS the workbook; binding_energies.E_b_void() derives the m=2 binding from the capillary term plus A_void_0*exp(-lambda*(m-1)) and never reads it. Note the same function's docstring records that `lambda` and `A_void_0` were once carried by the workbook and unread -- same bug, previously fixed. FIX THE CODE, then delete this entry; do not sweep E_b_v2 until then.

## Deferred observables

The planner will not propose levers for these until the entry is removed from `calibration_ledger.json` under `policy.deprioritized_observables`. They still appear in every score above - they are deferred, not ignored.

- **N_void** - SUPERSEDED REASON, KEPT AS RECORD: this entry previously claimed the cavity channel was structurally defective because the S14 vacancy triple moved N_void only 3.5x. That inference was WRONG and it was mine. The triple contains E_b_v2, which is wired to nothing, and f_cl_v, which is genuinely weak (2.5x over its whole box). V1 varied one lever at a time and found N_void responds by up to 4.8e4x -- to gamma_s and E_m_v, the two levers never varied in 19 stages. Unresponsiveness was an artefact of which levers had been tried, not a property of the channel. CURRENT REASON TO DEFER, WHICH IS A BLOCKER NOT A VERDICT: every setting that grows realistic cavities exceeds the V=2000 grid. gamma_s=1.5 and E_m_v=0.9 are mechanistically unrelated yet both halt at mean_n_v 833/819, occ_v 0.416/0.409, GRID+NOCONV. Those rows are invalid, so the large response is real but unmeasurable and cannot enter the reachability table. Re-enable when a stage runs at V ~ 20000, sampling the bracket INTERIORS (gamma_s 1.5-2.22, E_m_v 0.593-0.90) that no row has visited.
- **d_void** - Same blocker. d_void is a pure function of mean_n_v (d = 0.2824*m^(1/3) to four digits on every row), so a d_void near 2.65 on a grid-limited row is the CEILING, not agreement with the 2.12-2.9 band. Any d_void that looks in range at V=2000 must be treated as an artefact until re-run at V ~ 20000.

## Claimed stages

One row per machine. A machine claims a stage by running `plan.py --write` on it; the claim is not a lock, only a record of what that machine was last told to run.

| machine | stage | levers | rows | design |
|---|---|---|---|---|
| 0 (MacBook Pro) | **S20** | `L_hat`, `B_111`, `E_m_i` | 12 | `design/S20_calib.csv` |
| 1 (Matrix-PC) | **S22** | `s_i`, `eta`, `L_hat` | 20 | `design/S22_calib.csv` |
| 2 (Nasr Workstation) | **V3** | `gamma_s`, `E_m_v` | 12 | `design/V3_calib.csv` |

## Next stage

**V3** - worst residuals are d_111 (0.165x), N_111 (5.15x), N_100 (1.27x); sweeping gamma_s, E_m_v, which the ledger has not retired

Sweeping: `gamma_s`, `E_m_v`

Design: `design/V3_calib.csv` (12 rows)

Run it with:

```bash
python run_ensemble.py \
    --design design/V3_calib.csv \
    --conditions conditions_S8.json \
    --spec parameters_S4.json \
    --out results/V3_calib_machine2.jsonl \
    --machine 0 --of 1 \
    --equations bin_moment --i-discrete 100 --i-bin 36 \
    --v-discrete 5 --v-bin 20 --allow-mixed \
    --I 80000 --V 20000 --dose 15.0 --lnl 1 --rtol 1e-5 \
    --solver-mode full_system \
    --timeout-s 72000 --workers 12 --omp-threads 1
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

