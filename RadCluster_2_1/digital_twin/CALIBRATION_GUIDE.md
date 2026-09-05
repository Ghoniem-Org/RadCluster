# Calibration Guide - EUROFER97 digital twin

*Derived by `learn.py`; content last changed 2026-09-05 18:03:44Z (re-running with no new results leaves this file untouched). Do not hand-edit: edits are overwritten. Durable notes belong in `calibration_ledger.json` under `notes`, which is preserved across regenerations.*

## Goal

Find one parameter vector that puts all six observables inside their experimental ranges at EUROFER97 330 C / 15 dpa, neutron.

Condition: **EUROFER97, 330 C, 15.0 dpa, neutron** - targets in `targets_330C_15dpa.json`.

## Where the campaign stands

Best **extent-verified** row: **7700** from `W6a_V5k_machine0` - **2/6** observables in range, log-distance 0.65. This leaderboard ranks on the verified count first, so a row whose in-band claims have not survived a change of grid extent cannot take the top slot no matter how many bands it sits in.

| observable | model | target | range | ratio | in range |
|---|---|---|---|---|---|
| N_100 | 1.01e+22 | 4.97e+21 | 4.67e+21 - 9e+21 | 2.04 x | **no** |
| d_100 | 6.3 | 6.2 | 3.4 - 7 | 1.02 x | yes |
| N_111 | 1.44e+22 | 1.93e+21 | 1.73e+21 - 1.5e+22 | 7.44 x | yes |
| d_111 | 1.04 | 6.2 | 3.4 - 7 | 0.168 x | **no** |
| N_void | 1.38e+20 | 1.5e+21 | 3.6e+20 - 3.01e+21 | 0.0922 x | **no** |
| d_void | 3.1 | 2.6 | 2.12 - 2.9 | 1.19 x | **no** |

### Defect inventory (SIA content locked in loops)

Computed from the row own n<->d relation, so no lattice constant is assumed.

| character | model N x nbar | experiment N x nbar | ratio |
|---|---|---|---|
| <100> | 7.69e+24 | 3.65e+24 | 2.11 x |
| <111> | 2.58e+23 | 1.23e+24 | 0.211 x |
| **total** | 7.95e+24 | 4.87e+24 | **1.63 x** |

### Best row ignoring verification

Row **9305** from `B3_coal` scores **4/6** in range at log-distance 0.236 - better on the raw score than the verified leader above, but it is UNVERIFIED (no extent pair exists yet), so the ranking will not promote it.

| observable | model | target | range | ratio | in range |
|---|---|---|---|---|---|
| N_100 | 2.18e+21 | 4.97e+21 | 4.67e+21 - 9e+21 | 0.439 x | **no** |
| d_100 | 6.67 | 6.2 | 3.4 - 7 | 1.08 x | yes |
| N_111 | 3.6e+21 | 1.93e+21 | 1.73e+21 - 1.5e+22 | 1.86 x | yes |
| d_111 | 4.54 | 6.2 | 3.4 - 7 | 0.732 x | yes |
| N_void | 1.62e+21 | 1.5e+21 | 3.6e+20 - 3.01e+21 | 1.08 x | yes |
| d_void | 5.61 | 2.6 | 2.12 - 2.9 | 2.16 x | **no** |

Treat this as a lead, not a result: the verification gate exists because an unverified row can be a grid artefact (d_cavity = 0.2825*(0.37V)^(1/3) once took the top of this leaderboard for twelve stages). To promote it, re-run it at a second V extent.

## Can any lever still close this?

How far each observable must move from the best row, against the largest change any SINGLE lever has actually produced across its full tested span. Levers can compound, so a `no` is a strong hint that the residual is structural rather than a proof that it is.

| observable | ratio | needs | best single lever | that lever | single lever enough? |
|---|---|---|---|---|---|
| N_100 | 2.04 x | +104% | +16301094% | `B_111` | yes |
| d_100 | 1.02 x | +2% | +693% | `f_cl_v+E_b_v2+s_v` | yes |
| N_111 | 7.44 x | +644% | +113054% | `B_111` | yes |
| d_111 | 0.168 x | +494% | +113% | `loop_coal_pref` | **no** |
| N_void | 0.0922 x | +985% | +4111087% | `E_m_v` | yes |
| d_void | 1.19 x | +19% | +875% | `f_cl_v+E_m_v+gamma_s` | yes |

## What each lever does

A lever is **dead** when its full tested span moves every observable by less than 5%. Only pairs of rows differing in exactly one lever are used, and only rows that are real measurements (full dose, grid-clean, converged).

`inconclusive` means the span was attempted but no pair of valid rows bracketed it - usually because the rows piled at the grid ceiling or never reached dose. An inconclusive lever is an OPEN question, not a closed one.

| lever | tested span | verdict | valid pairs | peak response | moves | stages |
|---|---|---|---|---|---|---|
| `B_111` | B_111 0.2->0.6 | **live** (drives rows off-grid) | 23 | 16301094.4% | N_100, N_111, d_100, d_111, d_void | S16_calib_machine1, S20_calib_machine0, V13_calib_machine0 |
| `E_m_v` | E_m_v 0.5->1 | **live** (drives rows off-grid) | 18 | 4111086.8% | N_void, N_100, d_void, N_111, d_100 | S4_cavity, V12_calib_machine0, V1_calib_machine0, V2_calib_machine0, V5_calib_machine0, V6_calib_machine0, V7_calib_machine0, V8_calib_machine0 |
| `n_ref_conv` | n_ref_conv 10->100 | **live** (drives rows off-grid) | 8 | 1051551.1% | N_100, d_100, N_111, d_111 | M4_nref, S3_lifetime |
| `E_m_i` | E_m_i 0.25->0.6 | **live** (drives rows off-grid) | 22 | 41468.3% | N_100, N_111, d_100, N_void, d_111, d_void | S16_calib_machine1, S20_calib_machine0, V13_calib_machine0, V15_calib_machine0 |
| `L_hat` | L_hat 10->3e+03 | **live** (drives rows off-grid) | 24 | 30324.9% | N_100, N_111, d_100, d_111, d_void | S16_calib_machine1, S20_calib_machine0, S2_coarsen |
| `f_cl_v+E_m_v+gamma_s` | f_cl_v 0.317->0.55; E_m_v 0.594->0.8; gamma_s 2.22->2.33 | **live** | 1 | 29161.7% | N_void, d_void, d_100, N_111 | V9_calib_machine0 |
| `f_cl_v+E_m_v+gamma_s+dH2_abs_conv+rho_d` | f_cl_v 0.317->0.55; E_m_v 0.594->0.8; gamma_s 2.22->2.33; dH2_abs_conv 0.26->0.32; rho_d 1e+14->5e+14 | **live** | 1 | 3156.0% | N_void, N_111, d_void, N_100, d_111, d_100 | V11_calib_machine0 |
| `rho_d` | rho_d 1e+13->1e+16 | **live** (drives rows off-grid) | 15 | 1059.9% | N_100, N_111, d_100, N_void, d_void, d_111 | V12_calib_machine0, V1_calib_machine0, V5_calib_machine0, V6_calib_machine0, V7_calib_machine0, V8_calib_machine0, V9_calib_machine0 |
| `E_b_i2` | E_b_i2 0.55->0.875 | **live** (drives rows off-grid) | 6 | 991.2% | N_100, N_void, N_111, d_void, d_111, d_100 | S10_calib, S12_calib, S13_calib, S14_calib, S15_calib_machine0, S8_15dpa, V13_calib_machine0, V15_calib_machine0 |
| `f_cl_v+E_b_v2+s_v` | f_cl_v 0.317->0.7; E_b_v2 0.213->0.35; s_v 1.9->2.5 | **live** (drives rows off-grid) | 2 | 692.8% | d_100, N_void, N_111 | S14_calib |
| `f_cl_v` | f_cl_v 0.05->0.8 | **live** (drives rows off-grid) | 7 | 449.6% | N_void, N_100, N_111, d_100, d_void | S18_calib_machine2, V12_calib_machine0, V1_calib_machine0, V8_calib_machine0 |
| `s_i` | s_i 1->3 | **live** | 12 | 436.5% | d_100, N_111, d_111, N_100 | S17_calib_machine1 |
| `f_cl_i` | f_cl_i 0.02->0.25 | **live** (drives rows off-grid) | 11 | 426.8% | d_100, N_100, N_111, N_void, d_void, d_111 | B4_n100, B5_corner, S15_calib_machine0, S1_search, S2_coarsen, S3_lifetime, S8_15dpa, V10_calib_machine0 |
| `eta` | eta 0.05->0.249 | **live** (drives rows off-grid) | 3 | 397.2% | N_void, N_111, N_100, d_100, d_111 | S15_calib_machine0 |
| `gamma_s` | gamma_s 1.5->3 | **live** (drives rows off-grid) | 20 | 395.8% | d_100, N_100, N_void, N_111, d_void, d_111 | V10_calib_machine0, V1_calib_machine0, V2_calib_machine0, V5_calib_machine0, V6_calib_machine0, V7_calib_machine0, V8_calib_machine0 |
| `r_cap_i_b` | r_cap_i_b 4->8 | **live** | 1 | 192.9% | N_111, N_100, d_100 | S1_search |
| `dH2_abs_conv` | dH2_abs_conv 0.26->0.55 | **live** (drives rows off-grid) | 33 | 173.8% | N_void, N_111, N_100, d_100, d_void, d_111 | M3_absgate, S1_search, S3_lifetime, S4_cavity, S5_bracket, V10_calib_machine0, V13_calib_machine0, V15_calib_machine0, V9_calib_machine0 |
| `phi_max_junc` | phi_max_junc 0.05->0.6 | **live** | 9 | 163.9% | N_100, d_100, N_111, d_111 | B4_n100, S17_calib_machine1 |
| `i_mobile` | i_mobile 5->100 | **live** (drives rows off-grid) | 8 | 132.7% | N_100, N_111, d_100, d_111 | S17_calib_machine1, S1_search, S8_15dpa, V13_calib_machine0 |
| `loop_coal_pref` | loop_coal_pref 3->1e+06 | **live** (drives rows off-grid) | 13 | 112.9% | d_111, N_111, N_100, N_void, d_100 | B3_coal, B5_corner, S2_coarsen, S5_bracket |
| `loop_net_w_c` | loop_net_w_c 5->2.7e+03 | **live** | 32 | 96.0% | N_100, d_100, N_void, N_111, d_void, d_111 | V11_calib_machine0, V13_calib_machine0, V15_calib_machine0 |
| `Z_i` | Z_i 1->1.35 | **live** | 2 | 86.6% | N_111, N_100, d_100, N_void, d_void, d_111 | V5_calib_machine0, V6_calib_machine0 |
| `Z_v` | Z_v 0.9->1.35 | **live** | 3 | 65.2% | N_111, N_100, d_100 | V1_calib_machine0, V5_calib_machine0, V6_calib_machine0 |
| `E_b_hV_1` | E_b_hV_1 2.21->3 | **live** | 5 | 28.7% | N_100, d_100 | S18_calib_machine2, V11_calib_machine0 |
| `s_v` | s_v 2.5->3 | **live** (drives rows off-grid) | 2 | 23.2% | N_100, N_void, d_100 | V12_calib_machine0, V8_calib_machine0 |
| `loop_net_chi` | loop_net_chi 50->1e+03 | **live** | 7 | 22.3% | N_100, d_100 | V11_calib_machine0 |
| `E_a0_conv` | E_a0_conv 1.8->2.3 | **dead** (drives rows off-grid) | 2 | 2.5% | - | S10_calib, S12_calib, S13_calib, S14_calib, S8_15dpa |
| `dH2_conv` | dH2_conv 0.439->0.75 | **dead** | 3 | 0.0% | - | V9_calib_machine0 |
| `LOOP_COAL` | LOOP_COAL 0->1 | **inconclusive** (drives rows off-grid) | 0 | - | - | B3_coal |
| `Z_i+Z_p_i+Z_p_v+Z_gb_i+Z_gb_v` | Z_i 1.2->1.35; Z_p_i 1.2->1.45; Z_p_v 0.93->0.98; Z_gb_i 1.2->1.4; Z_gb_v 0.93->0.98 | **inconclusive** (drives rows off-grid) | 0 | - | - | S12_calib |
| `d_g+rho_p` | d_g 1e-06->4e-06; rho_p 1e+19->1e+20 | **inconclusive** (drives rows off-grid) | 0 | - | - | S13_calib |
| `E_b_v2` | E_b_v2 0.1->0.45 | **dead** | 5 | 0.0% | - | S18_calib_machine2, V1_calib_machine0 |
| `f_cl_v+E_m_v+rho_d` | f_cl_v 0.317->0.55; E_m_v 0.594->0.8; rho_d 1e+14->5e+14 | **inconclusive** | 0 | - | - | V10_calib_machine0 |

### Never varied

These columns exist in the design but have never taken more than one value, so the campaign has no evidence about them:

`E_m_h`, `Z_i_loop`, `Z_loop_model`, `lambda`, `loop_net_K_rec`, `r_cap_v_b`, `r_p`

## Stage history

Only stages that ran to 15.0 dpa are in scope; the rest are listed for provenance but cast no vote on any lever.

| stage | dose | scope | rows | valid | swept | best in-range |
|---|---|---|---|---|---|---|
| `AUTOVERIFY_7100` | 15 | in | 1 | 1 | - | 4 |
| `B1a_V2k` | 15 | in | 8 | 0 | - | - |
| `B1b_V8k` | 15 | in | 8 | 0 | - | - |
| `B2_I80k` | 15 | in | 4 | 0 | - | - |
| `B3_coal` | 15 | in | 7 | 6 | `LOOP_COAL`, `loop_coal_pref` | 4 |
| `B4_n100` | 15 | in | 6 | 6 | `f_cl_i`, `phi_max_junc` | 4 |
| `B5_corner` | 15 | in | 6 | 3 | `f_cl_i`, `loop_coal_pref` | 4 |
| `C1_pref1` | 15 | in | 3 | 3 | - | 2 |
| `F1a_V5k` | 15 | in | 3 | 0 | - | - |
| `F1b_V20k` | 15 | in | 3 | 0 | - | - |
| `F2a_V5k` | 15 | in | 3 | 0 | - | - |
| `F2b_V20k` | 15 | in | 3 | 0 | - | - |
| `F3a_V5k` | 15 | in | 3 | 0 | - | - |
| `F3b_V20k` | 15 | in | 3 | 0 | - | - |
| `F4a_V5k` | 15 | in | 3 | 0 | - | - |
| `F4b_V20k` | 15 | in | 3 | 0 | - | - |
| `G0_legacy` | 15 | in | 4 | 3 | - | 2 |
| `G1_flatD` | 15 | in | 2 | 1 | - | 0 |
| `G1a_vb45_machine0` | 15 | in | 1 | 1 | - | 1 |
| `G1b_V80k_machine0` | 15 | in | 1 | 1 | - | 1 |
| `G1c_V320k_machine0` | 15 | in | 1 | 1 | - | 1 |
| `G2_wolfer` | 15 | in | 4 | 3 | - | 0 |
| `G2a_vb45_machine0` | 15 | in | 1 | 0 | - | - |
| `G2b_V80k_machine0` | 15 | in | 1 | 1 | - | 3 |
| `G3_both` | 15 | in | 2 | 1 | - | 1 |
| `G3a_V20k_machine0` | 15 | in | 5 | 5 | - | 4 |
| `G3b_V80k_machine0` | 15 | in | 5 | 5 | - | 4 |
| `G4_I160k_machine0` | 15 | in | 1 | 1 | - | 4 |
| `G5a_disc_V1k` | 15 | in | 1 | 0 | - | - |
| `G5b_disc_V4k` | 15 | in | 1 | 0 | - | - |
| `G5c_bm_V1k` | 15 | in | 1 | 0 | - | - |
| `G5d_bm_V4k` | 15 | in | 1 | 0 | - | - |
| `L1a_I80k` | 15 | in | 3 | 3 | - | 4 |
| `L1b_I160k` | 15 | in | 3 | 3 | - | 4 |
| `M1_coal_off` | 15 | in | 4 | 3 | - | 2 |
| `M1_coal_on` | 15 | in | 4 | 4 | - | 2 |
| `M2_ic10` | 15 | in | 1 | 1 | - | 1 |
| `M2_ic20` | 15 | in | 1 | 0 | - | - |
| `M2_ic200` | 15 | in | 1 | 1 | - | 0 |
| `M2_ic50` | 15 | in | 1 | 1 | - | 1 |
| `M3_absgate` | 15 | in | 6 | 5 | `dH2_abs_conv` | 2 |
| `M4_nref` | 15 | in | 5 | 4 | `n_ref_conv` | 3 |
| `Q1_short` | 2 | out | 1 | 0 | - | - |
| `R1_repro_machine0` | 15 | in | 3 | 3 | - | 2 |
| `R2_regress` | 15 | in | 3 | 3 | - | 4 |
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
| `S1_search` | 15 | in | 5 | 5 | `f_cl_i`, `dH2_abs_conv`, `i_mobile`, `r_cap_i_b` | 0 |
| `S1_smoke` | 15 | in | 1 | 0 | - | - |
| `S20_calib_machine0` | 15 | in | 13 | 12 | `E_m_i`, `L_hat`, `B_111` | 3 |
| `S2_cap` | 1 | out | 6 | 0 | - | - |
| `S2_coarsen` | 15 | in | 8 | 0 | `f_cl_i`, `L_hat`, `loop_coal_pref` | - |
| `S2_nocap` | 1 | out | 10 | 0 | - | - |
| `S3_lifetime` | 15 | in | 8 | 6 | `f_cl_i`, `dH2_abs_conv`, `n_ref_conv` | 2 |
| `S3_nucleation` | 0.02 | out | 8 | 6 | `f_cl_i`, `d_g`, `rho_p+r_p`, `s_i` | 1 |
| `S4_bias` | 0.02 | out | 9 | 8 | `f_cl_i`, `Z_i`, `Z_i_loop`, `Z_p_i+Z_p_v`, `Z_gb_i+Z_gb_v`, `rho_d` | 0 |
| `S4_cavity` | 15 | in | 6 | 2 | `E_m_v`, `dH2_abs_conv` | 1 |
| `S4_smoke` | 0.002 | out | 2 | 0 | - | - |
| `S5_binmoment` | 0.005 | out | 2 | 2 | - | 1 |
| `S5_bracket` | 15 | in | 2 | 0 | `dH2_abs_conv`, `loop_coal_pref` | - |
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
| `T2_design_v2_machine0` | 15 | in | 154 | 0 | - | - |
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
| Mac.san.rr.com | 396 | 1.38 h | 0.04 h | 5.7 h | 1 |
| MacBook-Pro.local | 686 | 1.02 h | 0.02 h | 9.64 h | 13 |
| Nasr-Workstation | 569 | 3.35 h | 0.06 h | 19.79 h | 164 |

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

## Curated notes

- BIAS AXIS TESTED AND WEAK (B1_bias, 2026-09-03, 16 rows at 15 dpa, V=2000 and V=8000).  Z_i 1.35 -> 1.20 -> 1.05 moves N_void NON-MONOTONICALLY by <1.4x at both ceilings (6.20e20 -> 7.97e20 -> 7.00e20 at V=8000); relaxing Z_i, Z_p_i and Z_gb_i together to 1.05 (0.17x the workbook bias driving force, computed from the code's own k2: disloc 73%, prec 20%, gb 8% of the SIA sink) LOWERS N_void to 4.47e20.  Over the same rows f_cl_v 0.317 -> 0.70 gives x3.2-4.0.  The bias set is NOT the cavity-density lever; the cascade vacancy source is.  These 16 rows are classified INVALID (delta_FP 0.049-0.27, pile_111 = 1.000: I=20000 truncates the loops at 15 dpa) so they cast no vote in the lever table -- which is why this is a note.  Do not re-propose the bias axis on the strength of its `never varied` / `inconclusive` status without reading this first.
- CAVITY SIZE IS A GRID READOUT, MEASURED (B1a_V2k vs B1b_V8k, 2026-09-03).  The same 8 theta vectors at V=2000 and V=8000: mean_n_v scales x4.32-5.52 for a x4 ceiling, d_cavity 2.42-2.51 -> 3.94-4.12 nm, and every row that was 'in band' on size left the band.  kfrac_cav_2_4 collapses 0.95 -> 0.15.  Across all 411 at-dose rows ever run, occ_v is bimodal: sub-critical rows sit at median 0.0039 (grid-independent), super-critical rows at median 0.347 REGARDLESS of theta or V.  The campaign's own tail measurement explains it -- c(n) ~ n^-0.47, and with alpha < 2 the first moment diverges with the cutoff, so <n> is set by the ceiling.  CONSEQUENCE: d_void cannot be calibrated by any parameter in the current physics; tuning it is choosing V.  Score cavities on the windowed fractions (kfrac_cav_2_4, kfrac_cav_1_10) until the model gains a term that bounds growth above r*.  N_void is also grid-sensitive (x1.7-2.2 for x4 in V) -- lever RATIOS transfer, absolutes do not.
- REFERENCE VECTOR IS NOW B3_coal ROW 9305 (workbook default since 2026-09-04, at the author's direction).  11 cells changed in input_parameters.xlsx; previous workbook kept verbatim as input/input_parameters.BACKUP-pre9305-20260904.xlsx (sha 79f700b60e1c6b1e, the sha every campaign run through B5 recorded).  New sha e2d36d163357ed8e.  Changed: f_cl_i 0.09->0.124017, f_cl_v 0.317274->0.55, E_m_v 0.78->0.80, gamma_s 2.21795->2.30, E_b_i2 0.75->0.60, phi_max_junc 0.05->0.10, loop_net_w_c 2.482e-8->1.241e-9 m (100->5 b_111), rho_d 1e14->5e14, loop_coal_pref 30->20, LOOP_COAL stays 1.  WHY: best mean |log10(model/target)| of any full-dose row (0.194 vs 0.65 for the previous ledger best) and the FIRST vector ever to put d_111 in band (4.54 nm), via loop coalescence.  WHAT IT IS NOT: the row is INVALID by classify() -- delta_FP 0.148 vs a 1e-2 tolerance and pile_111 = 1.000 (the 1/2<111> population is against the top of an I=80000 grid), and its d_void 'agreement' is the V=20000 ceiling, not a fit.  It is the best AGREEMENT available, not a conserving row; treat absolute d_void and any 5/6 claim accordingly.  N_100 remains 2.3x below its band and four brackets (bias set, f_cl_i up, f_cl_i down, phi_max_junc) failed to recover it -- see the <100> conversion-channel note.  SIDE EFFECT: E_b_i2 = 0.60 was listed under 'Unaffordable lever values' from S14 (0 of 6 rows reached dose).  Row 9305 reaches 15 dpa at E_b_i2 = 0.60, so that entry is contradicted and should be dropped on the next learn.py pass that re-derives affordability.
- REFERENCE WORKBOOK WAS INCOMPLETE UNTIL 2026-09-04 (second fix, same day).  The row-9305 promotion wrote the 40 SAMPLED design parameters but not the FIXED block that run_ensemble.apply_theta writes at run time 'so a stale workbook cannot leak in' (run_ensemble.py, the `for f in spec['fixed']` loop).  Those keys are in no design CSV, so the promotion never saw them and my verification -- which compared design columns to the workbook -- passed while they were wrong.  Stale/missing: psucc_abs_pref 1 -> 2.0 (this is A_abs, the direct multiplier on the <100> absorption success gate, i.e. the workbook was running <100> conversion at HALF strength), gamma_a_conv 0.03 -> 0.02, loop_net_rho_max 1e16 -> 5e14, grow_boost_100 ABSENT -> 1.0.  n_ref_conv stays 30 (a sampled value beats the fixed default, and the design carries 30).  Workbook sha e2d36d163357ed8e -> b0db9433c84609f3.  SYMPTOM that found it: the notebook/simulation path could not reproduce the campaign -- N_100 came out 50x low at 0.004 dpa and 4 ORDERS low by 6 dpa, and the run froze at 0.4 dpa.  After the fix the same driver reproduces row 9305's ladder (N_111 within 1%, d_111 within 6%, cavities within 15%, N_100 the right order and trend) and runs to 40 dpa without freezing.  LESSON: promoting a calibration row means writing the sampled parameters AND the fixed block; verifying only the design columns is not sufficient.

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

