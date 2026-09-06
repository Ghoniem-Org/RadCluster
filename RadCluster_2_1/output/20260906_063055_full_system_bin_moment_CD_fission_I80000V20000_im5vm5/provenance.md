# RadCluster_2_1 run — 20260906_063055_full_system_bin_moment_CD_fission_I80000V20000_im5vm5

## (1) Material Data

_All input tables from the Excel workbook with parameter overrides applied._

### Production (Fission)

- C_i: 0.1093
- C_v: 0.1506
- G_He/G: 0.5−1
- eta: 0.248719
- f_cl_i: 0.124017
- f_cl_v: 0.65
- i_cascade: 20
- phi_dot: user
- s_i: 1
- s_v: 2.5
- v_cascade: 10

### Production (Fusion)

- C_i: 0.0553
- C_v: 0.1296
- G_He/G: ~10
- eta: 0.248719
- f_cl_i: 0.09
- f_cl_v: 0.317274
- i_cascade: 50
- phi_dot: user
- s_i: 1
- s_v: 2.5
- v_cascade: 20

### Energetics

- B2: 1.67e-29
- B3: 1.84e-58
- E_B_2i: 0.8
- E_B_2v: 0.22
- E_b_hV_1: 2.21315
- E_b_hV_2: 2
- E_f_i: 3.64
- E_f_v: 2
- E_m_2i: 0.42
- E_m_2v: 0.7
- E_m_h: 0.0735044
- E_m_i: 0.401809
- E_m_v: 0.8
- E_s_He: 2.35
- Omega: 1.18e-29
- a: 0.2867
- b_100: 0.2867
- b_111: 0.2482
- gamma_s: 2.3
- mu: 82
- nu: 0.29
- nu_h: 3e+12
- nu_i: 1e+13
- nu_v: 1e+13

### Diffusion

- B_rot: 2.627
- D0_h: 2.5e-08
- D0_i: 8.2e-08
- D0_v: 8.2e-08
- E_b_C_SIA: 0.6
- E_b_C_V: 0.45
- E_b_C_loop: 0.6
- E_b_Cr_SIA: 0.1
- E_b_Cr_V: 0.05
- E_b_Cr_loop: 0.1
- E_b_Mn_SIA: 0.2
- E_b_Mn_V: 0.1
- E_b_N_SIA: 0.6
- E_b_N_V: 0.4
- E_b_N_loop: 0.6
- E_b_W_V: 0.27
- E_m_1D: 0.34
- E_m_h: 0.0735044
- E_m_i: 0.401809
- E_m_v: 0.8
- L_hat: 1016.12
- c_C: 0.0005
- c_Cr: 0.094
- c_Mn: 0.0047
- c_N: 0.0002
- c_W: 0.0033
- i_mobile: 5
- nu0_1D: 6e+12
- s_1D: 1
- s_vc: 1
- v_mobile: 5

### Dissociation

- A_100: 3
- A_111: 3
- A_He_2: 0.55
- A_He_3: 0.4
- A_He_4: 0.75
- A_void_0: 1.2353
- A_void_1: 2.9064
- A_void_2: 3.4147
- A_void_3: 2.1504
- A_void_4: -0.159
- B2: 1.67e-29
- B3: 1.84e-58
- B_100: 0.3581
- B_111: 0.313688
- E_TM_0_4: 0.3
- E_TM_0_5: 0.1
- E_TM_0_6: 0
- E_TM_1_5: 1
- E_TM_1_6: 0.5
- E_TM_1_7: 0
- E_b_i2: 0.6
- E_b_v2: 0.2129
- E_f_v: 2
- alpha_He: 1.7
- gamma_s: 2.3
- gamma_sf: 0.6
- lambda: 0.715709
- lambda_v: 0.575
- m_He: 4.0026
- mu_He: 0.658
- n_tr: 25
- nu0_TM: 1e+12
- r_c: 0.333
- sigma_tr: 5

### Reactions

- A_1D: 2.632
- A_loop: 10.78
- A_sph: 7.818
- B_rot: 2.627
- C_floor: 1e-25
- E_a0_conv: 2.1
- G: 1e-07
- G_He_r: 1
- I: 80000
- I_bin: 18
- K_iv_pf: 21.77
- LOOP_COAL: 1
- LOOP_NETWORK_LOSS: 1
- L_He_max: mf
- T: 603.15
- T_star_conv_C: 450
- V: 20000
- VOID_NETWORK_LOSS: 1
- V_bin: 20
- Z_He: 1
- Z_gb_i: 1.4
- Z_gb_v: 0.93
- Z_i: 1.35
- Z_i_loop: 1.1813
- Z_loop_model: 0
- Z_p_i: 1.45
- Z_p_v: 0.93
- Z_v: 1
- absorb_boost_100: 1
- alpha_He: 1.7
- b0_fission: 0.01
- b0_fusion: 0.1
- ci1_seg: 1.11393e-11
- cv1_seg: 2.02149e-09
- dH2_abs_conv: 0.36
- dH2_conv: 0.439455
- dH_rev_conv: 0.3
- d_g: 1e-06
- gamma_a_conv: 0.02
- grow_boost_100: 1
- he_kinetics: quasi_steady_state
- i_discrete: 100
- loop_coal_pref: 20
- loop_net_K_rec: 1e-06
- loop_net_chi: 50
- loop_net_n_inc: nan
- loop_net_rho_max: 500000000000000
- loop_net_w_c: 1.241e-09
- loop_net_xi: 0
- n1_bin: 1
- n_group: 50
- n_j_min_frac: 0.6
- n_j_min_junc: 30
- n_loop_min: 4
- n_moments: 2
- n_ref_conv: 30
- nu0_conv: 10000000000000
- phi_max_junc: 0.1
- physics_option: bin_moment_CD_fission
- psucc_abs_pref: 2
- r_cap_i_b: 3
- r_cap_v_b: 1
- r_p: 2e-08
- r_ratio: 2
- rho_d: 500000000000000
- rho_net: 5e+14
- rho_p: 100000000000000000000
- shape_function: linear
- sigma_s_junc: 0.35
- solver_mode: full_system
- spectrum: fission
- v_discrete: 5
- void_net_chi: 1165
- void_net_m_inc: nan
- void_net_w_c: nan

### Derived (computed)

- A_1D: 0.524468
- A_loop: 10.7742
- A_sph: 7.79555
- B_rot: 2.627
- Cv_eq: 1.94301e-17
- D1D_base: 1.06694e-09
- Dh_Fe: 5.99502e-08
- Dh_eff: 5.99502e-08
- Di_Fe: 3.60942e-10
- Di_eff: 1.21823e-12
- Dv_Fe: 1.6992e-13
- Dv_eff: 9.47208e-15
- E_f_v: 2
- E_m_h: 0.0735044
- E_m_i: 0.401809
- E_m_v: 0.8
- E_s_He: 2.35
- G: 1e-07
- G_He: 1e-13
- G_He_r: 1
- L_hat: 1016.12
- Omega: 1.18e-29
- T: 603.15
- a_m: 2.867e-10
- b_111: 2.482e-10
- boundary_flux: absorption
- gamma_s: 2.3
- i_mobile: 5
- kBT: 0.0519754
- nu_h: 3e+12
- nu_i: 1e+13
- nu_v: 1e+13
- omega_h_eff: 7.29349e+11
- omega_i_eff: 1.48209e+07
- omega_v_eff: 115236
- r0: 1.41231e-10
- s_1D: 1
- spectrum: fission
- trap_SIA: 295.285
- trap_VAC: 16.9391
- trap_loop: 146.981
- v_mobile: 5

## (2) User Selections

### Run configuration

- C_floor: 1e-25
- G: 1e-07
- I: 80000
- I_bin: 18
- L_He_max: None
- T: 603.15
- V: 20000
- V_bin: 20
- alpha_He: 1.7
- boundary_flux: absorption
- excel_file: /Users/ghoni/Documents/GitHub/RadCluster/RadCluster_2_1/input/input_parameters.xlsx
- he_kinetics: quasi_steady_state
- i_discrete: 100
- i_mobile: 5
- physics_option: bin_moment_CD_fission
- shape_function: linear
- solver_mode: full_system
- spectrum: fission
- v_discrete: 5
- v_mobile: 5

## (3) Solver Configuration

### Solver settings

- atol: 1e-20
- log_time: True
- loop_conversion: 1
- n_points: 40
- rtol: 1e-05
- solver_method.concentration_threshold: 1e-22
- solver_method.linsol: gmres
- solver_method.preconditioner: Woodbury
- solver_method.window_pad: 20
- solver_method.window_width: 10
- t_span: [1e-06, 4e+08]
- timeout_s: 172800

## (4) Run Statistics

### Runtime and machine

- cpu_count: 16
- hostname: MacBook-Pro.local
- n_time_points: 10
- omp_num_threads: 12
- omp_threads_used: 12
- partial: False
- platform: macOS-26.6.2-arm64-arm-64bit
- process_rss_GB: 0.531
- processor: arm
- python: 3.11.7
- ram_available_GB: 84.54
- ram_total_GB: 128
- run_status: completed
- solver.ncfn: 2903
- solver.netf: 428
- solver.nfe: 43450
- solver.nli: 180551
- solver.nli_per_nni: 4.16
- solver.nlsetup: 6977
- solver.nni: 43447
- solver.npe: 3053
- solver.nps: 223199
- solver.steps: 27811
- timestamp: 20260906_063055
- wall_clock_s: 1317.37


---

# ⭐ THIS IS THE REFERENCE RUN  (supersedes 20260904_212202)

_Designated 2026-09-06. The full analysis behind it is tracked at
`digital_twin/CAMPAIGN_REPORT_2026-09-06_void_sweeping.md`; this section is the
run-local copy._

Produced by executing the notebook configuration with **PARAM_OVERRIDES empty** —
the workbook is the sole source of physics. Verified against the override run
that discovered the vector: **0.000e+00 relative difference** on every
observable at every one of the 37 dose points, so the shipped defaults
reproduce it exactly.

## What is new

The first vector in the campaign with **all six observables in band** at
330 °C / 15 dpa, and the first whose **cavity size is not a grid artefact**.

| parameter | value | change |
|---|---|---|
| `VOID_NETWORK_LOSS` | 1 | **new channel** |
| `void_net_chi` | 1165 | new |
| `f_cl_v` (fission) | 0.65 | was 0.55 |
| everything else | calibration row 9305 | unchanged |

| observable | model @15.72 dpa | band | margin |
|---|---|---|---|
| N_100 | 4.715e21 | ≥ 4.67e21 | +0.96% |
| d_100 | 5.939 | 3.4 – 7 | comfortable |
| N_111 | 3.498e21 | ≥ 1.73e21 | comfortable |
| d_111 | 4.637 | 3.4 – 7 | comfortable |
| N_void | 3.609e20 | ≥ 3.6e20 | +0.25% |
| d_void | 2.890 | ≤ 2.9 | +0.34% |

## The result that matters is NOT the 6/6

**Occupancy `mean_n_v/V` = 0.054, against 0.382 for the previous reference.**
Every prior run in this campaign sat near occ ≈ 0.35 (median 0.347 over 411
rows): the cavity distribution piled into the top bin and `d_cavity` was reading
the V = 20000 ceiling rather than physics. With sweeping active the distribution
turns over on its own, so cavity size is a prediction for the first time. That
holds robustly across the whole χ ≥ 700 family; the 6/6 does not.

## Why the existing parameters could not do this

Cavity growth never self-limits: `dm/dt = A_sph·m^(1/3)·drive`, and the drive
`ω_v c_v − ω_i c_i` is a **constant 30% of the vacancy capture rate from 1e-5 to
40 dpa**. So m ∝ t^1.5 unbounded — the model's own 40 dpa answer is ~7.7e7
vacancies (d ≈ 118 nm), and the old reference's 7307 was the ceiling.

- **`E_m_v` cancels.** It scales `k2_vac` (∝ D_v) and the cavity capture
  coefficient `K_v` (∝ ω_v) identically, so it drops out of the growth rate. A
  nucleation knob, not a size knob (N_voids ×41110, d_cavity +7.4% over its box).
- **The bias factors flip the loops.** Loops survive on a NEGATIVE residual,
  `Z_i_loop·ω_i·c_i − ω_v·c_v` = −4.69e-5, only −17% of its vacancy term.
  `Z_v` = 1.15 halves it: measured loop content ×65, d_111 4.17 → 13.9 nm,
  δ_FP 0.13. `Z_v` = 1.30 zeroes it outright.
- **A conservation bound floors the rest.** `S = S_I + ΔJ^d` closes on the old
  reference at S = 7.296e-5, S_I = 4.610e-5, ΔJ^d = 2.686e-5. Holding the loops
  fixed pins S_I, so **d_cavity ≥ 4.70 nm by any parameter choice**.

Sweeping escapes all three: it never touches c_v (so the loop drive survives),
it routes cavity vacancies straight into ΔJ^d (so S may fall below S_I), and
Λ ∝ d_cav ∝ m^(1/3) removes large cavities preferentially — the only
size-dependent negative feedback in the model.

## CAVEATS — read before using this vector

1. **The 6/6 sits on three band edges.** N_100 +0.96%, N_void +0.25%, d_void
   +0.34%. `f_cl_v` 0.65 → 0.66 (a 1.5% change) drops the score to 4/6.
   `learn.py`'s `worst_margin` would rate this ≈0.0025 and demote it — that
   tie-break exists to catch exactly this. The robust sibling is χ = 1150 /
   f_cl_v = 0.65 at 5/6 (d_void 2.904, missing by 0.14%).
2. **χ ≈ 1165 is a GLIDE-rate proxy, not an elastic capture radius.** It implies
   a swept volume ~1000× what network *climb* supplies (climb moves a line
   1.9 µm over the entire 40 dpa). Defensible only as glide-driven sweeping
   under stress, with χ absorbing the mobile-dislocation velocity. **Needs
   physical sign-off before publication.**
3. **χ and f_cl_v are COUPLED.** f_cl_v = 0.65 is stable only because sweeping
   is on: sweeping removes cavities, cutting the vacancy sink, raising c_v, and
   pushing the loop drive back negative. Turn the channel off while keeping
   f_cl_v = 0.65 and the loops run away — f_cl_v = 0.70 alone hung for 3h19m
   unconverged. Do not vary them independently.
4. **This is a 15 dpa calibration, not a dose-independent one.** At the 40 dpa
   endpoint the vector drifts out of band: N_voids 1.94e20 (below the 3.6e20
   floor), mean_n_v 487, N_loops_100 9.01e21 (at the 9e21 ceiling). The targets
   are 15 dpa data and the previous reference drifted too, but the 6/6 headline
   should not be read as holding at all doses.
5. **Conservation degrades with dose:** δ_FP 0.0604 at 15.72 dpa → 0.0696 at
   40 dpa, against the model's own 1e-2 gate. The previous reference was already
   0.057, so this is inherited rather than introduced by the channel — the OFF
   path is byte-identical and the χ = 1 case moved δ_FP by ×1.0013.
6. **Not extent-verified.** No second V extent exists for this vector, so
   `learn.py` will not promote it. That is the next step.

## Reproducing

Notebook `codes/Notebooks/RadCluster_2_1.ipynb` as shipped (grid I = 80000,
V = 20000, i_mobile 5, v_mobile 5, i_discrete 100, I_bin 18, V_bin 20,
rtol 1e-5, t_span (1e-6, 4e8), PARAM_OVERRIDES empty). Figures:
`python codes/make_dose_figures.py <run_dir> --in-place --annotated`.

`plots/plot_data.pkl` (37 MB) is deliberately NOT tracked — it is a regenerable
cache and a pickle of live objects, which will not load after a class refactor.
Everything irreplaceable in this directory is committed.
