# RadCluster_2_1 run — 20260906_145255_SMOKE_full_system_bin_moment_CD_fission_I200V100_im5vm5

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
- E_a0_conv: 2.1
- G: 1e-07
- G_He_r: 1
- I: 200
- I_bin: 6
- K_iv_pf: 21.77
- LOOP_COAL: 1
- LOOP_NETWORK_LOSS: 1
- L_He_max: mf
- T: 603.15
- T_star_conv_C: 450
- V: 100
- VOID_NETWORK_LOSS: 1
- V_bin: 6
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
- ci1_seg: 2.26227e-11
- cv1_seg: 4.08502e-09
- dH2_abs_conv: 0.36
- dH2_conv: 0.439455
- dH_rev_conv: 0.3
- d_g: 1e-06
- gamma_a_conv: 0.02
- grow_boost_100: 1
- he_kinetics: quasi_steady_state
- i_discrete: 10
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

- C_floor: 1e-15
- G: 1e-07
- I: 200
- I_bin: 6
- L_He_max: None
- T: 603.15
- V: 100
- V_bin: 6
- alpha_He: 1.7
- boundary_flux: absorption
- excel_file: /Users/ghoni/Documents/GitHub/RadCluster/RadCluster_2_1/input/input_parameters.xlsx
- he_kinetics: quasi_steady_state
- i_discrete: 10
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
- n_points: 10
- rtol: 1e-05
- solver_method.concentration_threshold: 1e-22
- solver_method.linsol: gmres
- solver_method.preconditioner: woodbury
- t_span: [1e-06, 100000]
- timeout_s: 600

## (4) Run Statistics

### Runtime and machine

- cpu_count: 16
- hostname: MacBook-Pro.local
- n_time_points: 10
- omp_num_threads: 2
- omp_threads_used: 2
- partial: False
- platform: macOS-26.6.2-arm64-arm-64bit
- process_rss_GB: 0.109
- processor: arm
- python: 3.11.7
- ram_available_GB: 83.89
- ram_total_GB: 128
- run_status: completed
- solver.ncfn: 4725
- solver.netf: 261
- solver.nfe: 56107
- solver.nli: 68573
- solver.nli_per_nni: 1.22
- solver.nlsetup: 9592
- solver.nni: 56104
- solver.npe: 4836
- solver.nps: 123636
- solver.steps: 34807
- timestamp: 20260906_145255
- wall_clock_s: 6.60016

