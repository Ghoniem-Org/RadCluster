# Campaign M — system size and binning

Production mobility (`i_mobile = v_mobile = 5`), fission cascade, `LOOP_COAL` and `prec_bw` at their defaults in BOTH arms.  All rows share one 45-point log grid to 20 dpa, so the comparison dose is the same point in every column (plan S3.2).

## M1 — cost against domain (fully discrete)

| rung | I = V | N_eq | dose reached | wall |
|---|---|---|---|---|
| M1_D1000 | 1000 | 2006 | 20 / 20 dpa | 5.2 min (done) |
| M1_D2000 | 2000 | 4006 | 0.03839 / 20 dpa | 61.1 min (done) |
| M1_D4000 | 4000 | — | nan / 20 dpa | nan min (running) |
| M1_D8000 | 8000 | — | _not run_ | — |
| M1_D16000 | 16000 | — | _not run_ | — |
| M1_D32000 | 32000 | — | _not run_ | — |
| M1_D10000 | 10000 | — | _not run_ | — |

## M2R — binning convergence at I = V = 1000

`i_discrete = v_discrete = 50` throughout; only the bin count moves, so a deviation below is closure error and nothing else.

| rung | I_bin=V_bin | r | N_eq | N_111 (m^-3) | d_111 (nm) | N_100 (m^-3) | d_100 (nm) | N_cav (m^-3) | d_cav (nm) |
|---|---|---|---|---|---|---|---|---|---|
| **M1_D1000 (exact)** | — | — | 2006 | 3.605e+21 | 3.060 | 5.381e+21 | 5.253 | 1.889e+20 | 1.535 |
| M2R_B2 | 2 | 4.472 | 114 | 3.741e+21 | 2.924 | 5.445e+21 | 5.222 | 1.885e+20 | 1.494 |
| M2R_B3 | 3 | 2.714 | 118 | 3.672e+21 | 2.893 | 5.271e+21 | 5.305 | 1.877e+20 | 1.503 |
| M2R_B4 | 4 | 2.115 | 122 | 3.656e+21 | 2.897 | 5.317e+21 | 4.969 | 1.845e+20 | 1.498 |
| M2R_B6 | 6 | 1.648 | 130 | 3.647e+21 | 2.883 | 5.270e+21 | 5.134 | 1.855e+20 | 1.511 |
| M2R_B8 | 8 | 1.454 | 138 | 3.646e+21 | 2.867 | 5.318e+21 | 5.213 | 1.875e+20 | 1.523 |
| M2R_B12 | 13 | 1.284 | 158 | 3.640e+21 | 2.913 | 5.369e+21 | 5.258 | 1.886e+20 | 1.534 |
| M2R_B16 | 17 | 1.206 | 174 | 3.639e+21 | 2.896 | 5.376e+21 | 5.266 | 1.890e+20 | 1.536 |
| M2R_B24 | 25 | 1.133 | 206 | 3.640e+21 | 2.862 | 5.370e+21 | 5.259 | 1.893e+20 | 1.537 |
| M2R_B32 | 33 | 1.098 | 238 | 3.640e+21 | 2.859 | 5.369e+21 | 5.257 | 1.893e+20 | 1.537 |

### M2R: deviation from the exact arm (%)

| rung | N_eq | r | N_111 | d_111 | N_100 | d_100 | N_cav | d_cav | max |
|---|---|---|---|---|---|---|---|---|---|
| M2R_B2 | 114 | 4.47 | +3.7 | -4.4 | +1.2 | -0.6 | -0.2 | -2.7 | **4.4** |
| M2R_B3 | 118 | 2.71 | +1.8 | -5.4 | -2.1 | +1.0 | -0.6 | -2.1 | **5.4** |
| M2R_B4 | 122 | 2.11 | +1.4 | -5.3 | -1.2 | -5.4 | -2.3 | -2.4 | **5.4** |
| M2R_B6 | 130 | 1.65 | +1.2 | -5.8 | -2.1 | -2.3 | -1.8 | -1.6 | **5.8** |
| M2R_B8 | 138 | 1.45 | +1.1 | -6.3 | -1.2 | -0.8 | -0.7 | -0.8 | **6.3** |
| M2R_B12 | 158 | 1.28 | +1.0 | -4.8 | -0.2 | +0.1 | -0.1 | -0.1 | **4.8** |
| M2R_B16 | 174 | 1.21 | +0.9 | -5.4 | -0.1 | +0.3 | +0.1 | +0.1 | **5.4** |
| M2R_B24 | 206 | 1.13 | +1.0 | -6.4 | -0.2 | +0.1 | +0.2 | +0.1 | **6.4** |
| M2R_B32 | 238 | 1.10 | +1.0 | -6.6 | -0.2 | +0.1 | +0.2 | +0.1 | **6.6** |

## M2 — binning convergence at I = V = 10000

`i_discrete = v_discrete = 50` throughout; only the bin count moves, so a deviation below is closure error and nothing else.

> **No exact arm: the discrete solve at this domain diverges.**  The rungs are listed so their mutual spread can be read, but the deviation columns are WITHHELD.  Scoring a ladder against its own finest rung is how one converges confidently to the wrong answer.

| rung | I_bin=V_bin | r | N_eq | N_111 (m^-3) | d_111 (nm) | N_100 (m^-3) | d_100 (nm) | N_cav (m^-3) | d_cav (nm) |
|---|---|---|---|---|---|---|---|---|---|
| M2_B4 | 4 | 3.761 | 122 | 4.038e+21 | 6.071 | 2.668e+18 | 2.064 | 1.869e+20 | 1.753 |
| M2_B6 | 6 | 2.418 | 130 | 3.852e+21 | 6.878 | 2.703e+18 | 2.042 | 1.913e+20 | 1.807 |
| M2_B8 | 8 | 1.939 | 138 | 3.862e+21 | 6.724 | 1.857e+18 | 1.904 | 2.043e+20 | 1.845 |
| M2_B12 | 13 | 1.555 | 158 | 3.834e+21 | 6.854 | 1.968e+18 | 1.953 | 2.039e+20 | 1.838 |
| M2_B16 | 17 | 1.393 | 174 | 3.829e+21 | 6.884 | 1.990e+18 | 1.960 | 2.039e+20 | 1.840 |
| M2_B24 | 25 | 1.247 | 206 | 3.832e+21 | 6.852 | 1.983e+18 | 1.956 | 2.041e+20 | 1.839 |
| M2_B32 | 33 | 1.180 | 238 | 3.828e+21 | 6.869 | 1.986e+18 | 1.957 | 2.041e+20 | 1.839 |

## What must be stated

1. **The domain is truncated.** V = 10000 sits below the mean cavity size at production dose, so d_cav reads the ceiling at the coarse rungs rather than the model. Both arms see the identical truncation, so it cancels out of the closure error — but nothing here may be compared to experiment.
2. **`delta_FP` runs ~0.05 throughout**, against the model's own 1e-2 gate.  A separate, disclosed defect: closure error and conservation error are different quantities.  Do not conflate them.
3. **The reference is a single trajectory, not a converged domain.**  M1 shows how cost grows with I; it does not show that I = V = 10000 is itself converged in domain.  A rung agreeing with it agrees with the exact solution ON THIS DOMAIN.
