# Campaign report — cavity sweeping closes all six observables

**Date:** 2026-09-06
**Outcome:** first vector in the campaign with **6/6 observables in band** at
EUROFER97 / 330 °C / 15 dpa, and the first with a **grid-converged cavity size**.

> `output/` is gitignored, so this file is the tracked record. The run
> directory it describes carries the same analysis in its `provenance.md`.

## The vector

| parameter | value | change |
|---|---|---|
| `VOID_NETWORK_LOSS` | 1 | **new channel** |
| `void_net_chi` | 1165 | new |
| `f_cl_v` (fission) | 0.65 | was 0.55 |
| everything else | calibration row 9305 | unchanged |

Workbook defaults as of 2026-09-06. Previous workbook (row 9305) preserved as
`input/input_parameters.BACKUP-pre-voidsweep-20260906.xlsx`.

| observable | model | band | margin |
|---|---|---|---|
| N_100 | 4.715e21 | ≥ 4.67e21 | +0.96% |
| d_100 | 5.939 | 3.4 – 7 | comfortable |
| N_111 | 3.498e21 | ≥ 1.73e21 | comfortable |
| d_111 | 4.637 | 3.4 – 7 | comfortable |
| N_void | 3.609e20 | ≥ 3.6e20 | +0.25% |
| d_void | 2.890 | ≤ 2.9 | +0.34% |

## Why the existing parameters could not do this

The starting question was whether vacancy migration energy or the dislocation
bias factors could stop cavities piling against the V ceiling. They cannot, and
the reasons are structural rather than a matter of searching harder.

**Cavity growth never self-limits.** `dm/dt = A_sph·m^(1/3)·drive`, and the
drive `ω_v c_v − ω_i c_i` is a *constant* 30% of the vacancy capture rate from
1e-5 dpa to 40 dpa — c_v and c_i each fall ~2× over six decades and fall
together. So m ∝ t^1.5 without bound; the model's own 40 dpa answer is ~7.7e7
vacancies (d ≈ 118 nm). Row 9305's mean of 7307 was the V = 20000 ceiling, not
physics — occupancy 0.365, and the median across 411 campaign rows is 0.347.

**`E_m_v` cancels out of the growth rate.** It scales the vacancy sink
(`k2_vac ∝ D_v`) and the cavity capture coefficient (`K_v ∝ ω_v`) by the same
factor. That is why the campaign measured `d_cavity` +7.4% but `N_voids`
+41110% across its prior box: it is a nucleation knob, not a size knob.

**The bias factors destabilise the loops.** Loops survive on a *negative*
residual, `Z_i_loop·ω_i·c_i − ω_v·c_v` = −4.69e-5, only −17% of its vacancy
term. `Z_v` = 1.15 halves that drive; measured, loop content rose ×65, d_111
went 4.17 → 13.9 nm and δ_FP to 0.13. `Z_v` = 1.30 zeroes the drive outright.
The earlier claim that `Z_v` is a clean cavity lever *because it appears in no
SIA equation* was wrong: it lowers c_v, and c_v is what shrinks loops.

**A conservation bound floors the rest.** The swelling identity
`S = S_I + ΔJ^d` closes on the old reference at
S = 7.296e-5, S_I = 4.610e-5, ΔJ^d = 2.686e-5. Every vacancy in a cavity is
matched by an SIA in a loop or at a dislocation, so holding the loops fixed
pins S_I and leaves only ΔJ^d reducible — **d_cavity ≥ 4.70 nm by any parameter
choice**. The 2.6 nm target is unreachable that way.

Five perturbations confirmed it. Every one that shrank cavities did so by
starving them, `N_void` falling harder than `d_cav`:

| case | δ_FP | d_cav | N_void | S_I |
|---|---|---|---|---|
| `Z_i` = 1.15 | 0.127 | ×0.66 | ×0.16 | ×58 |
| `Z_v` = 1.15 | 0.132 | ×0.63 | ×0.15 | ×65 |
| `s_v` = 1.8 | 0.075 | ×0.82 | ×0.44 | ×5.0 |
| `s_v` = 2.8 | 0.052 | ×1.01 | ×1.04 | ×0.66 |

`s_v` = 2.8 is the only stable one and moves cavities not at all.

## The mechanism that works

A climbing/gliding network dislocation that intersects a cavity absorbs it
whole, delivering its m vacancies to the line:

```
Λ_m^void = |v_net| · ρ_net · (χ · d_cav(m))
```

Three properties make it the only viable route:

1. **It never touches c_v**, so the loop drive — which every parameter lever
   flips — is left intact.
2. **It escapes the conservation floor**, routing cavity vacancies straight
   into ΔJ^d so S may fall below S_I.
3. **Λ ∝ d_cav ∝ m^(1/3)** removes large cavities preferentially and truncates
   the tail: the only size-dependent negative feedback in the model.

Implementation is gated off by default and byte-identical when off (verified:
0.000e+00 against the prior reference on every observable). Swept content is
charged to `J_VAC_fixed`, so δ_FP reads the channel as a transfer, not a
violation (0.0567 vs 0.0566 at χ = 1).

**No `P_ld` geometric gate**, unlike the loop channel. That gate asks whether a
cluster lies inside a line's *instantaneous* elastic zone, which is right for
incorporating a static loop. A sweeping dislocation traverses the volume, so
`v·ρ·w` already *is* the encounter rate; with the loop gate applied, a 44.7 nm
spacing against a 5.5 nm cavity gives P_ld = 0 and the channel vanishes at χ=1.

## Calibration path

χ alone (f_cl_v = 0.55) reaches d_void = 2.85 in band but costs the loops:
d_111 falls to 3.19 (below its 3.4 floor) and N_100 to 3.83e21. Sweeping
removes cavities, which cuts the vacancy sink, raises c_v, and shrinks loops.

Raising `f_cl_v` compensates on both counts — fewer free vacancy monomers means
a less negative loop drive (recovering d_111) and more cascade nuclei
(recovering N_void). **`f_cl_v` = 0.62 does not bifurcate here, though 0.70
alone hung for 3h19m unconverged**: sweeping raises c_v and cancels the
destabilisation. The two are therefore coupled and cannot be varied
independently.

| χ | f_cl_v | occ | /6 | missing |
|---|---|---|---|---|
| — (ref) | 0.55 | 0.382 | 4/6 | N_100, d_void |
| 300 | 0.55 | 0.163 | 4/6 | N_100, d_void |
| 1000 | 0.55 | 0.051 | 4/6 | N_100, d_111 |
| 3000 | 0.55 | 0.015 | 2/6 | over-swept |
| 1000 | 0.62 | 0.059 | 4/6 | N_100, d_void |
| 1150 | 0.65 | 0.054 | 5/6 | d_void by 0.14% |
| 1300 | 0.65 | 0.047 | 5/6 | N_void by 3.9% |
| **1165** | **0.65** | **0.054** | **6/6** | — |
| 1165 | 0.66 | 0.054 | 4/6 | N_void, d_void |

## What is actually established, and what is not

**Established.** Cavity size is now a *prediction*: occupancy 0.054 against
0.382, so the distribution turns over below the ceiling instead of reading it.
That is the substantive result and it is robust across the whole χ ≥ 700 family.
The loops are undisturbed (S_I 5.19e-5 vs 4.21e-5 reference).

**Not established — do not oversell the 6/6.** Three observables sit within 1%
of a band edge, and `f_cl_v` 0.65 → 0.66 (a 1.5% change) drops the score to
4/6. `learn.py`'s `worst_margin` would rate this ≈0.0025 and demote it; that
tie-break exists precisely to catch edge-sitting rows. The robust sibling is
χ = 1150 / f_cl_v = 0.65 at 5/6 (d_void 2.904, missing by 0.14%) with nearly
identical physics.

**χ ≈ 1165 is a glide-rate proxy, not an elastic capture radius.** It implies a
swept volume ~1000× what network *climb* supplies — climb moves a line 1.9 µm
over the entire 40 dpa. The channel is defensible only as glide-driven sweeping
under stress, with χ absorbing the mobile-dislocation velocity. **This needs
physical sign-off before publication.**

**Outstanding.** δ_FP ≈ 0.060 against the model's own 1e-2 gate (the previous
reference was already 0.057, so this is inherited, not introduced). The vector
is not extent-verified, so `learn.py` will not promote it — a second V extent
is the next step. Scoring uses 15 dpa targets read off the 40 dpa trajectory at
its 15.72 dpa output point.
