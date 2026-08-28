# EUROFER97 digital-twin calibration — campaign report

**Target:** six observables at 330 °C, 15 dpa, neutron (`targets_330C_15dpa.json`).
**Result: agreement was NOT achieved, and the reason is numerical, not parametric.**
The cavity axis has no converged solution in the usable parameter range, so two of
the six observables are not measurements at any affordable grid.

---

## 1. Headline

**Best VERIFIED result: 3 of 4 measurable observables in band.**

`V500_CONTROL_base_S1700` (stage V5, row 5500), extent-tested at
V = 2000 / 20000 / 200000 -- a **100x** range:

| observable | model | band | in band | extent-verified |
|---|---|---|---|---|
| N_100   | 8.284e21 | [4.67e21, 9e21]   | YES | YES |
| d_100   | 6.785    | [3.4, 7.0]        | YES | YES |
| N_111   | 1.417e22 | [1.73e21, 1.5e22] | YES | YES |
| d_111   | 1.050    | [3.4, 7.0]        | no  | YES |
| N_void  | 2.93e18  | [3.6e20, 3.01e21] | no  | **NO - never verifiable** |
| d_void  | 1.023    | [2.12, 2.9]       | no  | **NO - grid artefact** |

That row is the campaign's ORIGINAL anchor.  Twelve stages moved away from it
chasing cavity numbers that were functions of the grid.

**Leading unverified candidate: V15 row 7100** (`E_b_i2` 0.60, `dH2_abs_conv`
0.36, `loop_net_w_c` 5, `rho_d` 5e14) -- N_100 5.54e21, d_100 4.24,
N_111 9.06e21, all three in band, every value inside its declared prior, at the
`rho_d` where loop observables drift <= 4 % over a 16x grid range.  It is NOT
claimed: it has one extent point, and under the campaign's own rule that makes
it UNVERIFIED until its pair lands.

**Two of the six targets are not measurable at all** (S2, S6): `d_cavity` is
0.2825*(0.37V)^(1/3) to 0.04 %, and `N_void` has never survived a change of grid
extent in 1592 rows.  The campaign is therefore scored out of FOUR.

An earlier "4/6" in this campaign counted band membership on numbers that move
when the discretisation moves.  It is withdrawn.

## 2. The decisive finding: the cavity axis is pinned to the grid

At the best-scoring parameter point, the reported cavity size is reproduced to
**0.04 %** by a formula containing no physics:

```
d_cavity = 0.2825 * (0.37 V)^(1/3)
V =  20000  ->  predicted 5.505   observed 5.503
V =  80000  ->  predicted 8.739   observed 8.738
```

`mean_n_v` is a fixed 36–37 % of the grid ceiling at every V tried, i.e. the void
population sits against the ceiling and the "size" grows as V^(1/3) without limit
(13.9 nm at V=320k, 35 nm at V=5.1M).  There is no converged value to calibrate to.

Three warnings that did NOT fire:
* **`occ_v` is useless here.** One pinned point sat at occ_v = 0.042 and was wrong by 14x.
* **Refinement at fixed extent is not a convergence test.** v_bin 20 -> 47 moved
  `mean_n_v` by only -3.8 % at the anchor; extending V then moved it +316 %.
* **`topbin_v` is geometric, not physical** (see §3).

## 3. `topbin_v` and the `learn.py` fix (committed a92e8cf)

`topbin_v` is the vacancy CONTENT fraction in the top bin.  On a log grid that has
a floor set by the bin RATIO: for a tail c(n) ~ n^-alpha it tends to 1 - r^-(2-alpha).
Measured four times across a 34x span in mean void size and a 4x span in grid extent:

| point | mean_n_v | bin ratio r | topbin_v |
|---|---|---|---|
| G1a | 1566   | 1.193 | 0.1715 |
| G1b | 12411  | 1.196 | 0.1708 |
| G2a | 7115   | 1.193 | 0.1718 |
| G2b | 29599  | 1.196 | 0.1705 |

It tracks r and nothing else.  `TOPBIN_TOL = 0.02` would need ~630 vacancy bins and
is unreachable; it had no probe behind it, unlike every other bar in that file.

`run_ensemble.py` had already withdrawn truncation as a reject criterion
(`TRUNCATION_GATES = False`, 2026-08-05).  `learn.py` reintroduced it, so the ledger
and `CALIBRATION_GUIDE.md` — which `plan.py` plans from — reported a best-of-campaign
that excluded the campaign's best rows.  The fix makes the vacancy top bin a soft
`TOPBINV` flag while `pile` and the SIA-axis top bin keep disqualifying (1095 of 1226
flagged rows fail `pile`, so it rescues 105 rows, not the campaign).

Effect: **100 -> 182 valid rows, best 3/6 -> 4/6, 10 levers inconclusive -> live.**

## 4. What each lever does

Verified on clean one-factor-at-a-time pairs.

### Live and useful
| lever | effect | notes |
|---|---|---|
| `loop_net_w_c` | **clean N_100 knob** | 300 -> 2700 b_111 moves N_100 2.17e22 -> 8.73e20 while d_100 saturates at ~7.26 |
| `dH2_abs_conv` | **d_100 size knob** | 0.28 -> 0.32 moves d_100 -40.8 % with N_100 -6.5 %. ONLY in that window: 0.28 -> 0.55 costs 100 % of N_100 |
| `rho_d` | cavity size, but shared | 1e14 -> 5e14 shrinks voids and simultaneously breaks N_100, d_100, N_111 |

### Structurally inert
| lever | why |
|---|---|
| `loop_net_chi` | `P_ld = 0.5(1+tanh((chi*d - s_ld)/b))` saturates to 1.0 for every chi >= 50 — verified bit-identical at 50/150/400/1000 |
| `rho_d` above 5e14 | `loop_net_rho_max` caps `rho_net` at 5e14 and `rho_net` supersedes `rho_d`; 5e14 and 1e15 are bit-identical |
| `dH2_conv` | junction route runs at P_success = 0.064 vs absorption's 0.683 |
| `E_b_v2` | appears only in `create_excel.py`; never read by the physics |
| junction threshold | at `i_mobile` = 5, `n_j_min_eff` = max(2, min(30, ceil(0.6*5))) = 3, and size-comparability suppresses big-small junctions by 2e-6 |

### No authority over void SIZE (all move void NUMBER only)
| lever | d_void | N_void |
|---|---|---|
| `s_v` 2.50 -> 2.80   | +0.3 % | -18.9 % |
| `f_cl_v` 0.55 -> 0.45| +0.5 % | -10.7 % |
| `E_m_v` 0.80 -> 1.00 | +4.9 % | +48012 % |

These nulls are explained by the pinning of §2, not by kinetics: void size was set by
the grid ceiling in every one of those rows.

### Anomalous — do not trust
`E_b_hV_1` shows a non-monotonic single-point response (2.2131 -> mean_n_100 624.5,
**2.6 -> 528.4**, 3.0 -> 624.5) that reproduces bit-exactly on re-run, with the correct
parameters verified as reaching the workbook.  This looks like solution multistability
or a discrete gate, not a smooth lever.  OFAT gradients near it can report a branch
flip as a lever effect; the ledger currently marks it "live" on that basis.

## 5. Why d_111 cannot be fixed by calibration

d_111 = 1.04 nm (mean_n_111 ~ 18 SIA) against a 6.2 nm target, and it is 1.036 /
1.039 / 1.042 across every grid — the most stable number in the campaign.

* Not a TEM-visibility artifact: the visibility fields show 92 % of ⟨111⟩ loops above
  0.8 nm but only **1 % above 1.5 nm**.  The distribution is monodisperse at n ~ 18
  with no hidden tail.
* Not the junction threshold (see §4).
* ⟨111⟩ loops above `i_mobile` have no sink loss (`k2_SIA` = 0), no effective junction
  loss and negligible network loss — they simply never grow, because ⟨100⟩ loops win
  the mobile-SIA flux competition through the boosted absorption gate.

Levers that shift the split move d_111 by single-digit percent while moving d_100 by
tens of percent.  **6/6 requires a formulation change, not a parameter change.**

## 6. The cavity axis is pinned at EVERY usable rho_d (G3)

The same theta run at two grid extents, `rho_d` swept across its whole effective range
(`rho_net` saturates at `loop_net_rho_max` = 5e14, so 1e15 is bit-identical to 5e14):

| rho_d | mean_n_v V=20k | mean_n_v V=80k | ratio | verdict |
|---|---|---|---|---|
| 1e14 | 7104  | 29599 | 4.16 | PINNED |
| 3e14 | 3661  | 21830 | 5.96 | PINNED |
| 5e14 | 848.3 | 7646  | 9.01 | PINNED |
| 1e15 | 848.3 | 7646  | 9.01 | PINNED |

**There is no converged cavity solution anywhere in the usable range.**  The
rho_d = 5e14 row is the one that mattered: at V=20k it reports d_void = 2.67 nm,
comfortably inside the [2.12, 2.9] band, at occ_v = 0.042.  At V=80k the same theta
reports 5.56 nm.  That in-band cavity size was an artifact.

An earlier draft of this report described a TRADE-OFF — loops agreeing at low `rho_d`,
cavities at high `rho_d`.  That was wrong: the cavity column was never a measurement at
either end.  The correct statement is that **d_cavity and N_void are not measurable by
this model at 15 dpa on any affordable grid**, so the campaign is scored on four
observables, not six:

| measurable observable | status |
|---|---|
| d_100  | in band, grid-stable |
| N_111  | in band, 15 % grid sensitivity |
| N_100  | in band on 2 of 3 grids |
| d_111  | out — structural (§5) |

**Best defensible result: 3 of 4 measurable observables in band.**

## 6b. The SIA axis IS converged (G4) — there is only one defect

Doubling the SIA grid at fixed vacancy grid (I = 80000 -> 160000, i_bin 36 -> 40,
bin ratio held at r ~ 1.203):

| observable | I=80k | I=160k | change |
|---|---|---|---|
| N_100 | 6.882e21 | 6.884e21 | 0.0 % |
| d_100 | 5.900 | 5.901 | 0.0 % |
| N_111 | 1.200e22 | 1.199e22 | 0.0 % |
| d_111 | 1.039 | 1.039 | 0.0 % |
| mean_n_100 | 664.2 | 664.6 | 0.1 % |
| topbin_111 | 0.02703 | 1.0e-6 | -100 % |

`I = 80000` is already converged; the SIA axis is not a source of error.

Two consequences:

1. **N_100's grid sensitivity is not its own.**  N_100 moves 6.88e21 -> 9.50e21 when the
   VACANCY grid is extended, and 0.0 % when the SIA grid is doubled.  So its instability
   is downstream coupling to the unbounded cavity population, not SIA truncation.  The
   campaign has exactly ONE numerical defect, on the vacancy axis, and N_100 inherits it.
   Fixing the cavity channel should stabilise N_100 as a side effect.

2. **The SIA-axis top-bin gate also has a false positive.**  At I=80000 `topbin_111` =
   0.027 trips `TOPBIN_TOL` = 0.02, yet that row is identical to the I=160000 row where
   the flag clears (1e-6).  The gate rejected a converged row.  The §3 fix deliberately
   kept the SIA top bin as a disqualifier; this is one clear counter-example against
   that choice, but one case is not enough to loosen it — it should be re-probed on a
   proper I-ladder before changing.

## 6c. What each in-prior lever does at rho_d = 5e14 (V13/V15)

The three loop targets separate onto three levers, and the separation HELD when
stacked -- the first stacked prediction in this campaign to survive a run:

| lever | sets | effect on the others |
|---|---|---|
| `dH2_abs_conv` 0.34 -> 0.36 | **d_100** 8.72 -> 5.59 | N_100 -15 % |
| `E_b_i2` 0.75 -> 0.60       | **N_111** 1.59e23 -> 1.10e22 (14.5x) | N_100 crashes |
| `loop_net_w_c` 200 -> 5     | **N_100** back up into band | N_111 ~ unchanged |

`E_m_i` 0.40 -> 0.30 is a second route to the first pair (N_100 -79 %,
N_111 -53 %, d_100 unchanged).

**Objective (a) result.**  N_100 only reaches its band when `loop_net_w_c` is
about 5 b_111 or lower -- the loop -> network loss channel must be NEARLY OFF.
The w_c sweep is monotone (5 -> 5.54e21, 15 -> 3.27e21, 40 -> 1.43e21,
70 -> 8.35e20).  The channel exists to saturate loop density, but at the setting
where <100> density matches experiment it is barely active.  The campaign had
been running w_c = 100-2700, most of it outside the declared box [1, 200].

## 6d. The self-learning loop, redesigned (commits a92e8cf, 0a12d54, f1dc0af)

The campaign produced no improvement for twelve stages because it optimised an
objective that was one-third noise, ranked from a sample that excluded its own
best rows, and had no working way to tell the difference.  Four fixes:

| failure | fix |
|---|---|
| scored grid artefacts as results | `learn.py`: an observable is scored only if it survives a change of grid EXTENT (`EXTENT_TOL` = 10 %); rows rank on `n_in_range_verified` FIRST, so an unverified row can never displace a verified one |
| never generated its own verification | `verify.py`: finds the unverified rows the ledger rests on and emits/launches the same theta at 4x extent with the bin ratio held |
| ledger discarded its own best rows | the vacancy top bin became a soft `TOPBINV` flag (it is geometric: 0.1705/0.1708/0.1715/0.1718 across a 34x span in mean void size) |
| designs silently left their priors | `run_ensemble.py`: prior-box audit at the design chokepoint, stamped into provenance as `prior_violations` |

Pairing is on `theta_hash`, which is built from design columns only and is
therefore grid-independent, so pairs form ACROSS stage files: every grid ladder
retroactively certifies rows everywhere in the ledger.  `plan.py` inherits the
anchor fix for free -- it reads `led["best"]`, which now ranks by verified count.

THE INVARIANT for any verification row: `theta_hash` includes `theta_id`, so the
row must copy its parent verbatim except `row_id` and `cond_row_id`.  Renumber
`theta_id` and the pair never forms and the CPU is wasted silently.  (The
V12/V13/V15 designs renumber it and so cannot pair with their parents.)

Coverage after the change -- the diagnosis in one table:

| observable | rows ever extent-verified (of 1592) |
|---|---|
| d_100, d_111 | 30 |
| N_100, N_111 | 26 |
| d_cavity | 7 |
| **N_void** | **0** |

## 6e. The void campaign (W1-W7): density is reachable, size is not

Started from V500 (the verified 3/4 loop row) on the author's objective: land
cavity density AND size without losing the loops.

**The observable was wrong before the physics was.**  `d_cavity` derives from
`mean_n_v`, a number-weighted mean over a distribution the author identifies as
BIMODAL -- small He-stabilised bubbles plus a grown-void population.  A mean
across two modes reports the valley between them.  Measured: `mean_n_v` drifts
314-56000 % with grid extent while `N_voids`, a pure count, drifts 13-108 %, and
on one row the legacy d_cavity read 0.661 nm while a real population sat at
2.09 nm.  Size-resolved fields were added (S3 of this section).

**What each lever does to the cavity channel**

| lever | density | size |
|---|---|---|
| `f_cl_v` 0.32->0.70, `s_v` 2.5->1.8, `eta` 0.25->0.35 | +21 to +71 % | ZERO (d frozen 0.56-0.57) |
| `Z_v` 1.00->0.95 | zero (bit-identical) | zero -- extra vacancies go to RECOMBINATION |
| `rho_d` 1e14->2.5e14 | zero | zero, and costs the loops |
| `E_b_hV_1` 2.21->2.80 | zero (bit-identical) | zero -- not in the fission He path |
| **`E_m_v` 0.594->0.80** | **1e12 -> 1.1e21** | rises WITH density, cannot be separated |

Cascade supply was never the constraint: at 15 dpa with eta=0.25, f_cl_v=0.32 the
cascades make ~1e29 clustered vacancies/m^3 against a 1.5e21 target -- five
orders of magnitude of surplus.  SURVIVAL is the constraint.

**Result.**  Cavity DENSITY reaches the band (N_cav_tem = 1.11e21 at E_m_v 0.80
against a 1.5e21 target).  Cavity SIZE does not, and cannot: `E_m_v` is the only
lever with authority and it raises count and size together.

**And the cavity-active region is not converged.**  Three independent extent
pairs at E_m_v 0.75-0.78: N drifts 95-143 %, d drifts 55 %, all PINNED -- and
each LOSES N_111 at the larger grid, i.e. the unbounded cavity population drags
the loops through the vacancy balance.  So:

| region | cavities | loops |
|---|---|---|
| E_m_v = 0.594 (V500) | BOUNDED (ratio 1.00) but negligible, N_bub ~ 1e12 | 3/4 in band, 2.2 % over 100x in V |
| E_m_v >= 0.70 | present but UNBOUNDED | contaminated by the coupling |

**There is no setting with converged loops AND a converged, correct cavity
population.**

**Why -- and it is structural.**  The fission (Case 2) reduction assigns He as
`ell_bar * m^(2/3)` from ONE scalar `ell_bar = Q_tot / C_vac_tot`
(rate_kernels.cpp:258-273).  Two cavities of the same size therefore carry the
same He and the same binding, so a He-rich stable branch cannot coexist with a
He-poor growing branch.  That is precisely the bimodality the measurement shows,
and it is unreachable by any parameter.  It also explains the inert lever:
`E_b_hV_1` never enters that path (it uses `delta_He`/`beta_He_exp`).

The He SUPPLY is adequate -- 1 appm/dpa x 15 dpa = 1.27e24 He/m^3 is 847 He per
cavity at the experimental density, He/V ~ 1.09, right where an equilibrium
bubble sits.  The model cannot concentrate it into bubbles, not lacks it.

**Three retractions, all the same error.**  "Bounded at V=5000", "the bubble size
is in band at 2.85 nm", and "the shape is invariant to four digits" were each
drawn from rows sharing ONE grid and each overturned by the extent pair.  A
bounded window's diameter is a property of the window (d_cav_bub = 2.854-2.856
across a 5e7 range in count); a tail-inclusive diameter is grid-contaminated
(3.40 -> 5.29).  No cavity claim without its extent pair, however clean the
trend.

## 7. What would actually close this

### 7.1 The cavity channel does NOT need a size cap (G5 + author, 2026-08-20)

An earlier draft recommended adding a size cap to the cavity formulation.  That
was wrong on two counts.

**G5 -- the closure is not the problem.**  Same theta, same I, run in BOTH
solvers at two extents:

| solver | mean_n_v V=1000 -> V=4000 | ratio |
|---|---|---|
| discrete   | 13.54 -> 67.26 | 4.97 PINNED |
| bin_moment | 12.83 -> 65.04 | 5.07 PINNED |

The bin-moment closure reproduces discrete to 5.2 % on the cavity mean and 0.4 %
on the loops (d_100 11.539 vs 11.496; N_100 2.291e22 vs 2.26e22).  It is
faithful.  The unbounded growth is in the PHYSICS AT THESE PARAMETERS, not in
the discretisation -- both solvers do it.

**And the physics does not want a cap.**  Cavities in EUROFER97 at 330 C are
1-3 nm; at the target values they hold ~1.2e24 m^-3 of vacancies against
~4.9e24 m^-3 of SIA content in the loops -- a factor of 4 SMALLER.  There is no
inventory that needs saturating.  A model producing 8 nm cavities at 15 dpa has
the GROWTH BALANCE wrong; bolting on a cap would hide that, not fix it.

So the correct statement is not "the cavity channel lacks a mechanism" but
**"the cavity growth parameters are in the wrong regime, and the campaign never
searched the regime that matters."**

### 7.2 The untested lever: the dislocation bias ratio

Void growth is driven by the dislocation bias -- dislocations preferentially
absorb SIAs, leaving the vacancy excess that feeds cavities.  The campaign ran
`Z_i` = 1.35 against `Z_v` = 1.0, a bias of **1.35**, and never swept it with an
extent check.  Inside the declared boxes `Z_i` reaches 1.02 and `Z_v` reaches
1.05 -- a bias of **0.97**, which should not grow voids at all.

`E_m_v` and `rho_d` were tested (both fail to bound growth); the BIAS RATIO,
the strongest of the three, was not.  That is the gap.

**Lead:** row 382 (T3, V=5000) reports d_void = 2.35 nm against a ceiling
prediction of 3.47 -- **32 % BELOW** the formula, the only row in 1614 that does
not track it.  That is what genuinely self-limited growth looks like, and it is
the place to start.

### 7.3 The rest

1. **d_111 needs a formulation change.**  1.04-1.63 nm against 6.2, stable to
   three digits across every grid and lever.  <111> loops above `i_mobile` have
   no sink loss, no effective junction loss and negligible network loss -- they
   never grow because <100> takes the mobile-SIA flux through the boosted
   absorption gate.  Best single lever: +89 % against a +507 % requirement.
2. **Re-derive every lever sensitivity from the corrected ledger.**  The
   pre-2026-08-20 verdicts were computed partly against the two artefact
   observables and partly across spans outside the prior boxes.
3. **`loop_net_chi` is inert by construction** (`P_ld` saturates to 1.0 for
   every chi >= 50) and its box is [1, 60].  Objective (a) is carried entirely
   by `loop_net_w_c`, and N_100 only lands in band when w_c <~ 5 b_111 -- i.e.
   with the loop -> network channel NEARLY OFF.
4. **Do not generalise a grid property from one theta.**  Three times in this
   campaign a convergence claim measured at a single parameter point failed to
   transfer (v_bin refinement at the anchor; "loops are grid-exact"; "loops
   decouple at rho_d = 5e14" -- row 7100 drifts 16-49 % at that same rho_d).
   Extent-verify the theta you intend to report, not a neighbour.

## 8. Provenance

Stages this session: V11 (loop-network channel), V12 (cavity levers, stopped early),
G1/G2 (vacancy grid ladders), G3 (`rho_d` convergence scan), G4 (SIA grid ladder),
R1 (reproducibility replicate — all rows reproduce bit-exactly; the solver is
deterministic and the harness applies parameters correctly).

---

# Session 2026-08-21 — the ⟨111⟩ comparator, and two silent traps

## 9. `d_111` is a model deficiency, not a comparator artifact

§2 showed `d_cavity_nm` was a grid artifact because it was built from a
whole-distribution `mean_n_v`. `d_111_nm` and `d_100_nm` are built the same way
(`run_ensemble.observe`), so the same challenge applied to them. The visible
DENSITIES (`N_111_vis_*`) had been emitted since the f_100 work; the visible
DIAMETER never was. It is now.

**The cutoff was set from the data, not chosen.** `docs/Database/MicroData.xlsx`
carries measured loop size distributions. The 300 C / 15 dpa entries — the closest
condition to the 330 C / 15 dpa target — have a lower bin edge of **1 nm** and a
measured mean of **3.43 / 4.68 / 4.76 nm** across three areas. So d_min = 1.0.
Note this is NOT the 1.25 nm that happened to make `N_111` fit best; choosing the
window by which answer it produced would have been fitting the comparator.

**L1: rows 7400/7402/7403 at I = 80 000 and 160 000. Drift 0.0 % at every entry.**

| quantity | unwindowed | @0.8 | **@1.0** | @1.25 | @1.5 | band |
|---|---|---|---|---|---|---|
| `d_111` | 1.049 | 1.058 | **1.160** | 1.363 | 1.602 | [3.4, 7.0] ✗ |
| `d_100` | 6.683 | — | 6.693 | — | — | ✓ |

**The window does not rescue `d_111`.** The hypothesis was mine and it is refuted.
The ⟨111⟩ distribution is a dense pile sitting just ABOVE the cutoff — not a few
visible loops buried under invisible ones — so windowing barely moves the mean.
This is the opposite of the `d_cavity` case, and the reason it had to be measured
rather than argued.

**What the window does change is the density.** `N_111_vis_1` = 7.79e21,
extent-converged, and inside the ORIGINAL 1.1e22 ceiling. The 2026-08-19 raise to
1.5e22 was compensating for the wrong comparator and can be withdrawn.

### 9.1 The defect is a character PARTITION error, not a supply error

| | N (m⁻³) | d (nm) | SIA/loop | SIA content (m⁻³) |
|---|---|---|---|---|
| measured ⟨111⟩ | 1.93e21 | 6.20 | 636 | 1.23e24 |
| **model ⟨111⟩** | 1.43e22 | 1.05 | **18** | **2.60e23** |
| measured ⟨100⟩ | 4.97e21 | 6.20 | 735 | 3.65e24 |
| model ⟨100⟩ | 8.62e21 | 6.68 | 854 | 7.36e24 |

Total SIA content is 1.56× measured — the right order — and ⟨100⟩ is essentially
correct. The defect is confined to ⟨111⟩: **7.4× too many loops, each 35× too
small, holding 4.7× too little SIA content.** Interstitials that should build
⟨111⟩ loops end up in ⟨100⟩. This points at the ⟨111⟩ growth/dissociation balance
(`B_111`, `E_b_i2`, the derived `A_111`) — not at helium, and not at the grid.

## 10. Two silent traps found while setting up the Case 1 test

Both would have produced confident, wrong answers.

**(a) The spectrum sheet is selected at build time, the design writes to one sheet.**
`parameters_S4.json` pins `eta`, `f_cl_i`, `f_cl_v`, `s_i`, `s_v` to
`production_fission`, but `bin_moment_rates.py:711` and `rate_equations.py:155`
read `production_fusion` when the cascade is fusion. Any Case 1 run would have
**silently discarded the entire calibrated production vector**, substituted
workbook defaults, and still reported the design row in provenance. Fixed by
mirroring design writes to both sheets (a no-op for fission runs).

**(b) `cascade` switches two things at once.** It selects the He coupling case
AND the cascade source spectrum, and those sheets differ in four parameters no
design column controls:

| | fission | fusion |
|---|---|---|
| `i_cascade` (max SIA cluster) | 20 | 50 |
| `v_cascade` (max vacancy cluster) | 10 | 20 |
| `C_i`, `C_v` | 0.1093 / 0.1506 | 0.0553 / 0.1296 |

`v_cascade` is a direct void-nucleation lever. F1/F2 (Case 1, un-mirrored) returned
cavity densities 100–1000× over band and loops 25× over — attributing that to the
helium formulation would have been wrong. **F1/F2 are void as evidence about
helium and are not cited as such.** Fixed with a `mirror_production` condition flag,
verified on all four parameters. `G_He_r` is now also condition-settable, so the He
SUPPLY and the He FORMULATION can finally be moved independently — previously
`cascade` was the only way to change the supply, which made Case 1 vs Case 2
unfalsifiable in principle.

## 11. Status of the Case 1 question — OPEN

The corrected test (F3/F4: Case 1 with the fission cascade source held, at 1 and
10 appm/dpa) was launched and **stopped by the author before any rows completed**.
Whether a per-size `Q_m` bounds the cavity population is therefore **unanswered**.

What is established about Case 1 is structural, from the code alone: it carries
`Q_m(V)` and computes `ell_bar_m(m) = Q_m[m]/c_v[m]` per size class
(`rate_kernels.cpp:751-795`), so `ell(m)` is free rather than assigned
`ell_bar·m^(2/3)` from one scalar. That is the structure a critical radius needs,
and hence the structure bubble/void coexistence needs. It remains untested here.

Caveat for whoever runs it: in `bin_moment` mode Case 1 stores Q **per bin**, not
per size (`rate_kernels.cpp:2176`). Still size-resolved, but coarser than the
discrete branch.

### Helium supply is not the obstacle

At the workbook rate (`G_He_r` = 1 appm/dpa, Reactions sheet — authoritative, and
NOT the 0.75 spectrum default) 15 dpa gives **1.27e24 He/m³**. Holding the measured
cavity population as equilibrium bubbles needs 5.9e23 (He/V = 0.5) to 1.17e24
(He/V = 1.0) at the central target of 1.5e21 at 2.6 nm. The supply sits inside that
range. The measured BOR-60 population is quantitatively consistent with being
He-stabilized bubbles at the helium this model already generates; what the Case 2
reduction cannot do is LOCALIZE that helium.

## 12. Corrections issued this session

1. **"Windowing will close `d_111`"** — refuted by L1's own measurement (1.16 nm
   against a 3.4 nm floor). Withdrawn.
2. **"The model has no coordinate to represent bubbles vs voids"** — true of
   Case 2, false of the code. Case 1 has carried `Q_m` all along.
3. **He supply "9.55e23 He/m³"** — used the 0.75 spectrum default; the workbook
   sets `G_He_r` = 1, giving 1.27e24. The corrected figure strengthens the bubble
   hypothesis rather than weakening it.
4. **"Cavity density achieved"** (§6e) — too generous. Density enters band only at
   `E_m_v` ≥ 0.78, where `N_100` is 78–169 % over its ceiling. Reachable in
   isolation, never jointly.

## 13. Operational note

Killing a run's parent leaves its workers orphaned. 90 such workers from the
killed W-series were still resident (~36 GB, oldest > 24 h). **Kill children
first, then the parent.** Also: `pkill -f <pattern>` matches any process whose
command line contains the pattern — including a monitoring loop that names the
same file. It killed one this session.

## 14. Case 1 ANSWERED — per-size helium does not bound the cavities

§11 left this open. It is now closed, and the answer is negative.

**Correction to §10(b) first.** F3/F4 (mirrored) returned BIT-IDENTICAL results to
F1/F2 (un-mirrored). The reason: `bin_moment_rates.py:711` selects the production
sheet from `input_data.derived['spectrum']`, which is read from the workbook
Reactions sheet and is always `'fission'` -- NOT from `cascade`. The fusion
production sheet is never consulted, so `i_cascade`/`v_cascade`/`C_i`/`C_v` never
differed and there was no cascade-source confound.

Consequently **`cascade` changes only `he_mode`, and F1/F2 were already a clean
Case 1 test.** Declaring them "void as evidence about helium" in commit 23f746c
was wrong: the 100-1000x cavity densities were caused by Case 1 itself, not by
`v_cascade`. Having just found two genuine silent traps, a third was read into a
result that was in fact a physical finding. The `mirror_production` guard remains
valid for the case where a workbook sets `spectrum = fusion`, but it was inert
here and explains nothing about F1/F2.

**The measurement.** Same theta, same He supply (1 appm/dpa), only `he_mode`
differing:

| E_m_v | mode | N_voids | d_cav V5k → V20k | drift | N_100 | N_111 |
|---|---|---|---|---|---|---|
| 0.594 | Case 2 | 2.46e18 | 0.559 → 0.559 | **0 %** | 8.47e21 ✓ | 1.43e22 |
| 0.594 | Case 1 | 7.50e18 | 0.558 → 0.618 | 11 % | 7.14e21 ✓ | 1.90e22 ✗ |
| 0.700 | Case 2 | 2.11e19 | 1.581 → 3.291 | 108 % | 8.54e21 ✓ | 1.43e22 |
| 0.700 | Case 1 | 2.67e23 | 3.595 → 5.660 | 57 % | 4.19e22 ✗ | 1.35e23 ✗ |

Zero of six rows pass the criterion fixed in advance (`d_cav` in [2.12, 2.9] with
< 10 % extent drift and loops in band). Case 1 triples cavity nucleation at the
loop-best point while keeping cavities sub-nm and INTRODUCING grid drift where
Case 2 had none; above E_m_v = 0.70 it over-produces cavities 100x and loops 10x.
At 10 appm/dpa (F4) it is worse again: N_v ~ 3e24, N_111 ~ 5e24.

**Conclusion.** Giving helium a size coordinate is NOT sufficient for bubble/void
coexistence here. `ell(m) = Q_m/c_v[m]` is still a per-class MEAN, so cavities of
equal size share one gas loading and the saddle-node structure a conversion event
requires -- gas-rich stable below a critical radius, gas-poor growing above -- has
no room to form. The bimodal population the author describes needs the joint
(m, ell) state, not a size-resolved mean of it. Note also that the fixed criterion
was applied as written: Case 1 was not credited for being "closer" on any single
axis when it failed the conjunction.

Cost note: F3/F4 was CPU spent chasing a confound that did not exist. The check
that would have prevented it -- trace which sheet is ACTUALLY read, rather than
which one the switch appears to select -- takes one smoke test, and was run only
after the results came back identical.

## 15. ⟨111⟩ formulation review — the mechanism, and a fix that already exists

§9.1 localized the defect: 7.4× too many ⟨111⟩ loops at 1/35 the size, while
⟨100⟩ and total SIA content are right. This section identifies why.

### 15.1 The binding law is NOT the constraint — my nominated target was wrong

I proposed reviewing `B_111` / `E_b_i2` / `A_111`. At the loop-best theta, for a
⟨111⟩ loop at n = 18:

| channel | rate coefficient |
|---|---|
| capture `K_SIA_grow` | 1.269e+09 |
| thermal emission `G_SIA` | 4.830e-05 |

Thirteen orders of magnitude. Emission is irrelevant above n ≈ 10, so no setting
of the binding law can change `d_111`. **The dissociation balance is not the
mechanism and should not be pursued.**

### 15.2 The real mechanism: ⟨111⟩ loops cannot coarsen

`reaction_rates.py:507-513`:

```
if n > i_mobile:  D = 0            # sessile - no coalescence as projectile
elif n < 4:       D = Di_cluster_3D(n)
else:             D = D1D(n) / rot_factor
```

with `rot_factor = 1 + B_rot·L̂²`. Measured on the campaign theta:

| n | 1–3 | 4 | 5 | ≥ 6 |
|---|---|---|---|---|
| `D_SIA_eff` | 1.218e-12 | 6.645e-19 | 5.316e-19 | **0** |

D falls 1.8e6× between n = 3 and n = 4 and is exactly zero above `i_mobile` = 5.
So only n ≤ 3 genuinely moves, and **loop–loop coalescence is dead**: it needs a
mobile projectile, and no loop qualifies. A ⟨111⟩ loop can grow only by capturing
monomers one at a time.

Note `rot_factor` = 2.712e6 here, against the "≈ 6568" the source comment assumes
— because it scales as L̂² and the campaign runs `L_hat` = 1016, not 50. `L_hat`
is inside its prior [3.5, 3500], but the quadratic consequence for capture was
plainly not anticipated when that comment was written.

### 15.3 The flux competition is 85:15 against ⟨111⟩

Capture scales with loop radius, so the larger ⟨100⟩ population wins the mobile
SIA flux even though it is 40 % fewer in number:

| | per-loop coefficient | × gate | population-weighted share |
|---|---|---|---|
| ⟨111⟩ at n = 18 | 1.269e+09 | — | **15.0 %** |
| ⟨100⟩ at n = 854 | 8.742e+09 (6.89×) | ×1.367 → 9.41× | **85.0 %** |

This is rich-get-richer: bigger ⟨100⟩ captures more, so it gets bigger. And the
graph is a one-way funnel — three ⟨111⟩ → ⟨100⟩ edges (Marian junction, Dudarev
unary, absorption-as-partner) and **no back-conversion**. ⟨100⟩ is an absorbing
state that feeds on ⟨111⟩.

Two smells in that gate: `conv_psuccess_abs` = 1.367 is a *success probability
above unity* (the kernel runs at 137 % of the collision rate), and it sits there
because `dH2_abs_conv` = 0.26 is pinned at its prior FLOOR — a bound the spec
itself records as numerical, not physical, chosen to stop ⟨100⟩ running away.

Arithmetic on what a flux-share fix could buy: n must go 18 → 636 (35×). Giving
⟨111⟩ 100 % of the flux is only 6.7×. **Rebalancing the competition cannot close
`d_111`** — consistent with §7's measured +89 % against a +507 % requirement.

### 15.4 The fix already exists and has never been switched on

`reaction_rates.py:516-541` documents this exact defect (2026-08-16) and
implements the remedy:

> `D_SIA_eff` is hard-truncated to zero at `i_mobile`, so a loop that grows past
> the cutoff can never take part in a coalescence again... **That is why the
> ½⟨111⟩ mean size is pinned AT the cutoff**: `mean_n_111/i_mobile` = 1.32 ± 0.75
> at i_mobile = 40 and 0.84 ± 0.25 at i_mobile = 100 across 16 runs, i.e. the mean
> tracks the cutoff rather than the physics. Real ½⟨111⟩ loops in bcc Fe stay
> glissile to large size and coarsen by 1-D glide.

`D_loop_coal` continues the 1-D glide law past the cutoff, as a separate array so
large loops do not simultaneously become visible to the fixed-sink and cavity
channels. It is gated by `LOOP_COAL`, **default 0**.

`LOOP_COAL` appears in **none of the 79 designs, not in `parameters_S4.json`
(neither `parameters` nor `fixed`), and not in the workbook Reactions sheet.**
`re.get('LOOP_COAL', 0)` has therefore been 0 in every run of this campaign. The
mechanism the report has spent three sessions calling "structural" has a switch,
and the switch has never been flipped.

This also explains the shape of the failure rather than just its size: a mean that
tracks a numerical cutoff is exactly what §9 measured — a dense pile of loops just
above the TEM cutoff, which is why windowing could not rescue `d_111`.

### 15.5 Recommended next test

Enable `LOOP_COAL = 1` at the loop-best theta and re-measure `d_111` / `N_111`,
at two grid extents as usual. `loop_coal_pref` gives the capture-geometry
prefactor the 1-D-glide isotropic-equivalent estimate understates (the same
correction `absorb_boost_100` was created for on the ⟨100⟩ side, and likewise
left at 1.0).

Predictions, stated before the run so they can fail: ⟨111⟩ loop NUMBER falls and
mean SIZE rises at roughly fixed SIA content; `N_111` moves toward 1.93e21 from
1.43e22; `d_111` rises from 1.16 nm. Whether it reaches 3.4 nm is exactly the
open question. Risk to watch: coarsening frees mobile SIAs that ⟨100⟩ may absorb,
so `N_100`/`d_100` — currently the model's best-fitting observables — must be
re-checked, not assumed.

Both `LOOP_COAL` and `loop_coal_pref` should be added to the spec so the
dual-declaration audit (§ commit 4181488) can see them.
