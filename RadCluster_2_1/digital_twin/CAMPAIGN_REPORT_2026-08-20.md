# EUROFER97 digital-twin calibration — campaign report

**Target:** six observables at 330 °C, 15 dpa, neutron (`targets_330C_15dpa.json`).
**Result: agreement was NOT achieved, and the reason is numerical, not parametric.**
The cavity axis has no converged solution in the usable parameter range, so two of
the six observables are not measurements at any affordable grid.

---

## 1. Headline

| observable | band | best verified value | status |
|---|---|---|---|
| d_100  | [3.4, 7.0] nm            | 5.55 / 5.90 / 5.71 across 3 grids | **in band, grid-stable** |
| N_111  | [1.73e21, 1.5e22] m^-3   | 1.20–1.38e22                      | in band, moves 15 % with grid |
| N_100  | [4.67e21, 9e21] m^-3     | 6.88–9.50e21                      | in band on 2 of 3 grids |
| N_void | [3.6e20, 3.01e21] m^-3   | 4.93–7.53e20                      | **NOT MEASURABLE — unconverged axis** |
| d_111  | [3.4, 7.0] nm            | 1.04 (three digits, every grid)   | **out — structural** |
| d_cavity | [2.12, 2.9] nm         | 5.43 → 8.74 as V grows            | **NOT MEASURABLE — grid artifact** |

An earlier "4/6" claim in this campaign counted band membership on numbers that
move when the discretization moves.  It is withdrawn.

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

## 7. What would actually close this

1. **A size cap on the cavity channel.** The `dH2_abs_conv` spec note already
   anticipates the analogue for ⟨100⟩ loops ("restore the lower bound only after a
   size cap exists").  The same pathology is on the vacancy axis and needs the same fix.
2. **Test whether the pinning is the bin-moment closure or the physics** by running a
   reduced-V case in `discrete` mode against `bin_moment`.  Not yet done.
3. **`B_111` at rho_d >= 5e14** — zero clean pairs exist in 1565 rows.  It is the
   ledger's strongest lever on N_111 (5.15x) and N_111 is the binding constraint in the
   converged-cavity region, but its authority there is unmeasured.
4. **Replace `occ_v`/`topbin_v` with a real convergence test**: same theta at two grid
   EXTENTS, reject when `mean_n_v` scales with V.

## 8. Provenance

Stages this session: V11 (loop-network channel), V12 (cavity levers, stopped early),
G1/G2 (vacancy grid ladders), G3 (`rho_d` convergence scan), G4 (SIA grid ladder),
R1 (reproducibility replicate — all rows reproduce bit-exactly; the solver is
deterministic and the harness applies parameters correctly).
