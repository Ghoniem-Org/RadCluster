# Stock-Market Dynamics

**Graph-based cluster dynamics applied to S&P 500 market regimes.**

Stock-Market Dynamics is a standalone prototype that ports the Reaction Admissibility
Graph (RAG) cluster-dynamics method of
[RadCluster](https://github.com/Ghoniem/RadCluster) to an equity-market problem.
Instead of defect clusters binned by size, index constituents are binned into
**market-regime states** — momentum × volatility, extended in later revisions to
include turnover and size — and the capital distribution across those states evolves
under index entry (source), measured transitions and a herding flux (edges), and index
exit (sink).

The model forecasts the **capital distribution**, not returns. That distinction is the
central finding of the project and is stated plainly in Section 9: seven
pre-registered return-forecast tests were run and all seven failed walk-forward
validation. What survived is a regime census with calibrated uncertainty, and a
systematic quality/value/low-beta stock screen.

- **Branch** — this is an **orphan branch** of the RadCluster repository. It shares no
  history with `main` and contains none of the RadCluster source; checking it out
  replaces the working tree entirely.
- **Languages** — Python (`numpy`, `scipy`, `pandas`, `matplotlib`).
- **License** — MIT.
- **Method reference** — N. M. Ghoniem, *A Generalized Graph-Based Cluster Dynamics
  Framework for Irradiated Materials*, Journal of Nuclear Materials (submitted, 2026).
- **Project documentation** —
  [`stock-market-dynamics/REPORT.md`](stock-market-dynamics/REPORT.md) is the
  authoritative technical report (Revision 6);
  [`stock-market-dynamics/readme.md`](stock-market-dynamics/readme.md) is the model and
  pipeline summary. This file is the branch-level entry point.

---

## Contents

1. [Capabilities](#1-capabilities)
2. [Methodology](#2-methodology)
3. [Repository structure](#3-repository-structure)
4. [Installation](#4-installation)
5. [Quick start](#5-quick-start)
6. [Inputs and data](#6-inputs-and-data)
7. [Outputs](#7-outputs)
8. [Verification and validation](#8-verification-and-validation)
9. [Results](#9-results)
10. [Limitations](#10-limitations)
11. [Relationship to the other branches](#11-relationship-to-the-other-branches)
12. [License and contact](#12-license-and-contact)

---

## 1. Capabilities

| Capability | Where | Notes |
|---|---|---|
| Regime cluster state space | `src/build_regimes.py`, `src/rag.py` | 15 bins (momentum × volatility) in the prototype; 24 capital-weighted states (3 momentum × 2 volatility × 2 turnover × 2 size) in the validated engine |
| Master-equation evolution | `src/model.py`, `src/model4.py` | $d\mathbf{c}/dt = \mathbf{P} + \mathbf{S}\mathbf{J} - \mathbf{D}\mathbf{c}$, the same form RadCluster integrates |
| Measured transition kernel | `src/build_regimes.py` | month-to-month Markov transition matrix $T$ estimated from data, not assumed |
| Herding and sticky-drift kernels | `src/model4.py` | bilinear winner-chasing flux and a sticky-drift term, both calibrated — and both calibrate to **zero** on clean data |
| Market advection and spike damper | `src/model4.py` | mechanical advection plus a spike-persistence damper matched to observed post-peak decay |
| True historical membership | `src/build_panel_v2.py`, `src/fetch_shares_hist.py` | point-in-time S&P 500 constituents, avoiding backfilled-membership survivorship bias |
| Share-unit integrity repair | `src/fix_share_units.py` | segmentation and anchor-checked rescaling of filed share counts carrying $10^3$–$10^6$ unit errors |
| Conditional hindcasts | `src/hindcast*.py` | 2008–09 and 2020–21 concentration episodes, regime- and VIX-split variants |
| Scenario ensemble | `src/ensemble_forecast.py` | 200-scenario forward ensemble with explicit intervals |
| Fragility monitoring | `src/fragility_index.py` | top-lane share and HHI as percentiles of their own history |
| Portfolio risk | `src/portfolio_risk.py`, `src/paper_trade.py` | predicted tilt volatility, reversal statistics, paper-trading protocol |
| Systematic Buffett screen | `src/buffett_screen.py` | point-in-time quality/value/low-beta composite, annually rebalanced |
| Pre-registered negative tests | `src/comomentum.py`, `src/bubblechar.py`, `src/hf13f_test.py`, `src/dispersion_test.py`, `src/concentration_test.py` | each with its own write-up under `docs/` |

---

## 2. Methodology

### 2.1 The master equation

The state $\mathbf{c}(t)$ is the distribution of index capital across regime bins, and
evolves as

$$\frac{d\mathbf{c}}{dt} = \mathbf{P} + \mathbf{S}\mathbf{J} - \mathbf{D}\mathbf{c}$$

— the same three terms RadCluster integrates: a source, a stoichiometry-weighted sum
of fluxes, and a diagonal loss.

### 2.2 The RAG mapping

| Cluster dynamics (RadCluster) | Stock-Market Dynamics |
|---|---|
| cluster of size *n* (vertex) | regime bin: momentum × trailing volatility (× turnover × size in the validated engine) |
| concentration per size class | fraction, or capital share, of the index in each bin |
| cascade production (source) | index entry |
| growth, shrinkage, coalescence | measured Markov drift between bins; bilinear herding flux; mechanical market advection |
| rate kernels | transition matrix $T$, herding rate $h$, sticky-drift rate $\lambda$, spike-persistence damper |
| fixed sinks | index exit |
| provenance-stamped `output/` | `outputs/`, with the report and figures rebuilt from the run |

### 2.3 Calibration

Herding and sticky drift are free parameters, fitted rather than assumed. On the
repaired panel both calibrate to zero ($h_1 = \lambda_1 = 0$) — that is, once the data
are clean the measured drift plus advection account for the observed dynamics without
a behavioural term. The model keeps the kernels so the result can be re-tested, not
because they are needed.

---

## 3. Repository structure

```
/
├── README.md                       # this file — branch entry point
├── .gitignore                      # root-anchored; see Section 11
└── stock-market-dynamics/
    ├── REPORT.md                   # authoritative technical report (Revision 6)
    ├── readme.md                   # model and pipeline summary
    ├── src/                        # fetch → build → calibrate → hindcast → figures
    ├── data/                       # panels, transition matrices, test results
    ├── doc/                        # LaTeX manuscript and slides (+ build products)
    ├── docs/                       # one write-up per pre-registered test
    ├── figs/                       # figures
    └── outputs/                    # run outputs
```

---

## 4. Installation

```bash
git clone https://github.com/Ghoniem/RadCluster.git
cd RadCluster
git checkout Stock-Market-Dynamics
cd stock-market-dynamics
python -m pip install numpy scipy pandas matplotlib
```

Python 3.9 or newer. There is no compiled component: unlike `main`, this branch needs
no C++ toolchain and no SUNDIALS.

---

## 5. Quick start

The pipeline runs in order — fetch, build, calibrate, hindcast:

```bash
python3 src/fetch_data.py        # Yahoo Finance daily OHLCV (raw cache is not tracked)
python3 src/build_regimes.py     # bin constituents into regime states, estimate T
python3 src/hindcast.py          # baseline concentration hindcast
python3 src/hindcast_regime.py   # VIX-split, regime-dependent drift
python3 src/buffett_screen.py    # the quality/value/low-beta screen
python3 src/figures.py           # regenerate figures
```

Later revisions have their own entry points (`build_panel_v2.py`,
`fix_share_units.py`, `calibrate_rev5.py`, `hindcast_rev5.py`, `figures_rev5.py`);
`REPORT.md` says which revision produced which result.

---

## 6. Inputs and data

| Input | Source | Notes |
|---|---|---|
| Daily OHLCV | Yahoo Finance | 460 tickers, 2015–2021 in the prototype; extended to 2000–2022 in later revisions |
| Index membership | historical constituent lists | true point-in-time membership in `panel_v2.csv`; the prototype's backfilled list carries survivorship bias |
| Share counts | SEC company facts | repaired by `src/fix_share_units.py`; see Section 8 |
| Fundamentals | SEC company facts | net income, equity, assets for the quality and value scores |
| Volatility context | VIX | used for the regime-split drift test |

The raw price cache (~90 MB) is deliberately excluded; refetch with
`python3 src/fetch_data.py`.

---

## 7. Outputs

Results are written under `data/` and `outputs/` as CSV and JSON — transition
matrices (`Tcalm.npy`, `Tstress.npy`), monthly panels, per-test result files — with
figures under `figs/` and `doc/figs/`. The LaTeX manuscript and slides under `doc/`
are rebuilt from those artifacts, and `REPORT.md` is the narrative of record.

Note that this branch **tracks** its `.log`, `.aux`, `.out` and `.toc` files on
purpose: `data/fetch.log` and `data/regimes.log` are data, and the LaTeX products are
kept intentionally. See Section 11.

---

## 8. Verification and validation

**Conditional hindcasts.** The engine is validated against two concentration episodes.
On the repaired panel the hindcasts reach RMSE 0.0899 (2008–09) and 0.0714 (2020–21),
with the 2020–21 peak amplitude exact.

**Data integrity.** Filed share counts carried consistent $10^3$–$10^6\times$ unit
errors that a jump-based cleaner missed — in one case reading a utility at \$5,000
trillion, 400× the whole index, which flipped the measured top-lane share between
0.999 and 0.002 month to month. The repair segments each ticker at >20× split-adjusted
jumps and checks against two absolute anchors, market capitalisation in
[\$200 M, \$4 T] and average daily turnover in [0.02 %, 100 %]. Twenty-eight segments
were rescaled and seventeen dropped; **none were invented**, and the corrupt panel is
preserved at `data/panel_v2_preunitfix.csv` so the repair can be audited. Verification:
total binnable capitalisation \$10.2 T at 2007-12 against an index near \$13 T, and
\$41.3 T at 2021-12 against roughly \$40 T.

**Point-in-time discipline.** The screen rebalances each June on the latest fiscal year
ending at least 120 days earlier, the Fama–French end-date convention. Filing dates
are not used: SEC `filed` dates are unreliable before 2010, since restated facts carry
the restating filing's date.

---

## 9. Results

**What failed.** Seven pre-registered return-forecast tests were run against
2000–2022 S&P 500 data, and **all seven failed walk-forward validation** — comomentum
(a Lou & Polk non-replication), short interest, bubble characteristics, hedge-fund
crowding, dispersion, and the rest. Revision 6 removed their sections and figures from
the narrative rather than softening them. Each retains a write-up under
`stock-market-dynamics/docs/`.

**What survived.**

- *The regime capital-distribution engine*, as a census and scenario tool. A
  200-scenario ensemble twelve months out from 2022-12 puts the top-lane share at a
  median of 0.239 with a 90 % interval of [0.106, 0.503]. The interval width is the
  result: one-year concentration is dominated by driver and transition uncertainty, and
  a point forecast without intervals overstates knowledge roughly fourfold.
- *The systematic Buffett screen*, implementing the unlevered stock-selection leg of
  Frazzini, Kabiller & Pedersen's *Buffett's Alpha* decomposition — quality, value and
  betting-against-beta. The float leverage is noted and never applied, because a screen
  cannot replicate it.

**Standing verdict, unchanged since revision 4 and re-confirmed on clean data:** the
regime model does not forecast returns and must not drive capital deployment. Its
legitimate uses are regime census, conditional scenario analysis, and
concentration/fragility monitoring.

---

## 10. Limitations

`REPORT.md` holds the authoritative list. The ones to know before quoting any number:
the prototype universe is backfilled and therefore survivorship-biased, which is why
the validated engine moved to true historical membership; the share-count repair is a
repair, not a measurement, and is audited against the preserved pre-fix panel; and the
scenario ensemble describes the distribution of capital, not returns.

---

## 11. Relationship to the other branches

| Branch | Contents |
|---|---|
| `main` | The RadCluster suite — graph-based cluster dynamics for irradiated materials. |
| `Genealogy-Dynamics` | The same method applied to U.S. population genealogy. |
| `Stock-Market-Dynamics` | This branch. |

```bash
git checkout Stock-Market-Dynamics   # working tree becomes stock-market-dynamics/
git checkout main                    # RadCluster restored
```

Switching here deletes RadCluster's tracked files from the working tree — they remain
safe in `main` and in the object store — and leaves behind everything `main` did *not*
track: build output, run directories, virtual environments, several hundred megabytes
of it. The root `.gitignore` on this branch hides those leftovers so `git status`
stays readable.

That `.gitignore` is deliberately **narrower** than the one on `Genealogy-Dynamics`: it
does not ignore `*.log`, `*.aux`, `*.out` or `*.toc`, because this branch tracks
eighteen such files on purpose (Section 7). Do not copy the genealogy version over it.

Stage explicitly even so — `git add stock-market-dynamics/` rather than `git add -A` —
so an unrelated working tree can never be swept into this history.

---

## 12. License and contact

Released under the [MIT License](https://github.com/Ghoniem/RadCluster/blob/main/LICENSE),
© 2026 Nasr M. Ghoniem.

**Nasr M. Ghoniem** — Mechanical and Aerospace Engineering Department,
University of California, Los Angeles.
