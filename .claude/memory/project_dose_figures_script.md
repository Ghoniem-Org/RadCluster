---
name: project_dose_figures_script
description: "density_vs_dose.png / size_vs_dose.png come from codes/make_dose_figures.py, not from the notebook; the original script was lost in gitignored output/"
metadata: 
  node_type: memory
  type: project
  originSessionId: cd1032b8-b118-418f-a8f5-4745e69c19b5
  modified: 2026-09-05T19:29:09.904Z
---

`density_vs_dose.png` and `size_vs_dose.png` — the figures RadCluster_2_1 is
presented against — are **not produced by the notebook**. They come from
`RadCluster_2_1/codes/make_dose_figures.py` (committed 2026-09-05, `2f7a64c`).

**Why:** the script that originally drew the reference figures was written into
the gitignored `output/` tree and no longer exists. The only surviving relative,
`output/20260904_153005_f281fcf/make_figures.py`, is an earlier and different
variant — it plots the *campaign row's* `at_dose` trajectory from
`digital_twin/results/B3_coal.jsonl` rather than a run's own solution, and draws
a titled/marker style with no experimental bands. `docs/Database/SOURCES.md`
pointed readers at that dead path until it was corrected.

**How to apply:**
- Regenerate with `python codes/make_dose_figures.py <run_dir> --in-place`
  (add `--annotated` for the titled/marker variants). It reads the run's own
  `plots/plot_data.pkl`.
- It derives the six observables with the same expressions
  `digital_twin/run_ensemble.py` uses, so a figure and a ledger row for one
  vector agree by construction — keep them in sync if either changes.
- Two documented departures from the lost original: legend is 16 pt rather than
  12 (the project floor, see [[feedback_plot_fontsize]]; set `LEGEND_PT = 12` to
  match the originals exactly), and the density panel's y-floor is pinned at
  1e16 m^-3 instead of autoscaling — the un-nucleated <100> tail runs to ~1e13
  and costs four decades of vertical space.
- Anything generated into `output/` is gitignored and machine-local. Do not
  assume a figure or run directory survives; only the notebook config and this
  script are in version control. See [[project_reference_run]].
