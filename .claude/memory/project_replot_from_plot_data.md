---
name: project_replot_from_plot_data
description: "Any RadCluster figure can be redrawn from a run's plots/plot_data.pkl without re-integrating; results_y.npy is NOT sufficient because the <100> block is split off before it is written"
metadata:
  type: project
---

Every run directory carries `plots/plot_data.pkl` (~38 MB), written by
`visualization.save_plot_data`. It holds the **full** `results` dict plus
`input_data`, `rate_eq_obj` and `label` — so any plot function can be called
later with no re-integration:

```python
d = viz.load_plot_data(f'{run}/plots/plot_data.pkl')
viz.plot_mean_sizes_tem(d['results'], d['input_data'], d['rate_eq_obj'],
                        out_path=..., title=d['label'])
```

**Why:** `results_t.npy` / `results_y.npy` look like the raw-data artifacts but
are **not** replot-sufficient. `simulation.py` saves `results['y']`, which
`cpp_bridge` has already split the ⟨100⟩ block off from, so neither `y_sia100`
(per-size) nor `y_sia100_raw` (as solved) reaches those files — on the
2026-09-06 reference run `results_y.npy` is (189, 37) while the ⟨100⟩ block is
another 136 rows. Concluding "no ⟨100⟩ on disk ⇒ must re-run" from the `.npy`
files alone is wrong, and on that run would have cost a needless 22-minute
re-integration.

**Machine-local, though.** `.gitignore:41` excludes
`RadCluster_2_1/output/*/plots/plot_data.pkl` on purpose — 37 MB, 92 % of the
tracked reference directory, and a pickle of live `InputData`/`RateEquations`
objects that will not unpickle after a class refactor. So the pickle exists only
on the machine that ran the simulation: on a fresh clone the tracked reference
run has its PNGs, provenance and `.npy` files but no `plot_data.pkl`, and there
the figures genuinely are not regenerable without re-integrating. Check the file
is present before promising a cheap replot.

**How to apply:** before proposing a re-run to fix a figure, look in
`plots/plot_data.pkl`. Redraw only the affected figure rather than calling
`save_all_plots`, which rewrites all ~26 PNGs and re-pickles the 38 MB payload.
Note `save_all_plots` does *not* write `density_vs_dose.png` / `size_vs_dose.png`
— those come from [[project_dose_figures_script]]. See also
[[project_reference_run]].
