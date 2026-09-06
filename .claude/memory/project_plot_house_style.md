---
name: project_plot_house_style
description: "RadCluster's settled figure style as of 2026-09-06: no rendered titles, 22/20/16 pt type, lw=3, steelblue/darkviolet/tomato per population, driven by _PLOT_CONFIG in visualization.py"
metadata:
  type: project
---

The plotting style converged on 2026-09-06 (commit `1636750`). Match it in any
new figure code rather than reinventing sizes or colours.

**Titles are suppressed suite-wide.** `_PLOT_CONFIG['show_titles'] = False`;
every plot calls `_set_title(ax, ...)`, never `ax.set_title` directly. The
figures are captioned by the document or notebook that embeds them, so a baked-in
title duplicates the caption. `title=` kwargs are still plumbed through, so
`viz.set_plot_config({'show_titles': True})` restores them all, per-panel titles
included. The same `SHOW_TITLES` / `_set_title` pattern was extended to
`py_utils/loop_burgers_fraction.py` and `py_utils/size_distributions.py`, which
each carry their own module-level flag.

**Type scale** — set once in `plt.rcParams` at module load, not per call:

| | `visualization.py` | `codes/make_dose_figures.py` |
|---|---|---|
| axis labels / titles | 22 | 28 |
| tick labels | 20 | 24 |
| legend | 16 | 16 |
| dpi | 150 | 200 |

The dose figures run larger because they are the presentation headline figures.
`_INTERIOR_LEGEND_FONTSIZE = 12` is a deliberate exception for legends that sit
*inside* the axes, where 16 pt would cover the curves.

**Lines and colour.** `_LW = 3.0` for data, `_LW_THIN = 2.0` for reference or
overlay curves that must sit behind it. Colour encodes the population, suite-wide:
`_C_111 = 'steelblue'` (½⟨111⟩), `_C_100 = 'darkviolet'` (⟨100⟩), `'tomato'`
(voids/cavities). Labels come from `_LBL_111` / `_LBL_100`. In
`make_dose_figures.py` only, experimental points are coloured by *literature
source* with marker shape carrying the population (o = ½⟨111⟩, s = ⟨100⟩,
^ = cavity) — see [[project_dose_figures_script]].

**Axis knobs.** Plots are grouped into `concentration`, `scalar`, and `size_dist`;
each takes `xlim`/`ylim`/`xscale`/`yscale`, applied by `_apply_axis_config(ax, group)`
or `_apply_axis_config_to_fig(fig, group)` as the last step before `tight_layout`.
Override from a notebook with
`viz.set_plot_config({'concentration': {'ylim': (1e16, 1e22)}})` — partial dicts
merge. Several plots pin a deliberate bound *before* that call (e.g. `he_content`
floors y at 1e12, which clips free He out of view on purpose); read the comment
before "fixing" one.

Supersedes the older blanket "16 pt everywhere" note, which now describes only
the legend tier. See [[feedback_plot_fontsize]] and [[project_replot_from_plot_data]]
for redrawing existing runs in this style.
