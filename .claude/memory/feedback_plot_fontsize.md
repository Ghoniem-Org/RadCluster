---
name: 16pt fonts in plots
description: Plot text must be large and set via rcParams, never hardcoded; 16pt was the 2026-05-05 blanket rule, refined 2026-09-06 to a 22/20/16 scale
type: feedback
originSessionId: 72f427e0-249f-4dfb-a304-3919c5848d0f
---
Plot text must be set large, and set once via `plt.rcParams` at module load — never hardcoded per call.

**Superseded in part:** the original blanket rule was 16pt for everything. As of
2026-09-06 RadCluster uses a three-tier scale (labels 22, ticks 20, legend 16;
28/24/16 in the dose-figure script), so 16pt now describes the *legend* tier only.
See [[project_plot_house_style]] for the settled style. The rule below still holds:
never fall back to matplotlib's 10-12pt defaults, and never hardcode `fontsize=`.

**Why:** User stated this directly on 2026-05-05 as a permanent preference for plot readability. Default matplotlib sizes (10–12pt) are too small for the figures used in this project's notebooks and presentations.

**How to apply:**
- Set `plt.rcParams` for `axes.titlesize`, `axes.labelsize`, `xtick.labelsize`, `ytick.labelsize`, `legend.fontsize`, and `legend.title_fontsize` at module load time — at the 22/20/16 tiers, not a flat 16.
- Do not pass explicit `fontsize=` overrides on `legend()`, `set_title()`, `set_xlabel()`, etc. unless the user explicitly asks for a different size.
- When adding new plotting code, rely on the rcParams defaults rather than hard-coding sizes.
- The live wiring is `RadCluster_2_1/py_utils/visualization.py` (`_LABEL_FONTSIZE`/`_TICK_FONTSIZE`/`_PLOT_FONTSIZE` + the rcParams update at the top); new plots in any module should follow that pattern. The `RadCluster_1_0` path named here originally is archived.
