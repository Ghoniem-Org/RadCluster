---
name: project_reference_run
description: "The RadCluster_2_1 reference vector is the VOID_NETWORK_LOSS 6/6 vector (chi=1165, f_cl_v=0.65), workbook default since 2026-09-06; it supersedes row 9305"
metadata: 
  node_type: memory
  type: project
  originSessionId: cd1032b8-b118-418f-a8f5-4745e69c19b5
  modified: 2026-09-06T13:07:34.442Z
---

The **reference vector** for RadCluster_2_1, designated 2026-09-06, is the
first vector in the campaign to put **all six observables in band at 330 C /
15 dpa**:

    VOID_NETWORK_LOSS = 1        (cavity -> network-dislocation sweeping)
    void_net_chi      = 1165
    f_cl_v            = 0.65     (fission)
    everything else   = calibration row 9305

It is the **workbook default** in `input_parameters.xlsx` as of 2026-09-06.
The previous reference was row 9305 carried to 40 dpa
(`output/20260904_212202_..._I80000V20000_im5vm5`), which scored 4/6 and whose
cavity size was a **grid artefact** — see below.

| observable | model | band | margin |
|---|---|---|---|
| N_100 | 4.715e21 | >= 4.67e21 | +0.96% |
| d_100 | 5.939 | 3.4 - 7 | comfortable |
| N_111 | 3.498e21 | >= 1.73e21 | comfortable |
| d_111 | 4.637 | 3.4 - 7 | comfortable |
| N_void | 3.609e20 | >= 3.6e20 | +0.25% |
| d_void | 2.890 | <= 2.9 | +0.34% |

**The real result is not "6/6" — it is that cavity size became a prediction.**
Occupancy `mean_n_v/V` is 0.054, against 0.382 for row 9305. Every earlier run
in this campaign sat at occ ~0.35 (median 0.347 over 411 rows): the cavity
distribution piled into the top bin and `d_cavity` was reading the V ceiling,
not physics. With sweeping on, the distribution turns over on its own.

**How to apply:**
- Treat this as canonical for "the reference case". See
  [[project_void_network_loss]] for the mechanism and why the parameter
  levers could not do this.
- **Three observables sit within 1% of a band edge, and it is fragile.**
  f_cl_v 0.65 -> 0.66 (a 1.5% change) drops it to 4/6. `learn.py`'s
  `worst_margin` would score it ~0.0025 and demote it. Report it as
  "sits on three band edges", never as a settled calibration. The robust
  sibling is chi=1150 / f_cl_v=0.65 at 5/6 (d_void 2.904, misses by 0.14%)
  with nearly identical physics.
- **chi and f_cl_v are COUPLED and must not be varied independently.**
  f_cl_v = 0.65 is stable only because sweeping is on: sweeping removes
  cavities, cutting the vacancy sink and raising c_v, which pushes the loop
  drive back negative. Turn VOID_NETWORK_LOSS off while keeping f_cl_v = 0.65
  and the loops run away (f_cl_v = 0.70 alone hung for 3h19m unconverged).
- **chi ~ 1e3 is a GLIDE-rate proxy, not an elastic capture radius.** It means
  the swept volume is ~1000x what network climb supplies (climb moves a line
  1.9 um over the whole 40 dpa). Presenting it as a capture radius would be
  indefensible. This needs physical sign-off before publication.
- Still open: delta_FP ~ 0.060 against the model's own 1e-2 gate, and the
  vector is NOT extent-verified, so `learn.py` will not promote it
  ([[project_learn_py_verified_gate]]).
- Figures come from `codes/make_dose_figures.py`, not the notebook
  ([[project_dose_figures_script]]). `output/` is gitignored, so run
  directories live only on this machine.
