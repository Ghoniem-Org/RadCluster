---
name: project_void_network_loss
description: VOID_NETWORK_LOSS (cavity sweeping by dislocations) is the only mechanism that bounds cavity growth; the parameter levers structurally cannot
metadata: 
  node_type: memory
  type: project
  originSessionId: cd1032b8-b118-418f-a8f5-4745e69c19b5
  modified: 2026-09-06T13:07:54.819Z
---

`VOID_NETWORK_LOSS` — a moving network dislocation intersects a cavity and
absorbs it whole, delivering its m vacancies to the line:

    Lambda_m^void = |v_net| * rho_net * (chi * d_cav(m))

Added 2026-09-06 (commit on `main`), gated off by default, wired through
`reaction_rates.py` -> `cpp_bridge.py` -> `parameters.h` -> all three RHS
variants in `rate_kernels.cpp`, with the swept content charged to
`J_VAC_fixed`. OFF path verified byte-identical to the prior reference.

**Why no parameter could do this (all measured, not argued):**
- Cavity growth is `dm/dt ~ A_sph m^(1/3) * drive`, and the drive
  (`w_v c_v - w_i c_i`) is a **constant 30% of the vacancy capture rate from
  1e-5 to 40 dpa** — it never decays. So m ~ t^1.5, unbounded: the model's own
  answer at 40 dpa is ~7.7e7 vacancies (d ~ 118 nm). Row 9305's 7307 was
  purely the V=20000 ceiling.
- `E_m_v` **cancels**: it scales both `k2_vac` (∝ D_v) and the cavity capture
  coefficient `K_v` (∝ w_v), so it drops out of the growth rate. It is a
  nucleation knob (N_voids x41110 over its prior box, d_cavity +7.4%).
- The **bias factors destabilise the loops**. Loops survive on a NEGATIVE
  residual: `Z_i_loop w_i c_i - w_v c_v` = -4.69e-5, only -17% of its vacancy
  term. Z_v 1.0->1.15 halves that drive; measured loop content x65, d_111
  4.17 -> 13.9 nm, delta_FP 0.13. Z_v=1.30 zeroes it (bifurcation).
- A **conservation bound** blocks everything else: `S = S_I + Delta_J^d` pins
  cavity content to loop content, so with loops held fixed d_cavity cannot go
  below **4.70 nm** by any parameter choice. Sweeping escapes it by routing
  cavity vacancies straight into Delta_J^d.
- Because `Lambda ~ d_cav ~ m^(1/3)`, sweeping removes LARGE cavities
  preferentially and truncates the tail — the only size-dependent negative
  feedback in the model.

**How to apply:**
- Do not re-litigate `E_m_v` or the bias factors as cavity-size levers; the
  measurements are in and they cannot work. Five perturbations all shrank
  cavities only by starving them (N_void falling harder than d_cav).
- Sweeping **reduces** cavity density (it removes whole cavities); it does not
  subdivide content. Pair it with a nucleation increase, not stack it.
- `Lambda` is zero at build time because `ci1_seg`/`cv1_seg` are set by the
  segment refresh in `run_adaptive`, exactly like the loop channel. An
  all-zero array at build is expected, not a bug.
- Do NOT copy the loop channel's `P_ld` geometric gate onto cavities. It asks
  whether a cavity lies inside a line's instantaneous elastic zone; a sweeping
  dislocation traverses the volume, so `v*rho*w` IS the encounter rate and the
  gate double-counts, zeroing the channel at chi=1.
- See [[project_reference_run]] for the calibrated vector and its caveats.
