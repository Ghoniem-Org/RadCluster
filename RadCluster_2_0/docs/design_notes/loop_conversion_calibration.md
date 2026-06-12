# Phase 6 — Loop-Conversion Calibration: Harness & Key Finding

**Status:** Harness built and validated; calibration surfaced a physics decision
(below) that needs the author's call before the 6 knobs can be fit.
**Date:** 2026-06-12.

## What was built

`codes/Python_Testing/calibrate_loop_conversion.py` — a 1-D calibration harness:

- **Experimental target** — `f111_experimental("EUROFER")` loads the ½⟨111⟩-loop
  fraction f₁₁₁(T) for EUROFER-97 from the microstructure database via
  `loop_burgers_fraction.extract_f111_table`. The EUROFER trend is a sharp
  crossover: **f₁₁₁ ≈ 1.0 at 250–300 °C → ≈ 0 by 335–400 °C** (crossover
  ~330–350 °C).
- **Model evaluation** — `f111_model(T_C, t_end, overrides)` runs the C++ solver
  (conversion ON) at temperature T to a fixed dose and returns
  f₁₁₁ = S_I^{111} / (S_I^{111} + S_I^{100}).
- **Sweep + T\*** search + a model-vs-experiment plot.

Performance note: high-T runs are fast (~4 s, I=60); low-T runs are stiff and
slow (>75 s). For exploratory sweeps the harness accepts a `G` (dose-rate)
override to compress the time axis; production calibration should use the
experimental G offline.

## Key finding — the model over-converts (decision required)

With the default parameters, **the model drives f₁₁₁ → 0 at all temperatures**,
including 300 °C (confirmed at 300 °C and 400 °C). This does **not** match the
EUROFER data (f₁₁₁ ≈ 1 below ~330 °C).

**Why (not a bug):** the two adopted channels have different temperature
character —

- **Dudarev unary** transformation *is* gated by the thermodynamic driving force
  ΔF(n,T) (the `T_star_conv_C` knob), so it correctly switches off at low T.
- **Marian junction + absorption** are **purely kinetic** — temperature
  independent. Once any ⟨100⟩ nucleus forms by a ½⟨111⟩+½⟨111⟩ junction, the
  absorption channel (⟨100⟩_m + ½⟨111⟩_n → ⟨100⟩_{m+n}, driven by the *full*
  mobile ½⟨111⟩ diffusivity) efficiently sweeps up the mobile ½⟨111⟩ at **any**
  temperature.

So the Marian channels dominate and **mask** the Dudarev T-trend: calibrating
`T_star_conv_C` alone cannot reproduce f₁₁₁(T), because the junction/absorption
keep converting below T\*.

## The decision

Reproducing the experimental f₁₁₁(T) crossover requires the Marian kinetic
channels to **switch off where ⟨100⟩ is not yet thermodynamically stable**. Three
defensible options (author's call — this is the debated-mechanism question):

1. **Thermodynamically gate the Marian channels (recommended).** Multiply the
   junction branching φ_junc(n,n′) and the absorption rate by a stability gate
   tied to ΔF (e.g. nonzero only where the ⟨100⟩ *product* size has ΔF > 0, or a
   smooth `1−exp(−ΔF/k_BT)` factor). Physical reading: a junction/absorption only
   yields a *stable* ⟨100⟩ where ⟨100⟩ is the lower-energy configuration;
   otherwise the segment reverts to ½⟨111⟩ (Marian's metastability + Dudarev's
   thermodynamics combined). This makes the kinetic channels T-dependent and lets
   the single crossover temperature control f₁₁₁(T).
   *Implementation:* gate `phi_junc` and `K_100_absorb` in `reaction_rates.py`
   and the C++ `conv_phi_junc`/`K_100_absorb` by the existing `conversion_dF`.

2. **Suppress the kinetic channels (φ_max ≪ 1, weak absorption)** so the unary
   thermodynamic channel dominates and T\* controls f₁₁₁(T). Simpler, but loses
   the Marian growth-to-TEM-size mechanism.

3. **Keep both channels ungated** (current model) and accept that f₁₁₁ → 0
   wherever loops form — i.e. the model predicts ⟨100⟩ dominance at all
   irradiation temperatures, contrary to the low-T EUROFER data.

## Recommendation

Adopt **Option 1** (gate the Marian junction/absorption by the ΔF stability
factor). It keeps both mechanisms, ties the kinetic channels to ⟨100⟩ stability,
and reduces the calibration to the dominant knob (the crossover temperature) plus
the junction/absorption *magnitudes*. After gating, the EUROFER crossover
(~330–350 °C) sets `T_star_conv_C ≈ 340 °C` (vs the current placeholder 450 °C).

This is a model-physics change, so it is left for the author's decision rather
than applied unilaterally — it is precisely the ½⟨111⟩↔⟨100⟩ mechanism question
flagged at the start of this work.
