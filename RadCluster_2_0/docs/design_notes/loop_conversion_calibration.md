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

---

## Resolution (2026-06-12) — Marian two-step success probability + direct-rotation unary

The author identified the missing physics: the conversion reaction rate lacked a
**probability of success** carrying the **two-step barrier** through the
metastable ½⟨110⟩ intermediate (Marian Fig. 3: ½⟨111⟩ →[ΔH₁≈0.5 eV]→ ½⟨110⟩
→[ΔH₂≈1.0 eV]→ ⟨100⟩, with E₍₁₀₀₎ < E₍₁₁₁₎ < E₍₁₁₀₎).

**Implemented.** From the metastable ½⟨110⟩ the segment branches forward to
⟨100⟩ or reverts to ½⟨111⟩:
$$P_{\rm success}(T) = \frac{e^{-\Delta H_2/k_BT}}{e^{-\Delta H_2/k_BT}+e^{-\Delta H_{\rm rev}/k_BT}}$$
(`dH2_conv`≈1.0, `dH_rev_conv`≈0.30 eV). This scalar multiplies **both** the
junction yield φ_junc and the absorption rate K_100_absorb
(`reaction_rates.py` + C++ `conv_phi_junc`/`K_100_absorb`, passed as
`loop_conv_psuccess`).

**Diagnostic that pinned the over-conversion.** With the unary channel disabled
(`E_a0_conv` large) and the gate on, **f₁₁₁ = 1.0 at 300 °C *and* 450 °C** — the
gated junction/absorption do *not* over-convert. So the over-conversion was
**entirely the unary channel**: its barrier `E_a0 = 0.8 eV` (with ν₀=10¹³)
converts loops fast even at 300 °C.

**Physical fix — the unary is a *direct* rotation.** A single loop with no
junction partner cannot use the easy two-step path; it must rotate
½⟨111⟩→⟨100⟩ **directly**, which Marian gives as **> 2 eV**. So `E_a0_conv` was
raised from 0.8 to ~2.0–2.5 eV. This is fully consistent with Marian: the direct
path is hard (unary), the two-step path is easy *because the partner enables it*
(junction/absorption). With the high direct barrier the unary correctly switches
off at low T and turns on near the crossover.

**Behaviour (confirmed).** The model now reproduces the f₁₁₁(T) crossover: f₁₁₁
stays ≈ 1 at low T and drops toward 0 at high T, with the crossover temperature
set by `E_a0_conv` and the **dose** (a rate-vs-time competition,
ν₀·e^{−E_a/kT}·t ~ 1 — physically real; ion vs neutron irradiation differ for
this reason). At ~2 dpa, `E_a0_conv` ≈ 2.5 eV pushes the crossover above 450 °C;
the EUROFER data (~16 dpa, ~340 °C) imply `E_a0_conv` ≈ 2.0–2.3 eV. The precise
value needs a **dose-matched** calibration sweep (slow, real G; the harness
supports it).

**Defaults set:** `E_a0_conv = 2.5` (placeholder in Marian's >2 eV range),
`dH2_conv = 1.0`, `dH_rev_conv = 0.30` eV. Junction/absorption gated by
P_success; conservation unaffected.

### Remaining calibration step
Run `calibrate_loop_conversion.py` at the experimental dose (~16 dpa, real G) to
fit `E_a0_conv` (and optionally `dH_rev_conv`) so the modeled crossover lands at
~340–350 °C. This is a slow offline sweep (low-T full-domain runs are stiff);
the dominant knob is now `E_a0_conv`.
