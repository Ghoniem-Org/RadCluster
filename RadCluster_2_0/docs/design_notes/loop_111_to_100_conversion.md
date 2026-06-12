# Design Note — ½⟨111⟩ → ⟨100⟩ Loop Conversion in the EUROFER-97 RAG

**Status:** Approved for implementation (2026-06-11). Work proceeds on branch
**`CodeDevelopment`** (main merged in; all 27 `py_utils` sources present).
**Scope:** Add interstitial-loop Burgers-vector character (½⟨111⟩ vs ⟨100⟩) and
the ⟨111⟩→⟨100⟩ conversion physics to the RadCluster_2_0 reaction
admissibility graph (RAG).

> **Existing hooks found in code (reuse, do not re-derive):**
> `binding_energies.py:216` `E_b_loop_i(...)` already accepts
> `A_100=0.7160, B_100=0.3581` (Table 18) — *dormant* (line 254 uses only
> `A_111/B_111`); wire these into a ⟨100⟩ emission curve `E_b_loop_100`.
> Material constants for the loop-energy module are already in that file:
> `E_f_i=3.77` eV, `G_shear=82e9` Pa, `b_111=2.49e-10` m, `nu=0.29`,
> `gamma_sf=0.6` J/m², `Omega=1.18e-29` m³, `n_tr=25`, `sigma_tr=5`.

**Decisions taken:**

- ⟨100⟩ loops: **fully sessile** (`mobile_max = 0`).
- Conversion: **one-way** (½⟨111⟩ → ⟨100⟩).
- Mechanism: **two additive channels**
  1. **Marian–Wirth–Perlado (PRL 88, 255507, 2002)** — *kinetic*: ⟨100⟩ nuclei
     form by collision/junction of two comparable-size mobile ½⟨111⟩ loops,
     then grow by absorbing further mobile ½⟨111⟩ clusters.
  2. **Dudarev–Bullough–Derlet (PRL 100, 135503, 2008)** — *thermodynamic*:
     a single ½⟨111⟩ loop spontaneously reorients to ⟨100⟩ when the
     anisotropic elastic free energy (softening of c′ = (c₁₁−c₁₂)/2 toward the
     α–γ transition) makes ⟨100⟩ the stable configuration; this channel is
     explicitly *not* collision-driven (Arakawa *et al.*, PRL 96, 125506, 2006).

The two channels are complementary: Marian supplies the low-/moderate-T
collision route and the dominant route to TEM-visible ⟨100⟩; Dudarev supplies
the high-T thermodynamic route and the observed temperature trend
(½⟨111⟩-dominated below ~350 °C → ⟨100⟩-dominated above ~550 °C). Their
**relative weight is the calibration target** against the repository's empirical
`f₁₁₁(T, dose)` reduction
([`py_utils/loop_burgers_fraction.py`](../../py_utils/loop_burgers_fraction.py)).

---

## 1. Motivation and current state

In bcc Fe / ferritic–martensitic steels, interstitial loops occur with two
Burgers characters:

- **½a⟨111⟩** — glissile (1-D glide), mobile to large sizes (n ≳ 100), the
  dominant cascade product.
- **a⟨100⟩** — sessile, square {100} habit, favoured at high T / large size.

**Current code.** The EUROFER RAG carries a *single* SIA population, `bulk`
([`declaration.py:211`](../../py_utils/materials/eurofer97/declaration.py)),
that lumps 3-D-mobile small clusters with glissile ½⟨111⟩ loops; ⟨100⟩ is not
represented. A `FUTURE HOOK` at
[`declaration.py:205-210`](../../py_utils/materials/eurofer97/declaration.py)
anticipates the split, and the architecture already supports it:

- `INTER_POPULATION` is one of the ten edge classes
  ([`edge_classes.py:41`](../../py_utils/core/edge_classes.py)) and is
  **implemented** in the walker
  ([`graph_walker.py:188`](../../py_utils/core/graph_walker.py)) — it carries
  the Dudarev unary channel with **zero** core change.
- `COALESCENCE` already supports a `partner_population ≠ population`
  ([`graph_walker.py:196`](../../py_utils/core/graph_walker.py)) — it carries
  the Marian *absorption* channel with **zero** core change; only the Marian
  *junction* channel needs a one-line core addition (§4).

---

## 2. Adopted physics

### 2.1 Marian channel — collision / junction (kinetic)

**Premise (consistent with our cascade model):** cascades inject only ½⟨111⟩;
these are highly mobile in 1-D. Two sub-processes:

**(a) Nucleation by junction.** Two mobile ½⟨111⟩ loops of *comparable size*
collide and react (Marian Eq. 4),
$$\tfrac12[111] + \tfrac12[1\bar1\bar1] \rightarrow [100],$$
condensing into a single ⟨100⟩ loop of the combined SIA content. Marian's MD
shows the junction forms *instantaneously* on collision, driven by the energy
reduction of merging, and **only when the two loops are of approximately the
same size** — otherwise the smaller loop simply rotates into the ½⟨111⟩
orientation of the larger (ordinary coalescence). The atomic path is the
modified Eyre–Bullough two-step (Marian Eqs. 5–6) with barriers ΔH₁ ≈ 0.5 eV
(½⟨111⟩→½⟨110⟩) and ΔH₂ ≈ 1.0 eV (½⟨110⟩→⟨100⟩); the direct ⟨111⟩↔⟨100⟩
rotation (> 2.0 eV) is negligible.

→ **cross-population COALESCENCE** `bulk-111 × bulk-111 → bulk-100`, gated by a
size-comparability branching fraction φ_junc(n,n′). The *complement*
(1 − φ_junc) is ordinary ½⟨111⟩ coalescence (stays in `bulk-111`).

**(b) Growth by absorption.** Once formed, ⟨100⟩{100} loops are essentially
immobile (glide barrier > 2.5 eV) and act as **biased sinks for mobile
cascade-produced ½⟨111⟩ clusters** (Marian Eq. 7),
$$[100]_m + \tfrac12\langle111\rangle_n \rightarrow [100]_{m+n},$$
the absorbed ½⟨111⟩ cluster rotating into the ⟨100⟩ orientation. This is the
dominant route by which ⟨100⟩ grows to TEM-visible sizes.

→ **cross-population COALESCENCE** `bulk-100 × bulk-111 → bulk-100`. Because the
product population equals the primary population, this needs **no** core change.

### 2.2 Dudarev channel — unary thermodynamic transformation

A *single* ½⟨111⟩ loop of size n reorients to ⟨100⟩ when ⟨100⟩ becomes the
lower-free-energy configuration. The driving force is the anisotropic elastic
free-energy difference, which is strongly T-dependent because the shear modulus
c′ = (c₁₁−c₁₂)/2 softens toward the α–γ transition (T_c ≈ 912 °C = 1185 K,
spin-fluctuation driven). The ⟨100⟩[100] pure-edge prelogarithmic energy factor
*vanishes* as c′→0 (Dudarev Eq. 4), while all ½⟨111⟩ factors stay finite — so
⟨100⟩ wins at high T. Dudarev's three stability regions (for ~10 nm loops):

| Region | Approx. T | Behaviour |
|---|---|---|
| 1 | T ≲ 350 °C | ½⟨111⟩ unconditionally stable — **no conversion** |
| 2 | 350–550 °C | ½⟨111⟩{110} unstable; may reorient to ⟨100⟩ — **partial** |
| 3 | T ≳ 550 °C | ⟨100⟩ unconditionally stable — **conversion favoured** |

Size enters through the loop perimeter; boundaries are diffuse for a size
distribution, exactly the smeared transition seen in `loop_burgers_fraction.py`.

→ **INTER_POPULATION** `bulk-111 → bulk-100` (unary, size-fixed), gated by the
thermodynamic driving force ΔF(n,T).

### 2.3 Why both

Dudarev explicitly argues the spontaneous transformation is *not* primarily
collisional (⟨100⟩ is absent in non-magnetic bcc metals where ½⟨111⟩ loops are
equally mobile), while Marian shows collisions *do* nucleate ⟨100⟩ and dominate
its growth. Neither alone reproduces the full (T, dose) behaviour: Marian's
junction yield is weakly T-dependent (kinetic), whereas the experimental
½⟨111⟩→⟨100⟩ crossover is sharply T-driven (thermodynamic). Adding them lets the
junction channel carry the dose/size dependence and the unary channel carry the
temperature dependence; the calibration reconciles both against `f₁₁₁(T,dose)`.

---

## 3. Reaction energetics and rate kernels

All three new kernels and the modified self-coalescence kernel are below. New
physics parameters are collected in the table at the end of this section.

### 3.1 Loop free energies (shared by both channels)

Prismatic-loop free energy (Dudarev Eq. 5; equivalently Marian Eq. 3):
$$E_l^{X}(n,T) = P_X(n)\Big[\,\hat F_X(T)\,\ln\!\frac{4R^\*_X(n)}{e\,\delta}
            \;+\; F_\delta^{X}(T) \;+\; F_c^{X}\,\Big],\qquad X\in\{111,100\}$$

with loop geometry from n SIAs (platelet area $A = n\,\Omega/b_X$, equivalent
radius $R^\*_X=\sqrt{A/\pi}$, core cutoff δ ≈ 0.4 nm):

| | sides $N_X$ | $b_X$ | habit | perimeter $P_X(n)$ |
|---|---|---|---|---|
| ½⟨111⟩ | 6 (hexagon) | $\tfrac{\sqrt3}{2}a$ | {110} | $6\sqrt{2A/3\sqrt3}$ |
| ⟨100⟩ | 4 (square) | $a$ | {100} | $4\sqrt{A}$ |

**Prelogarithmic factors** $\hat F_X(T)$ (eV/Å), anisotropic elasticity:
- ⟨100⟩[100] is analytic (Dudarev Eq. 4):
  $$\hat F_{001}([100]) = \frac{a^2}{4\pi}(c_{11}{+}c_{12})
    \!\left[\frac{c_{44}(c_{11}{-}c_{12})}{c_{11}(c_{11}{+}c_{12}{+}2c_{44})}\right]^{1/2}
    \;\xrightarrow{c'\to0}\;0 .$$
- ½⟨111⟩ has no closed form; evaluate from anisotropic elasticity (Dudarev
  Eq. 3) or tabulate from their Fig. 2. Weak T-dependence relative to ⟨100⟩.

**T-dependence** enters through the elastic constants, dominated by
$$c'(T)=\tfrac12(c_{11}-c_{12}) \approx 56.8\,(1-T/T_c)^{1/2}\ \text{GPa},
  \qquad T_c = 1185\ \text{K},$$
with $c_{44}$, $(c_{11}+c_{12})$ taken from Dever (1972) interpolated toward
$T_c$.

**Zero-T core constants** (Dudarev best fit to experiment, eV/Å):
| term | 111[11̄2] | 111[1̄10] | 001[100] | 001[110] |
|---|---|---|---|---|
| $F_c$ (nonlinear core) | 0.46 | 0.47 | 0.33 | — |
| $F_\delta$ (core-traction) | 0.345 | 0.349 | 0.387 | 0.390 |

**Adopted: parametrised $\Delta f(T)$, two-term form.** We bypass the full
anisotropic-elasticity evaluation but keep the physics that only the ⟨100⟩
*prelogarithmic* part softens (its core terms do not). Splitting $f_{100}$ into a
T-independent core and a softening prelog tied to the analytic $\hat
F_{001}([100])\propto\sqrt{c_{11}-c_{12}}\propto(1-T/T_c)^{1/4}$ scaling
(Dudarev Eq. 4):
$$\boxed{\;f_{100}(T) = f_{100}^{\rm core} + f_{100}^{\rm pre}\,(1-T/T_c)^{1/4},
   \qquad \Delta f(T) = f_{111} - f_{100}(T)\;}\qquad
   \Delta F(n,T) = P(n)\,\Delta f(T),$$
with $T_c=1185$ K fixed and the exponent **pinned at $1/4$** by Eq. 4 — *not* a
free knob. The two per-length energies $f_{100}^{\rm core}=(F_\delta+F_c)_{100}$
and $f_{100}^{\rm pre}=\hat F_{001}\ln(4R^\*/e\delta)$ come from the §3.1 core
constants; $f_{111}$ likewise. The single fit target is the crossover
temperature $T^\*$ (where $\Delta f=0$), anchored to the Dudarev Fig. 4 regions
(sign change ≈ 350 °C, strongly positive by ≈ 550 °C). This two-term form
replaces the earlier single-exponent $f_{111}-f_{100}^{0}(1-T/T_c)^p$, which was
numerically ill-conditioned when $f_{111}\approx f_{100}^{0}$. The full
anisotropic route (Dudarev Eqs. 3–4) remains the high-fidelity upgrade behind the
same kernel.

### 3.2 Dudarev unary kernel `K_111to100[n]` (1-D)

Thermodynamic driving force (favourable when > 0):
$$\Delta F(n,T) = E_l^{111}(n,T) - E_l^{100}(n,T).$$

Thermally activated, one-way, gated by the driving force, with a
**size-dependent barrier**:
$$\boxed{\;\Gamma_{\rm uni}(n,T) = \nu_0\,
   \exp\!\Big(-\frac{E_a^{\rm uni}(n)}{k_BT}\Big)\,
   \max\!\Big[0,\;1-\exp\!\Big(-\frac{\Delta F(n,T)}{k_BT}\Big)\Big],\qquad
   E_a^{\rm uni}(n) = E_a^{0} + \gamma_a\,\frac{P(n)}{b_{111}}\;}$$

- $\Delta F\le 0$ (low T / small n): gating → 0 (Region 1). ✓
- $\Delta F\gg k_BT$ (high T): gating → 1, rate → $\nu_0 e^{-E_a^{\rm uni}(n)/k_BT}$. ✓
- **Size dependence.** Coherent reorientation proceeds segment-by-segment
  (Marian's propagating ½⟨110⟩→⟨100⟩ front), so the effective barrier grows with
  the number of dislocation segments to reorient, $P(n)/b_{111}\propto\sqrt n$.
  $E_a^{0}$ ≈ 0.5–1.0 eV (Marian's per-step ΔH₂); $\gamma_a$ tunes the size
  suppression. Net effect: the gating favours large n (more $\Delta F$) while the
  barrier suppresses it — yielding a **preferred conversion-size window** that
  reproduces Arakawa's spontaneous transformation of *small* loops while leaving
  large mobile ½⟨111⟩ to convert via the junction/absorption channels instead.
- Registered as a 1-D array over n; **nonzero only on the thermodynamic
  support** where ΔF(n,T) > 0 (`LoopEnergetics.conversion_mask`).
  **Phase-2 finding (2026-06-11):** with the default Dudarev constants +
  Table-18 binding, the crossover temperature *increases with size*
  (n=20→−13 °C, n=50→450 °C, n=200→611 °C), so the support is **small-loop
  biased** — small loops unary-convert (Arakawa) while large mobile ½⟨111⟩
  stay and convert via the Marian junction/absorption channels. This is
  physically consistent and reinforced by the size-dependent barrier
  $E_a^{\rm uni}(n)$, which further suppresses large-loop unary conversion.
  *Correction:* the conversion support is therefore **not** a single
  lower-`n_conv` cutoff; `bulk-100`'s `n_min` is a separate **loop-onset
  floor** `n_loop_min` (≈ 4, below which the loop-energy formula is invalid
  and no ⟨100⟩ loop exists), independent of the ΔF>0 mask. The size direction
  and absolute crossover are sensitive to the (approximate) elastic constants
  and are a Phase-6 calibration decision against `f₁₁₁(T,dose)`.

### 3.3 Marian junction kernel `K_111_junction[n,n′]` (2-D)

Collision rate = existing same-polarity ½⟨111⟩ coalescence kernel
$\mathcal K^{ii}_{n,n'}$ (Eq. 79, built from the SIA diffusivities), split by a
size-comparability branching fraction:
$$\boxed{\;K^{\rm junc}_{n,n'} = \varphi_{\rm junc}(n,n')\,\mathcal K^{ii}_{n,n'},
   \qquad
   \varphi_{\rm junc}(n,n') = \varphi_{\max}\,
   \exp\!\Big[-\frac{(\ln(n/n'))^2}{2\sigma_s^2}\Big]\,
   \Theta\!\big(\min(n,n')\ge n_{j,\min}\big)\;}$$

- Peaked at n = n′ (Marian "approximately the same size"); $\sigma_s$ sets the
  tolerance, $\varphi_{\max}\!\le\!1$ the peak yield, $n_{j,\min}$ the minimum
  size for a stable junction.
- Product is ⟨100⟩ of size n+n′ → deposited into `bulk-100`.

**Modified existing self-coalescence** (stays in `bulk-111`):
$$\mathcal K^{ii,\,\rm self}_{n,n'} = \big(1-\varphi_{\rm junc}(n,n')\big)\,
   \mathcal K^{ii}_{n,n'}.$$
This re-uses `_coalescence_kernel(...)` then multiplies by (1−φ) / φ for the two
edges, conserving total ½⟨111⟩ collision rate (no double counting).

### 3.4 Marian absorption kernel `K_100_absorb[m,n]` (2-D)

Capture of a mobile ½⟨111⟩ cluster (size n, diffusivity $D^{111}_n$) by a
sessile ⟨100⟩ loop (size m, $D^{100}_m\approx0$):
$$\boxed{\;K^{\rm abs}_{m,n} = \frac{8\pi}{\Omega^{2/3}}
   (\xi_m+\xi_n)\big(D^{100}_m+D^{111}_n\big)
   \;\approx\; \frac{8\pi}{\Omega^{2/3}}(\xi_m+\xi_n)\,D^{111}_n\;}$$
i.e. the same 2-D coalescence kernel form (`_coalescence_kernel`) evaluated with
the **cross-population** diffusivity pair — the mobile ½⟨111⟩ partner drives the
capture. Optionally scaled by an absorption-and-rotation probability
$\varphi_{\rm abs}\!\approx\!1$. Product ⟨100⟩ of size m+n → `bulk-100`.

### 3.5 ⟨100⟩ point-defect kernels (sessile families)

⟨100⟩ loops still exchange single point defects: `GROWTH` (absorb I₁),
`SHRINKAGE` (absorb V₁), `DISSOCIATION` (emit I₁), `SINK`, `RECOMBINATION` and
`ANNIHILATION` with vacancies. These reuse the loop-capture form
$A_{\rm loop}\,n^{1/2}\,\omega^{\rm eff}$ with **sessile** mobility (so they
never enter coalescence transport), and a **⟨100⟩-specific binding curve** for
emission, analogous to ½⟨111⟩'s $E_b^{\rm loop}(n)=A_{111}n^{+B_{111}}$ in
[`binding_energies.py:216`](../../py_utils/binding_energies.py):
$$E_b^{100}(n) = A_{100}\,n^{+B_{100}}\ \text{(small-n)} \;\xrightarrow{\tanh}\;
   E_b^{\rm cont,100}(n)\ \text{(continuum)} .$$
**Already in code:** `E_b_loop_i(...)` accepts `A_100=0.7160, B_100=0.3581`
(Table 18) but currently ignores them (line 254 uses only `A_111/B_111`). The
implementation adds a sibling `E_b_loop_100(n)` that activates the dormant
⟨100⟩ fit with the *same* `n_tr`/`sigma_tr` blend and continuum tail. ⟨100⟩ loops
are very stable, so emission is slow but nonzero. The Excel `Physical_Properties`
sheet already carries the `b_100` Burgers vector.

### 3.6 Physics parameters

**Genuinely new (9)** — add to the Excel sheets:

| Symbol | Meaning | First-cut value | Source / calibration | Sheet |
|---|---|---|---|---|
| $E_a^{0}$ | unary barrier offset | ≈ 0.5–1.0 eV | Marian ΔH₂ | Model_Parameters |
| $\gamma_a$ | unary barrier size slope [eV/segment] | ≈ 0.03 (cal.) | Arakawa onset + size window | Model_Parameters |
| $\nu_0$ | attempt frequency | ≈ 10¹³ s⁻¹ | Debye; **reuse existing lattice ν₀** | Physical_Properties |
| $f_{111}$ | per-length ½⟨111⟩ energy | ≈ 1.6 eV/Å | §3.1 (Dudarev constants) | Physical_Properties |
| $f_{100}^{\rm core}, f_{100}^{\rm pre}$ | ⟨100⟩ core + softening prelog | core ≈ 0.7, pre ≈ 1.0 eV/Å | §3.1 (Dudarev constants) | Physical_Properties |
| $T^\*$ | conversion crossover temperature | fit, ∈ [350, 550] °C | Dudarev Fig. 4 | Model_Parameters |
| $\varphi_{\max}$ | peak junction yield | 0.5–1.0 | Marian MD | Model_Parameters |
| $\sigma_s$ | log-size tolerance | 0.3–0.5 | Marian "comparable size" | Model_Parameters |
| $n_{j,\min}$ | min junction size | ≈ 30 (10–35) | Marian (junctions from n≈34–37) | Model_Parameters |

**Fixed / already in code (do not add as free knobs):**

| Symbol | Value | Where |
|---|---|---|
| $A_{100}, B_{100}$ | 0.7160, 0.3581 | **already in** `binding_energies.py:216` (Table 18), dormant — just wire in |
| $T_c$ | 1185 K | constant (Dudarev) |
| softening exponent | $1/4$ (pinned) | Dudarev Eq. 4 — *not* a fit knob (see §3.1) |
| δ | 0.4 nm | constant (Dudarev) |
| $E_f^i, G, b_{111}, \nu, \gamma_{\rm sf}, \Omega, n_{\rm tr}, \sigma_{\rm tr}$ | — | already in `binding_energies.py` |
| conversion support | *derived* — ΔF>0 mask (small-loop biased) | §3.2 (not free) |
| $n_{\rm loop\_min}$ | ≈ 4 (loop-onset floor = `bulk-100` `n_min`) | constant |
| $T^\*$, $n_{\rm ref}$ | calibration anchor for `LoopEnergetics` | 450 °C, 50 (default) |

Net: **9 new sheet parameters**; of these, 6
($E_a^0,\gamma_a,T^\*,\varphi_{\max},\sigma_s,n_{j,\min}$ — all bounded by the
papers) are genuine calibration knobs against `f₁₁₁(T,dose)`; the remaining 3
($f_{111}, f_{100}^{\rm core}, f_{100}^{\rm pre}$) are computed from the Dudarev
energetics.

---

## 4. Software structure

### 4.1 Populations (vertices)

| Population | Polarity | `n_min` | `mobile_max` | Role |
|---|---|---|---|---|
| `bulk-111` | SIA | 1 | `i_mobile` | 3-D-mobile small clusters + glissile ½⟨111⟩ loops; **SIA monomer pool**; **all cascade SIA source**; self-coalesces (now split by φ_junc) |
| `bulk-100` | SIA | `n_loop_min` (≈4) | **0** (sessile) | ⟨100⟩ loops; no glide, no self-coalescence, no cascade source; grows by junction + absorption + monomer capture; strong V/I sink. `n_min` is the loop-onset floor, *not* the conversion support (which is the ΔF>0 mask, see §3.2). |
| `bulk` (VAC) | VACANCY | 1 | `v_mobile` | voids/bubbles — **unchanged** |

Monomer population stays `bulk-111` only; ⟨100⟩ GROWTH/SHRINKAGE auto-draw from
the ½⟨111⟩ n=1 pool ([`graph_walker.py:125`](../../py_utils/core/graph_walker.py)).

*Possible refinement (not now):* Marian distinguishes ⟨100⟩{110} (forms first)
from ⟨100⟩{100} (stable, n>68). Folded here into one `bulk-100` whose energetics
use the lower envelope; a habit-plane sub-population could be added later.

### 4.2 Edges

| # | Edge | Class | Endpoints | Kernel | Core change? |
|---|---|---|---|---|---|
| 1 | unary transform (Dudarev) | `INTER_POPULATION` | `bulk-111 → bulk-100` | `K_111to100[n]` (§3.2) | none |
| 2 | junction (Marian a) | `COALESCENCE` | `bulk-111 × bulk-111 → bulk-100` | `K_111_junction[n,n′]` (§3.3) | **yes — `product_population`** |
| 3 | absorption (Marian b) | `COALESCENCE` | `bulk-100 × bulk-111 → bulk-100` | `K_100_absorb[m,n]` (§3.4) | none (product = population) |
| — | ½⟨111⟩ self-coal (modify) | `COALESCENCE` | `bulk-111 × bulk-111 → bulk-111` | `(1−φ_junc)·𝒦ⁱⁱ` (§3.3) | none |
| 4–9 | ⟨100⟩ sessile families | GROWTH/SHRINKAGE/DISSOCIATION/SINK/RECOMBINATION/ANNIHILATION | `bulk-100` (+ vacancies) | §3.5 | none |

`bulk-111` retains its full existing P1–P8 edge set (only its self-coalescence
kernel changes). `bulk-100` gets **no** SOURCE and **no** self-COALESCENCE.

### 4.3 State layout

Add a third discrete (or bin-moment) block in `build_eurofer_rag`
([`declaration.py:448`](../../py_utils/materials/eurofer97/declaration.py)):

```python
layout.add_discrete("SIA111", I, population=sia_111)
layout.add_discrete("SIA100", I, population=sia_100)   # NEW
layout.add_discrete("VAC",    V, population=vac_bulk)
```

`N_eq` grows by ≈ I (or P·I_bin). He and conservation blocks unchanged. ⟨100⟩
sizes below `n_conv` carry zero conversion kernels, so the block is allocated
full-length I.

---

## 5. The one Layer-1 core change

Required *only* by the Marian junction edge (#2); the unary (#1) and absorption
(#3) edges use existing machinery.

1. **`Edge` / `COALESCENCE` spec.** Allow an optional `product_population` on a
   `COALESCENCE` edge and permit `changes_population = True` for that case.
   Files: [`core/rag.py:38`](../../py_utils/core/rag.py),
   [`core/edge_classes.py:127`](../../py_utils/core/edge_classes.py).

2. **`_c_coalescence` walker.** When `product_population` is set, redirect **only
   the gain deposit** `da[ij] += …` to the product block's view; leave the loss
   terms `da[i] -= …`, `db[j] -= …` on the source/partner blocks unchanged. The
   `same`-pair ½ factor and the `2.0*rate` gain factor are preserved verbatim —
   junction reuses the validated self-coalescence numerics, only the gain lands
   in `bulk-100`. File:
   [`core/graph_walker.py:196`](../../py_utils/core/graph_walker.py).

**Conservation is preserved structurally.** All three channels conserve signed
defect q = χn: the junction product has size n+n′, the absorption product m+n,
and the unary transfer relabels population at fixed size — all same polarity. The
`δ_FP` diagnostic (CLAUDE.md §8) is unchanged and serves as the correctness gate.

---

## 6. C++ mirror (deferred — the real implementation lift)

After the Python `GraphWalker` reference is validated, mirror in
`cpp_utils/rate_equations.cpp`, `cpp_utils/parameters.h`,
[`py_utils/cpp_bridge.py`](../../py_utils/cpp_bridge.py):

- second SIA block in the state vector and RHS;
- the unary inter-population term and the two cross-population coalescence
  deposits (junction → `bulk-100`, absorption → `bulk-100`);
- the φ_junc split of the ½⟨111⟩ self-coalescence;
- **bin-moment mode:** the junction and absorption edges couple two bin-moment
  SIA systems → reconstruct → transfer → project between the ½⟨111⟩ and ⟨100⟩
  bin grids.

---

## 7. Proposed build order

Files to touch are listed per step (all paths under `RadCluster_2_0/`).

1. **Layer-1 core** — ✅ **DONE (2026-06-11)**. `COALESCENCE` `product_population`
   + walker gain-redirect; 5/5 conservation tests pass.
   - `py_utils/core/edge_classes.py` (documented optional `product_population`)
   - `py_utils/core/rag.py` (`Edge.__post_init__`: same-polarity product validation)
   - `py_utils/core/graph_walker.py:196` (`_c_coalescence`: gain → `dprod`)
   - `codes/Python_Testing/check_coalescence_product_population.py` (new test)
   - ⚠️ **Pre-existing bug fixed here (Python reference only):** the
     same-population coalescence gain was `2.0*rate` (with a `*0.5` already
     applied) → it *created* SIA content (`d/dt Σn·cₙ = +0.76` in a no-over-top
     probe). Corrected to `rate`.
   - ✅ **C++ audited — correct, no factor-2.** `cpp_utils/.../rate_kernels.cpp`
     `K_ii_coal(n,np)` uses an *asymmetric single-diffusivity* convention
     (returns `…·D_np`, the projectile's D only): both orderings a→b and b→a
     fire separately and sum to the full `(D_a+D_b)` rate, so no symmetry factor
     is needed and nothing is over-produced. The Python reference used a
     *symmetric* `(D_a+D_b)` kernel and tried to compensate with `0.5`/`2.0` —
     the `2.0` was wrong. **Net: C++ production results are valid; the fix
     brings the Python reference into agreement with the C++.**
2. **Energetics module** — ✅ **DONE (2026-06-11)**. `E_l^{111/100}(n,T)`,
   two-term `Δf(T)` (§3.1, exponent pinned ¼), `ΔF(n,T)`, single-knob
   calibration to `(T*, n_ref)`, `conversion_mask`. 7/7 self-test checks pass.
   - **new** `py_utils/loop_energetics.py` (`LoopEnergetics` dataclass;
     `python py_utils/loop_energetics.py` runs the self-test)
   - `py_utils/binding_energies.py` (**new** `E_b_loop_100(n)` — activates the
     dormant `A_100=0.7160, B_100=0.3581`)
   - Finding: crossover T increases with size ⇒ small-loop-biased unary support
     (see §3.2). C++ port (Phase 7) must mirror this.
3. **Kernels** — ✅ **DONE (2026-06-11)**. 1-D ingredients in
   `py_utils/reaction_rates.py`: `Gamma_uni[n]` (§3.2, uses `LoopEnergetics` +
   size-dependent barrier, floored at `n_loop_min`), `phi_junc[n,n']` (§3.3),
   sessile ⟨100⟩ `K_100_grow/K_100_shrink/G_100/k2_100` (§3.5, `G_100` via
   `E_b_loop_100`). 18/18 tests pass
   (`codes/Python_Testing/check_loop_conversion_kernels.py`).
   Params read from `inp.reactions` with defaults (Excel rows added in Phase 5).
   - The 2-D `K_111_junction = φ·𝒦ⁱⁱ`, `(1−φ)·𝒦ⁱⁱ` self-coal split, and
     `K_100_absorb` are assembled in `declaration.py` (Phase 4) from `phi_junc`
     + `D_SIA_eff`, matching the existing 2-D-in-declaration architecture.
4. **Layer-2 declaration** — ✅ **DONE (2026-06-11)**. FUTURE HOOK realized:
   `bulk-111`/`bulk-100` populations; `_absorb_kernel` helper; the φ-split
   (`K_111_self`/`K_111_junction`), `K_100_absorb`, `Gamma_uni`, and sessile
   ⟨100⟩ kernels registered; conversion + ⟨100⟩-family edges added;
   `add_discrete("SIA100", …)` (and bin-moment block). Now **3 populations,
   24 edges**; discrete N_eq grows by I. Integration test
   `check_loop_conversion_integration.py` (structural + composed conservation,
   total SIA content flat to 0). `check_eurofer_rag.py` updated for the new
   population name + SIA100 block.
   - `py_utils/materials/eurofer97/declaration.py`
5. **Inputs** — 9 new sheet parameters (§3.6).
   - `py_utils/input_data.py`, `py_utils/create_excel.py`, `input/*.xlsx`
6. **Validate** against the Python walker: conservation (δ_FP, δ_He) **and** the
   `f₁₁₁(T)` / `f₁₁₁(dose)` trend from `loop_burgers_fraction.py`; calibrate the
   six §3.6 knobs here.
   - **new** `codes/Python_Testing/check_loop_conversion.py`
7. **C++ mirror** last, gated on the Python reference (§6).
   - `cpp_utils/parameters.h`, `cpp_utils/rate_equations.{cpp,h}`, `py_utils/cpp_bridge.py`

---

## 8. Summary

- **Two additive channels:** Marian collision (junction nucleation + ½⟨111⟩
  absorption growth, kinetic, carries dose/size dependence) and Dudarev unary
  transformation (thermodynamic, carries the T trend).
- **Graph topology:** two SIA populations; one `INTER_POPULATION` edge (unary);
  two cross-population `COALESCENCE` edges (junction, absorption); the existing
  self-coalescence kernel split by φ_junc; ⟨100⟩ sessile point-defect families.
- **Only one Layer-1 change:** add `product_population` to `COALESCENCE` (gain
  redirect) — needed solely by the junction edge.
- **All debated physics is localised** to three kernels + the loop-energy module,
  so the Marian/Dudarev weighting can be re-calibrated against `f₁₁₁(T,dose)`
  without touching the graph or the solver. `n_conv` is *derived* from
  ΔF(n,T) > 0, not a free parameter.

### Key references
- J. Marian, B. D. Wirth, J. M. Perlado, *Phys. Rev. Lett.* **88**, 255507 (2002).
- S. L. Dudarev, R. Bullough, P. M. Derlet, *Phys. Rev. Lett.* **100**, 135503 (2008).
- K. Arakawa *et al.*, *Phys. Rev. Lett.* **96**, 125506 (2006) (in-situ transformation).
- B. L. Eyre, R. Bullough, *Philos. Mag.* **12**, 31 (1965) (original shear reactions).
