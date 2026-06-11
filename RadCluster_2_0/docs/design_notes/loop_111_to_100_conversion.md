# Design Note — ½⟨111⟩ → ⟨100⟩ Loop Conversion in the EUROFER-97 RAG

**Status:** Draft for review — *no code written yet.*
**Scope:** Add interstitial-loop Burgers-vector character (½⟨111⟩ vs ⟨100⟩) and
the ⟨111⟩→⟨100⟩ conversion physics to the RadCluster_2_0 reaction
admissibility graph (RAG).

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

**Adopted: parametrised $\Delta f(T)$.** We bypass the full
anisotropic-elasticity evaluation and parametrise the per-unit-length
free-energy difference, tying the ⟨100⟩ softening to the same $c'(T)$ law that
makes $\hat F_{001}([100])\to0$:
$$\boxed{\;\Delta f(T) = f_{111} - f_{100}^{0}\,(1-T/T_c)^{p}\;}\qquad
  \Delta F(n,T) = P(n)\,\Delta f(T),$$
with $T_c=1185$ K fixed, $f_{111}$ and $f_{100}^{0}$ the (≈ zero-T) per-length
energies of the two characters (eV/Å, scale set by the §3.1 core constants), and
the exponent $p$ chosen so $\Delta f$ changes sign near 350 °C and is strongly
positive by 550 °C (Dudarev Fig. 4 regions). Conversion is favoured for
$T>T^\*$ where $f_{111}=f_{100}^{0}(1-T^\*/T_c)^{p}$. The full anisotropic route
(Dudarev Eqs. 3–4) remains the high-fidelity upgrade behind the same kernel.

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
- Registered as a 1-D array over n; zero below the size where ΔF first turns
  positive — this fixes `n_conv` (the `n_min` of `bulk-100`) *self-consistently*
  from the energetics rather than as a free parameter.

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
emission, analogous to ½⟨111⟩'s $E_b^{\rm loop}(n)=A_{111}n^{-B_{111}}$:
$$E_b^{100}(n) = A_{100}\,n^{-B_{100}},$$
fit to the ⟨100⟩ formation-energy curve of §3.1 (⟨100⟩ loops are very stable, so
emission is slow but nonzero). The Excel `Physical_Properties` sheet already
carries the `b_100` Burgers vector.

### 3.6 New physics parameters

| Symbol | Meaning | First-cut value | Source / calibration |
|---|---|---|---|
| $E_a^{0}$ | unary barrier (offset) | ≈ 0.5–1.0 eV | Marian ΔH₂ |
| $\gamma_a$ | unary barrier size slope | tune | Arakawa onset + size window |
| $\nu_0$ | attempt frequency | ≈ 10¹³ s⁻¹ | Debye |
| $f_{111}, f_{100}^{0}$ | per-length loop energies | O(0.3–0.5) eV/Å | §3.1 core constants |
| $p$ | $\Delta f(T)$ softening exponent | tune | Dudarev Fig. 4 (sign change ≈350 °C) |
| $T_c$ | α–γ / spin-fluct. softening | 1185 K | Dudarev |
| δ | core cutoff | ≈ 0.4 nm | Dudarev |
| $\varphi_{\max}$ | peak junction yield | 0.5–1.0 | MD / calibration |
| $\sigma_s$ | log-size tolerance | 0.3–0.5 | Marian "comparable size" |
| $n_{j,\min}$ | min junction size | O(10) SIAs | Marian (junctions from n≈34–37) |
| $A_{100}, B_{100}$ | ⟨100⟩ binding fit | from §3.1 curve | — |
| $n_{\rm conv}$ | ⟨100⟩ onset size | *derived* from ΔF>0 | §3.2 (not free) |

---

## 4. Software structure

### 4.1 Populations (vertices)

| Population | Polarity | `n_min` | `mobile_max` | Role |
|---|---|---|---|---|
| `bulk-111` | SIA | 1 | `i_mobile` | 3-D-mobile small clusters + glissile ½⟨111⟩ loops; **SIA monomer pool**; **all cascade SIA source**; self-coalesces (now split by φ_junc) |
| `bulk-100` | SIA | `n_conv` | **0** (sessile) | ⟨100⟩ loops; no glide, no self-coalescence, no cascade source; grows by junction + absorption + monomer capture; strong V/I sink |
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

1. **Layer-1:** `COALESCENCE` `product_population` + walker gain-redirect
   (small, unit-testable against `δ_FP`).
2. **Energetics module:** `E_l^{111/100}(n,T)`, `ΔF(n,T)`, T-dependent elastic
   constants (or the parametrised `Δf(T)` first cut) — new helper alongside
   `binding_energies.py`.
3. **Kernels:** `K_111to100` (§3.2), `φ_junc` → `K_111_junction` and modified
   self-coal (§3.3), `K_100_absorb` (§3.4), `E_b^{100}` (§3.5) + Table entries.
4. **Layer-2 declaration:** two SIA populations, edges 1–9, monomer wiring.
5. **Validate** against the Python walker: conservation (δ_FP, δ_He) **and** the
   empirical `f₁₁₁(T)` / `f₁₁₁(dose)` trend from `loop_burgers_fraction.py` —
   calibrate the §3.6 parameters here.
6. **C++ mirror** last, gated on the Python reference.

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
