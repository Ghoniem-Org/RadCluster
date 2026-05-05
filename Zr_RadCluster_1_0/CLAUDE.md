# RadCluster_1_0 — Physics and Solver Reference

## Source Document

All equations cite:
> Ghoniem, N.M. (2026), *"A Cluster Dynamics Model for Radiation Damage Evolution
> in Ferritic-Martensitic Steels"* (`docs/Formulation/rate_equations.tex`)

---

## 1. State Vector

### Naming convention

| Symbol | Meaning |
|---|---|
| $I$ | max SIA cluster size (code: `I`, was `N`) |
| $V$ | max vacancy cluster size (code: `V`, was `M`) |
| $i_{\rm mobile}$ | max mobile SIA size (code: `i_mobile`, was `n_max_i`) |
| $v_{\rm mobile}$ | max mobile vacancy size (code: `v_mobile`, was `m_max_v`) |
| $i_{\rm discrete}$ | max discrete SIA size (individually tracked) |
| $v_{\rm discrete}$ | max discrete vacancy size (individually tracked) |
| $I_{\rm bin}$ | number of SIA bin-moment equations beyond $i_{\rm discrete}$ |
| $V_{\rm bin}$ | number of VAC bin-moment equations beyond $v_{\rm discrete}$ |
| $i_{\rm cascade}$ | max SIA cluster size from cascade (was `m1`) |
| $v_{\rm cascade}$ | max vacancy cluster size from cascade (was `n1`) |

### full_CD modes  (Eqs. ME\_SIA, ME\_vac, ME\_He)

| Segment | Indices | Variable |
|---|---|---|
| SIA clusters | $0 \ldots I-1$ | $c_i$, $i = 1,\ldots,I$ |
| Vacancy/He grid | $I \ldots I + V(\ell_{\max}+1) - 1$ | $c_{v,\ell}$, $v=1,\ldots,V$, $\ell=0,\ldots,\ell_{\max}(v)$ |
| Free He | last | $c_h$ |

Total equations: $I + V(\ell_{\max}+1) + 1$.

**fission** physics option uses He-reduction **Case 2** (Eq. Ebv\_eff, decoupled):
decoupled scalar $\ell_{\rm tot}$ → reduced system size $I + V + 2$.

**fusion** physics option uses He-reduction **Case 1** (Eq. case1\_cm, mean-field):
mean-field $\bar{\ell}(v)$ per void class → reduced system size $I + 2V + 1$.

### bin_moment_CD modes  (Chapter 9 — hybrid discrete + bin-moment)

Sizes $1 \ldots i_{\rm discrete}$ tracked individually (one ODE per size).
Sizes $i_{\rm discrete}+1 \ldots I$ grouped into $I_{\rm bin}$ logarithmic bins,
each with $P$ moments (Eq. bin\_moments\_def):

$$\mu_k^{(0)} = \sum_{n \in \mathcal{B}_k} c_n \qquad \text{(zeroth moment — bin density)}$$

$$\mu_k^{(1)} = \sum_{n \in \mathcal{B}_k} n\,c_n \qquad \text{(first moment — bin content)}$$

$$\mu_k^{(2)} = \sum_{n \in \mathcal{B}_k} n^2\,c_n \qquad \text{(second moment — lognormal only)}$$

Same for vacancies: $v_{\rm discrete}$ discrete + $V_{\rm bin}$ bins with $P$ moments each.

The number of moments per bin $P$ is determined by the `shape_function` parameter:

| `shape_function` | $P$ | Closure | Truncation error |
|---|---|---|---|
| `"constant"` | 1 | Piecewise-constant (Eq. pc\_ansatz) | $O((r-1)^2)$ |
| `"linear"` | 2 | Hat-function / dual-basis (Eq. shape\_ansatz) | $O((r-1)^3)$ |
| `"lognormal"` | 3 | Log-normal shape (Eq. lognormal\_ansatz) | $O((r-1)^4)$ |

System size: $N_{\rm eq} = i_{\rm discrete} + P \cdot I_{\rm bin} + v_{\rm discrete} + P \cdot V_{\rm bin} + n_{\rm He}$.

When $I_{\rm bin} = 0$ and $i_{\rm discrete} = I$, all equations are discrete → recovers full\_CD.

---

## 2. Cascade Production Spectrum  (Eqs. G\_eta–Ci, Tables 2, 5)

Survival-corrected production rate (Eq. G\_eta):

$$G = \eta\, G_{\rm NRT}$$

Power-law cluster production fraction (Eq. eps\_i\_model):

$$\epsilon_m^{(i)} = C_i\, m^{-s_i}, \qquad m = 2,\ldots,m_1$$

Normalisation constant (Eq. Ci):

$$C_i = \frac{f_i^{\rm cl}}{\displaystyle\sum_{m=2}^{m_1} m^{1-s_i}}, \qquad
  C_v = \frac{f_v^{\rm cl}}{\displaystyle\sum_{n=2}^{n_1} n^{1-s_v}}$$

Production rates (Eq. Pmi):

$$P_m^{(i)} = \epsilon_m^{(i)}\, G, \qquad
  P_1^{(i)} = G\!\left(1 - \sum_{m=2}^{m_1} m\,\epsilon_m^{(i)}\right) \qquad [\text{atom fraction s}^{-1}]$$

| Parameter | Fission | Fusion |
|---|---|---|
| $\eta$ (survival) | 0.30 | 0.28 |
| $f_i^{\rm cl}$ | 0.58 | 0.65 |
| $s_i$ | 1.6 | 1.5 |
| $m_1$ | 20 | 50 |
| $f_v^{\rm cl}$ | 0.15 | 0.20 |
| $s_v$ | 2.5 | 2.3 |
| $n_1$ | 10 | 20 |

---

## 3. Defect Diffusion and Solute Trapping  (Eqs. Di–D3D\_eff\_alloy)

Point-defect diffusivity (Eq. Di–Dh):

$$D_\alpha = a^2\,\nu_\alpha \exp\!\left(-\frac{E_m^\alpha}{k_{\rm B} T}\right)$$

Jump frequencies (Eqs. omega\_i–omega\_v):

$$\omega_\alpha = \nu_\alpha \exp(-E_m^\alpha / k_{\rm B} T) = D_\alpha / a^2 \qquad [\text{s}^{-1}]$$

EUROFER solute-trapping effective frequency (Eq. omega\_eff):

$$\omega_\alpha^{\rm eff} = \frac{\nu_\alpha\,\exp(-E_m^\alpha / k_{\rm B} T)}{1 + \displaystyle\sum_s z_s\,c_s \exp\!\left(\dfrac{E_b^{\alpha\text{-}s}}{k_{\rm B} T}\right)}$$

This applies to SIA, vacancy, and He. SIA clusters are further modified via loop
trapping (Eq. D1D\_eff\_Eurofer).

SIA cluster 1D glide (Eq. D1D\_boxed):

$$D_n^{\rm 1D} = \frac{3a^2\,\nu_0}{2\,n^{s}} \exp\!\left(-\frac{E_m^{\rm 1D}}{k_{\rm B} T}\right), \qquad s \approx 0.5\text{–}1.0$$

Mean free path in EUROFER (Eq. L\_eff): $L_{\rm eff} = D_{n,\rm eff}^{\rm 1D} / \nu_{R,\rm eff}$.

Mobility cutoffs: $n_{\max}^i = 100$ SIAs, $m_{\max}^v = 5$ vacancies.

---

## 4. Reaction Rate Constants  (Eqs. K3D–P8, Tables 28, 30)

Geometric prefactors (Eq. constants):

$$A_{\rm sph}  = (48\pi^2)^{1/3} \approx 7.818, \qquad
  A_{\rm loop} = 8\sqrt{\pi/\!\sqrt{3}} \approx 10.78$$

$$A_{\rm 1D}   = \frac{9}{8\pi^{2/3}} \approx 2.632, \qquad
  B_{\rm rot}  = \frac{4}{\pi}\!\left(\frac{8\pi}{3}\right)^{\!1/3} \approx 2.627$$

**P1: Vacancy–SIA recombination** (Eq. P1):

$$\mathcal{K}_{iv} = 4\sqrt{3}\,\pi\,(\omega_i^{\rm eff} + \omega_v^{\rm eff}) \approx 21.8\,\omega_i^{\rm eff}$$

**P2: Point-defect absorption by spherical cavity** (Eq. P2):

$$\mathcal{K}_{\alpha,m}^{\rm cav} = A_{\rm sph}\, m^{1/3}\, \omega_\alpha^{\rm eff}$$

**P3: Point-defect absorption by dislocation loop** (Eqs. P3\_i, P3\_v):

$$\mathcal{K}_{i,n}^{\rm loop} = A_{\rm loop}\, n^{1/2}\, Z_i^{\rm loop}\, \omega_i^{\rm eff}, \qquad
  \mathcal{K}_{v,n}^{\rm loop} = A_{\rm loop}\, n^{1/2}\, \omega_v^{\rm eff}$$

**P4: Fixed sinks** (Eqs. D\_total, P4d–P4p):

$$\mathcal{D}_\alpha = \mathcal{D}_\alpha^d + \mathcal{D}_\alpha^{gb} + \mathcal{D}_\alpha^p \qquad [\text{s}^{-1}]$$

$$\mathcal{D}_\alpha^d = Z_\alpha^d\,\rho_d\,a^2\,\omega_\alpha^{\rm eff}, \qquad
  \mathcal{D}_\alpha^{gb} = \frac{4\pi^2}{d_g^2}\,a^2\,\omega_\alpha^{\rm eff}, \qquad
  \mathcal{D}_\alpha^p = \frac{8\pi\,r_p}{a}\,c_p^N\,\omega_\alpha^{\rm eff}$$

**P5: Thermal emission** (Eqs. P5v, P5i, P5h):

Vacancy emission from cavity $(m,\ell)$:

$$\varepsilon_v(m,\ell) = A_{\rm sph}\,(m{-}1)^{1/3}\,\omega_v^{\rm eff}\,\exp\!\bigl(-E_b^v(m,\ell)/k_{\rm B}T\bigr)$$

SIA emission from loop of $n$ SIAs:

$$\varepsilon_i(n) = A_{\rm loop}\,(n{-}1)^{1/2}\,\omega_i^{\rm eff}\,\exp\!\bigl(-E_b^i(n)/k_{\rm B}T\bigr)$$

He emission from bubble $(m,\ell)$ — uses $m^{1/3}$ (He removal doesn't change vacancy count):

$$\varepsilon_h(m,\ell) = A_{\rm sph}\,m^{1/3}\,\omega_h^{\rm eff}\,\exp\!\bigl(-E_B^{\rm He}(m,\ell)/k_{\rm B}T\bigr)$$

**P6: SIA cluster–cavity interaction (mixed 1D/3D)** (Eq. P6):

$$\mathcal{K}_{n,m}^{\rm eff} = \frac{A_{\rm sph}\, m^{1/3}\, \omega_n^{\rm 1D}}{1 + B_{\rm rot}\,\hat{L}^2\, m^{-1/3}},
  \qquad \hat{L} = L/a$$

Full cavity rate constant with 3D/1D transition (Eq. K\_cav):

$$\mathcal{K}_{n,m}^{(\rm cav)} = \begin{cases}
  A_{\rm sph}\,m^{1/3}\,\omega_i^{\rm eff},
    & 1 \le n \le 3 \text{ (3D)}, \\[4pt]
  \dfrac{A_{\rm sph}\,m^{1/3}\,\omega_n^{\rm 1D}}{1 + B_{\rm rot}\,\hat{L}^2\,m^{-1/3}},
    & 4 \le n \le n_{\max}^i \text{ (1D/3D)}, \\[4pt]
  0,
    & n > n_{\max}^i \text{ (sessile)}.
\end{cases}$$

**P7: Trap mutation** (Eq. P7):

$$\Gamma_{\rm TM}(m, \ell) = \nu_0 \exp\!\left(-\frac{E_{\rm TM}(m,\ell)}{k_{\rm B} T}\right)$$

**P8: Radiation re-solution** (Eq. P8):

$$\Gamma_{\rm res}(m, \ell) = b_0\, \ell\, \dot{\phi}$$

---

## 5. Master Equations  (Eqs. ME\_SIA, ME\_vac, ME\_He)

All concentrations are atomic fractions; rate constants carry units of s$^{-1}$.
Equations below are the **general (coalescence) form**: every pair of mobile clusters
that can react does so, governed by $\mathcal{K}^{ii}$, $\mathcal{K}^{vv}$, and $\mathcal{K}^{vi}$.

General cluster–cluster rate constant (3D, spherical):

$$\mathcal{K}_{n,n'}^{ii} = 8\pi\,(\xi_n + \xi_{n'})\,(\omega_n^{\rm eff} + \omega_{n'}^{\rm eff}),
  \qquad \xi_n = \left(\frac{3n}{8\pi}\right)^{\!1/3}$$

(analogous forms for $\mathcal{K}^{vv}$ and $\mathcal{K}^{vi}$).

The **mono-defect form** ($n_{\max}^i = 1$, $m_{\max}^v = 1$) is recovered by restricting
all coalescence sums to their $n'=1$ or $m'=1$ terms, reducing
$\mathcal{K}^{ii} \to \mathcal{K}_{i,n}^{\rm loop}$,
$\mathcal{K}^{vv} \to \mathcal{K}_{\alpha,m}^{\rm cav}$,
$\mathcal{K}^{vi} \to \mathcal{K}_{iv}$ (Eqs. monodef\_i–monodef\_iv).

---

**SIA clusters** (Eq. ME\_SIA), $n = 1,\ldots,N$:

$$\frac{dc_n}{dt} =
  \underbrace{G_n}_{\text{production}}$$
$$+ \underbrace{\frac{1}{2}\sum_{n'=1}^{\min(n-1,\,n_{\max}^i)}
  \mathcal{K}_{n',n-n'}^{ii}\,c_{n'}\,c_{n-n'}}_{\text{SIA–SIA coalescence gain } (I_{n'} + I_{n-n'} \to I_n)}$$
$$- \underbrace{c_n \sum_{n'=1}^{n_{\max}^i}
  \mathcal{K}_{n,n'}^{ii}\,c_{n'}}_{\text{SIA–SIA coalescence loss } (I_n + I_{n'} \to I_{n+n'})}$$
$$+ \underbrace{\varepsilon_i(n{+}1)\,c_{n+1}
  - \varepsilon_i(n)\,c_n\,\delta_{n \ge 2}}_{\text{thermal emission (SIA from loops)}}$$
$$+ \underbrace{\sum_{m'=1}^{m_{\max}^v}
  \mathcal{K}_{m',n+m'}^{vi}\,c_{m',0}\,c_{n+m'}}_{\text{V–I annihilation gain } (V_{m'} + I_{n+m'} \to I_n)}$$
$$- \underbrace{\sum_{m'=1}^{m_{\max}^v}
  \mathcal{K}_{m',n}^{vi}\,c_{m',0}\,c_n}_{\text{V–I annihilation loss } (V_{m'} + I_n \to I_{n-m'})}$$
$$- \underbrace{\sum_m \mathcal{K}_{n,m}^{(\rm cav)}\,c_n\,\bar{c}_m}_{\text{SIA cluster–cavity absorption}}$$
$$- \underbrace{\mathcal{D}_i(n)\,c_n}_{\text{fixed-sink loss (mobile only: } \mathcal{D}_i = 0 \text{ for } n > n_{\max}^i\text{)}}$$
$$+ \underbrace{\delta_{n,1}\sum_{m,\ell}\Gamma_{\rm TM}(m,\ell)\,c_{m,\ell}}_{\text{trap mutation (SIA source at } n=1\text{)}}$$

---

**Vacancy/He clusters** (Eq. ME\_vac), $m = 1,\ldots,M$, $\ell = 0,\ldots,\ell_{\max}(m)$:

$$\frac{dc_{m,\ell}}{dt} =
  \underbrace{G_{m,\ell}}_{\text{production}}$$
$$+ \underbrace{\frac{1}{2}\sum_{m'=1}^{\min(m-1,\,m_{\max}^v)}
  \mathcal{K}_{m',m-m'}^{vv}\,c_{m',0}\,c_{m-m',\ell}}_{\text{V–V coalescence gain } (V_{m'} + V_{m-m'} \to V_m)}$$
$$- \underbrace{c_{m,\ell}\sum_{m'=1}^{m_{\max}^v}
  \mathcal{K}_{m,m'}^{vv}\,c_{m',0}}_{\text{V–V coalescence loss } (V_m + V_{m'} \to V_{m+m'})}$$
$$+ \underbrace{\varepsilon_v(m{+}1,\ell)\,c_{m+1,\ell}
  - \varepsilon_v(m,\ell)\,c_{m,\ell}}_{\text{thermal emission (vacancy from cavities)}}$$
$$+ \underbrace{\sum_{n=1}^{n_{\max}^i}
  \mathcal{K}_{n,m+n}^{(\rm cav)}\,c_n\,c_{m+n,\ell}}_{\text{SIA-induced shrinkage gain } (I_n + V_{m+n} \to V_m)}$$
$$- \underbrace{\sum_{n=1}^{n_{\max}^i}
  \mathcal{K}_{n,m}^{(\rm cav)}\,c_n\,c_{m,\ell}}_{\text{SIA-induced shrinkage loss } (I_n + V_m \to V_{m-n})}$$
$$+ \underbrace{A_{\rm sph}\,m^{1/3}\,\omega_h^{\rm eff}\,c_h\,c_{m,\ell-1}\,\delta_{\ell \ge 1}
  - A_{\rm sph}\,m^{1/3}\,\omega_h^{\rm eff}\,c_h\,c_{m,\ell}\,\delta_{\ell < \ell_{\max}}}_{\text{He absorption by cavities}}$$
$$+ \underbrace{\varepsilon_h(m,\ell{+}1)\,c_{m,\ell+1}
  - \varepsilon_h(m,\ell)\,c_{m,\ell}\,\delta_{\ell \ge 1}}_{\text{He emission from cavities}}$$
$$+ \underbrace{\Gamma_{\rm TM}(m{-}1,\ell)\,c_{m-1,\ell}
  - \Gamma_{\rm TM}(m,\ell)\,c_{m,\ell}}_{\text{trap mutation}}$$
$$+ \underbrace{b_0\,\dot{\phi}\,(\ell{+}1)\,c_{m,\ell+1}
  - b_0\,\dot{\phi}\,\ell\,c_{m,\ell}}_{\text{radiation re-solution}}$$
$$- \underbrace{\mathcal{D}_v\,c_{m,\ell}\,\delta_{m \le m_{\max}^v}}_{\text{fixed-sink loss (mobile clusters only)}}$$

---

**Free He** (Eq. ME\_He):

$$\frac{dc_h}{dt} =
  \underbrace{G_{\rm He}}_{\text{production}}$$
$$- \underbrace{\sum_{m=1}^{M} A_{\rm sph}\,m^{1/3}\,\omega_h^{\rm eff}\,c_h
  \sum_{\ell=0}^{\ell_{\max}(m)-1} c_{m,\ell}}_{\text{He absorption by cavities}}$$
$$- \underbrace{\mathcal{D}_h\,c_h}_{\text{fixed-sink loss}}$$
$$+ \underbrace{\sum_{m=1}^{M}\sum_{\ell=1}^{\ell_{\max}(m)}
  \varepsilon_h(m,\ell)\,c_{m,\ell}}_{\text{He emission from cavities}}$$
$$+ \underbrace{\sum_{m=1}^{M}\sum_{\ell=1}^{\ell_{\max}(m)}
  b_0\,\dot{\phi}\,\ell\,c_{m,\ell}}_{\text{radiation re-solution}}$$

> **Note:** He–He coalescence is excluded because interstitial He–He binding
> energy in bcc Fe is negligible ($\lesssim 0.05$ eV).

---

## 6. He-Vacancy State-Space Reduction  (Section 8)

### Case 1 — Fusion (fast He equilibration, Eq. case1\_cm)

Track marginal $\bar{c}_m = \sum_\ell c_{m,\ell}$ and mean loading (Eq. lbar\_def):

$$\bar{\ell}(m,t) = \frac{\sum_\ell \ell\,c_{m,\ell}}{\bar{c}_m}$$

Auxiliary ODE for He content per class $N_{\rm He}(m) = \bar{\ell}(m)\,\bar{c}_m$ (Eq. case1\_NHe).
System size: $N + 2M + 1$.

### Case 2 — Fission (decoupled, Eq. Ebv\_eff)

Track only $\bar{c}_m$ and scalar $\ell_{\rm tot} = \sum_m \ell(m)\,\bar{c}_m$.
He-corrected effective vacancy binding energy (Eq. Ebv\_eff):

$$E_b^{v,\rm eff}(m) = E_f^v + \frac{2\gamma_s\,\Omega}{r_m} - \frac{\ell(m)\,k_{\rm B}T}{m}$$

He allocation proportional to capture rate (Eq. allocation):

$$\ell(m) = \ell_{\rm tot}\,\frac{m^{1/3}\,\bar{c}_m}{\displaystyle\sum_{m'} m'^{1/3}\,\bar{c}_{m'}}$$

System size: $N + M + 2$.

---

## 7. Size-Bin Moment Reduction  (Chapter 9, Eqs. bin\_boundaries–delta\_FP\_binned)

Bin boundaries (Eq. bin\_boundaries): first $n_0 = n_{\max}^i$ sizes resolved individually;
remaining follow geometric progression (Eq. bin\_ratio):

$$s_{k+1} = \lfloor r\, s_k \rfloor + 1, \qquad r = 1.5\text{–}2.0$$

Number of grouped bins (Eq. NB):

$$N_B^{\rm group} = \left\lceil \frac{\ln(N/n_{\max}^i)}{\ln r} \right\rceil, \qquad
  N_B = n_{\max}^i + N_B^{\rm group}$$

Bin moments (Eq. bin\_moments\_def):

$$\mu_k^{(0)} = \sum_{n \in \mathcal{B}_k} c_n, \qquad
  \mu_k^{(1)} = \sum_{n \in \mathcal{B}_k} n\,c_n$$

### Bin-integrated rate equations

The bin-moment equations are obtained by summing the corrected master
equation (Eq. ME\_SIA, Section 5) over all sizes $n \in \mathcal{B}_k$.
Every rate process from ME\_SIA must appear in the bin projection.

---

**Zeroth-moment equation** (Eq. mu0\_eq):

$$\frac{d\mu_k^{(0)}}{dt} =
  \underbrace{\sum_{n \in \mathcal{B}_k} G_n}_{\text{production}}$$
$$+ \underbrace{J_{k-1 \to k} - J_{k \to k+1}}_{\text{inter-bin flux (advection in size space)}}$$
$$+ \underbrace{\sum_{n \in \mathcal{B}_k}\bigl[
    \frac{1}{2}\!\!\sum_{n'}\mathcal{K}_{n',n-n'}^{ii}\,c_{n'}\,c_{n-n'}
    - c_n\!\sum_{n'}\mathcal{K}_{n,n'}^{ii}\,c_{n'}\bigr]}_{\text{i–i coalescence (gain + loss)}}$$
$$+ \underbrace{\sum_{n \in \mathcal{B}_k}\bigl[
    \varepsilon_i(n{+}1)\,c_{n+1} - \varepsilon_i(n)\,c_n\,\delta_{n\ge 2}
    \bigr]}_{\text{thermal SIA emission}}$$
$$+ \underbrace{\sum_{n \in \mathcal{B}_k}\bigl[
    \sum_{m'}\mathcal{K}_{m',n+m'}^{vi}\,c_{m',0}\,c_{n+m'}
    - \sum_{m'}\mathcal{K}_{m',n}^{vi}\,c_{m',0}\,c_n
    \bigr]}_{\text{V–I annihilation (gain + loss); } n=1 \text{ uses } \mathcal{K}_{iv}}$$
$$- \underbrace{\sum_{n \in \mathcal{B}_k}\sum_m
    \mathcal{K}_{n,m}^{(\rm cav)}\,c_n\,\bar{c}_m}_{\text{SIA cluster–cavity absorption}}$$
$$- \underbrace{\mathcal{D}_i\,\mu_k^{(0)}}_{\text{fixed sinks}}$$
$$+ \underbrace{\delta_{1 \in \mathcal{B}_k}\sum_{m,\ell}
    \Gamma_{\rm TM}(m,\ell)\,c_{m,\ell}}_{\text{trap mutation (only bin containing } n=1\text{)}}$$

The first-moment equation $d\mu_k^{(1)}/dt$ is the same sum weighted by $n$.

---

### Vacancy and He equations (unchanged)

The vacancy cluster equations and He equations in the bin-moment system
are **identical** to the full per-size equations (Eq. ME\_vac, Eq. ME\_He)
with all corrected rate processes from Section 5:
- P1 recombination at $m=1$ uses $\mathcal{K}_{iv}$
- SIA-induced shrinkage includes both gain ($V_{m+n} \to V_m$) and loss
  ($V_m \to V_{m-n}$) for all mobile $n$
- 1D/3D mixed cavity absorption for $n = 4,\ldots,n_{\max}^i$
- He capture does NOT change the vacancy size class $m$

Only the SIA population is bin-grouped; the vacancy and He populations
retain full per-size resolution.

---

### Closure and flux

Three intra-bin shape functions are available (selected via `shape_function`):

**Piecewise-constant** (`"constant"`, $P=1$, Eq. pc\_ansatz):

$$c_n \approx \frac{\mu_k^{(0)}}{\Delta_k} \quad \text{for } n \in \mathcal{B}_k$$

**Hat-function / linear** (`"linear"`, $P=2$, Eqs. shape\_ansatz–hat\_moment1):

$$c_n = \phi_{k,0}(n)\,\mu_k^{(0)} + \phi_{k,1}(n)\,\mu_k^{(1)}$$

where $\phi_{k,0}(n) = (S_2 - S_1 n)/\Delta$, $\phi_{k,1}(n) = (\Delta_k n - S_1)/\Delta$.

**Log-normal** (`"lognormal"`, $P=3$, Eqs. lognormal\_ansatz–lognormal\_params):

$$c_n \propto \frac{1}{n}\exp\!\left[-\frac{(\ln n - m_k)^2}{2\sigma_k^2}\right]$$

where $\sigma_k^2 = \ln(\mu_k^{(2)}\mu_k^{(0)}/{\mu_k^{(1)}}^2)$,
$m_k = \ln(\mu_k^{(1)}/\mu_k^{(0)}) - \sigma_k^2/2$,
and the amplitude is set so $\sum_{n \in \mathcal{B}_k} c_n = \mu_k^{(0)}$.
Falls back to linear when $\sigma_k^2 \le 0$ (monodisperse limit).

Inter-bin upwind flux (Eq. flux\_upwind).

Conservation diagnostic (Eq. delta\_FP\_binned): $\delta_{\rm FP}^{\rm bin}$.

---

## 8. Post-Processing  (Eqs. swelling\_identity–delta\_He)

Define $S_I(t) = \sum_{n=1}^N n\,c_n(t)$ (SIA inventory) and
$S(t) = \sum_{m=1}^M m\,\bar{c}_m(t)$ (swelling) (Eq. SI\_S).

Swelling identity (Eq. swelling\_identity):

$$S(t) = S_I(t) + \Delta J^d(t)$$

where $\Delta J^d(t) = \int_0^t [Z_i^d\,\rho_d\,a^2 \sum_n n\,\omega_n^{\rm eff}\,c_n - Z_v^d\,\rho_d\,a^2\,\omega_v^{\rm eff}\,c_{1,0}]\,dt'$ is the cumulative net bias flux (Eq. DeltaJd).

Frenkel pair conservation diagnostic (Eq. delta\_FP):

$$\delta_{\rm FP}(t) = \frac{\bigl|S(t) - S_I(t) - \Delta J^d(t)\bigr|}{S(t) + S_I(t) + \Delta J^d(t)}$$

He conservation (exact, Eq. dSHe\_exact):

$$\frac{dS_{\rm He}}{dt} = G_{\rm He} - \mathcal{D}_h\,c_h - \mathcal{D}_v\sum_{m=1}^{m_{\max}^v}\sum_\ell \ell\,c_{m,\ell}$$

He is lost only through (i) free He absorbed at fixed sinks ($\mathcal{D}_h\,c_h$, small)
and (ii) He carried away when mobile vacancy clusters ($m \le m_{\max}^v$) reach fixed sinks
(also small, since mobile clusters carry little He). Immobile voids ($m > m_{\max}^v$) do not
diffuse to sinks and their He is retained.

He conservation diagnostic (Eq. delta\_He):

$$\delta_{\rm He}(t) = \frac{\displaystyle\left|c_h + \sum_{m,\ell}\ell\,c_{m,\ell} - \int_0^t G_{\rm He}\,dt' - \bigl[c_h(0) + \sum_{m,\ell}\ell\,c_{m,\ell}(0)\bigr]\right|}{\displaystyle c_h + \sum_{m,\ell}\ell\,c_{m,\ell}}$$

Both $\delta_{\rm FP}$ and $\delta_{\rm He}$ should remain below $\sim 10^{-6}$; values above $10^{-3}$ indicate a coding error.

---

## 9. Solver Modes

| Mode | Description |
|---|---|
| `full_system` | Full system, SUNDIALS CVODE BDF, dense/band/GMRES linear solver. |
| `active_window` | Two independent sliding windows (SIA + VAC) + OpenMP-parallel RHS. Thread count is auto-picked from `N_eq` (overridable via `OMP_NUM_THREADS`); when OpenMP is unavailable or the auto-pick lands on 1, the same code path runs serial transparently. |

Legacy aliases `cpp_full` → `full_system` and `sliding_OpenMP` → `active_window` are accepted silently for back-compatibility with older configs and saved output directories.

Physics is selected by two orthogonal axes:

| Axis | Values | Meaning |
|---|---|---|
| `equations` | `discrete` \| `bin_moment` | Per-size ODEs (Eqs. ME\_SIA/ME\_vac/ME\_He) vs. Chapter 9 bin-moment grouping. |
| `cascade`   | `fission`  \| `fusion`     | He coupling case + cascade spectrum. `fission` → Case 2 (decoupled, Eq. 175). `fusion` → Case 1 (mean-field, Eq. 174). |

The four combinations map to the canonical `physics_option` strings used internally and on disk (the historical `full_CD` token is preserved in the combined string):

| equations | cascade | physics_option |
|---|---|---|
| `discrete`   | `fission` | `full_CD_fission` |
| `discrete`   | `fusion`  | `full_CD_fusion` |
| `bin_moment` | `fission` | `bin_moment_CD_fission` |
| `bin_moment` | `fusion`  | `bin_moment_CD_fusion` |

`RadClusterSimulation` accepts either the new `(equations, cascade)` pair or the legacy single `physics_option=` kwarg. Legacy `equations='full_CD'` is silently aliased to `'discrete'`. Helpers `make_physics_option(eq, cas)` / `split_physics_option(po)` are exported from `RadCluster_1_0.py_utils`.

### Preconditioner options (for GMRES linear solver, `linsol=2`)

| `prec_type` | Name | Description |
|---|---|---|
| 0 | Jacobi | Diagonal scaling $P = \text{diag}(I - \gamma J)$ (legacy) |
| 1 | **Woodbury** | Bordered-banded SMW preconditioner (default for GMRES) |

The Woodbury preconditioner exploits the Jacobian structure
$J = T + U V^T$ where $T$ is banded (half-bandwidth
$b = \max(2 i_{\rm mobile}, 2 v_{\rm mobile}) + 1$) and $U V^T$ is a
rank-$r$ correction ($r = i_{\rm mobile} + v_{\rm mobile}$) from
mobile species coupling. Uses LAPACK `dgbtrf`/`dgbtrs` for the band
and `dgetrf`/`dgetrs` for the $r \times r$ Schur complement.

**Default selection:** Woodbury is used only for `full_system` mode
(`window_mode=0`) with GMRES. For sliding-window modes (3, 4), the
active system is kept small (50--200 unknowns) so Jacobi+GMRES
converges efficiently; Woodbury's 58-RHS setup cost is
counterproductive at that scale.

Parameters: `prec_type` (0/1), `prec_bw` (auto), `prec_rank` (auto).
See `Docs/Formulation/Jacobian_Preconditioner.tex` for derivation.

---

## 10. Binding Energies  (Eqs. v\_binding–Eb\_blended, Tables 18–19)

Void capillary (Eqs. v\_void\_binding–v\_binding\_cap):

$$E_b^v(m) = E_f^v - A_{\rm void}\bigl[m^{2/3} - (m-1)^{2/3}\bigr],
  \qquad A_{\rm void} = 4\pi\,\gamma_s\,r_0^2$$

Bubble vacancy binding (Eq. EbV\_bubble):

$$E_b^v(m,\ell) = E_b^v(m) + \text{He gas-pressure correction via virial EOS}$$

He virial EOS (Eqs. 64–65):

$$P_{\rm He} V = N_{\rm He}\,k_{\rm B} T\bigl[1 + B_2(N_{\rm He}/V) + B_3(N_{\rm He}/V)^2\bigr]$$

$$B_2 = 1.67\times10^{-29}\ \text{m}^3/\text{atom}, \qquad
  B_3 = 1.84\times10^{-58}\ \text{m}^6/\text{atom}^2$$

He binding to bubble (Eq. EbHe\_blended, Table 19):

$$E_b^{\rm He}(m, \ell) = E_s^{\rm He} + P_{\rm He}(m,\ell)\,\Omega
  - f_{\rm blend}\, A^{\rm He}(\alpha)\exp\!\bigl(-\mu(\ell - \alpha)\bigr)$$

Interstitial loop (Eqs. Eb\_smalln\_fit–Eb\_blended):

$$E_b^{\rm loop}(n) = A_{111}\, n^{-B_{111}} \quad (n \le n_{\rm tr}),
  \quad \text{blended to continuum via } \tanh\!\left(\frac{n - n_{\rm tr}}{\sigma_{\rm tr}}\right)$$

---

## 11. Parameter Sources

| Table | Content |
|---|---|
| Table 2 | Cascade production parameters (fission/fusion) |
| Table 5 | bcc Fe lattice, energetics (updated DFT values) |
| Table 8 | He EOS virial coefficients |
| Table 16 | EUROFER dissolved solute concentrations and trapping energies |
| Table 18 | Atomistic fitting amplitudes $A(m)$ for void correction |
| Table 19 | He binding amplitudes $A^{\rm He}(\alpha)$ |
| Table 25 | Geometric rate constant prefactors |
| Table 26 | Dislocation bias factors and sink parameters |
| Table 27 | Trap mutation barriers $E_{\rm TM}(m, \ell)$ |
| Table 28 | Full rate constant table (3D reactions) |
| Table 29 | Model parameters (solver settings) |
| Table 30 | Rate constant table (1D/mixed reactions) |
