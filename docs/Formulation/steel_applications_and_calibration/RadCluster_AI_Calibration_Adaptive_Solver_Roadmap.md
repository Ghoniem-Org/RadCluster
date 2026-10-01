# RadCluster-AI: Parameter Identification and Adaptive Solver Control for Irradiated-Microstructure Evolution

## 1. Problem Statement

RadCluster simulates microstructure evolution in irradiated materials up to high dose levels, typically on the order of \(\sim 100\) dpa, over a range of temperatures. The model includes several numerical modes, including direct ODE integration and hybrid formulations in which small clusters are treated discretely while larger clusters are represented by bin-moment or moment-based equations.

The code also supports multiple linear algebra and time-integration strategies, including sparse matrix methods, GMRES with preconditioners, and KLU-type sparse direct solvers. In regimes with strong cluster coalescence, helium accumulation, or cavity evolution, the Jacobian/matrix structure may evolve dynamically from nearly arrow-tridiagonal or bordered diagonal forms into more general sparse structures.

The central difficulty is that the physical and numerical problems are coupled:

\[
\text{physical parameters} \rightarrow \text{microstructure evolution} \rightarrow \text{matrix structure} \rightarrow \text{solver performance}.
\]

Thus, RadCluster faces two distinct but interacting AI/ML problems:

\[
\boxed{\text{A. Identify uncertain physical parameters from experimental targets.}}
\]

\[
\boxed{\text{B. Adapt the numerical method to the evolving microstructure and matrix structure.}}
\]

The experimental database contains measurements such as total SIA loop density, cavity density, loop size, cavity size, and size distributions at selected doses and temperatures. These quantities can be treated as calibration targets or ground truth, subject to experimental uncertainty.

---

## 2. Physical-Parameter Identification

Let the critical physical parameter vector be

\[
\boldsymbol{\theta}=(\theta_1,\theta_2,\ldots,\theta_p),
\]

where \(p\approx 10\text{--}12\). These parameters may include migration energies, binding energies, capture efficiencies, loop/cavity bias factors, helium-vacancy association parameters, sink strengths, coalescence coefficients, or reduced cavity-physics parameters.

The RadCluster forward model can be written abstractly as

\[
\mathbf{y}_{\rm sim}=\mathcal{M}(\boldsymbol{\theta},T,d;\boldsymbol{\eta}),
\]

where:

- \(\boldsymbol{\theta}\) = physical parameters,
- \(T\) = irradiation temperature,
- \(d\) = dose,
- \(\boldsymbol{\eta}\) = numerical settings,
- \(\mathbf{y}_{\rm sim}\) = predicted microstructural observables.

Typical target observables are

\[
N_{\rm loop}, \quad \bar{R}_{\rm loop}, \quad N_{\rm cav}, \quad \bar{R}_{\rm cav}, \quad p_{\rm loop}(R), \quad p_{\rm cav}(R).
\]

Because densities may vary by orders of magnitude, errors in number density should usually be measured in log-space:

\[
e_N=\log_{10}N_{\rm sim}-\log_{10}N_{\rm exp}.
\]

A general calibration loss can be written as

\[
\mathcal{L}(\boldsymbol{\theta})=
\sum_j w_j
\left[
\log y_{j,{\rm sim}}(\boldsymbol{\theta})-
\log y_{j,{\rm exp}}
\right]^2
+
w_{\rm dist}D\left(p_{\rm sim}(R),p_{\rm exp}(R)\right),
\]

where \(D\) is a distribution-distance metric.

---

## 3. Fitting Size Distributions

The availability of size distributions is highly valuable. Calibration should not rely only on mean sizes and number densities.

Instead of fitting only

\[
N, \quad \bar{R},
\]

one should fit the full distribution

\[
p_{\rm sim}(R;T,d)
\quad \text{against} \quad
p_{\rm exp}(R;T,d).
\]

Useful metrics include:

### 3.1 Wasserstein Distance

\[
D_W=W_1(p_{\rm sim},p_{\rm exp}).
\]

This is physically meaningful because it measures how much distributional mass must be shifted in size space to transform the simulated distribution into the experimental one.

### 3.2 Moment-Weighted Loss

\[
\mathcal{L}_{\rm dist}
=
w_0\left[\log N_{\rm sim}-\log N_{\rm exp}\right]^2
+
w_1\left[\bar{R}_{\rm sim}-\bar{R}_{\rm exp}\right]^2
+
w_2\left[\sigma_{R,{\rm sim}}-\sigma_{R,{\rm exp}}\right]^2
+
w_3D_W.
\]

### 3.3 Histogram Likelihood

If TEM or other experimental measurements produce binned histograms, the observed counts can be modeled as Poisson variables:

\[
n_k^{\rm exp}\sim \text{Poisson}(n_k^{\rm sim}).
\]

The corresponding log-likelihood is

\[
\log p(\mathbf{n}^{\rm exp}|\boldsymbol{\theta})
=
\sum_k
\left[
 n_k^{\rm exp}\log n_k^{\rm sim}
 -n_k^{\rm sim}
 -\log(n_k^{\rm exp}!)
\right].
\]

This is statistically more defensible than fitting smoothed curves when the experimental data are histogram counts.

---

## 4. Bayesian Calibration

The most scientifically defensible parameter-identification approach is Bayesian calibration.

Assign physical priors to the uncertain parameters:

\[
\theta_i\sim p_i(\theta_i),
\]

with ranges constrained by physics, atomistic calculations, literature data, or expert judgment.

The posterior distribution is

\[
p(\boldsymbol{\theta}|\mathbf{y}_{\rm exp})
\propto
p(\mathbf{y}_{\rm exp}|\boldsymbol{\theta})p(\boldsymbol{\theta}).
\]

This is preferable to seeking only a single best-fit parameter set, because the RadCluster inverse problem is likely to be partially non-identifiable. Many different parameter combinations may reproduce the same loop and cavity density data.

Bayesian calibration can reveal:

- which parameters are strongly constrained by the data,
- which parameters are weakly constrained,
- which parameters are correlated,
- which combinations of parameters are identifiable,
- where additional experiments would be most informative.

For example, one may find that a combined activation quantity such as

\[
E_m+E_b
\]

is identifiable, while \(E_m\) and \(E_b\) separately are not.

### Practical Bayesian/Optimization Methods

Candidate methods include:

| Method | Use |
|---|---|
| Latin hypercube or Sobol sampling | Initial exploration of parameter space |
| Gaussian-process Bayesian optimization | Expensive simulations with limited runs |
| Sequential Monte Carlo | Approximate posterior exploration |
| Ensemble Kalman inversion | Efficient approximate calibration |
| CMA-ES or differential evolution | Robust global optimization |
| Multi-fidelity Bayesian calibration | Combine reduced and full RadCluster models |

A practical sequence is:

\[
\text{Sobol sampling}
\rightarrow
\text{sensitivity screening}
\rightarrow
\text{surrogate modeling}
\rightarrow
\text{Bayesian calibration}.
\]

---

## 5. Sensitivity and Identifiability Analysis

Before final calibration, one should ask:

\[
\boxed{\text{Which of the 10--12 parameters can actually be inferred from the available data?}}
\]

Local sensitivities are

\[
S_i=\frac{\partial \mathbf{y}}{\partial \theta_i}.
\]

Global sensitivity analysis can be performed using Sobol indices:

\[
S_i, \quad S_{ij}, \quad S_i^{\rm total}.
\]

Parameters can then be classified into:

| Parameter class | Meaning |
|---|---|
| Identifiable | Strongly affects measured targets |
| Correlated | Only certain combinations are identifiable |
| Sloppy | Weakly constrained by available data |
| Regime-inactive | Not important in the studied dose/temperature range |

This step prevents overfitting sparse experimental data with too many adjustable parameters.

---

## 6. Surrogate Modeling of RadCluster Outputs

Because RadCluster simulations may be computationally expensive, build a surrogate model:

\[
\widehat{\mathcal{M}}:
(\boldsymbol{\theta},T,d)
\mapsto
\left[N_{\rm loop},\bar{R}_{\rm loop},N_{\rm cav},\bar{R}_{\rm cav},p(R)\right].
\]

Candidate surrogate models include:

| Surrogate | Strength |
|---|---|
| Gaussian process | Good uncertainty estimates; useful for limited data |
| Random forest | Robust for irregular responses |
| Gradient-boosted trees | Strong tabular performance |
| Polynomial chaos expansion | Good if response is smooth |
| Active-subspace model | Useful if only a few parameter combinations matter |
| Deep neural network | Useful only after many simulations |
| Neural operator | Useful for full time histories or distributions |

For a first implementation, Gaussian processes or gradient-boosted trees are preferable to deep neural networks.

The surrogate should also predict computational cost:

\[
t_{\rm CPU}=\mathcal{C}(\boldsymbol{\theta},T,d,\boldsymbol{\eta}).
\]

The cost surrogate is essential because some parameter choices will lead to severe stiffness, strong coalescence, or difficult sparse matrix structures.

---

## 7. Active Learning and Cost-Aware Simulation Design

Rather than sample the 10--12 dimensional parameter space blindly, use active learning to select the next RadCluster simulations.

At each stage, choose

\[
\boldsymbol{\theta}_{n+1}
=
\arg\max_{\boldsymbol{\theta}}
\left[
\text{expected information gain}
-
\lambda\,\text{expected computational cost}
\right].
\]

A useful cost-aware acquisition function is

\[
A(\boldsymbol{\theta})
=
\frac{
\text{expected improvement or uncertainty reduction}
}{
\widehat{t}_{\rm CPU}(\boldsymbol{\theta})
}.
\]

This avoids spending excessive computational effort in parameter regions that are expensive but not informative.

---

## 8. Adaptive Numerical Method Selection

The numerical problem can be formulated as a policy-learning problem.

At each time or dose step, define a numerical state vector:

\[
\mathbf{s}_n=
\left[
N_{\rm species},
N_{\rm nonzero},
\text{bandwidth},
\text{arrow strength},
\text{border size},
\kappa_{\rm est},
\rho_{\rm fill},
r_{\rm coalescence},
r_{\rm emission},
r_{\rm absorption},
\Delta t,
n_{\rm GMRES},
n_{\rm rejected},
\ldots
\right].
\]

The solver action is

\[
a_n\in
\{
\text{direct ODE},
\text{hybrid moments},
\text{GMRES+Jacobi},
\text{GMRES+ILU},
\text{GMRES+block preconditioner},
\text{KLU sparse},
\text{dense fallback},
\text{physics reduction},
\text{adaptive rebinning}
\}.
\]

The objective is to minimize computational cost while preserving accuracy:

\[
J=t_{\rm wall}+\lambda_{\rm fail}I_{\rm fail}+\lambda_{\rm err}\epsilon_{\rm numerical}.
\]

This defines a solver-policy problem:

\[
\pi(\mathbf{s}_n)\rightarrow a_n.
\]

---

## 9. Supervised Solver Classification

The recommended first approach is supervised learning, not reinforcement learning.

Generate a solver-performance database from RadCluster runs:

\[
\{\mathbf{s}_n,a_n^{\rm best},t(a_n),\epsilon(a_n)\}.
\]

At selected steps, run several candidate solvers or preconditioners and record which one is fastest and stable.

Useful features include:

| Feature class | Examples |
|---|---|
| Matrix size | Number of species, number of equations, number of active clusters |
| Sparsity | nnz, fill ratio, bandwidth, border size |
| Stiffness | Jacobian norm, estimated condition number, spectral radius estimate |
| Solver history | GMRES iterations, rejected time steps, KLU factorization time |
| Physics activity | Coalescence rate, absorption rate, emission rate, cavity-growth rate |
| Hybrid structure | Number of discrete clusters, number of moment bins |

Candidate classifiers:

| Model | Advantage |
|---|---|
| Decision tree | Interpretable rules |
| Random forest | Robust nonlinear classification |
| Gradient-boosted trees | Excellent performance on tabular data |
| Logistic classifier | Simple baseline |
| Neural network | Useful later with large data |

A good first choice is a gradient-boosted tree model, followed by extraction of simplified decision-tree rules for implementation inside RadCluster.

---

## 10. Contextual Bandit Solver Selection

Solver selection can also be formulated as a contextual bandit problem:

\[
\text{context}=\mathbf{s}_n,
\quad
\text{action}=a_n,
\quad
\text{reward}=-t_{\rm wall}-\lambda I_{\rm fail}.
\]

This is safer and more practical than full reinforcement learning because it does not require long-horizon credit assignment.

A safe implementation is:

1. Maintain a conservative default solver, such as KLU or a verified sparse direct method.
2. Allow the bandit model to choose among safe candidate solvers.
3. If convergence is poor or the solver fails, revert to the conservative fallback.
4. Record the failure or poor performance as new training data.

This creates an adaptive solver-selection mechanism that improves with use but remains robust.

---

## 11. Matrix-Structure Classification

The Jacobian sparsity pattern can be treated as a graph:

\[
G_J=(V,E),
\]

where vertices are cluster species or moment variables and edges are nonzero couplings.

Graph and sparsity features include:

\[
\text{degree distribution},\quad
\text{bandwidth},\quad
\text{connected components},\quad
\text{border size},\quad
\text{treewidth estimate},\quad
\text{fill-in estimate}.
\]

The matrix can then be classified as:

\[
G_J\mapsto \text{matrix class}.
\]

Possible classes and corresponding solver strategies are:

| Matrix structure | Likely solver strategy |
|---|---|
| Tridiagonal or banded | Banded direct solver |
| Arrow tridiagonal | Sherman--Morrison/Woodbury or bordered solver |
| Bordered diagonal | Block elimination / Schur complement |
| Sparse irregular | KLU / sparse LU |
| Large sparse, mildly stiff | GMRES + ILU/block preconditioner |
| Coalescence-dominated dense coupling | Physics reduction / grouping / moment closure |

This is one of the most promising AI-assisted numerical components because it is structure-aware rather than purely black-box.

---

## 12. Adaptive Discrete/Moment Partitioning

RadCluster already has a hybrid representation:

\[
\text{small clusters: discrete equations},
\]

\[
\text{large clusters: moment or bin-moment equations}.
\]

The adaptive question is where to place the transition between discrete and moment representations.

Let \(n_c\) be a cutoff cluster size or population boundary. Instead of prescribing it manually, define

\[
n_c=\pi(\mathbf{s},T,d,\boldsymbol{\theta}).
\]

Possible criteria are:

\[
\epsilon_{\rm closure}(n_c)<\epsilon_{\rm tol},
\]

\[
\frac{t_{\rm discrete}(n_c)}{t_{\rm moment}(n_c)}>\tau,
\]

\[
\text{distribution smoothness above }n_c\text{ is sufficient}.
\]

Useful features include:

\[
\left[
\frac{\partial C_n}{\partial n},
\frac{\partial^2 C_n}{\partial n^2},
\text{tail mass},
\text{tail variance},
\text{coalescence activity},
\text{gas loading},
\text{cavity pressure regime}
\right].
\]

A practical rule is:

\[
n_c=\min n
\quad \text{such that}\quad
\left|\frac{C_{n+1}-C_n}{C_n}\right|<\epsilon_1
\]

and

\[
\left|\frac{C_{n+2}-2C_{n+1}+C_n}{C_n}\right|<\epsilon_2.
\]

AI can then learn how \(\epsilon_1\), \(\epsilon_2\), and \(n_c\) should vary with the physical and numerical state.

---

## 13. Physics-Constrained Neural Networks

Neural networks can be useful, but they should be used carefully. They should not replace RadCluster initially.

Possible roles include:

### 13.1 Neural Surrogate for Outputs

\[
(\boldsymbol{\theta},T,d)\mapsto \mathbf{y}_{\rm targets}.
\]

This is useful only after enough RadCluster simulations have been accumulated.

### 13.2 Neural Correction Model

Use RadCluster as the baseline and learn only a correction:

\[
\mathbf{y}_{\rm exp}
=
\mathbf{y}_{\rm RadCluster}(\boldsymbol{\theta})
+
\Delta_{\rm ML}(\boldsymbol{\theta},T,d).
\]

This is safer than replacing the physics model.

### 13.3 Latent Distribution Model

Represent a size distribution using a low-dimensional basis:

\[
p(R,t)\approx \sum_{k=1}^m z_k(t)\phi_k(R).
\]

Then learn the evolution of the latent variables \(z_k(t)\), while enforcing physical constraints such as:

\[
N\geq 0,\quad R\geq 0,
\]

and, where applicable,

\[
\text{vacancy conservation},\quad
\text{interstitial conservation},\quad
\text{helium conservation},\quad
\text{nonnegative cluster populations}.
\]

---

## 14. Digital Twin Formulation

RadCluster can be framed as a digital twin for irradiation microstructure evolution.

The microstructural state is

\[
\mathbf{x}(d,T)=
\left[
C_{\rm SIA}(n),
C_{\rm V}(n),
C_{\rm HeV}(m,\ell),
N_{\rm loop},
N_{\rm cav},
p_{\rm cav}(R),
\ldots
\right].
\]

The uncertain parameters are \(\boldsymbol{\theta}\). Experimental observations are

\[
\mathbf{z}_k=H(\mathbf{x}_k)+\boldsymbol{\epsilon}_k.
\]

Possible data-assimilation approaches include:

| Method | Use |
|---|---|
| Ensemble Kalman filter | Sequential state/parameter updates |
| Ensemble smoother | Calibration using all dose/temperature data |
| Particle filter | Nonlinear/non-Gaussian inference, but expensive |
| Variational data assimilation | Optimize full trajectory |

This is particularly attractive when measurements are available at multiple doses, for example

\[
0.1,\ 1,\ 10,\ 50,\ 100\ \text{dpa}.
\]

The calibration then fits the evolution trajectory, not only isolated endpoint data.

---

## 15. Multi-Objective Optimization

The RadCluster calibration problem is naturally multi-objective.

Targets may include:

\[
N_{\rm loop},\quad
R_{\rm loop},\quad
N_{\rm cav},\quad
R_{\rm cav},\quad
p(R),\quad
t_{\rm CPU}.
\]

A multi-objective formulation is

\[
\min_{\boldsymbol{\theta},\boldsymbol{\eta}}
\left[
\mathcal{L}_{\rm loop},
\mathcal{L}_{\rm cavity},
\mathcal{L}_{\rm distribution},
t_{\rm CPU}
\right].
\]

This produces a Pareto set:

\[
\mathcal{P}=\left\{
\boldsymbol{\theta}:
\text{no other parameter set improves all objectives}
\right\}.
\]

This is useful when one parameter set fits loop density well but fails cavity size, while another does the opposite.

---

## 16. Joint Physical and Numerical Optimization

The full problem may be written as

\[
\min_{\boldsymbol{\theta},\boldsymbol{\eta}}
\mathcal{L}
\left[
\mathcal{M}(\boldsymbol{\theta};\boldsymbol{\eta}),
\mathbf{y}_{\rm exp}
\right]
+
\lambda t_{\rm CPU}(\boldsymbol{\theta},\boldsymbol{\eta})
+
\mu \mathcal{R}(\boldsymbol{\theta}),
\]

where:

- \(\boldsymbol{\theta}\) are physical parameters,
- \(\boldsymbol{\eta}\) are numerical settings,
- \(\mathcal{R}(\boldsymbol{\theta})\) penalizes unphysical parameter combinations.

However, physical inference should not be corrupted by numerical approximation. Therefore, a two-tier strategy is recommended:

\[
\boxed{\text{Tier 1: Calibrate physics using verified numerical tolerances.}}
\]

\[
\boxed{\text{Tier 2: Learn the fastest solver policy that preserves Tier-1 accuracy.}}
\]

The optimizer should not be allowed to choose a faster solver if that solver changes the inferred physical conclusions.

---

## 17. Recommended RadCluster-AI Architecture

A practical implementation can be organized into four modules.

### 17.1 Experiment/Target Database

Store experimental measurements and uncertainties:

\[
T,\quad d,\quad \dot{d},\quad \text{alloy},\quad \text{measurement type},\quad N,\quad \bar{R},\quad p(R),\quad \sigma_{\rm exp}.
\]

Uncertainty estimates are essential:

\[
\sigma_N,\quad \sigma_R,\quad \sigma_{p(R)}.
\]

Without experimental uncertainties, calibration results can be misleading.

### 17.2 Simulation Campaign Manager

Given a parameter set \(\boldsymbol{\theta}\), the campaign manager runs RadCluster and records:

\[
\mathbf{y}_{\rm sim},\quad
t_{\rm wall},\quad
\text{solver trace},\quad
\text{matrix trace}.
\]

The trace should include:

\[
N_{\rm eq}(t),\quad
\text{nnz}(t),\quad
\text{bandwidth}(t),\quad
\text{GMRES iterations}(t),\quad
\text{KLU factorization time}(t),\quad
\Delta t(t),\quad
\text{coalescence events}(t).
\]

The numerical trace is as valuable as the physical output.

### 17.3 Bayesian Calibration Engine

The calibration engine estimates

\[
p(\boldsymbol{\theta}|\mathbf{y}_{\rm exp}).
\]

Recommended sequence:

1. Sobol or Latin hypercube sampling.
2. Global sensitivity analysis.
3. Surrogate construction.
4. Bayesian optimization or ensemble inversion.
5. Posterior uncertainty quantification.

### 17.4 Adaptive Solver Policy

Use matrix and solver trace data to learn

\[
\pi:\mathbf{s}_n\mapsto a_n.
\]

The first version can be rule-based:

```text
if matrix is banded and coalescence is inactive:
    use banded/direct structure solver
elif arrow or bordered structure is detected:
    use bordered/block solver or Woodbury correction
elif sparse irregular and fill is modest:
    use KLU
elif GMRES iteration count is below threshold:
    use GMRES + current preconditioner
else:
    switch to KLU or physics reduction
```

Later, the rule-based logic can be replaced or augmented with learned classifiers or contextual bandit models.

---

## 18. Recommended Initial Workflow

The most practical implementation sequence is:

### Step 1: Select the Critical Physical Parameters

Choose the 10--12 most important uncertain parameters and assign physical bounds:

\[
\theta_i^{\min}\leq \theta_i\leq \theta_i^{\max}.
\]

Avoid including many weak or poorly constrained parameters at the first stage.

### Step 2: Define the Target Vector

For each experiment, define

\[
\mathbf{y}=\left[
\log N_{\rm SIA},
\log \bar{R}_{\rm SIA},
\log N_{\rm cav},
\log \bar{R}_{\rm cav},
\text{distribution metrics}
\right].
\]

### Step 3: Run an Initial Sobol Campaign

Use about

\[
N_{\rm runs}\sim 10p \text{ to } 50p,
\]

where \(p=10\text{--}12\). This gives roughly 100--600 initial simulations, depending on computational cost.

### Step 4: Record Solver and Matrix Traces

For every simulation, save both physical outputs and numerical metadata.

### Step 5: Perform Sensitivity Analysis

Identify which parameters influence which observables.

### Step 6: Build Two Surrogates

Physics surrogate:

\[
\widehat{\mathcal{M}}(\boldsymbol{\theta},T,d).
\]

Cost surrogate:

\[
\widehat{\mathcal{C}}(\boldsymbol{\theta},T,d,\boldsymbol{\eta}).
\]

### Step 7: Calibrate Physical Parameters

Use Bayesian optimization, ensemble Kalman inversion, or surrogate-assisted posterior sampling to infer plausible physical parameter sets.

### Step 8: Train Adaptive Solver Selector

Use the accumulated solver traces to train a solver classifier or contextual bandit policy.

---

## 19. AI Methods Suitable for RadCluster

### 19.1 For Physical Parameter Identification

1. Bayesian calibration
2. Gaussian-process surrogate modeling
3. Active learning / Bayesian experimental design
4. Global sensitivity analysis
5. Multi-fidelity modeling
6. Multi-objective optimization
7. Ensemble Kalman inversion
8. Physics-informed surrogate correction
9. Distribution-aware inverse modeling

### 19.2 For Adaptive Numerical Methods

1. Contextual bandit solver selection
2. Supervised solver-policy classification
3. Matrix-graph classification
4. Cost-aware adaptive time integration
5. Adaptive discrete/moment partitioning
6. Learned preconditioner selection
7. Runtime anomaly detection
8. Reinforcement learning for solver control, but only after simpler methods are successful

---

## 20. Key Warnings

### 20.1 Do Not Overfit Sparse Data

With only a few doses and temperatures, the inverse problem is likely underdetermined.

Expect parameter sloppiness:

\[
\boldsymbol{\theta}_1\neq \boldsymbol{\theta}_2
\quad \text{but} \quad
\mathcal{M}(\boldsymbol{\theta}_1)
\approx
\mathcal{M}(\boldsymbol{\theta}_2).
\]

The best result is not necessarily a single best parameter set:

\[
\boldsymbol{\theta}_{\rm best}.
\]

A better result is a posterior distribution:

\[
p(\boldsymbol{\theta}|\text{data}).
\]

### 20.2 Do Not Let Numerical Error Masquerade as Physics

Solver settings must not alter inferred physical conclusions. Calibration should be done with verified tolerances first. Adaptive solver learning should then be trained to reproduce the verified solution at lower cost.

### 20.3 Do Not Replace RadCluster with a Neural Network Too Early

The most defensible AI role is not model replacement. It is:

\[
\boxed{\text{calibration + uncertainty quantification + adaptive solver intelligence}.}
\]

---

## 21. Recommended Conceptual Framing

A strong conceptual description is:

\[
\boxed{
\text{Physics-constrained, uncertainty-aware, adaptive simulation intelligence for irradiated microstructure evolution.}
}
\]

The AI system supports RadCluster by performing three functions:

\[
\boxed{\text{1. Inferring uncertain physical parameters from sparse irradiation data.}}
\]

\[
\boxed{\text{2. Quantifying uncertainty, sensitivity, and identifiability.}}
\]

\[
\boxed{\text{3. Learning when to switch numerical representations, solvers, and preconditioners.}}
\]

A suitable project or paper title would be:

**RadCluster-AI: Bayesian Calibration and Adaptive Solver Control for Cluster-Dynamics Simulation of Irradiated Microstructures**

or, more compactly:

**AI-Assisted RadCluster for Uncertainty-Aware Microstructure Evolution and Adaptive Numerical Acceleration**

---

## 22. Bottom-Line Recommendation

The most robust RadCluster-AI strategy is not a single black-box neural network. It is a layered architecture:

\[
\boxed{
\text{Bayesian parameter calibration}
+
\text{surrogate modeling}
+
\text{active learning}
+
\text{adaptive solver policy learning}.
}
\]

The immediate implementation should begin with:

1. A carefully defined 10--12 parameter vector with physical bounds.
2. A target database including densities, mean sizes, and size distributions with uncertainties.
3. A Sobol or Latin-hypercube simulation campaign.
4. Full recording of physical outputs, matrix structure, and solver traces.
5. Sensitivity and identifiability analysis.
6. Physics and cost surrogates.
7. Bayesian calibration.
8. Supervised or bandit-based solver adaptation.

This approach preserves the mechanistic value of RadCluster while adding AI where it is most useful: uncertainty quantification, parameter inference, cost-aware simulation design, and adaptive numerical acceleration.
