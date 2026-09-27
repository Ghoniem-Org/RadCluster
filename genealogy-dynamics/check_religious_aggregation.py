#!/usr/bin/env python3
"""Aggregation check: with disaffiliation off, the group-aggregated religious
rhs must equal the two-sex rhs of the aggregated state.

sum_g rhs_religious[(d,sex,g)] == rhs_twosex[(d,sex)] for every (d,sex),
when rho_g = 0 for all groups (random religious pairing: the pairing
P(g1,g2) = pi^M_{g1} pi^F_{g2} factorizes, so the group sums collapse onto
the two-sex weights). Any rho > 0 introduces genuine religious assortment
that reshapes the marginal pedigree-class distribution of mothers -- that
is the point of the religious dimension, not a bug.

This validates the religious mating flux, the S*J birth accumulation, the
source split, and the state indexing in one shot.
"""
import os

os.environ.setdefault("GENEALOGY_MAX_DEPTH", "14")

import numpy as np

from genealogy.populations import N_CLUSTERS, state_index
from genealogy.religion import (G, GROUPS, N_STATE_R, rstate_index, RHO_VEC,
                                immigrant_composition)
from genealogy.kernels import SEX_RATIO_AT_BIRTH, IMMIGRANT_MALE_FRACTION
from genealogy.model import rhs_religious, rhs as rhs_twosex

S = SEX_RATIO_AT_BIRTH
SIGMA = IMMIGRANT_MALE_FRACTION
n = N_CLUSTERS

rng = np.random.default_rng(7)
c = np.abs(rng.normal(0.5, 0.3, N_STATE_R))   # religious state
beta = np.linspace(0.035, 0.022, n)
I, mu, eps, alpha = 1.0, 0.009, 0.00075, 0.6
iot = immigrant_composition(2000.0)

r = rhs_religious(2000.0, c, I=I, beta=beta, alpha=alpha, mu=mu, eps=eps,
                  sigma=SIGMA, s=S, rule="deepest_inheritance",
                  rho=np.zeros(G),   # rho=0: pairing factorizes -> identity
                  iota=iot, delta=0.0)
R = r.reshape(G, 2, n).sum(axis=0)            # (2, n): [sex, d], aggregated
r_agg = np.concatenate([R[0], R[1]])          # two-sex block order

# two-sex IC = group aggregation of the religious IC
C = c.reshape(G, 2, n).sum(axis=0)
c2 = np.concatenate([C[0], C[1]])
r2 = rhs_twosex(2000.0, c2, I=I, beta_pat=beta, beta_mat=beta, alpha=alpha,
                mu=mu, eps=eps, sigma=SIGMA, s=S, rule="deepest_inheritance")

dev = np.max(np.abs(r_agg - r2))
scale = np.max(np.abs(r2))
print(f"max |agg(religious rhs) - twosex rhs| = {dev:.3e}  (scale {scale:.3e})")
assert dev < 1e-9 * max(scale, 1.0), "AGGREGATION CHECK FAILED"
print("aggregation check PASSED")

# also with delta>0 the *total* balance still holds (disaffiliation internal)
from genealogy.kernels import mating_flux_religious
from genealogy.model import total_balance_religious
r3 = rhs_religious(2000.0, c, I=I, beta=beta, alpha=alpha, mu=mu, eps=eps,
                   sigma=SIGMA, s=S, rule="deepest_inheritance", rho=RHO_VEC,
                   iota=iot, delta=0.008)
_, B = mating_flux_religious(c, beta, alpha, RHO_VEC, "deepest_inheritance")
pred = total_balance_religious(I, B, c.sum(), mu, eps)
print(f"total-balance dev with disaffiliation on: {abs(r3.sum()-pred):.3e}")
assert abs(r3.sum() - pred) < 1e-9
print("disaffiliation conservation check PASSED")
