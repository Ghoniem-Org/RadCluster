"""Resume the stochastic colonial ensembles one config at a time.

Saves each completed ensemble to disk so a killed process only loses the
in-flight config. Idempotent: skips configs whose .npz already exists.
Then run aggregate_stochastic.py to build summaries and figures.

Usage: python3 resume_stochastic.py [outdir]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_stochastic_colonial as R

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "outputs",
    "20260927_161916_stochastic")
os.makedirs(OUTDIR, exist_ok=True)

CONFIGS = [
    ("shallowest_inheritance", 0.8, 300, 20260926, 1.0),
    ("shallowest_inheritance", 0.0, 300, 20260927, 1.0),
    ("deepest_inheritance", 0.8, 300, 20260928, 1.0),
    ("deepest_inheritance", 0.0, 300, 20260929, 1.0),
    ("shallowest_inheritance", 0.8, 30, 20260930, 1000.0),
]

t_start = time.time()
for rule, alpha, n_traj, seed, lam in CONFIGS:
    tag = f"{rule}_a{alpha}" + ("_meanfield" if lam != 1.0 else "")
    fn = os.path.join(OUTDIR, f"ensemble_{tag}.npz")
    if os.path.exists(fn):
        print(f"skip (exists): {tag}", flush=True)
        continue
    print(f"ensemble start: {tag} N={n_traj} lam={lam} seed={seed}", flush=True)
    if (rule, alpha) not in R.ODE:
        R.run_ode(rule, alpha)
        print(f"  ODE done: {rule} alpha={alpha}", flush=True)
    ens = R.run_ensemble(rule, alpha, n_traj, seed, lam=lam)
    np.savez_compressed(
        fn,
        years=ens["years"], total=ens["total"], clusters=ens["clusters"],
        n_leaps=np.array(ens["n_leaps"]),
        rule=np.array(rule), alpha=np.array(alpha),
        seed=np.array(seed), lam=np.array(lam))
    print(f"  done: {ens['n_leaps']} leaps, "
          f"{(time.time() - t_start) / 60:.1f} min elapsed -> {fn}", flush=True)
print("ALL ENSEMBLES COMPLETE", flush=True)
