"""Aggregate saved stochastic ensembles into summaries, figures, provenance.

Loads ensemble_*.npz from OUTDIR (written by resume_stochastic.py),
recomputes ODE baselines, then executes the summary/figure/provenance
block of run_stochastic_colonial.main() verbatim against the loaded data.

Usage: python3 aggregate_stochastic.py [outdir]
"""
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_stochastic_colonial as R

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "outputs",
    "20260927_161916_stochastic")

CONFIGS = [
    ("shallowest_inheritance", 0.8, 300, 20260926, 1.0),
    ("shallowest_inheritance", 0.0, 300, 20260927, 1.0),
    ("deepest_inheritance", 0.8, 300, 20260928, 1.0),
    ("deepest_inheritance", 0.0, 300, 20260929, 1.0),
    ("shallowest_inheritance", 0.8, 30, 20260930, 1000.0),
]

results = {}
mf = None
for rule, alpha, n_traj, seed, lam in CONFIGS:
    tag = f"{rule}_a{alpha}" + ("_meanfield" if lam != 1.0 else "")
    fn = os.path.join(OUTDIR, f"ensemble_{tag}.npz")
    if not os.path.exists(fn):
        raise SystemExit(f"missing ensemble file: {fn} -- run resume_stochastic.py first")
    z = np.load(fn, allow_pickle=True)
    ens = dict(years=z["years"], total=z["total"], clusters=z["clusters"],
               n_leaps=int(z["n_leaps"]), rule=str(z["rule"]),
               alpha=float(z["alpha"]), seed=int(z["seed"]),
               lam=float(z["lam"]))
    if (rule, alpha) not in R.ODE:
        print(f"ODE baseline: {rule} alpha={alpha}", flush=True)
        R.run_ode(rule, alpha)
    if lam == 1.0:
        results[(rule, alpha)] = (ens, R.summarize(rule, alpha, ens))
        print(f"summarized: {tag}", flush=True)
    else:
        mf = ens

# Execute main()'s summary/figure/provenance block verbatim, except the
# mean-field rerun (mf is already loaded from disk).
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "run_stochastic_colonial.py")).read().splitlines(keepends=True)
import textwrap
start = next(i for i, l in enumerate(src)
             if "# ---- summary table" in l)
end = next(i for i, l in enumerate(src)
           if "total wall time" in l)
block = "".join(src[start:end])
block = textwrap.dedent(block)
block = block.replace(
    'mf = run_ensemble("shallowest_inheritance", 0.8, 30, 20260930, lam=1000.0)',
    '# mf loaded from disk (resume path)')

ns = dict(R.__dict__)  # module globals: np, pd, N, RULES, ...
ns["plt"] = plt
ns["matplotlib"] = matplotlib
ns.update(dict(results=results, mf=mf, OUTDIR=OUTDIR, CONFIGS=CONFIGS,
               t_start=time.time()))
exec(compile(block, "stochastic_summary_block", "exec"), ns)
print("AGGREGATION COMPLETE", flush=True)
