"""Cluster-dynamics model for S&P 500 regime distribution.

State: n(t) in R^15, n_k = fraction of index stocks in regime bin k
        bins = momentum quintile (0..4) x volatility tercile (0..2), bin = mq*3+vt

Master equation (monthly discrete time):
    n(t+1) = T' n(t) + h * H(n(t)) + s - d*n(t)
  T      : empirical drift (Markov transition matrix, measured 2015-2019)
  H_k    : herding flux = n_k * (q_k - qbar(n)), q_k = momentum quintile of bin k
           zero-sum by construction; the bilinear "mating" analog.
           h = 0 -> pure drift (null); h > 0 -> capital chases winners (ratchet)
  s, d   : measured index entry/exit (source/sink), small
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
N_BIN = 15

def bin_mom_q(b):
    return b // 3

def load_regimes():
    reg = pd.read_csv(os.path.join(DATA, 'regimes.csv'), parse_dates=['date'])
    return reg

def monthly_distributions(reg):
    """Pivot to (month x bin) fractions."""
    ct = reg.groupby(['date', 'bin']).size().unstack(fill_value=0)
    for b in range(N_BIN):
        if b not in ct.columns:
            ct[b] = 0
    ct = ct[list(range(N_BIN))]
    return ct.div(ct.sum(axis=1), axis=0).sort_index()

def transition_matrix(reg, d0, d1):
    """Empirical monthly T[i,j] = P(bin i -> bin j) over [d0, d1]."""
    sub = reg[(reg.date >= d0) & (reg.date <= d1)].copy()
    piv = sub.pivot(index='ticker', columns='date', values='bin')
    dates = sorted(piv.columns)
    counts = np.zeros((N_BIN, N_BIN))
    for a, b in zip(dates[:-1], dates[1:]):
        pa = piv[a].dropna().astype(int)
        pb = piv[b].dropna().astype(int)
        both = pa.index.intersection(pb.index)
        for t in both:
            counts[pa[t], pb[t]] += 1
    row = counts.sum(axis=1, keepdims=True)
    T = np.divide(counts, row, out=np.zeros_like(counts), where=row > 0)
    # empty rows -> self-loop
    empty = (row[:, 0] == 0)
    T[empty, :] = 0
    T[empty, np.arange(N_BIN)[empty]] = 1.0
    return T

def entry_exit_rates(reg, d0, d1):
    """Measured monthly entry/exit fractions of the panel."""
    sub = reg[(reg.date >= d0) & (reg.date <= d1)].copy()
    piv = sub.pivot(index='ticker', columns='date', values='bin')
    dates = sorted(piv.columns)
    ins, outs, tot = [], [], []
    prev = set(piv[dates[0]].dropna().index)
    for d in dates[1:]:
        cur = set(piv[d].dropna().index)
        ins.append(len(cur - prev)); outs.append(len(prev - cur)); tot.append(len(prev))
        prev = cur
    tot = np.array(tot, dtype=float)
    return float(np.mean(np.array(ins) / tot)), float(np.mean(np.array(outs) / tot))

def herding_flux(n):
    q = np.array([bin_mom_q(b) for b in range(N_BIN)], dtype=float)
    qbar = float(np.dot(n, q))
    return n * (q - qbar) / 4.0  # scaled: |H_k| <= n_k/2, a perturbation to drift

def step(n, T, h, s_rate, d_rate):
    n = T.T @ n + h * herding_flux(n)
    n = n + s_rate / N_BIN - d_rate * n
    n = np.clip(n, 0, None)
    s = n.sum()
    return n / s if s > 0 else np.full(N_BIN, 1.0 / N_BIN)

def simulate(n0, T, h, s_rate, d_rate, months):
    traj = [n0.copy()]
    n = n0.copy()
    for _ in range(months):
        n = step(n, T, h, s_rate, d_rate)
        traj.append(n.copy())
    return np.array(traj)

def rmse(a, b):
    return float(np.sqrt(np.mean((a - b) ** 2)))

def transition_matrix_by_regime(reg, vix_monthly, d0, d1, thresh=20.0):
    """Separate T for calm (VIX<=thresh) vs stress (VIX>thresh) origin months.

    vix_monthly: Series indexed by 'YYYY-MM' string -> mean VIX.
    """
    sub = reg[(reg.date >= d0) & (reg.date <= d1)].copy()
    piv = sub.pivot(index='ticker', columns='date', values='bin')
    dates = sorted(piv.columns)
    c_calm = np.zeros((N_BIN, N_BIN))
    c_stress = np.zeros((N_BIN, N_BIN))
    for a, b in zip(dates[:-1], dates[1:]):
        ym = a.strftime('%Y-%m') if hasattr(a, 'strftime') else str(a)[:7]
        c = c_stress if vix_monthly.get(ym, 0) > thresh else c_calm
        pa = piv[a].dropna().astype(int)
        pb = piv[b].dropna().astype(int)
        for t in pa.index.intersection(pb.index):
            c[pa[t], pb[t]] += 1

    def norm(counts):
        row = counts.sum(axis=1, keepdims=True)
        T = np.divide(counts, row, out=np.zeros_like(counts), where=row > 0)
        empty = (row[:, 0] == 0)
        T[empty, np.arange(N_BIN)[empty]] = 1.0
        return T
    return norm(c_calm), norm(c_stress), int(c_calm.sum()), int(c_stress.sum())

def simulate_regime(n0, T_calm, T_stress, regimes, h, s_rate, d_rate):
    """regimes: list of bool, True=stress, one per step."""
    traj = [n0.copy()]
    n = n0.copy()
    for r in regimes:
        n = step(n, T_stress if r else T_calm, h, s_rate, d_rate)
        traj.append(n.copy())
    return np.array(traj)

def monthly_vix():
    v = pd.read_csv(os.path.join(DATA, 'vix.csv'), parse_dates=['date'])
    v['ym'] = v['date'].dt.strftime('%Y-%m')
    return v.groupby('ym')['close'].mean()

def calibrate_h(reg, T, s_rate, d_rate, d0, d1, hs=(0, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0)):
    dist = monthly_distributions(reg)
    dist.index = pd.to_datetime(dist.index)
    d0, d1 = pd.Timestamp(d0), pd.Timestamp(d1)
    dates = [d for d in dist.index if d0 <= d <= d1]
    actual = dist.loc[dates].values
    n0 = actual[0]
    out = []
    for h in hs:
        traj = simulate(n0, T, h, s_rate, d_rate, len(dates) - 1)
        out.append((h, rmse(traj, actual)))
    return out, dates, actual
