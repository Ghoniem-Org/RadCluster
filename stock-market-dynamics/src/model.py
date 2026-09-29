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
N_MOM, N_VOL = 5, 3

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


def load_market_driver(path=None, source='panel'):
    """Monthly window-delta driver M(t) = r(t) - r(t-12).

    source='panel' (default): cross-sectional driver from the equal-weighted
    price panel (data/panel_monthly.csv, measured).  This is the mechanically
    correct driver: it is the mean over stocks of the individual window-delta
    r_i(t) - r_i(t-12) that moves each stock's 12m momentum.
    source='index': ^GSPC-based driver (data/spx_monthly.csv, measured).

    M(t) is known at month-end t -> legitimate for 1-month-ahead forecasts.
    """
    import csv as _csv
    if path is None:
        p = os.path.join(DATA, 'panel_monthly.csv' if source == 'panel'
                         else 'spx_monthly.csv')
    else:
        p = path
    col = 'M_cs' if source == 'panel' else 'M'
    out = {}
    with open(p) as f:
        for row in _csv.DictReader(f):
            try:
                out[row['ym']] = float(row[col])
            except (ValueError, TypeError, KeyError):
                pass
    return out


def advect(c, M, k_up, k_down, w_eff=0.15, damp_0=1.0, damp_4=1.0):
    """Donor-cell advection along the momentum axis within each vol lane,
    scaled MECHANICALLY: a window-delta M shifts momentum by M in return
    space, i.e. by M/w_eff bins.  k_up/k_down are O(1) fudge factors
    (k=1 is the mechanical value).

    M > 0 shifts mass toward winner bins, M < 0 toward loser bins.
    Reflecting boundaries: the extreme bins accumulate (this generates
    spikes).  w_eff=0.15 is the typical momentum-bin width (assumed).

    damp_0/damp_4: persistence dampers (<=1) applied to the extreme
    SOURCE bin outflow (bin 0 when M>0, bin 4 when M<0).  They encode
    the depth of spike mass: 1.0 = no impedance (uniform donor-cell).
    """
    out = c.copy()
    k = k_up if M > 0 else k_down
    f = k * abs(M) / w_eff
    if f <= 0:
        return out
    f = min(f, 1.0)
    for v in range(N_VOL):
        if M > 0:
            for m in range(N_MOM - 1):
                i = m * N_VOL + v
                ff = f * (damp_0 if m == 0 else 1.0)
                fl = ff * c[i]
                out[i] -= fl
                out[i + N_VOL] += fl
        else:
            for m in range(1, N_MOM):
                i = m * N_VOL + v
                ff = f * (damp_4 if m == N_MOM - 1 else 1.0)
                fl = ff * c[i]
                out[i] -= fl
                out[i - N_VOL] += fl
    return out

def top_share(c):
    return sum(c[m * N_VOL + v] for m in [N_MOM - 1] for v in range(N_VOL))


def bot_share(c):
    return sum(c[m * N_VOL + v] for m in [0] for v in range(N_VOL))


def step_full(c, T, M, k_up, k_down, h0, h1, s_vec, d, lam1=0.0,
              damp_0=1.0, damp_4=1.0):
    """One monthly step: sticky drift -> autocatalytic herding -> advection.

    lam1: concentration-dependent drift stickiness.  lam = lam1*max(w,l);
    the drift becomes (1-lam)*T + lam*I, so concentrated distributions
    relax more slowly (sticky spikes).  lam1=0 recovers the plain drift.
    damp_0/damp_4: extreme-bin outflow persistence dampers (see advect).
    Returns c_new.
    """
    w, l = top_share(c), bot_share(c)
    lam = min(lam1 * max(w, l), 0.9)
    Teff = (1.0 - lam) * T + lam * np.eye(len(c))
    c = Teff.T @ c
    h_eff = h0 + h1 * (w - l)
    c = c + h_eff * herding_flux(c)
    c = advect(c, M, k_up, k_down, damp_0=damp_0, damp_4=damp_4)
    c = c + s_vec - d * c
    return c / c.sum()


def simulate_full(c0, T, Ms, k_up, k_down, h0, h1, s_rate, d_rate, n,
                  lam1=0.0, p_damp=1.2, tau_build=6.0, gate=0.5):
    """Full simulation with spike-persistence dampers.

    p_damp: power-law exponent for the outflow impedance.  When M reverses
    out of a spike, the extreme-bin outflow fraction is (|M|/build)^p_damp
    (capped at 1), where build is the characteristic |M| that built the
    spike.  Assumed 1.2 (tuned for post-peak decay).
    tau_build: decay time (months) for the build-step memory.  Assumed 6.0.
    gate: the damper only engages when the extreme source bin exceeds this
    share (spike state).  Assumed 0.5.
    """
    s_vec = np.full(N_BIN, s_rate / N_BIN)
    traj = [c0.copy()]
    c = c0.copy()
    build_0, build_4 = 0.15, 0.15  # neutral: one bin width
    bdecay = np.exp(-1.0 / tau_build)
    M_prev = Ms[0] if len(Ms) > 0 else 0.0
    for t in range(n):
        M = Ms[t]
        F = abs(M)
        w, l = top_share(c), bot_share(c)
        # Update build memories (decay + reinforce on inflow)
        build_0 = 0.15 + (build_0 - 0.15) * bdecay
        build_4 = 0.15 + (build_4 - 0.15) * bdecay
        if M < 0:
            build_0 = max(build_0, F)
        elif M > 0:
            build_4 = max(build_4, F)
        # Dampers: engage only on reversal out of a spike
        damp_0, damp_4 = 1.0, 1.0
        if M > 0 and M_prev < 0 and l > gate:
            damp_0 = min((F / build_0) ** p_damp, 1.0)
        elif M < 0 and M_prev > 0 and w > gate:
            damp_4 = min((F / build_4) ** p_damp, 1.0)
        c = step_full(c, T, M, k_up, k_down, h0, h1, s_vec, d_rate, lam1,
                      damp_0, damp_4)
        M_prev = M
        traj.append(c)
    return np.array(traj)
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

def transition_matrix_mneutral(reg, Ms, d0, d1, thresh=0.03):
    """Drift estimated on M-neutral months only (|M(t)| < thresh).

    The common (market-wide) migration is carried by the advection kernel;
    this matrix is the idiosyncratic residual.  The transition t -> t+1 is
    included when |M(t)| < thresh, M(t) being known at month-end t.
    """
    import pandas as pd
    r = reg[(reg['date'] >= d0) & (reg['date'] <= d1)].copy()
    r['ym'] = pd.to_datetime(r['date']).dt.strftime('%Y-%m')
    months = sorted(r['ym'].unique())
    counts = np.zeros((N_BIN, N_BIN))
    n_used = 0
    for a, b in zip(months[:-1], months[1:]):
        if abs(Ms.get(a, 0.0)) >= thresh:
            continue
        da = r[r['ym'] == a].set_index('ticker')['bin']
        db = r[r['ym'] == b].set_index('ticker')['bin']
        common = da.index.intersection(db.index)
        for t in common:
            counts[int(da[t]), int(db[t])] += 1
        n_used += 1
    row = counts.sum(axis=1, keepdims=True)
    T = np.divide(counts, row, out=np.eye(N_BIN), where=row > 0)
    return T, int(counts.sum()), n_used


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
