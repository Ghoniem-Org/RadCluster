"""Rev-4 model: 24-bin capital-weighted state.

Bins: 3 momentum x 2 vol x 2 turnover x 2 size = 24.
  b = ((m*2 + v)*2 + t)*2 + s,  m in 0..2, v,t,s in 0..1
State c in R^24: cap-weighted bin shares (measured from PIT shares).
Kernels (ported from rev-3, generalized):
  1. sticky drift:      T_eff = (1-lam) T + lam I, lam = lam1 * max(w,l)
  2. autocatalytic herding (signed): h_eff = h0 + h1*(w-l), flux along mom axis
  3. fractional advection: exact conservative remap of momentum bins under
     return-space shift delta = kappa*M (unequal bin widths handled exactly)
  4. spike-persistence damper on extreme-lane outflow at reversal (assumed params)
Source/sink: measured monthly entry cap-share blended with trailing entry-bin mix.
Labels: measured / estimated / calibrated / assumed per parameter (see PARAMS).
"""
import numpy as np
import pandas as pd

N_MOM, N_VOL, N_TO, N_SZ = 3, 2, 2, 2
N_BIN = N_MOM * N_VOL * N_TO * N_SZ

MOM_EDGES = np.array([-1.0, 0.0, 0.30, 3.0])     # fixed absolute (measured choice)
VOL_EDGES = np.array([0.0, 0.28, 5.0])            # fixed absolute
TO_EDGES = np.array([0.0, 0.008, 10.0])           # fixed absolute turnover
SZ_EDGES = np.array([0.0, 1e10, 1e16])            # $10B nominal (measured choice)

def bin_of(m, v, t, s):
    return ((m * 2 + v) * 2 + t) * 2 + s

def coords(b):
    s = b % 2; t = (b // 2) % 2; v = (b // 4) % 2; m = b // 8
    return m, v, t, s

def mom_of(b):
    return b // 8

def assign_bins(df):
    """df needs mom, vol, turnover, mcap. Returns bin array (-1 if unbinnable)."""
    m = np.digitize(df['mom'].values, MOM_EDGES) - 1
    v = np.digitize(df['vol'].values, VOL_EDGES) - 1
    t = np.digitize(df['turnover'].values, TO_EDGES) - 1
    s = np.digitize(df['mcap'].values, SZ_EDGES) - 1
    ok = (df['mom'].notna() & df['vol'].notna() & df['turnover'].notna()
          & df['mcap'].notna()).values
    m = np.clip(m, 0, 2); v = np.clip(v, 0, 1); t = np.clip(t, 0, 1); s = np.clip(s, 0, 1)
    b = ((m * 2 + v) * 2 + t) * 2 + s
    return np.where(ok, b, -1)

# ---------------- state ----------------
def cap_state(month_df):
    """Cap-weighted 24-bin shares + coverage diagnostics."""
    d = month_df[(month_df['bin'] >= 0) & (month_df['mcap'] > 0)].copy()
    tot = d['mcap'].sum()
    c = np.zeros(N_BIN)
    if tot > 0:
        for b, g in d.groupby('bin'):
            c[int(b)] = g['mcap'].sum() / tot
    return c, {'n_bin': len(d), 'cap_cov': float(tot),
               'n_members': int((month_df['has_price'] == 1).sum())}

# ---------------- measured drift ----------------
def transition_matrix(panel, d0, d1, vix=None, stress_thresh=30.0, regime=None):
    """Cap-weighted T[i,j] pooled over [d0,d1].
    regime: None (pooled), 'calm', or 'stress' (by VIX at origin month)."""
    sub = panel[(panel['month'] >= d0) & (panel['month'] <= d1)
                & (panel['bin'] >= 0) & (panel['mcap'] > 0)].copy()
    sub = sub.sort_values(['ticker', 'month'])
    sub['bin_next'] = sub.groupby('ticker')['bin'].shift(-1)
    sub['mcap_next'] = sub.groupby('ticker')['mcap'].shift(-1)
    sub['month_next'] = sub.groupby('ticker')['month'].shift(-1)
    # consecutive months only
    sub['dm'] = (sub['month_next'].dt.year * 12 + sub['month_next'].dt.month
                 - (sub['month'].dt.year * 12 + sub['month'].dt.month))
    tr = sub[(sub['dm'] == 1) & (sub['bin_next'] >= 0)].copy()
    if regime is not None and vix is not None:
        tr['ym'] = tr['month'].dt.strftime('%Y-%m')
        tr['stress'] = tr['ym'].map(vix).fillna(0) > stress_thresh
        tr = tr[tr['stress'] == (regime == 'stress')]
    counts = np.zeros((N_BIN, N_BIN))
    for (a, b), g in tr.groupby(['bin', 'bin_next']):
        counts[int(a), int(b)] += g['mcap'].sum()
    row = counts.sum(axis=1, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        T = np.divide(counts, row, out=np.zeros_like(counts), where=row > 0)
    empty = (row[:, 0] == 0)
    T[empty, np.arange(N_BIN)[empty]] = 1.0
    return T, float(counts.sum()), int(tr['month'].nunique())

def entry_stats(panel, d0, d1):
    """Measured monthly entry cap-share and entry-bin mix over [d0,d1]."""
    sub = panel[(panel['month'] >= d0) & (panel['month'] <= d1)].copy()
    months = sorted(sub['month'].unique())
    erates, mixes = [], []
    prev = None
    for m in months:
        d = sub[sub['month'] == m]
        cur = set(d[d['has_price'] == 1]['ticker'])
        if prev is not None:
            new = cur - prev
            if cur:
                nd = d[d['ticker'].isin(new) & (d['bin'] >= 0) & (d['mcap'] > 0)]
                tot = d[(d['bin'] >= 0) & (d['mcap'] > 0)]['mcap'].sum()
                erates.append(nd['mcap'].sum() / tot if tot > 0 else 0.0)
                mix = np.zeros(N_BIN)
                for b, g in nd.groupby('bin'):
                    mix[int(b)] = g['mcap'].sum()
                if mix.sum() > 0:
                    mixes.append(mix / mix.sum())
        prev = cur
    erate = float(np.mean(erates)) if erates else 0.0
    emix = np.mean(mixes, axis=0) if mixes else np.full(N_BIN, 1/N_BIN)
    return erate, emix / emix.sum()

# ---------------- kernels ----------------
def herding_flux(c):
    q = np.array([mom_of(b) for b in range(N_BIN)], dtype=float)
    qbar = float(np.dot(c, q))
    return c * (q - qbar) / 2.0

def top_share(c):
    return float(c[[b for b in range(N_BIN) if mom_of(b) == 2]].sum())

def bot_share(c):
    return float(c[[b for b in range(N_BIN) if mom_of(b) == 0]].sum())

def advect_matrix(delta):
    """Exact conservative 3x3 momentum remap for return-space shift delta.
    Uniform within-bin prior (mechanical). Reflecting at [-1, 3]."""
    e = MOM_EDGES
    A = np.zeros((3, 3))  # A[dest, src]
    for i in range(3):
        lo, hi = e[i] + delta, e[i + 1] + delta
        lo_c, hi_c = max(lo, e[0]), min(hi, e[-1])
        w = e[i + 1] - e[i]
        if hi_c <= lo_c or w <= 0:
            A[i, i] = 1.0
            continue
        for j in range(3):
            ov = max(0.0, min(hi_c, e[j + 1]) - max(lo_c, e[j]))
            A[j, i] = ov / w
        s = A[:, i].sum()
        A[:, i] /= s if s > 0 else 1.0
    return A

def advect(c, M, kappa=1.0, damp_lo=1.0, damp_hi=1.0):
    """Apply fractional advection within each (v,t,s) lane.
    damp_lo/damp_hi: outflow impedance for momentum lane 0 / 2 (<=1)."""
    A = advect_matrix(kappa * M)
    out = np.zeros_like(c)
    for v in range(N_VOL):
        for tt in range(N_TO):
            for s in range(N_SZ):
                lane = [bin_of(m, v, tt, s) for m in range(3)]
                cl = c[lane]
                if M != 0:
                    # damper: scale outflow from extreme lanes
                    B = A.copy()
                    if damp_lo < 1.0:  # lane m=0
                        keep = B[0, 0] + (1 - damp_lo) * (1 - B[0, 0])
                        B[0, 0] = keep
                        B[1:, 0] *= damp_lo
                    if damp_hi < 1.0:  # lane m=2
                        keep = B[2, 2] + (1 - damp_hi) * (1 - B[2, 2])
                        B[2, 2] = keep
                        B[:2, 2] *= damp_hi
                    cl = B @ cl
                out[lane] = cl
    return out

def step_full(c, T, M, P, damp_lo=1.0, damp_hi=1.0, erate=0.0, emix=None):
    """One monthly step. P: dict of params (kappa,h0,h1,lam1)."""
    w, l = top_share(c), bot_share(c)
    lam = min(P['lam1'] * max(w, l), 0.9)
    Teff = (1.0 - lam) * T + lam * np.eye(N_BIN)
    c = Teff.T @ c
    h_eff = P['h0'] + P['h1'] * (w - l)
    c = c + h_eff * herding_flux(c)
    c = advect(c, M, kappa=P['kappa'], damp_lo=damp_lo, damp_hi=damp_hi)
    if erate > 0 and emix is not None:
        c = (1 - erate) * c + erate * emix
    c = np.clip(c, 0, None)
    s = c.sum()
    return c / s if s > 0 else np.full(N_BIN, 1.0 / N_BIN)

def simulate(c0, T, Ms, P, n, erate=0.0, emix=None,
             p_damp=1.2, tau_build=6.0, gate=0.5):
    """Full simulation with spike-persistence damper (rev-3 logic, 24-bin)."""
    traj = [c0.copy()]
    c = c0.copy()
    bdec = np.exp(-1.0 / tau_build)
    build_lo, build_hi = 0.15, 0.15
    M_prev = Ms[0] if len(Ms) else 0.0
    for t in range(n):
        M = Ms[t]
        F = abs(M)
        w, l = top_share(c), bot_share(c)
        build_lo = 0.15 + (build_lo - 0.15) * bdec
        build_hi = 0.15 + (build_hi - 0.15) * bdec
        if M < 0:
            build_lo = max(build_lo, F)
        elif M > 0:
            build_hi = max(build_hi, F)
        dlo, dhi = 1.0, 1.0
        if M > 0 and M_prev < 0 and l > gate:
            dlo = min((F / build_lo) ** p_damp, 1.0)
        elif M < 0 and M_prev > 0 and w > gate:
            dhi = min((F / build_hi) ** p_damp, 1.0)
        c = step_full(c, T, M, P, damp_lo=dlo, damp_hi=dhi,
                      erate=erate, emix=emix)
        M_prev = M
        traj.append(c)
    return np.array(traj)

# ---------------- market driver + forecast ----------------
def market_driver(panel):
    """M(t) = cross-sectional mean of per-stock mom change, known at month-end t.
    Measured."""
    d = panel[(panel['has_price'] == 1) & panel['mom'].notna()].copy()
    d = d.sort_values(['ticker', 'month'])
    d['mom_prev'] = d.groupby('ticker')['mom'].shift(1)
    d['dm'] = (d['month'].dt.year * 12 + d['month'].dt.month
               - (d.groupby('ticker')['month'].shift(1).dt.year * 12
                  + d.groupby('ticker')['month'].shift(1).dt.month))
    d = d[d['dm'] == 1]
    d['dmom'] = d['mom'] - d['mom_prev']
    M = d.groupby(d['month'].dt.strftime('%Y-%m'))['dmom'].mean()
    return M.to_dict()

def mech_forecast(spy_monthly_ret, t_ym):
    """Mechanical M(t+1) nowcast from the momentum definition.
    mom_i(t+1)-mom_i(t) ~= r_i(t)-r_i(t-12), so M(t+1) ~= R(t)-R(t-12),
    both known at month-end t. Zero lookahead. Measured (not fitted).
    spy_monthly_ret: dict ym -> SPY monthly total return.
    Returns (Mhat, info)."""
    yms = sorted(spy_monthly_ret.keys())
    if t_ym not in yms:
        return 0.0, {'method': 'mech', 'note': 'ym missing'}
    j = yms.index(t_ym)
    if j < 12:
        return 0.0, {'method': 'mech', 'note': 'insufficient history'}
    Mhat = float(spy_monthly_ret[t_ym] - spy_monthly_ret[yms[j - 12]])
    return Mhat, {'method': 'mech'}

def ar1_forecast(M_hist, t_ym):
    """Expanding-window AR(1) forecast of M at t+1 using M <= t. No lookahead."""
    yms = sorted(k for k in M_hist if k <= t_ym)
    if len(yms) < 12:
        return 0.0, {}
    y = np.array([M_hist[k] for k in yms])
    X = y[:-1]
    Y = y[1:]
    Xc = X - X.mean()
    phi = float(np.dot(Xc, Y - Y.mean()) / np.dot(Xc, Xc)) if np.dot(Xc, Xc) > 0 else 0.0
    phi = float(np.clip(phi, -0.95, 0.95))
    a = float(Y.mean() - phi * X.mean())
    return float(a + phi * y[-1]), {'a': a, 'phi': phi, 'n': len(y)}

PARAM_LABELS = {
    'kappa': 'mechanical (=1: a pp of window-delta is a pp of momentum)',
    'h0': 'assumed 0 (no baseline herding)',
    'h1': 'calibrated on 2000-2007 conditional hindcast',
    'lam1': 'calibrated on 2000-2007 conditional hindcast',
    'p_damp': 'assumed/tuned on 2000-2007 post-peak decay',
    'tau_build': 'assumed 6 months',
    'gate': 'assumed 0.5 extreme-lane share',
    'gamma': 'assumed 2.0 (tilt aggressiveness)',
    'tcost': 'assumed 5 bps one-way',
}

def rmse(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return float(np.sqrt(np.mean((a - b) ** 2)))
