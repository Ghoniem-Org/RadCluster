"""Rev-5 risk metric 1: ensemble distribution forecast with uncertainty intervals.

M=200 scenarios, 12-month horizon from 2022-12 (end of panel).
Per scenario:
  - Tcalm/Tstress drawn from Dirichlet posteriors. Mean = cap-weighted T from
    params_rev5.json (measured 2000-2007); precision = Kish effective sample
    size per origin bin, N_eff = (sum w)^2/sum(w^2) over cap-weighted
    transitions (estimated; the Dirichlet form is assumed).
  - Market driver M(t): AR(1) fit on M history through 2022-12 with
    bootstrapped residuals (estimated).
  - VIX regime (stress = VIX>30): 2-state Markov chain fit 2000-2022 (measured).
Model: rev-5 calibrated kernels (h1=0, lam1=0; advection + damper active).
Output: per-horizon per-bin 5/50/95 percentiles; top/bottom lane and HHI
intervals. Saves outputs/ensemble_forecast.json.
"""
import os, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, '..', 'outputs')
from model4 import (N_BIN, assign_bins, cap_state, simulate, market_driver,
                    top_share, bot_share)

M_SCEN = 200
HORIZON = 12
SEED = 20260928


def kish_neff(panel, d0='2000-01-01', d1='2007-12-31'):
    """Per-origin-bin Kish effective sample size of cap-weighted transitions."""
    sub = panel[(panel['month'] >= d0) & (panel['month'] <= d1)
                & (panel['bin'] >= 0) & (panel['mcap'] > 0)].copy()
    sub = sub.sort_values(['ticker', 'month'])
    sub['bin_next'] = sub.groupby('ticker')['bin'].shift(-1)
    sub['month_next'] = sub.groupby('ticker')['month'].shift(-1)
    sub['dm'] = (sub['month_next'].dt.year * 12 + sub['month_next'].dt.month
                 - (sub['month'].dt.year * 12 + sub['month'].dt.month))
    tr = sub[(sub['dm'] == 1) & (sub['bin_next'] >= 0)].copy()
    neff = np.zeros(N_BIN)
    for b, g in tr.groupby('bin'):
        w = g['mcap'].values.astype(float)
        w = w[np.isfinite(w) & (w > 0)]
        if len(w) and w.sum() > 0:
            neff[int(b)] = (w.sum() ** 2) / (np.dot(w, w))
    return neff


def fit_ar1_M(M_hist, end_ym='2022-12'):
    yms = sorted(k for k in M_hist if k <= end_ym)
    y = np.array([M_hist[k] for k in yms])
    X, Y = y[:-1], y[1:]
    Xc = X - X.mean()
    phi = float(np.dot(Xc, Y - Y.mean()) / np.dot(Xc, Xc)) if np.dot(Xc, Xc) > 0 else 0.0
    phi = float(np.clip(phi, -0.95, 0.95))
    a = float(Y.mean() - phi * X.mean())
    resid = Y - (a + phi * X)
    return a, phi, resid, y[-1]


def fit_vix_markov(end_ym='2022-12'):
    v = pd.read_csv(os.path.join(DATA, 'raw', 'VIX.csv'), parse_dates=['date'])
    v['ym'] = v['date'].dt.strftime('%Y-%m')
    vx = v.groupby('ym')['close'].mean()
    yms = sorted(k for k in vx.index if k <= end_ym)
    s = np.array([1 if vx[k] > 30.0 else 0 for k in yms])
    Tm = np.zeros((2, 2))
    for i, j in zip(s[:-1], s[1:]):
        Tm[i, j] += 1
    Tm = Tm / Tm.sum(axis=1, keepdims=True)
    return Tm, int(s[-1])


def main():
    rng = np.random.default_rng(SEED)
    panel = pd.read_csv(os.path.join(DATA, 'panel_v2.csv'),
                        parse_dates=['month'], low_memory=False)
    panel['bin'] = assign_bins(panel)
    P = json.load(open(os.path.join(DATA, 'params_rev5.json')))
    Tcalm0 = np.array(P['_Tcalm'])
    Tstress0 = np.array(P['_Tstress'])
    erate, emix = P['_erate'], np.array(P['_emix'])

    neff = kish_neff(panel)
    print('Kish N_eff: median=%.0f min=%.0f max=%.0f' % (
        np.median(neff[neff > 0]), neff[neff > 0].min(), neff.max()), flush=True)

    M_hist = market_driver(panel)
    a, phi, resid, M_last = fit_ar1_M(M_hist)
    print('M AR(1): a=%.4f phi=%.3f resid_sd=%.4f' % (a, phi, resid.std()), flush=True)
    Tvix, s0 = fit_vix_markov()
    print('VIX Markov:\n', np.round(Tvix, 3), ' start=%d' % s0, flush=True)

    months = sorted(panel['month'].unique())
    c0, _ = cap_state(panel[panel['month'] == months[-1]])
    print('forecast origin:', months[-1].date(), ' top_lane=%.3f' % top_share(c0), flush=True)
    horizons = pd.date_range(months[-1] + pd.offsets.MonthEnd(1), periods=HORIZON, freq='M')
    hlabels = [d.strftime('%Y-%m') for d in horizons]

    # Dirichlet posteriors
    def draw_T(T0):
        T = np.zeros_like(T0)
        for i in range(N_BIN):
            alpha = neff[i] * T0[i] + 0.5
            alpha = np.clip(alpha, 0.5, None)
            T[i] = rng.dirichlet(alpha)
        return T

    paths = np.zeros((M_SCEN, HORIZON + 1, N_BIN))
    paths[:, 0, :] = c0
    for s in range(M_SCEN):
        Tc, Ts = draw_T(Tcalm0), draw_T(Tstress0)
        # M path
        Mp = np.zeros(HORIZON)
        m = M_last
        for h in range(HORIZON):
            m = a + phi * m + rng.choice(resid)
            Mp[h] = m
        # regime path
        reg = np.zeros(HORIZON, dtype=int)
        st = s0
        for h in range(HORIZON):
            st = int(rng.choice([0, 1], p=Tvix[st]))
            reg[h] = st
        # simulate with per-month T
        cc = c0.copy()
        traj = [cc.copy()]
        from model4 import step_full
        bdec = np.exp(-1.0 / P['tau_build'])
        build_lo, build_hi = 0.15, 0.15
        M_prev = Mp[0]
        for h in range(HORIZON):
            M = Mp[h]; F = abs(M)
            w, l = top_share(cc), bot_share(cc)
            build_lo = 0.15 + (build_lo - 0.15) * bdec
            build_hi = 0.15 + (build_hi - 0.15) * bdec
            if M < 0:
                build_lo = max(build_lo, F)
            elif M > 0:
                build_hi = max(build_hi, F)
            dlo, dhi = 1.0, 1.0
            if M > 0 and M_prev < 0 and l > P['gate']:
                dlo = min((F / build_lo) ** P['p_damp'], 1.0)
            elif M < 0 and M_prev > 0 and w > P['gate']:
                dhi = min((F / build_hi) ** P['p_damp'], 1.0)
            T = Ts if reg[h] == 1 else Tc
            cc = step_full(cc, T, M, P, damp_lo=dlo, damp_hi=dhi,
                           erate=erate, emix=emix)
            M_prev = M
            traj.append(cc.copy())
        paths[s] = np.array(traj)
        if (s + 1) % 50 == 0:
            print(f'  scenario {s+1}/{M_SCEN}', flush=True)

    # percentiles
    top = paths[:, :, [b for b in range(N_BIN) if b // 8 == 2]].sum(axis=2)
    bot = paths[:, :, [b for b in range(N_BIN) if b // 8 == 0]].sum(axis=2)
    hhi = (paths ** 2).sum(axis=2)

    def pct(x):
        return [float(np.percentile(x, q)) for q in (5, 50, 95)]

    out = {
        'origin': months[-1].strftime('%Y-%m-%d'),
        'horizons': hlabels,
        'M_scenarios': M_SCEN,
        'top_lane': [pct(top[:, h]) for h in range(1, HORIZON + 1)],
        'bot_lane': [pct(bot[:, h]) for h in range(1, HORIZON + 1)],
        'hhi': [pct(hhi[:, h]) for h in range(1, HORIZON + 1)],
        'bins': {str(b): [pct(paths[:, h, b]) for h in range(1, HORIZON + 1)]
                 for b in range(N_BIN)},
        'note': ('Ensemble of 200 Dirichlet-T + AR(1)-M + VIX-Markov scenarios. '
                 'T posterior: Dirichlet(Kish N_eff * T_rev5 + 0.5). '
                 'Labels: T measured 2000-2007; Dirichlet/EffN assumed; '
                 'M AR(1) estimated; VIX Markov measured.'),
    }
    json.dump(out, open(os.path.join(OUT, 'ensemble_forecast.json'), 'w'), indent=1)
    # headline numbers
    t12 = [pct(top[:, 12]), pct(bot[:, 12]), pct(hhi[:, 12])]
    print('12-mo ahead 5/50/95:', flush=True)
    print('  top lane: %.3f / %.3f / %.3f' % tuple(t12[0]), flush=True)
    print('  bot lane: %.3f / %.3f / %.3f' % tuple(t12[1]), flush=True)
    print('  HHI:      %.4f / %.4f / %.4f' % tuple(t12[2]), flush=True)
    print('wrote outputs/ensemble_forecast.json', flush=True)


if __name__ == '__main__':
    main()
