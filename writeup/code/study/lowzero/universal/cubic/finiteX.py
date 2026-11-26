# Discovery: finite-X models (no fitted parameter): eps' from the Roberts/Taniguchi-Thorne secondary term
import os, sys, json, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mpmath as mp
pc = json.load(open(os.path.join(HERE, 'pred_cache.json'))); edges = np.array(pc['edges']); h = np.diff(edges)
kap = float(4 * mp.zeta(mp.mpf(1) / 3) / (5 * mp.gamma(mp.mpf(2) / 3)**3 * mp.zeta(mp.mpf(5) / 3))); z3 = float(mp.zeta(3))
KC = {-1: 4 * math.sqrt(3) * z3 * kap, 1: 12 * z3 * kap}
def epsp(sgn, X): return KC[sgn] * (10 / 33) * (2**(11 / 6) - 1) * X**(-1 / 6)    # ln p-weighted, dyadic [X, 2X)
cTT = {-1: -1.02494734415798, 1: -1.77526087516441}
ym13 = 1.5 * np.diff(edges**(2 / 3)) / h
def jack(num, den, P):
    n = len(num); N = num.sum(0); D = den.sum(0); m = N / D; mj = (N - num) / (D - den); dm = mj - mj.mean(0)
    cov = (n - 1) / n * dm.T @ dm; sig = np.sqrt(np.diag(cov)); w = 1 / sig**2
    s = np.sum(P * m * w) / np.sum(P * P * w); sj = (mj * (P * w)).sum(1) / np.sum(P * P * w)
    return m, cov, sig, s, math.sqrt((n - 1) / n * np.sum((sj - sj.mean())**2))
def coarse(x, g=4): return x.reshape(x.shape[:-1] + (x.shape[-1] // g, g)).sum(-1)
for sgn in (-1, 1):
    k = '%+d' % sgn; z = np.load(os.path.join(HERE, 'cubic_disc_%s.npz' % k)); num, den, D = z['num'], z['den'], np.abs(z['D'])
    e = epsp(sgn, 5000.0); P0 = np.array(pc[k]['pred'])
    models = {'M0 limit': P0, 'M1 limit + TT share': P0 + (1 / (1 + e) - 1) * cTT[sgn] * ym13, 'M2 limit/(1+eps)': P0 / (1 + e)}
    print('D %s 0: eps\' = %.4f, 1/(1+eps\') = %.3f' % ('<' if sgn < 0 else '>', e, 1 / (1 + e)))
    for sub, sel in (('all', np.ones(len(D), bool)), ('|D|<7500', D < 7500), ('|D|>=7500', D >= 7500)):
        for name, P in models.items():
            m, cov, sig, s, se = jack(num[sel], den[sel], P)
            mc, covc, _, _, _ = jack(coarse(num[sel]), coarse(den[sel]), coarse(P) / 4)
            r = mc - coarse(P) / 4; c2 = float(r @ np.linalg.solve(covc, r)); d2 = float(np.sum((m - P)**2 / sig**2))
            print('   %-9s %-20s s = %.3f +- %.3f  chi2 coarse %6.1f  chi2 fine diag %6.1f' % (sub, name, s, se, c2, d2))
        if sub == 'all':
            for name in ('M1 limit + TT share',):
                P = models[name]; print('      first 6 bins: data', np.round(m[:6], 2), ' M0', np.round(P0[:6], 2), ' M1', np.round(P[:6], 2))
    # s by y-range for M0
    for lo, hi in ((0, 0.2), (0.2, 0.35), (0.35, 1.0), (1.0, 2.0), (2.0, 3.0)):
        sel = (edges[:-1] >= lo - 1e-9) & (edges[1:] <= hi + 1e-9)
        m, cov, sig, s, se = jack(num[:, sel], den[:, sel], P0[sel])
        print('   M0 amplitude on y in (%.2f, %.2f]: %.3f +- %.3f' % (lo, hi, s, se))
