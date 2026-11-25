"""Stress tests of the round-146 reading (own code): splits, trimming, height bands, round 142's heights."""
import os, numpy as np, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from myfit import load, design, ols, PR
from altlib import clo, factor, first
from mpmath import polylog
R = os.path.dirname(os.path.abspath(__file__))+'/../../../../'
LI2 = {p: float(polylog(2, 1.0 / p)) for p in PR}
A123 = {2: (1.2939, -3.30, 5.4), 3: (1.0715, -3.06, 3.5), 5: (0.7984, -2.67, 2.9), 7: (0.6462, -2.38, 2.6),
        11: (0.4789, -2.04, 2.7), 13: (0.4270, -1.84, 1.9)}   # derivation.md sect. 5 compact form


def compact(p, L):
    if p in A123:
        a1, a2, a3 = A123[p]
        return -0.5 * LI2[p] + a1 / L + a2 / L ** 2 + a3 / L ** 3
    return clo(p, int(round(L * 1000)))[0]


def fit(rows, yfun, preds, nboot=1000, seed=1):
    v = np.array([yfun(r) for r in rows])
    X, nn = design(rows)
    b = ols(X, v)[nn:]
    rng = np.random.default_rng(seed); n = len(rows); B = []
    for _ in range(nboot):
        ix = rng.integers(0, n, n); Xb = X[ix].copy(); Xb[:, 1:] -= Xb[:, 1:].mean(0)
        B.append(ols(Xb, v[ix])[nn:])
    B = np.array(B)
    P = {k: ols(X, np.array(pv))[nn:] for k, pv in preds.items()}
    return b, B, P


def amp(bb, z0, one, ps=slice(0, 3)):
    return np.sum((bb[ps] - z0[ps]) * (one[ps] - z0[ps])) / np.sum((one[ps] - z0[ps]) ** 2)


def report(name, rows, yfun, predfun, nboot=1000):
    preds = {k: [f(r) for r in rows] for k, f in predfun.items()}
    b, B, P = fit(rows, yfun, preds, nboot)
    se = B.std(0, ddof=1)
    z3 = (b - P['3pt']) / se
    out = '%-26s n=%3d fit ' % (name, len(rows)) + ' '.join('%+.4f' % x for x in b[:4]) + ' se ' + ' '.join('%.4f' % x for x in se[:4])
    out += ' | z3pt ' + ' '.join('%+.1f' % x for x in z3[:4])
    if 'clo' in P:
        zc = (b - P['clo']) / se
        out += ' | zclo ' + ' '.join('%+.1f' % x for x in zc[:4])
        g = amp(b, P['clo'], P['3pt']); gs = np.std([amp(x, P['clo'], P['3pt']) for x in B], ddof=1)
        out += ' | gC %.2f+-%.2f' % (g, gs)
    if 'lead' in P:
        g = amp(b, P['lead'], P['3pt']); gs = np.std([amp(x, P['lead'], P['3pt']) for x in B], ddof=1)
        out += ' | g %.3f+-%.3f' % (g, gs)
    print(out, flush=True)
    return b, se, P


def y_of(r, lo=0, hi=1e9, trim=None):
    m = (r['g'] >= lo) & (r['g'] <= hi)
    y = r['la'][m] - np.log(r['Lj'][m])
    if trim is not None:
        y = y[y > trim]
    return np.var(y, ddof=1)


def cl_pred(r, lo=600, hi=1000):
    m = (r['g'] >= lo) & (r['g'] <= hi)
    Le = 1 / np.mean(1 / r['Lj'][m])
    return sum(clo(p, int(round(Le * 1000)))[0] for p in factor(r['q']))


def cp_pred(r, lo=600, hi=1000):
    m = (r['g'] >= lo) & (r['g'] <= hi)
    Le = 1 / np.mean(1 / r['Lj'][m])
    return sum(compact(p, Le) for p in factor(r['q']))


if __name__ == '__main__':
    rows = [r for r in load() if r['D'] not in (197, -395)]
    base = {'3pt': lambda r: r['pred']['cum2_3pt'], 'clo': cl_pred, 'lead': lambda r: r['pred']['cum2_lead']}
    report('all (frozen 3pt)', rows, y_of, base)
    report('all (compact formula)', rows, y_of, {'3pt': cp_pred, 'clo': cl_pred, 'lead': lambda r: r['pred']['cum2_lead']})
    for s in ['C', 'Bp']:
        report('set ' + s, [r for r in rows if r['set'] == s], y_of, base)
    report('D > 0', [r for r in rows if r['D'] > 0], y_of, base)
    report('D < 0', [r for r in rows if r['D'] < 0], y_of, base)
    for s, cut in [('C', 250), ('Bp', 550)]:
        report('set %s |D| <= %d' % (s, cut), [r for r in rows if r['set'] == s and r['q'] <= cut], y_of, base)
        report('set %s |D| > %d' % (s, cut), [r for r in rows if r['set'] == s and r['q'] > cut], y_of, base)
    ally = np.concatenate([r['la'] - np.log(r['Lj']) for r in rows])
    for qq in [0.001, 0.01]:
        t = np.quantile(ally, qq)
        print('trim pooled lowest %.1f%% (y > %.3f):' % (100 * qq, t), end=' ')
        report('', rows, lambda r, t=t: y_of(r, trim=t), base)
    # height bands with the compact-formula prediction
    for lo, hi in [(600, 800), (800, 1000)]:
        report('band %d-%d (compact)' % (lo, hi), rows, lambda r, lo=lo, hi=hi: y_of(r, lo, hi),
               {'3pt': lambda r, lo=lo, hi=hi: cp_pred(r, lo, hi), 'clo': lambda r, lo=lo, hi=hi: cl_pred(r, lo, hi),
                'lead': lambda r: r['pred']['cum2_lead']})
