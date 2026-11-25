"""Round 142's stored values (heights 100-600): the same reading with the retro predictions and the wave closure,
per set (R120, Dn, Fn) and pooled; remove-the-wave map. Own code."""
import os, numpy as np, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from myfit import design, ols, PR
from altlib import clo, factor
from sympy import jacobi_symbol
R = os.path.dirname(os.path.abspath(__file__))+'/../../../../'
PRED = json.load(open(R + 'writeup/code/study/lowzero/varlp/retro_pred_r142.json'))


def kron(D, p):
    if p == 2:
        return 0 if D % 2 == 0 else (1 if D % 8 in (1, 7) else -1)
    return int(jacobi_symbol(D % p, p)) if D % p else 0


rows = []
for s in ['R120', 'Dn', 'Fn']:
    V = np.load(R + 'writeup/data/study/r142/vals_%s.npz' % s)
    for e in PRED[s]:
        D = e['D']
        if D == 1153: continue
        g = V['D%d__gamma' % D]; la = V['D%d__logabs' % D]
        m = (g >= 100) & (g <= 600)
        Lj = np.log(abs(D) * g[m] / (2 * np.pi))
        rows.append(dict(set=s, D=D, q=abs(D), g=g[m], la=la[m], Lj=Lj, pred=e))
print('L-functions', len(rows))


def run(rr, name, nboot=1000):
    v = np.array([np.var(r['la'] - np.log(r['Lj']), ddof=1) for r in rr])
    X, nn = design(rr)
    b = ols(X, v)[nn:]
    rng = np.random.default_rng(20261126); B = []
    for _ in range(nboot):
        ix = rng.integers(0, len(rr), len(rr)); Xb = X[ix].copy(); Xb[:, 1:] -= Xb[:, 1:].mean(0)
        B.append(ols(Xb, v[ix])[nn:])
    B = np.array(B); se = B.std(0, ddof=1)
    p3 = ols(X, np.array([r['pred']['cum2_3pt'] for r in rr]))[nn:]
    pl = ols(X, np.array([r['pred']['cum2_lead'] for r in rr]))[nn:]
    pc = ols(X, np.array([sum(clo(p, int(round(1000 / np.mean(1 / r['Lj']))))[0] for p in factor(r['q'])) for r in rr]))[nn:]
    Le = 1 / np.mean(np.concatenate([1 / r['Lj'] for r in rr]))
    amp = lambda bb, z0: np.sum((bb[:3] - z0[:3]) * (p3[:3] - z0[:3])) / np.sum((p3[:3] - z0[:3]) ** 2)
    print('%-6s n=%d <L>=%.2f fit %s se %s' % (name, len(rr), Le, ' '.join('%+.4f' % x for x in b[:4]), ' '.join('%.4f' % x for x in se[:4])))
    print('        3pt %s z %s | closure %s z %s | fit-3pt %s | 3pt-clo %s' % (
        ' '.join('%+.4f' % x for x in p3[:4]), ' '.join('%+.1f' % x for x in ((b - p3) / se)[:4]),
        ' '.join('%+.4f' % x for x in pc[:4]), ' '.join('%+.1f' % x for x in ((b - pc) / se)[:4]),
        ' '.join('%+.4f' % x for x in (b - p3)[:4]), ' '.join('%+.4f' % x for x in (p3 - pc)[:4])))
    print('        gamma(lead->3pt) %.3f +- %.3f   gC(closure->3pt) %.2f +- %.2f' % (
        amp(b, pl), np.std([amp(x, pl) for x in B], ddof=1), amp(b, pc), np.std([amp(x, pc) for x in B], ddof=1)))
    return b, B, p3, pc


run(rows, 'all')
for s in ['R120', 'Dn', 'Fn']:
    run([r for r in rows if r['set'] == s], s)
print('remove-the-wave map, round 142 (all sets):')
for p in [2, 3, 5, 7]:
    out = []
    for r in rows:
        ch = kron(r['D'], p)
        if ch == 0: continue
        y = r['la'] - np.log(r['Lj']); th = r['g'] * np.log(p); L = r['Lj']
        k = np.arange(1, 41)[:, None]
        c = (ch ** k) * p ** (-k / 2.0) * np.cos(k * th)
        A = (c / k).sum(0); eps = -(2 * np.log(p) / L) * c.sum(0)
        w = 1 / (1 + eps); w = w / w.sum(); Y2 = y - A - np.log1p(eps)
        m2 = np.sum(w * Y2); v2 = np.sum(w * (Y2 - m2) ** 2) * len(y) / (len(y) - 1)
        out.append((v2 - np.var(y, ddof=1), np.mean(L)))
    out = np.array(out)
    print('  p=%d n=%d dVar %+.5f +- %.5f  <L> %.2f  closure %+.5f  departure %+.5f' % (
        p, len(out), out[:, 0].mean(), out[:, 0].std(ddof=1) / np.sqrt(len(out)), out[:, 1].mean(),
        clo(p, int(round(out[:, 1].mean() * 1000)))[0], out[:, 0].mean() - clo(p, int(round(out[:, 1].mean() * 1000)))[0]))
