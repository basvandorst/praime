"""Replication analysis: residual Re/Im(G - pred) on the undrawn genus-2 curves (rep/*.npz), by root number, profile in T,
window-profile correlation with rounds 151+152, the 2 | N subgroup; then everything pooled (rounds 151 + 152 + replication)."""
import json, math, sys, glob
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, load, jk_mean
rep = []
for fn in sorted(glob.glob('rep/*.npz')):
    z = np.load(fn); info = json.loads(str(z['info']))
    if info.get('status') != 'COMPLETE': print('NOT COMPLETE', fn, info.get('status')); continue
    rep.append(dict(key=info['key'], rnd=0, cls='g2', Q=info['Q'], mu=info['mu'], eps=info['eps'][0], c=info['c'], T2=150, g=z['g'], info=info))
old = [r for r in load() if r['cls'] == 'g2']
Ts = [10, 20, 30, 60, 90, 120, 150]
def resid(r, a, b):
    g = r['g']; m = (g > a) & (g <= b)
    return np.exp(2j * theta(g[m], r['Q'], r['mu'])).sum() + r['eps'] * r['c'] * (b - a) / (2 * math.pi)
for lab, S in (('replication', rep), ('rounds 151+152', old), ('all pooled', old + rep)):
    x = np.array([resid(r, 0, 150) for r in S]); ep = np.array([r['eps'] for r in S])
    m, s = jk_mean(x.real); mi, si = jk_mean(x.imag)
    mp, sp = jk_mean(x.real[ep > 0]); mm, sm = jk_mean(x.real[ep < 0])
    print('%-15s n=%2d  Re %+5.2f +- %4.2f  Im %+5.2f +- %4.2f  | eps+1 n=%2d %+5.2f +- %4.2f  eps-1 n=%2d %+5.2f +- %4.2f | mean c eps+1 %.2f eps-1 %.2f | sd %.1f' % (
        lab, len(S), m, s, mi, si, (ep > 0).sum(), mp, sp, (ep < 0).sum(), mm, sm, np.mean([r['c'] for r in S if r['eps'] > 0]), np.mean([r['c'] for r in S if r['eps'] < 0]), x.real.std(ddof=1)))
    print('      profile ' + ' '.join('T%d %+.2f(%.2f)' % (T, *jk_mean([resid(r, 0, T).real for r in S])) for T in Ts))
    # (20,150] increment
    print('      increment (20,150]: %+.2f +- %.2f' % jk_mean([resid(r, 20, 150).real for r in S]))
    # slopes (registered estimator)
    X = np.array([((resid(r, 0, 150) - r['eps'] * r['c'] * 150 / (2 * math.pi) * 0).real, 0) for r in S])
    G = np.array([resid(r, 0, 150) - r['eps'] * r['c'] * 150 / (2 * math.pi) for r in S]); P = np.array([-r['eps'] * r['c'] * 150 / (2 * math.pi) for r in S])
    for e, nm in ((None, 'all'), (1, 'eps+1'), (-1, 'eps-1')):
        sel = np.ones(len(S), bool) if e is None else (ep == e)
        a = (G[sel].real * P[sel]); b = P[sel] ** 2; sx, sy = a.sum(), b.sum(); loo = (sx - a) / (sy - b); n = sel.sum()
        print('      slope %-6s %.3f +- %.3f' % (nm, sx / sy, math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))), end='')
    print()
edges = np.arange(10, 151, 10)
def W(S):
    return np.array([[resid(r, edges[j], edges[j + 1]) for j in range(len(edges) - 1)] for r in S])
if len(rep) >= 10:
    Wr = W(rep).mean(0); Wo = W(old).mean(0)
    print('window profile correlation replication vs rounds 151+152: Re %.3f  Im %.3f' % (np.corrcoef(Wr.real, Wo.real)[0, 1], np.corrcoef(Wr.imag, Wo.imag)[0, 1]))
    # std of a complex array is the complex spread (a real number): take Re and Im separately (fixed 25 Nov 2026)
    Wx = W(rep); n = len(rep); mm = Wx.mean(0); sr = Wx.real.std(0, ddof=1) / math.sqrt(n); si = Wx.imag.std(0, ddof=1) / math.sqrt(n)
    print('replication window chi2 vs 0: Re %.1f  Im %.1f (14 windows, 14 dof)' % (np.sum((mm.real / sr) ** 2), np.sum((mm.imag / si) ** 2)))
two = [r for r in rep if int(round(r['Q'])) % 2 == 0]
if two:
    x = np.array([resid(r, 0, 150).real for r in two]); print('replication 2 | N: n=%d mean %+.2f (values %s)' % (len(two), x.mean(), np.round(x, 1).tolist()))
