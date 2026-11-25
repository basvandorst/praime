"""H5 (Sato-Tate moments per curve) + robustness of the slope (bootstrap, robust estimators, null noise from far
twists r < 1) + covariates + comparison with the other degree-4 classes. Own code; a_p from own counts (cown, p <= 1e5)
for the 51 sampled curves and from lfunan (p <= 3e4) for the other 49."""
import sys, os, math, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from g2data import load, slope
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
rng = np.random.default_rng(20261125)
allrows = load(('g2', 'sym3', 'rs')); G2 = [r for r in allrows if r['cls'] == 'g2']
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
PR3 = primes_upto(30000)
for r in G2:
    fn = 'cown/%s.npz' % r['key']
    if os.path.exists(fn):
        z = np.load(fn); p = z['p']; N = int(z['N']); good = np.array([N % int(q) != 0 for q in p]); e1 = z['e1'][good & (p > 100)]
        e2 = z['e2'][good & (p > 30) & (p <= 1000)]; src = 'own<=1e5'
    else:
        L = pari('lfungenus2([%s,%s])' % (r['spec'][1], r['spec'][2])); N = int(pari.lfunparams(L)[0]); an = pari.lfunan(L, 30000)
        p = PR3[(PR3 > 100) & np.array([N % int(q) != 0 for q in PR3])]; e1 = np.array([float(an[int(q) - 1]) for q in p]) / np.sqrt(p)
        e2 = np.array([]); src = 'lfunan<=3e4'
    n = len(e1); r['N'] = N; r['src'] = src
    r['M2'], r['M4'], r['M6'] = [float(np.mean(e1 ** k)) for k in (2, 4, 6)]
    r['zM2'] = (r['M2'] - 1) / math.sqrt(2 / n); r['zM4'] = (r['M4'] - 3) / math.sqrt((84 - 9) / n); r['zM6'] = (r['M6'] - 14) / math.sqrt((594 - 196) / n)
    r['E2'] = float(np.mean(e2)) if len(e2) else float('nan'); r['E22'] = float(np.mean(e2 ** 2)) if len(e2) else float('nan')
    r['pred'] = -r['eps'] * r['c'] * 150 / (2 * math.pi); r['ratio'] = r['G'].real / r['pred'] if r['pred'] != 0 else float('nan')
print('H5: Sato-Tate moments of a_p/sqrt(p) (USp(4): M2 = 1, M4 = 3, M6 = 14; E[e2] = 1, E[e2^2] = 2)')
zz = np.array([[r['zM2'], r['zM4'], r['zM6']] for r in G2])
print('  z-scores over 100 curves: M2 mean %+.2f sd %.2f max|z| %.2f ; M4 mean %+.2f sd %.2f max|z| %.2f ; M6 mean %+.2f sd %.2f max|z| %.2f' % (
    zz[:, 0].mean(), zz[:, 0].std(), np.abs(zz[:, 0]).max(), zz[:, 1].mean(), zz[:, 1].std(), np.abs(zz[:, 1]).max(), zz[:, 2].mean(), zz[:, 2].std(), np.abs(zz[:, 2]).max()))
E2 = np.array([r['E2'] for r in G2 if not math.isnan(r['E2'])]); E22 = np.array([r['E22'] for r in G2 if not math.isnan(r['E22'])])
print('  e2 moments (30 < p <= 1000, 51 curves): E[e2] %.3f +- %.3f (range %.3f..%.3f), E[e2^2] %.3f (range %.3f..%.3f)' % (E2.mean(), E2.std() / math.sqrt(len(E2)), E2.min(), E2.max(), E22.mean(), E22.min(), E22.max()))
low = sorted([r for r in G2 if r['eps'] > 0], key=lambda r: r['ratio'])[:8]
print('  the 8 lowest G/pred among eps = +1:')
for r in low: print('    %-13s c %.3f ratio %+.3f  M2 %.3f (z %+.1f) M4 %.3f (z %+.1f) M6 %.2f (z %+.1f) E[e2] %s  [%s]' % (r['key'], r['c'], r['ratio'], r['M2'], r['zM2'], r['M4'], r['zM4'], r['M6'], r['zM6'], ('%.3f' % r['E2']) if not math.isnan(r['E2']) else '-', r['src']))
pos = [r for r in G2 if r['eps'] > 0]
w = np.array([r['pred'] ** 2 for r in pos]); rat = np.array([r['ratio'] for r in pos])
for nm in ('zM4', 'zM2', 'zM6'):
    x = np.array([r[nm] for r in pos]); X = np.vstack([np.ones_like(x), x]).T; W = np.diag(w)
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ rat)
    # jackknife se of the weighted slope on the covariate
    bl = []
    for i in range(len(x)):
        m = np.arange(len(x)) != i; Xm = X[m]; Wm = np.diag(w[m]); bl.append(np.linalg.solve(Xm.T @ Wm @ Xm, Xm.T @ Wm @ rat[m])[1])
    bl = np.array(bl); se = math.sqrt((len(x) - 1) / len(x) * np.sum((bl - bl.mean()) ** 2))
    print('  eps+1: weighted regression of G/pred on %s: slope %+.4f +- %.4f (intercept %.3f)' % (nm, beta[1], se, beta[0]))
# robustness of the slope
print('\nrobustness (all 100, eps+1):')
for nm, S in (('all 100', G2), ('eps+1', pos)):
    Gs = np.array([r['G'] for r in S]); P = np.array([r['pred'] for r in S], complex)
    b0, se0 = slope(Gs, P)
    bs = []
    for _ in range(4000):
        i = rng.integers(0, len(S), len(S)); x = (Gs[i] * P[i].conj()).real.sum(); y = (np.abs(P[i]) ** 2).sum(); bs.append(x / y)
    bs = np.array(bs)
    ratios = np.array([(g * p.conjugate()).real / abs(p) ** 2 for g, p in zip(Gs, P)]); ww = np.abs(P) ** 2
    order = np.argsort(ratios); cw = np.cumsum(ww[order]) / ww.sum(); wmed = ratios[order][np.searchsorted(cw, 0.5)]
    # drop the 5 heaviest curves
    hv = np.argsort(-ww)[:5]; m = np.ones(len(S), bool); m[hv] = False; bdrop, sedrop = slope(Gs[m], P[m])
    # Huber-type: clip standardized residuals at 2 sd of the per-curve scatter (sd 11)
    res = (Gs - P).real; resc = np.clip(res, -22, 22); bh = 1 + (resc * P.real).sum() / (P.real ** 2).sum()
    print('  %-8s slope %.4f jk %.4f | bootstrap sd %.4f, 2.5-97.5%% [%.3f, %.3f], P(b >= 1) = %.4f | weighted median ratio %.3f | without 5 heaviest %.3f +- %.3f | clipped residuals %.3f' % (
        nm, b0, se0, bs.std(), *np.percentile(bs, [2.5, 97.5]), np.mean(bs >= 1), wmed, bdrop, sedrop, bh))
# null noise: twists r < 1 far from small rationals, u = -ln r in [6, 9]
us = rng.uniform(6, 9, 200)
sd_null = []; sd_obs = []
for r in G2:
    Z = np.array([(r['ph'] * np.exp(1j * u * r['gg'])).sum() for u in us])
    r['nullvar'] = float(np.mean(Z.real ** 2))
P = np.array([r['pred'] for r in G2]); V = np.array([r['nullvar'] for r in G2])
se_null = math.sqrt((P ** 2 * V).sum()) / (P ** 2).sum()
res = np.array([(r['G'] - r['pred']).real for r in G2])
print('\nnull noise from r = e^{-u}, u in [6, 9] (200 values): mean Var(Re) per curve %.1f (sd %.1f); observed Re residual variance at r = 1: %.1f' % (V.mean(), math.sqrt(V.mean()), res.var()))
print('  slope error implied by the null noise: %.4f (jackknife 0.0198); deficit / null se = %.2f' % (se_null, (1 - 0.9286) / se_null))
# other degree-4 classes pooled
S34 = [r for r in allrows if r['cls'] in ('sym3', 'rs')]
b34, s34 = slope([r['G'] for r in S34], [-r['eps'] * r['c'] * r['T2'] / (2 * math.pi) for r in S34])
bs3, ss3 = slope([r['G'] for r in S34 if r['cls'] == 'sym3' and r['eps'] > 0], [-r['eps'] * r['c'] * r['T2'] / (2 * math.pi) for r in S34 if r['cls'] == 'sym3' and r['eps'] > 0])
print('\nSym^3 + RS pooled (n=%d): %.3f +- %.3f ; Sym^3 eps+1 (n=%d): %.3f +- %.3f ; genus 2 minus (Sym^3+RS): %.3f +- %.3f (%.1f se)' % (
    len(S34), b34, s34, sum(1 for r in S34 if r['cls'] == 'sym3' and r['eps'] > 0), bs3, ss3, 0.9286 - b34, math.hypot(0.0198, s34), (0.9286 - b34) / math.hypot(0.0198, s34)))
b4, s4 = slope([r['G'] for r in allrows], [-r['eps'] * r['c'] * r['T2'] / (2 * math.pi) for r in allrows])
print('all degree 4 (n=%d): %.3f +- %.3f' % (len(allrows), b4, s4))
# covariates on eps+1: 2 | N, rank parity (m0), number of bad primes, size of the p = 2 factor
def sub(cond, nm):
    A = [r for r in pos if cond(r)]; B = [r for r in pos if not cond(r)]
    ba = slope([r['G'] for r in A], [r['pred'] for r in A]) if len(A) > 2 else (float('nan'), float('nan'))
    bb = slope([r['G'] for r in B], [r['pred'] for r in B]) if len(B) > 2 else (float('nan'), float('nan'))
    print('  eps+1 %-26s yes n=%2d %.3f +- %.3f | no n=%2d %.3f +- %.3f' % (nm, len(A), *ba, len(B), *bb))
print('\ncovariates (eps = +1 curves):')
sub(lambda r: r['N'] % 2 == 0, '2 | N')
sub(lambda r: r['m0'] == 2, 'central zero of order 2')
sub(lambda r: len([q for q in range(2, 3000) if r['N'] % q == 0 and all(q % d for d in range(2, int(q ** 0.5) + 1))]) >= 2, '>= 2 bad primes')
sub(lambda r: r['c'] > 3, 'c > 3')
sub(lambda r: r['N'] > 1500, 'N > 1500')
sub(lambda r: len(r['spec'][1]) and 'x^6' not in r['spec'][1], 'deg f = 5')
print('  m0 values eps+1:', sorted(set(r['m0'] for r in pos)), ' counts', {m: sum(1 for r in pos if r['m0'] == m) for m in set(r['m0'] for r in pos)})
json.dump({r['key']: dict(N=r['N'], eps=r['eps'], c=r['c'], ratio=r['ratio'], M2=r['M2'], M4=r['M4'], M6=r['M6'], zM4=r['zM4'], E2=r['E2'], m0=r['m0'], nullvar=r['nullvar'], src=r['src']) for r in G2}, open('h5_table.json', 'w'), indent=1)
