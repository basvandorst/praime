"""Diagnostics: (1) why the regression-matched off-diagonal part differs from the group mean at the mean L
(decomposition: L-averaging vs other bad primes vs base L-dependence); (2) L-dependence of the p = 2 off part;
(3) noise levels from round 123 for the power of the fresh tests. usage: diag.py"""
import os, sys, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bklib import Grid, dS_diag, dsig_off, sigma2_full, LGRID, W, L_nodes

HERE = os.path.dirname(os.path.abspath(__file__)); LZ = os.path.abspath(os.path.join(HERE, '..'))
PR = (2, 3, 5, 7, 11, 13)
idx = lambda l: int(round(l / W)) - 1
def pf(n):
    n = abs(n); ps = []; p = 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: ps.append(n)
    return ps
G = Grid(eps_max=7.6, h=0.0025, pmax=10 ** 5)   # 1e5 primes suffice (checked: Phi_off differs by < 3e-6)

print('(2) p = 2 off part at l = 7 and 7.5 as a function of fixed L')
for L in (7.5, 8.0, 8.5, 9.0, 9.38, 10.0, 10.5, 11.0):
    o = dsig_off(G, L, 2)
    print('  L %.2f: off(7) %+.4f off(7.5) %+.4f   diag(7) %+.4f' % (L, o[idx(7)], o[idx(7.5)], dS_diag(2, L, 7)))

rows = json.load(open(os.path.join(LZ, 'r120.json')))['rows']
two = [r for r in rows if r['D'] % 2 == 0]
print('\n(1) decomposition for p = 2 (n = %d with 2 | D) at l = 7 and 7.5' % len(two))
for l in (7.0, 7.5):
    a = np.mean([dsig_off(G, r['L'], 2)[idx(l)] for r in two])                       # single prime, mean L
    b = np.mean([sum(w * dsig_off(G, L, 2)[idx(l)] for L, w in zip(*L_nodes(abs(r['D']), 100, 600, 8))) for r in two])  # single, L-avg
    c = np.mean([dsig_off(G, r['L'], 2, rest=tuple(q for q in pf(r['D']) if q != 2))[idx(l)] for r in two])               # other bad primes, mean L
    d = np.mean([sum(w * dsig_off(G, L, 2, rest=tuple(q for q in pf(r['D']) if q != 2))[idx(l)] for L, w in zip(*L_nodes(abs(r['D']), 100, 600, 8))) for r in two])
    print('  l %.1f: single/meanL %+.4f | single/L-avg %+.4f | with other bad primes/meanL %+.4f | both %+.4f' % (l, a, b, c, d))

# (3) noise: s(l) of round 123 regressions, and the expected se for the fresh designs
r123 = json.load(open(os.path.join(LZ, 'r123.json')))
X = np.array([[1] + [int(r['D'] % p == 0) for p in PR] for r in rows], float); XtXi = np.linalg.inv(X.T @ X)
LS = [0.5 * i for i in range(1, 17)]
s = np.array([np.mean([r123['se']['%d_%s' % (p, l)] / math.sqrt(XtXi[i + 1, i + 1]) for i, p in enumerate(PR)]) for l in LS])
print('\n(3) residual SD s(l) of round 123 per-l regressions:', ' '.join('%.4f' % v for v in s))
pj = json.load(open(os.path.join(HERE, 'predict.json')))
for T in ('B', 'C'):
    Dl = pj['fresh'][T].get('D') or [r['D'] for r in rows]
    Xn = np.array([[1] + [int(D % p == 0) for p in PR] for D in Dl], float); Xi = np.linalg.inv(Xn.T @ Xn)
    full = np.array(pj['fresh'][T]['full']); dg = np.array(pj['fresh'][T]['diag'])
    off = np.array([full[idx(l)] - dg[idx(l)] for l in LS])
    se = np.array([[s[j] * math.sqrt(Xi[i + 1, i + 1]) for i in range(6)] for j in range(16)])
    print('  Test %s: expected se (from s of round 123) p=2: %s' % (T, ' '.join('%.4f' % v for v in se[:, 0])))
    print('          expected chi2 excess of the diagonal-only law, sum (off/se)^2 over 96 cells: %.0f' % ((off / se) ** 2).sum())
    print('          expected z of gamma (slope of resid on off, se-naive, independent cells): %.1f;  /sqrt(96/15) for correlation: %.1f' % (
        math.sqrt(((off / se) ** 2).sum()), math.sqrt(((off / se) ** 2).sum() / (96 / 15))))
# round 123 regression-matched: same numbers, as a calibration of the above
# complex: r121 per-character residual SD at l = 1 against s_real(1)
r121 = json.load(open(os.path.join(LZ, 'r121.json')))['rows']
y = np.array([r['sigma2_1'] for r in r121]); x = np.array([r['x'] for r in r121]); A = np.column_stack([np.ones(len(y)), x])
res = y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
print('\n  complex (r121) residual SD at l = 1: %.4f; real (r123) s(1) = %.4f' % (res.std(ddof=2), s[1]))
NsA = [r['N'] for r in r121]
XA = np.array([[1] + [int(N % p == 0) for p in (2, 3, 5, 7)] for N in NsA], float); XiA = np.linalg.inv(XA.T @ XA)
fullA = np.array(pj['fresh']['A']['full']); dgA = np.array(pj['fresh']['A']['diag'])
offA = np.array([fullA[idx(l)] - dgA[idx(l)] for l in LS])
sA = s * res.std(ddof=2) / s[1]
seA = np.array([[sA[j] * math.sqrt(XiA[i + 1, i + 1]) for i in range(4)] for j in range(16)])
print('  Test A: expected se p=2: %s' % ' '.join('%.4f' % v for v in seA[:, 0]))
print('          expected chi2 excess of the diagonal-only law over 64 cells: %.0f; z(gamma) naive %.1f, corr-adjusted %.1f' % (
    ((offA / seA) ** 2).sum(), math.sqrt(((offA / seA) ** 2).sum()), math.sqrt(((offA / seA) ** 2).sum() / (64 / 10))))
