"""Alternatives to the three-point prediction, matched through the registered design on the fresh data.
 closure : wave closure (explicit-formula wave of each bad prime + independence), all orders, per L-function
           sum over its bad primes at its zero-weighted L.
 first   : first order -Li2/2 + (ln p/L) l(2+l)   (= LDM)
 M-only  : 3pt with Delta V3,P replaced by the far-field model M (uses derivation's L^2*residual table: NOT own)
 noT3    : 3pt with Delta V3 = 0 (GUE T3 or T3 without bad-prime factors) (uses derivation's tables: NOT own)
Also the mean shift under the closure, and the 'remove the wave' map on the real values."""
import os, numpy as np, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from myfit import load, stat, design, ols, PR
from mpmath import polylog
from functools import lru_cache


def factor(n):
    f = []; d = 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: f.append(n)
    return f


@lru_cache(maxsize=None)
def clo(p, Lr, n=1 << 13):
    L = Lr / 1000.
    K = max(4, int(60 / np.log2(p)) + 2)
    th = 2 * np.pi * np.arange(n) / n
    k = np.arange(1, K + 1)[:, None]
    c = p ** (-k / 2.0) * np.cos(k * th)
    A = (c / k).sum(0)
    eps = -(2 * np.log(p) / L) * c.sum(0)
    w = 1 + eps
    X = A + np.log(w)
    m = np.mean(w * X)
    return -(np.mean(w * X * X) - m * m), -m


def first(p, L):
    l = np.log(p / (p - 1.0))
    return -0.5 * float(polylog(2, 1.0 / p)) + np.log(p) / L * l * (2 + l)


if __name__ == "__main__":
  pass
rows = load()
kept = [r for r in rows if r['D'] not in (197, -395)]
v = np.array([stat(r, 'lnL')[0] for r in kept]); mu = np.array([stat(r, 'lnL')[1] for r in kept])
X, nn = design(kept)
b = ols(X, v); bm = ols(X, mu)
rng = np.random.default_rng(20261126); n = len(kept); B = []; Bm = []
for _ in range(2000):
    ix = rng.integers(0, n, n); Xb = X[ix].copy(); Xb[:, 1:] -= Xb[:, 1:].mean(0)
    B.append(ols(Xb, v[ix])); Bm.append(ols(Xb, mu[ix]))
B = np.array(B); Bm = np.array(Bm); se = B.std(0, ddof=1); sem = Bm.std(0, ddof=1)
Cv = np.cov(B[:, nn:nn + 4].T)

pc = []; pcm = []; pf = []
for r in kept:
    Le = 1.0 / r['pred']['mean_invL']
    fs = factor(r['q'])
    pc.append(sum(clo(p, int(round(Le * 1000)))[0] for p in fs))
    pcm.append(sum(clo(p, int(round(Le * 1000)))[1] for p in fs))
    pf.append(sum(first(p, Le) for p in fs))
cc = ols(X, np.array(pc))[nn:]; ccm = ols(X, np.array(pcm))[nn:]; cf = ols(X, np.array(pf))[nn:]
c3 = ols(X, np.array([r['pred']['cum2_3pt'] for r in kept]))[nn:]
c3m = ols(X, np.array([r['pred']['cum1_3pt'] for r in kept]))[nn:]
cld = ols(X, np.array([r['pred']['cum2_ldm'] for r in kept]))[nn:]
Leff = 1 / np.mean([r['pred']['mean_invL'] for r in kept])
print('zero-weighted L of kept: %.3f' % Leff)
# derivation's table: L^2 * (dV3P - intM) and dV3P at L = 10, 11 (p = 2,3,5,7) -> interpolate to Leff
tab = {2: ((10, .1723, -5.45), (11, .1533, -5.44)), 3: ((10, .1770, -6.33), (11, .1587, -6.35)),
       5: ((10, .1499, -6.49), (11, .1359, -6.49)), 7: ((10, .1262, -6.22), (11, .1152, -6.21))}
fb = b[nn:]
print('p                 ' + ' '.join('%9d' % p for p in PR))
print('fit               ' + ' '.join('%+9.5f' % x for x in fb))
print('se                ' + ' '.join('%9.5f' % x for x in se[nn:]))
for name, c in [('3pt (frozen)', c3), ('LDM (frozen)', cld), ('first order', cf), ('wave closure', cc)]:
    z = (fb - c) / se[nn:]
    d = (fb - c)[:4]
    print('%-17s' % name + ' '.join('%+9.5f' % x for x in c) + '   z ' + ' '.join('%+.1f' % x for x in z)
          + '   chi2(4) %.1f' % (d @ np.linalg.solve(Cv, d)))
alts = {}
for nm in ['M-only', 'noT3']:
    c = c3[:4].copy()
    for i, p in enumerate([2, 3, 5, 7]):
        (L0, P0, r0), (L1, P1, r1) = tab[p]
        t = (Leff - L0) / (L1 - L0)
        P = P0 + t * (P1 - P0); res = (r0 + t * (r1 - r0)) / Leff ** 2
        Q = np.log(p) / Leff * np.log(p / (p - 1.)) ** 2
        c[i] = c3[i] - res if nm == 'M-only' else c3[i] - P - Q
    z = (fb[:4] - c) / se[nn:nn + 4]; d = fb[:4] - c
    print('%-17s' % nm + ' '.join('%+9.5f' % x for x in c) + '   z ' + ' '.join('%+.1f' % x for x in z)
          + '   chi2(4) %.0f' % (d @ np.linalg.solve(Cv, d)))


def amp(bb, z0, one):
    return np.sum((bb - z0) * (one - z0)) / np.sum((one - z0) ** 2)
for nm, z0 in [('closure', cc), ('first', cf)]:
    g = amp(fb[:3], z0[:3], c3[:3]); gs = np.std([amp(x[nn:nn + 3], z0[:3], c3[:3]) for x in B], ddof=1)
    g4 = amp(fb[:4], z0[:4], c3[:4]); gs4 = np.std([amp(x[nn:nn + 4], z0[:4], c3[:4]) for x in B], ddof=1)
    print('amplitude (0 = %s, 1 = 3pt): p=2,3,5: %.3f +- %.3f ; p=2..7: %.3f +- %.3f' % (nm, g, gs, g4, gs4))
print('\nmean shift: fit ' + ' '.join('%+.5f' % x for x in bm[nn:]))
print('            se  ' + ' '.join('%.5f' % x for x in sem[nn:]))
print('       3pt/2pt  ' + ' '.join('%+.5f' % x for x in c3m) + '  z ' + ' '.join('%+.1f' % x for x in (bm[nn:] - c3m) / sem[nn:]))
print('       closure  ' + ' '.join('%+.5f' % x for x in ccm) + '  z ' + ' '.join('%+.1f' % x for x in (bm[nn:] - ccm) / sem[nn:]))
fo = np.array([np.log(p) * np.log(p / (p - 1.)) / Leff for p in PR])
print('  first order   ' + ' '.join('%+.5f' % x for x in fo))

# remove-the-wave map on the real values, p good (chi(p) = +-1): predicted contrast if p were missing
print('\nremove-the-wave map on real data (mean over L-functions with p good, unweighted):')
from sympy import jacobi_symbol
def kron(D, p):
    if p == 2:
        if D % 2 == 0: return 0
        return 1 if D % 8 in (1, 7) else -1
    return int(jacobi_symbol(D % p, p)) if D % p else 0
for p in [2, 3, 5, 7]:
    out = []
    for r in kept:
        ch = kron(r['D'], p)
        if ch == 0: continue
        y = r['la'] - np.log(r['Lj']); th = r['g'] * np.log(p); L = r['Lj']
        K = 40; k = np.arange(1, K + 1)[:, None]
        c = (ch ** k) * p ** (-k / 2.0) * np.cos(k * th)
        A = (c / k).sum(0); eps = -(2 * np.log(p) / L) * c.sum(0)
        w = 1 / (1 + eps); w = w / w.sum()
        Y2 = y - A - np.log1p(eps)
        m2 = np.sum(w * Y2); v2 = np.sum(w * (Y2 - m2) ** 2) * len(y) / (len(y) - 1)
        out.append((v2 - np.var(y, ddof=1), m2 - y.mean(), np.mean(L)))
    out = np.array(out)
    print('  p=%d n=%d  dVar %+.5f +- %.5f   dmean %+.5f   <L> %.2f   (closure at <L>: %+.5f)' % (
        p, len(out), out[:, 0].mean(), out[:, 0].std(ddof=1) / np.sqrt(len(out)), out[:, 1].mean(), out[:, 2].mean(),
        clo(p, int(round(out[:, 2].mean() * 1000)))[0]))
