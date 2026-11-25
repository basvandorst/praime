"""(1) smooth-weight (V^2) diagonal models for k = 1 and k = 2: what the controls say about them;
(2) control variates built from singles instead of pairs; fitted beta; (3) pairs split by gcd."""
import json, math, os, sys
import numpy as np
from scipy.special import erfc
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from stats import jk, slope, P, R, tri, pairs, TT, PP, SS, m
from analysis import gl, laurent_W, ST
import mpmath as mp
EG = float(mp.euler)


def szeta_pow(k):
    zs = [1.0] + [(-1) ** n * ST[n] / math.factorial(n) for n in range(k)]
    zk = [1.0] + [0.0] * k
    for _ in range(k):
        new = [0.0] * (k + 1)
        for a in range(k + 1):
            for b in range(k + 1 - a):
                new[a + b] += zk[a] * zs[b]
        zk = new
    return zk


def diagk(Ds, g, t, W):
    k = len(Ds); Zk = szeta_pow(k)
    c = [sum(Zk[a] * g[n - a] for a in range(n + 1)) for n in range(k + 1)]
    ell = sum(np.log(abs(D) * t / (2 * math.pi)) for D in Ds) / 2
    Wc = [1.0] + list(W)
    E = [ell ** n / math.factorial(n) for n in range(k + 1)]
    tot = 0.0
    for a in range(k + 1):
        for b in range(k + 1 - a):
            tot = tot + c[a] * E[b] * Wc[k - a - b]
    return 2 * tot


def main():
    M1 = m['single_trap']; M2 = m['pair_trap']; M3 = m['triple_trap']
    print('=== (1) diagonal models with weight W, controls ===')
    for lab, W in [('sharp', [0, 0, 0])] + [('V^2 delta %.1f' % d, laurent_W(lambda u, d=d: (0.5 * erfc(d * u / 2)) ** 2)) for d in (1.0, 2.0, 4.0)]:
        # k = 1: g(sigma) = h(sigma/2) = c(sigma): from hm1 = [h(0), h'(0)] with h(z) = c(2z): c'(0) = hm1[1]/2
        r1 = []
        for D in R:
            h0, h1 = SS[D]['hm1']; g = [h0, h1 / 2]
            Id = gl(lambda t: diagk((D,), g, t, W), 100, 1000); r1.append(M1[R.index(D)] / Id)
        r2 = []
        for i, p in enumerate(pairs):
            Id = gl(lambda t: diagk(p, PP[p]['g'], t, W), 100, 1000); r2.append(M2[i] / Id)
        r1 = np.array(r1); r2 = np.array(r2)
        print('  %-14s singles data/H_W %.4f +- %.4f | pairs data/H_W %.4f +- %.4f' % (
            lab, r1.mean(), r1.std(ddof=1) / math.sqrt(56), *jk(lambda k: r2[k].mean(), pairs)))
    print('=== (2) control variates ===')
    pf = np.array([TT[t]['full'] for t in tri]); r3 = M3 / pf
    rs = {D: M1[i] / SS[D]['P1'] for i, D in enumerate(R)}
    rp = {p: M2[i] / PP[p]['full'] for i, p in enumerate(pairs)}
    cs = np.array([rs[a] + rs[b] + rs[c] - 3 for a, b, c in tri])
    cpair = np.array([rp[(a, b)] + rp[(a, c)] + rp[(b, c)] - 3 for a, b, c in tri])
    print('  raw                         %.4f +- %.4f' % jk(lambda k: r3[k].mean(), tri))
    print('  pair CV (registered)        %.4f +- %.4f' % jk(lambda k: (r3[k] - cpair[k]).mean(), tri))
    print('  singles CV r3 - sum(r1-1)   %.4f +- %.4f' % jk(lambda k: (r3[k] - cs[k]).mean(), tri))
    # pair CV after removing the pair offset estimated on all 1540 pairs: E[r3] = theta_CV + 3 (rbar2 - 1)
    r2all = np.array([rp[p] for p in pairs])
    def est(k):
        keepc = set()
        sub = [tri[i] for i in k]
        chars = set(c for t in sub for c in t)
        kp = [i for i, p in enumerate(pairs) if p[0] in chars and p[1] in chars]
        return (r3[k] - cpair[k]).mean() + 3 * (r2all[kp].mean() - 1)
    print('  pair CV + 3(mean pair - 1)  %.4f +- %.4f' % jk(est, tri))
    # 'triple relative to pairs relative to singles': r3 - sum(r2 - 1) + sum(r1 - 1)  (the genuine three-body excess)
    cpair_s = np.array([(rp[(a, b)] - rs[a] - rs[b] + 1) + (rp[(a, c)] - rs[a] - rs[c] + 1) + (rp[(b, c)] - rs[b] - rs[c] + 1) for a, b, c in tri])
    print('  r3 - sum(r2 - r1 - r1 + 1) - sum(r1 - 1) [= r3 - sum(r2-1) + sum(r1-1)]  %.4f +- %.4f' % jk(lambda k: (r3[k] - cpair[k] + cs[k]).mean(), tri))
    # fitted beta (OLS of r3 - 1 on cpair), jackknifed
    def bfit(k):
        x = cpair[k]; y = r3[k]
        b = slope(x, y); return b
    print('  fitted beta (r3 on cpair)   %.3f +- %.3f' % jk(bfit, tri))
    def cvb(k):
        b = slope(cpair[k], r3[k]); return (r3[k] - b * cpair[k]).mean()
    print('  CV with fitted beta         %.4f +- %.4f' % jk(cvb, tri))
    print('=== (3) pairs by gcd of conductors ===')
    from math import gcd
    g = np.array([gcd(abs(a), abs(b)) > 1 for a, b in pairs])
    r2 = r2all
    for lab, sel in (('coprime', ~g), ('gcd > 1', g)):
        k = np.nonzero(sel)[0]; sub = [pairs[i] for i in k]
        print('  %-8s n %4d mean %.4f +- %.4f' % (lab, len(k), *jk(lambda kk: r2[k][kk].mean(), sub)))
    # triples with 2 shared: which pairs / which primes
    nsh = np.array([(gcd(abs(t[0]), abs(t[1])) > 1) + (gcd(abs(t[0]), abs(t[2])) > 1) + (gcd(abs(t[1]), abs(t[2])) > 1) for t in tri])
    k = np.nonzero(nsh == 2)[0]
    even = np.array([sum(D % 2 == 0 for D in tri[i]) for i in k])
    for e in sorted(set(even)):
        kk = k[even == e]
        print('  n shared 2, %d even D: n %3d raw %.4f CV %.4f' % (e, len(kk), r3[kk].mean(), (r3 - cpair)[kk].mean()))


if __name__ == '__main__':
    main()
