# Reviewer's own zero-side statistics: theta_E by mpmath, one-point sums G(ln r) at many r over bands,
# X12 by the banded direct double sum for ALL pairs (sigma 0.5, 1.0), saved for the regressions.
import os, time, sys
import numpy as np
import mpmath as mp
from scipy.special import erf

SC = os.path.dirname(os.path.abspath(__file__))+'/'
R = np.load(SC + 'rv_curves.npz')
ids, N = R['ids'], R['N']
Z = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/study/r144/fresh_zeros.npz')
nc = len(ids)
t0 = time.time()
mp.mp.dps = 20
zs, ph = [], []
for k in range(nc):
    z = np.sort(Z[ids[k]])
    lnN = float(mp.log(N[k])); l2p = float(mp.log(2 * mp.pi))
    th = np.array([float(t / 2 * lnN - t * l2p + mp.im(mp.loggamma(1 + 1j * mp.mpf(t)))) for t in z])
    th = np.mod(2 * th, 2 * np.pi)
    zs.append(z); ph.append(np.exp(1j * th))
print('theta done', time.time() - t0)

# rationals r = a/b, coprime, a*b <= 120 (lines) -> one-point sums on bands
from math import gcd
RAT = sorted({(a, b) for a in range(1, 121) for b in range(1, 121) if a * b <= 120 and gcd(a, b) == 1},
             key=lambda x: (x[0] * x[1], x))
BANDS = {'b20_400': (20, 400), 'b20_200': (20, 200), 'b200_400': (200, 400), 'w20_396': (20, 396),
         'w20_208': (20, 208), 'w208_396': (208, 396)}
xi = np.array([np.log(a / b) for a, b in RAT])
CTRL = np.array([0.35, -0.35, 0.55, -0.55, 0.9, -0.9, 1.3, -1.3])
G = {bn: np.zeros((nc, len(RAT)), complex) for bn in BANDS}
GC = {bn: np.zeros((nc, len(CTRL)), complex) for bn in BANDS}
for k in range(nc):
    for bn, (lo, hi) in BANDS.items():
        m = (zs[k] > lo) & (zs[k] <= hi)
        g = zs[k][m]; e = ph[k][m]
        G[bn][k] = (e[None, :] * np.exp(-1j * np.outer(xi, g))).sum(1)
        GC[bn][k] = (e[None, :] * np.exp(-1j * np.outer(CTRL, g))).sum(1)
print('one-point done', time.time() - t0)


def Xpair(i, j, sig, lo=20.0, hi=396.0):
    zi, zj, ei, ej = zs[i], zs[j], ph[i], ph[j]
    D = 7.0 * sig
    a = np.searchsorted(zj, zi - D); b = np.searchsorted(zj, zi + D)
    cnt = b - a
    rep = np.repeat(np.arange(len(zi)), cnt)
    start = np.cumsum(cnt) - cnt
    j2 = a[rep] + (np.arange(cnt.sum()) - start[rep])
    g1 = zi[rep]; g2 = zj[j2]
    m = 0.5 * (g1 + g2)
    s = sig / np.sqrt(2)
    W = 0.5 * (erf((hi - m) / s) - erf((lo - m) / s))
    K = np.exp(-(g2 - g1) ** 2 / (2 * sig * sig))
    return np.sum(ei[rep] * np.conj(ej[j2]) * K * W) / (hi - lo)


if __name__ == '__main__':
    I, J = np.triu_indices(nc, 1)
    out = {}
    for sig in (0.5, 1.0):
        for (lo, hi, tag) in ((20.0, 396.0, ''), (20.0, 208.0, '_lo'), (208.0, 396.0, '_hi')):
            if sig == 1.0 and tag:
                continue
            X = np.array([Xpair(i, j, sig, lo, hi) for i, j in zip(I, J)])
            out['X_%s%s' % (sig, tag)] = X
            print('X sigma', sig, tag, 'done', time.time() - t0)
    np.savez(SC + 'rv_zeros.npz', ids=ids, I=I, J=J, rat=np.array(RAT), ctrl=CTRL,
             **{'G_' + k: v for k, v in G.items()}, **{'GC_' + k: v for k, v in GC.items()}, **out)
    print('saved', time.time() - t0)
