"""Quadrature check: |L(1/2+it, chi_D)| at the midpoints t = 100.025 + k/20 (k = 0..17999) by PARI (own call,
realprecision 30), for the three characters of the triple with the largest trapezoid-Simpson discrepancy; the
h = 0.025 trapezoid of |L1 L2 L3|^2 vs the stored h = 0.05 one; Euler-Maclaurin endpoint term; spot values vs mpmath."""
import sys, os, math, json
import numpy as np, cypari2, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
pari = cypari2.Pari(); pari.allocatemem(400 * 10**6); pari.set_real_precision(30)
Ds = [int(a) for a in sys.argv[1:4]]
Z = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/study/r149/reserved_values.npz')
mid = {}
for D in Ds:
    f = os.path.join(HERE, 'mid_%d.npy' % D)
    if os.path.exists(f):
        mid[D] = np.load(f); continue
    pari('Li = lfuninit(lfuncreate(%d), [1/2, 0, 1001])' % D)
    v = pari('vector(18000, k, abs(lfun(Li, 1/2 + I*(100 + (2*k-1)/40))))')
    mid[D] = np.array([float(x) for x in v]); np.save(f, mid[D])
    print('D', D, 'done', flush=True)
h = 0.05
y = np.ones(18001); ym = np.ones(18000)
for D in Ds:
    y *= Z['absL_%d' % D] ** 2; ym *= mid[D] ** 2
T1 = h * (y.sum() - 0.5 * (y[0] + y[-1]))
T2 = 0.5 * (T1 + h * ym.sum())
# Euler-Maclaurin endpoint correction -h^2/12 (f'(b) - f'(a)) estimated with the fine grid
hf = h / 2
yf = np.empty(36001); yf[0::2] = y; yf[1::2] = ym
fa = (-3 * yf[0] + 4 * yf[1] - yf[2]) / (2 * hf); fb = (3 * yf[-1] - 4 * yf[-2] + yf[-3]) / (2 * hf)
EM = -(hf ** 2 / 12) * (fb - fa)
print('triple %s: trap h=0.05 %.10g  trap h=0.025 %.10g  rel diff %.2e ; EM endpoint term (h=0.025) rel %.2e' % (Ds, T1, T2, T1 / T2 - 1, EM / T2))
# spot check of a few midpoint values by an independent Hurwitz sum (mpmath)
mp.mp.dps = 20
rng = np.random.default_rng(3)
for D in Ds:
    for k in rng.integers(0, 18000, 1):
        t = 100 + (2 * (k + 1) - 1) / 40
        q = abs(D); s = mp.mpc(0.5, t)
        val = abs(sum(int(pari.kronecker(D, a)) * mp.zeta(s, mp.mpf(a) / q) for a in range(1, q + 1)) * mp.mpf(q) ** (-s))
        print('  D %d t %.3f PARI %.12f Hurwitz %.12f' % (D, t, mid[D][k], float(val)))
