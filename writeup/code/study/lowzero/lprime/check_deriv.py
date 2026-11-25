"""Independent check of lprime_values.py on a sample of B: |L'(rho)| from
 (a) mpmath: L(s, chi_D) = q^-s sum_{a=1}^{q} chi_D(a) zeta(s, a/q), L' by mpmath's Hurwitz derivative (no PARI);
 (b) PARI values only (lfun of order 0, realprecision 50): 4-point central difference along t.
usage: check_deriv.py VALDIR ZEROFILE"""
import sys, os, math, numpy as np, mpmath as mp, cypari2
mp.mp.dps = 25
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6); pari('default(realprecision, 50)')
vdir, zf = sys.argv[1], sys.argv[2]


def kron(D, n):
    return int(pari('kronecker(%d, %d)' % (D, n)))


def L_and_dL(D, s):
    q = abs(D); v = mp.mpf(0); d = mp.mpf(0); lq = mp.log(q)
    for a in range(1, q):
        c = kron(D, a)
        if c == 0: continue
        h = mp.zeta(s, mp.mpf(a) / q); hd = mp.zeta(s, mp.mpf(a) / q, 1)
        v += c * h; d += c * (hd - lq * h)
    return q ** (-s) * v, q ** (-s) * d


rows = []
for D in (401, -404, 520, -643):
    fn = os.path.join(vdir, 'D%d.npz' % D)
    if not os.path.exists(fn): continue
    z = np.load(fn); g = z['gamma']; pick = [0, len(g) // 2, len(g) - 1]
    Li = pari('lfuninit(lfuncreate(%d), [1/2, 0, 605])' % D)
    for j in pick:
        x = g[j]; s = mp.mpc(0.5, x)
        Lv, dL = L_and_dL(D, s)
        h = 1e-4; f = lambda u: complex(pari.lfun(Li, pari('1/2 + I*(%.17g)' % u)))
        fd = (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h)   # dL/dt = i L'
        rows.append((D, x, float(z['logabs'][j]), math.log(abs(complex(dL))), math.log(abs(fd)), abs(complex(Lv))))
        print('D=%5d gamma=%9.4f  log|L\'| PARI(lfun,1) %.12f  mpmath-Hurwitz %.12f  PARI finite-diff %.12f   |L(rho)| mpmath %.1e' % rows[-1], flush=True)
r = np.array([x[2:] for x in rows])
print('max |PARI - mpmath| in log|L\'|: %.2e;  max |PARI - finite difference|: %.2e;  n = %d' % (np.abs(r[:, 0] - r[:, 1]).max(), np.abs(r[:, 0] - r[:, 2]).max(), len(r)))
