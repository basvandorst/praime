"""Independent L'(1/2 + i gamma, chi_D) through the Hurwitz zeta function (mpmath), 5 zeros each of 10 L-functions,
against the stored log|L'| values."""
import os, numpy as np, mpmath as mp, sys, time
from sympy import jacobi_symbol
mp.mp.dps = 25
R = os.path.dirname(os.path.abspath(__file__))+'/../../../../'


def kron(D, n):
    # Kronecker symbol (D/n), n >= 1
    res = 1
    while n % 2 == 0:
        n //= 2
        if D % 2 == 0: return 0
        res *= 1 if D % 8 in (1, 7) else -1
    if n == 1: return res
    return res * int(jacobi_symbol(D % n, n)) if np.gcd(D, n) == 1 else 0


def lprime(D, g):
    q = abs(D); s = mp.mpc(0.5, g)
    tot = mp.mpc(0); totL = mp.mpc(0)
    for a in range(1, q):
        c = kron(D, a)
        if c == 0: continue
        totL += c * mp.zeta(s, mp.mpf(a) / q)
        tot += c * mp.zeta(s, mp.mpf(a) / q, 1)
    Lv = q ** (-s) * totL
    Ld = q ** (-s) * tot - mp.log(q) * Lv
    return Ld, Lv


Ds = [int(x) for x in sys.argv[1:]]
rng = np.random.default_rng(7)
for D in Ds:
    f = 'C' if abs(D) <= 400 else 'Bp'
    V = np.load(R + 'writeup/data/study/r146/vals_%s.npz' % f)
    g = V['D%d__gamma' % D]; la = V['D%d__logabs' % D]
    ix = np.sort(rng.choice(len(g), 5, replace=False))
    t0 = time.time(); out = []
    for i in ix:
        Ld, Lv = lprime(D, g[i])
        out.append((g[i], float(mp.log(abs(Ld))), la[i], float(abs(Lv) / abs(Ld))))
    print('D %5d' % D, ' '.join('[%.3f: %.10f vs %.10f d=%.1e |L/L\'|=%.0e]' % (a, b, c, b - c, e) for a, b, c, e in out),
          '%.0fs' % (time.time() - t0), flush=True)
