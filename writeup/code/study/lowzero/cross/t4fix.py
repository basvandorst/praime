"""(T4) repeated with a composite Gauss-Legendre rule (the 20-node rule of theory.py is not converged at w = 4.5):
(1/2pi^2) int_{-w}^{w} (w - |e|) Re D_psi(-i e) de, D_psi = (L'/L)'(1+z, psi) - B_psi(z), against the closed form."""
import os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from crosslib import *
import cypari2
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
P6 = primes_upto(10 ** 6); PS = primes_upto(10 ** 4)
def setup(D1, D2):
    d = int(pari.quaddisc(D1 * D2)); Lf = pari.lfuncreate(d)
    rem = [(int(p), int(pari.kronecker(d, int(p)))) for p in pari.factor(abs(D1 * D2))[0] if d % int(p) != 0]
    x = (kron_table([D1], P6)[0] * kron_table([D2], P6)[0]).astype(float); return d, Lf, rem, x
def logabsL(Lf, rem, s):
    v = complex(pari.lfun(Lf, s)); out = math.log(abs(v))
    for p, c in rem: out += math.log(abs(1 - c * complex(p) ** (-complex(s))))
    return out
def closed(D1, D2, w, Lf, rem, x):
    xs = x[:len(PS)]; s = 0.0
    for p, c in zip(PS, xs):
        if c == 0: continue
        for m in range(2, 60):
            if float(p) ** m > 1e16: break
            s += (1 / m - 1 / m ** 2) * (c ** m) * float(p) ** -m * (1 - math.cos(m * w * math.log(p)))
    return (logabsL(Lf, rem, 1) - 0.5 * (logabsL(Lf, rem, complex(1, w)) + logabsL(Lf, rem, complex(1, -w))) - s) / math.pi ** 2
def Dpsi(Lf, rem, x, s):
    v0 = complex(pari.lfun(Lf, s)); v1 = complex(pari.lfun(Lf, s, 1)); v2 = complex(pari.lfun(Lf, s, 2)); s = complex(s)
    out = v2 / v0 - (v1 / v0) ** 2
    for p, c in rem: out += -c * math.log(p) ** 2 * p ** (-s) / (1 - c * p ** (-s)) ** 2
    lp = np.log(P6.astype(float)); u = x * P6.astype(float) ** (-s); out -= np.sum(lp ** 2 * u ** 2 / (1 - u) ** 2)
    return out
xg, wg = np.polynomial.legendre.leggauss(16)
for D1, D2 in ((401, -403), (-995, 997)):
    d, Lf, rem, x = setup(D1, D2); w = 4.5; npc = 6; tot = 0.0
    for j in range(npc):
        a, b = w * j / npc, w * (j + 1) / npc; e = (a + b) / 2 + (b - a) / 2 * xg; we = (b - a) / 2 * wg
        vals = np.array([Dpsi(Lf, rem, x, complex(1, -t)).real for t in e]); tot += np.sum(we * (w - e) * vals)
    rec = 2 * tot / (2 * math.pi ** 2); cl = closed(D1, D2, w, Lf, rem, x)
    print('D1=%d D2=%d w=%.1f: recipe integral (6 x 16 nodes) %+.7f, closed form %+.7f, diff %+.1e' % (D1, D2, w, rec, cl, rec - cl), flush=True)
