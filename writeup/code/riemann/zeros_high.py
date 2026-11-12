"""H4 at height 2.68e11: the zeros near zero #10^12 from the prime powers up to X = c sqrt(T/2 pi).

Same fit as primes_to_zeros.py (smoothed steps Kc(ln X (t - gamma)), least squares), scored
against Odlyzko's table (zeros 10^12 + 1 ... 10^12 + 10^4, accurate to 1e-8). The heights are
t = B + delta with B = 267653395647; B ln n mod 2 pi is taken in mpmath, theta(B + delta) - theta(B)
by its Taylor series, and theta(B)/pi + 1 mod 1 with mpmath.siegeltheta.
For H5 it also prints the error of each zero against its gap to the nearest neighbour.

usage: python3 zeros_high.py C1,C2,... K0,K0,...   (windows of 60 zeros starting at index K0)
"""
import os
import re
import sys
import time

import mpmath as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import primes_to_zeros as pz  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
txt = open(os.path.join(HERE, '..', '..', 'data', 'riemann', 'odlyzko_zeros3.txt')).read().splitlines()
B = int(re.search(r'gamma - (\d+)', ' '.join(txt[:6])).group(1))
Z = np.array([float(s) for s in txt if re.fullmatch(r'\s*\d+\.\d+\s*', s)])
mp.mp.dps = 40
OFF = float(mp.frac(mp.siegeltheta(B) / mp.pi + 1))


def dtheta(d):
    Bf = float(B)
    return 0.5 * np.log(Bf / (2 * np.pi)) * d + d * d / (4 * Bf) - d ** 3 / (12 * Bf * Bf)


class HighPrimes:
    def __init__(self, X):
        n, l = pz.lam_table(int(X))
        self.U = np.log(X)
        self.ln = np.log(n)
        self.coef = l * pz.wfun(self.ln / self.U) / (np.sqrt(n) * self.ln)
        two = 2 * mp.pi
        self.phi = np.array([float(mp.fmod(B * mp.log(int(m)), two)) for m in n])
        self.count = len(n)

    def N(self, d):
        s = np.zeros_like(d)
        for i in range(0, len(d), 2000):
            s[i:i + 2000] = -np.sin(self.phi[None, :] + np.outer(d[i:i + 2000], self.ln)) @ self.coef
        return OFF + dtheta(d) / np.pi + s / np.pi


def window(P, k0, nz=60, mg=25):
    g = Z[k0:k0 + nz]
    dd = 2 * np.pi / np.log((B + g[0]) / (2 * np.pi))
    t = np.arange(g[0] - (mg + 6) * dd, g[-1] + (mg + 6) * dd, min(dd / 10, np.pi / (4 * P.U)))
    N = P.N(t)
    est = pz.crossings(t, N)
    p0 = est[(est > g[0] - mg * dd) & (est < g[-1] + mg * dd)]
    sm = (t > g[0] - (mg - 5) * dd) & (t < g[-1] + (mg - 5) * dd)
    p = pz.fit(t[sm], N[sm], p0, P.U)
    e1 = np.array([np.min(np.abs(p - z)) for z in g])
    e0 = np.array([np.min(np.abs(est - z)) for z in g])
    gaps = np.minimum(Z[k0 + 1:k0 + nz + 1] - g, g - Z[k0 - 1:k0 + nz - 1]) / dd
    return e0, e1, gaps


if __name__ == '__main__':
    root = np.sqrt((B + 1300) / (2 * np.pi))
    cs = [float(x) for x in sys.argv[1].split(',')]
    starts = [int(s) for s in sys.argv[2].split(',')]
    for c in cs:
        allE, allG = [], []
        t0 = time.time()
        P = HighPrimes(c * root)
        rows = []
        for k0 in starts:
            e0, e1, gaps = window(P, k0)
            allE += list(e1); allG += list(gaps)
            j = int(np.argmax(e1))
            rows.append(f"zeros {k0}-{k0 + 59}: max {e1.max():.1e}, rms {np.sqrt(np.mean(e1 ** 2)):.1e} "
                        f"(crossings {e0.max():.1e}); worst zero's nearest gap {gaps[j]:.2f} mean spacings")
        print(f"c = {c}: X = {c * root:.0f}, {P.count} prime powers [{time.time() - t0:.0f} s]\n   "
              + "\n   ".join(rows), flush=True)
        E, Gp = np.array(allE), np.array(allG)
        for lo, hi in [(0, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 9)]:
            m = (Gp >= lo) & (Gp < hi)
            if m.any():
                print(f"nearest gap {lo}-{hi} spacings: {m.sum()} zeros, median error {np.median(E[m]):.1e}, max {E[m].max():.1e}")
