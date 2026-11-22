"""Own discrete leading-order law P(a,b,c) ~ S({0,a,a+b,a+b+c}) e^(-(a+b+c)/lam), even a,b,c >= 2, a+b+c <= G.
Singular series from primes <= G only (primes > G give a constant factor, irrelevant after normalisation). CMI in bits."""
import sys, math
import numpy as np
import sympy

def law_cmi(lam, G):
    ev = np.arange(2, G + 1, 2)
    A, B, C = np.meshgrid(ev, ev, ev, indexing='ij'); m = (A + B + C) <= G
    A, B, C = A[m], B[m], C[m]; del m
    pts = [np.zeros_like(A), A, A + B, A + B + C]
    logS = np.zeros(len(A))
    for p in sympy.primerange(3, G + 1):
        r = np.stack([x % p for x in pts], 1); r.sort(1)
        nu = 1 + np.count_nonzero(np.diff(r, axis=1), axis=1)
        val = 1 - nu / p
        logS += np.where(val > 0, np.log(np.maximum(val, 1e-300)), -np.inf) - 4 * math.log(1 - 1 / p)
    w = np.exp(logS - (A + B + C) / lam); w /= w.sum()
    n = G // 2 + 1
    T = np.zeros((n, n, n)); T[A // 2, B // 2, C // 2] = w
    ab = T.sum(2); bc = T.sum(0); bb = ab.sum(0)
    def Hh(x):
        x = x[x > 0]; return -np.sum(x * np.log2(x))
    return Hh(ab) + Hh(bc) - Hh(bb) - Hh(T), float(np.exp(-G / lam))

MULTS = (12, 20)
if __name__ == "__main__":
    for lam in map(float, sys.argv[1].split(',')):
        for mult in MULTS:
            G = int(mult * lam) // 2 * 2
            j, _ = law_cmi(lam, G)
            print('lam %.4f G %d (%d lam) CMI bits %.6f' % (lam, G, mult, j), flush=True)
