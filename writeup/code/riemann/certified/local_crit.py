"""The local criterion: one odd n, no other prime known.

Kernel support a = 1.9 < n ln(1 + 2/n): only n and the even neighbours n +- 1 lie inside.
Lambda(n +- 1) is ln 2 if n +- 1 is a power of 2, else 0 (no prime needed to know this).
Lhat(n) = (zero side - even neighbours) / K(0); certified iff delta / K(0) < (1/2) ln n.
"""
import math, sys, time
import numpy as np
from kbm import K, K0
from cert import sample
from sympy import isprime, perfect_power

import os
S = os.environ.get('WORK', os.getcwd()) + '/'
ZFILE = os.environ.get('ZFILE', S + 'z2m.f64')
G = np.load(ZFILE) if ZFILE.endswith('.npy') else np.fromfile(ZFILE)
EPS = float(os.environ.get('EPS', '1e-5'))

def pow2(m):
    return m > 0 and m & (m - 1) == 0

def local(n, a, beta, m, c):
    T = c * beta * n / a
    z, d, info = sample(float(n), G, T, a, beta, m, EPS)
    for e in (n - 1, n + 1):
        if pow2(e):
            z -= math.log(2) * math.sqrt(n / e) * float(K(n * math.log(e / n), a, beta, m))
    k0 = K0(beta, m)
    return z / k0, d / k0, T, info

def best_T(n, a=1.9):
    """Smallest certified T over a grid of (beta, m, c)."""
    best = None
    for m in (2, 3, 4, 5, 6):
        for beta in np.arange(4.0, 30.01, 1.0):
            for c in (1.02, 1.05, 1.1, 1.2, 1.35, 1.5, 1.75, 2.0, 2.5):
                T = c * beta * n / a
                if T > G[-1] or (best and T >= best[0]):
                    continue
                L, d, T, info = local(n, a, beta, m, c)
                if d < 0.5 * math.log(n):
                    if best is None or T < best[0]:
                        best = (T, m, beta, c, d, L)
                    break
    return best

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'best':
        for n in [1009, 1011, 10007, 10009, 100003, 100005]:
            t0 = time.time(); b = best_T(n)
            T, m, beta, c, d, L = b
            print(f'n = {n}: certified with T = {T:.0f} = {T/n:.3f} n = {T/(math.pi*n):.3f} pi n (m = {m}, beta = {beta}, c = {c}); delta = {d:.3f} < {0.5*math.log(n):.3f}; Lhat = {L:.4f}, Lambda = {math.log(n) if isprime(n) else 0:.4f} [{time.time()-t0:.0f} s]')
    elif mode == 'range':
        lo, hi, m, beta, c = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
        t0 = time.time(); stats = dict(cert=0, uncert=0, wrong=0, skip=0); worst_err = 0; worst_d = 0; worst_ratio = 0
        for n in range(lo | 1, hi, 2):
            if perfect_power(n):
                stats['skip'] += 1; continue
            L, d, T, info = local(n, 1.9, beta, m, c)
            truth = math.log(n) if isprime(n) else 0.0
            err = abs(L - truth); worst_err = max(worst_err, err); worst_d = max(worst_d, d)
            assert err <= d, (n, L, truth, d)   # the bound must hold
            if d < 0.5 * math.log(n):
                stats['cert'] += 1
                if (L > 0.5 * math.log(n)) != (truth > 0): stats['wrong'] += 1
            else:
                stats['uncert'] += 1
        print(f'odd n in [{lo}, {hi}), m = {m}, beta = {beta}, T = {c} beta n / 1.9 = {c*beta/1.9:.3f} n: {stats}; largest |Lhat - Lambda| = {worst_err:.2e}, largest delta = {worst_d:.3f} [{time.time()-t0:.0f} s]')
