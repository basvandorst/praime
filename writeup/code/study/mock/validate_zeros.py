#!/usr/bin/env python3
"""validate_zeros.py -- the explicit formula used by mock_sum.py, checked against true zeros (PARI lfunzeros).

For every non-principal character chi mod q (q = 101, 151), test function f with fhat = T_b * T_b (T_b the triangle
(1-|u|/b)_+, b = s/2; support of fhat = [-s, s]; f(x) = (b sinc^2(pi b x))^2 decays like x^-4):
   zero side:  W(chi) = sum over zeros 1/2 + i gamma, |gamma| <= T, of f(gamma log q / 2 pi)
   prime side: Wbar_a + P(chi),  Wbar_a = (1/L) int f(x) [log(q/pi) + Re psi(1/4 + a/2 + i pi x/L)] dx,
               P(chi) = -(2/L) sum_{n <= q^s} Lambda(n) n^-1/2 fhat(log n/L) Re chi(n)      (L = log q, a = parity)
Independent of mock_sum.py: primes, characters (chareval) and zeros all from PARI; digamma from scipy.
Also checks Wbar_odd - Wbar_even against pi*fhat(0)/L-type size (printed).
"""
import math, sys
import numpy as np
from scipy import integrate, special
import cypari2

pari = cypari2.Pari()
pari.allocatemem(500 * 10 ** 6)


def fhat_TT(u, s):
    b = s / 2.0
    t = abs(u) / b
    if t <= 1:
        return b * (2.0 / 3 - t * t + t ** 3 / 2)
    if t <= 2:
        return b * (2 - t) ** 3 / 6
    return 0.0


def f_TT(x, s):
    b = s / 2.0
    y = math.pi * b * x
    F = b if abs(y) < 1e-12 else b * (math.sin(y) / y) ** 2
    return F * F


def wbar(q, s, a):
    L = math.log(q)
    def g(x):
        return f_TT(x, s) * (math.log(q / math.pi) + special.psi(0.25 + a / 2 + 1j * math.pi * x / L).real)
    tot = 0.0
    edges = [0, 1, 2, 5, 10, 20, 50, 100, 400, 2000]
    for lo, hi in zip(edges[:-1], edges[1:]):
        v, _ = integrate.quad(g, lo, hi, limit=400, epsabs=1e-13, epsrel=1e-12)
        tot += v
    return 2 * tot / L


def main():
    T = 150.0
    worst = 0.0
    for q in (101, 151):
        L = math.log(q)
        G = pari(f'znstar({q},1)')
        for s in (1.0, 2.0):
            X = int(q ** s)
            pp = [int(p) for p in pari.primes([2, X])]
            plist = []  # (n, Lambda, fhat)
            for p in pp:
                if p == q:
                    continue
                k = 1
                while p ** k <= X:
                    plist.append((p ** k, math.log(p), fhat_TT(k * math.log(p) / L, s)))
                    k += 1
            wb = {0: wbar(q, s, 0), 1: wbar(q, s, 1)}
            diffs = []
            for k in range(1, q - 1):
                chi = pari(f'[{k}]')
                a = 0 if k % 2 == 0 else 1
                Ps = 0.0
                for n, lam, fh in plist:
                    if fh == 0.0:
                        continue
                    e = float(pari.chareval(G, chi, n))
                    Ps += lam / math.sqrt(n) * fh * math.cos(2 * math.pi * e)
                Ps *= -2.0 / L
                Lf = pari.lfuncreate(pari(f'[znstar({q},1),[{k}]]'))
                zs = [float(z) for z in pari.lfunzeros(Lf, pari(f'[-{T},{T}]'))]
                W = sum(f_TT(g * L / (2 * math.pi), s) for g in zs)
                d = W - (wb[a] + Ps)
                diffs.append(d)
            diffs = np.array(diffs)
            worst = max(worst, np.max(np.abs(diffs)))
            print(f'q={q} s={s}: {q-2} characters, Wbar_even={wb[0]:.6f} Wbar_odd={wb[1]:.6f} '
                  f'(fhat(0)={fhat_TT(0, s):.4f}); max |zeros - (Wbar+P)| = {np.max(np.abs(diffs)):.2e}, '
                  f'mean = {np.mean(diffs):+.1e}', flush=True)
    tail = 2 * (s / 2) ** 2 / (math.pi * s / 2) ** 4 / (3 * (T * math.log(101) / (2 * math.pi)) ** 3)
    print(f'tail bound of the truncated zero sum (per unit density, x^-4 envelope): {tail:.1e}')
    print('PASS' if worst < 1e-6 else 'FAIL')


if __name__ == '__main__':
    main()
