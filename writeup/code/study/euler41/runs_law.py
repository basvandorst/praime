"""Round 103 law: expected number of n in [N0, N1) with f(n), ..., f(n+r-1) all prime, f(n) = n^2 + n + 41
(Bateman-Horn for the pattern {0..r-1}): S_r * int dn / prod_i ln f(n+i), S_r = prod_p (1 - nu_p/p)/(1 - 1/p)^r,
nu_p = |{n mod p : p | prod_{i<r} f(n+i)}|. Primes p up to PMAX exactly; beyond, the factor is taken as 1 (its log
averages zero since nu_p has mean r). usage: runs_law.py"""
import math, json, sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import common as C

PMAX = 2_000_000; N0, N1 = 10 ** 3, 10 ** 7; RS = range(1, 11)


def roots(p):
    # n^2 + n + 41 = 0 mod p  <=>  (2n+1)^2 = -163 mod p
    if p == 2:
        return [n for n in range(2) if (n * n + n + 41) % 2 == 0]
    d = (-163) % p
    if d == 0:
        return [((p - 1) * pow(2, -1, p)) % p]
    if pow(d, (p - 1) // 2, p) != 1:
        return []
    s = int(sympy_sqrt(d, p)); inv2 = pow(2, -1, p)
    return sorted({((s - 1) * inv2) % p, ((-s - 1) * inv2) % p})


def sympy_sqrt(d, p):
    from sympy.ntheory import sqrt_mod
    return sqrt_mod(d, p)


if __name__ == '__main__':
    P = C.small_primes(PMAX).tolist()
    logS = {r: 0.0 for r in RS}
    for p in P:
        R = roots(p)
        for r in RS:
            nu = len({(x - i) % p for x in R for i in range(r)}) if R else 0
            logS[r] += math.log1p(-nu / p) - r * math.log1p(-1 / p)
    u = np.linspace(N0, N1, 400001); out = {}
    for r in RS:
        dens = np.prod([1 / np.log(u + i) ** 0 for i in range(1)], axis=0)
        integrand = np.ones_like(u)
        for i in range(r):
            integrand = integrand / np.log((u + i) ** 2 + (u + i) + 41)
        I = float(np.trapezoid(integrand, u))
        out[str(r)] = {'S_r': math.exp(logS[r]), 'expected': math.exp(logS[r]) * I}
        print(r, 'S_r %.4f' % math.exp(logS[r]), 'expected count %.1f' % (math.exp(logS[r]) * I), flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs_law.json'), 'w'), indent=1)
