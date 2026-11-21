"""Round 2c: the scatter of prime counts along the quadratics 4n^2 + bn + c (Ulam's diagonals), from Bateman-Horn
constants of polynomial pairs. A Montgomery-Soundararajan-type law for a curve instead of an interval.

X_n = 1[f(n) prime], q_n = C / ln f(n)                                   (Bateman-Horn, conjecture)
E[X_n X_{n+d}] = C_d / (ln f(n) ln f(n+d)),  C_d = BH constant of the pair {f(x), f(x+d)}   (same conjecture)
Var(count) = sum_n q_n (1 - q_n) + 2 sum_{d>=1} (C_d - C^2) sum_n 1/(ln f(n) ln f(n+d))
E[(obs - E)^2 / E] = Var / E,  so the theory for G.chi2_BH is the mean over the curves of Var/E.

C_d / C^2 = prod_p (1 - omega_d(p)/p) / (1 - w(p)/p)^2,
  w(p) = #roots of f mod p, omega_d(p) = #(R u (R - d)) = 2w - #(R n (R - d)).
For p <= P0 everything by brute force (periodic in d). For p > P0 (p odd, p does not divide 4):
  generic omega = 2w;  p | d -> omega = w;  p | 16 d^2 - D, p not dividing dD -> omega = 3 (w = 2, one root shared).
Primes p > P1 = 30000 with p | 16d^2 - D change C_d by a factor (1-3/p)/(1-4/p) < 1 + 4e-5 and are ignored.
The generic product runs to 2e6 (the tail sum w^2/p^2 beyond it is below 1e-7).
"""
import math, os, json, sys
import numpy as np
from sympy.ntheory import sqrt_mod
from common import small_primes, WINDOWS
import probes as PR

P0, P1, PG = 200, 30000, 2_000_000
_PR = small_primes(PG)


def legendre_vec(D, ps):
    d = np.mod(D, ps); e = (ps - 1) // 2; r = np.ones_like(ps); base = d.copy()
    while np.any(e > 0):
        m = (e & 1) == 1
        r[m] = (r[m] * base[m]) % ps[m]
        base = (base * base) % ps
        e >>= 1
    return np.where(d == 0, 0, np.where(r == 1, 1, -1))


def ratio_Cd(f, dmax):
    """array R[d] = C_d / C^2 for d = 0..dmax (R[0] unused)."""
    a, b, c = f
    D = b * b - 4 * a * c
    d = np.arange(dmax + 1, dtype=np.int64)
    logR = np.zeros(dmax + 1)
    # small primes: brute force on d mod p
    for p in _PR[_PR <= P0]:
        p = int(p)
        roots = [x for x in range(p) if (a * x * x + b * x + c) % p == 0]
        w = len(roots)
        tab = np.empty(p)
        for t in range(p):
            om = len(set(roots) | set((r - t) % p for r in roots))
            tab[t] = math.log((1 - om / p) / (1 - w / p) ** 2) if om < p else -np.inf
        logR += tab[d % p]
    # generic part for p > P0
    ps = _PR[_PR > P0].astype(np.int64)
    w = 1 + legendre_vec(D, ps)
    gen = np.log1p(-2 * w / ps) - 2 * np.log1p(-w / ps)
    logR += gen.sum()
    # p | d  (p <= dmax): omega = w instead of 2w
    for p, wp, g in zip(ps[ps <= dmax], w[ps <= dmax], gen[ps <= dmax]):
        p = int(p)
        corr = math.log(1 / (1 - wp / p)) - g if wp < p else 0.0
        logR[p::p] += corr
    # p | 16 d^2 - D, p not dividing dD: omega = 3 instead of 4 (needs w = 2)
    for p, wp, g in zip(ps[ps <= P1], w[ps <= P1], gen[ps <= P1]):
        p = int(p)
        if wp != 2:
            continue
        x = (D * pow(16, -1, p)) % p
        rts = sqrt_mod(x, p, all_roots=True)
        corr = math.log((1 - 3 / p) / (1 - 2 / p) ** 2) - g
        for r in rts:
            if r % p == 0:
                continue
            logR[int(r)::p] += corr
    return np.exp(logR)


def curve_var_over_E(f, n, C, R):
    """Var/E of the count of primes among f(n) for the given consecutive n (array)."""
    a, b, c = f
    v = (a * n * n + b * n + c).astype(np.float64)
    wgt = 1.0 / np.log(v)
    q = C * wgt
    E = q.sum()
    N = len(n)
    var = float(np.sum(q * (1 - q)))
    if N > 1:
        cor = np.correlate(wgt, wgt, mode='full')[N:]      # lags 1..N-1: sum_n w_n w_{n+d}
        var += 2 * float(np.sum((R[1:N] * C * C - C * C) * cor))
    return var / E


def main():
    out = {}
    cache = {}
    curves = PR.curves()
    for win, (L, W) in WINDOWS.items():
        out[win] = {}
        for blk in range(-1, 10):
            Lb, Wb = (L, W) if blk < 0 else (L + blk * (W // 10), W // 10)
            vals = []
            for f in curves:
                n, v = PR._curve_ns(f, Lb, Wb)
                if len(n) == 0:
                    continue
                key = f
                if key not in cache or len(cache[key][1]) < len(n) + 1:
                    cache[key] = (PR.bh(f), ratio_Cd(f, max(len(n), 4000)))
                C, R = cache[key]
                vals.append(curve_var_over_E(f, n, C, R))
            out[win][str(blk)] = {'G.chi2_BH': float(np.mean(vals))}
        print(win, out[win]['-1'], flush=True)
    D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
    json.dump(out, open(os.path.join(D, 'theory_geo.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
