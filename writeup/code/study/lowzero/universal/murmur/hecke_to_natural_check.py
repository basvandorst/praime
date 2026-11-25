"""Reverse direction: from the Hecke-twisted family lines (fact 2054 / Theorem D, harmonic weights, prime level)
to the natural-weight density at prime level.

Natural weights = harmonic weights x L(1, Sym^2 f) / zeta^{(N)}(2) = harmonic x sum_t lambda_f(t^2)/t, so
   A_nat(rho) = sum_{t,m >= 1, t^2/m^2 = rho} (1/t) * mu(m) c_m(t^2) / (phi(m) t)        [Theorem D amplitude at l = t^2]
must equal the line amplitude derived from the one-point law with lines + Plancherel (local_factors_check.py):
   A_prime((s/d)^2) = (6/pi^2) s^{-2} prod_{p|d} p^2/((p-1)^2 (p+1))   (d squarefree, (s,d) = 1; 0 if d not squarefree),
(the factor at p = N is ignored: N large).  Checked numerically for small s, d (sum over g = gcd-multiplier up to G)."""
from math import gcd, pi
import sympy as sp


def ramanujan(m, n):
    # c_m(n) = sum_{d | gcd(m,n)} mu(m/d) d
    g = gcd(m, n)
    return sum(sp.mobius(m // d) * d for d in sp.divisors(g))


def phi(m):
    return int(sp.totient(m))


def A_nat(s, d, G=3000):
    tot = 0.0
    for g in range(1, G + 1):
        t, m = s * g, d * g
        mu = sp.mobius(m)
        if mu == 0:
            continue
        tot += float(mu * ramanujan(m, t * t)) / (phi(m) * t * t)
    return tot


def A_prime(s, d):
    if any(e > 1 for e in sp.factorint(d).values()):
        return 0.0
    v = 6 / pi ** 2 / s ** 2
    for p in sp.factorint(d):
        v *= p * p / ((p - 1) ** 2 * (p + 1))
    return v


print(" s  d | from Theorem D (l = t^2), harmonic->natural | from one-point law + Plancherel | ratio")
for s, d in ((1, 1), (2, 1), (3, 1), (1, 2), (1, 3), (1, 6), (2, 3), (3, 2), (5, 6), (1, 4), (4, 1), (1, 5), (7, 10)):
    if gcd(s, d) != 1:
        continue
    a, b = A_nat(s, d), A_prime(s, d)
    print(f"{s:2d} {d:2d} | {a:.6f} | {b:.6f} | {a / b if b else float('nan'):.5f}")
