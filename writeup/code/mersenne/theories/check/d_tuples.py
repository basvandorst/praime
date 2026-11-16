# D second method: enumerate p = 23 mod 24 (p = 24j - 1) with p + offsets all prime,
# by residue filtering mod primes 5..997 then gmpy2.is_prime on survivors (no sieve array).
# Also the Hardy-Littlewood prediction C7 * sum_j 1/ln(24j+9)^7, and E2.
import numpy as np
import gmpy2
from math import log, isqrt
from sympy import primerange, isprime

JMAX = 41_666_666
PATS = {223: (0, 2, 8, 12, 14, 18, 20),   # 24j-1, +1, +7, +11, +13, +17, +19 (24j+5 composite)
        251: (0, 2, 6, 8, 12, 18, 20)}    # 24j-1, +1, +5, +7, +11, +17, +19 (24j+13 composite)
FILT = [r for r in primerange(5, 1000)]

found = {223: [], 251: []}
CH = 2_000_000
for v, offs in PATS.items():
    # p+off composite for the excluded slot: 223 needs 24j+5 = p+6 composite; 251 needs 24j+13 = p+14 composite
    excl = 6 if v == 223 else 14
    for j0 in range(1, JMAX + 1, CH):
        j = np.arange(j0, min(j0 + CH, JMAX + 1), dtype=np.int64)
        p = 24 * j - 1
        for r in FILT:
            pr = p % r
            keep = np.ones(len(p), dtype=bool)
            for o in offs:
                keep &= ((pr + o) % r != 0) | (p + o == r)
            j = j[keep]
            p = p[keep]
        for jj in j:
            jj = int(jj)
            pp = 24 * jj - 1
            if all(gmpy2.is_prime(pp + o) for o in offs) and not gmpy2.is_prime(pp + excl):
                found[v].append(jj)
    print("byte", v, "count", len(found[v]), "first j", found[v][:1],
          "<= 416666:", sum(1 for x in found[v] if x <= 416666))
    print("   all j:", found[v])

# Hardy-Littlewood: forms 24j + o
FORMS = {223: (-1, 1, 7, 11, 13, 17, 19), 251: (-1, 1, 5, 7, 11, 17, 19)}
rlist = list(primerange(2, 5_000_000))


def C7(forms):
    lg = 0.0
    for r in rlist:
        if r in (2, 3):
            nu = 0
        else:
            nu = len({o % r for o in forms})
        lg += log(1 - nu / r) - 7 * log(1 - 1 / r)
    # tail beyond 5e6: nu = 7, log factor ~ -21/r^2 ; sum_{r>R} 1/r^2 ~ 1/(R ln R)
    R = 5e6
    lg += -21.0 / (R * log(R))
    return np.exp(lg)


def Ssum(J):
    s = 0.0
    for j0 in range(1, J + 1, 5_000_000):
        j = np.arange(j0, min(j0 + 5_000_000, J + 1), dtype=np.float64)
        s += np.sum(np.log(24 * j + 9) ** -7.0)
    return s


S9 = Ssum(JMAX)
S7 = Ssum(416_666)
for v, f in FORMS.items():
    c = C7(f)
    print(f"HL byte {v}: C7 = {c:.4f}  pred to 1e9 = {c * S9:.3f}  pred to 1e7 = {c * S7:.3f}")

# E2
lim = isqrt(10**9 + 50)
N = sum(1 for p in primerange(5, lim + 1) if isprime(p * p - 2))
print("E2 primes 5 <= p <=", lim, "with p^2-2 prime:", N)
print("10^9+7, 10^9+9 prime:", isprime(10**9 + 7), isprime(10**9 + 9))
