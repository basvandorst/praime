"""Puzzle 9: independent checks for writeup/9.md (plain Python + numpy + sympy, no C).

1. Every square / perfect power S(k) = p_1 + ... + p_k for k <= 10^6, from a direct sieve.
   The same for the odd-prime sums T(k) = S(k) - 2 (Seidov's variant).
2. The large records (McCranie 709838, Resta 126789311423, Seidov 55560480) with Lucy's
   prime-sum recursion, done twice: exactly mod 2^64 (uint64 wrap-around) and approximately
   in float64. A sum that agrees with n^2 mod 2^64 and lies within 2^63 of n^2 equals n^2.
Run from the repo root: .venv/bin/python writeup/code/p9_verify.py   (about a minute)
"""
from math import isqrt
import numpy as np
from sympy import isprime, integer_nthroot, primerange

# 1. Direct sieve: all perfect powers among S(k) and T(k), k <= 10^6.
K = 10**6
s, found_S, found_T = 0, [], []
for k, p in enumerate(primerange(2, 15485864), start=1):   # p_1000000 = 15485863
    s += p
    for m in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        if integer_nthroot(s, m)[1]:
            found_S.append((k, s, m))
        if k > 1 and integer_nthroot(s - 2, m)[1]:
            found_T.append((k, s - 2, m))
assert k == K
assert [(k, m) for k, _, m in found_S] == [(9, 2), (2474, 2), (6694, 2), (7785, 2), (709838, 2)], found_S
assert [(k, m) for k, _, m in found_T] == [(3, 3), (24, 2), (28, 2), (32, 2), (46, 2), (296, 2),
                                           (5327, 2), (12227, 2)], found_T
print("k <= 10^6: squares of S at k = 9, 2474, 6694, 7785, 709838; no cube or higher power")
print("k <= 10^6: T(3) = 8 = 2^3 and the seven squares of Seidov's list; nothing else")


# 2. Lucy's recursion for pi(x) and sum_{p <= x} p, vectorised with numpy.
def lucy(n):
    r = isqrt(n)
    V = [n // i for i in range(1, r + 1)]
    V += list(range(V[-1] - 1, 0, -1))
    V = np.array(V, dtype=np.int64)
    cnt = len(V)
    Vu = V.astype(np.uint64)
    C = V - 1                                                    # count of 2..v
    half = np.where(V % 2 == 0, (Vu // 2) * (Vu + 1), Vu * ((Vu + 1) // 2))
    Sm = half - np.uint64(1)                                     # sum of 2..v, mod 2^64
    Sf = V.astype(np.float64) * (V.astype(np.float64) + 1) / 2 - 1   # the same, in float64
    asc = V[::-1]

    def idx(x):
        return np.where(x <= r, cnt - x, n // np.maximum(x, 1) - 1)

    def one(x):
        return cnt - x if x <= r else n // x - 1

    for p in range(2, r + 1):
        if C[one(p)] == C[one(p - 1)]:
            continue                                             # p is not prime
        cp, smp, sfp = C[one(p - 1)], Sm[one(p - 1)], Sf[one(p - 1)]
        m = cnt - np.searchsorted(asc, p * p)                    # entries with v >= p^2
        j = idx(V[:m] // p)
        C[:m] -= C[j] - cp
        Sm[:m] -= np.uint64(p) * (Sm[j] - smp)
        Sf[:m] -= p * (Sf[j] - sfp)
    return int(C[0]), int(Sm[0]), float(Sf[0])


def check(x, k, total):
    """pi(x) = k, x prime, and sum_{p <= x} p = total."""
    pi, sm, sf = lucy(x)
    assert isprime(x) and pi == k, (x, pi, k)
    assert sm == total % 2**64 and abs(sf - total) < 2**63, (x, sm, sf)


check(10729219, 709838, 1916357**2)
print("McCranie: pi(10729219) = 709838, sum = 1916357^2")
check(1097877707, 55560480, 172287403**2 + 2)
print("Seidov:   pi(1097877707) = 55560480, sum of odd primes = 172287403^2")
check(3531577135439, 126789311423, 468726713734**2)
assert 468726713734**2 == 219704732167875184222756
print("Resta:    pi(3531577135439) = 126789311423, sum = 468726713734^2 = 219704732167875184222756")
