"""Puzzle 1121: R_m = phi(F_m)/phi(F_m - 2) for Fermat numbers F_m = 2^(2^m) + 1.

Claim: R_m = 2 for m <= 4, R_5 = 261735/131072, and R_m > 2 for every m >= 6 (whether F_m is prime or not).
So "F_m prime <=> R_m = 2" holds exactly when no Fermat prime beyond 65537 exists; R_m = 2 forces m <= 4.

Proof of R_m > 2 for m >= 6 (checked here):
  F_m - 2 = F_0 F_1 ... F_(m-1) (pairwise coprime), so phi(F_m - 2)/(F_m - 2) <= [2^31/(2^32 - 1)] * phi(F_5)/F_5 = D.
  R_m = F_m/(F_m - 2) * N / (phi(F_m - 2)/(F_m - 2)) > N/D with N = phi(F_m)/F_m >= 1 - sum 1/q over primes q | F_m.
  Every prime q | F_m is 1 mod 2^(m+2) (Lucas), so the j-th smallest is >= j 2^(m+2) + 1 and there are at most 2^m of them:
  sum 1/q < (1 + m ln 2)/2^(m+2).  R_m > 2 iff N > 2D.  For m >= 11 the bound alone suffices; for m = 6..10 we
  find all prime factors below 2^32 by testing every q = k 2^(m+2) + 1 < 2^32, and bound the rest by (#rest)/2^32.
"""
from fractions import Fraction as Fr
from math import log
import numpy as np
from sympy import isprime, factorint, totient

F = lambda m: 2**(2**m) + 1
# small m exactly (F_0..F_5 factor instantly)
for m in range(6):
    R = Fr(int(totient(F(m))), int(totient(F(m) - 2)))
    print(f"m = {m}: R = {R} = {float(R):.6f}")
    if m <= 4: assert R == 2 and isprime(F(m))
assert Fr(int(totient(F(5))), int(totient(F(5) - 2))) == Fr(261735, 131072)
phi5 = 640 * 6700416
assert phi5 == int(totient(F(5)))
twoD = 2 * Fr(2**31, 2**32 - 1) * Fr(phi5, F(5))      # R_m > 2 iff N > 2D (and N > 2D is implied by 1 - S > 2D)
thr = 1 - twoD                                        # need sum 1/q < thr
print(f"need sum 1/q over primes q | F_m below {float(thr):.6e}")
# general m >= 11
for m in range(11, 200):
    assert (1 + m * log(2)) / 2**(m + 2) < float(thr)
print("m >= 11: (1 + m ln 2)/2^(m+2) < threshold, so R_m > 2 (the bound decreases with m)")
# m = 6..10: all prime factors q < 2^32, q = k 2^(m+2) + 1
B = 2**32
for m in range(6, 11):
    step = 2**(m + 2)
    k = np.arange(1, B // step + 1, dtype=np.uint64); q = k * np.uint64(step) + np.uint64(1)
    q = q[q < B]
    x = np.full(len(q), 2, dtype=np.uint64)
    for _ in range(m): x = (x * x) % q                 # 2^(2^m) mod q, products < 2^64
    found = [int(v) for v in q[x == q - 1] if isprime(int(v))]
    S = sum(Fr(1, v) for v in found)
    rest_bits = 2**m + 1 - sum(log(v, 2) for v in found)   # log2 of the cofactor bound
    tail = Fr(int(rest_bits / 32) + 1, B)              # at most rest_bits/32 further factors, each >= 2^32
    assert S + tail < thr
    print(f"m = {m}: prime factors below 2^32: {found}; sum 1/q < {float(S + tail):.3e} -> R_m > 2")
# the page's values for m = 6, 7 (complete factorisations by PARI)
import cypari2
pari = cypari2.Pari()
def phi_F(i):
    ph = 1
    for p in [int(t) for t in pari.factor(F(i))[0]]:
        assert isprime(p); ph *= p - 1
    return ph
for m, val in ((6, Fr(18764930006193, 9367860543488)),
               (7, Fr(9615364552842273224968401793233, 4800163769175502470780823797760))):
    den = 1
    for i in range(m): den *= phi_F(i)                 # phi(F_m - 2) = prod phi(F_i), i < m
    R = Fr(phi_F(m), den)
    assert R == val and R > 2
    print(f"m = {m}: R = {float(R):.6f} (page value confirmed)")
print("all claims verified")
