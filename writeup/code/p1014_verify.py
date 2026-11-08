"""Puzzle 1014: for consecutive primes p < q < r, none of p | q^2 + r^2, q | p^2 + r^2, r | p^2 + q^2 holds
for 13 < p <= 4e18.  (about 2 s)

1. Mod p, q^2 + r^2 = g^2 + (g+h)^2 (g = q - p, h = r - q); mod q, p^2 + r^2 = g^2 + h^2; mod r, p^2 + q^2 = (g+h)^2 + h^2.
2. Direct check of the three conditions on every triple with p < 1.1e7.
3. Above that, gaps are <= 1476 (complete below 4e18), so the three residues are positive and below 1.09e7 < p.
"""
import numpy as np
from math import isqrt
N = 11_000_000
sv = np.ones(N + 400, bool); sv[:2] = False
for i in range(2, isqrt(N + 400) + 1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0].astype(np.int64)
p, q, r = P[:-2], P[1:-1], P[2:]
k = p < N; p, q, r = p[k], q[k], r[k]
A = (q*q + r*r) % p == 0; B = (p*p + r*r) % q == 0; C = (p*p + q*q) % r == 0
print("a) p | q^2 + r^2 for p =", p[A].tolist())
print("b) q | p^2 + r^2 for p =", p[B].tolist())
print("c) r | p^2 + q^2 for p =", p[C].tolist())
assert p[A].tolist() == [2, 5, 13] and not B.any() and not C.any()
g, h = q - p, r - q
x, y, z = g*g + (g+h)**2, g*g + h*h, (g+h)**2 + h*h
assert np.all(((q*q + r*r) - x) % p == 0) and np.all(((p*p + r*r) - y) % q == 0) and np.all(((p*p + q*q) - z) % r == 0)
# last p with x >= p, y >= q, z >= r (the page's thresholds 2477, 293, 2203)
print("largest p with x >= p:", p[x >= p].max(), "| largest q with y >= q:", q[y >= q].max(), "| largest r with z >= r:", r[z >= r].max())
G = 1476                                             # largest gap below 4e18
assert G*G + (2*G)**2 < N and 2*G*G < N
print("all claims verified")
