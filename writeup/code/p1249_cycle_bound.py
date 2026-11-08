"""Puzzles 1249/1250: cycles of p -> GPD(a p + b) and the product of their cofactors.

On a cycle p_1 .. p_L (p_(L+1) = p_1) write a p_i + b = s_i p_(i+1). Then
    S = prod s_i = prod (a p_i + b)/p_(i+1) = prod (a + b/p_i),
an integer with a^L < S <= (a + b/m)^L, m = min p_i.  S is odd (a p + b is odd), and for a = 3, b = 2 also prime
to 3, so S >= a^L + 1 (a = 2) or S >= a^L + 2 (a = 3).  If every prime below X ends in a known cycle, a new
cycle has m > X, so a^L ((1 + b/(a X))^L - 1) >= 1 or 2; this gives a lower bound on its length L.
"""
from fractions import Fraction as Fr
from math import prod
cycles = {(2, 3): [[5, 13, 29, 61], [43, 89, 181, 73, 149]],
          (3, 2): [[5, 17, 53, 23, 71, 43, 131, 79, 239, 719, 127, 383, 1151, 691, 83, 251, 151, 13, 41], [167, 503, 1511, 907, 389]]}
def gpd(n):
    g, d = 1, 2
    while d * d <= n:
        while n % d == 0: g, n = d, n // d
        d += 1
    return max(g, n) if n > 1 else g
for (a, b), cs in cycles.items():
    for c in cs:
        L = len(c)
        assert all(gpd(a * c[i] + b) == c[(i + 1) % L] for i in range(L))           # it is a cycle
        s = [(a * c[i] + b) // c[(i + 1) % L] for i in range(L)]
        assert prod(s) == prod(Fr(a * p + b, p) for p in c)                          # S = prod (a + b/p)
        print(f"a={a}, b={b}: cycle of length {L}, cofactors {s}, S = {prod(s)} = {a}^{L} + {prod(s) - a**L}")
def minlen(a, b, X, gap):
    L = 1
    while a**L * ((1 + Fr(b, a * X))**L - 1) < gap: L += 1
    return L
for (a, b, gap) in ((2, 3, 1), (3, 2, 2)):
    for X, src in ((2 * 10**8, "our check"), (10**10 if a == 2 else 61 * 10**11, "page search")):
        print(f"a={a}, b={b}: if every prime < {X:.3g} ends in a known cycle ({src}), a new cycle has length >= {minlen(a, b, X, gap)}")
print("all checks passed")
