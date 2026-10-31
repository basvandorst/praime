"""Puzzle 221 (A153136): greedy primes such that no sum of at least two terms is prime.
Checks a(1..10) (a(10) = 42939980423, Donovan Johnson 2008, OEIS): every new term p makes p + S composite for every
nonempty subset sum S of the earlier terms, and for a(2..9) that no smaller prime works (full greedy recomputation).
Then the heuristic: P(p works) ~ (1 - 2/ln p)^E, E = number of nonempty subsets with even sum (odd p + odd S is even),
first solution where the expected count reaches 1; tested on a(10), applied to a(11)."""
import math
from itertools import combinations
from sympy import isprime, nextprime
seq = [2, 7, 13, 43, 103, 1627, 25349, 315743, 7338823, 42939980423]
def sums(t):
    s = {0}
    for x in t: s |= {y + x for y in s}
    return sorted(s - {0})
for k in range(1, len(seq)):
    S = sums(seq[:k]); p = seq[k]
    assert isprime(p) and all(not isprime(p + x) for x in S), k
    if k <= 8:                                   # greedy minimality, cheap up to a(9)
        q = nextprime(seq[k - 1])
        while q < p:
            assert any(isprime(q + x) for x in S), (k, q); q = nextprime(q)
print("a(1..10) satisfy the rule; a(2..9) are the greedy minimum")
def first(prev, lo):
    E = sum(1 for x in sums(prev) if x % 2 == 0)
    x, acc = lo, 0.0
    while acc < 1:                               # integrate density (1/ln x) * (1 - 2/ln x)^E
        lx = math.log(x); acc += x * 0.01 / lx * (1 - 2 / lx) ** E; x *= 1.01
    return E, x
E, x = first(seq[:9], seq[8]); print(f"a(10): {E} even subset sums; heuristic first solution near {x:.2g} (real 4.29e10)")
E, x = first(seq, seq[9]); print(f"a(11): {E} even subset sums; heuristic first solution near {x:.2g}")
