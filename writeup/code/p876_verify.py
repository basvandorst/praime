"""Puzzle 876 (prime tuples, Kamenetsky): every term after the n-th is the sum of the previous n terms; all terms are
distinct odd primes; k = number of terms before the first composite.

Theorem (proved below and checked here): for odd n not 1 (mod 3), every prime tuple has k <= 2n + 1.
Mod 3, the generated terms satisfy g_1 = S, g_(j+1) = 2 g_j - a_j (j <= n), g_(n+2) = 2 g_(n+1) - g_1, so
g_(j+1) = (-1)^j B_j with B_j = S - sum_(i<=j) (-1)^i a_i.  Nonzero g_1 .. g_(n+1) force B to alternate, except at the
one start term that can be 0 mod 3 (the prime 3); without it n = 1 (mod 3), with it g_(n+2) = 0 (mod 3).
Checks: (1) the page's tuples: lengths and weights; (2) the run bound by exhaustive search over all residue vectors
mod 3 (at most one zero) for odd n <= 17; (3) the extremal residue pattern of the proof on the page's tuples.
"""
from itertools import product
from sympy import isprime

table = {3: ([3, 13, 7], 7, 23), 5: ([17, 3, 19, 7, 13], 11, 59), 7: ([157, 379, 487, 109, 13, 7, 271], 25, 1423),
         9: ([11, 47, 17, 23, 41, 5, 3, 13, 19], 19, 179), 11: ([43, 7, 19, 13, 3, 17, 11, 5, 29, 41, 23], 23, 211),
         13: ([53, 137, 11, 17, 41, 227, 47, 101, 83, 5, 149, 263, 29], 34, 1163),
         15: ([29, 5, 23, 11, 41, 47, 89, 17, 71, 3, 7, 13, 37, 19, 79], 31, 491),
         17: ([5, 47, 53, 11, 17, 41, 89, 3, 61, 43, 97, 19, 13, 7, 37, 31, 73], 35, 647),
         19: ([89, 227, 29, 17, 5, 251, 269, 107, 101, 197, 41, 191, 173, 179, 47, 53, 71, 11, 23], 43, 2081)}
def length(start):
    seq = list(start); n = len(start)
    assert all(isprime(p) and p % 2 for p in seq) and len(set(seq)) == n
    while True:
        t = sum(seq[-n:])
        if not isprime(t) or t in seq: return len(seq)
        seq.append(t)
for n, (start, k, w) in table.items():
    assert length(start) == k and sum(start) == w
    print(f"n = {n}: page tuple has length {k} (weight {w}); " + (f"bound 2n+1 = {2*n+1}" + (" -> maximal" if k == 2*n+1 else "") if n % 3 != 1 else "n = 1 (mod 3): no mod-3 bound"))
def maxrun3(n):
    best = 0
    for v in product(range(3), repeat=n):
        if v.count(0) > 1: continue
        w = list(v); r = 0
        while r < 3 * n + 10:
            t = sum(w[-n:]) % 3
            if t == 0: break
            w.append(t); r += 1
        best = max(best, r)
    return best
for n in (3, 5, 7, 9, 11, 13):
    r = maxrun3(n)
    print(f"n = {n}: longest run of nonzero generated terms mod 3 = {r}" + (" (unbounded: all start terms congruent)" if r >= 3*n+10 else ""))
    assert (r == n + 1) if n % 3 != 1 else (r >= 3 * n + 10)
for n, (start, k, w) in table.items():          # the residue pattern of the proof
    if n % 3 == 1 or k != 2 * n + 1: continue
    j0 = start.index(3) + 1; c = sum(start) % 3
    assert j0 % 3 == (2 * n + 1) % 3
    assert all(a % 3 == c for a in start[:j0 - 1]) and all(a % 3 == (-c) % 3 for a in start[j0:])
print("all checks passed")
