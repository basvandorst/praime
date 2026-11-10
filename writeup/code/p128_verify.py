# Puzzle 128, independent check with Python integers (math.isqrt, a dict of prefix sums; no numpy).
# (a) no square among p_1^2 + ... + p_n^2 with p_1 = 2 or 3 and p_n <= 1e6.
# (c) adjacent windows with equal power sums, first window ending at p_j, j <= 3000: 259 solutions for first
#     powers (as in p128_search.py), none for squares; the page's non-adjacent example is recomputed.
# Run from the repo root: .venv/bin/python3 writeup/code/p128_verify.py
from math import isqrt
from sympy import primerange
P = list(primerange(2, 10**6))
for start in (0, 1):
    s = 0; n = 0
    for idx, p in enumerate(P[start:]):
        s += p * p
        if idx >= 1:
            n += 1; assert isqrt(s) ** 2 != s, (start, p)
    print(f"(a) p_1 = {P[start]}: {n} sums with last prime <= {P[-1]}, no square")
def adjacent(k, J):
    Q = [0]
    for p in P: Q.append(Q[-1] + p ** k)
    where = {q: l for l, q in enumerate(Q)}
    assert 2 * Q[J] <= Q[-1]
    return [(P[i], P[j - 1], P[j], P[where[2 * Q[j] - Q[i]] - 1]) for j in range(1, J + 1) for i in range(j) if 2 * Q[j] - Q[i] in where]
s1 = adjacent(1, 3000); s2 = adjacent(2, 3000)
assert len(s1) == 259 and s2 == []
print(f"(c) j <= 3000: first powers {len(s1)} solutions, e.g. {s1[1]} (first..last of each window); squares {len(s2)}")
A = sum(p * p for p in primerange(131, 648)); B = sum(p * p for p in primerange(941, 1034))
assert A == B
print(f"page: 131^2 + ... + 647^2 = 941^2 + ... + 1033^2 = {A} (not adjacent)")
assert sum(p * p for p in primerange(41, 174)) == 586 ** 2
print("all checks passed")
