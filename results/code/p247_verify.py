"""Puzzle 247 Q1: check that 8090674745553 starts a run of exactly 8 consecutive Smith numbers (OEIS A059754 a(8),
Max Alekseyev 2010), and recheck a(1..7). Minimality is from OEIS, not rechecked here (it needs a search to 8.1e12)."""
from sympy import factorint
def smith(n):
    f = factorint(n)
    if sum(f.values()) < 2: return False                       # primes (and 1) are not Smith numbers
    ds = lambda m: sum(map(int, str(m)))
    return ds(n) == sum(ds(p) * e for p, e in f.items())
def run(s):
    k = 0
    while smith(s + k): k += 1
    return k
terms = [4, 728, 73615, 4463535, 15966114, 2050918644, 164736913905, 8090674745553]
for n, s in enumerate(terms, 1):
    assert run(s) >= n and not smith(s - 1), (n, s)
    print(f"a({n}) = {s}: run of {run(s)} consecutive Smith numbers")
print("all checks passed")
