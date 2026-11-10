# Puzzle 1015: independent check in plain Python (sympy.factorint, recursive search).
# (1) The page's table of smallest starts per length (below 10^8): the longest sequence from each start has the
#     stated length; for lengths <= 12 the start is also the smallest among the semiprimes below 400 that are not
#     themselves a concatenation a&b of two primes (such a number continues the sequence of a*b).
# (2) Our starts for length 22: longest sequence exactly 22; one such sequence is printed.
# Run from the repo root: .venv/bin/python3 writeup/code/p1015_verify.py
import sys
from functools import lru_cache
from sympy import factorint, isprime
sys.setrecursionlimit(10000)
def pq(n):
    f = factorint(n)
    if sum(f.values()) != 2: return None
    ps = sorted(f); return (ps[0], ps[0]) if len(ps) == 1 else (ps[0], ps[1])
@lru_cache(maxsize=None)
def best(n):
    """longest sequence starting at n, as a tuple; () if n is not a semiprime"""
    f = pq(n)
    if f is None: return ()
    p, q = f
    a = best(int(str(p) + str(q))); b = best(int(str(q) + str(p))) if p != q else ()
    return (n,) + (a if len(a) >= len(b) else b)
table = {1: 6, 2: 4, 3: 74, 4: 121, 5: 123, 6: 34, 7: 15, 8: 95, 9: 26, 10: 10, 11: 391, 12: 178, 13: 1253, 14: 6274,
         15: 453, 16: 145, 17: 311507, 18: 78331, 19: 1354427, 20: 35239, 21: 28907, 23: 16064229}
for L, s in table.items(): assert len(best(s)) == L, (L, s, len(best(s)))
print("page table: all 22 starts have the stated length")
def has_predecessor(n):
    t = str(n)
    return any(t[i] != "0" and isprime(int(t[:i])) and isprime(int(t[i:])) for i in range(1, len(t)))
first = {}
for n in range(4, 400):
    L = len(best(n))
    if L and L not in first and not has_predecessor(n): first[L] = n
assert all(first[L] == table[L] for L in range(1, 13)), first
print("smallest start for lengths 1..12 (all n < 400):", [first[L] for L in range(1, 13)])
for s in (414428717, 978652727):
    seq = best(s); assert len(seq) == 22 and not has_predecessor(s)
    print(f"length 22 from {s}:")
    for x in seq: p, q = pq(x); print(f"   {x} = {p} * {q}")
print("all checks passed")
