# Puzzle 58 / 54: runs of k consecutive primes p with 4p^2 + 1 prime.  Reads the candidate runs printed by
# p58_search.c (data/p58_candidates.txt: consecutive primes passing the small-prime filter on 4p^2 + 1), tests
# 4p^2 + 1 exactly (PARI isprime, 25-digit numbers) for each member, and reports the runs with >= K true members
# in increasing order of the start.  Also checks that the members are consecutive primes.
# Usage: p58_verify.py K           Run from the repo root.
import re, sys
from sympy import isprime, nextprime
import cypari2
pari = cypari2.Pari()
K = int(sys.argv[1]) if len(sys.argv) > 1 else 8
runs = []
for line in open("writeup/data/p58_candidates.txt"):
    m = re.match(r"run (\d+) start (\d+):(.*)\| next (\d+)", line)
    if not m: continue
    ps = list(map(int, m.group(3).split())); nxt = int(m.group(4))
    assert all(isprime(p) for p in ps) and all(nextprime(a) == b for a, b in zip(ps, ps[1:])) and nextprime(ps[-1]) == nxt
    good = [pari.isprime(4 * p * p + 1) == 1 for p in ps]
    # longest block of consecutive True
    best, cur, start = 0, 0, None
    for i, g in enumerate(good):
        cur = cur + 1 if g else 0
        if cur > best: best, start = cur, ps[i - cur + 1]
    if best >= K: runs.append((start, best))
runs.sort()
for start, k in runs[:10]: print(f"k = {k} run of consecutive primes p with 4p^2+1 prime starting at {start}")
print(f"{len(runs)} true runs with k >= {K}; earliest start:", runs[0][0] if runs else None)
