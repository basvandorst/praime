# Puzzle 58 / 54: runs of k consecutive primes p with 4p^2 + 1 prime.
# Mode 1 (default): check data/p58_runs.txt: every listed run is a maximal run (members consecutive primes with
#   4p^2 + 1 prime, proven with PARI; the primes before and after fail).
# Mode 2: p58_verify.py K candidate_file: read the candidate runs printed by p58_search.c (consecutive primes passing
#   the small-prime filter on 4p^2 + 1), test 4p^2 + 1 exactly for each member, report the true runs with >= K members.
import re, sys
from sympy import isprime, nextprime, prevprime
import cypari2
pari = cypari2.Pari()
def good(p): return pari.isprime(4 * p * p + 1) == 1
if len(sys.argv) > 2:
    K = int(sys.argv[1]); runs = []
    for line in open(sys.argv[2]):
        m = re.match(r"run (\d+) start (\d+):(.*)\| next (\d+)", line)
        if not m: continue
        ps = list(map(int, m.group(3).split())); nxt = int(m.group(4))
        assert all(isprime(p) for p in ps) and all(nextprime(a) == b for a, b in zip(ps, ps[1:])) and nextprime(ps[-1]) == nxt
        g = [good(p) for p in ps]; best, cur, start = 0, 0, None
        for i, x in enumerate(g):
            cur = cur + 1 if x else 0
            if cur > best: best, start = cur, ps[i - cur + 1]
        if best >= K: runs.append((start, best))
    runs.sort()
    for start, k in runs[:10]: print(f"k = {k} run of consecutive primes p with 4p^2+1 prime starting at {start}")
    print(f"{len(runs)} true runs with k >= {K}; earliest start:", runs[0][0] if runs else None)
else:
    for line in open("writeup/data/p58_runs.txt"):
        m = re.match(r"k=(\d+) start=(\d+) members=(.*) next_prime=(\d+)", line)
        if not m: continue
        k, s, ps, nxt = int(m.group(1)), int(m.group(2)), list(map(int, m.group(3).split())), int(m.group(4))
        assert len(ps) == k and ps[0] == s and all(isprime(p) and good(p) for p in ps)
        assert all(nextprime(a) == b for a, b in zip(ps, ps[1:])) and nextprime(ps[-1]) == nxt
        assert not good(prevprime(s)) and not good(nxt)
        print(f"k = {k} maximal run from {s} to {ps[-1]}: members consecutive primes, 4p^2+1 prime (PARI), neighbours fail")
    print("all checks passed")
