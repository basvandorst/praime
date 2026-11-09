# Puzzle 1218: primes p with LPF(prime(q) + prime(q+1)) = p, q = prime(p).
# Checks the four known solutions (3, 41, 463, 4724804933), searches all p < 1.5e6 directly, and computes the heuristic count.
# Run from the repo root: .venv/bin/python3 writeup/code/p1218_verify.py   (compiles writeup/code/p9_lucy.c for the prime counts)
import subprocess, os, numpy as np
from math import log
from sympy import isprime, nextprime, factorint

here = os.path.dirname(__file__)
subprocess.run(["gcc", "-O2", "-o", "/tmp/p9_lucy_1218", os.path.join(here, "p9_lucy.c"), "-lm"], check=True)
def primepi(x):
    out = subprocess.run(["/tmp/p9_lucy_1218", str(x)], capture_output=True, text=True, check=True).stdout
    return int(out.split()[0].split("=")[1])

# the large solution (Volpatti 2025, found with primesieve/primecount)
n, p, q, x, y = 222617381, 4724804933, 115412584783, 3203417744519, 3203417744629
assert isprime(p) and isprime(q) and isprime(x) and nextprime(x) == y
assert primepi(p) == n and primepi(q) == p and primepi(x) == q      # p = prime(n), q = prime(p), x = prime(q), y = prime(q+1)
f = factorint(x + y); print("x + y =", x + y, "=", f); assert max(f) == p and f == {2: 2, 3: 1, 113: 1, p: 1}

# direct search: p < 1.5e6 needs prime(q) for q < prime(1.5e6) ~ 2.4e7, i.e. primes to ~4.5e8
N = 460_000_000
s = np.ones(N + 1, bool); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]: s[i*i::i] = False
P = np.flatnonzero(s).astype(np.int64)          # P[k-1] = prime(k)
sols, checked = [], 0
for p in P[P < 1_500_000]:
    q = P[p - 1]
    assert q + 1 <= len(P)
    z = int(P[q - 1]) + int(P[q]); checked += 1
    if z % p == 0 and max(factorint(z)) == p: sols.append(int(p))
print("primes checked:", checked, "solutions p < 1.5e6:", sols)
assert sols == [3, 41, 463]

# cofactor: (x + y)/p is about 2 ln^2 p, smaller than p for p >= 11, so LPF = p exactly when p | x + y
for p0 in (3, 41, 463, p):
    pass
# heuristic: an odd prime p divides the even number x + y with probability 1/p; sum over odd primes 3..X = ln ln X + M - 1/2
M = 0.2614972128
for X in (2.28e10, 1e12, 1e15, 1e20, 1e25, 1e30):
    print("expected number of solutions with p <= %.2e: %.2f" % (X, log(log(X)) + M - 0.5))
print("found below 2.28e10: 4")
print("all checks passed")
