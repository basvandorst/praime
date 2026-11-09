# Puzzle 1021, independent check of data/p1021_table.txt: every row (n, type, p, x, factors): p prime, the next
# prime(s) follow, the sum equals x, x is squarefree with exactly n prime factors (sympy, plain integers; primality
# is BPSW above 2^64 here, the proofs are in p1021_verify.py via PARI).  Then minimality for n <= 7 by brute force
# over consecutive primes (omega from a smallest-prime-factor sieve).  Run from the repo root.
import numpy as np
from sympy import isprime, nextprime, factorint
rows = [l.split() for l in open("writeup/data/p1021_table.txt") if not l.startswith("#")]
tab = {}
for n, typ, p, x, *fs in rows:
    n, p, x = int(n), int(p), int(x); fs = list(map(int, fs))
    q = nextprime(p); s = p + q + 1 if typ == 'A' else p + q + nextprime(q)
    assert isprime(p) and s == x, (n, typ)
    f = factorint(x); assert sorted(f) == fs and all(e == 1 for e in f.values()) and len(f) == n, (n, typ)
    tab[(n, typ)] = p
print(f"{len(rows)} rows checked (n up to {max(n for n, _ in tab)}): primes, consecutiveness, sums and factorisations agree")
# brute force n <= 7: least p over consecutive primes below 1e7
N = 10**7 + 100; spf = np.zeros(3 * N, np.int32)
for i in range(2, int((3 * N) ** 0.5) + 1):
    if spf[i] == 0: spf[i*i::i][spf[i*i::i] == 0] = i
def omega_sqfree(x):
    k = 0
    while x > 1:
        p = spf[x] or x; x //= p
        if x % p == 0: return -1
        k += 1
    return k
s = np.ones(N, bool); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]: s[i*i::i] = False
P = np.flatnonzero(s)
least = {}
for i in range(len(P) - 2):
    for typ, x in (('A', int(P[i] + P[i+1] + 1)), ('B', int(P[i] + P[i+1] + P[i+2]))):
        k = omega_sqfree(x)
        if k > 0 and (k, typ) not in least: least[(k, typ)] = int(P[i])
for n in range(1, 8):
    for typ in 'AB': assert least[(n, typ)] == tab[(n, typ)], (n, typ, least[(n, typ)], tab[(n, typ)])
print("brute force over consecutive primes below 1e7 confirms the least p for n <= 7, both types")
print("all checks passed")
