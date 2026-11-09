"""Puzzles 1034 / 1091, Python version of p1034_search.c (same algorithm, for smaller ranges).
A chain p, p + 6k, p + 6k + 6(k+1), ... of primes is a run of consecutive n with f(n) = 3n^2 + 3n + c prime.
f has no root mod q (q >= 5) iff D = 9 - 12c is a non-residue mod q; the Hardy-Littlewood constant
C(c) = 3 prod_{q>=5} (1 - w_q/q)/(1 - 1/q), w_q = 1 + (D|q), measures the prime density of f.
For c in [-CMAX, CMAX], odd and prime to 3, with (D|q) = -1 for q = 5..23 and C(c) >= CMIN, scan n = 0..N and
report runs of >= RMIN consecutive prime values.      Usage: p1034_search.py CMAX N CMIN RMIN"""
import sys
from math import prod
from gmpy2 import mpz, is_prime, legendre
from sympy import primerange
CMAX, N, CMIN, RMIN = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
qs = list(primerange(5, 2000)); q7 = qs[:7]
ncand = nscan = 0; best = (0, None)
for c in range(-CMAX, CMAX + 1):
    if c % 2 == 0 or c % 3 == 0: continue
    D = 9 - 12 * c
    if any(legendre(D % q, q) != -1 for q in q7): continue
    ncand += 1
    C = 3.0
    for q in qs:
        w = 1 + int(legendre(D % q, q)) if D % q else 1
        C *= (1 - w / q) / (1 - 1 / q)
    if C < CMIN: continue
    nscan += 1
    run = 0
    for n in range(N + 1):
        v = 3 * n * n + 3 * n + c
        pr = v > 1 and is_prime(v)
        if pr: run += 1
        if not pr or n == N:
            if run >= RMIN:
                n0 = (n if pr else n - 1) - run + 1
                print(f"run {run}: c={c} C={C:.2f} n0={n0} p={3*n0*n0+3*n0+c} step={6*(n0+1)}", flush=True)
                if run > best[0]: best = (run, c)
            run = 0
print(f"c in [-{CMAX}, {CMAX}]: {ncand} pass the residue filter, {nscan} scanned (C >= {CMIN}), longest run {best}")
