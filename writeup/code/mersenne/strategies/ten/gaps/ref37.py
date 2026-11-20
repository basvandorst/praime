"""Independent reference for prime37.html: gaps between consecutive primes up to 1000 and 10000 (sympy)."""
import json
import sys
from collections import Counter
from sympy import primerange

out = {}
for X in (1000, 10000):
    P = list(primerange(2, X + 1))
    G = [b - a for a, b in zip(P, P[1:])]
    rec, best = [], 0
    for i, g in enumerate(G):
        if g > best:
            best = g
            rec.append([i + 1, P[i], P[i + 1], g])
    cnt = Counter(G)
    out[X] = {"primes": len(P), "gaps": len(G), "sum": sum(G), "max": max(G),
              "count": {str(k): v for k, v in sorted(cnt.items())}, "records": rec, "list": G}
    print("up to %d: %d primes, %d gaps, sum %d (= %d - 2), mean %.3f, largest %d" % (
        X, len(P), len(G), sum(G), P[-1], sum(G) / len(G), max(G)))
    print("  how often each gap occurs:", dict(sorted(cnt.items())))
    print("  most common gap:", max(cnt, key=cnt.get), "; record gaps (n, p_n, p_n+1, gap):", rec)
json.dump(out, open("ref37.json", "w"))
X = 10000
P = list(primerange(2, X + 1))
for i in range(len(P) - 1):
    print("prime %d-%d: %d -> %d, gap %d" % (i + 1, i + 2, P[i], P[i + 1], P[i + 1] - P[i]))
