# Puzzle 1104: validation of every claim in writeup/1104.md.
# 1. The theorem S_n = mu(d) phi(p-1)/phi(d) (mod p), d = (p-1)/gcd(n, p-1), checked for all p < 2000 and n <= 60.
# 2. The parity rule S_n = (number of odd primitive roots) (mod 2), same range.
# 3. The table data/p1104_table.txt (least p with S_n prime, n = 1..60): primitive roots found by the direct test, S_n recomputed,
#    the listed S_n proved prime (PARI isprime, APR-CL), every smaller p shown composite (sympy.isprime) or excluded by a rule.
# 4. Q1: p1104_q1.c lists p <= 10^6 with p-1 squarefree and S_1, S_2 prime; S_3..S_5 tested here. Only p = 43 has five prime sums.
# 5. The heuristic expectations, calibrated on the S_1 and S_1 & S_2 counts of step 4.
# Run from the repo root (about 5 minutes): .venv/bin/python3 writeup/code/p1104_verify.py
import os, subprocess
from math import gcd, log, prod
from sympy import primerange, factorint, isprime, totient
from sympy.functions.combinatorial.numbers import mobius
import cypari2
pari = cypari2.Pari(); pari.allocatemem(256*2**20)
here = os.path.dirname(__file__)

def proots(p):
    m = p - 1; qs = list(factorint(m))
    return [g for g in range(1, p) if all(pow(g, m//q, p) != 1 for q in qs)]
def squarefree(m): return all(e == 1 for e in factorint(m).values())

# 1 and 2
for p in primerange(3, 2000):
    R = proots(p); m = p - 1; odd = sum(r & 1 for r in R) % 2
    for n in range(1, 61):
        S = sum(pow(r, n, p) for r in R) % p
        d = m//gcd(n, m)
        expect = (int(mobius(d))*int(totient(m))//int(totient(d))) % p
        assert S == expect, (p, n)
        assert (S == 0) == (not squarefree(d))
        assert sum(pow(r, n) for r in R) % 2 == odd
print("theorem and parity rule checked for p < 2000, n <= 60")

# 3
table = {}
for line in open(os.path.join(here, "..", "data", "p1104_table.txt")):
    n, p, S = line.split(); table[int(n)] = (int(p), int(S))
assert sorted(table) == list(range(1, 61))
pmax = max(p for p, S in table.values())
roots = {p: proots(p) for p in primerange(3, pmax + 1)}
for n in range(1, 61):
    pn, Sn = table[n]
    for p in primerange(3, pn + 1):
        R = roots[p]; m = p - 1; d = m//gcd(n, m)
        if p < pn:
            if not squarefree(d): continue                       # p | S_n, proved
            if p > 3 and sum(r & 1 for r in R) % 2 == 0: continue   # S_n even and > 2, proved
            assert not isprime(sum(pow(r, n) for r in R)), (n, p)
        else:
            S = sum(pow(r, n) for r in R); assert S == Sn
            assert pari.isprime(pari(S)) == 1, (n, p)
    print("n = %2d: least p = %6d, S_n has %3d digits, proved prime; all smaller p composite or excluded" % (n, pn, len(str(Sn))), flush=True)
gusev = {1: 3, 2: 5, 3: 23, 4: 5, 5: 11, 6: 67, 7: 103, 8: 463, 9: 719, 10: 1499, 11: 1907, 12: 2131, 13: 1087, 14: 43, 15: 647, 16: 31,
         17: 3343, 18: 4603, 19: 11959, 20: 3343, 21: 6863, 22: 179, 23: 2347, 24: 4603}
assert all(table[n][0] == gusev[n] for n in gusev), "Gusev's table reproduced"

# 4
subprocess.run(["gcc", "-O2", "-o", "/tmp/p1104_q1", os.path.join(here, "p1104_q1.c")], check=True)
out = subprocess.run(["/tmp/p1104_q1", "1000000"], capture_output=True, text=True, check=True)
cands = [tuple(map(int, l.split())) for l in out.stdout.split("\n") if l.strip()]
stats = out.stderr.strip(); print(stats)
assert "p-1 squarefree 29396, S1 prime 1291, S1 and S2 prime 69" in stats
runs = {}
for p, S1, S2 in cands:
    R = proots(p); S = [sum(pow(r, n) for r in R) for n in range(1, 8)]
    assert S[0] == S1 and S[1] == S2
    k = 0
    while k < 7 and isprime(S[k]): k += 1
    runs[p] = k
print("p <= 10^6 with S_1, S_2 prime:", len(cands), "; with S_1..S_3 prime:", sorted(p for p, k in runs.items() if k >= 3))
assert [p for p, k in runs.items() if k >= 4] == [23, 43] and [p for p, k in runs.items() if k >= 5] == [43] and runs[43] == 5

# 5 heuristic. S_n ~ phi(p-1) p^n/(n+1). Conditions shared by all n: p-1 squarefree (else p | S_1), the number of odd roots odd (probability 1/2).
# Given that, S_n is prime with probability about 2/ln S_n. Joint corrections: S_n mod 3 depends only on n mod 2, S_n mod 5 only on n mod 4.
def expected(kmax, X, lo=47):
    E = 0.0
    for p in primerange(lo, X):
        m = p - 1
        if not squarefree(m): continue
        ph = int(totient(m)); pr = 1.0
        for n in range(1, kmax + 1): pr *= min(1.0, 2/log(ph*p**n/(n + 1)))
        odd_n = len([n for n in range(1, kmax + 1) if n % 2]); even_n = kmax - odd_n
        boost3 = (3/2)**max(0, odd_n - 1)*(3/2)**max(0, even_n - 1)            # mod 3: one condition per parity class instead of one per n
        classes5 = {n % 4 for n in range(1, kmax + 1)}; boost5 = (5/4)**(kmax - len(classes5))
        E += 0.5*pr*boost3*boost5
    return E
print("calibration, p in [47, 10^6): expected S_1 prime %.0f (found 1291), S_1 and S_2 prime %.1f (found 69)" % (expected(1, 10**6), expected(2, 10**6)))
for k in (5, 6):
    print("k = %d: expected number of p in [47, 10^6) with S_1..S_k prime: %.3f" % (k, expected(k, 10**6)))
# tail: the sum behaves like c_k X/ln^(k+1) X; integrate the per-prime term with density 1/ln x beyond 10^6
import numpy as np
for k in (5, 6):
    odd_n = (k + 1)//2; even_n = k - odd_n
    boost = (3/2)**(odd_n - 1)*(3/2)**(even_n - 1)*(5/4)**(k - len({n % 4 for n in range(1, k + 1)}))
    const = 0.5*0.3745*boost*prod(2/(n + 1) for n in range(1, k + 1))      # 0.3745 = density of primes with p-1 squarefree (29396/78497)
    for X in (1e8, 1e9, 1e10, 1e11, 1e12):
        xs = np.logspace(6, log(X, 10), 2000); f = const/np.log(xs)**(k + 1)
        E = np.trapezoid(f, xs)
        print("k = %d: expected further solutions with 10^6 < p < %.0e: %.2f" % (k, X, E))
print("all checks passed")
