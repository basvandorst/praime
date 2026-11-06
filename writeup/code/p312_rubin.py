"""Puzzle 312 (Rubin): primes with every sum of three consecutive terms prime.

Q2: Rubin's rule (P1 = 3, P2 = 5, then the first prime not used yet with P(i-2) + P(i-1) + P(i)
prime) is run under both readings of 'not appeared yet' (in the term sequence only / in either
sequence) until a sum repeats a prime that already appeared.
Q1: the repaired greedy rule (the new term AND the new sum must both be unused primes) is run
for STEPS steps; the odd primes below the largest term that are still missing are counted.
Usage: p312_rubin.py [STEPS]
"""
import sys, bisect
from sympy import isprime, primerange, nextprime

def rubin(either):
    seq, sums = [3, 5], []
    S, T = {3, 5}, {}
    for i in range(2, 10**6):
        q = 7
        while q in S or (either and q in T) or not isprime(seq[-2] + seq[-1] + q):
            q = nextprime(q)
        s = seq[-2] + seq[-1] + q
        if s in S or s in T:
            where = f"term P({seq.index(s) + 1})" if s in S else f"sum S({T[s]})"
            return i + 1, q, s, where, seq, sums
        seq.append(q); sums.append(s); S.add(q); T[s] = i - 1

for either in (False, True):
    i, q, s, where, seq, sums = rubin(either)
    print(f"reading '{'either sequence' if either else 'term sequence only'}': at P({i}) = {q} the sum {s} repeats {where}")
    if either:
        print("  terms:", seq[:30]); print("  sums: ", sums[:30])

STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
LIM = 60 * STEPS
unused = list(primerange(7, LIM))          # sorted list of primes not yet used
used = {3, 5}
seq = [3, 5]
for i in range(STEPS):
    a = seq[-2] + seq[-1]
    for j, q in enumerate(unused):
        s = a + q
        if s not in used and isprime(s):
            break
    else:
        raise SystemExit("candidate list exhausted; raise LIM")
    del unused[j]
    k = bisect.bisect_left(unused, s)
    if k < len(unused) and unused[k] == s:
        del unused[k]
    used.add(q); used.add(s); seq.append(q)
frontier = max(seq)
assert frontier < LIM
missing = [p for p in unused if p < frontier]
print(f"repaired rule, {STEPS} steps: largest term {frontier}, {len(missing)} odd primes below it not yet placed, "
      f"smallest missing {missing[0] if missing else None}: every odd prime below {missing[0] if missing else frontier} appears exactly once")
