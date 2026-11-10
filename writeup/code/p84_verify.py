# Puzzle 84: sets of non-primes whose subset sums are all non-prime (TWA): the lexicographically earliest set
# (A052349), the odd variant (A128687) and the even-after-1 variant (A128688).  The OEIS terms (Noe 2007, Johnson
# 2010) extend the page's lists.
# (1) Every subset sum of each OEIS set is non-prime (sympy).
# (2) Minimality recomputed from scratch up to BOUND[kind]: with the prefix a(1..n) fixed, a(n+1) is the smallest
#     m > a(n) of the right parity with m + s non-prime for every subset sum s of the prefix.  (An m < a(n) that
#     fits the prefix would also fit every shorter prefix and would have been taken earlier, so m > a(n) loses
#     nothing.)  One pass over m with a segmented sieve; peak memory 570 MB, 400 s.
# (3) Nash's congruence construction, run greedily, gives a valid but much larger next term.
# Run from the repo root: (ulimit -v 1500000; .venv/bin/python3 writeup/code/p84_verify.py)
import sys
import numpy as np
from sympy import isprime, primerange
from sympy.ntheory.modular import crt

seqs = {
 "any":  ("A052349", [1, 8, 24, 25, 86, 1260, 1890, 14136, 197400, 10467660, 1231572090, 682616834970]),
 "odd":  ("A128687", [1, 9, 15, 39, 105, 731, 2805, 59535, 2658795, 78060135, 67554224745]),
 "even": ("A128688", [1, 8, 24, 86, 90, 780, 5940, 52350, 278460, 40768260, 6847205430, 5027286840810])}
page = {"any": 10, "odd": 8, "even": 7}                                    # terms printed on the page
BOUND = {"any": 1231572090, "odd": 78060135, "even": 6847205430}           # minimality is recomputed up to here
if len(sys.argv) > 1: BOUND = {k: min(v, int(float(sys.argv[1]))) for k, v in BOUND.items()}   # reduced run
SEG = 1 << 25                                                              # candidates per segment
DENSE = 6                                                                  # full-array passes before going sparse

def subset_sums(S):
    sums = {0}
    for x in S: sums |= {s + x for s in sums}
    return sorted(sums)

for kind, (A, S) in seqs.items():
    ss = subset_sums(S)
    assert all(not isprime(s) for s in ss), kind
    print(f"{kind} ({A}): {len(S)} terms, all {len(ss)} subset sums non-prime", flush=True)

BASE = np.array(list(primerange(2, 90000)), dtype=np.int64)                # 90000^2 > 6.9e9 + largest subset sum

def window(lo, n):
    """bool array w with w[i] = (lo + i is prime), 0 <= i < n"""
    w = np.ones(n, dtype=bool)
    for p in BASE:
        p = int(p)
        if p * p >= lo + n: break
        start = max(p * p, (lo + p - 1) // p * p) - lo
        w[start::p] = False
    if lo < 2: w[:2 - lo] = False
    return w

def survivors(lo, n, ss, kind):
    """all m in [lo, lo + n) of the right parity with m + s non-prime for every s in ss (ss sorted, ss[0] = 0)"""
    w = window(lo, n + ss[-1])
    cand = np.ones(n, dtype=bool)
    if kind != "any": cand[(lo + (kind == "even")) % 2::2] = False         # drop the wrong parity
    for s in ss[:DENSE]: cand &= ~w[s:s + n]
    idx = np.flatnonzero(cand)
    for s in ss[DENSE:]:
        if len(idx) == 0: break
        idx = idx[~w[idx + s]]
    del w, cand
    return [lo + int(i) for i in idx]

for kind, (A, S) in seqs.items():
    confirmed = 1
    for n in range(1, len(S)):
        if S[n] > BOUND[kind]: break
        ss = subset_sums(S[:n]); found = None; lo = S[n - 1] + 1
        while lo <= S[n] and found is None:
            got = survivors(lo, min(SEG, S[n] + 1 - lo), ss, kind)
            if got: found = got[0]
            lo += SEG
        assert found == S[n], (kind, n + 1, found, S[n])
        confirmed = n + 1
        print(f"  {kind}: a({n + 1}) = {S[n]} is the least admissible m > a({n})", flush=True)
    print(f"{kind} ({A}): a(1..{confirmed}) confirmed minimal; page had {page[kind]} terms", flush=True)

# (3) Nash's construction for the next odd term: m = 1 (mod 2) kills the odd sums; then each odd prime p takes the
# residue that kills the most remaining sums.  The CRT solution is valid; compare its size with the true terms.
def nash(prefix, parity):
    T = set(subset_sums(prefix)); mods = [2]; res = [parity]
    T = {s for s in T if (s + parity) % 2}                                 # m + s even (and > 2) is composite
    for p in primerange(3, 10**6):
        if not T: break
        best = max(range(p), key=lambda r: sum((s + r) % p == 0 for s in T))
        T = {s for s in T if (s + best) % p}
        mods.append(p); res.append(best)
    m, M = crt(mods, res); m = int(m); M = int(M)
    while m <= max(prefix) or any(isprime(m + s) for s in subset_sums(prefix)): m += M
    return m, len(mods)
for kind, parity, k in (("odd", 1, 8), ("even", 0, 7)):
    S = seqs[kind][1]
    m, npr = nash(S[:k], parity)
    assert all(not isprime(m + s) for s in subset_sums(S[:k]))
    print(f"Nash greedy after the page's {k} {kind} terms: {m} ({len(str(m))} digits, {npr} moduli) against the true a({k + 1}) = {S[k]}")
print("all checks passed")
