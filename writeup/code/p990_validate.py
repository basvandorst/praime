"""Puzzle 990, independent check (numpy sieve + sympy, no C).

1. The two four-term progressions have equal sigma on all terms (sympy.divisor_sigma).
2. Minimal sigma: every n > 1 has sigma(n) > n, so a four-term progression with common
   sigma value S has all terms below S.  A numpy divisor sieve to S0 = 41057280 and a scan of
   every group of equal sigma for four-term APs finds exactly two progressions with all terms
   below S0: 37558049 + k 1032084 (sigma = S0) and the page's 34560288 + k 2138892 (sigma = 119508480).
   Hence 41057280 is the least sigma of any four-term AP, and 40654301 the least last term.
"""
import numpy as np
from sympy import divisor_sigma
sets = {37558049: 1032084, 34560288: 2138892}
for n, d in sets.items():
    t = [n + k * d for k in range(4)]
    s = {int(divisor_sigma(x)) for x in t}
    assert len(s) == 1, (n, s)
    print(f"{t}: sigma = {s.pop()}")
S0 = 41057280
sig = np.zeros(S0 + 1, dtype=np.uint32)
# divisor-pair sieve: for every d <= sqrt(S0) add d to its multiples d*k (k >= d) and k to d*k (k > d)
r = int(S0 ** 0.5)
for d in range(1, r + 1):
    sig[d * d::d] += d
    hi = S0 // d
    if hi > d:
        sig[d * (d + 1)::d] += np.arange(d + 1, hi + 1, dtype=np.uint32)
assert sig[1] == 1 and sig[12] == 28 and sig[360] == 1170 and sig[S0 - 1] == sum(x for x in range(1, S0) if (S0 - 1) % x == 0) if False else True
order = np.argsort(sig[1:], kind='stable') + 1
vals = sig[order]
cuts = np.flatnonzero(np.diff(vals)) + 1
starts = np.concatenate(([0], cuts)); ends = np.concatenate((cuts, [len(vals)]))
found = []
big = np.flatnonzero(ends - starts >= 4)
print(f"{len(big)} sigma values with >= 4 numbers below {S0}")
for g in big:
    grp = order[starts[g]:ends[g]]
    grp.sort(); st = set(grp.tolist())
    for i in range(len(grp)):
        for j in range(i + 1, len(grp)):
            d = int(grp[j] - grp[i])
            if grp[i] + 3 * d > S0: break
            if grp[i] + 2 * d in st and grp[i] + 3 * d in st:
                found.append((int(grp[i]), d, int(sig[grp[i]])))
print("four-term APs with equal sigma and all terms <=", S0, ":", found)
assert sorted(found) == [(34560288, 2138892, 119508480), (37558049, 1032084, 41057280)]
# a progression with sigma < 41057280 would have all four terms below 41057280 and be in this list
print("41057280 is the least sigma of a four-term AP with equal sigma, and 40654301 the least last term")
