"""Puzzle 326: terms n with n = pi(product of the factorials of the digits of n): 0, 3, 34, 52, 2800414, ... .
For every digit count k <= KMAX (argument): every multiset of k digits whose factorial product P lies in the window where
pi(P) has k digits is tested; the leading digits of pi(P) are first estimated by li(P) with a generous error (2 sqrt(P) ln P)
to skip multisets that cannot match, and pi(P) is computed exactly (sympy.primepi) for the rest.  Prints every term found.
"""
import sys, math, time
from itertools import combinations_with_replacement
from sympy import primepi
from math import factorial
from mpmath import li
KMAX = int(sys.argv[1])
F = [factorial(d) for d in range(10)]
def bounds(N):
    l = math.log(N); ll = math.log(l)
    return N * (l + ll - 1), N * (l + ll)
found = []
for k in range(1, KMAX + 1):
    lo = bounds(max(10**(k-1), 6))[0] if k > 1 else 0
    hi = bounds(10**k)[1]
    t0 = time.time()
    byP = {}
    for ms in combinations_with_replacement(range(10), k):
        P = 1
        for d in ms: P *= F[d]
        if lo <= P <= hi:
            byP.setdefault(P, []).append(ms)
    tested = 0
    for P, mss in sorted(byP.items()):
        # leading digits of pi(P) from li(P): |pi(P) - li(P)| < sqrt(P) ln P for P <= 1e13 (checked numerically), so the
        # first ~4 digits of pi(P) are those of li(P) unless a carry; keep multisets consistent with either neighbour
        if P < 100:
            est, err = int(primepi(P)), 0
        else:
            est = int(li(P)); err = int(2 * math.sqrt(P) * math.log(P)) + 1
        a, b = str(max(est - err, 0)), str(est + err)
        if len(a) == len(b) == k:
            pref = 0
            while pref < k and a[pref] == b[pref]: pref += 1
            lead = a[:pref]
        elif len(a) <= k <= len(b):
            lead = ''
        else:
            continue
        ok = []
        for ms in mss:
            cnt = {d: ms.count(d) for d in range(10)}
            if all(cnt[int(c)] >= lead.count(c) for c in set(lead)):
                ok.append(ms)
        if not ok: continue
        tested += 1
        n = int(primepi(P)); sn = str(n)
        for ms in ok:
            if len(sn) == k and sorted(sn) == sorted(str(d) for d in ms):
                found.append((n, ms, P)); print("FOUND", n, ms, P, flush=True)
    print(f"k={k}: {len(byP)} distinct products in window, {tested} needed pi(P), {time.time()-t0:.1f} s", flush=True)
print("all found:", found)
