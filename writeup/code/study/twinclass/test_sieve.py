"""Check tcsieve against sympy on two small ranges (all statistics, all classes)."""
import sys; sys.path.insert(0, "/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/twinclass")
import numpy as np, sympy, math
from tcsieve import count_window
for L0, W in ((1000003, 200000), (10**7 - 5, 123457)):
    P, WP, X = count_window(L0, W, Q=2310, nsub=3, seg_odd=1 << 14, verbose=False)
    edges = [L0 + (W * s) // 3 for s in range(4)]
    P2 = np.zeros_like(P); WP2 = np.zeros_like(WP); X2 = {h: np.zeros_like(P) for h in X}
    for p in sympy.primerange(L0, L0 + W):
        s = sum(p >= e for e in edges[1:]); c = p % 2310
        P2[s, c] += 1; WP2[s, c] += 1 / math.log(p)
        for h in X:
            if sympy.isprime(p + h): X2[h][s, c] += 1
    print(L0, W, 'P ok', (P == P2).all(), 'WP maxdiff', abs(WP - WP2).max(), 'X ok', all((X[h] == X2[h]).all() for h in X), 'primes', P.sum(), 'twins', X[2].sum())
