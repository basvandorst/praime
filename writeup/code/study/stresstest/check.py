"""Independent recount: gmpy2.next_prime chains from 1e9 and 1e12 (first 200k primes), compare (g1,g2) counts with the sieve-based H."""
import sys, gmpy2, numpy as np
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import common as C
for L0 in (10**9, 10**12):
    ps = []; p = gmpy2.next_prime(L0 - 1)
    while len(ps) < 200002:
        ps.append(int(p)); p = gmpy2.next_prime(p)
    P = np.array(ps, dtype=np.int64)
    S = np.nonzero(C.primes_mask(L0, int(P[-1] - L0 + 1)))[0].astype(np.int64) + L0
    assert np.array_equal(S, P), L0
    g1 = P[1:-1] - P[:-2]; g2 = P[2:] - P[1:-1]
    print(L0, 'sieve == next_prime chain for', len(P), 'primes; cell(6,6)', int(np.sum((g1 == 6) & (g2 == 6))), 'max gap', int(max(g1.max(), g2.max())))
