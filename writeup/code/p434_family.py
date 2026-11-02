"""Puzzle 434: the one-unknown family x = Q*p, p = c*r - 1 (p, r prime). Then y = sigma(Q)*c*r and
F = L*r + sigma(M) + Q with M = sigma(Q)*c, L = sigma(M) - M - Q*c, so F = 0 forces
r = (sigma(M) + Q) / (M + Q*c - sigma(M)). Usage: p434_family.py N QMAX (M <= N, Q <= QMAX)."""
import numpy as np, sys, math, gmpy2
from sympy import divisor_sigma
N = int(sys.argv[1]); Qmax = int(sys.argv[2])
S = math.isqrt(N)
s = np.zeros(N + 1, dtype=np.uint32)
for d in range(1, S + 1):
    s[d::d] += np.uint32(d)
for k in range(1, S + 1):
    hi = N // k
    if hi <= S: break
    s[k * (S + 1)::k] += np.arange(S + 1, hi + 1, dtype=np.uint32)
cand = tested = 0
hits = []
for Q in range(1, Qmax + 1):
    sQ = int(s[Q])
    cmax = N // sQ
    c = np.arange(2, cmax + 1, 2, dtype=np.int64)
    M = sQ * c
    sM = s[M].astype(np.int64)
    den = M + Q * c - sM
    num = sM + Q
    ok = (den > 0) & (num % np.maximum(den, 1) == 0)
    for cc, n_, d_ in zip(c[ok].tolist(), num[ok].tolist(), den[ok].tolist()):
        r = n_ // d_
        cand += 1
        if r < 3 or not gmpy2.is_prime(r):
            continue
        p = cc * r - 1
        if p < 3 or not gmpy2.is_prime(p) or math.gcd(Q, p) != 1:
            continue
        x = Q * p
        tested += 1
        y = int(divisor_sigma(x))
        if int(divisor_sigma(y)) == y + x and Q > 1:
            hits.append((x, Q, cc, r))
            print('HIT x=%d = %d * %d, p+1 = %d * %d' % (x, Q, p, cc, r), flush=True)
print('Q <=', Qmax, 'candidates', cand, 'prime pairs tested', tested, 'hits', len(hits))
