"""Puzzle 309: one-unknown family x = Q*p, p - 1 = c1*r1, p + 1 = c2*r2 (p, r1, r2 prime):
c1*c2*(f(x) - x) = N - D*p with D = c1*c2*Q + c1*phi(M2) - c2*sigma(M1), M1 = phi(Q)*c1, M2 = sigma(Q)*c2,
N = c2*sigma(M1)*(c1 - 1) + c1*phi(M2)*(c2 - 1), so f(x) = x forces p = N/D. Usage: p309_family.py QMAX CMAX."""
import sys, math
import numpy as np, gmpy2
from sympy import totient, divisor_sigma
Qmax, Cmax = int(sys.argv[1]), int(sys.argv[2])
L = 5 * Qmax * Cmax
sig = np.zeros(L + 1, dtype=np.int64)
S = math.isqrt(L)
for d in range(1, S + 1):
    sig[d::d] += d
for k in range(1, S + 1):
    hi = L // k
    if hi <= S: break
    sig[k * (S + 1)::k] += np.arange(S + 1, hi + 1, dtype=np.int64)
ph = np.arange(L + 1, dtype=np.int64)
sv = np.ones(L + 1, dtype=np.bool_); sv[:2] = False
for i in range(2, S + 1):
    if sv[i]: sv[i*i::i] = False
for p in np.nonzero(sv)[0]:
    p = int(p); v = ph[p::p]; v -= v // p
cands = tested = 0
best = []
c2 = np.arange(2, Cmax + 1, 2, dtype=np.int64)
for Q in range(1, Qmax + 1):
    fQ, sQ = int(ph[Q]), int(sig[Q])
    M2 = sQ * c2
    keep = M2 <= L
    c2k, M2k = c2[keep], M2[keep]
    phM2 = ph[M2k]
    for c1 in range(2, Cmax + 1, 2):
        M1 = fQ * c1
        if M1 > L: break
        sM1 = int(sig[M1])
        g = np.gcd(c2k, c1)
        D = c1 * c2k * Q + c1 * phM2 - c2k * sM1
        Nn = c2k * sM1 * (c1 - 1) + c1 * phM2 * (c2k - 1)
        ok = (g == 2) & ((c1 * c2k) % 3 == 0) & (D > 0)
        ok &= (Nn % np.where(D > 0, D, 1) == 0)
        for cc2, n_, d_ in zip(c2k[ok].tolist(), Nn[ok].tolist(), D[ok].tolist()):
            p = n_ // d_
            cands += 1
            if (p - 1) % c1 or (p + 1) % cc2:
                continue
            r1, r2 = (p - 1) // c1, (p + 1) // cc2
            if not (gmpy2.is_prime(p) and gmpy2.is_prime(r1) and gmpy2.is_prime(r2)):
                continue
            x = Q * p
            tested += 1
            f = int(divisor_sigma(totient(x))) - int(totient(divisor_sigma(x)))
            print(('HIT ' if f == x else 'prime triple, f != x ') + 'x=%d Q=%d c1=%d c2=%d p=%d' % (x, Q, c1, cc2, p), flush=True)
print('Q <= %d, c <= %d: integer candidates %d, prime triples %d' % (Qmax, Cmax, cands, tested))
