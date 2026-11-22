"""Own segmented sieve: per residue class mod 210, prime counts Y_a and the full gap table N[a, class]
(class k = gap 2k for k < 20, class 20 = gaps >= 40), consecutive pairs with both primes inside [L, L+W).
usage: count.py L W out.npz"""
import sys, math
import numpy as np

L, W, out = int(float(sys.argv[1])), int(float(sys.argv[2])), sys.argv[3]
Q = 210
U = [a for a in range(1, Q) if math.gcd(a, Q) == 1]
idx = np.full(Q, -1, np.int64); idx[U] = np.arange(48)
H = L + W
r = int(math.isqrt(H)) + 1
s = np.ones(r + 1, bool); s[:2] = False
for i in range(2, int(r ** 0.5) + 1):
    if s[i]:
        s[i * i::i] = False
small = np.nonzero(s)[0][1:]          # odd primes up to sqrt(H)
SEG = 50_000_000
N = np.zeros((48, 20), np.int64); last = None
for lo in range(L, H, SEG):
    hi = min(lo + SEG, H)
    base = lo | 1                     # first odd >= lo
    n = (hi - base + 1) // 2          # odd numbers base, base+2, ... < hi
    m = np.ones(n, bool)
    for p in small:
        p = int(p)
        if p * p >= hi:
            break
        st = max(p * p, (base + p - 1) // p * p)
        if st % 2 == 0:
            st += p
        m[(st - base) // 2::p] = False
    P = base + 2 * np.nonzero(m)[0].astype(np.int64)
    if lo <= 2 < hi:
        P = np.concatenate([[2], P])
    if last is not None:
        P = np.concatenate([[last], P])
    last = int(P[-1])
    a = idx[P[:-1] % Q]; g = np.diff(P)
    assert (a >= 0).all()
    N += np.bincount(a * 20 + np.minimum(g // 2, 20) - 1, minlength=960).reshape(48, 20)
    del m, P, a, g
Y = N.sum(1)
np.savez(out, N=N, Y=Y, U=np.array(U), L=L, W=W)
print('primes counted (pairs)', int(Y.sum()), 'twins', int(N[:, 0].sum()), 'gap6', int(N[:, 2].sum()), 'ge40', int(N[:, 19].sum()))
