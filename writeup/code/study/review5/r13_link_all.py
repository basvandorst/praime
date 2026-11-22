"""What does the Hardy-Littlewood law predict for the spread of twin counts across the 15 admissible classes mod 210
over a whole window (the round-13 statistic Var/mean, ddof=1)?  Not computed in the study; round 31 says the measured
spreads (0.40 / 0.18 / 1.00 / 0.53) are 'consistent with a rigidity of about 0.5-0.6'.
Per class a: t restricted to m = a + 210 j has density r_a = 14 S2/lam^2 and pair density 14 S4(210 k)/lam^4, so
v_a = 1 - n r_a + 2 sum_k (n-k) 14 S4(210k) / (lam^4 n r_a); spread = (15 v_a - v(W)) / 14."""
import math, sys
import numpy as np
from lib import primes_upto

S2 = 1.3203236316937391
PR = primes_upto(20_000_000).astype(np.float64); P5 = PR[PR >= 5]
G4 = 8 * 27 / 16 * math.exp(float(np.sum(np.log1p(-4 / P5) - 4 * np.log1p(-1 / P5))))


def S4_mult210(n):
    """S4(210 k) for k = 0..n (exact closed form, k=0 set to 0)."""
    k = np.arange(n + 1, dtype=np.int64)
    S = np.full(n + 1, G4 * 3.0 * (5 / 3))                  # p = 5, 7 always divide d
    for p in primes_upto(n):
        p = int(p)
        if p < 11:
            continue
        S[p::p] *= (1 - 2 / p) / (1 - 4 / p)                # p | k
    f = lambda p: (1 - 3 / p) / (1 - 4 / p)
    for sgn in (-2, 2):                                     # p | 210 k + sgn
        cof = np.abs(210 * k + sgn)
        cof[0] = 1
        while True:
            ev = (cof % 2 == 0) & (cof > 0)
            if not ev.any():
                break
            cof[ev] //= 2
        for p in primes_upto(int(math.isqrt(210 * n + 2)) + 1):
            p = int(p)
            if p < 11:
                continue
            k0 = (-sgn * pow(210, -1, p)) % p
            S[k0::p] *= f(p)
            idx = np.arange(k0, n + 1, p)
            for _ in range(40):
                q = cof[idx] % p == 0
                if not q.any():
                    break
                cof[idx[q]] //= p
        big = cof > 1                                       # remaining cofactor is a prime > sqrt
        S[big] *= (1 - 3 / cof[big]) / (1 - 4 / cof[big])
    S[0] = 0.0
    return S


for name, L, W, vW in (('train', 1e7, 9e7, 0.33), ('A', 1e10, 5e7, None), ('B', 1e12, 5e7, None), ('C', 1e15, 2e7, None)):
    n = int(W // 210)
    S = S4_mult210(n)
    kk = np.arange(1, n, dtype=np.float64)
    cum = float(np.sum((n - kk) * S[1:n]))
    for lam in (math.log(L), math.log(L + W / 2), math.log(L + W)):
        ra = 14 * S2 / lam ** 2
        va = 1 - n * ra + 2 * 14 * cum / (lam ** 4 * n * ra)
        for v in ((0.3, 0.45) if vW is None else (vW - 0.05, vW + 0.05)):
            print(name, 'lam %.2f' % lam, 'v_a %.3f' % va, 'v(W)=%.2f' % v, 'predicted spread %.3f' % ((15 * va - v) / 14), flush=True)
