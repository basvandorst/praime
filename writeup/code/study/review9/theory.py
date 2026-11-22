"""Independent HL prediction for the twin counts per admissible class mod Q.
F(c) = sum_{k=-(n-1)}^{n-1} w(|k|) (rho(c + Qk) - 1), d = c + Qk != 0, w(k) = sum_j r_j r_{j+k},
rho(d) = S({0,2,d,d+2})/S({0,2})^2 restricted to p not dividing Q: prod_p (1 - nu_p/p)/(1 - 2/p)^2,
nu_p = #{0,2,d,d+2} mod p.  Var X_a = sum r(1-r) + F(0); Cov(X_a, X_b) = F(b - a)."""
import math
import numpy as np

PMAX = 10 ** 7


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


PR = primes_upto(PMAX)


def gen(p):
    return np.log1p(-4.0 / p) - 2 * np.log1p(-2.0 / p)


def setup(Q, L0, W):
    qf = [int(p) for p in PR[:20] if Q % p == 0]
    P = PR[(Q % PR) != 0]
    base = float(np.sum(gen(P.astype(float)))) - 4.0 / (PMAX * math.log(PMAX))   # tail p > PMAX
    c = math.prod(1 / (1 - 1 / p) ** 2 for p in qf) * math.exp(
        float(np.sum(np.log1p(-2.0 / P[P > 2]) - 2 * np.log1p(-1.0 / P[P > 2]))) )
    n = W // Q
    m = L0 + Q * np.arange(n) + Q / 2.0
    r = c / (np.log(m) * np.log(m + 2))
    N = 1 << int(math.ceil(math.log2(2 * n)))
    f = np.fft.rfft(r, N)
    w = np.fft.irfft(f * np.conj(f), N)[:n]
    return dict(Q=Q, qf=qf, P=P, base=base, c=c, n=n, r=r, w=w)


def F(S, c):
    Q, n, P, base = S['Q'], S['n'], S['P'], S['base']
    k = np.arange(-(n - 1), n, dtype=np.int64)
    d = c + Q * k
    wt = S['w'][np.abs(k)].copy()
    if c == 0:
        wt[n - 1] = 0.0; d[n - 1] = 1                         # drop d = 0
    L = np.full(d.size, base)
    lim = int(math.isqrt(int(np.abs(d).max()) + 2)) + 1
    for sh, nu in ((0, 2), (-2, 3), (2, 3)):                 # p | d (nu 2), p | d-2, p | d+2 (nu 3)
        cof = np.abs(d + sh)
        cof[cof == 0] = 1
        for q in [2] + [int(x) for x in S['qf']]:
            while True:
                m = (cof % q == 0) & (cof > 1)
                if not m.any(): break
                cof[m] //= q
        for p in P[P <= lim]:
            p = int(p)
            i0 = (-(c + sh) * pow(Q, -1, p) + (n - 1)) % p    # index i = k + n - 1 with p | d + sh
            idx = np.arange(i0, d.size, p)
            idx = idx[cof[idx] % p == 0]
            if idx.size == 0: continue
            L[idx] += math.log((1 - nu / p) / (1 - 2 / p) ** 2) - gen(p)
            while idx.size:
                cof[idx] //= p
                idx = idx[cof[idx] % p == 0]
        big = cof > 1
        pb = cof[big].astype(float)
        L[big] += np.log1p(-nu / pb) - 2 * np.log1p(-2 / pb) - gen(pb)
    rho = np.exp(L)
    return float(np.sum(wt * (rho - 1.0))), rho, d
