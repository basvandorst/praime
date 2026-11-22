"""Section 10 check: kappa3/mean of prime counts in aligned intervals [n, n+H), n = 0 mod H, from HL sums.
Own singular series: S0(X) = 2^k (5/4)^k prod_{p != 2,5} (1 - nu_p/p) / (1 - 1/p)^k."""
import math, sys, json
from itertools import combinations
import numpy as np
from sieve import small_primes, mask

PR = small_primes(10 ** 7)


def tail(k, H):
    p = PR[(PR >= H)].astype(float)
    p = p[(p != 2) & (p != 5)]
    return math.exp(np.sum(np.log1p(-k / p) - k * np.log1p(-1 / p)))


def S0(sets, H):
    sets = np.asarray(sets, dtype=np.int64); n, k = sets.shape
    logS = np.full(n, k * math.log(2 * 1.25))
    for p in PR[PR < H]:
        p = int(p)
        if p in (2, 5):
            continue
        r = np.sort(sets % p, axis=1)
        nu = 1 + np.count_nonzero(np.diff(r, axis=1), axis=1)
        f = 1 - nu / p
        logS += np.log(np.maximum(f, 1e-300)) - k * math.log(1 - 1 / p)
        logS[f <= 0] = -np.inf
    return np.exp(logS) * tail(k, H)


def sums(H):
    U = np.array([a for a in range(H) if math.gcd(a, 10) == 1])
    s1 = S0(U[:, None], H).sum()
    i, j = np.triu_indices(len(U), 1)
    s2 = S0(np.stack([U[i], U[j]], 1), H).sum()
    s3 = 0.0
    it = combinations(U.tolist(), 3)
    while True:
        ch = [x for _, x in zip(range(200000), it)]
        if not ch:
            break
        s3 += S0(np.array(ch), H).sum()
    return s1, s2, s3


def f(lam, s):
    s1, s2, s3 = s
    m1 = s1 / lam; m2 = m1 + 2 * s2 / lam ** 2; m3 = m1 + 6 * s2 / lam ** 2 + 6 * s3 / lam ** 3
    return (m3 - 3 * m1 * m2 + 2 * m1 ** 3) / m1


def block_theory(L, W, s):
    x = np.linspace(L, L + W, 2001)
    return float(np.mean(f(np.log(x), s)))


if __name__ == '__main__':
    out = {}
    L, W = 10 ** 7, 9 * 10 ** 7
    M = mask(L, W)
    for H in (50, 500):
        s = sums(H)
        th, me = [], []
        for b in range(10):
            Lb = L + b * W // 10
            c = M[b * W // 10:(b + 1) * W // 10].reshape(-1, H).sum(1).astype(float)
            mid = Lb + np.arange(len(c)) * H + H / 2
            e = H / np.log(mid); r = c - e
            me.append(float(np.mean((r - r.mean()) ** 3) / e.mean()))
            th.append(block_theory(Lb, W // 10, s))
        d = np.array(me) - np.array(th)
        out[H] = {'sums': [float(v) for v in s], 'meas': float(np.mean(me)), 'law': float(np.mean(th)),
                  't9': float(d.mean() / (d.std(ddof=1) / math.sqrt(10))),
                  's1_over_H(should be 1)': float(s[0] / H)}
        print(H, out[H], flush=True)
    json.dump(out, open('/home/bas/scratch/study-review2/s10_skew.json', 'w'), indent=1)
