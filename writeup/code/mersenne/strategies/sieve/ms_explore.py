"""Exploration: per band, 1/k hit rate and clustering of factors (number of k<=K hits per p)."""
import numpy as np, gmpy2, math
N = 500_000
P, F = np.loadtxt("ms_first.txt", dtype=np.int64, unpack=True)
pairs = np.loadtxt("ms_pairs.txt", dtype=np.int64)
n = len(P); idx = {int(p): i for i, p in enumerate(P)}
small = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
K = 64
probs = [[] for _ in range(n)]   # admissible 1/k for k in 3..K (k=1 handled exactly)
k1 = np.zeros(n, dtype=bool)
for k in range(1, K + 1):
    qv = 2 * k * P + 1
    m8 = (qv % 8 == 1) | (qv % 8 == 7)
    proper = np.ones(n, dtype=bool)
    sm = P < 64
    proper[sm] = [int(x) < (1 << int(p)) - 1 for x, p in zip(qv[sm], P[sm])]
    cand = m8 & proper
    for r in small:
        cand &= ~((qv % r == 0) & (qv != r))
    for i in np.nonzero(cand)[0]:
        if gmpy2.is_prime(int(qv[i])):
            if k == 1: k1[i] = True
            else: probs[i].append(1.0 / k)
hits = np.zeros(n, dtype=np.int64)    # hits with 2<=k<=K
for p, k in pairs:
    if 2 <= k <= K: hits[idx[int(p)]] += 1
bands = [(2, 10 ** 4), (10 ** 4, 10 ** 5), (10 ** 5, N + 1)]
for lo, hi in bands:
    b = np.nonzero((P >= lo) & (P < hi) & ~k1)[0]   # primes without the k=1 factor
    exp_h = sum(sum(probs[i]) for i in b); var_h = sum(sum(x * (1 - x) for x in probs[i]) for i in b)
    obs_h = int(hits[b].sum())
    # Poisson-binomial distribution of hits per p, predicted counts for 0,1,2,>=3
    pred = np.zeros(4)
    for i in b:
        d = np.zeros(len(probs[i]) + 1); d[0] = 1
        for x in probs[i]:
            d[1:] = d[1:] * (1 - x) + d[:-1] * x; d[0] *= (1 - x)
        pred[0] += d[0]; pred[1] += d[1] if len(d) > 1 else 0
        pred[2] += d[2] if len(d) > 2 else 0; pred[3] += d[3:].sum() if len(d) > 3 else 0
    obs = [int(np.sum(hits[b] == 0)), int(np.sum(hits[b] == 1)), int(np.sum(hits[b] == 2)), int(np.sum(hits[b] >= 3))]
    print(f"[{lo},{hi}) no k=1 factor: n={len(b)} hits k=3..{K}: obs {obs_h} exp {exp_h:.1f} z={(obs_h-exp_h)/math.sqrt(var_h):+.2f};"
          f" #p with 0/1/2/3+ hits obs {obs} pred {[round(x,1) for x in pred]}")
