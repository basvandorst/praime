"""Is 'independence across q' (joint CDF of T vs product of per-q CDFs within 0.0045) informative? Same metric on
the empirical first-order chain (own code, sampler from check74_75 logic re-implemented here) and on data subsamples."""
import math, json, numpy as np
P = np.load('p1e9.npy'); Q = 30030; QS = (3, 5, 7, 11, 13); KMAX = 150
CLS = np.array([a for a in range(Q) if math.gcd(a, Q) == 1]); IDX = np.full(Q, -1); IDX[CLS] = np.arange(len(CLS))
g = np.diff(P); cl = IDX[P[1:] % Q]; o = np.argsort(cl, kind='stable'); gs = g[o]; st = np.searchsorted(cl[o], np.arange(len(CLS))); c = np.bincount(cl)
rng = np.random.default_rng(9); n = 200000
cur = CLS[rng.integers(0, len(CLS), n)]; D = np.zeros(n, np.int64)
seen = {q: np.zeros((n, q), bool) for q in QS}; f = {q: np.full(n, KMAX + 1) for q in QS}
for q in QS: seen[q][:, 0] = True
for k in range(1, KMAX + 1):
    r = IDX[cur]; gg = gs[st[r] + (rng.random(n) * c[r]).astype(np.int64)]; D += gg; cur = (cur - gg) % Q
    for q in QS:
        seen[q][np.arange(n), D % q] = True; f[q][(seen[q].sum(1) == q - 1) & (f[q] > KMAX)] = k
T = np.max(np.stack([f[q] for q in QS]), 0)
def metric(idx):
    j = np.array([np.mean(T[idx] <= k) for k in range(KMAX + 1)]); p = np.array([np.prod([np.mean(f[q][idx] <= k) for q in QS]) for k in range(KMAX + 1)])
    return float(np.max(np.abs(j - p))), float((j - p)[np.argmax(np.abs(j - p))])
print('chain all', metric(np.arange(n)))
print('chain 20k subsets', [metric(np.arange(i, i + 20000)) for i in range(0, n, 20000)][:5])
