"""K by an independent sum (primes to 1e8 by segmented sieve + integral tail), brute-force mean over even g of S ln S,
measured gap entropy (plug-in + Miller-Madow) of [1e7, 1e8), and intercept fits of the round-52 differences."""
import math, json
import numpy as np
import mpmath

def seg_primes(lo, hi, base):
    m = np.ones(hi - lo, bool)
    for p in base:
        if p * p >= hi: break
        s = max(p * p, ((lo + p - 1) // p) * p)
        m[s - lo::p] = False
    if lo <= 1: m[:2 - lo] = False
    return np.nonzero(m)[0] + lo

base = seg_primes(0, 10 ** 4 + 1, [])
base = seg_primes(0, 10 ** 4 + 1, base)
# small base up to 1e4 via the sieve above needs the primes to 100; redo cleanly
def simple(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for p in range(2, int(n ** .5) + 1):
        if s[p]: s[p * p::p] = False
    return np.nonzero(s)[0]
base = simple(10 ** 4 + 10)
# K: sum over odd p <= 1e8 by segments; tail sum_{p>X} ln(1+1/(p-2))/(p-1) ~ int_X^inf dt/(t^2 ln t) ~ 1/(X ln X)
K_sum = 0.0; lnC2 = 0.0; X = 10 ** 8
for lo in range(0, X, 10 ** 7):
    P = seg_primes(lo, lo + 10 ** 7, base).astype(np.float64)
    P = P[P > 2]
    K_sum += float(np.sum(np.log1p(1 / (P - 2)) / (P - 1)))
    lnC2 += float(np.sum(np.log1p(-1 / (P - 1) ** 2)))
tail = 1 / (X * math.log(X))
K = math.log(2) + lnC2 + K_sum + tail
print('K', K, 'C2', math.exp(lnC2), '(K-ln2)/ln2 bits', (K - math.log(2)) / math.log(2))
# brute force mean_{g even <= G} S ln S / 2
G = 2 * 10 ** 6
k = np.arange(1, G // 2 + 1); lnS = np.full(len(k), math.log(2) + lnC2)
for p in simple(G // 2 + 1)[1:]:
    lnS[p - 1::p] += math.log((p - 1) / (p - 2))
S = np.exp(lnS)
print('brute K (g<=2e6)', 0.5 * float(np.mean(S * lnS)), 'mean S', float(np.mean(S)))
# measured entropy over [1e7, 1e8): gaps with lower prime in the window
cnt = {}; last = None; lam_sum = 0.0; n = 0; asym_sum = 0.0
for lo in range(10 ** 7, 10 ** 8 + 10 ** 7, 10 ** 7):
    hi = lo + 10 ** 7
    P = seg_primes(lo, hi, base).astype(np.int64)
    if lo >= 10 ** 8:
        P = P[:1]   # only the first prime above 1e8 closes the last gap
    if last is not None: P = np.concatenate([[last], P])
    g = np.diff(P); low = P[:-1]
    v, c = np.unique(g, return_counts=True)
    for a, b in zip(v.tolist(), c.tolist()): cnt[a] = cnt.get(a, 0) + b
    lam_sum += float(np.sum(np.log(low.astype(np.float64)))); n += len(g)
    asym_sum += float(np.sum(np.log2(math.e * np.log(low.astype(np.float64)))))
    last = int(P[-1])
c = np.array(list(cnt.values()), float); N = c.sum(); p = c / N
H = -np.sum(p * np.log2(p)) + (len(c) - 1) / (2 * N * math.log(2))
print('gaps', int(N), 'mean lam', lam_sum / n, 'H data', H, 'asym at mean lam', math.log2(math.e * lam_sum / n) - K / math.log(2),
      'mean of asym over primes', asym_sum / n - K / math.log(2))
r = json.load(open('/home/bas/praime/writeup/data/study/r52.json'))
lam = np.array([r[w]['mean_lambda'] for w in ('train', '1e9-2e9', 'A', 'B', 'C')])
d = np.array([r[w]['H_data'] - r[w]['H_asym_at_mean_lambda'] for w in ('train', '1e9-2e9', 'A', 'B', 'C')])
dd = np.array([r[w]['diff'] for w in ('train', '1e9-2e9', 'A', 'B', 'C')])
for name, X_ in (('1/lam', 1 / lam), ('ln lam/lam', np.log(lam) / lam)):
    A_ = np.vstack([np.ones(5), X_]).T
    print('fit a + b*' + name, 'data-asym:', np.linalg.lstsq(A_, d, rcond=None)[0], ' data-discrete:', np.linalg.lstsq(A_, dd, rcond=None)[0])
