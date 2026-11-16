# Item G, stage 2: exact bin means of Delta(n) and of -S(n) (mod 12 strands), sampled pred(n).
import numpy as np, math, time, os, sys
from sympy import primerange, isprime
d = os.path.dirname(os.path.abspath(__file__))
P = {a: np.load(os.path.join(d, 'G_p%d.npy' % a)) for a in (1, 5, 7, 11)}
bits = np.load(os.path.join(d, 'G_bits.npy'))
edges = np.unique(np.round(np.exp(np.linspace(math.log(1e6), math.log(1e9), 61)))).astype(np.int64)
print('edges:', len(edges), edges[:3], edges[-3:])
nb = len(edges) - 1

def isp_odd(x):
    i = x >> 1
    return (bits[i >> 3] >> (7 - (i & 7)).astype(np.uint8)) & 1

# spot-check the bit array and class arrays
rng = np.random.default_rng(12345)
xs = rng.integers(0, 5 * 10**8, 2000) * 2 + 1
assert all(int(v) == isprime(int(x)) for v, x in zip(isp_odd(xs), xs))
for a in P:
    assert all(isprime(int(x)) and int(x) % 12 == a for x in P[a][rng.integers(0, len(P[a]), 200)])

# G_ab(t) = #{(p, q): p = a, q = b mod 12, p + q < t}
def G(a, b, t):
    pa = P[a][:np.searchsorted(P[a], t)]
    keys = (np.uint32(t - 1) - pa)[::-1]          # ascending, >= 0
    return int(np.searchsorted(P[b], keys, side='right').sum())

# check G against brute force on a small t
t = 20000
pr = [p for p in primerange(5, t)]
for a, b in [(1, 11), (5, 7), (1, 1), (11, 11)]:
    bf = sum(1 for p in pr if p % 12 == a for q in pr if q % 12 == b and p + q < t)
    assert bf == G(a, b, t), (a, b)

SPEC = {0: [((1, 11), 2), ((5, 7), -2)],          # (1,11)+(11,1) - (5,7) - (7,5); symmetric sums
        2: [((1, 1), 1), ((7, 7), -1)],
        6: [((1, 5), 2), ((7, 11), -2)],
        10: [((5, 5), 1), ((11, 11), -1)]}
pairs = sorted({pq for c in SPEC for pq, _ in SPEC[c]})
t0 = time.time()
Gt = {pq: np.array([G(pq[0], pq[1], int(t)) for t in edges], np.int64) for pq in pairs}
print('G done in %.0f s' % (time.time() - t0), flush=True)
# check symmetry (a,b) vs (b,a) on bin sums for one pair
gb = np.array([G(11, 1, int(t)) for t in edges[-3:]])
assert np.all(gb == Gt[(1, 11)][-3:])

def ncount(c, lo, hi):        # number of n = c mod 12 in [lo, hi)
    f = lambda x: (x - c + 11) // 12      # count of n = c mod 12 in [0, x)
    return f(hi) - f(lo)

meanD = {c: np.zeros(nb) for c in SPEC}
for c in SPEC:
    tot = sum(w * Gt[pq] for pq, w in SPEC[c])
    for i in range(nb):
        meanD[c][i] = (tot[i + 1] - tot[i]) / ncount(c, edges[i], edges[i + 1])

# exact mean of S(n) per bin: for n = c mod 12 and s >= 5 prime (s^2 = 1 mod 12), r = n - s^2 = c - 1 mod 12
S_primes = np.array(list(primerange(5, math.isqrt(10**9) + 1)), np.int64)
rcls = {0: 11, 2: 1, 6: 5, 10: None}
meanS = {c: np.zeros(nb) for c in SPEC}
for c in (0, 2, 6):
    Pr = P[rcls[c]]
    for i in range(nb):
        lo, hi = edges[i], edges[i + 1]
        sq = S_primes[S_primes ** 2 < hi] ** 2
        # r in [lo - s^2, hi - s^2), r >= 1  (s^2 < n  <=>  r >= 1)
        cnt = np.searchsorted(Pr, np.maximum(hi - sq, 0), 'left') - np.searchsorted(Pr, np.maximum(lo - sq, 0), 'left')
        meanS[c][i] = cnt.sum() / ncount(c, lo, hi)
# n = 10 mod 12: S(n) = 0 identically (r = 9 mod 12); S(n-2) uses r = 7 mod 12
meanS2_10 = np.zeros(nb)
for i in range(nb):
    lo, hi = edges[i] - 2, edges[i + 1] - 2
    sq = S_primes[S_primes ** 2 < hi] ** 2
    cnt = np.searchsorted(P[7], np.maximum(hi - sq, 0), 'left') - np.searchsorted(P[7], np.maximum(lo - sq, 0), 'left')
    meanS2_10[i] = cnt.sum() / ncount(10, edges[i], edges[i + 1])

# sampled pred(n) = S(n) * SIG(n) / Q(n); 2000 n per bin and class
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
def S_of(n):
    out = np.zeros(len(n), np.int64)
    for s in S_primes:
        x = n - s * s
        ok = x >= 1
        if not ok.any(): break
        out[ok] += isp_odd(x[ok])
    return out
def pi_(a, n): return np.searchsorted(P[a], n, 'right').astype(np.int64)
def SIG(c, n):
    if c == 0: return 2 * (pi_(1, n) + pi_(11, n)) - 2 * (pi_(5, n) + pi_(7, n))
    if c == 2: return 2 * pi_(1, n) - 2 * pi_(7, n)
    if c == 6: return 2 * (pi_(1, n) + pi_(5, n)) - 2 * (pi_(7, n) + pi_(11, n))
    return 2 * pi_(5, n) - 2 * pi_(11, n)
predm = {c: np.zeros(nb) for c in SPEC}; preds = {c: np.zeros(nb) for c in SPEC}
Sm_samp = np.zeros(nb)
t0 = time.time()
for c in SPEC:
    ns = []
    for i in range(nb):
        lo, hi = int(edges[i]), int(edges[i + 1])
        first = lo + ((c - lo) % 12)
        k = (hi - 1 - first) // 12 + 1
        ns.append(first + 12 * rng.integers(0, k, NS))
    n = np.concatenate(ns).astype(np.int64)
    S = S_of(n - 2) if c == 10 else S_of(n)
    Q = np.searchsorted(S_primes ** 2, n, 'right')
    pr = S * SIG(c, n) / Q
    pr = pr.reshape(nb, NS);
    predm[c] = pr.mean(1); preds[c] = pr.std(1, ddof=1) / math.sqrt(NS)
    if c == 0:
        Sm_samp = S.reshape(nb, NS).mean(1)
print('pred sampling done in %.0f s' % (time.time() - t0), flush=True)
# exact check of a few Delta values per n against direct pair counting (small n)
np.savez(os.path.join(d, 'G_results.npz'), edges=edges, **{'D%d' % c: meanD[c] for c in SPEC},
         **{'S%d' % c: meanS[c] for c in SPEC}, S2_10=meanS2_10, **{'pm%d' % c: predm[c] for c in SPEC},
         **{'ps%d' % c: preds[c] for c in SPEC})

want = [1_000_000, 100_000_000, 501_187_234, 891_250_938]
for w in want:
    i = int(np.nonzero(edges == w)[0][0])
    print('bin [%d, %d):' % (edges[i], edges[i + 1]))
    for c in SPEC:
        extra = '  mean -S(n) = %.2f' % (-meanS[c][i]) if c != 10 else '  mean -S(n-2) = %.2f' % (-meanS2_10[i])
        print('   n=%2d mod 12: mean Delta = %9.3f   pred (sample of %d) = %9.2f +- %.2f%s' % (
            c, meanD[c][i], NS, predm[c][i], preds[c][i], extra))
    print('   sampled mean -S (n=0) = %.2f (exact %.2f)' % (-Sm_samp[i], -meanS[0][i]))
for c in SPEC:
    r = np.corrcoef(meanD[c], predm[c])[0, 1]
    print('corr over %d bins, n=%2d: Delta vs sampled pred = %.4f   mean |Delta - pred| = %.2f' % (nb, c, r, np.mean(np.abs(meanD[c] - predm[c]))))
