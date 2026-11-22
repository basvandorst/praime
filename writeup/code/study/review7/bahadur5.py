"""Own computation of the order-4 and order-5 Bahadur terms (ln P) for a gap g, as defined in r8.py:
s_X = S(B u X)/S(B), B = {0, g}, X subset of the live interior; K_k = joint cumulant via set partitions of the
s-values; term_k = (-1)^k sum_X K_k / lam^k / prod_{i in X} (1 - s_i/lam).
Own singular series: per-subset product over primes p <= g (residue counting with np.sort) times the closed
constant for p > g (all points distinct mod p), the constant from primes up to 2e7 (tail < 1e-9)."""
import sys, math, json, time
from itertools import combinations
import numpy as np

LAMS = np.array([16.0, 18.0, 20.0, 22.0, 24.0, 27.0])


def primes_upto(N):
    s = np.ones(N + 1, dtype=bool); s[:2] = False
    for p in range(2, math.isqrt(N) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0]


def bell_partitions(items):
    if not items:
        yield []
        return
    f, rest = items[0], items[1:]
    for p in bell_partitions(rest):
        yield [[f]] + p
        for i in range(len(p)):
            yield p[:i] + [[f] + p[i]] + p[i + 1:]


g = int(sys.argv[1])
PR = primes_upto(20_000_000)
small = PR[PR <= g]; big = PR[PR > g].astype(float)
B = np.array([0, g])


def const(kx):          # prod over p > g of (1 - (2+kx)/p) / ((1 - 2/p) (1 - 1/p)^kx)
    return math.exp(float(np.sum(np.log1p(-(2 + kx) / big) - np.log1p(-2 / big) - kx * np.log1p(-1 / big))))


CONST = {k: const(k) for k in range(1, 6)}


def s_of(pts):
    """pts: (N, kx) interior points -> S(B u X)/S(B)."""
    N, kx = pts.shape
    full = np.concatenate([np.tile(B, (N, 1)), pts], 1)
    L = np.zeros(N)
    for p in small:
        p = int(p)
        R = np.sort(full % p, axis=1)
        nu = 1 + np.count_nonzero(np.diff(R, axis=1), axis=1)
        nuB = len({0, g % p})
        with np.errstate(divide='ignore'):
            L += np.log(np.maximum(1 - nu / p, 0)) - math.log(1 - nuB / p) - kx * math.log(1 - 1 / p)
    return np.exp(L) * CONST[kx]


interior = np.arange(1, g)
s1all = s_of(interior[:, None])
live = interior[s1all > 0]; s1 = s1all[s1all > 0]; n = len(live)
print('g', g, 'live', n, flush=True)
t0 = time.time()
s2 = np.zeros((n, n)); s3 = np.zeros((n, n, n)); s4 = np.zeros((n, n, n, n), dtype=np.float64)
I2 = np.array(list(combinations(range(n), 2))); v = s_of(live[I2])
for a, b in [(0, 1), (1, 0)]:
    s2[I2[:, a], I2[:, b]] = v
from itertools import permutations
I3 = np.array(list(combinations(range(n), 3))); v = s_of(live[I3])
for pm in permutations(range(3)):
    s3[tuple(I3[:, k] for k in pm)] = v
I4 = np.array(list(combinations(range(n), 4)))
v4 = np.concatenate([s_of(live[I4[i:i + 100000]]) for i in range(0, len(I4), 100000)])
for pm in permutations(range(4)):
    s4[tuple(I4[:, k] for k in pm)] = v4
print('tables', round(time.time() - t0), flush=True)
look = {1: lambda c: s1[c[0]], 2: lambda c: s2[c[0], c[1]], 3: lambda c: s3[c[0], c[1], c[2]],
        4: lambda c: s4[c[0], c[1], c[2], c[3]]}


def cumulant(cols, sk, k):
    K = np.zeros(len(cols[0]))
    for part in bell_partitions(list(range(k))):
        m = len(part)
        pr = np.full(len(cols[0]), (-1.0) ** (m - 1) * math.factorial(m - 1))
        for blk in part:
            pr = pr * (sk if len(blk) == k else look[len(blk)]([cols[j] for j in sorted(blk)]))
        K += pr
    return K


def add(cols, K, k, t):
    for i, lam in enumerate(LAMS):
        q = 1 - s1 / lam
        den = np.ones(len(K))
        for c in cols:
            den = den * q[c]
        t[i] += (-1) ** k * float(np.sum(K / lam ** k / den))


t4 = np.zeros(len(LAMS))
cols = [I4[:, j] for j in range(4)]
add(cols, cumulant(cols, v4, 4), 4, t4)
t5 = np.zeros(len(LAMS))
it = combinations(range(n), 5)
while True:
    Q = np.array([x for _, x in zip(range(150000), it)], dtype=np.int64)
    if len(Q) == 0:
        break
    cols = [Q[:, j] for j in range(5)]
    add(cols, cumulant(cols, s_of(live[Q]), 5), 5, t5)
print('secs', round(time.time() - t0))
ref = json.load(open('/home/bas/praime/writeup/data/study/r38.json')).get(str(g))
print('mine t4', np.round(t4, 6).tolist())
print('mine t5', np.round(t5, 6).tolist())
if ref:
    print('r38  t4', np.round(ref['order4_term_ln'], 6).tolist())
    print('r38  t5', np.round(ref['order5_term_ln'], 6).tolist())
    print('max rel diff t4 %.2e t5 %.2e' % (np.max(np.abs(t4 / np.array(ref['order4_term_ln']) - 1)),
                                           np.max(np.abs(t5 / np.array(ref['order5_term_ln']) - 1))))
