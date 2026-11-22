"""(a) Synthetic check of the second-order polymer formula for a hard-core gas of pair polymers.
(b) Own Bahadur form (orders 2-4) and the overlap correction against own exact inclusion-exclusion, gaps 52 and 70."""
import math, itertools, random, sys, json
from functools import lru_cache
import numpy as np
from sympy import primerange
sys.path.insert(0, '/home/bas/scratch/study-review4')
from hl_exact import U_sums, tail_T

# ---------- (a) synthetic: Z = sum over sets of pairwise-disjoint pairs (matchings) of prod w ----------
def logZ_matchings(n, w):
    @lru_cache(None)
    def Z(mask):                       # vertices still free = bits in mask
        if mask == 0:
            return 1.0
        i = (mask & -mask).bit_length() - 1
        rest = mask & ~(1 << i)
        tot = Z(rest)                  # i unmatched
        m = rest
        while m:
            j = (m & -m).bit_length() - 1
            tot += w[i][j] * Z(rest & ~(1 << j))
            m &= m - 1
        return tot
    return math.log(Z((1 << n) - 1))


def second_order(n, w, variant):
    first = sum(w[i][j] for i in range(n) for j in range(i + 1, n))
    sq = sum(w[i][j] ** 2 for i in range(n) for j in range(i + 1, n))
    ov = sum(w[i][j] * w[i][k] for i in range(n) for j in range(n) for k in range(n) if len({i, j, k}) == 3)
    if variant == 'claimed':
        return first - 0.5 * (sq + ov)
    if variant == 'no_self':           # B = B' left out
        return first - 0.5 * ov
    if variant == 'first':
        return first


random.seed(1)
n = 8
base = [[0.0] * n for _ in range(n)]
for i in range(n):
    for j in range(i + 1, n):
        base[i][j] = base[j][i] = random.uniform(-1, 1)
print('(a) synthetic hard-core pair gas, n = 8, activities eps * U(-1,1):')
for eps in (0.08, 0.04, 0.02, 0.01):
    w = [[eps * base[i][j] for j in range(n)] for i in range(n)]
    ex = logZ_matchings(n, w)
    print('  eps %.2f  err first-order %.2e   err claimed 2nd-order %.2e   err without B=B\' %.2e' % (
        eps, second_order(n, w, 'first') - ex, second_order(n, w, 'claimed') - ex, second_order(n, w, 'no_self') - ex))

# ---------- (b) Bahadur form on prime data ----------
def S(H, g):
    H = sorted(set(H)); k = len(H)
    v = 1.0
    for p in primerange(2, g + 1):
        nu = len({h % p for h in H})
        if nu == p:
            return 0.0
        v *= (1 - nu / p) / (1 - 1 / p) ** k
    return v * TT[(k, g)]


TT = {}
res = {}
for g in (52, 70):
    for k in range(2, 8):
        TT[(k, g)] = float(tail_T(k, g))
    U, _ = U_sums(g)
    Uf = [float(u) for u in U]
    B = [0, g]
    SB = S(B, g)
    I = [m for m in range(1, g) if S(B + [m], g) > 0]
    n = len(I)
    rho = np.array([S(B + [m], g) / SB for m in I])
    s2 = np.zeros((n, n))
    for a, b in itertools.combinations(range(n), 2):
        s2[a, b] = s2[b, a] = S(B + [I[a], I[b]], g) / SB
    s3 = {}
    for t in itertools.combinations(range(n), 3):
        s3[t] = S(B + [I[x] for x in t], g) / SB
    s4 = {}
    for t in itertools.combinations(range(n), 4):
        s4[t] = S(B + [I[x] for x in t], g) / SB
    for lam in (17.8, 23.0, 34.5):
        exact = sum((-1) ** j * Uf[j] / lam ** (j + 1) for j in range(len(Uf)))
        p = rho / lam
        m2 = s2 / lam ** 2
        L1 = float(np.sum(np.log1p(-p)))
        k2 = {}
        for a, b in itertools.combinations(range(n), 2):
            k2[(a, b)] = (m2[a, b] - p[a] * p[b]) / ((1 - p[a]) * (1 - p[b]))
        L2 = L1 + sum(k2.values())
        L3 = L2
        for (a, b, c), v in s3.items():
            k3 = v / lam ** 3 - m2[a, b] * p[c] - m2[a, c] * p[b] - m2[b, c] * p[a] + 2 * p[a] * p[b] * p[c]
            L3 -= k3 / ((1 - p[a]) * (1 - p[b]) * (1 - p[c]))
        L4 = L3
        for (a, b, c, d), v in s4.items():
            m3 = lambda x, y, z: s3[tuple(sorted((x, y, z)))] / lam ** 3
            k4 = (v / lam ** 4 - m3(a, b, c) * p[d] - m3(a, b, d) * p[c] - m3(a, c, d) * p[b] - m3(b, c, d) * p[a]
                  - m2[a, b] * m2[c, d] - m2[a, c] * m2[b, d] - m2[a, d] * m2[b, c]
                  + 2 * (m2[a, b] * p[c] * p[d] + m2[a, c] * p[b] * p[d] + m2[a, d] * p[b] * p[c]
                         + m2[b, c] * p[a] * p[d] + m2[b, d] * p[a] * p[c] + m2[c, d] * p[a] * p[b])
                  - 6 * p[a] * p[b] * p[c] * p[d])
            L4 += k4 / ((1 - p[a]) * (1 - p[b]) * (1 - p[c]) * (1 - p[d]))
        # polymer second order for pair polymers
        sq = sum(v * v for v in k2.values())
        K = np.zeros((n, n))
        for (a, b), v in k2.items():
            K[a, b] = K[b, a] = v
        ov = sum(K[i, j] * K[i, k] for i in range(n) for j in range(n) for k in range(n) if j != k and j != i and k != i)
        corr = -0.5 * (sq + ov)
        sh = lambda L: SB / lam * math.exp(L) / exact - 1
        row = {'n_live': n, 'exact': exact, 'o2': sh(L2), 'o3': sh(L3), 'o4': sh(L4), 'o3+corr': sh(L3 + corr), 'o4+corr': sh(L4 + corr),
               'corr': corr}
        res['%d@%.1f' % (g, lam)] = row
        print('(b) gap %d lam %.1f live %d: exact %.6e  rel. bias o2 %+.3f%% o3 %+.3f%% o4 %+.3f%%  o3+corr %+.3f%%  o4+corr %+.3f%%' % (
            g, lam, n, exact, 100 * row['o2'], 100 * row['o3'], 100 * row['o4'], 100 * row['o3+corr'], 100 * row['o4+corr']), flush=True)
json.dump(res, open('/home/bas/scratch/study-review4/polymer.json', 'w'), indent=1)
