# own Bahadur order-k terms for a gap g: ln P(no prime inside | 0, g prime) term of order k
#   = (-1)^k sum_{k-sets} kappa_k(s)/(lam^k prod q),  s_X = S({0,g} u X)/S({0,g}), q_i = 1 - s_i/lam
import numpy as np, math, sys, json, time
from itertools import combinations
from math import comb

def primes_upto(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]

def setup(g):
    SP = [int(p) for p in primes_upto(g)]            # p <= g: residues may collide
    BIG = primes_upto(20_000_000); BIG = BIG[BIG > g].astype(float)
    # C_j = prod_{p>g} (1-(2+j)/p)/(1-2/p) * (1-1/p)^(-j)
    Cj = {j: math.exp(float(np.sum(np.log1p(-(2+j)/BIG) - np.log1p(-2/BIG) - j*np.log1p(-1/BIG)))) for j in range(1, 6)}
    return SP, Cj

def logS_small(pts, SP):
    """sum over p<=g of log[(1-nu_p/p)(1-1/p)^-k] for rows of pts (N,k); -inf if inadmissible"""
    N, k = pts.shape
    out = np.zeros(N)
    for p in SP:
        r = np.sort(pts % p, axis=1)
        nu = 1 + np.count_nonzero(np.diff(r, axis=1), axis=1)
        with np.errstate(divide='ignore'):
            out += np.log(1 - nu/p) - k*math.log(1 - 1/p)
    return out

def s_values(X, g, SP, Cj, base):
    """s_X for rows X (N,j) of interior points"""
    N, j = X.shape
    pts = np.concatenate([np.zeros((N, 1), np.int64), np.full((N, 1), g, np.int64), X], 1)
    return np.exp(logS_small(pts, SP) - base) * Cj[j]

def set_partitions(s):
    if not s:
        yield []; return
    first, rest = s[0], s[1:]
    for p in set_partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i+1:]
        yield [[first]] + p

def combidx(sorted_cols):
    """combinatorial number system index of sorted index tuples (columns ascending)"""
    idx = np.zeros(sorted_cols.shape[0], np.int64)
    for t in range(sorted_cols.shape[1]):
        c = sorted_cols[:, t]
        idx += np.array([comb(int(v), t+1) for v in range(c.max()+1)], np.int64)[c]
    return idx

def run(g, K, lams, chunk=400000):
    t0 = time.time()
    SP, Cj = setup(g)
    base = logS_small(np.array([[0, g]], np.int64), SP)[0]
    I = np.arange(1, g, dtype=np.int64)
    s1 = s_values(I[:, None], g, SP, Cj, base)
    live = I[s1 > 0]; s1 = s1[s1 > 0]; n = len(live)
    tables = {1: s1}
    for j in range(2, K):
        Q = np.array(list(combinations(range(n), j)), np.int64)
        v = np.concatenate([s_values(live[Q[i:i+chunk]], g, SP, Cj, base) for i in range(0, len(Q), chunk)])
        tables[j] = v      # stored in lexicographic order of combinations == need combidx order
        # reorder to combinadic index
        order = np.empty(len(Q), np.int64); order[combidx(Q)] = np.arange(len(Q)); tables[j] = v[order]
        del Q
    parts = list(set_partitions(list(range(K))))
    coef = [(-1)**(len(p)-1)*math.factorial(len(p)-1) for p in parts]
    blocks = sorted({tuple(sorted(b)) for p in parts for b in p if len(b) < K}, key=lambda b: (len(b), b))
    tot = np.zeros(len(lams)); it = combinations(range(n), K)
    while True:
        Q = np.array([x for _, x in zip(range(chunk), it)], np.int64)
        if len(Q) == 0: break
        m = {}
        for b in blocks:
            cols = Q[:, list(b)]
            m[b] = tables[1][cols[:, 0]] if len(b) == 1 else tables[len(b)][combidx(cols)]
        m[tuple(range(K))] = s_values(live[Q], g, SP, Cj, base)
        kap = np.zeros(len(Q))
        for p, c in zip(parts, coef):
            pr = np.full(len(Q), float(c))
            for b in p: pr *= m[tuple(sorted(b))]
            kap += pr
        for li, lam in enumerate(lams):
            q = 1 - s1/lam
            den = np.ones(len(Q))
            for t in range(K): den *= q[Q[:, t]]
            tot[li] += (-1)**K * float(np.sum(kap/den)) / lam**K
    return n, tot, time.time()-t0

if __name__ == '__main__':
    g, K = int(sys.argv[1]), int(sys.argv[2])
    lams = [float(x) for x in sys.argv[3].split(',')]
    n, tot, dt = run(g, K, lams)
    print(json.dumps({'g': g, 'order': K, 'live': n, 'lams': lams, 'term': tot.tolist(), 'secs': round(dt)}), flush=True)
