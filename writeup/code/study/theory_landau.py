"""Round 5: consecutive n with n^2 + 1 prime, from Bateman-Horn singular series of shifts (preregister/round5.txt).

S_f(B) = prod_p (1 - omega_B(p)/p)/(1 - 1/p)^k,  omega_B(p) = |{(r - b) mod p : r root of x^2+1 mod p, b in B}|.
Computation: p <= Ps (Ps >= span and >= 2k+2) by the residues; Ps < p <= 2e7 by the generic factor
g_k(p) = (1 - k w_p/p)/(1 - 1/p)^k, w_p = 1 + chi_4(p), through cumulative sums; plus, for primes p > Ps dividing
d^2 + 4 for a difference d of B (then two shifted roots coincide: r - b = r' - b' with r - r' = 2i, d^2 = -4),
the exact factor with omega = 2k - (number of such pairs). Primes p > Ps cannot divide a difference d <= span.

usage: theory_landau.py check | build | values
"""
import sys, os, math, json, pickle, time
from itertools import combinations
import numpy as np
from common import small_primes

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'data', 'study')
CACHE = os.environ.get('STUDY_CACHE5', '/tmp/study_landau.pkl')
X = 20_000_000
_P = small_primes(X).astype(np.int64)
_Pf = _P.astype(np.float64)
_W = np.where(_P == 2, 1, np.where(_P % 4 == 1, 2, 0)).astype(np.float64)
KMAX = 7
_CS = {k: np.concatenate(([0.0], np.cumsum(np.where(_P > 60, np.log(np.abs(1 - k * _W / _Pf)) - k * np.log1p(-1 / _Pf), 0.0))))
       for k in range(1, KMAX + 1)}          # p <= 60 is always done by the residues (Ps >= 60)
SPF = None
FWIN = {'F1': (1_000_000, 2_000_000), 'F2': (100_000_000, 1_000_000), 'F3': (10_000_000_000, 200_000)}
LAMG = np.linspace(2 * math.log(1e6) - 0.01, math.log((1e10 + 2e5) ** 2 + 1) + 0.01, 48)
HMAX = 400


def roots(p):
    if p == 2:
        return [1]
    if p % 4 == 3:
        return []
    return [x for x in range(p) if (x * x + 1) % p == 0] if p < 2000 else None


_ROOTS = {int(p): roots(int(p)) for p in _P[_P < 2000]}


def _spf(n):
    global SPF
    if SPF is None or len(SPF) <= n:
        N = max(n + 1, 2_000_000)
        s = np.zeros(N, dtype=np.int64)
        for p in small_primes(int(N ** 0.5) + 1):
            p = int(p)
            m = s[p::p] == 0
            s[p::p][m] = p
        s[s == 0] = np.arange(N)[s == 0]
        SPF = s
    return SPF


def prime_factors(n):
    spf = _spf(n)
    out = set()
    while n > 1:
        p = int(spf[n]); out.add(p); n //= p
    return out


def Sf_many(arr):
    arr = np.asarray(arr, dtype=np.int64)
    n, k = arr.shape
    span = int((arr.max(1) - arr.min(1)).max())
    Ps = max(span, 2 * k + 2, 60)
    logS = np.zeros(n)
    dead = np.zeros(n, dtype=bool)
    for p in _P[_P <= Ps]:
        p = int(p)
        R = _ROOTS[p]
        if not R:
            logS += -k * math.log(1 - 1 / p)
            continue
        vals = np.concatenate([(r - arr) % p for r in R], axis=1)
        vals.sort(axis=1)
        om = 1 + np.count_nonzero(np.diff(vals, axis=1), axis=1)
        f = 1 - om / p
        dead |= f <= 0
        logS += np.log(np.where(f > 0, f, 1.0)) - k * math.log(1 - 1 / p)
    i0 = int(np.searchsorted(_P, Ps, side='right'))
    logS += _CS[k][-1] - _CS[k][i0]
    # collisions at large p: p | d^2 + 4
    for row in range(n):
        if dead[row]:
            continue
        b = arr[row]
        cnt = {}
        for x, y in combinations(b.tolist(), 2):
            for p in prime_factors((y - x) ** 2 + 4):
                if p > Ps:
                    cnt[p] = cnt.get(p, 0) + 1
        for p, c in cnt.items():
            if p <= X:
                logS[row] += math.log((1 - (2 * k - c) / p) / (1 - 2 * k / p))
    S = np.exp(logS)
    S[dead] = 0.0
    return S


class FPattern:
    """Bahadur form (as theory2.Pattern) with S_f; B = {0, g}, interior points 0 < m < g."""

    def __init__(self, B, order):
        B = sorted(B)
        s = B[-1]
        self.SB = float(Sf_many(np.array([B]))[0])
        self.k = len(B)
        self.s1 = np.zeros(0); self.pairs = self.trip = self.quad = None
        if self.SB == 0:
            return
        I = np.array([m for m in range(1, s) if m not in B], dtype=np.int64)
        Barr = np.array(B, dtype=np.int64)
        Ssub = lambda cols: Sf_many(np.concatenate([np.tile(Barr, (cols.shape[0], 1)), cols], 1)) / self.SB
        if len(I) == 0:
            return
        s1 = Ssub(I[:, None]); live = I[s1 > 0]; s1 = s1[s1 > 0]; n = len(live)
        self.s1 = s1
        if order >= 2 and n >= 2:
            i, j = np.triu_indices(n, 1)
            s2v = Ssub(np.stack([live[i], live[j]], 1))
            self.s2 = np.zeros((n, n)); self.s2[i, j] = s2v; self.s2[j, i] = s2v
            self.pairs = (i, j, s2v)
        if order >= 3 and n >= 3:
            T = np.array(list(combinations(range(n), 3)), dtype=np.int64)
            a, b, c = T[:, 0], T[:, 1], T[:, 2]
            s3v = Ssub(np.stack([live[a], live[b], live[c]], 1))
            self.trip = (a, b, c, s3v)
            if order >= 4 and n >= 4:
                s3 = np.zeros((n, n, n))
                for (x, y, z) in ((a, b, c), (a, c, b), (b, a, c), (b, c, a), (c, a, b), (c, b, a)):
                    s3[x, y, z] = s3v
                Q = np.array(list(combinations(range(n), 4)), dtype=np.int64)
                a, b, c, d = Q[:, 0], Q[:, 1], Q[:, 2], Q[:, 3]
                self.quad = (a, b, c, d, Ssub(np.stack([live[a], live[b], live[c], live[d]], 1)), s3)

    def log_avoid(self, lam, order):
        import theory2 as T2
        P = T2.Pattern.__new__(T2.Pattern)
        P.SB = self.SB; P.s1 = self.s1; P.pairs = self.pairs; P.trip = self.trip; P.quad = self.quad
        if self.pairs is not None:
            P.s2 = self.s2
        return P.log_avoid(lam, order)


def order_for(g):
    return 4 if g <= 60 else 3


def build():
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    t0 = time.time()
    for g in range(2, HMAX + 1, 2):
        if g in cache:
            continue
        o = order_for(g)
        pt = FPattern([0, g], o)
        la = np.array([pt.log_avoid(l, o) for l in LAMG]) if pt.SB > 0 else None
        cache[g] = (pt.SB, la)
        if time.time() - t0 > 450:
            break
    pickle.dump(cache, open(CACHE, 'wb'))
    print('cached', len(cache), 'of', HMAX // 2, flush=True)


def share(cache, g, N0, W):
    SB, la = cache[g]
    if SB == 0:
        return 0.0
    S0 = float(Sf_many(np.array([[0]]))[0])
    n = np.linspace(N0, N0 + W, 2001)
    lam = np.log(n * n + 1)
    w = 1 / lam
    f = SB / S0 / lam * np.exp(np.interp(lam, LAMG, la))
    return float(np.trapezoid(w * f, n) / np.trapezoid(w, n))


GS = [2, 4, 6, 8, 10, 12, 20, 30, 50]
QS = [3, 5, 7, 13]


def values(N0, W, cache):
    g = {h: share(cache, h, N0, W) for h in range(2, HMAX + 1, 2)}
    tot = sum(g.values())
    t = {'norm': tot}
    for h in GS:
        t['L.gap%d' % h] = g[h]
    t['L.ge60'] = sum(v for h, v in g.items() if h >= 60) / tot
    t['L.ge100'] = sum(v for h, v in g.items() if h >= 100) / tot
    for q in QS:
        t['L.same%d' % q] = sum(v for h, v in g.items() if h % q == 0) / tot
    return t


def main():
    if sys.argv[1] == 'check':
        print('S_f({0}) =', Sf_many(np.array([[0]]))[0], ' (Bateman-Horn constant of n^2+1: 1.3728134628)')
        print('S_f({0,2}) =', Sf_many(np.array([[0, 2]]))[0])
        return
    if sys.argv[1] == 'build':
        build(); return
    cache = pickle.load(open(CACHE, 'rb'))
    res = {}
    for win, (N0, W) in FWIN.items():
        res[win] = {'-1': values(N0, W, cache)}
        for b in range(10):
            res[win][str(b)] = values(N0 + b * (W // 10), W // 10, cache)
        print(win, {k: round(v, 5) for k, v in res[win]['-1'].items()}, flush=True)
    json.dump(res, open(os.path.join(OUT, 'theory_landau.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
