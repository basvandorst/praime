"""Round 14 (preregister/round14.txt): consecutive Sophie Germain primes from the singular series of linear forms.

Tokens: the form x + b is the token 2b, the form 2x + 2m + 1 the token 2m + 1. For an odd prime q both have the root
-t/2 (mod q), so nu_q = number of distinct tokens mod q; for q = 2 only the even tokens have roots, at (t/2) mod 2.
Beyond the span of the tokens no two roots coincide, so the generic factor G_k of theory2 applies.
usage: r14.py data | build | theory
"""
import sys, os, math, json, pickle, time
from itertools import combinations
import numpy as np
import theory2 as T2
import common as C

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
CACHE = os.environ.get('STUDY_CACHE14', '/tmp/study_r14.pkl')
WIN = dict(C.WINDOWS); WIN['C'] = (10 ** 15, 2 * 10 ** 7)
LAMG = np.linspace(math.log(1e7) - 0.01, math.log(1e15 + 2e7) + 0.01, 60)
HS = list(range(6, 301, 6))
W1 = [6, 12, 18, 24, 30, 36, 42, 48, 60, 90, 120]
LN2 = math.log(2)


def S_tok(arr):
    arr = np.asarray(arr, dtype=np.int64)
    n, k = arr.shape
    span = int((arr.max(1) - arr.min(1)).max())
    P = max(span, 3)
    logS = np.zeros(n); dead = np.zeros(n, dtype=bool)
    # q = 2: roots (t/2) mod 2 for the even tokens only
    ev = (arr % 2 == 0)
    r2 = np.where(ev, (arr // 2) % 2, -1)
    nu2 = (np.any(r2 == 0, axis=1)).astype(int) + (np.any(r2 == 1, axis=1)).astype(int)
    f = 1 - nu2 / 2
    dead |= f <= 0
    logS += np.log(np.where(f > 0, f, 1.0)) - k * math.log(0.5)
    for q in C.small_primes(P)[1:]:
        q = int(q)
        r = np.sort(arr % q, axis=1)
        nu = 1 + np.count_nonzero(np.diff(r, axis=1), axis=1)
        f = 1 - nu / q
        dead |= f <= 0
        logS += np.log(np.where(f > 0, f, 1.0)) - k * math.log(1 - 1 / q)
    S = np.exp(logS) * T2.G(k, P)
    S[dead] = 0.0
    return S


def S_sets(sets):
    out = np.zeros(len(sets)); by = {}
    for i, s in enumerate(sets):
        by.setdefault(len(s), []).append(i)
    for k, idx in by.items():
        out[idx] = S_tok(np.array([sorted(sets[i]) for i in idx], dtype=np.int64))
    return out


class SGPattern:
    """B = tokens of the two SG primes at 0 and h; interior events {2m, 2m+1}; moments s, and the numbers of new even
    (x-forms, weight 1/lam) and odd (2x-forms, weight 1/(lam + ln 2)) tokens."""

    def __init__(self, h, order=3):
        B = frozenset([0, 1, 2 * h, 2 * h + 1])
        self.SB = float(S_tok(np.array([sorted(B)]))[0])
        ev = [frozenset([2 * m, 2 * m + 1]) for m in range(1, h)]
        U1 = [B | e for e in ev]
        s1 = S_sets(U1) / self.SB
        live = [i for i in range(len(ev)) if s1[i] > 0]
        ev = [ev[i] for i in live]; self.s1 = s1[live]
        cnt = lambda u: (sum(1 for t in u - B if t % 2 == 0), sum(1 for t in u - B if t % 2 == 1))
        self.n1 = np.array([cnt(e) for e in ev]).reshape(-1, 2)
        n = len(ev); self.pairs = self.trip = None
        if n >= 2:
            i, j = np.triu_indices(n, 1)
            U = [B | ev[a] | ev[b] for a, b in zip(i, j)]
            v = S_sets(U) / self.SB
            nn = np.array([cnt(u) for u in U])
            self.s2 = np.zeros((n, n)); self.s2[i, j] = v; self.s2[j, i] = v
            self.n2 = np.zeros((n, n, 2), dtype=int); self.n2[i, j] = nn; self.n2[j, i] = nn
            self.pairs = (i, j)
        if order >= 3 and n >= 3:
            T = np.array(list(combinations(range(n), 3)), dtype=np.int64)
            U = [B | ev[a] | ev[b] | ev[c] for a, b, c in T]
            self.trip = (T[:, 0], T[:, 1], T[:, 2], S_sets(U) / self.SB, np.array([cnt(u) for u in U]))

    def log_avoid(self, lam, order):
        l2 = lam + LN2
        w = lambda nn: lam ** nn[..., 0] * l2 ** nn[..., 1]
        p = self.s1 / w(self.n1)
        out = float(np.sum(np.log1p(-p)))
        if order >= 2 and self.pairs is not None:
            m2 = self.s2 / w(self.n2)
            i, j = self.pairs
            out += float(np.sum((m2[i, j] - p[i] * p[j]) / ((1 - p[i]) * (1 - p[j]))))
            if order >= 3 and self.trip is not None:
                a, b, c, s3, n3 = self.trip
                k3 = s3 / w(n3) - m2[a, b] * p[c] - m2[a, c] * p[b] - m2[b, c] * p[a] + 2 * p[a] * p[b] * p[c]
                out -= float(np.sum(k3 / ((1 - p[a]) * (1 - p[b]) * (1 - p[c]))))
        return out


def build():
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    t0 = time.time()
    for h in HS:
        if h in cache:
            continue
        pt = SGPattern(h, 3)
        cache[h] = (pt.SB, np.array([pt.log_avoid(l, 3) for l in LAMG]), np.array([pt.log_avoid(l, 1) for l in LAMG]))
        if time.time() - t0 > 450:
            break
    pickle.dump(cache, open(CACHE, 'wb'))
    print('cached', len(cache), 'of', len(HS), flush=True)


def share(cache, h, L, W, k=1):
    SB, la3, la1 = cache[h]
    S0 = float(S_tok(np.array([[0, 1]]))[0])
    u = np.linspace(math.log(L), math.log(L + W), 2001)
    wgt = np.exp(u) / (u * (u + LN2))
    la = la3 if k == 1 else la1
    f = SB / S0 / (u * (u + LN2)) * np.exp(np.interp(u, LAMG, la))
    return float(np.trapezoid(wgt * f, u) / np.trapezoid(wgt, u))


def values(L, W, cache, k=1):
    g = {h: share(cache, h, L, W, k) for h in HS}
    t = {('W1.gap%d' % h): g[h] for h in W1}
    tot = sum(g.values())
    t['W2.le300'] = tot
    for q in (5, 7):
        t['W3.same%d' % q] = sum(v for h, v in g.items() if h % q == 0) / tot
    return t


def theory():
    cache = pickle.load(open(CACHE, 'rb'))
    res = {}
    for win, (L, W) in WIN.items():
        res[win] = {'-1': values(L, W, cache), 'order1': values(L, W, cache, 2)}
        for b in range(10):
            res[win][str(b)] = values(L + b * (W // 10), W // 10, cache)
        print(win, {k: round(v, 5) for k, v in res[win]['-1'].items()}, flush=True)
    json.dump(res, open(os.path.join(OUT, 'theory_r14.json'), 'w'), indent=1)


def data():
    with open(os.path.join(OUT, 'results_r14.jsonl'), 'w') as f:
        for win, (L, W) in WIN.items():
            M = C.primes_mask(L, W)
            M2 = C.primes_mask(2 * L, 2 * W + 2)
            P = np.nonzero(M)[0].astype(np.int64) + L
            sg = P[M2[2 * P + 1 - 2 * L]]
            sg = sg[sg > 5]
            for b in range(-1, 10):
                lo, hi = (L, L + W) if b < 0 else (L + b * (W // 10), L + (b + 1) * (W // 10))
                T = sg[(sg >= lo) & (sg < hi)]
                hgap = np.diff(T)
                st = {'n_sg': int(len(T)), 'n_gaps': int(len(hgap))}
                for x in W1:
                    st['W1.gap%d' % x] = float(np.mean(hgap == x))
                small = hgap <= 300
                st['W2.le300'] = float(np.mean(small))
                for q in (5, 7):
                    st['W3.same%d' % q] = float(np.mean(hgap[small] % q == 0))
                f.write(json.dumps({'window': win, 'block': b, 'stats': st}) + '\n')
            print(win, len(sg), flush=True)


if __name__ == '__main__':
    {'data': data, 'build': build, 'theory': theory}[sys.argv[1]]()
