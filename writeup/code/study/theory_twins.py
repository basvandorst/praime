"""Round 4: the law of consecutive twin primes from the singular series (preregister/round4.txt).

EventPattern: a boundary set B and interior events e (sets of offsets). X_e = 1[all of e prime]. Hardy-Littlewood:
    E prod_{e in A} X_e  (given B prime) = S(B u U) / (S(B) lam^{|U minus B|}),  U = union of the events in A.
The Bahadur form of round 2 then gives ln P(no event | B) = sum ln(1-p_e) + sum c~ - sum k3~ (orders 1..3).

usage: theory_twins.py build | values
"""
import sys, os, math, json, pickle, time
from itertools import combinations
import numpy as np
import theory2 as T2
import common as C

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'data', 'study')
CACHE = os.environ.get('STUDY_CACHE4', '/tmp/study_twins.pkl')
WIN = dict(C.WINDOWS); WIN['C'] = (10 ** 15, 2 * 10 ** 7)
LAMG = np.linspace(math.log(1e7) - 0.01, math.log(1e15 + 2e7) + 0.01, 60)
HS = list(range(6, 301, 6))
W1 = [6, 12, 18, 24, 30, 36, 42, 48, 60, 90, 120, 150, 180, 210]


def S_sets(sets, P):
    """singular series of a list of python sets (any sizes), grouped by size."""
    out = np.zeros(len(sets))
    by = {}
    for i, s in enumerate(sets):
        by.setdefault(len(s), []).append(i)
    for k, idx in by.items():
        arr = np.array([sorted(sets[i]) for i in idx], dtype=np.int64)
        out[idx] = T2.S_many(arr, P=P)
    return out


class EventPattern:
    def __init__(self, B, events, order=3):
        self.B = frozenset(B)
        span = max(max(B), max((max(e) for e in events), default=0)) - min(B)
        P = max(span, 2)
        self.SB = T2.S1(sorted(B))
        U1 = [self.B | frozenset(e) for e in events]
        s1 = S_sets(U1, P) / self.SB
        live = [i for i in range(len(events)) if s1[i] > 0]
        self.s1 = s1[live]
        ev = [frozenset(events[i]) for i in live]
        self.n1 = np.array([len(e - self.B) for e in ev])
        n = len(ev)
        self.pairs = self.trip = None
        if order >= 2 and n >= 2:
            i, j = np.triu_indices(n, 1)
            U = [self.B | ev[a] | ev[b] for a, b in zip(i, j)]
            s2v = S_sets(U, P) / self.SB
            n2v = np.array([len(u - self.B) for u in U])
            self.s2 = np.zeros((n, n)); self.s2[i, j] = s2v; self.s2[j, i] = s2v
            self.n2 = np.zeros((n, n), dtype=int); self.n2[i, j] = n2v; self.n2[j, i] = n2v
            self.pairs = (i, j)
        if order >= 3 and n >= 3:
            T = np.array(list(combinations(range(n), 3)), dtype=np.int64)
            U = [self.B | ev[a] | ev[b] | ev[c] for a, b, c in T]
            self.trip = (T[:, 0], T[:, 1], T[:, 2], S_sets(U, P) / self.SB,
                         np.array([len(u - self.B) for u in U]))

    def log_avoid(self, lam, order):
        p = self.s1 / lam ** self.n1
        out = float(np.sum(np.log1p(-p)))
        if order >= 2 and self.pairs is not None:
            m2 = self.s2 / lam ** self.n2
            i, j = self.pairs
            out += float(np.sum((m2[i, j] - p[i] * p[j]) / ((1 - p[i]) * (1 - p[j]))))
            if order >= 3 and self.trip is not None:
                a, b, c, s3, n3 = self.trip
                k3 = s3 / lam ** n3 - m2[a, b] * p[c] - m2[a, c] * p[b] - m2[b, c] * p[a] + 2 * p[a] * p[b] * p[c]
                out -= float(np.sum(k3 / ((1 - p[a]) * (1 - p[b]) * (1 - p[c]))))
        return out


def twin_pattern(h, order):
    return EventPattern([0, 2, h, h + 2], [(m, m + 2) for m in range(1, h)], order)


def build():
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    t0 = time.time()
    for h in HS:
        if h in cache:
            continue
        pt = twin_pattern(h, 3)
        cache[h] = (pt.SB, [np.array([pt.log_avoid(l, o) for l in LAMG]) for o in (1, 2, 3)])
        if time.time() - t0 > 450:
            break
    pickle.dump(cache, open(CACHE, 'wb'))
    print('cached', len(cache), 'of', len(HS), flush=True)


def share(cache, h, L, W, order=3):
    """share per twin pair: P(B) / (S(0,2)/lam^2) = S(B)/S(0,2) / lam^2 * avoid, weight 1/ln^2 t."""
    SB, las = cache[h]
    la = las[order - 1]
    S02 = T2.S1([0, 2])
    u = np.linspace(math.log(L), math.log(L + W), 2001)
    w = np.exp(u) / u ** 2
    f = SB / S02 / u ** 2 * np.exp(np.interp(u, LAMG, la))
    return float(np.trapezoid(w * f, u) / np.trapezoid(w, u))


def values(L, W, cache, order=3):
    g = {h: share(cache, h, L, W, order) for h in HS}
    t = {('W1.gap%d' % h): g[h] for h in W1}
    tot = sum(g.values())
    t['W2.le300'] = tot
    for q in (5, 7):
        t['W3.same%d' % q] = sum(v for h, v in g.items() if h % q == 0) / tot
    return t


def main():
    if sys.argv[1] == 'build':
        build(); return
    cache = pickle.load(open(CACHE, 'rb'))
    assert all(h in cache for h in HS)
    res = {}
    for win, (L, W) in WIN.items():
        res[win] = {'-1': values(L, W, cache)}
        res[win]['order2'] = values(L, W, cache, 2)
        res[win]['order1'] = values(L, W, cache, 1)
        for b in range(10):
            res[win][str(b)] = values(L + b * (W // 10), W // 10, cache)
        print(win, {k: round(v, 5) for k, v in res[win]['-1'].items()}, flush=True)
    json.dump(res, open(os.path.join(OUT, 'theory_twins.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
