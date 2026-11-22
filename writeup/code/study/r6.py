"""Round 6 (preregister/round6.txt): second differences of consecutive primes, and bits per prime.

usage: r6.py build | theory | stats | bits
Uses the round-3 pattern cache (theory_r3.py: same lam grid, same frozen orders) and adds the patterns (0, h, 2h + k).
"""
import sys, os, math, json, pickle, time
import numpy as np
import theory2 as T2
import theory_r3 as R3
import common as C
import probes as PR

OUT = R3.OUT
KS = [-6, -2, 0, 2, 6]
HM = 150


def pats():
    return [(0, h, 2 * h + k) for k in KS for h in range(2, HM + 1, 2) if h + k >= 2]


def build():
    cache = pickle.load(open(R3.CACHE, 'rb'))
    t0 = time.time()
    todo = sorted([B for B in pats() if B not in cache], key=lambda B: B[-1])
    for B in todo:
        o = R3.order_for(B[-1])
        pt = T2.Pattern(list(B), o)
        cache[B] = (pt.SB, np.array([pt.log_avoid(l, o) for l in R3.LAMG]) if pt.SB > 0 else None)
        if time.time() - t0 > 450:
            break
    pickle.dump(cache, open(R3.CACHE, 'wb'))
    print('cached', sum(B in cache for B in pats()), 'of', len(pats()), flush=True)


def tvals(L, W, cache):
    t = {}
    for k in KS:
        t['N.d2_%+d' % k] = sum(R3.share(cache, (0, h, 2 * h + k), L, W) for h in range(2, HM + 1, 2) if h + k >= 2)
    return t


def theory():
    cache = pickle.load(open(R3.CACHE, 'rb'))
    res = {}
    for win, (L, W) in R3.WIN.items():
        res[win] = {'-1': tvals(L, W, cache)}
        for b in range(10):
            res[win][str(b)] = tvals(L + b * (W // 10), W // 10, cache)
        print(win, {k: round(v, 6) for k, v in res[win]['-1'].items()}, flush=True)
    json.dump(res, open(os.path.join(OUT, 'theory_r6.json'), 'w'), indent=1)


def H(c):
    p = c / c.sum(); p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))


def sstats(M, L):
    P = np.nonzero(M)[0].astype(np.int64) + L
    g = np.diff(P)
    d = g[1:] - g[:-1]
    st = {'n_gaps': int(len(g))}
    for k in KS:
        st['N.d2_%+d' % k] = float(np.mean(d == k))
    c = np.minimum(g // 2, 20) - 1                       # 0..19
    pair = np.bincount(c[:-1] * 20 + c[1:], minlength=400).reshape(20, 20)
    trip = np.bincount((c[:-2] * 20 + c[1:-1]) * 20 + c[2:], minlength=8000).astype(float)
    st['_pair'] = pair.tolist()
    n2, n3 = pair.sum(), trip.sum()
    # conditional entropies with Miller-Madow: H(X|Y) = H(X,Y) - H(Y), each with +(K-1)/(2N ln 2)
    mm = lambda cnt: H(cnt) + (np.count_nonzero(cnt) - 1) / (2 * cnt.sum() * math.log(2))
    st['H1_mm'] = mm(pair.ravel().astype(float)) - mm(pair.sum(1).astype(float))
    st['H2_mm'] = mm(trip) - mm(pair.ravel().astype(float))
    return st


def stats():
    with open(os.path.join(OUT, 'results_r6.jsonl'), 'w') as f:
        for win, (L, W) in R3.WIN.items():
            M = C.primes_mask(L, W)
            for b in range(-1, 10):
                Lb, Wb = (L, W) if b < 0 else (L + b * (W // 10), W // 10)
                f.write(json.dumps({'window': win, 'block': b,
                                    'stats': sstats(M[Lb - L:Lb - L + Wb], Lb)}) + '\n')
            print(win, flush=True)


def class_law(L, W, weight):
    """class probabilities (c = 1..20) of the law P(h) proportional to weight(h, lam), window average."""
    hs = np.arange(2, R3.HMAX + 1, 2)
    u = np.linspace(math.log(L), math.log(L + W), 401)
    w = np.exp(u) / u
    out = np.zeros(20)
    for ui, wi in zip(u, w):
        p = np.array([weight(h, ui) for h in hs]); p /= p.sum()
        cl = np.zeros(20)
        np.add.at(cl, np.minimum(hs // 2, 20) - 1, p)
        out += wi * cl
    return out / out.sum()


def bits():
    th = json.load(open(os.path.join(OUT, 'theory_r3.json')))
    rows = [json.loads(l) for l in open(os.path.join(OUT, 'results_r6.jsonl'))]
    Sp = {h: T2.S1([0, h]) for h in range(2, R3.HMAX + 1, 2)}
    res = {}
    for win, (L, W) in R3.WIN.items():
        whole = [r for r in rows if r['window'] == win and r['block'] == -1][0]['stats']
        pair = np.array(whole['_pair'], float); N = pair.sum()
        nxt = pair.sum(0)                                    # counts of the predicted symbol
        J = np.clip(np.array(th[win]['-1']['_joint']), 0, None); J /= J.sum()
        M0 = class_law(L, W, lambda h, lam: math.exp(-h / lam))
        M1 = class_law(L, W, lambda h, lam: Sp[h] * math.exp(-h / lam))
        M2 = J.sum(0)
        M3 = J / J.sum(1, keepdims=True)
        ce = lambda q: float(-np.sum(nxt * np.log2(q)) / N)
        r = {'M0 coins': ce(M0), 'M1 crude Gallagher': ce(M1), 'M2 cluster law, i.i.d.': ce(M2),
             'M3 cluster law, Markov': float(-np.sum(np.where(pair > 0, pair * np.log2(np.where(M3 > 0, M3, 1)), 0)) / N),
             'measured H(c_n+1 | c_n)': whole['H1_mm'], 'measured H(c_n+1 | c_n, c_n-1)': whole['H2_mm']}
        res[win] = r
        print(win, {k: round(v, 4) for k, v in r.items()}, flush=True)
    json.dump(res, open(os.path.join(OUT, 'bits_r6.json'), 'w'), indent=1)


if __name__ == '__main__':
    {'build': build, 'theory': theory, 'stats': stats, 'bits': bits}[sys.argv[1]]()
