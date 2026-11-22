"""Round 8 (preregister/round8.txt): the fifth-order Bahadur term for single gaps with span <= 120, order 4 for
spans 122..200, and the gap tail at train recomputed.

kappa_k(A) = sum over set partitions pi of A of (-1)^(|pi|-1) (|pi|-1)! prod_{blocks} mu(block); every partition
of a k-set carries lam^-k, so kappa_k = K_k / lam^k with K_k built from s_X = S(B u X)/S(B) alone.
The connected Bahadur term of order k is (-1)^k kappa_k / prod (1 - p).
usage: r8.py build | tail
"""
import sys, os, math, json, pickle, time
from itertools import combinations
import numpy as np
import theory2 as T2
import theory_r3 as R3

OUT = R3.OUT
CACHE8 = os.environ.get('STUDY_CACHE8', '/tmp/study_r8.pkl')
LAMT = np.linspace(math.log(1e7) - 0.01, math.log(1e8) + 0.01, 7)      # train range only


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


PARTS5 = list(set_partitions(list(range(5))))
assert len(PARTS5) == 52


def tables(B):
    """live interior and dense s-tables up to order 4 (s_X = S(B u X)/S(B))."""
    s = B[-1]
    SB = T2.S1(B)
    I = np.array([m for m in range(1, s) if m not in B], dtype=np.int64)
    Barr = np.array(B, dtype=np.int64)
    Ssub = lambda cols: T2.S_many(np.concatenate([np.tile(Barr, (cols.shape[0], 1)), cols], 1), P=max(s, 2)) / SB
    s1 = Ssub(I[:, None]); live = I[s1 > 0]; s1 = s1[s1 > 0]; n = len(live)
    i, j = np.triu_indices(n, 1)
    s2 = np.zeros((n, n)); v = Ssub(np.stack([live[i], live[j]], 1)); s2[i, j] = v; s2[j, i] = v
    T = np.array(list(combinations(range(n), 3)), dtype=np.int64)
    v3 = Ssub(live[T])
    s3 = np.zeros((n, n, n))
    from itertools import permutations
    for perm in permutations(range(3)):
        s3[T[:, perm[0]], T[:, perm[1]], T[:, perm[2]]] = v3
    return SB, live, s1, s2, s3, Ssub, n


def order4_term(live, s1, s2, s3, Ssub, n, lams):
    Q = np.array(list(combinations(range(n), 4)), dtype=np.int64)
    s4v = Ssub(live[Q])
    a, b, c, d = Q.T
    K4 = (s4v - s3[a, b, c] * s1[d] - s3[a, b, d] * s1[c] - s3[a, c, d] * s1[b] - s3[b, c, d] * s1[a]
          - s2[a, b] * s2[c, d] - s2[a, c] * s2[b, d] - s2[a, d] * s2[b, c]
          + 2 * (s2[a, b] * s1[c] * s1[d] + s2[a, c] * s1[b] * s1[d] + s2[a, d] * s1[b] * s1[c]
                 + s2[b, c] * s1[a] * s1[d] + s2[b, d] * s1[a] * s1[c] + s2[c, d] * s1[a] * s1[b])
          - 6 * s1[a] * s1[b] * s1[c] * s1[d])
    out = []
    for lam in lams:
        q = 1 - s1 / lam
        out.append(float(np.sum(K4 / lam ** 4 / (q[a] * q[b] * q[c] * q[d]))))
    return np.array(out), s4v, Q


def order5_term(live, s1, s2, s3, s4v, Q4, Ssub, n, lams):
    from itertools import permutations
    s4 = np.zeros((n,) * 4, dtype=np.float64)
    for perm in permutations(range(4)):
        s4[tuple(Q4[:, k] for k in perm)] = s4v
    sfun = {1: lambda idx: s1[idx[0]], 2: lambda idx: s2[idx[0], idx[1]], 3: lambda idx: s3[idx[0], idx[1], idx[2]],
            4: lambda idx: s4[idx[0], idx[1], idx[2], idx[3]]}
    tot = np.zeros(len(lams))
    combs = combinations(range(n), 5)
    while True:
        chunk = np.array(list(next(combs) for _ in range(400000)) if False else [x for _, x in zip(range(400000), combs)],
                         dtype=np.int64)
        if len(chunk) == 0:
            break
        s5 = Ssub(live[chunk])
        cols = [chunk[:, k] for k in range(5)]
        K5 = np.zeros(len(chunk))
        for part in PARTS5:
            m = len(part)
            coef = (-1) ** (m - 1) * math.factorial(m - 1)
            prod = np.ones(len(chunk))
            for blk in part:
                if len(blk) == 5:
                    prod = prod * s5
                else:
                    prod = prod * sfun[len(blk)]([cols[k] for k in blk])
            K5 += coef * prod
        for li, lam in enumerate(lams):
            q = 1 - s1 / lam
            den = q[cols[0]] * q[cols[1]] * q[cols[2]] * q[cols[3]] * q[cols[4]]
            tot[li] += float(np.sum(K5 / lam ** 5 / den))
    return -tot                                       # sign (-1)^5


def build():
    cache = pickle.load(open(CACHE8, 'rb')) if os.path.exists(CACHE8) else {}
    t0 = time.time()
    for h in list(range(50, 121, 2)) + list(range(122, 201, 2)):
        if h in cache:
            continue
        SB, live, s1, s2, s3, Ssub, n = tables([0, h])
        if SB == 0:
            cache[h] = None; continue
        t4, s4v, Q4 = order4_term(live, s1, s2, s3, Ssub, n, LAMT)
        t5 = order5_term(live, s1, s2, s3, s4v, Q4, Ssub, n, LAMT) if h <= 120 else np.zeros(len(LAMT))
        cache[h] = (t4, t5)
        print(h, n, np.round(t4, 5), np.round(t5, 5), round(time.time() - t0), flush=True)
        pickle.dump(cache, open(CACHE8, 'wb'))
        if time.time() - t0 > 480:
            break
    print('done', sum(h in cache for h in list(range(50, 121, 2)) + list(range(122, 201, 2))), 'of 76')


def tail_shares(c3, c8, L, W, variant):
    u = np.linspace(math.log(L), math.log(L + W), 2001); w = np.exp(u) / u
    g = {}
    for h in range(2, R3.HMAX + 1, 2):
        SB, la = c3[(0, h)]
        if SB == 0:
            g[h] = 0.0; continue
        lav = np.interp(u, R3.LAMG, la)
        if variant == 'round8' and h in c8 and c8[h] is not None:
            t4, t5 = c8[h]
            lav = lav + np.interp(u, LAMT, t4 if h > 120 else t5)
        g[h] = float(np.trapezoid(w * SB / u * np.exp(lav), u) / np.trapezoid(w, u))
    tot = sum(g.values())
    return {('T1.tail%d' % t): sum(v for h, v in g.items() if h >= t) / tot for t in (54, 72, 90)}


def blocks():
    from scipy.stats import t as st
    c3 = pickle.load(open(R3.CACHE, 'rb')); c8 = pickle.load(open(CACHE8, 'rb'))
    L, W = R3.WIN['train']
    rows = [json.loads(l) for l in open(os.path.join(OUT, 'results_r3.jsonl'))]
    bl = sorted([r for r in rows if r['window'] == 'train' and r['block'] >= 0], key=lambda r: r['block'])
    out = {}
    for variant in ('round3', 'round8'):
        th = [tail_shares(c3, c8, L + b * (W // 10), W // 10, variant) for b in range(10)]
        out[variant] = {}
        for k in th[0]:
            r = np.array([bl[b]['stats'][k] - th[b][k] for b in range(10)])
            z = r.mean() / (r.std(ddof=1) / math.sqrt(10))
            out[variant][k] = {'t9': float(z), 'p': float(2 * st.sf(abs(z), 9))}
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(OUT, 'r8_blocks.json'), 'w'), indent=1)


def tail():
    """train tail shares: round-3 law, with order 5 added for spans 50..120 and order 4 replacing order 3 for 122..200."""
    c3 = pickle.load(open(R3.CACHE, 'rb')); c8 = pickle.load(open(CACHE8, 'rb'))
    L, W = R3.WIN['train']
    u = np.linspace(math.log(L), math.log(L + W), 2001); w = np.exp(u) / u
    res = {}
    for variant in ('round3', 'round8'):
        g = {}
        for h in range(2, R3.HMAX + 1, 2):
            SB, la = c3[(0, h)]
            if SB == 0:
                g[h] = 0.0; continue
            lav = np.interp(u, R3.LAMG, la)
            if variant == 'round8' and h in c8 and c8[h] is not None:
                t4, t5 = c8[h]
                if h > 120:
                    lav = lav + np.interp(u, LAMT, t4)          # round 3 had order 3 here: add the order-4 term
                else:
                    lav = lav + np.interp(u, LAMT, t5)          # round 3 had order 4 here: add the order-5 term
            f = SB / u * np.exp(lav)
            g[h] = float(np.trapezoid(w * f, u) / np.trapezoid(w, u))
        tot = sum(g.values())
        res[variant] = {('tail%d' % t): sum(v for h, v in g.items() if h >= t) / tot for t in (54, 72, 90)}
        res[variant]['norm'] = tot
    v3 = json.load(open(os.path.join(OUT, 'verdicts_r3.json')))['train']
    res['primes'] = {('tail%d' % t): v3['T1.tail%d' % t]['x'] for t in (54, 72, 90)}
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(OUT, 'r8_tail.json'), 'w'), indent=1)


if __name__ == '__main__':
    {'build': build, 'tail': tail, 'blocks': blocks}[sys.argv[1]]()
