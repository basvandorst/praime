"""Round 38 (diagnostic, preregister/round38.txt): order 4 and order 5 of the Bahadur series (ln P terms) for long
gaps, on the lam grid of round 32. Order 5 sums the 52-partition cumulant over all 5-subsets of the live interior
points; the order-4 values are kept in a colex-indexed vector (C(n,4) floats) instead of a dense n^4 table, so gaps
with up to ~100 live points fit in memory.
usage: r38.py build g1,g2,...   (resumable; stops after ~8 minutes)   |   r38.py report   |   r38.py u
(u: mean U of the first occurrences of the gaps 180..300 with order 5 added; needs STUDY_CACHE3 and STUDY_CACHE32)"""
import os, sys, json, math, pickle, time
from itertools import combinations
from math import comb
import numpy as np
import r8

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
CACHE = os.environ.get('STUDY_CACHE38', '/tmp/study_r38.pkl')
LAMS = np.array([16.0, 18.0, 20.0, 22.0, 24.0, 27.0])
CH = 300000


def colex(cols):
    r = np.zeros(len(cols[0]), dtype=np.int64)
    for k, c in enumerate(cols):
        r += np.array([comb(int(v), k + 1) for v in range(c.max() + 1)], dtype=np.int64)[c]
    return r


def terms(g):
    SB, live, s1, s2, s3, Ssub, n = r8.tables([0, g])
    q = [1 - s1 / lam for lam in LAMS]
    # order 4, and the colex vector of S over 4-subsets
    s4c = np.zeros(comb(n, 4))
    t4 = np.zeros(len(LAMS))
    it = combinations(range(n), 4)
    while True:
        Q = np.array([x for _, x in zip(range(CH), it)], dtype=np.int64)
        if len(Q) == 0:
            break
        a, b, c, d = Q.T
        s4v = Ssub(live[Q])
        s4c[colex([a, b, c, d])] = s4v
        K4 = (s4v - s3[a, b, c] * s1[d] - s3[a, b, d] * s1[c] - s3[a, c, d] * s1[b] - s3[b, c, d] * s1[a]
              - s2[a, b] * s2[c, d] - s2[a, c] * s2[b, d] - s2[a, d] * s2[b, c]
              + 2 * (s2[a, b] * s1[c] * s1[d] + s2[a, c] * s1[b] * s1[d] + s2[a, d] * s1[b] * s1[c]
                     + s2[b, c] * s1[a] * s1[d] + s2[b, d] * s1[a] * s1[c] + s2[c, d] * s1[a] * s1[b])
              - 6 * s1[a] * s1[b] * s1[c] * s1[d])
        for i, lam in enumerate(LAMS):
            t4[i] += float(np.sum(K4 / lam ** 4 / (q[i][a] * q[i][b] * q[i][c] * q[i][d])))
    # order 5
    look = {1: lambda x: s1[x[0]], 2: lambda x: s2[x[0], x[1]], 3: lambda x: s3[x[0], x[1], x[2]],
            4: lambda x: s4c[colex(x)]}
    t5 = np.zeros(len(LAMS))
    it = combinations(range(n), 5)
    while True:
        Q = np.array([x for _, x in zip(range(CH), it)], dtype=np.int64)
        if len(Q) == 0:
            break
        cols = [Q[:, k] for k in range(5)]
        s5 = Ssub(live[Q])
        K5 = np.zeros(len(Q))
        for part in r8.PARTS5:
            m = len(part)
            prod = np.full(len(Q), (-1.0) ** (m - 1) * math.factorial(m - 1))
            for blk in part:
                blk = sorted(blk)
                prod *= s5 if len(blk) == 5 else look[len(blk)]([cols[k] for k in blk])
            K5 += prod
        for i, lam in enumerate(LAMS):
            den = q[i][cols[0]] * q[i][cols[1]] * q[i][cols[2]] * q[i][cols[3]] * q[i][cols[4]]
            t5[i] -= float(np.sum(K5 / lam ** 5 / den))        # sign (-1)^5, as r8.order5_term
    return n, t4, t5


def klass(g):
    return 'six' if g % 6 == 0 else ('ten' if g % 10 == 0 else 'other')


def u_mode(cache):
    from scipy.stats import kstest
    import theory2 as T2
    import theory_r3 as R3
    c3 = pickle.load(open(R3.CACHE, 'rb')); c4 = pickle.load(open(os.environ['STUDY_CACHE32'], 'rb'))
    GS = list(range(180, 301, 2))
    P = {int(g): v for g, v in json.load(open(os.path.join(OUT, 'oes_first_occurrence.json')))['P'].items() if int(g) in GS}
    ratio = {g: v[2] / v[1] for g, v in cache.items()}
    fit = {}
    for k in ('six', 'ten', 'other'):
        gs = np.array([g for g in ratio if klass(g) == k], dtype=float)
        R = np.array([ratio[int(g)] for g in gs])
        fit[k] = [np.polyfit(gs, R[:, i], 1) for i in range(len(LAMS))]
    res = {'g': [], 'class': [], 'order5': [], 'U3': [], 'U4': [], 'U5': [], 'U5_nearest': [], 'U_geometric': []}
    for g in GS:
        SB, la = c3[(0, g)]
        pt = T2.Pattern([0, g], 2)
        oc = np.array([T2.overlap_correction(pt, l) for l in R3.LAMG])
        u = np.linspace(math.log(1e7), math.log(P[g]), 2001)
        base = np.interp(u, R3.LAMG, la + oc)
        t4g = np.asarray(c4[g])
        if g in cache:
            r = ratio[g]; rn = r; src = 'computed'
        else:
            r = np.array([np.polyval(fit[klass(g)][i], g) for i in range(len(LAMS))])
            same = [h for h in ratio if klass(h) == klass(g)]
            rn = ratio[min(same, key=lambda h: abs(h - g))]; src = 'trend'
        t4 = np.interp(u, LAMS, t4g); t5 = np.interp(u, LAMS, r * t4g); t5n = np.interp(u, LAMS, rn * t4g)
        rr = np.interp(u, LAMS, r)
        I = lambda lav: float(np.trapezoid(SB / u * np.exp(lav) * np.exp(u) / u, u))
        res['g'].append(g); res['class'].append(klass(g)); res['order5'].append(src)
        res['U3'].append(I(base)); res['U4'].append(I(base + t4)); res['U5'].append(I(base + t4 + t5))
        res['U5_nearest'].append(I(base + t4 + t5n)); res['U_geometric'].append(I(base + t4 + t5 / (1 - rr)))
    g = np.array(res['g'])
    for key in ('U3', 'U4', 'U5', 'U5_nearest', 'U_geometric'):
        v = np.array(res[key]); ks = kstest(v, 'expon')
        summ = {'n': len(v), 'mean': float(v.mean()), 'se': float(v.std(ddof=1) / math.sqrt(len(v))), 'ks_p': float(ks.pvalue)}
        for name, m in (('six', g % 6 == 0), ('not_six', g % 6 != 0)):
            summ['mean_' + name] = float(v[m].mean()); summ['se_' + name] = float(v[m].std(ddof=1) / math.sqrt(m.sum()))
        res[key + '_summary'] = summ
        print(key, {k: round(x, 3) for k, x in summ.items()}, flush=True)
    res['n_computed'] = int(sum(s == 'computed' for s in res['order5']))
    json.dump(res, open(os.path.join(OUT, 'r38_u.json'), 'w'), indent=1)


if __name__ == '__main__':
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    if sys.argv[1] == 'build':
        t0 = time.time()
        for g in map(int, sys.argv[2].split(',')):
            if g in cache:
                continue
            if time.time() - t0 > 240:
                break
            t1 = time.time()
            cache[g] = terms(g)
            pickle.dump(cache, open(CACHE, 'wb'))
            print(g, 'live', cache[g][0], 'secs', round(time.time() - t1),
                  'ratio', np.round(cache[g][2] / cache[g][1], 3), flush=True)
    elif sys.argv[1] == 'u':
        u_mode(cache)
    else:
        out = {str(g): {'live': int(v[0]), 'lam': LAMS.tolist(), 'order4_term_ln': v[1].tolist(),
                        'order5_term_ln': v[2].tolist(), 'ratio5_4': (v[2] / v[1]).tolist()} for g, v in sorted(cache.items())}
        json.dump(out, open(os.path.join(OUT, 'r38.json'), 'w'), indent=1)
        for g, d in out.items():
            print(g, d['live'], np.round(d['ratio5_4'], 3))
