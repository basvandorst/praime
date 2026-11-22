"""Round 32 (preregister/round32.txt): U_g for gaps 180..300 with the order-4 term and the overlap correction.
usage: r32b.py build | test"""
import os, sys, json, math, pickle, time
from itertools import combinations
import numpy as np
from scipy.stats import kstest
import theory2 as T2
import theory_r3 as R3
import r8

OUT = R3.OUT
CACHE = os.environ.get('STUDY_CACHE32', '/tmp/study_r32.pkl')
LAMS = np.array([16.0, 18.0, 20.0, 22.0, 24.0, 27.0])
GS = list(range(180, 301, 2))


def order4(g):
    SB, live, s1, s2, s3, Ssub, n = r8.tables([0, g])
    tot = np.zeros(len(LAMS))
    it = combinations(range(n), 4)
    while True:
        Q = np.array([x for _, x in zip(range(500000), it)], dtype=np.int64)
        if len(Q) == 0:
            break
        s4v = Ssub(live[Q]); a, b, c, d = Q.T
        K4 = (s4v - s3[a, b, c] * s1[d] - s3[a, b, d] * s1[c] - s3[a, c, d] * s1[b] - s3[b, c, d] * s1[a]
              - s2[a, b] * s2[c, d] - s2[a, c] * s2[b, d] - s2[a, d] * s2[b, c]
              + 2 * (s2[a, b] * s1[c] * s1[d] + s2[a, c] * s1[b] * s1[d] + s2[a, d] * s1[b] * s1[c]
                     + s2[b, c] * s1[a] * s1[d] + s2[b, d] * s1[a] * s1[c] + s2[c, d] * s1[a] * s1[b])
              - 6 * s1[a] * s1[b] * s1[c] * s1[d])
        for li, lam in enumerate(LAMS):
            q = 1 - s1 / lam
            tot[li] += float(np.sum(K4 / lam ** 4 / (q[a] * q[b] * q[c] * q[d])))
    return tot


if sys.argv[1] == 'build':
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    t0 = time.time()
    for g in GS:
        if g in cache:
            continue
        cache[g] = order4(g)
        pickle.dump(cache, open(CACHE, 'wb'))
        if time.time() - t0 > 480:
            break
    print('cached', sum(g in cache for g in GS), 'of', len(GS), round(time.time() - t0), 's', flush=True)
else:
    c3 = pickle.load(open(R3.CACHE, 'rb')); c4 = pickle.load(open(CACHE, 'rb'))
    P = {int(g): v for g, v in json.load(open(os.path.join(OUT, 'oes_first_occurrence.json')))['P'].items() if int(g) in GS}
    res = {'g': [], 'U3': [], 'U4': []}
    for g in GS:
        SB, la = c3[(0, g)]
        pt = T2.Pattern([0, g], 2)
        oc = np.array([T2.overlap_correction(pt, l) for l in R3.LAMG])
        u = np.linspace(math.log(1e7), math.log(P[g]), 2001)
        base = np.interp(u, R3.LAMG, la + oc)
        t4 = np.interp(u, LAMS, c4[g])
        for key, lav in (('U3', base), ('U4', base + t4)):
            res[key].append(float(np.trapezoid(SB / u * np.exp(lav) * np.exp(u) / u, u)))
        res['g'].append(g)
    for key in ('U3', 'U4'):
        v = np.array(res[key]); ks = kstest(v, 'expon')
        res[key + '_summary'] = {'n': len(v), 'mean': float(v.mean()), 'se': float(v.std(ddof=1) / math.sqrt(len(v))),
                                 'median': float(np.median(v)), 'ks_p': float(ks.pvalue)}
        print(key, res[key + '_summary'], flush=True)
    json.dump(res, open(os.path.join(OUT, 'r32b.json'), 'w'), indent=1)
