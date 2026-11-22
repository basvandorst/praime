"""Round 29 (preregister/round29.txt): the law (order 3 + overlap correction) against the counts of gaps 460-700 below 4e18.
usage: r29.py build g1 g2 ... | compare"""
import os, sys, json, math, pickle, time
import numpy as np
import theory2 as T2
import common as C
from r26 import LAMG

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
CACHE = os.environ.get('STUDY_CACHE29', '/tmp/study_r29.pkl')
GS = [460, 482, 500, 520, 538, 556, 580, 610, 640, 670, 700]

if sys.argv[1] == 'build':
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    for g in [int(x) for x in sys.argv[2:]]:
        if g in cache:
            continue
        t0 = time.time()
        pt = T2.Pattern([0, g], 3)
        cache[g] = (pt.SB, np.array([pt.log_avoid(l, 3) + T2.overlap_correction(pt, l) for l in LAMG]))
        pickle.dump(cache, open(CACHE, 'wb'))
        print(g, 'live', len(pt.s1), round(time.time() - t0), 's', flush=True)
else:
    cache = pickle.load(open(CACHE, 'rb'))
    data = json.load(open(os.path.join(OUT, 'oes_gap_counts_4e18.json')))['N']
    M = C.primes_mask(2, 10 ** 7 - 2); P = np.append(np.nonzero(M)[0] + 2, 10000019); g0 = np.diff(P)
    u = np.linspace(math.log(1e7), math.log(4e18), 20001)
    out = {}
    for g in GS:
        SB, la = cache[g]
        sh = SB / u * np.exp(np.interp(u, LAMG, la))
        law = float(np.trapezoid(sh * np.exp(u) / u, u)) + int(np.sum(g0 == g))
        d = data[str(g)]
        out[g] = {'data': d, 'law': law, 'rel': (d - law) / d, 'g_over_lnx': g / math.log(4e18)}
        print(g, d, '%.4g' % law, '%+.3f' % out[g]['rel'], 'g/lnx %.1f' % out[g]['g_over_lnx'], flush=True)
    json.dump(out, open(os.path.join(OUT, 'r29.json'), 'w'), indent=1)
