"""Round 26 (preregister/round26.txt): the cluster law (+ polymer overlap correction) against the counts of the gaps
72..454 below 4e18. usage: r26.py build | compare"""
import os, sys, json, math, pickle, time
import numpy as np
import theory2 as T2
import common as C

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
CACHE = os.environ.get('STUDY_CACHE26', '/tmp/study_r26.pkl')
LAMG = np.linspace(math.log(1e7) - 0.01, math.log(4e18) + 0.01, 40)
GS = list(range(72, 455, 2))


def build():
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    t0 = time.time()
    for g in GS:
        if g in cache:
            continue
        o = 4 if g <= 120 else 3
        pt = T2.Pattern([0, g], o)
        if pt.SB == 0:
            cache[g] = (0.0, None); continue
        cache[g] = (pt.SB, np.array([pt.log_avoid(l, o) + T2.overlap_correction(pt, l) for l in LAMG]))
        if time.time() - t0 > 450:
            break
    pickle.dump(cache, open(CACHE, 'wb'))
    print('cached', sum(g in cache for g in GS), 'of', len(GS), flush=True)


def compare():
    cache = pickle.load(open(CACHE, 'rb'))
    data = json.load(open(os.path.join(OUT, 'oes_gap_counts_4e18.json')))['N']
    M = C.primes_mask(2, 10 ** 7 - 2); P = np.append(np.nonzero(M)[0] + 2, 10000019); g0 = np.diff(P)
    u = np.linspace(math.log(1e7), math.log(4e18), 20001)
    out = {}
    for g in GS:
        SB, la = cache[g]
        if SB == 0:
            continue
        sh = SB / u * np.exp(np.interp(u, LAMG, la))
        law = float(np.trapezoid(sh * np.exp(u) / u, u)) + int(np.sum(g0 == g))
        d = data[str(g)]
        out[g] = {'data': d, 'law': law, 'rel': (d - law) / d}
    for g in (72, 100, 150, 200, 250, 300, 350, 400, 454):
        if g in out:
            print(g, out[g]['data'], '%.4g' % out[g]['law'], '%+.4f' % out[g]['rel'], flush=True)
    json.dump(out, open(os.path.join(OUT, 'r26.json'), 'w'), indent=1)


if __name__ == '__main__':
    {'build': build, 'compare': compare}[sys.argv[1]]()
