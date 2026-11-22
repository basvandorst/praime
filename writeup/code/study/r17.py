"""Round 17 (preregister/round17.txt): the frozen cluster law at 1e18, primes by gmpy2.next_prime (no sieve).
usage: r17.py data | build | test"""
import sys, os, math, json, pickle, time
import numpy as np
import gmpy2
from scipy.stats import t as student_t
import theory2 as T2

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
CACHE = os.environ.get('STUDY_CACHE17', '/tmp/study_r17.pkl')
L, W, K = 10 ** 18, 2 * 10 ** 6, 10
LAM = math.log(L + W / 2)
HMAX = int(12 * LAM + 40) // 2 * 2
GAPS = [2, 4, 6, 8, 10, 12, 18, 30]
QS = [3, 4, 5, 7, 8]
TH = [int(round(k * LAM / 2)) * 2 for k in (3, 4, 5)]


def order_for(span):
    return 4 if span <= 120 else (3 if span <= 300 else 2)


def data():
    ps = []
    p = gmpy2.next_prime(L - 1)
    while p < L + W:
        ps.append(int(p)); p = gmpy2.next_prime(p)
    ps.append(int(p))                      # the first prime beyond the window closes the last gap
    np.save(os.path.join(os.environ.get('STUDY_SCRATCH', '/tmp'), 'r17_primes.npy'), np.array(ps, dtype=object))
    print(len(ps) - 1, 'primes in the window', flush=True)


def build():
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    t0 = time.time()
    pats = [(0, h) for h in range(2, HMAX + 1, 2)] + [(0, h, 2 * h) for h in range(2, 151, 2)]
    for B in sorted(pats, key=lambda b: b[-1]):
        if B in cache:
            continue
        o = order_for(B[-1])
        pt = T2.Pattern(list(B), o)
        cache[B] = pt.SB / LAM ** (len(B) - 1) * math.exp(pt.log_avoid(LAM, o)) if pt.SB > 0 else 0.0
        if time.time() - t0 > 450:
            break
    pickle.dump(cache, open(CACHE, 'wb'))
    print('cached', sum(B in cache for B in pats), 'of', len(pats), flush=True)


def H_mm(counts):
    c = counts[counts > 0]; n = c.sum(); p = c / n
    return float(-np.sum(p * np.log2(p)) + (len(c) - 1) / (2 * n * math.log(2)))


def stats(P):
    g = np.diff(P)
    d = {}
    for h in GAPS:
        d['S.gap%d' % h] = float(np.mean(g == h))
    d['N.eqgap'] = float(np.mean(g[1:] == g[:-1]))
    for q in QS:
        r = P % q
        d['A.same%d' % q] = float(np.mean(r[1:-1] == r[:-2])) if False else float(np.mean(r[1:] == r[:-1]))
    cl = np.minimum(g // 2, 20)
    d['T4.H_class'] = H_mm(np.bincount(cl, minlength=21)[1:].astype(float))
    for t in TH:
        d['T1.tail%d' % t] = float(np.mean(g >= t))
    return d


def test():
    cache = pickle.load(open(CACHE, 'rb'))
    gl = {h: cache[(0, h)] for h in range(2, HMAX + 1, 2)}
    tot = sum(gl.values())
    T = {'norm': tot}
    for h in GAPS:
        T['S.gap%d' % h] = gl[h]
    T['N.eqgap'] = sum(cache[(0, h, 2 * h)] for h in range(2, 151, 2))
    for q in QS:
        T['A.same%d' % q] = sum(v for h, v in gl.items() if h % q == 0) / tot
    cl = np.zeros(20)
    for h, v in gl.items():
        cl[min(h // 2, 20) - 1] += v / tot
    T['T4.H_class'] = float(-np.sum(cl * np.log2(cl)))
    for t in TH:
        T['T1.tail%d' % t] = sum(v for h, v in gl.items() if h >= t) / tot
    P = np.array([int(x) - L for x in np.load(os.path.join(os.environ.get('STUDY_SCRATCH', '/tmp'), 'r17_primes.npy'),
                                              allow_pickle=True)], dtype=np.int64)
    whole = stats(P)
    blocks = []
    for b in range(K):
        sel = (P >= b * W // K) & (P < (b + 1) * W // K)
        idx = np.nonzero(sel)[0]
        blocks.append(stats(P[idx[0]:idx[-1] + 2]))        # include the gap after the last prime of the block
    # estimator bias of the class entropy at the block size, from the law
    rng = np.random.default_rng(17)
    nb = int(np.mean([np.sum((P >= b * W // K) & (P < (b + 1) * W // K)) for b in range(K)]))
    bias = float(np.mean([H_mm(rng.multinomial(nb, cl).astype(float)) for _ in range(200)]) - T['T4.H_class'])
    res = {'primes': int(len(P) - 1), 'lam': LAM, 'norm_law': tot, 'H_class_bias_block': bias, 'stats': {}}
    n = len(P) - 1
    for s in whole:
        r = np.array([bl[s] - T[s] - (bias if s == 'T4.H_class' else 0.0) for bl in blocks])
        z = float(r.mean() / (r.std(ddof=1) / math.sqrt(K)))
        zb = (whole[s] - T[s]) / math.sqrt(T[s] * (1 - T[s]) / n) if not s.startswith('T4') else None
        res['stats'][s] = {'x': whole[s], 'T': T[s], 'rel': (whole[s] - T[s]) / T[s], 't9': z,
                           'p': float(2 * student_t.sf(abs(z), K - 1)), 'z_binom': zb}
        print('%-12s %.6g %.6g %+.4f t9 %.2f zb %s' % (s, whole[s], T[s], (whole[s] - T[s]) / T[s], z,
                                                     '%.2f' % zb if zb is not None else ''), flush=True)
    json.dump(res, open(os.path.join(OUT, 'r17.json'), 'w'), indent=1)


if __name__ == '__main__':
    {'data': data, 'build': build, 'test': test}[sys.argv[1]]()
