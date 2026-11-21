"""Round 2 predictions for every window and block, from theory2.Pattern (HL + Bahadur cluster expansion).

usage: theory_r2.py build   -> computes ln(avoidance) of every pattern on a common lam grid, cache file
       theory_r2.py values  -> integrates over every window/block, writes data/study/theory_r2.json
Frozen choices (preregister/round2.txt): order 4 for span <= 120, order 3 above; single gaps to h = 12 lam + 40.
"""
import sys, os, math, json, pickle, time
import numpy as np
import theory2 as T2
import common as C
import probes as PR

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'data', 'study')
CACHE = os.environ.get('STUDY_CACHE', '/tmp/study_patterns_r2.pkl')
LAMG = np.linspace(math.log(1e7) - 0.01, math.log(1e12 + 5e7) + 0.01, 64)
HMAX = int(12 * LAMG[-1] + 40) // 2 * 2          # 370
JCAP = 38                                          # joint law exact for h1, h2 <= 38 (classes 1..19)


def order_for(span):
    return 4 if span <= 120 else 3


def patterns():
    pats = [(0, h) for h in range(2, HMAX + 1, 2)]
    pats += [(0, h, 2 * h) for h in range(2, HMAX // 2 + 1, 2)]
    pats += [(0, a, a + b) for a in range(2, JCAP + 1, 2) for b in range(2, JCAP + 1, 2) if a != b]
    return pats


def build():
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    t0 = time.time()
    for B in patterns():
        if B in cache:
            continue
        o = order_for(B[-1])
        pt = T2.Pattern(list(B), o)
        la = np.array([pt.log_avoid(l, o) for l in LAMG]) if pt.SB > 0 else None
        cache[B] = (pt.SB, la)
        if time.time() - t0 > 400:
            break
    pickle.dump(cache, open(CACHE, 'wb'))
    print('cached', len(cache), 'of', len(patterns()))


def share(cache, B, L, W):
    SB, la = cache[B]
    if SB == 0:
        return 0.0
    k = len(B)
    return T2.window_avg(L, W, lambda u: SB / u ** (k - 1) * np.exp(np.interp(u, LAMG, la)))


def mi_bits(J):
    J = J / J.sum()
    px, py = J.sum(1), J.sum(0)
    nz = J > 0
    return float(np.sum(J[nz] * np.log2(J[nz] / np.outer(px, py)[nz])))


def values(L, W, cache):
    t = {}
    g = {h: share(cache, (0, h), L, W) for h in range(2, HMAX + 1, 2)}
    tot = sum(g.values())
    t['norm_gap_law'] = tot                       # should be close to 1
    t['mean_gap_law'] = sum(h * v for h, v in g.items())
    for h in PR.GAPS:
        t['S.gap%d' % h] = g[h]
    t['N.eqgap'] = sum(share(cache, (0, h, 2 * h), L, W) for h in range(2, HMAX // 2 + 1, 2))
    for q in PR.SAME_Q:
        t['A.same%d' % q] = sum(v for h, v in g.items() if h % q == 0) / tot
    for H in PR.VAR_H:
        t['S.var%d_pairs' % H] = T2.var_ratio_exact(H, L, W)
        t['S.var%d' % H] = T2.var_ratio_aligned(H, L, W)     # round 2b: intervals on a grid aligned mod 10
    # joint law of two consecutive gaps -> classes c = min(h/2, 20); class 20 = all h >= 40 (from marginals)
    K = 20
    J = np.zeros((K + 1, K + 1))
    for a in range(2, JCAP + 1, 2):
        for b in range(2, JCAP + 1, 2):
            J[a // 2, b // 2] = share(cache, (0, a, a + b), L, W) if a != b else share(cache, (0, a, 2 * a), L, W)
    marg = np.zeros(K + 1)
    for h, v in g.items():
        marg[min(h // 2, K)] += v / tot
    J = J / tot                                   # same normalisation as the marginal
    for c in range(1, K):
        J[c, K] = marg[c] - J[c, 1:K].sum()
        J[K, c] = marg[c] - J[1:K, c].sum()
    J[K, K] = 1 - J[1:K + 1, 1:K + 1].sum() + J[K, K]
    t['joint_negative_cells'] = int((J < 0).sum())
    J = np.clip(J[1:, 1:], 0, None)
    t['I.migap_lag1'] = mi_bits(J)
    t['D.gain1_gap'] = t['I.migap_lag1']           # expected cross-validated gain = MI (minus a ~1e-4 penalty)
    return t


def main():
    if sys.argv[1] == 'build':
        build(); return
    cache = pickle.load(open(CACHE, 'rb'))
    missing = [B for B in patterns() if B not in cache]
    assert not missing, '%d patterns missing' % len(missing)
    res = {}
    for win, (L, W) in C.WINDOWS.items():
        K = 10
        res[win] = {'-1': values(L, W, cache)}
        for b in range(K):
            res[win][str(b)] = values(L + b * (W // K), W // K, cache)
        print(win, {k: round(v, 6) for k, v in res[win]['-1'].items()}, flush=True)
    json.dump(res, open(os.path.join(OUT, 'theory_r2.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
