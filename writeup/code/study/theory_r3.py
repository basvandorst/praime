"""Round 3 predictions (preregister/round3.txt): the frozen round-2 cluster law on a lam grid that reaches 1e15,
plus tails, four primes in a row in arithmetic progression (CPAP-4), entropies and the 361-cell joint law.

usage: theory_r3.py build | values
Rule changes against round 2, both on terms below 1e-7: the equal-gap sum stops at h = 150 (span 300), and the
CPAP-4 sum at h = 150; only h = 0 mod 6 gives CPAP-4 (otherwise {0,h,2h,3h} covers every class mod 3).
"""
import sys, os, math, json, pickle, time
import numpy as np
import theory2 as T2
import common as C
import probes as PR

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'data', 'study')
CACHE = os.environ.get('STUDY_CACHE3', '/tmp/study_patterns_r3.pkl')
WIN = dict(C.WINDOWS); WIN['C'] = (10 ** 15, 2 * 10 ** 7)
LAMG = np.linspace(math.log(1e7) - 0.01, math.log(1e15 + 2e7) + 0.01, 80)
HMAX = int(12 * LAMG[-1] + 40) // 2 * 2
JCAP = 38
EQMAX = 150


def order_for(span):
    return 4 if span <= 120 else 3


def patterns():
    pats = [(0, h) for h in range(2, HMAX + 1, 2)]
    pats += [(0, h, 2 * h) for h in range(2, EQMAX + 1, 2)]
    pats += [(0, a, a + b) for a in range(2, JCAP + 1, 2) for b in range(2, JCAP + 1, 2) if a != b]
    pats += [(0, h, 2 * h, 3 * h) for h in range(6, EQMAX + 1, 6)]
    return pats


def build():
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    t0 = time.time()
    todo = [B for B in patterns() if B not in cache]
    todo.sort(key=lambda B: B[-1])
    for B in todo:
        o = order_for(B[-1])
        pt = T2.Pattern(list(B), o)
        la = np.array([pt.log_avoid(l, o) for l in LAMG]) if pt.SB > 0 else None
        cache[B] = (pt.SB, la)
        if time.time() - t0 > 450:
            break
    pickle.dump(cache, open(CACHE, 'wb'))
    print('cached', len([B for B in patterns() if B in cache]), 'of', len(patterns()), flush=True)


def share(cache, B, L, W):
    SB, la = cache[B]
    if SB == 0:
        return 0.0
    k = len(B)
    return T2.window_avg(L, W, lambda u: SB / u ** (k - 1) * np.exp(np.interp(u, LAMG, la)))


def tail_thresholds(win):
    L, W = WIN[win]
    lm = math.log(L + W / 2)
    return [int(round(k * lm / 2)) * 2 for k in (3, 4, 5)]


def joint(cache, g, tot, L, W):
    K = 20
    J = np.zeros((K + 1, K + 1))
    for a in range(2, JCAP + 1, 2):
        for b in range(2, JCAP + 1, 2):
            J[a // 2, b // 2] = share(cache, (0, a, a + b) if a != b else (0, a, 2 * a), L, W)
    marg = np.zeros(K + 1)
    for h, v in g.items():
        marg[min(h // 2, K)] += v / tot
    J = J / tot
    for c in range(1, K):
        J[c, K] = marg[c] - J[c, 1:K].sum()
        J[K, c] = marg[c] - J[1:K, c].sum()
    J[K, K] = 1 - J[1:K + 1, 1:K + 1].sum() + J[K, K]
    return J[1:, 1:], marg[1:]


def H_bits(p):
    p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))


def values(win, L, W, cache):
    t = {}
    g = {h: share(cache, (0, h), L, W) for h in range(2, HMAX + 1, 2)}
    tot = sum(g.values())
    t['norm_gap_law'] = tot
    for h in PR.GAPS:
        t['S.gap%d' % h] = g[h]
    t['N.eqgap'] = sum(share(cache, (0, h, 2 * h), L, W) for h in range(2, EQMAX + 1, 2))
    for q in PR.SAME_Q:
        t['A.same%d' % q] = sum(v for h, v in g.items() if h % q == 0) / tot
    for H in PR.VAR_H:
        t['S.var%d' % H] = T2.var_ratio_aligned(H, L, W)
    J, marg = joint(cache, g, tot, L, W)
    Jc = np.clip(J, 0, None)
    px, py = Jc.sum(1), Jc.sum(0)
    nz = Jc > 0
    t['I.migap_lag1'] = float(np.sum(Jc[nz] * np.log2(Jc[nz] / np.outer(px, py)[nz])))
    t['T4.H_class'] = H_bits(marg)
    t['T4.H_cond'] = H_bits(Jc.ravel()) - H_bits(px)
    for th in tail_thresholds(win):
        t['T1.tail%d' % th] = sum(v for h, v in g.items() if h >= th) / tot
    t['T3.cpap4'] = sum(share(cache, (0, h, 2 * h, 3 * h), L, W) for h in range(6, EQMAX + 1, 6))
    t['_joint'] = J.tolist()
    return t


def main():
    if sys.argv[1] == 'build':
        build(); return
    cache = pickle.load(open(CACHE, 'rb'))
    missing = [B for B in patterns() if B not in cache]
    assert not missing, '%d patterns missing' % len(missing)
    res = {}
    for win, (L, W) in WIN.items():
        res[win] = {'-1': values(win, L, W, cache)}
        for b in range(10):
            res[win][str(b)] = values(win, L + b * (W // 10), W // 10, cache)
        print(win, {k: round(v, 6) for k, v in res[win]['-1'].items() if not k.startswith('_')}, flush=True)
    json.dump(res, open(os.path.join(OUT, 'theory_r3.json'), 'w'))


if __name__ == '__main__':
    main()
