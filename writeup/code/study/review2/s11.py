"""Section 11 check: own S_G, exact inclusion-exclusion for P(d = 2), P(d = 4), and own nearest-neighbour counts
in box G1 (a, b in [1e5, 1e5 + 2000)), with a mirror (b + ai) analysis of the block test."""
import math, json, sys, time
from itertools import combinations
import numpy as np
import gmpy2
from sieve import small_primes

PR = small_primes(2 * 10 ** 6)
KAPPA = 4 / math.pi


def sqrtm1(p):
    for g in range(2, p):
        r = pow(g, (p - 1) // 4, p)
        if r * r % p == p - 1:
            return r


def SG(pts):
    """pts: list of (x, y). prod over Gaussian primes pi of (1 - nu/N)/(1 - 1/N)^k."""
    k = len(pts); logS = 0.0
    nu = len({(x + y) % 2 for x, y in pts})
    if nu >= 2:
        return 0.0
    logS += math.log(1 - nu / 2) - k * math.log(0.5)
    for p in PR[1:]:
        p = int(p)
        if p % 4 == 1:
            if p > 200:     # all residues distinct beyond the span: generic factor
                logS += 2 * (math.log(1 - k / p) - k * math.log(1 - 1 / p)); continue
            r = sqrtm1(p)
            for s in (r, p - r):
                nu = len({(x + y * s) % p for x, y in pts})
                if nu >= p:
                    return 0.0
                logS += math.log(1 - nu / p) - k * math.log(1 - 1 / p)
        else:
            if p * p > 4 * 10 ** 12:
                break
            nu = len({(x % p, y % p) for x, y in pts})
            if nu >= p * p:
                return 0.0
            logS += math.log(1 - nu / p ** 2) - k * math.log(1 - 1 / p ** 2)
    return math.exp(logS)


def avoid_exact(pts, lam):
    """P(no Gaussian prime at pts | 0 prime) = sum_A (-1)^|A| S(0 u A)/S(0) lam^-|A|."""
    S0 = SG([(0, 0)])
    tot = 0.0
    for k in range(len(pts) + 1):
        for A in combinations(pts, k):
            tot += (-1) ** k * SG([(0, 0)] + list(A)) / S0 / lam ** k
    return tot


def nn_box(A, B, K, margin=12):
    lo, loB, n = A - margin, B - margin, K + 2 * margin
    M = np.zeros((n, n), dtype=bool)
    for i in range(n):
        a = lo + i; a2 = a * a
        M[i] = [gmpy2.is_prime(a2 + (loB + j) ** 2) for j in range(n)]
    offs = sorted([(x, y) for x in range(-margin, margin + 1) for y in range(-margin, margin + 1)
                   if 0 < x * x + y * y <= 82], key=lambda v: v[0] ** 2 + v[1] ** 2)
    core = M[margin:margin + K, margin:margin + K]
    nn = np.zeros((K, K), dtype=np.int64)
    for x, y in offs:
        sh = M[margin + x:margin + x + K, margin + y:margin + y + K]
        upd = core & sh & (nn == 0)
        nn[upd] = x * x + y * y
    assert not np.any(core & (nn == 0))
    return core, nn


if __name__ == '__main__':
    import sys
    BOX = sys.argv[1] if len(sys.argv) > 1 else 'G1'
    out = {}
    out['pair_const_0_1+i'] = SG([(0, 0), (1, 1)])
    A, B, K = {'G1': (10 ** 5, 10 ** 5, 2000), 'G2': (5 * 10 ** 6, 5 * 10 ** 6, 1000)}[BOX]
    lam = math.log((A + K / 2) ** 2 + (B + K / 2) ** 2) / KAPPA
    n2 = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    n4 = [(2, 0), (-2, 0), (0, 2), (0, -2)]
    ge4 = avoid_exact(n2, lam); ge8 = avoid_exact(n2 + n4, lam)
    out['exact'] = {'P(d=2)': 1 - ge4, 'P(d=4)': ge4 - ge8, 'P(d>=8)': ge8, 'lam_eff': lam,
                       'density': KAPPA / math.log((A + K / 2) ** 2 + (B + K / 2) ** 2)}
    print(out, flush=True)
    t0 = time.time()
    core, nn = nn_box(A, B, K)
    print('sieve', round(time.time() - t0), int(core.sum()), flush=True)
    a = np.arange(K)[:, None] + A; b = np.arange(K)[None, :] + B
    aa = np.broadcast_to(a, (K, K)); bb = np.broadcast_to(b, (K, K))
    # mirror symmetry of the data: nn(a, b) == nn(b, a)
    out['count'] = int(core.sum())
    out['mirror_identical'] = bool(np.array_equal(core, core.T) and np.array_equal(nn, nn.T))
    T = json.load(open('/home/bas/praime/writeup/data/study/theory_r11.json'))[BOX]
    stats = {}
    for d in (2, 4, 8, 10):
        law = [T[str(k)]['d=%d' % d] for k in range(10)]
        # (i) as in the study: blocks = strips of b, all points
        m1 = [float(np.mean(nn[:, k * (K // 10):(k + 1) * (K // 10)][core[:, k * (K // 10):(k + 1) * (K // 10)]] == d)) for k in range(10)]
        r = np.array(m1) - np.array(law); t1 = r.mean() / (r.std(ddof=1) / math.sqrt(10))
        # (ii) one point per mirror pair (a < b), 10 contiguous groups of equal size in b
        sel = core & (aa < bb)
        bv = bb[sel]; dv = nn[sel]; o = np.argsort(bv, kind='stable'); dv = dv[o]
        g = np.array_split(dv == d, 10); m2 = np.array([x.mean() for x in g])
        r2 = m2 - T['-1']['d=%d' % d]; t2 = r2.mean() / (r2.std(ddof=1) / math.sqrt(10))
        # (iii) the same 10 strips, but the a < b half only, combined with weights
        stats[d] = {'data_all': float(np.mean(nn[core] == d)), 'law': T['-1']['d=%d' % d],
                    't9_strips_all': float(t1), 'blocks_all': m1,
                    'data_dedup': float(np.mean(dv == d)), 't9_dedup_equal_groups': float(t2),
                    'sd_blocks_all': float(np.std(m1, ddof=1)), 'sd_blocks_dedup': float(np.std(m2, ddof=1)),
                    'binom_sd_block_all': float(math.sqrt(0.2 * 0.8 / (core.sum() / 10))),
                    'binom_sd_block_dedup': float(math.sqrt(0.2 * 0.8 / (sel.sum() / 10)))}
        print(d, {k: v for k, v in stats[d].items() if k != 'blocks_all'}, flush=True)
    out['stats'] = stats
    # correlation of strip k and strip l block means for d = 10 through shared mirror squares
    json.dump(out, open('/home/bas/scratch/study-review2/s11_%s.json' % BOX, 'w'), indent=1)
