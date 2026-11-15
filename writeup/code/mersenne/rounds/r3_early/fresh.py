"""Independent recheck of the 8 train-chosen single-feature rules (rules.json) on
fresh survivors 12000 < k < 20000 (exponent 19937), with separate code:
reduction by GMP division (%), residues by separate % p, i0 found by comparing
the unreduced value with M, features computed with numpy from per-iteration
scalars (one float and ten small residues per iteration; no big integers kept).

usage: fresh.py survivors.txt kmin kmax out.json [k1,k2,...]   (optional explicit k list,
used to cross-check against llrun.py on the test range)
"""
import sys, json, math, time
import numpy as np
import gmpy2
from gmpy2 import mpz
from multiprocessing import Pool

PR = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
EXP = {107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937}
RULES = json.load(open('rules.json'))
JH = {'j32': lambda k: 32, 'j64': lambda k: 64, 'j128': lambda k: 128, 'j256': lambda k: 256,
      'k/16': lambda k: k // 16, 'k/8': lambda k: k // 8, 'k/4': lambda k: k // 4, 'k/2': lambda k: k // 2}


def trajectory(k, kind, jmax):
    M = mpz(2) ** k - 1
    v = mpz(4) if kind == 'LL' else mpz(3)
    i0 = None
    hz = np.zeros(jmax + 1)
    res = np.zeros((jmax + 1, len(PR)), dtype=np.int64)
    for i in range(1, jmax + 1):
        w = v * v - 2 if kind == 'LL' else v * v
        if i0 is None and w >= M:
            i0 = i
        v = w % M
        hz[i] = (gmpy2.popcount(v) - k / 2) / (math.sqrt(k) / 2)
        res[i] = [int(v % p) for p in PR]
    return i0, hz, res


def feature(name, k, i0, hz, res, j):
    h = hz[i0:j + 1]; r = res[i0:j + 1]; n = len(h)
    if name == 'a_zero5_z':
        c = (r[:, PR.index(5)] == 0).sum()
        return (c - n / 5) / math.sqrt(n * 0.2 * 0.8)
    if name == 'a_hw_lag1z':
        return float((h[1:] * h[:-1]).sum() / math.sqrt(n - 1))
    if name == 'a_hw_meanz':
        return float(h.sum() / math.sqrt(n))
    if name == 'a_hw_sd':
        return float(h.std())
    if name == 'v_mod17':
        return res[j, PR.index(17)] / 17
    if name == 'a_chi2_z':
        chi, df = 0.0, 0
        for c, p in enumerate(PR):
            cnt = np.bincount(r[:, c], minlength=p)
            chi += ((cnt - n / p) ** 2).sum() / (n / p)
            df += p - 1
        return (chi - df) / math.sqrt(2 * df)
    raise KeyError(name)


def work(k):
    out = {'k': k}
    for kind in ('LL', 'PRP'):
        L = k - 2 if kind == 'LL' else k      # last index of the run
        jmax = max(k // 2, min(256, L))
        i0, hz, res = trajectory(k, kind, jmax)
        out[kind + '_i0'] = i0
        for hname, rule in RULES.items():
            seq, fname = rule['feature'].split(':')
            if seq != kind:
                continue
            j = JH[hname](k)
            if j > L:                         # same cap as analysis.py (k < 258 only)
                j = k // 2
            out[hname] = rule['sign'] * feature(fname, k, i0, hz, res, j)
    return out


if __name__ == '__main__':
    fn, kmin, kmax, outfn = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    if len(sys.argv) > 5:
        ks = [int(x) for x in sys.argv[5].split(',')]
    else:
        ks = [int(x) for x in open(fn) if kmin < int(x) < kmax]
    ks.sort(reverse=True)
    t0 = time.time()
    with Pool(2) as pool:
        rows = sorted(pool.imap_unordered(work, ks, chunksize=1), key=lambda r: r['k'])
    json.dump(rows, open(outfn, 'w'))
    print('runs', len(rows), 'wall %.0f s' % (time.time() - t0))
    if len(sys.argv) > 5:
        sys.exit()
    lab = np.array([r['k'] in EXP for r in rows])
    print('fresh survivors %d, exponents %s' % (len(rows), [r['k'] for r in rows if r['k'] in EXP]))
    print('horizon | rule (oriented, train threshold) | exponent kept | rejection of non-exponents | share of non-exponents scoring below the exponent')
    for hname, rule in RULES.items():
        s = np.array([r[hname] for r in rows])
        thr = rule['threshold']
        kept = bool((s[lab] >= thr).all())
        rej = (s[~lab] < thr).mean()
        below = (s[~lab] < s[lab].min()).mean()
        print('%-5s | %+d*%s >= %.4f | %s | %.3f | %.3f' % (hname, rule['sign'], rule['feature'], thr, kept, rej, below))
