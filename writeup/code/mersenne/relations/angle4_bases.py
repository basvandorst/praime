"""Angle 4: digit features in bases 2..36 (terms with >= 2 digits in that base).

pal_b    number of palindromes                              (hi)
dist_b   mean of (#distinct digits / min(#digits, b))        (lo = structure)
few_b    number of terms with >= 4 digits and <= 2 distinct digits   (hi)
run_b    mean longest run of one repeated digit              (hi)
dsum_b   mean digit sum / (#digits * (b-1)/2)                (two-sided)
"""
import json, sys, time
import numpy as np

BASES = range(2, 37)


def digits_matrix(a, b):
    a = np.array(a, dtype=np.int64)
    rows = []
    x = a.copy()
    while np.any(x > 0):
        rows.append(np.where(x > 0, x % b, -1))
        x //= b
    D = np.array(rows).T  # LSB first, -1 padding
    nd = (D >= 0).sum(1)
    return D, nd


def base_feats(a, b):
    D, nd = digits_matrix(a, b)
    keep = nd >= 2
    D, nd = D[keep], nd[keep]
    pal = 0
    dist = []
    few = 0
    run = []
    dsum = []
    for row, m in zip(D, nd):
        d = row[:m]
        if np.array_equal(d, d[::-1]):
            pal += 1
        u = len(set(d.tolist()))
        dist.append(u / min(m, b))
        if m >= 4 and u <= 2:
            few += 1
        best = cur = 1
        for i in range(1, m):
            cur = cur + 1 if d[i] == d[i - 1] else 1
            best = max(best, cur)
        run.append(best)
        dsum.append(d.sum() / (m * (b - 1) / 2))
    return {'pal': pal, 'dist': float(np.mean(dist)), 'few': few, 'run': float(np.mean(run)),
            'dsum': float(np.mean(dsum))}


FEATS = {'pal': 'hi', 'dist': 'lo', 'few': 'hi', 'run': 'hi', 'dsum': 'two'}


def all_stats(a):
    out = {}
    for b in BASES:
        f = base_feats(a, b)
        for k in FEATS:
            out['%s_%d' % (k, b)] = f[k]
    return out


SIDES = {'%s_%d' % (k, b): s for b in BASES for k, s in FEATS.items()}

if __name__ == '__main__':
    D = json.load(open('nulls.json'))
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    t0 = time.time()
    A = D['real']
    real = all_stats(A)
    out = {s: {'side': SIDES[s], 'real': real[s], 'a': [], 'b': []} for s in SIDES}
    for fam in ('a', 'b'):
        for seq in D[fam][:nmax]:
            st = all_stats(seq)
            for s in SIDES:
                out[s][fam].append(st[s])
    json.dump(out, open('angle4.json', 'w'))
    print('done', time.time() - t0)
