import numpy as np
from common import *
d = r122_zeros(); Ds = sorted(d)
dall = real_zeros()
print('n', len(Ds), 'subset of 181:', set(Ds) <= set(dall), 'every third offset0:', Ds == sorted(dall)[0::3], len(sorted(dall)[0::3]))
ls = np.array([0.5, 1.0, 1.5, 2.0])
res = {}
for mode in ['dense', 'trim', 'rand']:
    S = []; Ls = []; nz = []
    for D in Ds:
        g = d[D]
        kw = dict(frac=0.1) if mode == 'trim' else (dict(mode='rand', rng=np.random.default_rng(abs(D)), nrand=2000) if mode == 'rand' else {})
        s, L, n = unfold_stats(g, abs(D), 900, 1500, ls, **kw); S.append(s); Ls.append(L); nz.append(n)
    S = np.array(S); Ls = np.array(Ls)
    x = np.array([sum(dS(1, p, Ls[k]) for p in primes_of(D)) for k, D in enumerate(Ds)])
    b, se = ols(np.column_stack([np.ones(len(Ds)), x]), S[:, 1])
    print(mode, 'beta %.3f +- %.3f a %.4f' % (b[1], se[1], b[0]), 'mean L %.2f' % Ls.mean(), 'nz mean %.0f min %d max %d' % (np.mean(nz), min(nz), max(nz)))
    if mode == 'dense':
        y = S[:, 1]; rng = np.random.default_rng(5); bb = []
        for _ in range(2000):
            i = rng.integers(0, len(y), len(y)); bb.append(ols(np.column_stack([np.ones(len(y)), x[i]]), y[i])[0][1])
        print('bootstrap sd %.3f' % np.std(bb))
        two = np.array([abs(D) % 2 == 0 for D in Ds])
        print('mean pred drop for 2|D: %.4f' % np.mean([dS(1, 2, Ls[k]) for k in range(len(Ds)) if two[k]]))
        print('check zero range', min(g.min() for g in d.values()), max(g.max() for g in d.values()))
        # x with only primes <= 13 vs all
        x13 = np.array([sum(dS(1, p, Ls[k]) for p in primes_of(D) if p <= 13) for k, D in enumerate(Ds)])
        b2, se2 = ols(np.column_stack([np.ones(len(Ds)), x13]), y); print('only p<=13: beta %.3f +- %.3f' % (b2[1], se2[1]))
        sign = np.array([D > 0 for D in Ds], float)
        b3, se3 = ols(np.column_stack([np.ones(len(Ds)), x, sign, Ls]), y); print('with sign, L: beta %.3f +- %.3f' % (b3[1], se3[1]))
        for j, l in enumerate(ls):
            xl = np.array([sum(dS(l, p, Ls[k]) for p in primes_of(D)) for k, D in enumerate(Ds)])
            bl, sel = ols(np.column_stack([np.ones(len(Ds)), xl]), S[:, j]); print('  l=%.1f beta %.3f +- %.3f' % (l, bl[1], sel[1]))
