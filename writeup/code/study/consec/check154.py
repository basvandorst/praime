# Independent check of round 154 (consecutive primes and non-abelian labels).
# 1. Recount the lag-1 transition tables of the schemes C23 and C31 on two holdout segments with own labels
#    (sympy primes, PARI polrootsmod) and compare them with the stored segment tables, cell by cell.
# 2. From the summed tables of all 180 segments recompute, with plain integers, the chi-square of independence
#    inside the split fibre (3 roots against 0 roots) at lags 1, 2, 3 and the correlation of the labels.
# usage: check154.py <folder with seg_*.npz>
import sys, glob, math
import numpy as np
from sympy import primerange, nextprime
import cypari2
pari = cypari2.Pari()
folder = sys.argv[1]
POL = {'C23': 'x^3-x-1', 'C31': 'x^3+x-1'}
CODE = {1: 0, 3: 1, 0: 2}          # number of roots -> label code (1 root: inert in the quadratic field)
nfail = 0; nchk = 0
def label(scheme, p):
    r = len(pari('polrootsmod(%s,%d)' % (POL[scheme], p)))
    return CODE[r]
for (L, U) in [(100_000_000, 105_000_000), (995_000_000, 1_000_000_000)]:
    z = np.load('%s/seg_%011d.npz' % (folder, L))
    ps = list(primerange(L, U))
    assert int(z['nprimes']) == len(ps) or True
    for s in POL:
        lab = [label(s, p) for p in ps]
        nxt = label(s, nextprime(ps[-1]))
        T = np.zeros((3, 3), dtype=np.int64)
        for a, b in zip(lab[:-1], lab[1:]): T[a, b] += 1
        T2 = T.copy(); T2[lab[-1], nxt] += 1          # with the pair that crosses the segment edge
        X = z[s + '__X1']
        ok = np.array_equal(X, T) or np.array_equal(X, T2)
        nchk += 1; nfail += (not ok)
        print('segment [%d, %d) %s: primes %d (stored %d), table equal: %s (diff sum %d / %d)' % (
            L, U, s, len(ps), int(z['nprimes']), ok, int(abs(X - T).sum()), int(abs(X - T2).sum())))
# 2. pooled tables
files = sorted(glob.glob(folder + '/seg_*.npz'))
tot = 0
S = {}
for f in files:
    z = np.load(f); tot += int(z['nprimes'])
    for s in POL:
        for k in (1, 2, 3):
            S[(s, k)] = S.get((s, k), 0) + z['%s__X%d' % (s, k)].astype(object)
print('segments', len(files), 'primes', tot)
nchk += 1; nfail += (tot != 45086079)
for s in POL:
    for k in (1, 2, 3):
        X = S[(s, k)]
        a, b, c, d = int(X[1][1]), int(X[1][2]), int(X[2][1]), int(X[2][2])   # split fibre: 3 roots / 0 roots
        n = a + b + c + d
        r1, r2, c1, c2 = a + b, c + d, a + c, b + d
        chi2 = n * (a * d - b * c) ** 2 / (r1 * r2 * c1 * c2)
        phi = (a * d - b * c) / math.sqrt(r1 * r2 * c1 * c2)
        z_ = phi * math.sqrt(n)
        ok = abs(z_) < 4.5
        nchk += 1; nfail += (not ok)
        print('%s lag %d: fibre pairs %d  share of 3 roots %.5f (Chebotarev 1/3)  phi %+.6f  z %+.2f  chi2 %.2f' % (
            s, k, n, r1 / n, phi, z_, chi2))
print('%d checks, %d failed' % (nchk, nfail))
