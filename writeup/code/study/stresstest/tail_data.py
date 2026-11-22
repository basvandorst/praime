import sys, math
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import numpy as np, common as C
from scipy.stats import chi2 as CHI2
th = np.load('/home/bas/scratch/stress-B/tail.npz')
BINS = [(2, 36), (38, 52), (54, 70)]
def bi(g):
    return np.select([g <= 36, g <= 52, g <= 70], [0, 1, 2], 3)
def cells(A):
    C6 = np.zeros((3, 3, 2))
    for a in range(2, 71, 2):
        for b in range(2, 71, 2):
            C6[int(bi(a)), int(bi(b)), int(a % 6 == 0 and b % 6 == 0)] += A[a // 2 - 1, b // 2 - 1]
    return C6
for tag, L0 in (('1e9', 10 ** 9), ('1e12', 10 ** 12)):
    P = np.nonzero(C.primes_mask(L0, 10 ** 8 + 2000))[0].astype(np.int64) + L0
    n = int(np.sum(P < L0 + 10 ** 8)); g1 = P[1:n + 1] - P[:n]; g2 = P[2:n + 2] - P[1:n + 1]
    b1, b2 = bi(g1), bi(g2); s6 = ((g1 % 6 == 0) & (g2 % 6 == 0)).astype(int)
    O = np.zeros((3, 3, 2)); m = (b1 < 3) & (b2 < 3)
    np.add.at(O, (b1[m], b2[m], s6[m]), 1)
    print('==', tag, 'N', n)
    for nm in ('T', 'Tc'):
        E = n * cells(th[tag + '_' + nm]); keep = np.ones((3, 3, 2), bool); keep[0, 0, :] = False   # (A,A) was tested already
        z = (O - E) / np.sqrt(E); c = float(np.sum(z[keep] ** 2)); df = int(keep.sum())
        print(nm, 'chi2 %.2f df %d p %.4f' % (c, df, CHI2.sf(c, df)))
        print(' z not-both-6 :\n', np.round(z[:, :, 0], 1), '\n z both-6 :\n', np.round(z[:, :, 1], 1))
    print(' observed counts (not-both-6 / both-6):\n', O[:, :, 0].astype(int), '\n', O[:, :, 1].astype(int))
    for nm in ('T', 'Tc'):
        E = n * cells(th[tag + '_' + nm]); keep = np.ones((3, 3, 2), bool); keep[0, 0, :] = False
        print(tag, nm, 'aggregate tail mass: obs %d exp %.0f rel %+.4f z %+.2f' % (O[keep].sum(), E[keep].sum(), O[keep].sum() / E[keep].sum() - 1, (O[keep].sum() - E[keep].sum()) / math.sqrt(E[keep].sum())))
