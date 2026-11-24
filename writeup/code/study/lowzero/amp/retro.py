"""RETROSPECTIVE (all data already known). The amplitude gamma of the bad-prime off-diagonal part
 (1) in test B (frozen predictions, r134_B.json): reproduce 0.932, then split by conductor, sign of D, prime, length;
 (2) in the 181 real L-functions of rounds 117/120 (100 <= |D| <= 400, the sample of round 123 and of test C's
     discriminants) at the same heights 100-600, and in its two halves 100-300, 300-600.
usage: retro.py  -> retro_out.txt (stdout), retro.npz"""
import os, sys, json, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'bk'))
from common import *
import gridcache

t0 = time.time()
out = {}
# ---------------- (1) test B ----------------
Ds, S, F, Dg = load_B(); X = design(Ds); n = len(Ds)
g, se = gamma_boot(S, F, Dg, X, nb=1000)
print('[1] test B: n %d, gamma %.4f +- %.4f (bootstrap, 1000), pooled per-L-function slope %.4f' % (n, g, se, pooled(S, F, Dg)))
absD = np.array([abs(D) for D in Ds]); sg = np.array([D > 0 for D in Ds])
splits = {'|D|<=550': absD <= 550, '|D|>550': absD > 550, 'D>0': sg, 'D<0': ~sg}
for nm, m in splits.items():
    ix = np.nonzero(m)[0]
    gg, ss = gamma_boot(S[ix], F[ix], Dg[ix], X[ix], nb=500)
    print('    %-9s n %3d  gamma %.3f +- %.3f   pooled %.3f' % (nm, len(ix), gg, ss, pooled(S[ix], F[ix], Dg[ix])))
for i, p in enumerate(PR6):
    gg, ss = gamma_boot(S, F, Dg, X, nb=500, rows=slice(i, i + 1))
    print('    p = %2d   gamma %.3f +- %.3f' % (p, gg, ss))
for nm, cols in (('l<=2', [j for j, l in enumerate(LS) if l <= 2]), ('2<l<=5', [j for j, l in enumerate(LS) if 2 < l <= 5]), ('l>5', [j for j, l in enumerate(LS) if l > 5])):
    gg, ss = gamma_boot(S, F, Dg, X, nb=500, cols=cols)
    print('    %-8s gamma %.3f +- %.3f' % (nm, gg, ss))
# where the statistic draws its weight: share of sum off^2 per prime and per l-range
A = np.linalg.pinv(X.T @ X) @ X.T; off = (A @ F)[1:] - (A @ Dg)[1:]; res = (A @ S)[1:] - (A @ Dg)[1:]
w = off ** 2 / (off ** 2).sum()
print('    weight of each prime in gamma (share of sum off^2): ' + ' '.join('%d: %.2f' % (p, w[i].sum()) for i, p in enumerate(PR6)))
print('    weight of l > 5 for p = 2: %.2f of the total' % w[0, [j for j, l in enumerate(LS) if l > 5]].sum())
print('    p=2 res/off per l: ' + ' '.join('%.2f' % (res[0, j] / off[0, j]) if abs(off[0, j]) > 1e-3 else '  . ' for j in range(16)))
out['B'] = dict(gamma=g, se=se)

# ---------------- (2) rounds 117/120 sample at 100-600 ----------------
G = gridcache.get()
Z = load_r120_zeros(); D2 = sorted([int(k[1:]) for k in Z], key=lambda D: (abs(D), D)); n2 = len(D2)
X2 = design(D2)
res2 = {}
for (a, b) in ((100, 600), (100, 300), (300, 600)):
    S2 = np.array([sigma2_set(Z['D%d' % D], abs(D), a, b) for D in D2])
    F2, Dg2 = predict(G, D2, a, b)
    gg, ss = gamma_boot(S2, F2, Dg2, X2, nb=1000)
    Lm = np.mean([float(np.dot(*L_nodes(abs(D), a, b, 8))) for D in D2])
    print('[2] rounds 117/120 sample (n %d, 100 <= |D| <= 400), heights %d-%d, mean L %.3f: gamma %.4f +- %.4f, pooled %.4f, level max |mean(S - full)| %.4f' % (
        n2, a, b, Lm, gg, ss, pooled(S2, F2, Dg2), np.abs((S2 - F2).mean(0)).max()))
    res2[(a, b)] = (S2, F2, Dg2)
    np.savez(os.path.join(HERE, 'r120_%d_%d.npz' % (a, b)), S=S2, F=F2, Dg=Dg2, D=np.array(D2))
S2, F2, Dg2 = res2[(100, 600)]
absD2 = np.array([abs(D) for D in D2]); sg2 = np.array([D > 0 for D in D2])
for nm, m in {'|D|<=250': absD2 <= 250, '|D|>250': absD2 > 250, 'D>0': sg2, 'D<0': ~sg2}.items():
    ix = np.nonzero(m)[0]
    gg, ss = gamma_boot(S2[ix], F2[ix], Dg2[ix], X2[ix], nb=500)
    print('    %-9s n %3d  gamma %.3f +- %.3f   pooled %.3f' % (nm, len(ix), gg, ss, pooled(S2[ix], F2[ix], Dg2[ix])))
for i, p in enumerate(PR6[:3]):
    gg, ss = gamma_boot(S2, F2, Dg2, X2, nb=500, rows=slice(i, i + 1))
    print('    p = %2d   gamma %.3f +- %.3f' % (p, gg, ss))
for nm, cols in (('l<=4', [j for j, l in enumerate(LS) if l <= 4]), ('l>4', [j for j, l in enumerate(LS) if l > 4])):
    gg, ss = gamma_boot(S2, F2, Dg2, X2, nb=500, cols=cols)
    print('    %-8s gamma %.3f +- %.3f' % (nm, gg, ss))
# cross-check against round 123's published per-prime fits (same zeros, same statistic)
r123 = json.load(open(os.path.join(LOW, 'r123.json')))
A2 = np.linalg.pinv(X2.T @ X2) @ X2.T; fit2 = (A2 @ S2)[1:]
pub = np.array([[r123['fit']['%d_%s' % (p, float(l))] for l in LS] for p in PR6])
print('    check: our per-prime fits against r123.json: max |diff| %.2e' % np.abs(fit2 - pub).max())
print('done %.0fs' % (time.time() - t0))
