import numpy as np, sys
from common import *
ls = np.arange(0.5, 8.01, 0.5)
d = real_zeros()
print('n L-functions', len(d))
# check: all fundamental discriminants 100..400
import sympy
def isfund(D):
    if D % 4 == 1: return sympy.ntheory.factor_.core(abs(D)) == abs(D)
    if D % 4 == 0:
        m = D // 4
        return m % 4 in (2, 3) and sympy.ntheory.factor_.core(abs(m)) == abs(m)
    return False
fund = [D for D in range(-400, 401) if 100 <= abs(D) <= 400 and isfund(D)]
print('fundamental discs', len(fund), 'missing', sorted(set(fund) - set(d)), 'extra', sorted(set(d) - set(fund)))
rng = np.random.default_rng(1)
Ds = sorted(d)
res = {k: [] for k in ['dense', 'dense_lc', 'rand', 'trim', 'const']}
Ls = []; nz = []
for D in Ds:
    g = d[D]; q = abs(D)
    s, L, n = unfold_stats(g, q, 100, 600, ls); res['dense'].append(s); Ls.append(L); nz.append(n)
    res['dense_lc'].append(unfold_stats(g, q, 100, 600, ls, center='l')[0])
    res['rand'].append(unfold_stats(g, q, 100, 600, ls, mode='rand', rng=rng, nrand=2000)[0])
    res['trim'].append(unfold_stats(g, q, 100, 600, ls, frac=0.1)[0])
    # constant term -chi(-1)/8 + 1: sign: chi(-1)=sign(D)
    res['const'].append(unfold_stats(g, q, 100, 600, ls, const=7/8 if D > 0 else 9/8)[0])
np.savez('/home/bas/scratch/review-logic3/real.npz', Ds=np.array(Ds), ls=ls, L=np.array(Ls), nz=np.array(nz),
         **{k: np.array(v) for k, v in res.items()})
print('done')
