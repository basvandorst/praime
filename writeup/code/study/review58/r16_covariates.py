"""Genus-2 residual Re(G - pred) at T = 150 against covariates: conductor, c, root number/rank, reduction types, model degree,
moments; group means with jackknife errors and simple slopes."""
import json, math, sys
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import jk_mean
T = {o['key']: o for o in json.load(open('r2_table.json'))}
A = {}
for line in open('r4_arith.jsonl'):
    o = json.loads(line); A[o['key']] = o
keys = list(T)
res = np.array([T[k]['res'][0] for k in keys]); resi = np.array([T[k]['res'][1] for k in keys])
def grp(lab, mask):
    a = res[mask]; b = res[~mask]
    ma, sa = jk_mean(a) if len(a) > 2 else (a.mean(), float('nan')); mb, sb = jk_mean(b) if len(b) > 2 else (b.mean(), float('nan'))
    print('%-40s yes n=%2d %+5.2f +- %4.2f | no n=%2d %+5.2f +- %4.2f' % (lab, len(a), ma, sa, len(b), mb, sb))
N = np.array([T[k]['N'] for k in keys]); eps = np.array([T[k]['eps'] for k in keys]); c = np.array([T[k]['c'] for k in keys])
bad = [A[k]['bad'] for k in keys]
grp('prime conductor', np.array([len(b) == 1 and b[0][1] == 1 for b in bad]))
grp('2 | N', np.array([any(p == 2 for p, e in b) for b in bad]))
grp('3 | N', np.array([any(p == 3 for p, e in b) for b in bad]))
grp('some exponent >= 2', np.array([any(e >= 2 for p, e in b) for b in bad]))
grp('N > 1500', N > 1500)
grp('eps = +1', eps > 0)
grp('c > 1', c > 1)
deg5 = np.array(['x^6' not in T[k]['spec'][1] for k in keys])
grp('deg f = 5 (rational Weierstrass point at infinity)', deg5)
hzero = np.array([T[k]['spec'][2] in ('0', 'x^0') for k in keys])
grp('h in {0, 1}', hzero)
m4 = np.array([A[k]['m4'] for k in keys]); grp('m4 > 2.93 (median-ish)', m4 > np.median(m4))
for lab, x in (('log N', np.log(N)), ('c', c), ('eps', eps), ('m4', m4), ('m2', np.array([A[k]['m2'] for k in keys])), ('mean lambda_p', np.array([A[k]['m1'] for k in keys]))):
    X = np.c_[np.ones_like(x), x]; b = np.linalg.lstsq(X, res, rcond=None)[0]
    loo = np.array([np.linalg.lstsq(np.delete(X, i, 0), np.delete(res, i), rcond=None)[0] for i in range(len(res))]); n = len(res)
    se = np.sqrt((n - 1) / n * np.sum((loo - loo.mean(0)) ** 2, 0))
    print('slope of residual on %-14s %+.3f +- %.3f (z %+.1f)' % (lab, b[1], se[1], b[1] / se[1]))
