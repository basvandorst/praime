"""Profinite CG(P) at lam = 24.16 against round 139's fresh windows (retrospective: the fresh numbers were known; the model
has no parameter fitted to them).  Per cell (q, k): all-direction mean of W at 1e4, 1e5, 1e6 and the slope
(W(1e6) - W(1e4))/ln 100; error = window scatter/sqrt(n).  Fresh numbers and the y = 17 comparator from
writeup/data/study/r139/rep_fresh.json (read only)."""
import json, numpy as np, math, glob
F = json.load(open('../../../../writeup/data/study/r139/rep_fresh.json'))['cells']
FO = json.load(open('firstorder.json'))['cells']
HS = ('10000', '100000', '1000000')
res = {}
for fn in sorted(glob.glob('cg4_*_24.16.json')):
    P = int(fn.split('_')[1]); J = json.load(open(fn)); n = len(J)
    print('CG(%d), %d windows, lam 24.16' % (P, n))
    print('  cell   CG W(1e4) W(1e5) W(1e6)   CG slope          fresh slope        LR17      CG/LR  fresh/LR  first-order/LR  z(fresh-CG)')
    for cell in sorted(next(iter(J.values())).keys(), key=lambda c: tuple(map(int, c.split('_')))):
        W = {H: np.array([np.mean(v[cell]['blk'][H]) for v in J.values()]) for H in HS}
        sl = (W['1000000'] - W['10000']) / math.log(100); s, e = sl.mean(), sl.std(ddof=1) / math.sqrt(n)
        f = F[cell]['all']; fs, fe = f['slope']; lr = f['lr']
        z = (fs - s) / math.hypot(fe, e)
        res['%d_%s' % (P, cell)] = dict(n=n, W=[float(W[H].mean()) for H in HS], slope=[float(s), float(e)], fresh=[fs, fe], lr=lr, z=z)
        print('  %-5s  %.4f %.4f %.4f   %+.5f+-%.5f   %+.5f+-%.5f   %+.5f   %.2f   %.2f      %.2f           %+.1f' % (
            cell, *[W[H].mean() for H in HS], s, e, fs, fe, lr, s / lr, fs / lr, FO[cell]['s_inf'] / lr, z))
json.dump(res, open('cg4sum.json', 'w'), indent=1)
