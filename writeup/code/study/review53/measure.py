"""Independent measurement of M1, M2 (all pairs), M3 (registered triples) from reserved_values.npz.
Trapezoid and Simpson on the full grid, Simpson on the even sub-grid (h = 0.1) as a quadrature check,
partial integrals over thirds [100,400], [400,700], [700,1000], top-1% contribution share.
Writes meas.npz in this directory."""
import json, itertools, numpy as np, os
HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/study/r149/reserved_values.npz')
S = json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/moment3/sets.json'))
R = S['reserved']; tri = [tuple(t) for t in S['triples']]; pairs = [tuple(p) for p in S['pairs']]
assert len(R) == 56 and len(tri) == 2000 and len(pairs) == 1540
V = {D: Z['absL_%d' % D].astype(float) ** 2 for D in R}
h = 0.05
n = 18001
t = 100 + h * np.arange(n)


def trap(y):
    return h * (y.sum() - 0.5 * (y[0] + y[-1]))


def simp(y, hh=h):
    return hh / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())


def parts(y):
    # thirds: k = 0..6000, 6000..12000, 12000..18000
    return [trap(y[a:b + 1]) for a, b in ((0, 6000), (6000, 12000), (12000, 18000))]


out = {}
for name, tuples in (('single', [(D,) for D in R]), ('pair', pairs), ('triple', tri)):
    M_tr, M_si, M_si2, M_parts, top1, mx = [], [], [], [], [], []
    for tp in tuples:
        y = np.ones(n)
        for D in tp:
            y = y * V[D]
        M_tr.append(trap(y)); M_si.append(simp(y)); M_si2.append(simp(y[::2], 2 * h))
        M_parts.append(parts(y))
        srt = np.sort(y)[::-1]
        top1.append(srt[:180].sum() / y.sum())
        mx.append(y.max() / y.mean())
    out[name + '_trap'] = np.array(M_tr); out[name + '_simp'] = np.array(M_si); out[name + '_simp2'] = np.array(M_si2)
    out[name + '_parts'] = np.array(M_parts); out[name + '_top1'] = np.array(top1); out[name + '_maxmean'] = np.array(mx)
# fourth moments of singles
out['single4_trap'] = np.array([trap(V[D] ** 2) for D in R])
out['single4_parts'] = np.array([parts(V[D] ** 2) for D in R])
np.savez(os.path.join(HERE, 'meas.npz'), **out)
for name in ('single', 'pair', 'triple'):
    a, b, c = out[name + '_trap'], out[name + '_simp'], out[name + '_simp2']
    print('%-6s trap vs simp: max rel %.2e  median %.2e ; simp(h) vs simp(2h): max rel %.2e' % (
        name, np.max(np.abs(a / b - 1)), np.median(np.abs(a / b - 1)), np.max(np.abs(b / c - 1))))
    print('       top-1%% of grid points carry: median %.3f  max %.3f ; max/mean median %.0f max %.0f' % (
        np.median(out[name + '_top1']), out[name + '_top1'].max(), np.median(out[name + '_maxmean']), out[name + '_maxmean'].max()))
