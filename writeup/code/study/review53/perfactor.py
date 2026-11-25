"""Per-factor AFE diagonal ('diagonal of the product of the k per-factor approximate functional equations, sharp
cutoffs sqrt(X_j)') for registered triples and pairs, against the recipe P(t) (frozen hm, verified by myrecipe.py),
at t = 100, 316.23, 1000; integral over [100, 1000] of the ratio via quadratic interpolation in ln t weighted by P.
usage: perfactor.py NTRI NPAIR SEED"""
import sys, os, math, json, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from diagbrute import chi_vec, per_factor
from analysis import Pdens, MONOS3, MONOS2, gl
from stats import TT, PP, tri, pairs


def per_factor2(Ds, t):
    X = [abs(D) * t / (2 * math.pi) for D in Ds]
    x = [int(math.floor(math.sqrt(XX))) for XX in X]
    ch = [chi_vec(D, xx) for D, xx in zip(Ds, x)]
    m1 = np.arange(1, x[0] + 1); m2 = np.arange(1, x[1] + 1)
    Mm = (m1[:, None] * m2[None, :]).ravel(); W = (ch[0][1:][:, None] * ch[1][1:][None, :]).ravel()
    c = np.bincount(Mm, weights=W)
    n = np.arange(len(c)); n[0] = 1
    return 4 * np.sum(c[1:] ** 2 / n[1:])


ntri, npair, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(seed)
ts = np.array([100.0, 316.2277660168, 1000.0])
for kind, items, n in (('pair', pairs, npair), ('triple', tri, ntri)):
    idx = rng.choice(len(items), n, replace=False)
    rat_all = []
    for i in idx:
        Ds = items[i]
        if kind == 'triple':
            Pt = Pdens(Ds, TT[Ds]['hm'], ts, MONOS3); pfv = np.array([per_factor(Ds, t) for t in ts])
            Pint = TT[Ds]['full']
        else:
            Pt = Pdens(Ds, PP[Ds]['hm'], ts, MONOS2); pfv = np.array([per_factor2(Ds, t) for t in ts])
            Pint = PP[Ds]['full']
        rat = pfv / Pt
        co = np.polyfit(np.log(ts), rat, 2)
        mon = MONOS3 if kind == 'triple' else MONOS2
        hm = TT[Ds]['hm'] if kind == 'triple' else PP[Ds]['hm']
        I = gl(lambda tt: np.polyval(co, np.log(tt)) * Pdens(Ds, hm, tt, mon), 100, 1000) / Pint
        rat_all.append(I)
        print('%-6s %-24s ratio at t: %s  integral ratio %.4f' % (kind, Ds, np.round(rat, 4).tolist(), I), flush=True)
    rat_all = np.array(rat_all)
    print('%s: per-factor diagonal / recipe over [100,1000]: mean %.4f sd %.4f min %.4f max %.4f (n %d)' % (
        kind, rat_all.mean(), rat_all.std(ddof=1), rat_all.min(), rat_all.max(), len(rat_all)), flush=True)
