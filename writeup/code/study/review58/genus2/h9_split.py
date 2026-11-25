"""Exploratory: split the law's constant into the all-integers part 1 (dual index m = 1, fact 2053) and the arithmetic
part c - 1: Re G = b1 (-eps T/2pi) + b2 (-eps (c - 1) T/2pi), least squares, jackknife over curves; and a power law
Re G ~ b (-eps T/2pi) |c|^beta sign(c) (grid in beta). Genus 2 (100) and, for comparison, Sym^3 + RS (95)."""
import sys, os, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from g2data import load
for cls in (('g2',), ('sym3', 'rs')):
    rows = load(cls); y = np.array([r['G'].real for r in rows]); U = np.array([-r['eps'] * r['T2'] / (2 * math.pi) for r in rows]); c = np.array([r['c'] for r in rows])
    X = np.vstack([U, U * (c - 1)]).T; n = len(y)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    loo = np.array([np.linalg.lstsq(X[np.arange(n) != i], y[np.arange(n) != i], rcond=None)[0] for i in range(n)])
    se = np.sqrt((n - 1) / n * ((loo - loo.mean(0)) ** 2).sum(0))
    print('%-9s n=%3d  b1 (part 1) %.3f +- %.3f   b2 (part c - 1) %.3f +- %.3f' % ('+'.join(cls), n, beta[0], se[0], beta[1], se[1]))
    best = None
    for bt in np.arange(0.6, 1.41, 0.05):
        Z = U * np.sign(c) * np.abs(c) ** bt; b = (Z * y).sum() / (Z * Z).sum(); rss = ((y - b * Z) ** 2).sum()
        if best is None or rss < best[2]: best = (bt, b, rss)
    print('          power law: best beta %.2f (b %.3f, rss %.0f); beta = 1 rss %.0f' % (best[0], best[1], best[2], ((y - ((U * c * y).sum() / ((U * c) ** 2).sum()) * U * c) ** 2).sum()))
