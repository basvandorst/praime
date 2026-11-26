# True characters: smallest eigenvalue of the truncated form (float), and zeta with signs allowed.
import numpy as np, datum as boot6
for (kappa, q, name, pts) in [(1, 4, 'chi_-4', [(3,-1),(5,1),(7,-1)]), (1, 3, 'chi_-3', [(2,-1),(5,-1),(7,1)]), (0, 5, 'chi_5', [(2,-1),(3,-1),(7,-1)])]:
    for X in (4.5, 6, 9):
        F = boot6.Form(X, 20, False, kappa, q)
        p = [x for x, s in pts if x <= X]; s = [s for x, s in pts if x <= X]
        print(name, 'X', X, 'true character lam_min %.3e' % F.mineig(p, s), ' empty %.3e' % F.mineig([], []))
# zeta datum, signs allowed: one and two points at X = 3 and 4
for X in (3.0, 4.0):
    F = boot6.Form(X, 12, True, 0, 1.0)
    g = np.arange(1.04, X + 1e-9, 0.02); res = []
    for i, x in enumerate(g):
        for s1 in (1, -1):
            res.append((F.mineig([x], [s1]), x, s1))
            for y in g[i:]:
                for s2 in (1, -1):
                    res.append((F.mineig([x, y], [s1, s2]), x, s1, y, s2))
    res.sort(reverse=True)
    print('zeta datum with signs, X', X, 'top:', [tuple(round(float(v), 4) for v in r) for r in res[:6]])
