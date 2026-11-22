"""Independent check of r39.rho: the ratios of singular series from theory2.S1 (general code, all primes), with the
factors p | 210 removed by their closed values."""
import numpy as np, theory2 as T2, r39
loc = lambda nu_big, nu_parts, p: (1 - nu_big / p) / np.prod([(1 - v / p) for v in nu_parts])
def small(Sbig, parts, p_sets):
    f = 1.0
    for p in (2, 3, 5, 7):
        f *= loc(len({x % p for x in p_sets[0]}), [len({x % p for x in s}) for s in p_sets[1:]], p) * (1 - 1 / p) ** (len(p_sets[0]) - sum(len(s) for s in p_sets[1:]))
    return f
worst = 0
for g in (2, 6, 22):
    n = 400
    R = {kind: r39.rho(n, g, kind) for kind in ('XX', 'XYp', 'XYm', 'YY')}
    for k in (1, 11, 13, 37, 121, 143, 391):
        d = 210 * k
        cases = {'XX': ([0, g, d, d + g], [[0, g], [d, d + g]]), 'XYp': ([0, g, d], [[0, g], [d]]),
                 'XYm': ([-d, 0, g], [[0, g], [-d]]), 'YY': ([0, d], [[0], [d]])}
        for kind, (big, parts) in cases.items():
            if kind == 'YY' and g != 2:
                continue
            Sb = T2.S1(sorted(big)); Sp = np.prod([T2.S1(sorted(s)) if len(s) > 1 else 1.0 for s in parts])
            # the S of a single point is 1; remove p | 210 factors (both points of a pair in the same class mod 210)
            sm = small(Sb, parts, [big] + parts)
            val = Sb / Sp / sm
            worst = max(worst, abs(val / R[kind][k] - 1))
print('largest relative difference', worst)
