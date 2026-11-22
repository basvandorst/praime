"""Size of successive Bahadur terms for large gaps: ln-rate change from order 2 -> 3, and the pair-overlap correction."""
import sys, math, numpy as np
sys.path.insert(0, '/home/bas/scratch/study-review4')
import u_check as UC
for g, lams in ((200, (19.0, 23.0)), (300, (21.4, 23.0)), (420, (25.5,))):
    SB, n, L3 = UC.rate_fn(g)
    # recompute pieces
    import itertools
    for lam in lams:
        # order-2 log-avoid and correction via a second call path
        T = {k: float(UC.tail_T(k, g)) for k in (2, 3, 4)}
        I = np.arange(1, g)
        s1 = UC.S_rows(np.stack([np.zeros_like(I), np.full_like(I, g), I], 1), g, T) / SB
        live = I[s1 > 0]; rho = s1[s1 > 0]; m = len(live)
        i, j = np.triu_indices(m, 1)
        s2v = UC.S_rows(np.stack([np.zeros_like(i), np.full_like(i, g), live[i], live[j]], 1), g, T) / SB
        p = rho / lam
        L1 = np.sum(np.log1p(-p))
        c = (s2v / lam ** 2 - p[i] * p[j]) / ((1 - p[i]) * (1 - p[j]))
        L2 = L1 + c.sum()
        K = np.zeros((m, m)); K[i, j] = c; K[j, i] = c
        r = K.sum(1)
        corr = -0.5 * (np.sum(c * c) + np.sum(r * r - (K * K).sum(1)))
        l3 = L3(lam)
        print('g %d lam %.1f g/lam %.1f live %d: sum p %.2f; L1 %.3f, order-2 term %+.4f, order-3 term %+.4f, pair-overlap corr %+.4f' % (
            g, lam, g / lam, m, p.sum(), L1, L2 - L1, l3 - L2, corr), flush=True)
