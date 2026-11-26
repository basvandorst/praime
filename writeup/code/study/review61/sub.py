# review 162: are the registered readings carried by the first bins? correlation and delta chi2 on sub-ranges of the HI bins (own data, corrected last bin)
import json, math, numpy as np
R = '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/review_162/'
PJ = json.load(open(R + 'pred_own_T1000.json'))
for r2 in (1, 2):
    z = np.load(R + 'mine_r%d.npz' % r2); n192 = z['n192'].astype(float); d192 = z['d192'].astype(float); n = len(n192)
    fi = [0] + [16 * j + 8 for j in range(24)]
    nc = np.stack([n192[:, fi[k]:fi[k + 1]].sum(1) for k in range(24)], 1); dc = np.stack([d192[:, fi[k]:fi[k + 1]].sum(1) for k in range(24)], 1)
    P3 = np.array(PJ['sig'][str(r2)]['pred_reg']); P2 = np.array(PJ['sig'][str(r2)]['deg2_reg'])
    for a0, a1 in ((4, 24), (4, 12), (12, 24), (8, 24), (16, 24)):
        num, den = nc[:, a0:a1], dc[:, a0:a1]; N = num.sum(0); Dn = den.sum(0); m = N / Dn; mj = (N[None] - num) / (Dn[None] - den)
        C = (n - 1) / n * (mj - mj.mean(0)).T @ (mj - mj.mean(0)); Ci = np.linalg.inv(C)
        def g(P): a = P @ Ci @ m / (P @ Ci @ P); r = m - a * P; return a, r @ Ci @ r
        a3, c3 = g(P3[a0:a1]); a2, c2 = g(P2[a0:a1]); cr = np.array([np.corrcoef(mj[i], P3[a0:a1])[0, 1] for i in range(n)])
        print('r2=%d coarse bins %2d..%2d (u %.3f..%.3f): corr deg3 %.3f +- %.3f, deg2 %.3f; GLS amp %.3f, chi2 %.1f on %d; dchi2(deg2-deg3) %.0f'
              % (r2, a0, a1 - 1, (2 * a0 - 1) / 24, (2 * a1 - 1) / 24, np.corrcoef(m, P3[a0:a1])[0, 1], math.sqrt((n - 1) / n * ((cr - cr.mean())**2).sum()), np.corrcoef(m, P2[a0:a1])[0, 1], a3, c3, a1 - a0, c2 - c3))
