# review 162: (i) the last-bin bug: registered readings with the last HI bin as registered (u <= 47/24) and as computed (u <= 2);
# (ii) small y: data / prediction in the LO bins for discovery and holdout (is m at small y a function of y alone?)
import json, math, numpy as np
Q = '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/quartic/data/'
F = json.load(open('/home/bas/praime/writeup/code/study/lowzero/universal/quartic/predictions_quartic.json'))
def jk(num, den):
    n = len(num); N = num.sum(0); Dn = den.sum(0); m = N / Dn; mj = (N[None] - num) / (Dn[None] - den); dm = mj - mj.mean(0); return m, (n - 1) / n * dm.T @ dm, mj
for r2 in (1, 2):
    S = F['sig'][str(r2)]; P = np.array(S['pred']); P2 = np.array(S['comparators']['deg2_kernel'])
    for rng in ('hold', 'disc'):
        z = np.load(Q + '%s_r%d.npz' % (rng, r2)); num, den = z['num'], z['den']
        fi = [0] + [8 * j + 4 for j in range(24)]
        for lab, last in (('as registered (u<=47/24)', 188), ('as run (u<=2)', 192)):
            fj = fi[:-1] + [last]
            nc = np.stack([num[:, fj[k]:fj[k + 1]].sum(1) for k in range(24)], 1); dc = np.stack([den[:, fj[k]:fj[k + 1]].sum(1) for k in range(24)], 1)
            m, C, mj = jk(nc[:, 4:24], dc[:, 4:24]); Ci = np.linalg.inv(C); p = P[4:24]; w = 1 / np.diag(C)
            a = p @ Ci @ m / (p @ Ci @ p); c2 = (m - a * p) @ Ci @ (m - a * p); q = P2[4:24]; b = q @ Ci @ m / (q @ Ci @ q); c22 = (m - b * q) @ Ci @ (m - b * q)
            sd = (p * m * w).sum() / (p * p * w).sum()
            print('r2=%d %s %-26s GLS amp %.3f chi2 %.1f, dchi2 deg2 %.0f, diag amp %.3f, corr %.3f' % (r2, rng, lab, a, c2, c22 - c2, sd, np.corrcoef(m, p)[0, 1]))
        mL, CL, _ = jk(nc[:, 1:4], dc[:, 1:4]); print('      LO bins u in (1/24,7/24]: data %s, prediction %s, ratio %s, |D| ~ %d' % (np.round(mL, 2), np.round(P[1:4], 2), np.round(mL / P[1:4], 2), int(np.median(np.abs(z['D'])))))
