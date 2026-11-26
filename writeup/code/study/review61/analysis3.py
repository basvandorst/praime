# review 162, task 3: fine structure. Fine bins of u (1/192): data against the scaled own prediction; where are the large residuals;
# is there a spike at the predicted places u = k/3, k/6 (cube lines) or at the quadratic-resolvent places y = n^2/4; smooth-only fits.
import json, math, numpy as np
R = '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/review_162/'
PJ = json.load(open(R + 'pred_own_T1000.json'))
ue = np.linspace(0, 2, 385); uc = (ue[:-1] + ue[1:]) / 2
for r2 in (1, 2):
    z = np.load(R + 'mine_r%d.npz' % r2); n = len(z['D']); num = z['n192'].astype(float); den = z['d192'].astype(float)
    N = num.sum(0); Dn = den.sum(0); m = N / Dn; mj = (N[None] - num) / (Dn[None] - den); sig = np.sqrt((n - 1) / n * ((mj - mj.mean(0))**2).sum(0))
    P = np.array(PJ['sig'][str(r2)]['pred_u192']); hi = slice(56, 376); w = 1 / sig[hi]**2
    s = float(np.sum(P[hi] * m[hi] * w) / np.sum(P[hi]**2 * w)); zz = (m[hi] - s * P[hi]) / sig[hi]
    print('r2=%d: 320 fine bins, diag amplitude %.3f, chi2 %.0f; |z|>4 in %d bins; prediction range on fine bins %.1f .. %.1f' % (r2, s, np.sum(zz**2), np.sum(np.abs(zz) > 4), P[hi].min(), P[hi].max()))
    top = np.argsort(-np.abs(zz))[:8]
    print('   largest residuals (u, data, s*pred, z):', '; '.join('%.3f %.2f %.2f %+.1f' % (uc[56 + k], m[56 + k], s * P[56 + k], zz[k]) for k in sorted(top)))
    for lab, u0 in (('u=1/3', 1 / 3), ('u=1/2', 0.5), ('u=2/3', 2 / 3), ('u=1', 1.0), ('u=4/3', 4 / 3), ('u=5/3', 5 / 3), ('y=1/4', 0.25**(1 / 3)), ('y=9/4', 2.25**(1 / 3)), ('y=4', 4**(1 / 3))):
        k0 = int(round(u0 * 192)); ks = range(k0 - 3, k0 + 3)
        print('   %-6s bins %s | data %s | s*pred %s | sigma %.2f' % (lab, ' '.join('%.3f' % uc[k] for k in ks), ' '.join('%6.2f' % m[k] for k in ks), ' '.join('%6.2f' % (s * P[k]) for k in ks), sig[k0]))
    # smooth-only vs smooth + prediction shape (amplitude free)
    Y = uc[hi]**3; B = np.stack([np.ones_like(Y), Y**-0.5, Y**(-1 / 3), Y**-0.25], 1); sw = np.sqrt(w)
    def chi(X): c, *_ = np.linalg.lstsq(X * sw[:, None], m[hi] * sw, rcond=None); r = (m[hi] - X @ c) * sw; return float(r @ r), c
    c0, _ = chi(B); c1, cc = chi(np.column_stack([B, P[hi]])); c2, _ = chi(np.column_stack([B, np.array(PJ['sig'][str(r2)]['deg2_u192'])[hi]]))
    print('   smooth only (1, y^-1/2, y^-1/3, y^-1/4): chi2 %.0f; + deg-3 prediction (free amplitude %.3f): %.0f; + deg-2 comparator instead: %.0f' % (c0, cc[-1], c1, c2))
