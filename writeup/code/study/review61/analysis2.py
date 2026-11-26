# review 162, task 2: generic oscillation fit with free alpha: m = b0 + b1 y^-1/2 + b2 y^-1/3 + y^-beta (a1 cos + a2 sin)(2 pi c y^alpha),
# alpha, c, beta on grids (c through the number of cycles K over the range), linear parameters by weighted least squares (diagonal
# jackknife errors). Fine bins of u (width 1/192, u in (7/24, 47/24]) and the 20 bins in y and in sqrt(y). Calibration: the same fit
# on the noise-free degree-3 prediction and on the degree-2 comparator. Error of alpha: 20-group jackknife over fields.
import sys, json, math, time, numpy as np
R = '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/review_162/'
PJ = json.load(open(R + 'pred_own_T1000.json')); t0 = time.time()
ylo, yhi = (7 / 24)**3, (47 / 24)**3
AL = np.round(np.arange(0.20, 0.7001, 0.01), 3); KS = np.arange(0.5, 25.0001, 0.05); BE = np.round(np.arange(-0.1, 0.5001, 0.1), 2)
def design_points(edges, nsub):
    # sub-points uniform in y inside each bin (prime density ~ uniform in p within a bin)
    e = np.asarray(edges); f = (np.arange(nsub) + 0.5) / nsub; return e[:-1, None] + (e[1:] - e[:-1])[:, None] * f[None, :]
def profile(m, w, Y):
    """Y: (nb, nsub) sub-points; returns chi2[alpha, beta] minimised over K, and argmin K."""
    sw = np.sqrt(w); Bk = np.stack([np.ones_like(Y), Y**-0.5, Y**(-1 / 3)], -1).mean(1) * sw[:, None]
    Q, _ = np.linalg.qr(Bk); mw = m * sw; mp = mw - Q @ (Q.T @ mw); base = float(mp @ mp)
    out = np.full((len(AL), len(BE)), np.inf); arg = np.zeros((len(AL), len(BE)), int)
    for ia, al in enumerate(AL):
        span = yhi**al - ylo**al; ph = 2 * math.pi * (KS / span)[:, None, None] * (Y**al)[None]     # (nK, nb, nsub)
        cs, sn = np.cos(ph), np.sin(ph)
        for ib, be in enumerate(BE):
            amp = Y**-be
            O1 = (cs * amp[None]).mean(2) * sw[None]; O2 = (sn * amp[None]).mean(2) * sw[None]      # (nK, nb)
            O1 = O1 - (O1 @ Q) @ Q.T; O2 = O2 - (O2 @ Q) @ Q.T
            a11 = (O1 * O1).sum(1); a22 = (O2 * O2).sum(1); a12 = (O1 * O2).sum(1); b1 = O1 @ mp; b2 = O2 @ mp
            det = a11 * a22 - a12**2; red = (a22 * b1**2 - 2 * a12 * b1 * b2 + a11 * b2**2) / det
            c2 = base - red; k = int(np.argmin(c2)); out[ia, ib] = c2[k]; arg[ia, ib] = k
    return out, arg
def summarize(ch, arg):
    prof = ch.min(1); ia = int(np.argmin(prof))
    if 0 < ia < len(AL) - 1:   # parabolic refinement
        y0, y1, y2 = prof[ia - 1:ia + 2]; d = (y0 - 2 * y1 + y2); ah = AL[ia] + 0.01 * 0.5 * (y0 - y2) / d if d > 0 else AL[ia]
    else: ah = AL[ia]
    ib = int(np.argmin(ch[ia])); K = KS[arg[ia, ib]]
    i3 = int(np.argmin(np.abs(AL - 1 / 3))); i2 = int(np.argmin(np.abs(AL - 0.5)))
    return dict(alpha=float(ah), chi_min=float(prof[ia]), beta=float(BE[ib]), cycles=float(K), chi_a033=float(prof[i3]), chi_a050=float(prof[i2]),
                prof=prof.tolist())
res = {}
for r2 in (1, 2, 0):
    z = np.load(R + 'mine_r%d.npz' % r2); S = PJ['sig'][str(r2)]; n = len(z['D']); o = {}
    sets = {'u192': (z['n192'][:, 56:376].astype(float), z['d192'][:, 56:376].astype(float), np.array(PJ['edges']['u192'])[56:377], 56, 376, 2),
            'y': (z['ny'].astype(float), z['dy'].astype(float), np.array(PJ['edges']['y']), 0, 20, 40),
            's': (z['ns'].astype(float), z['ds'].astype(float), np.array(PJ['edges']['s']), 0, 20, 40)}
    grp = np.random.default_rng(7).permutation(n) % 20
    for key, (num, den, e, a0, a1, nsub) in sets.items():
        N = num.sum(0); Dn = den.sum(0); m = N / Dn; mj = (N[None] - num) / (Dn[None] - den); sig = np.sqrt((n - 1) / n * ((mj - mj.mean(0))**2).sum(0)); w = 1 / sig**2
        Y = design_points(e, nsub); dof = len(m) - 8
        r = {}
        for lab, vec in (('data', m), ('pred3', np.array(S['pred_' + key])[a0:a1]), ('pred2', np.array(S['deg2_' + key])[a0:a1])):
            ch, arg = profile(vec, w, Y); r[lab] = summarize(ch, arg)
        if key == 'u192' and r2 in (1, 2):
            al_g = []
            for g in range(20):
                keep = grp != g; mg = num[keep].sum(0) / den[keep].sum(0); ch, arg = profile(mg, w, Y); al_g.append(summarize(ch, arg)['alpha'])
            al_g = np.array(al_g); r['alpha_jk_err'] = float(math.sqrt(19 / 20 * np.sum((al_g - al_g.mean())**2))); r['alpha_groups'] = al_g.tolist()
        o[key] = r
        d = r['data']
        print('r2=%d %-4s (%d bins, dof %d): data alpha %.3f (cycles %.2f, beta %.1f) chi2 min %.1f | chi2 at 1/3 %.1f, at 1/2 %.1f; '
              'calibration: deg-3 prediction gives alpha %.3f, deg-2 comparator %.3f%s  [%.0f s]'
              % (r2, key, len(m), dof, d['alpha'], d['cycles'], d['beta'], d['chi_min'], d['chi_a033'], d['chi_a050'], r['pred3']['alpha'], r['pred2']['alpha'],
                 ('; alpha group-jackknife error %.3f' % r['alpha_jk_err']) if 'alpha_jk_err' in r else '', time.time() - t0), flush=True)
    res[str(r2)] = o
json.dump(res, open(R + 'analysis2.json', 'w'))
