# review 162: statistics on own counts. (A) registered bins, y bins, sqrt(y) bins: amplitude, correlation, delta chi2 deg-2 vs deg-3
# (jackknife over fields); (C) conditional amplitudes vs cut-off L0; (B) generic oscillation fit with free alpha on fine bins.
import sys, json, math, numpy as np
R = '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/review_162/'
PJ = json.load(open(R + 'pred_own_T1000.json'))
ylo, yhi = (7 / 24)**3, (47 / 24)**3
def jk(num, den):
    n = len(num); N = num.sum(0); Dn = den.sum(0); m = N / Dn; mj = (N[None] - num) / (Dn[None] - den)
    dm = mj - mj.mean(0); return m, (n - 1) / n * dm.T @ dm, mj
def jkerr(vals): n = len(vals); return math.sqrt((n - 1) / n * np.sum((vals - vals.mean())**2))
def corr_rows(M, P):   # correlation of each row of M with P
    Mc = M - M.mean(1, keepdims=True); Pc = P - P.mean(); return (Mc @ Pc) / np.sqrt((Mc**2).sum(1) * (Pc**2).sum())
def gls(m, Ci, P):
    a = float(P @ Ci @ m / (P @ Ci @ P)); r = m - a * P; return a, 1 / math.sqrt(P @ Ci @ P), float(r @ Ci @ r)
res = {}
for r2 in (1, 2, 0):
    z = np.load(R + 'mine_r%d.npz' % r2); S = PJ['sig'][str(r2)]; n = len(z['D']); out = {}
    n192, d192 = z['n192'].astype(float), z['d192'].astype(float)
    # registered bins from the 1/192 bins: edges (2j+1)/24 = (16j+8)/192; HI = bins 4..23 of the coarse grid
    fi = [0] + [16 * j + 8 for j in range(24)]
    numR = np.stack([n192[:, fi[k]:fi[k + 1]].sum(1) for k in range(24)], 1); denR = np.stack([d192[:, fi[k]:fi[k + 1]].sum(1) for k in range(24)], 1)
    sets = {'reg': (numR[:, 4:24], denR[:, 4:24], 4, 24), 'y': (z['ny'].astype(float), z['dy'].astype(float), 0, 20), 's': (z['ns'].astype(float), z['ds'].astype(float), 0, 20)}
    for key, (num, den, a0, a1) in sets.items():
        m, C, mj = jk(num, den); sig = np.sqrt(np.diag(C)); w = 1 / sig**2; Ci = np.linalg.inv(C) if n >= 200 else np.diag(w)
        P3 = np.array(S['pred_' + key])[a0:a1]; P2 = np.array(S['deg2_' + key])[a0:a1]; PL = np.array(S['line_' + key])[a0:a1]
        a3, a3e, c3 = gls(m, Ci, P3); a2, a2e, c2 = gls(m, Ci, P2); aL, _, cL = gls(m, Ci, PL)
        sd = float(np.sum(P3 * m * w) / np.sum(P3**2 * w)); sdj = (mj * (P3 * w)[None]).sum(1) / np.sum(P3**2 * w)
        cd3 = float(np.sum((m - sd * P3)**2 * w)); sd2 = float(np.sum(P2 * m * w) / np.sum(P2**2 * w)); cd2 = float(np.sum((m - sd2 * P2)**2 * w))
        cr3 = float(np.corrcoef(m, P3)[0, 1]); cr2 = float(np.corrcoef(m, P2)[0, 1]); cj3 = corr_rows(mj, P3); cj2 = corr_rows(mj, P2)
        out[key] = dict(gls3=a3, gls3e=a3e, chi3=c3, gls2=a2, chi2=c2, dchi_deg2=c2 - c3, dchi_line=cL - c3, diag3=sd, diag3e=jkerr(sdj),
                        chi3_diag=cd3, chi2_diag=cd2, dchi_deg2_diag=cd2 - cd3, corr3=cr3, corr3e=jkerr(cj3), corr2=cr2, corr2e=jkerr(cj2),
                        zero=float(m @ Ci @ m) - c3)
        print('r2=%d %-3s: corr deg3 %.3f+-%.3f, deg2 %.3f+-%.3f | GLS amp %.3f+-%.3f chi2 %.1f; dchi2(deg2-deg3) %.0f (diag cov %.0f), (line-deg3) %.0f, (zero-deg3) %.0f | diag amp %.3f+-%.3f'
              % (r2, key, cr3, jkerr(cj3), cr2, jkerr(cj2), a3, a3e, c3, c2 - c3, cd2 - cd3, cL - c3, float(m @ Ci @ m) - c3, sd, jkerr(sdj)), flush=True)
        if key == 'reg':
            out['m_reg'] = m.tolist(); out['sig_reg'] = sig.tolist()
            # (C) conditional amplitudes with observed local types at l <= L0
            for L0 in (5, 7, 13, 29, 47):
                Pc = np.array(S['cond%d_reg' % L0])[4:24]; s_ = float(np.sum(Pc * m * w) / np.sum(Pc**2 * w)); sj = (mj * (Pc * w)[None]).sum(1) / np.sum(Pc**2 * w)
                ag, age, cg = gls(m, Ci, Pc); out['cond%d' % L0] = dict(diag=s_, diage=jkerr(sj), gls=ag, glse=age, chi=cg, corr=float(np.corrcoef(m, Pc)[0, 1]))
                print('      cond L0=%2d: diag amp %.3f+-%.3f, GLS %.3f+-%.3f (chi2 %.1f), corr %.3f' % (L0, s_, jkerr(sj), ag, age, cg, np.corrcoef(m, Pc)[0, 1]))
    res[str(r2)] = out
json.dump(res, open(R + 'analysis1.json', 'w'), indent=0)
