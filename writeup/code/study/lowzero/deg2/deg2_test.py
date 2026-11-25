"""Degree-2 round: number variance of the zeros of elliptic-curve L-functions against the frozen per-curve predictions
(predictions JSON of deg2_predict.py). Self-contained (numpy, scipy).
  1. completeness rule (registered): s(t) = #{0 < gamma <= t} - theta_E(t)/pi + r/2 at 202 heights, its running mean over
     21 heights, the largest unfolded gap; excluded if max|running mean| > 1.0 or max|s| > 2.0 (after the repair of
     deg2_zeros.py); a set flagged only by a gap > 3.0 is kept if the repair searched it (reported).
  2. Sigma2_i(l), l = 0.5 .. 10: u = theta_E(gamma)/pi for the zeros in the band, windows x = arange(u_0, u_end - l, 0.01),
     counts by searchsorted, np.var; edge check with x in [u_0 + 3, u_end - 3 - l].
  3. readings T1-T8 and decisions D1-D7, S_* (see round_draft.txt), bootstrap over curves (numpy default_rng(SEED), 2000 resamples).
usage: deg2_test.py PRED.json ZERODIR OUT.json [--sha SHA256] [--synthetic MODEL|noise] [--seed S]
  --synthetic full|k1|deg1|offzeta|gue: replace the data by MODEL_i + noise drawn from the noise model in
  noise_model.json (dry run / power); --synthetic noise: residuals only (needs zeros)."""
import sys, os, json, math, hashlib, argparse
import numpy as np
from scipy.special import loggamma
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument('pred'); ap.add_argument('zdir'); ap.add_argument('out')
ap.add_argument('--sha', default=None); ap.add_argument('--synthetic', default=None)
ap.add_argument('--noise', default=None); ap.add_argument('--seed', type=int, default=20261125)
ap.add_argument('--nboot', type=int, default=2000)
A = ap.parse_args()
txt = open(A.pred).read(); sha = hashlib.sha256(txt.encode()).hexdigest()
if A.sha and A.sha != sha: sys.exit('predictions file does not match the registered sha256: %s' % sha)
P = json.loads(txt); LS = P['l']; lo, hi = P['band']; CU = P['curves']; GUE = np.array(P['gue'])


def theta_E(N, t):
    t = np.asarray(t, float)
    return t / 2 * math.log(N) - t * math.log(2 * math.pi) + np.imag(loggamma(1 + 1j * t))


def completeness(g, N, r):
    t = np.concatenate([[lo], np.linspace(lo, hi, 200), [hi]])
    s = np.searchsorted(g, t, side='right') - theta_E(N, t) / math.pi + r / 2
    sm = np.convolve(s, np.ones(21) / 21, mode='valid')
    b = g[(g >= lo) & (g <= hi)]; d = np.diff(theta_E(N, b) / math.pi)
    return float(np.abs(s).max()), float(np.abs(sm).max()), float(d.max())


def sigma2(g, N, edge=0.0):
    b = g[(g >= lo) & (g <= hi)]; u = theta_E(N, b) / math.pi; v = []
    for l in LS:
        x = np.arange(u[0] + edge, u[-1] - edge - l, 0.01); c = np.searchsorted(u, x + l) - np.searchsorted(u, x)
        v.append(np.var(c))
    return np.array(v), len(b)


ALTS = ('k1', 'deg1', 'offzeta', 'k1p2', 'k1p3', 'w1p2', 'w1p3')
M = {k: np.array([c[k] for c in CU]) for k in ('full',) + ALTS}
info = []; keep = []
if A.synthetic in (None, 'noise'):
    S = []; SE = []; NZ = []
    for i, c in enumerate(CU):
        f = os.path.join(A.zdir, c['id'] + '.npy')
        g = np.sort(np.load(f)); r = json.load(open(os.path.join(A.zdir, c['id'] + '.log.json')))['r']
        ms, msm, gap = completeness(g, c['N'], r)
        excl = msm > 1.0 or ms > 2.0
        info.append({'id': c['id'], 'max_s': ms, 'max_smooth': msm, 'gap': gap, 'excluded': bool(excl)})
        v, nz = sigma2(g, c['N']); ve, _ = sigma2(g, c['N'], edge=3.0)
        S.append(v); SE.append(ve); NZ.append(nz); keep.append(not excl)
    S = np.array(S); SE = np.array(SE); NZ = np.array(NZ)
else:
    NM = json.load(open(A.noise)); C = np.array(NM['cov']); n0 = NM['zeros_per_curve']
    rng0 = np.random.default_rng(A.seed + 7); NZ = np.array([c.get('nz_expected', n0) for c in CU], float)
    E = rng0.multivariate_normal(np.zeros(len(LS)), C, size=len(CU)) * np.sqrt(n0 / NZ)[:, None]
    S = M[A.synthetic] + E if A.synthetic != 'gue' else GUE[None, :] + E
    SE = S.copy(); keep = [True] * len(CU)
keep = np.array(keep); idx = np.nonzero(keep)[0]
S = S[idx]; SE = SE[idx]; M = {k: v[idx] for k, v in M.items()}; CU = [CU[i] for i in idx]; n, nl = S.shape
if A.synthetic == 'noise':                       # dry run in noise-only mode: data = law + measured residual (permuted)
    rng0 = np.random.default_rng(A.seed + 11); R = S - M['full']; S = M['full'] + R[rng0.permutation(n)]

# ------------------------------------------------------------ design for the per-prime contrasts
G2 = ['0', '2', 'b']; G3 = ['0', '2', '3', 'b']                      # baselines |a_2| = 1, |a_3| = 1
X = np.array([[1.0] + [float(c['g2'] == g) for g in G2] + [float(c['g3'] == g) for g in G3] for c in CU])
XtXi = np.linalg.pinv(X.T @ X)


def ols(Y):
    c = XtXi @ X.T @ Y; return c, Y - X @ c


def readings(ix):
    Sx = S[ix]; Mx = {k: v[ix] for k, v in M.items()}; Xx = X[ix]; Xi = np.linalg.pinv(Xx.T @ Xx)
    out = {}
    res = Sx - Mx['full']; out['level'] = res.mean(0)
    w = 1 / np.maximum(res.var(0, ddof=1), 1e-12)
    def slope(y, x, centre=True):
        if centre: y = y - y.mean(0); x = x - x.mean(0)
        return float((w * (x * y).sum(0)).sum() / (w * (x * x).sum(0)).sum())
    out['gamma_A'] = slope(Sx, Mx['full'])
    for k in ALTS:
        out['gamma_' + k] = slope(Sx - Mx[k], Mx['full'] - Mx[k])
    out['lessrigid'] = (Sx[LR[ix]] - GUE).mean(0) if LR[ix].any() else np.full(nl, np.nan)
    out['gamma_gue_level'] = slope(Sx - GUE, Mx['full'] - GUE, centre=False)
    cf = Xi @ Xx.T @ Sx; out['coef'] = cf[1:]
    for k in ('full', 'k1', 'deg1', 'offzeta'):
        out['pc_' + k] = (Xi @ Xx.T @ Mx[k])[1:]
    return out


LR = (M['full'][:, LS.index(6.0)] - GUE[LS.index(6.0)]) > 0.02      # curves predicted less rigid than GUE (frozen)
R0 = readings(np.arange(n))
rng = np.random.default_rng(A.seed); B = [readings(rng.integers(0, n, n)) for _ in range(A.nboot)]
def bse(key): return np.nanstd(np.array([b[key] for b in B], dtype=float), axis=0, ddof=1)
# Hotelling T^2 of the level
D = S - M['full']; m = D.mean(0); Cd = np.cov(D.T)
if n > nl + 2:
    T2 = float(n * m @ np.linalg.solve(Cd, m)); F = T2 * (n - nl) / ((n - 1) * nl); pT2 = float(stats.f.sf(F, nl, n - nl))
else:
    T2 = float('nan'); pT2 = float('nan')                  # too few curves (development set): not defined
# chi2 of the per-prime coefficients against each hypothesis, covariance = (X'X)^-1_gg' (x) C_res
cf, res = ols(S); Cres = np.cov(res.T, ddof=X.shape[1]) * (n - 1) / (n - X.shape[1])
def chi2(rows, key):
    d = (R0['coef'][rows] - R0['pc_' + key][rows]).ravel() if key != 'gue' else R0['coef'][rows].ravel()
    Cr = Cres if n - X.shape[1] > nl + 2 else np.diag(np.diag(Cres))   # development set: diagonal only
    V = np.kron(XtXi[1:, 1:][np.ix_(rows, rows)], Cr); return float(d @ np.linalg.solve(V, d))
p2 = [0, 1, 2]; p3 = [3, 4, 5, 6]
chi = {f'chi2_p{p}_{k}': chi2(rows, k) for p, rows in ((2, p2), (3, p3)) for k in ('full', 'k1', 'deg1', 'offzeta', 'gue')}
sec = {k: (float(R0[k]), float(bse(k))) for k in ['gamma_A'] + ['gamma_' + a for a in ALTS] + ['gamma_gue_level']}
lr_pred = (M['full'][LR] - GUE).mean(0) if LR.any() else np.full(nl, np.nan)
# the curves are fixed by the frozen predictions, so the error of their mean is the measurement noise only
lr_se = (S[LR] - M['full'][LR]).std(0, ddof=1) / math.sqrt(max(LR.sum(), 1)) if LR.sum() > 1 else np.full(nl, np.nan)
i6 = LS.index(6.0); lr_noise6 = 0.0042                 # registered: per-curve noise at l = 6 (development residuals)
lr_count = int(((S[LR, i6] - GUE[i6]) > 3 * lr_noise6).sum()); lr_count_pred = int(((M['full'][LR, i6] - GUE[i6]) > 3 * lr_noise6).sum())
lev_se = bse('level'); coef_se = bse('coef')
i_less = [LS.index(l) for l in (6.0, 7.0, 8.0, 9.0, 10.0)]
out = {'n': n, 'excluded': [x for x in info if x['excluded']], 'completeness': info, 'zeros': NZ.tolist() if hasattr(NZ, 'tolist') else NZ,
       'sha256': sha, 'l': LS, 'level': R0['level'].tolist(), 'level_se': lev_se.tolist(), 'T2': T2, 'p_T2': pT2,
       'gammas': sec, 'chi2': chi, 'cells_p2': 3 * nl, 'cells_p3': 4 * nl,
       'coef': R0['coef'].tolist(), 'coef_se': coef_se.tolist(), 'pred_coef': {k: R0['pc_' + k].tolist() for k in ('full', 'k1', 'deg1', 'offzeta')},
       'groups': {'g2': {g: int(sum(c['g2'] == g for c in CU)) for g in ('0', '1', '2', 'b')},
                  'g3': {g: int(sum(c['g3'] == g for c in CU)) for g in ('0', '1', '2', '3', 'b')}},
       'sigma2': S.tolist(), 'edge_level': (SE - M['full']).mean(0).tolist(), 'synthetic': A.synthetic,
       'lessrigid': {'n': int(LR.sum()), 'measured': R0['lessrigid'].tolist(), 'se': lr_se.tolist(), 'law': lr_pred.tolist(),
                     'count_above_gue_3sd_l6': lr_count, 'count_law': lr_count_pred, 'noise_sd_l6': float(lr_noise6)}}
# registered decisions (thresholds in round_draft.txt)
g = sec
crit = {'chi2_0.001_60': float(stats.chi2.ppf(0.999, 3 * nl)), 'chi2_0.001_80': float(stats.chi2.ppf(0.999, 4 * nl))}
def passes(k): return bool(abs(g[k][0] - 1) <= 4 * g[k][1])
def rejects(k): return bool(g[k][0] >= 5 * g[k][1] and g[k][0] >= 0.5)
lv_ok = bool(pT2 >= 0.001 and np.all(np.abs(R0['level']) <= 4 * lev_se))
i4, i6 = LS.index(4.0), LS.index(6.0)
dec = {
    'D1_level_law': lv_ok,
    'D1_level_offset_within_allowance': bool((not lv_ok) and np.abs(R0['level']).max() <= 0.002),
    'D2_weights_law': passes('gamma_deg1') and passes('gamma_A'),
    'D2_deg1_rejected': rejects('gamma_deg1'),
    'D3_doubled_orbit_law': passes('gamma_k1p2') and passes('gamma_k1'),
    'D3_k1p2_rejected': rejects('gamma_k1p2'),
    'D3_k1_rejected': rejects('gamma_k1'),
    'D4_first_traversal_law': passes('gamma_w1p2'),
    'D4_w1p2_rejected': rejects('gamma_w1p2'),
    'D5_p2_shape_law': bool(chi['chi2_p2_full'] <= crit['chi2_0.001_60']),
    'D5_p2_gue_rejected': bool(chi['chi2_p2_gue'] > crit['chi2_0.001_60']),
    'D5_p2_deg1_rejected': bool(chi['chi2_p2_deg1'] > crit['chi2_0.001_60']),
    'D6_less_rigid_than_GUE': bool(LR.sum() >= 5 and all(R0['lessrigid'][i] >= 5 * lr_se[i] and abs(R0['lessrigid'][i] - lr_pred[i]) <= 4 * lr_se[i]
                                                        for i in (i4, i6))),
    'D7_sign_a2_2': bool(np.all(R0['coef'][1, i_less] > 0)),
    'S_k1p3_law': passes('gamma_k1p3'), 'S_k1p3_rejected': rejects('gamma_k1p3'),
    'S_w1p3_law': passes('gamma_w1p3'), 'S_w1p3_rejected': rejects('gamma_w1p3'),
    'S_offzeta_law': passes('gamma_offzeta'), 'S_offzeta_rejected': rejects('gamma_offzeta'),
    'S_p3_shape_law': bool(chi['chi2_p3_full'] <= crit['chi2_0.001_80'])}
out['decisions'] = dec; out['crit'] = crit
json.dump(out, open(A.out, 'w'))
print('n = %d curves (excluded %d), zeros in band %d (mean %.0f), sha256 %s' % (n, len(out['excluded']), int(np.sum(NZ)), float(np.mean(NZ)), sha[:16]))
print('groups', out['groups'])
print('T1 level, data - law, l = 0.5..10:'); print('   ' + ' '.join('%+.4f' % v for v in R0['level']))
print('   se ' + ' '.join('%.4f' % v for v in lev_se)); print('   Hotelling T2 %.1f on %d, p = %.3g' % (T2, nl, pT2))
for k, (v, s) in sec.items(): print('%-18s %.4f +- %.4f' % (k, v, s))
for k, v in chi.items(): print('%-20s %.1f' % (k, v))
print('p = 2 coefficients (|a_2| = 0, 2, 2|N, each minus |a_2| = 1) at l = 2, 4, 6, 8, 10: fit / law / k1 / deg1')
for j, gname in enumerate(G2):
    ii = [LS.index(l) for l in (2.0, 4.0, 6.0, 8.0, 10.0)]
    print('  %-3s fit  ' % gname + ' '.join('%+.4f' % R0['coef'][j, i] for i in ii) + '   se ' + ' '.join('%.4f' % coef_se[j, i] for i in ii))
    for k in ('full', 'k1', 'deg1'):
        print('      %-5s' % k + ' '.join('%+.4f' % R0['pc_' + k][j, i] for i in ii))
print('less rigid than GUE: %d curves predicted above GUE by > 0.02 at l = 6; mean (Sigma2 - GUE) at l = 2, 4, 6, 8, 10:' % LR.sum())
ii = [LS.index(l) for l in (2.0, 4.0, 6.0, 8.0, 10.0)]
print('   measured ' + ' '.join('%+.4f' % R0['lessrigid'][i] for i in ii) + '   se ' + ' '.join('%.4f' % lr_se[i] for i in ii))
print('   law      ' + ' '.join('%+.4f' % lr_pred[i] for i in ii))
print('   curves with Sigma2(6) - GUE(6) > 3 x per-curve noise (%.4f): %d of %d (law: %d)' % (lr_noise6, lr_count, LR.sum(), lr_count_pred))
print('edge check (windows 3 spacings from the band ends), data - law: ' + ' '.join('%+.4f' % v for v in out['edge_level']))
print('decisions', dec)
