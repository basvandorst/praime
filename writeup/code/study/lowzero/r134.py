"""Round 134 (preregister/round134.txt): the number variance Sigma2(l), l = 0.5 .. 8, of Dirichlet L-function zeros
against the frozen predictions of bk/predictions_per_lfunction.json: the full Bogomolny-Keating law (diagonal and
off-diagonal term, all bad primes) and the diagonal law alone. Per l, the regression of round 123 on the data and on
the predictions; readings (1)-(6) of the registration.
usage: r134.py A | B | C     (A: the 88 complex characters of round 121; B, C: zeros from r134z.py)"""
import os, sys, json, math
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); ZD = os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros')
T = sys.argv[1]
P = json.load(open(os.path.join(H, 'bk', 'predictions_per_lfunction.json')))
LG = P['l']; LS = [0.5 * i for i in range(1, 17)]; IDX = [LG.index(l) for l in LS]
lo, hi = P[T]['height']
FULL = np.array(P[T]['full'])[:, IDX]; DIAG = np.array(P[T]['diag'])[:, IDX]
def sigma2(g, q):
    """Sigma2(l) on the primary grid for one side: g = positive ordinates in [lo, hi]; also the number of windows per l."""
    u = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e)); v, w = [], []
    for l in LS:
        x = np.arange(u[0], u[-1] - l, 0.01); cnt = np.searchsorted(u, x + l) - np.searchsorted(u, x)
        v.append(np.var(cnt)); w.append(len(x))
    return np.array(v), np.array(w)
if T == 'A':
    PR = (2, 3, 5, 7); Z = {}
    for i in range(4):
        z = np.load(os.path.join(ZD, 'r121_part%d.npz' % i)); Z.update({k: z[k] for k in z.files})
    keys = [r['key'] for r in json.load(open(os.path.join(H, 'r121.json')))['rows']]
    Q = [int(k.split('_')[0][1:]) for k in keys]; assert Q == P['A']['N']
    S = []
    for k, q in zip(keys, Q):
        vs, ws = [], []
        for side in (1, -1):
            g = np.sort(side * Z[k][side * Z[k] > 0]); g = g[(g >= lo) & (g <= hi)]; v, w = sigma2(g, q); vs.append(v * w); ws.append(w)
        S.append(sum(vs) / sum(ws))
else:
    PR = (2, 3, 5, 7, 11, 13); Q = P[T]['D']; S = []
    Z = dict(np.load(os.path.join(ZD, 'r134_%s.npz' % T)))
    for D in Q:
        g = np.sort(Z['D%d' % D]); g = g[(g >= lo) & (g <= hi)]; S.append(sigma2(g, abs(D))[0])
S = np.array(S); n = len(Q); nc = len(PR) * len(LS)
X = np.array([[1] + [int(q % p == 0) for p in PR] for q in Q], float)
def coef(X, Y):
    """OLS of every column of Y on X: coefficients of the indicators (primes x lengths) and their OLS standard errors."""
    XtXi = np.linalg.pinv(X.T @ X); c = XtXi @ X.T @ Y; res = Y - X @ c
    s = np.sqrt((res ** 2).sum(0) / (len(Y) - X.shape[1])); return c[1:], np.sqrt(np.diag(XtXi))[1:, None] * s[None, :]
def stats(ix):
    fit, se = coef(X[ix], S[ix]); pf, _ = coef(X[ix], FULL[ix]); pd, _ = coef(X[ix], DIAG[ix])
    off = pf - pd; res = fit - pd
    gam = float((res * off).sum() / (off ** 2).sum())
    y = (S - DIAG)[ix]; x = (FULL - DIAG)[ix]; y = y - y.mean(0); x = x - x.mean(0)
    return fit, se, pf, pd, gam, float((x * y).sum() / (x * x).sum())
all_ = np.arange(n); fit, se, pf, pd, gam, pooled = stats(all_); off = pf - pd; res = fit - pd
rng = np.random.default_rng(134); bg, bp = [], []
for _ in range(2000):
    b = stats(rng.integers(0, n, n)); bg.append(b[4]); bp.append(b[5])
bg = np.array(bg); bp = np.array(bp)
chi_f = float((((fit - pf) / se) ** 2).sum()); chi_d = float((((fit - pd) / se) ** 2).sum())
rho = float((fit[:, 0] * pd[:, 0]).sum() / (pd[:, 0] ** 2).sum()); rho_pred = float((pf[:, 0] * pd[:, 0]).sum() / (pd[:, 0] ** 2).sum())
i2 = [LS.index(l) for l in (6, 6.5, 7, 7.5, 8)]; r2 = res[0, i2]; o2 = off[0, i2]
lev = (S - FULL).mean(0); levd = (S - DIAG).mean(0)
rho_lo, rho_hi = (0.35, 0.75) if T == 'A' else (0.40, 0.72)
R = {'1_full_passes': bool(0.7 <= gam <= 1.3), '1_diag_fails': bool(np.percentile(bg, 2.5) > 0.3),
     '2_full_passes': bool(chi_f <= nc), '2_diag_fails': bool(chi_d >= 4 * chi_f),
     '3_full_passes': bool(rho_lo <= rho <= rho_hi), '3_diag_fails': bool(rho <= 0.80),
     '4_signs': bool(all(np.sign(r2) == np.array([-1, 1, -1, 1, -1]))),
     '5_level': bool(np.abs(lev).max() <= 0.004), '6_pooled': bool(0.7 <= pooled <= 1.3)}
print('test %s: %d L-functions, groups %s' % (T, n, [int(X[:, i + 1].sum()) for i in range(len(PR))]))
print('(1) gamma %.3f, bootstrap se %.3f, 2.5%% point %.3f, 97.5%% point %.3f' % (gam, bg.std(), np.percentile(bg, 2.5), np.percentile(bg, 97.5)))
print('(2) chi2 full %.1f, diagonal %.1f on %d cells' % (chi_f, chi_d, nc))
print('    per prime, full / diagonal: ' + '  '.join('%d: %.1f / %.1f' % (p, (((fit - pf) / se) ** 2)[i].sum(), (((fit - pd) / se) ** 2)[i].sum()) for i, p in enumerate(PR)))
print('(3) rho at l = 0.5: %.3f (full law %.3f, diagonal law 1)' % (rho, rho_pred))
print('(4) p = 2, fit - diagonal at l = 6 .. 8: ' + ' '.join('%+.4f' % v for v in r2) + '   predicted ' + ' '.join('%+.4f' % v for v in o2) + '   se ' + ' '.join('%.4f' % v for v in se[0, i2]))
print('(5) mean Sigma2 - full: max |.| %.4f; mean Sigma2 - diagonal: max |.| %.4f' % (np.abs(lev).max(), np.abs(levd).max()))
print('    per l: ' + ' '.join('%+.4f' % v for v in lev))
print('(6) pooled per-L-function slope %.3f, bootstrap se %.3f' % (pooled, bp.std()))
print('readings', R)
json.dump({'test': T, 'n': n, 'primes': PR, 'l': LS, 'gamma': gam, 'gamma_se': float(bg.std()), 'gamma_q025': float(np.percentile(bg, 2.5)),
           'gamma_q975': float(np.percentile(bg, 97.5)), 'chi2_full': chi_f, 'chi2_diag': chi_d, 'cells': nc, 'rho': rho, 'rho_pred': rho_pred,
           'p2_res': r2.tolist(), 'p2_pred': o2.tolist(), 'p2_se': se[0, i2].tolist(), 'level_full': lev.tolist(), 'level_diag': levd.tolist(),
           'pooled': pooled, 'pooled_se': float(bp.std()), 'fit': fit.tolist(), 'se': se.tolist(), 'pred_full': pf.tolist(), 'pred_diag': pd.tolist(),
           'sigma2': S.tolist(), 'readings': R}, open(os.path.join(H, 'r134_%s.json' % T), 'w'))
