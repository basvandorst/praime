"""Round 134, calibrated reading of the registered data (after the forty-second review; not pre-registered).
The registered chi2 used the OLS errors of Sigma2 itself and treated the 16 lengths as independent. Here:
 (a) GLS chi2 of the fitted per-prime coefficients against each law, with the covariance
     (residual covariance of Sigma2 - law across the lengths) x (X'X)^-1, and its null distribution;
 (b) the amplitude gamma of the off-diagonal term with a null simulation: data = full law + resampled residual rows;
 (c) the absolute level: chi2 of mean(Sigma2 - law) over the 16 lengths with the covariance of the mean;
 (d) gamma split by prime, by length and by sign of D.
usage: r134b.py A | B"""
import os, sys, json, math
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); T = sys.argv[1]
P = json.load(open(os.path.join(H, 'bk', 'predictions_per_lfunction.json'))); R = json.load(open(os.path.join(H, 'r134_%s.json' % T)))
LS = R['l']; IDX = [P['l'].index(l) for l in LS]; PR = R['primes']
Q = P[T]['N'] if T == 'A' else P[T]['D']
S = np.array(R['sigma2']); F = np.array(P[T]['full'])[:, IDX]; Dg = np.array(P[T]['diag'])[:, IDX]
X = np.array([[1] + [int(q % p == 0) for p in PR] for q in Q], float); n, k = X.shape
XtXi = np.linalg.inv(X.T @ X); A = XtXi @ X.T
def gls(law):
    c = (A @ (S - law))[1:]; res = (S - law) - X @ (A @ (S - law)); C = res.T @ res / (n - k)
    V = np.kron(XtXi[1:, 1:], C); v = c.ravel()                        # order: prime, then length
    return float(v @ np.linalg.pinv(V) @ v), res
def gamma(Sx, rows=slice(None), cols=slice(None), ix=None):
    ix = np.arange(n) if ix is None else ix; Ai = np.linalg.pinv(X[ix].T @ X[ix]) @ X[ix].T
    fit = (Ai @ Sx[ix])[1:]; pf = (Ai @ F[ix])[1:]; pd = (Ai @ Dg[ix])[1:]; off = (pf - pd)[rows, cols]; res = (fit - pd)[rows, cols]
    return float((res * off).sum() / (off ** 2).sum())
cf, resF = gls(F); cd, resD = gls(Dg)
rng = np.random.default_rng(1342); gam = gamma(S)
boot = np.array([gamma(S, ix=rng.integers(0, n, n)) for _ in range(2000)])
null_g, null_c, null_reg = [], [], []
fit0, se0 = np.array(R['fit']), np.array(R['se'])
for _ in range(1000):
    Sn = F + resF[rng.integers(0, n, n)]                                # the full law is true; residual rows resampled
    null_g.append(gamma(Sn))
    c = (A @ (Sn - F))[1:]; r = (Sn - F) - X @ (A @ (Sn - F)); V = np.kron(XtXi[1:, 1:], r.T @ r / (n - k)); null_c.append(float(c.ravel() @ np.linalg.pinv(V) @ c.ravel()))
    rs = Sn - X @ (A @ Sn); s = np.sqrt((rs ** 2).sum(0) / (n - k)); se = np.sqrt(np.diag(XtXi))[1:, None] * s[None, :]
    null_reg.append(float(((((A @ Sn)[1:] - (A @ F)[1:]) / se) ** 2).sum()))
null_g = np.array(null_g); null_c = np.array(null_c); null_reg = np.array(null_reg)
pval = float((np.abs(null_g - null_g.mean()) >= abs(gam - null_g.mean())).mean())
def level(law):
    d = S - law; m = d.mean(0); C = np.cov(d.T) / n; return float(m @ np.linalg.pinv(C) @ m)
lf, ld = level(F), level(Dg)
ev = np.linalg.eigvalsh(np.corrcoef(resF.T)); neff = float(ev.sum() ** 2 / (ev ** 2).sum())
split = {}
for i, p in enumerate(PR[:3]):
    b = np.array([gamma(S, rows=slice(i, i + 1), ix=rng.integers(0, n, n)) for _ in range(500)]); split['p=%d' % p] = [gamma(S, rows=slice(i, i + 1)), float(b.std())]
for nm, cols in (('l<=4', [j for j, l in enumerate(LS) if l <= 4]), ('l>4', [j for j, l in enumerate(LS) if l > 4]), ('whole l', [j for j, l in enumerate(LS) if l == int(l)]), ('half l', [j for j, l in enumerate(LS) if l != int(l)])):
    b = np.array([gamma(S, cols=cols, ix=rng.integers(0, n, n)) for _ in range(500)]); split[nm] = [gamma(S, cols=cols), float(b.std())]
print('test %s: %d L-functions, %d cells' % (T, n, len(PR) * len(LS)))
print('(a) GLS chi2: full %.1f, diagonal %.1f; null of the full law: mean %.1f, 95%% point %.1f, 99%% point %.1f' % (cf, cd, null_c.mean(), np.percentile(null_c, 95), np.percentile(null_c, 99)))
print('    registered chi2 under the full law: null mean %.1f, 95%% point %.1f, 99%% point %.1f (registered threshold %d); effective number of independent lengths %.1f of %d' % (null_reg.mean(), np.percentile(null_reg, 95), np.percentile(null_reg, 99), len(PR) * len(LS), neff, len(LS)))
print('(b) gamma %.4f, bootstrap se %.4f, 95%% interval %.3f-%.3f; null (full law true): %.3f +- %.4f, sets at or below the observed value %d of 1000, two-sided p %.3f, z %.2f' % (gam, boot.std(), np.percentile(boot, 2.5), np.percentile(boot, 97.5), null_g.mean(), null_g.std(), int((null_g <= gam).sum()), pval, (gam - null_g.mean()) / null_g.std()))
print('(c) level chi2 on %d lengths: full %.1f, diagonal %.1f' % (len(LS), lf, ld))
print('(d) gamma by split: ' + '  '.join('%s %.3f +- %.3f' % (kk, v[0], v[1]) for kk, v in split.items()))
json.dump({'test': T, 'gls_full': cf, 'gls_diag': cd, 'gls_null_mean': float(null_c.mean()), 'gls_null_q95': float(np.percentile(null_c, 95)), 'gls_null_q99': float(np.percentile(null_c, 99)),
           'reg_null_mean': float(null_reg.mean()), 'reg_null_q95': float(np.percentile(null_reg, 95)), 'reg_null_q99': float(np.percentile(null_reg, 99)), 'neff': neff,
           'gamma': gam, 'gamma_se': float(boot.std()), 'gamma_q025': float(np.percentile(boot, 2.5)), 'gamma_q975': float(np.percentile(boot, 97.5)),
           'null_gamma_mean': float(null_g.mean()), 'null_gamma_sd': float(null_g.std()), 'null_below': int((null_g <= gam).sum()), 'p_two_sided': pval,
           'level_full': lf, 'level_diag': ld, 'split': split}, open(os.path.join(H, 'r134b_%s.json' % T), 'w'), indent=0)
