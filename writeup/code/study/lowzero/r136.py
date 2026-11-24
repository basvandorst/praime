"""Round 136 (preregister/round136.txt): is the amplitude of the bad-prime off-diagonal term 1 (Bogomolny-Keating) or
0.932 (round 134, test B)? The statistic of r134.py on a fresh test, against the frozen predictions of
amp/predictions_alt.json; readings R1-R5, the level reading and the completeness of the zero sets.
usage: r136.py C | Dn | Bp     (zeros: lzeros/r134_C.npz, lzeros/r136_Dn.npz, lzeros/r136_Bp.npz)"""
import os, sys, json, math
import numpy as np
from scipy.special import loggamma
H = os.path.dirname(os.path.abspath(__file__)); ZD = os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros'); T = sys.argv[1]
P = json.load(open(os.path.join(H, 'amp', 'predictions_alt.json'))); LS = P['l']; lo, hi = P[T]['height']; Q = P[T]['D']
F = np.array(P[T]['full']); Dg = np.array(P[T]['diag']); FL = np.array(P[T]['full_LP3']); PR = (2, 3, 5, 7, 11, 13)
Z = dict(np.load(sys.argv[2] if len(sys.argv) > 2 else os.path.join(ZD, ('r134_C.npz' if T == 'C' else 'r136_%s.npz' % T))))   # argv[2]: a synthetic set, for dry runs
def sigma2(g, q):
    u = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e)); v = []
    for l in LS:
        x = np.arange(u[0], u[-1] - l, 0.01); v.append(np.var(np.searchsorted(u, x + l) - np.searchsorted(u, x)))
    return v
def theta(t, q, a): return np.imag(loggamma(0.25 + a / 2 + 0.5j * t)) + t / 2 * math.log(q / math.pi)
def smean(g, D):
    """Mean over the band of N(t) - theta_chi(t)/pi (trapezoid on a grid of 4000 points), N counted from height 0."""
    q = abs(D); a = 0 if D > 0 else 1; t = np.linspace(lo, hi, 4000)
    return float(np.mean(np.searchsorted(g, t) - (theta(t, q, a) - theta(0.0, q, a)) / math.pi))
S, comp = [], []
for D in Q:
    g = np.sort(Z['D%d' % D]); comp.append(smean(g, D)); g = g[(g >= lo) & (g <= hi)]; S.append(sigma2(g, abs(D)))
S = np.array(S); comp = np.array(comp); n = len(Q); nc = len(PR) * len(LS)
X = np.array([[1] + [int(q % p == 0) for p in PR] for q in Q], float)
def stats(Sx, ix=None, rows=slice(None), cols=slice(None), law=None):
    ix = np.arange(n) if ix is None else ix; Xi = X[ix]; A = np.linalg.pinv(Xi.T @ Xi) @ Xi.T; Fx = F if law is None else law
    fit = (A @ Sx[ix])[1:]; pf = (A @ Fx[ix])[1:]; pd = (A @ Dg[ix])[1:]; off = (pf - pd)[rows, cols]; res = (fit - pd)[rows, cols]
    r1 = float((res * off).sum() / (off ** 2).sum())
    y = (Sx - Dg)[ix]; x = (Fx - Dg)[ix]; yc = y - y.mean(0); xc = x - x.mean(0); r2 = float((xc * yc).sum() / (xc * xc).sum())
    Hm = Xi @ A; hy = Hm @ y; hx = Hm @ x; by = hy - hy.mean(0); bx = hx - hx.mean(0); wy = y - hy; wx = x - hx
    return r1, r2, float((bx * by).sum() / (bx * bx).sum()), float((wx * wy).sum() / (wx * wx).sum()), fit, pf, pd
r1, r2, bet, wit, fit, pf, pd = stats(S)
rng = np.random.default_rng(136); bs = np.array([stats(S, ix=rng.integers(0, n, n))[:4] for _ in range(2000)])
A0 = np.linalg.pinv(X.T @ X) @ X.T; resF = (S - F) - X @ (A0 @ (S - F)); rng2 = np.random.default_rng(1361)
nul = np.array([stats(F + resF[rng2.integers(0, n, n)])[:4] for _ in range(1000)]); dbw = nul[:, 2] - nul[:, 3]
rs = S - X @ (A0 @ S); se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X)))[1:, None] * np.sqrt((rs ** 2).sum(0) / (n - X.shape[1]))[None, :]
rho = float((fit[:, 0] * pd[:, 0]).sum() / (pd[:, 0] ** 2).sum()); i2 = [LS.index(l) for l in (6, 6.5, 7, 7.5, 8)]
def split(mask=None, rows=slice(None), cols=slice(None)):
    ix0 = np.arange(n) if mask is None else np.nonzero(mask)[0]; r = np.random.default_rng(136)
    v = stats(S, ix=ix0, rows=rows, cols=cols)[0]; b = [stats(S, ix=ix0[r.integers(0, len(ix0), len(ix0))], rows=rows, cols=cols)[0] for _ in range(500)]
    return [v, float(np.std(b))]
Dq = np.array(Q); R5 = {'D>0': split(Dq > 0), 'D<0': split(Dq < 0), 'l<=4': split(cols=[j for j, l in enumerate(LS) if l <= 4]),
                      'l>4': split(cols=[j for j, l in enumerate(LS) if l > 4]), 'p=2': split(rows=slice(0, 1))}
def level(law): d = S - law; m = d.mean(0); return float(m @ np.linalg.pinv(np.cov(d.T) / n) @ m), float(np.abs(m).max())
lv = {'full': level(F), 'diag': level(Dg), 'LP3': level(FL), 'AMP': level(Dg + 0.932 * (F - Dg))}
lp3 = stats(S, law=FL)
print('test %s: %d L-functions, heights %d-%d, %d zeros in the band; completeness: mean of N - theta/pi over the band %.4f +- %.4f, largest |.| %.3f' % (T, n, lo, hi, sum(int(((np.sort(Z['D%d' % D]) >= lo) & (np.sort(Z['D%d' % D]) <= hi)).sum()) for D in Q), comp.mean(), comp.std() / math.sqrt(n), np.abs(comp).max()))
print('R1 gamma   %.4f +- %.4f (95%% %.3f-%.3f)' % (r1, bs[:, 0].std(), np.percentile(bs[:, 0], 2.5), np.percentile(bs[:, 0], 97.5)))
print('R2 pooled  %.4f +- %.4f   (decision line 0.966)' % (r2, bs[:, 1].std()))
print('R3 between %.4f +- %.4f, within %.4f +- %.4f; difference %.4f, null sd %.4f (%.1f sd)' % (bet, bs[:, 2].std(), wit, bs[:, 3].std(), bet - wit, dbw.std(), (bet - wit) / dbw.std()))
print('   null (full law true): gamma %.3f +- %.4f, pooled %.3f +- %.4f' % (nul[:, 0].mean(), nul[:, 0].std(), nul[:, 1].mean(), nul[:, 1].std()))
print('R4 rho(0.5) %.3f; p = 2 at l = 6..8: ' % rho + ' '.join('%+.4f' % v for v in (fit - pd)[0, i2]) + '  BK ' + ' '.join('%+.4f' % v for v in (pf - pd)[0, i2]) + '  se ' + ' '.join('%.4f' % v for v in se[0, i2]))
print('R5 ' + '  '.join('%s %.3f +- %.3f' % (k, v[0], v[1]) for k, v in R5.items()))
print('level chi2 on 16 lengths (max |mean difference|): ' + '  '.join('%s %.1f (%.4f)' % (k, v[0], v[1]) for k, v in lv.items()))
print('against the off-diagonal part at L + 0.3: gamma %.4f, pooled %.4f' % (lp3[0], lp3[1]))
json.dump({'test': T, 'n': n, 'R1': r1, 'R1_se': float(bs[:, 0].std()), 'R1_q': [float(np.percentile(bs[:, 0], 2.5)), float(np.percentile(bs[:, 0], 97.5))], 'R2': r2, 'R2_se': float(bs[:, 1].std()),
           'between': bet, 'between_se': float(bs[:, 2].std()), 'within': wit, 'within_se': float(bs[:, 3].std()), 'null_sd_between_minus_within': float(dbw.std()),
           'null_gamma': [float(nul[:, 0].mean()), float(nul[:, 0].std())], 'null_pooled': [float(nul[:, 1].mean()), float(nul[:, 1].std())], 'rho': rho,
           'p2_res': (fit - pd)[0, i2].tolist(), 'p2_BK': (pf - pd)[0, i2].tolist(), 'p2_se': se[0, i2].tolist(), 'R5': R5, 'level': lv, 'LP3': [lp3[0], lp3[1]],
           'completeness_mean': float(comp.mean()), 'completeness_max': float(np.abs(comp).max()), 'sigma2': S.tolist()}, open(os.path.join(H, 'r136_%s.json' % T), 'w'))
