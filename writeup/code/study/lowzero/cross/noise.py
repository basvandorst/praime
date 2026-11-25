"""NOISE ONLY on the already-used real zero sets (R120 + B + Dn, 546 L-functions, heights 100-600).
Prints: completeness; per-function variance of dn (single-function level, known); the dispersion of the pair
cross-covariances AFTER removing the regression on the predictors (residual sd); jackknife (delete-one-L-function)
standard errors of the registered amplitude readings for the old x old design and for subsets that mimic the fresh
designs (186 x 186, 186 x 360); off-line single-function projections (frequencies ln y, y not a prime power) as the
noise of the per-function orbit amplitude. NO coefficient, NO amplitude value, NO comparison with a prediction is printed.
usage: noise.py -> noise_out.txt; saves cache_C.npz (pair covariances; not to be read for coefficients)"""
import os, sys, json, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from crosslib import *
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..')); LOW = os.path.join(REPO, 'writeup/code/study/lowzero'); ZD = os.path.join(REPO, 'writeup/data/riemann/lzeros')
PJ = json.load(open(os.path.join(LOW, 'triple', 'predictions_kappa3.json')))
FILES = {'R120': ['r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz'], 'B': ['r134_B.npz'], 'Dn': ['r136_Dn.npz']}
T1, T2, H = 100.0, 600.0, 0.02
W = 0.3 * np.arange(1, 17)
Ds, G = [], []
for s in ('R120', 'B', 'Dn'):
    Z = {}
    for f in FILES[s]:
        z = np.load(os.path.join(ZD, f)); Z.update({k: z[k] for k in z.files})
    for D in PJ[s]['D']: Ds.append(D); G.append(np.sort(Z['D%d' % D]))
N = len(Ds); print('%d L-functions' % N)
# completeness (round 138 rule)
bad = []
for D, g in zip(Ds, G):
    q = abs(D); a = 0 if D > 0 else 1; t = np.linspace(T1, T2, 200)
    s = np.searchsorted(g, t) - (theta(t, q, a) - theta(0.0, q, a)) / math.pi
    b = g[(g >= T1) & (g <= T2)]; u = b / (2 * math.pi) * np.log(q * b / (2 * math.pi * math.e))
    if np.abs(s).max() > 1.5 or np.diff(u).max() > 3.0: bad.append((D, float(np.abs(s).max()), float(np.diff(u).max())))
print('completeness: %d flagged %s' % (len(bad), bad))
# pair covariances
CACHE = os.path.join(HERE, 'cache_C.npz')
if os.path.exists(CACHE):
    z = np.load(CACHE); C = z['C']; OFF = z['OFF']; YOFF = z['YOFF']
else:
    C = np.zeros((len(W), N, N)); YOFF = np.array([1.6, 2.5, 3.5, 4.5, 6.0, 6.5, 10.0, 12.0, 14.0, 15.0])
    OFF = np.zeros((len(W), N, len(YOFF)), complex); t0 = time.time()
    for k, w in enumerate(W):
        M = int(math.floor((T2 - w - T1) / H + 1e-9)) + 1; tg = T1 + H * np.arange(M)
        Dl = np.empty((N, M))
        for i, (D, g) in enumerate(zip(Ds, G)): Dl[i] = dn_series(g, abs(D), 0 if D > 0 else 1, w, tg)
        C[k] = Dl @ Dl.T / M
        E = np.exp(-1j * np.outer(tg + w / 2, np.log(YOFF)))          # M x len(YOFF)
        OFF[k] = 2 * (Dl @ E) / M
        print('  w=%.1f done (%.0fs)' % (w, time.time() - t0), flush=True)
    np.savez(CACHE, C=C, OFF=OFF, YOFF=YOFF, D=np.array(Ds), W=W)
# single-function variance (known level; for the correlation coefficient)
V = np.array([np.diag(C[k]) for k in range(len(W))])
print('mean variance of dn over the 546 sets, per w: ' + ' '.join('%.3f' % v for v in V.mean(1)))
# off-line projections: relative noise of a would-be orbit at y (amplitude scale (2/pi) y^-1/2 |sin(w ln y/2)|)
print('off-line projections (y not a prime power): rms of the complex projection / line scale (2/pi) y^-1/2 |sin(w ln y/2)|, sqrt(mean over functions of |a|^2/2)')
for j, y in enumerate(YOFF):
    sc = 2 / math.pi / math.sqrt(y) * np.abs(np.sin(W * math.log(y) / 2))
    rel = np.sqrt((np.abs(OFF[:, :, j]) ** 2).mean(1) / 2) / sc
    print('   y=%5.1f ' % y + ' '.join('%.3f' % r for r in rel))
# design
P = primes_upto(200000); X = kron_table(Ds, P)
iu = np.triu_indices(N, 1); I, J = iu
PR = (2, 3, 5, 7)
def AT(p, w, M_tg):
    """finite-T exact diagonal coefficients on the actual grid: (A, B) with the factor (1 + <cos 2m(t+w/2) ln p>)"""
    lp = math.log(p); A = 0.0; B = 0.0; m = 1
    while 1.0 / (m * m * p ** m) > 1e-14:
        f = 1 + np.mean(np.cos(2 * m * (M_tg + w / 2) * lp))
        t = 2 / math.pi ** 2 / (m * m * p ** m) * math.sin(m * w * lp / 2) ** 2 * f
        if m % 2: A += t
        else: B += t
        m += 1
    return A, B
ATAB = np.zeros((len(W), len(PR), 2))
for k, w in enumerate(W):
    M = int(math.floor((T2 - w - T1) / H + 1e-9)) + 1; tg = T1 + H * np.arange(M)
    for a, p in enumerate(PR): ATAB[k, a] = AT(p, w, tg)
REST = pair_cov_full(X[:, 4:], W, P[4:])       # primes > 7, T -> infinity form
psi = {p: (X[:, a].astype(float)[I] * X[:, a].astype(float)[J]) for a, p in enumerate(PR)}
bb = {p: (psi[p] ** 2) for p in PR}
def design(k, sel):
    Z = np.column_stack([np.ones(sel.sum())] + [psi[p][sel] for p in PR] + [bb[p][sel] for p in PR] + [REST[k][I[sel], J[sel]]])
    pred = sum(ATAB[k, a, 0] * psi[p][sel] + ATAB[k, a, 1] * bb[p][sel] for a, p in enumerate(PR)) + REST[k][I[sel], J[sel]]
    return Z, pred
def readings_and_jack(sel, groups):
    """sel: boolean mask over pairs; groups: list of index arrays of L-functions to delete one at a time (jackknife)
    returns dict of jackknife se for the amplitudes and the residual sd per w (no estimates are returned)."""
    Ii, Jj = I[sel], J[sel]; nf = N
    stats = {}
    COEF_full = []; COEF_jk = []; FULL_full = []; FULL_jk = []; ressd = []
    allf = np.unique(np.concatenate([Ii, Jj]))
    for k in range(len(W)):
        Z, pred = design(k, sel); y = C[k][Ii, Jj]
        ZtZ = Z.T @ Z; Zty = Z.T @ y; c = np.linalg.solve(ZtZ, Zty); e = y - Z @ c
        ressd.append(e.std())
        # full-law slope pieces: centred x = pred, y
        xs = pred; Sx = xs.sum(); Sy = y.sum(); Sxx = (xs * xs).sum(); Sxy = (xs * y).sum(); n = len(y)
        cj = []; fj = []
        for f in allf:
            m = (Ii == f) | (Jj == f); Zm = Z[m]; ym = y[m]; xm = xs[m]
            cj.append(np.linalg.solve(ZtZ - Zm.T @ Zm, Zty - Zm.T @ ym))
            fj.append(np.array([n - m.sum(), Sx - xm.sum(), Sy - ym.sum(), Sxx - (xm * xm).sum(), Sxy - (xm * ym).sum()]))
        COEF_full.append(c); COEF_jk.append(np.array(cj)); FULL_full.append(np.array([n, Sx, Sy, Sxx, Sxy])); FULL_jk.append(np.array(fj))
    COEF_full = np.array(COEF_full); COEF_jk = np.array(COEF_jk); FULL_full = np.array(FULL_full); FULL_jk = np.array(FULL_jk)
    nj = len(allf)
    def jse(full_val, jk_vals):
        return float(np.sqrt((nj - 1) / nj * ((jk_vals - jk_vals.mean()) ** 2).sum()))
    out = {}
    for a, p in enumerate(PR):
        A = ATAB[:, a, 0]; Bt = ATAB[:, a, 1]
        g = lambda co: (co[..., 1 + a] * A[:, None] if co.ndim == 3 else co[:, 1 + a] * A).sum(0) / (A ** 2).sum()
        out['gamma_%d' % p] = jse(None, (COEF_jk[:, :, 1 + a] * A[:, None]).sum(0) / (A ** 2).sum())
        out['gammaB_%d' % p] = jse(None, (COEF_jk[:, :, 5 + a] * Bt[:, None]).sum(0) / (Bt ** 2).sum())
        out['se_s_%d(w)' % p] = [jse(None, COEF_jk[k, :, 1 + a]) for k in range(len(W))]
    out['gamma_rest'] = jse(None, COEF_jk[:, :, 9].mean(0))
    # pooled full-law slope over pairs and w with an intercept per w
    def slope(F):
        n, Sx, Sy, Sxx, Sxy = [F[..., i] for i in range(5)]
        cov = Sxy - Sx * Sy / n; var = Sxx - Sx * Sx / n
        return cov.sum(0) / var.sum(0)
    out['gamma_full'] = jse(None, slope(FULL_jk))
    out['ressd'] = ressd; out['npairs'] = int(sel.sum()); out['nfun'] = int(nj)
    out['pred_sd_full'] = [float(design(k, sel)[1].std()) for k in range(len(W))]
    return out
def show(name, o):
    print('\n[%s] %d pairs of %d L-functions' % (name, o['npairs'], o['nfun']))
    print('  residual sd of the pair covariances after the regression, per w: ' + ' '.join('%.4f' % v for v in o['ressd']))
    print('  sd over pairs of the full-law prediction, per w:                 ' + ' '.join('%.4f' % v for v in o['pred_sd_full']))
    for p in PR:
        print('  p=%d jackknife se of the per-w coefficient: ' % p + ' '.join('%.5f' % v for v in o['se_s_%d(w)' % p]))
    print('  jackknife se of the amplitudes: ' + ', '.join('%s %.4f' % (k, o[k]) for k in ['gamma_%d' % p for p in PR] + ['gammaB_%d' % p for p in PR] + ['gamma_rest', 'gamma_full']))
    print('  expected z against "independent" (1/se): ' + ', '.join('%s %.0f' % (k, 1 / o[k]) for k in ['gamma_%d' % p for p in PR] + ['gammaB_2', 'gammaB_3', 'gamma_rest', 'gamma_full']))
allsel = np.ones(len(I), bool)
t0 = time.time(); o = readings_and_jack(allsel, None); show('old x old (all 546)', o); print('  (%.0fs)' % (time.time() - t0), flush=True)
rng = np.random.default_rng(140)
for rep in range(2):
    new = np.zeros(N, bool); new[rng.choice(N, 186, replace=False)] = True
    s_nn = new[I] & new[J]; s_no = new[I] != new[J]
    show('subset %d: 186 x 186 (mimics new x new)' % rep, readings_and_jack(s_nn, None))
    show('subset %d: 186 x 360 (mimics new x old)' % rep, readings_and_jack(s_no, None))
json.dump({'note': 'noise only', 'W': W.tolist()}, open(os.path.join(HERE, 'noise_meta.json'), 'w'))
