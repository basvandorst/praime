"""Round 138 (preregister/round138.txt): the third cumulant kappa3(l) of the zero count of Dirichlet L-functions against
the frozen predictions of triple/predictions_kappa3.json (full law; GUE; Q0 = triple-orbit term kept for bad primes;
twopt = bad primes act through the two-point function only; diagonly = bad primes act only through the triple-orbit
term). Completeness rule, readings T1-T6 and decisions D1-D4 as registered.
usage: r138.py TEST [ZEROFILE]   TEST in Fn Bp R120 B Dn C; default zero files: lzeros/r138_Fn.npz, r138_Bp.npz,
r117/r120 parts, r134_B.npz, r136_Dn.npz, r134_C.npz (ZEROFILE: a synthetic set, for dry runs)"""
import os, sys, json, math
import numpy as np
from scipy import stats
from scipy.special import loggamma
H = os.path.dirname(os.path.abspath(__file__)); ZD = os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros'); T = sys.argv[1]
PJ = json.load(open(os.path.join(H, 'triple', 'predictions_kappa3.json'))); P = PJ[T]; LS = PJ['l']; lo, hi = P['height']; Q = P['D']
GUE = np.array(PJ['gue_kappa3']); GUEF = np.array(PJ['gue_kf3']); PR = (2, 3, 5, 7, 11, 13)
FILES = {'Fn': ['r138_Fn.npz'], 'Bp': ['r138_Bp.npz'], 'R120': ['r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz'], 'B': ['r134_B.npz'], 'Dn': ['r136_Dn.npz'], 'C': ['r134_C.npz']}
Z = {}
for f in ([sys.argv[2]] if len(sys.argv) > 2 else [os.path.join(ZD, f) for f in FILES[T]]):
    z = np.load(f); Z.update({k: z[k] for k in z.files})
def theta(t, q, a): return np.imag(loggamma(0.25 + a / 2 + 0.5j * t)) + t / 2 * math.log(q / math.pi)
def complete(g, D):
    """Registered rule: |N - theta/pi| <= 1.5 on 200 heights of the band (and its ends) and largest unfolded gap <= 3.0."""
    q = abs(D); a = 0 if D > 0 else 1; t = np.linspace(lo, hi, 200); s = np.searchsorted(g, t) - (theta(t, q, a) - theta(0.0, q, a)) / math.pi
    b = g[(g >= lo) & (g <= hi)]; u = b / (2 * math.pi) * np.log(q * b / (2 * math.pi * math.e)); return float(np.abs(s).max()), float(np.diff(u).max())
def cumulants(g, q, edge=0.0):
    u = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e)); k3, kf = [], []
    for l in LS:
        x = np.arange(u[0] + edge, u[-1] - edge - l, 0.01); c = (np.searchsorted(u, x + l) - np.searchsorted(u, x)).astype(float); m = c.mean()
        k = float(((c - m) ** 3).mean()); k3.append(k); kf.append(k - 3 * c.var() + 2 * m)
    return k3, kf
keep, excl, K3, KF, K3e = [], [], [], [], []
for i, D in enumerate(Q):
    g = np.sort(Z['D%d' % D]); s, gap = complete(g, D)
    if s > 1.5 or gap > 3.0: excl.append((D, round(s, 2), round(gap, 2))); continue
    b = g[(g >= lo) & (g <= hi)]; a, f = cumulants(b, abs(D)); K3.append(a); KF.append(f); K3e.append(cumulants(b, abs(D), 3.0)[0]); keep.append(i)
K3 = np.array(K3); KF = np.array(KF); K3e = np.array(K3e); n = len(keep); Dk = [Q[i] for i in keep]
PRD = {k: np.array(P[k])[keep] for k in ('full', 'Q0', 'twopt', 'diagonly', 'kf3')}; KF1 = PRD['kf3'] - (PRD['full'] - PRD['twopt'])
X = np.array([[1] + [int(abs(D) % p == 0) for p in PR] for D in Dk], float)
def reg(ix):
    Xi = X[ix]; A = np.linalg.pinv(Xi.T @ Xi); Hm = A @ Xi.T; return A, Hm
def gam(fit, pf, palt, w): d = pf - palt; return float((w * (fit - palt) * d).sum() / (w * d * d).sum())
def readings(ix):
    A, Hm = reg(ix); S = K3[ix]; fit = (Hm @ S)[1:]; res = S - X[ix] @ (Hm @ S); Cres = res.T @ res / (len(ix) - X.shape[1])
    se = np.sqrt(np.diag(A))[1:, None] * np.sqrt(np.diag(Cres))[None, :]; w = 1 / se[:3] ** 2
    pc = {k: (Hm @ PRD[k][ix])[1:] for k in ('full', 'Q0', 'twopt', 'diagonly')}
    d = S - PRD['full'][ix]; m = d.mean(0); C = np.cov(d.T); Ci = np.linalg.pinv(C)
    a = PRD['full'][ix].mean(0) - GUE; b = S.mean(0) - GUE; gA = float(a @ Ci @ b / (a @ Ci @ a))
    g = {k: gam(fit[:3], pc['full'][:3], pc[k][:3], w) for k in ('Q0', 'twopt', 'diagonly')}
    wl = 1 / np.diag(Cres); x = PRD['full'][ix] - PRD['Q0'][ix]; y = S - PRD['Q0'][ix]; x = x - x.mean(0); y = y - y.mean(0); gQp = float((wl * x * y).sum() / (wl * x * x).sum())
    return dict(m=m, C=C, gA=gA, gQ=g['Q0'], g2=g['twopt'], goff=g['diagonly'], gQp=gQp, fit=fit, se=se, pc=pc, A=A, Cres=Cres)
r = readings(np.arange(n)); rng = np.random.default_rng(1371); B = [readings(rng.integers(0, n, n)) for _ in range(2000)]
bse = {k: float(np.std([b[k] for b in B])) for k in ('gA', 'gQ', 'g2', 'goff', 'gQp')}; mse = np.std([b['m'] for b in B], axis=0)
T2h = float(n * r['m'] @ np.linalg.pinv(r['C']) @ r['m']); Fst = T2h * (n - 16) / ((n - 1) * 16); pT = float(stats.f.sf(Fst, 16, n - 16))
crit = float(stats.chi2.ppf(0.999, 16)); chi = {}
for j, p in enumerate((2, 3, 5)):
    Ci = np.linalg.pinv(r['A'][j + 1, j + 1] * r['Cres'])
    chi[p] = {k: float((r['fit'][j] - alt) @ Ci @ (r['fit'][j] - alt)) for k, alt in (('full', r['pc']['full'][j]), ('Q0', r['pc']['Q0'][j]), ('twopt', r['pc']['twopt'][j]), ('diagonly', r['pc']['diagonly'][j]), ('GUE', 0 * r['fit'][j]))}
i6 = [LS.index(l) for l in (3, 3.5, 4, 4.5, 5.5)]
mf = (KF - PRD['kf3']).mean(0); sef = (KF - PRD['kf3']).std(0, ddof=1) / math.sqrt(n); A0, H0 = reg(np.arange(n)); cff = (H0 @ KF)[1:]; rf = KF - X @ (H0 @ KF)
sef_c = np.sqrt(np.diag(A0))[1:, None] * np.sqrt((rf ** 2).sum(0) / (n - X.shape[1]))[None, :]; cfull = (H0 @ PRD['kf3'])[1:]; c2pt = (H0 @ KF1)[1:]
Dec = {'T1_pass': bool(pT >= 0.001 and np.all(np.abs(r['m']) <= 4 * mse)), 'T1_max_abs_m': float(np.abs(r['m']).max()),
       'D1_full_passes': bool(abs(r['gQ'] - 1) <= 4 * bse['gQ'] and chi[2]['full'] < crit), 'D1_Q0_rejected': bool(r['gQ'] >= 0.5 and r['gQ'] >= 5 * bse['gQ']),
       'D1b_full_passes': bool(abs(r['goff'] - 1) <= 4 * bse['goff']), 'D1b_diag_rejected': bool(r['goff'] >= 5 * bse['goff']),
       'D2_twopt_rejected': bool(r['g2'] >= 5 * bse['g2']), 'D3_GUE_rejected': bool(r['gA'] >= 5 * bse['gA'] and chi[2]['GUE'] > crit),
       'T6_signs': bool(np.all(r['fit'][0, i6] > 0))}
print('test %s: %d L-functions kept of %d, heights %d-%d; excluded by the completeness rule: %s' % (T, n, len(Q), lo, hi, excl if excl else 'none'))
print('groups (2,3,5,7,11,13): %s' % [int(X[:, j + 1].sum()) for j in range(6)])
print('T1 level: mean(kappa3 - full) per l  ' + ' '.join('%+.5f' % v for v in r['m']))
print('          bootstrap se               ' + ' '.join('%.5f' % v for v in mse))
print('          mean kappa3                ' + ' '.join('%+.4f' % v for v in K3.mean(0)))
print('          full law / GUE             ' + ' '.join('%+.4f' % v for v in PRD['full'].mean(0)) + ' / ' + ' '.join('%+.4f' % v for v in GUE))
print('          Hotelling T2 %.1f, p %.4g; max |m| %.5f, max |m|/se %.1f' % (T2h, pT, np.abs(r['m']).max(), np.abs(r['m'] / mse).max()))
print('T2 gamma_A %.4f +- %.4f' % (r['gA'], bse['gA']))
for p in (2, 3, 5): print('T3 chi2_%d on 16 lengths (0.1%% point %.1f): ' % (p, crit) + '  '.join('%s %.1f' % (k, v) for k, v in chi[p].items()))
for j, p in enumerate((2, 3, 5)):
    print('   p=%d fit  %s' % (p, ' '.join('%+.4f' % v for v in r['fit'][j]))); print('       full %s' % ' '.join('%+.4f' % v for v in r['pc']['full'][j])); print('       se   %s' % ' '.join('%.4f' % v for v in r['se'][j]))
print('T4 gamma_Q %.4f +- %.4f (pooled per L-function %.4f +- %.4f); gamma_off %.4f +- %.4f; gamma_2pt %.5f +- %.5f' % (r['gQ'], bse['gQ'], r['gQp'], bse['gQp'], r['goff'], bse['goff'], r['g2'], bse['g2']))
print('T5 k[3] level: mean(k[3] - law) ' + ' '.join('%+.4f' % v for v in mf) + ' | se ' + ' '.join('%.4f' % v for v in sef))
print('   k[3] - GUE: data ' + ' '.join('%+.3f' % v for v in KF.mean(0) - GUEF) + ' | law ' + ' '.join('%+.3f' % v for v in PRD['kf3'].mean(0) - GUEF))
for j, p in enumerate((2, 3, 5)): print('   k[3] coefficient p=%d: max |z| against the full law %.1f, against two-point only %.1f' % (p, np.abs((cff[j] - cfull[j]) / sef_c[j]).max(), np.abs((cff[j] - c2pt[j]) / sef_c[j]).max()))
print('T6 p = 2 coefficient at l = 3, 3.5, 4, 4.5, 5.5: ' + ' '.join('%+.4f' % v for v in r['fit'][0, i6]))
print('edge check (windows 3 spacings inside the band): mean(kappa3 - full) ' + ' '.join('%+.5f' % v for v in (K3e - PRD['full']).mean(0)))
print('decisions', Dec)
json.dump({'test': T, 'n': n, 'excluded': excl, 'm': r['m'].tolist(), 'm_se': mse.tolist(), 'T2': T2h, 'T2_p': pT, 'gamma_A': r['gA'], 'gamma_Q': r['gQ'], 'gamma_Q_pooled': r['gQp'], 'gamma_off': r['goff'],
           'gamma_2pt': r['g2'], 'se': bse, 'chi2': {str(k): v for k, v in chi.items()}, 'fit': r['fit'].tolist(), 'fit_se': r['se'].tolist(), 'pred_full': r['pc']['full'].tolist(), 'pred_Q0': r['pc']['Q0'].tolist(),
           'pred_diagonly': r['pc']['diagonly'].tolist(), 'kappa3_mean': K3.mean(0).tolist(), 'kf3_minus_law': mf.tolist(), 'kf3_se': sef.tolist(), 'edge': (K3e - PRD['full']).mean(0).tolist(), 'decisions': Dec,
           'kappa3': K3.tolist()}, open(os.path.join(H, 'r138_%s.json' % T) if len(sys.argv) <= 2 else os.devnull, 'w'))
