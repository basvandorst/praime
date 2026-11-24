"""RETROSPECTIVE table: gamma measured in test B and in the rounds 117/120 sample (same heights, same statistic), by
split, against what each alternative predicts on the same split:
   BK  = 1 (the Bogomolny-Keating law), AMP = 0.932 (plain amplitude), LP3 = the off-diagonal change at L + 0.3.
Then the joint reading: null covariance of (between, within) slopes per sample, the common amplitude by GLS from the
four numbers, and the chi2 of each hypothesis.  usage: table.py -> table_out.txt"""
import os, sys, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'bk'))
from common import *
import gridcache

G = gridcache.get()
rng = np.random.default_rng(2024)
samples = {}
DB, SB, FB, DgB = load_B(); a = np.load(os.path.join(HERE, 'alt_B.npz')); samples['B'] = dict(D=DB, S=SB, F=FB, Dg=DgB, Fl=a['Flp3'])
z = np.load(os.path.join(HERE, 'r120_100_600.npz')); a = np.load(os.path.join(HERE, 'alt_R120.npz'))
pos = {int(D): i for i, D in enumerate(a['D'])}; order = [pos[int(D)] for D in z['D']]
assert np.abs(a['F'][order] - z['F']).max() < 1e-12
samples['R120'] = dict(D=[int(D) for D in z['D']], S=z['S'], F=z['F'], Dg=z['Dg'], Fl=a['Flp3'][order])
# band-specific data and predictions for B (zeros of r134_B.npz) and R120 (saved by retro.py); LP3 per band
zB = np.load(os.path.join(LZ, 'r134_B.npz'))
bands = {}
for (lo, hi) in ((100, 300), (300, 600)):
    SBb = np.array([sigma2_set(zB['D%d' % D], abs(D), lo, hi) for D in DB]); FBb, DgBb = predict(G, DB, lo, hi); FlBb, _ = predict(G, DB, lo, hi, off_dL=0.3)
    bands[('B', lo, hi)] = (SBb, FBb, DgBb, FlBb)
    zz = np.load(os.path.join(HERE, 'r120_%d_%d.npz' % (lo, hi))); FlRb, _ = predict(G, list(zz['D']), lo, hi, off_dL=0.3)
    bands[('R120', lo, hi)] = (zz['S'], zz['F'], zz['Dg'], FlRb)

def row(label, S, F, Dg, Fl, X, ix=None, **kw):
    ix = np.arange(len(S)) if ix is None else np.asarray(ix)
    g = gamma(S[ix], F[ix], Dg[ix], X[ix], **kw)
    bs = [gamma(S[j], F[j], Dg[j], X[j], **kw) for j in (ix[rng.integers(0, len(ix), len(ix))] for _ in range(500))]
    lp3 = gamma(Fl[ix], F[ix], Dg[ix], X[ix], **kw)
    print('  %-22s %.3f +- %.3f   | BK 1.000  AMP 0.932  LP3 %.3f   | z(BK) %+.1f  z(AMP) %+.1f  z(LP3) %+.1f' % (
        label, g, np.std(bs), lp3, (g - 1) / np.std(bs), (g - 0.932) / np.std(bs), (g - lp3) / np.std(bs)))
    return g, float(np.std(bs)), lp3

for nm in ('B', 'R120'):
    s = samples[nm]; X = design(s['D']); D = np.array(s['D']); n = len(D)
    print('%s (n = %d), gamma measured +- bootstrap se | predicted | z of the measured value under each:' % (nm, n))
    row('all', s['S'], s['F'], s['Dg'], s['Fl'], X)
    row('D > 0', s['S'], s['F'], s['Dg'], s['Fl'], X, ix=np.nonzero(D > 0)[0])
    row('D < 0', s['S'], s['F'], s['Dg'], s['Fl'], X, ix=np.nonzero(D < 0)[0])
    for lo, hi in ((100, 300), (300, 600)):
        Sb, Fb, Dgb, Flb = bands[(nm, lo, hi)]; row('heights %d-%d' % (lo, hi), Sb, Fb, Dgb, Flb, X)
    row('l <= 4', s['S'], s['F'], s['Dg'], s['Fl'], X, cols=[j for j, l in enumerate(LS) if l <= 4])
    row('l > 4', s['S'], s['F'], s['Dg'], s['Fl'], X, cols=[j for j, l in enumerate(LS) if l > 4])
    row('p = 2', s['S'], s['F'], s['Dg'], s['Fl'], X, rows=slice(0, 1))
    row('p = 3', s['S'], s['F'], s['Dg'], s['Fl'], X, rows=slice(1, 2))
    ab = np.abs(D); med = np.median(ab)
    row('|D| <= %d' % med, s['S'], s['F'], s['Dg'], s['Fl'], X, ix=np.nonzero(ab <= med)[0])
    row('|D| > %d' % med, s['S'], s['F'], s['Dg'], s['Fl'], X, ix=np.nonzero(ab > med)[0])

# ---------------- joint reading ----------------
def bw(S, F, Dg, X):
    H = X @ np.linalg.pinv(X.T @ X) @ X.T; x = F - Dg; y = S - Dg
    xb = H @ x; xw = x - xb; yb = H @ y; yw = y - yb
    xbc = xb - xb.mean(0); ybc = yb - yb.mean(0)
    return np.array([(xbc * ybc).sum() / (xbc ** 2).sum(), (xw * yw).sum() / (xw ** 2).sum()])

est, covs, preds_lp3 = [], [], []
print('\nJoint reading (between = indicator part, within = the rest of the per-L-function off pattern):')
for nm in ('B', 'R120'):
    s = samples[nm]; X = design(s['D']); n = len(s['D']); A = np.linalg.pinv(X.T @ X) @ X.T
    res = (s['S'] - s['F']) - X @ (A @ (s['S'] - s['F']))
    sims = np.array([bw(s['F'] + res[rng.integers(0, n, n)], s['F'], s['Dg'], X) for _ in range(1000)])
    C = np.cov(sims.T); e = bw(s['S'], s['F'], s['Dg'], X); lp = bw(s['Fl'], s['F'], s['Dg'], X)
    print('  %-5s between %.3f, within %.3f; null sd %.3f, %.3f, null correlation %+.2f; null mean %.3f %.3f; between - within %+.3f (z %+.1f); LP3 predicts %.3f %.3f' % (
        nm, e[0], e[1], math.sqrt(C[0, 0]), math.sqrt(C[1, 1]), C[0, 1] / math.sqrt(C[0, 0] * C[1, 1]), sims[:, 0].mean(), sims[:, 1].mean(),
        e[0] - e[1], (e[0] - e[1]) / math.sqrt(C[0, 0] + C[1, 1] - 2 * C[0, 1]), lp[0], lp[1]))
    est.append(e); covs.append(C); preds_lp3.append(lp)
y = np.concatenate(est); V = np.zeros((4, 4)); V[:2, :2] = covs[0]; V[2:, 2:] = covs[1]; Vi = np.linalg.inv(V)
one = np.ones(4); ahat = float(one @ Vi @ y / (one @ Vi @ one)); se = math.sqrt(1 / (one @ Vi @ one))
chi = lambda m: float((y - m) @ Vi @ (y - m))
print('  common amplitude (GLS, 4 numbers): %.3f +- %.3f; chi2 about it %.1f on 3 dof' % (ahat, se, chi(ahat * one)))
print('  chi2 on 4 dof: BK (all 1) %.1f;  AMP (all 0.932) %.1f;  LP3 (its own predictions) %.1f' % (chi(one), chi(0.932 * one), chi(np.concatenate(preds_lp3))))
from scipy.stats import chi2 as C2
for nm, v, d in (('BK', chi(one), 4), ('AMP', chi(0.932 * one), 4), ('LP3', chi(np.concatenate(preds_lp3)), 4), ('common amplitude', chi(ahat * one), 3)):
    print('    %-17s p = %.2g' % (nm, C2.sf(v, d)))
