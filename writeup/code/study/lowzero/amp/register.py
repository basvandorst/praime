"""Numbers for the draft registration: for each fresh sample (C: the rounds 117/120 discriminants at 600-1000;
Dn: 700 < |D| <= 1000 at 100-600; Bp: test B's discriminants at 600-1000) the readings gamma, between, within, pooled,
gamma by sign, by l range and for p = 2, under BK (all 1), AMP (0.932) and LP3 (off-diagonal change at L + 0.3), with
the expected standard errors from the noise model (residual rows of rounds 117/120 at 100-600, scaled by
sqrt(zeros in reference / zeros in the new band)), and the error rates of the decision rule.
Uses only alt_*.npz (predictions) and r120_100_600.npz (noise model). usage: register.py -> register_out.txt"""
import os, sys, json, math
import numpy as np
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'bk'))
from common import *
def nzeros(D, lo, hi):
    k = abs(D); N = lambda t: t / (2 * math.pi) * math.log(k * t / (2 * math.pi * math.e))
    return N(hi) - N(lo)
SAMPLES = {nm: ([int(d) for d in np.load(os.path.join(HERE, 'alt_%s.npz' % nm))['D']], lo, hi)
           for nm, lo, hi in (('C', 600, 1000), ('Dn', 100, 600), ('Bp', 600, 1000))}

def bw(S, F, Dg, X):
    H = X @ np.linalg.pinv(X.T @ X) @ X.T; x = F - Dg; y = S - Dg
    xb = H @ x; xw = x - xb; yb = H @ y; yw = y - yb
    xbc = xb - xb.mean(0); ybc = yb - yb.mean(0)
    return (xbc * ybc).sum() / (xbc ** 2).sum(), (xw * yw).sum() / (xw ** 2).sum()

def readings(S, F, Dg, X, D):
    D = np.asarray(D); out = {}
    out['gamma'] = gamma(S, F, Dg, X)
    out['between'], out['within'] = bw(S, F, Dg, X)
    out['pooled'] = pooled(S, F, Dg)
    out['gamma D>0'] = gamma(S, F, Dg, X, ix=np.nonzero(D > 0)[0]); out['gamma D<0'] = gamma(S, F, Dg, X, ix=np.nonzero(D < 0)[0])
    out['gamma l<=4'] = gamma(S, F, Dg, X, cols=[j for j, l in enumerate(LS) if l <= 4]); out['gamma l>4'] = gamma(S, F, Dg, X, cols=[j for j, l in enumerate(LS) if l > 4])
    out['gamma p=2'] = gamma(S, F, Dg, X, rows=slice(0, 1))
    return out

z = np.load(os.path.join(HERE, 'r120_100_600.npz')); Dref = [int(d) for d in z['D']]; Xr = design(Dref)
Ar = np.linalg.pinv(Xr.T @ Xr) @ Xr.T; resid = (z['S'] - z['F']) - Xr @ (Ar @ (z['S'] - z['F']))
nref = np.mean([nzeros(D, 100, 600) for D in Dref])
REG = {}
for nm in ('C', 'Dn', 'Bp'):
    Ds, lo, hi = SAMPLES[nm]; a = np.load(os.path.join(HERE, 'alt_%s.npz' % nm)); F, Dg, Fl = a['F'], a['Dg'], a['Flp3']; X = design(Ds); n = len(Ds)
    nz = np.mean([nzeros(D, lo, hi) for D in Ds]); sc = math.sqrt(nref / nz); rng = np.random.default_rng(77)
    sims = [readings(F + sc * resid[rng.integers(0, len(resid), n)], F, Dg, X, Ds) for _ in range(600)]
    keys = list(sims[0].keys()); se = {k: float(np.std([s[k] for s in sims])) for k in keys}
    amp = readings(Dg + 0.932 * (F - Dg), F, Dg, X, Ds); lp3 = readings(Fl, F, Dg, X, Ds)
    Lm = np.mean([float(np.dot(*L_nodes(abs(D), lo, hi, 8))) for D in Ds])
    print('sample %s: n %d, heights %d-%d, |D| %d-%d, mean L %.3f, mean zeros per set %.0f, groups 2..13 %s' % (
        nm, n, lo, hi, min(abs(d) for d in Ds), max(abs(d) for d in Ds), Lm, nz, [int(X[:, i + 1].sum()) for i in range(6)]))
    print('   %-11s %7s %7s %7s %7s   (BK | AMP | LP3 | expected se)' % ('reading', 'BK', 'AMP', 'LP3', 'se'))
    for k in keys:
        print('   %-11s %7.3f %7.3f %7.3f %7.3f' % (k, 1.0, amp[k], lp3[k], se[k]))
    cor = np.corrcoef([s['between'] for s in sims], [s['within'] for s in sims])[0, 1]
    print('   null correlation between/within %.2f; sd of (between - within) %.3f' % (cor, np.std([s['between'] - s['within'] for s in sims])))
    REG[nm] = dict(n=n, L=Lm, zeros=nz, se=se, amp=amp, lp3=lp3)
# decision rule on gamma: threshold halfway between BK and AMP
for combo in (('C',), ('C', 'Dn'), ('C', 'Dn', 'Bp')):
    w = np.array([1 / REG[c]['se']['gamma'] ** 2 for c in combo]); s = 1 / math.sqrt(w.sum())
    lp = float((w * np.array([REG[c]['lp3']['gamma'] for c in combo])).sum() / w.sum())
    thr = (1 + 0.932) / 2
    print('combined %-10s se %.4f; threshold %.3f: P(below | BK) %.3f, P(above | AMP) %.3f, P(above | LP3 %.3f) %.3f; BK-AMP separation %.1f sigma' % (
        '+'.join(combo), s, thr, norm.cdf((thr - 1) / s), norm.sf((thr - 0.932) / s), lp, norm.sf((thr - lp) / s), (1 - 0.932) / s))
json.dump(REG, open(os.path.join(HERE, 'register.json'), 'w'), indent=1, default=float)
