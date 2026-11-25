"""Retrospective on the fifty-seventh review's count files (read-only): twin weighted sums ln p ln(p+2) and theta sums by
p mod 420 on the grid X_k = floor(exp(k/2048)).  Folded to the 15 twin classes mod 210: spread of the weighted twin sums,
the linear part kappa (theta_b + theta_{b+2}) with kappa = S phi/A (round 155, Lambda weights), and the remainder; each
normalised by mean * w_eff, w_eff = sum w^2 / sum w (Poisson variance of a weighted sum), against the count-based
Hardy-Littlewood predictions of fourpoint.py.  -> rc_retro.json"""
import sys; sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/twinclass')
import json, math
import numpy as np
from scipy.stats import chi2
import fourpoint as F
RV = '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/review155/'
Q = 210; adm = F.admissible(Q, 2); K = len(adm); C = F.constants(Q, 2); kap = C['S'] * F.phi(Q) / C['A']
out = {}
for fname, lo, hi in (('rc_disc.npz', 1e8, 1e9), ('rc_disc.npz', 1e9, 2e9), ('rc_full.npz', 2e9, 1e10), ('rc_full.npz', 1e10, 2e10), ('rc_full.npz', 2e9, 2e10)):
    d = np.load(RV + fname); X = d['X'].astype(float)
    i0 = int(np.searchsorted(X, lo)); i1 = int(np.searchsorted(X, hi))      # cells (X[i-1], X[i]] for i in (i0, i1]
    sl = slice(i0 + 1, i1 + 1)
    L0, L1 = X[i0], X[i1]
    tw = d['twc'][sl].sum(0); th = d['thc'][sl].sum(0)
    fold = lambda a: np.bincount(np.arange(420) % Q, weights=a, minlength=Q)
    T = fold(tw)[adm]; TH = fold(th)
    Z = kap * (TH[adm] + TH[[(b + 2) % Q for b in adm]]); Z = Z - Z.mean() + T.mean()
    R = T - Z
    N2 = d['N2'][sl]; P2 = d['P2'][sl]
    weff = float(np.sum(P2 ** 2 / np.maximum(N2, 1)) / P2.sum())
    m = T.mean() * weff
    u = {'X': float(np.var(T, ddof=1) / m), 'Z': float(np.var(Z, ddof=1) / m), 'R': float(np.var(R, ddof=1) / m),
         'beta': float(np.cov(T, Z)[0, 1] / np.var(Z, ddof=1)), 'twins': int(N2.sum())}
    D = F.decompose(Q, 2, int(L0), int(L1 - L0))
    u.update(pred_X=D['vX'], pred_Z=D['vZ'], pred_R=D['vR'], ratio_X=u['X'] / D['vX'], ratio_R=u['R'] / D['vR'],
             p_lower_R=float(chi2.cdf(14 * u['R'] / D['vR'], 14)), L0=float(L0), L1=float(L1))
    out['%s (%.0e, %.0e]' % (fname, lo, hi)] = u
    print('(%.3g, %.3g] twins %d: X %.3f/%.3f (%.2f)  Z %.3f/%.3f  R %.3f/%.3f (%.2f, p_low %.3f)  beta %.2f' % (
        L0, L1, u['twins'], u['X'], D['vX'], u['ratio_X'], u['Z'], D['vZ'], u['R'], D['vR'], u['ratio_R'], u['p_lower_R'], u['beta']), flush=True)
json.dump(out, open('rc_retro.json', 'w'), indent=1)
