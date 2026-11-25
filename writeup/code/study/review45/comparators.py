import numpy as np, json, sys
from scipy import stats
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from analysis import *
SCR = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
fl = json.load(open(SCR + 'fastlaw_out.json'))
Dsel = [D for D in Ds if D != 1153 and str(D) in fl]
a = arrays(Dsel); n = len(Dsel); print('n =', n)
def bandW(D, lo, hi):
    t = np.linspace(lo, hi, 2001); return np.trapezoid(np.log(abs(D)*t/(2*np.pi)), t)
def comb(D, key):
    r = np.array(fl[str(D)][key]); return r
def full_from(r): return r[0] + r[2] + 3*r[1] - 2*ls
# full band reconstructed from the two sub-bands
FB = []
for D in Dsel:
    wl, wh = bandW(D, 100, 350), bandW(D, 350, 600)
    FB.append((wl*comb(D, 'lo') + wh*comb(D, 'hi'))/(wl + wh))
FB = np.array(FB)   # n x 5 x 16: kfP, sigma2, kQ, kQ11, kQall
print('reconstructed full vs JSON full: max |diff| %.1e' % np.max(np.abs(np.array([full_from(r) for r in FB]) - a['full'])))
print('kQ mine vs JSON kQ: max |diff| %.1e' % np.max(np.abs(FB[:, 2] - a['kQ'])))
print('Q0 check: full + (kQall - kQ) vs JSON Q0: max |diff| %.1e' % np.max(np.abs(a['full'] + FB[:, 4] - FB[:, 2] - a['Q0'])))
full = a['full']; y = a['k3']; X = a['X']
alts = {
  'GUE': np.tile(gue, (n, 1)),
  'Q0 (orbit kept)': a['Q0'],
  '(i-a) whole Q-part sign flipped': full - 2*a['kQ'],
  '(i-b) bad-prime orbit change flipped': 2*a['Q0'] - full,
  '(ii) diagonal: m=n=1 only': full - FB[:, 2] + FB[:, 3],
  '(iii) law at L-0.3': np.array([full_from(comb(D, 'Lm')) for D in Dsel]),
  '(iii) law at L+0.3': np.array([full_from(comb(D, 'Lp')) for D in Dsel]),
  'H_diag': a['diagonly'], 'H_2pt': a['twopt'],
  'GUE + diagonal Q-term (random-phase)': gue[None, :] + a['kQ'],
}
d = y - full; C = np.cov(d.T)/n; Ci = np.linalg.inv(C)
B, A, Cres, res = ols(X, y); pf = ols(X, full)[0]
idx = [(r, j) for r in (1, 2, 3) for j in range(16)]
V = np.array([[A[r, r2]*Cres[j, j2] for (r2, j2) in idx] for (r, j) in idx]); Vi = np.linalg.inv(V)
obs = np.array([B[r, j] - pf[r, j] for (r, j) in idx])
print('%-40s | level: gamma (law 0, alt 1), z_alt | contrasts p=2,3,5: gamma, z_alt | T2 of data vs alt, chi2_2/3/5 vs alt' % 'alternative')
for nm, alt in alts.items():
    Dl = alt.mean(0) - full.mean(0)
    g = Dl @ Ci @ d.mean(0)/(Dl @ Ci @ Dl); s = 1/np.sqrt(Dl @ Ci @ Dl)
    pa = ols(X, alt)[0]; dc = np.array([pa[r, j] - pf[r, j] for (r, j) in idx])
    gc = dc @ Vi @ obs/(dc @ Vi @ dc); sc = 1/np.sqrt(dc @ Vi @ dc)
    da = y - alt; ma = da.mean(0); T2a = n*ma @ np.linalg.solve(np.cov(da.T), ma)
    ch = []
    for r in (1, 2, 3):
        Vp = A[r, r]*Cres; e = B[r] - pa[r]; ch.append(e @ np.linalg.solve(Vp, e))
    print('%-40s | %+.3f +- %.3f (alt at %.1f se) | %+.3f +- %.3f (alt at %.1f se) | T2 %.0f, chi2 %.0f / %.0f / %.0f' % (
        nm, g, s, abs(1 - g)/s, gc, sc, abs(1 - gc)/sc, T2a, *ch))
    if nm.startswith('(iii)') or nm.startswith('(ii)'):
        print('      level alt-law x1e4:', np.round(Dl*1e4, 1).tolist())
# (iv) smooth 3-parameter function per prime
print('(iv) per-prime coefficient vs l: law (0 parameters) against a quadratic in l (3 parameters, GLS)')
for r, p in ((1, 2), (2, 3), (3, 5)):
    Vp = A[r, r]*Cres; Vpi = np.linalg.inv(Vp)
    M = np.vstack([np.ones(16), ls, ls**2]).T
    c = np.linalg.solve(M.T @ Vpi @ M, M.T @ Vpi @ B[r]); e = B[r] - M @ c
    chi_q = e @ Vpi @ e; el = B[r] - pf[r]; chi_l = el @ Vpi @ el
    # same quadratic fitted to the law's coefficients: how much of the law's shape is non-smooth
    cl = np.linalg.solve(M.T @ Vpi @ M, M.T @ Vpi @ pf[r]); el2 = pf[r] - M @ cl
    print('  p=%d: chi2 law %.1f (16 dof) | quadratic %.1f (13 dof, p=%.2g) | law shape vs its own best quadratic: chi2 %.0f' % (
        p, chi_l, chi_q, 1 - stats.chi2.cdf(chi_q, 13), el2 @ Vpi @ el2))
