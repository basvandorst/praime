"""Final numbers for round_draft.txt:
 - per fresh sample (C, Dn, Bp): rho at l = 0.5 and the p = 2 off-diagonal cells at l = 6..8 under BK, AMP and LP3;
 - error rates of the decision rules on gamma and on the pooled slope;
 - the descriptive conductor trend of gamma (four subsamples of B and R120, retrospective) and what it implies for C, Dn, Bp;
 - writes predictions_alt.json (per-L-function predictions on the 16 primary lengths: BK full, diag, LP3 full) and its sha256.
usage: regfinal.py -> regfinal_out.txt"""
import os, sys, json, math, hashlib
import numpy as np
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'bk'))
from common import *

REG = json.load(open(os.path.join(HERE, 'register.json')))
out = {'l': LS}
for nm, lo, hi in (('C', 600, 1000), ('Dn', 100, 600), ('Bp', 600, 1000)):
    a = np.load(os.path.join(HERE, 'alt_%s.npz' % nm)); Ds = [int(d) for d in a['D']]; X = design(Ds); A = np.linalg.pinv(X.T @ X) @ X.T
    F, Dg, Fl = a['F'], a['Dg'], a['Flp3']
    cf = (A @ F)[1:]; cd = (A @ Dg)[1:]; cl = (A @ Fl)[1:]; ca = cd + 0.932 * (cf - cd)
    rho = lambda c: float((c[:, 0] * cd[:, 0]).sum() / (cd[:, 0] ** 2).sum())
    i2 = [LS.index(l) for l in (6, 6.5, 7, 7.5, 8)]
    print('%s (heights %d-%d): rho(l=0.5) BK %.3f AMP %.3f LP3 %.3f | p=2 off cells l=6..8: BK %s | AMP %s | LP3 %s' % (
        nm, lo, hi, rho(cf), rho(ca), rho(cl), ' '.join('%+.4f' % (cf - cd)[0, j] for j in i2), ' '.join('%+.4f' % (ca - cd)[0, j] for j in i2),
        ' '.join('%+.4f' % (cl - cd)[0, j] for j in i2)))
    out[nm] = dict(height=[lo, hi], D=Ds, full=F.tolist(), diag=Dg.tolist(), full_LP3=Fl.tolist())
thr = (1 + 0.932) / 2
for key in ('gamma', 'pooled'):
    for combo in (('C',), ('C', 'Dn'), ('C', 'Dn', 'Bp')):
        w = np.array([1 / REG[c]['se'][key] ** 2 for c in combo]); s = 1 / math.sqrt(w.sum())
        lp = float((w * np.array([REG[c]['lp3'][key] for c in combo])).sum() / w.sum())
        print('rule on %-6s %-9s se %.4f: P(below %.3f | BK) %.4f, P(above | AMP) %.4f, P(above | LP3 %.3f) %.4f' % (
            key, '+'.join(combo), s, thr, norm.cdf((thr - 1) / s), norm.sf((thr - 0.932) / s), lp, norm.sf((thr - lp) / s)))
# descriptive conductor trend (retrospective numbers from table_out.txt)
sub = [(175.0, 1.001, 0.030), (325.0, 1.017, 0.043), (475.0, 0.956, 0.041), (625.0, 0.917, 0.039)]
# mean |D| of each half, recomputed
P = json.load(open(os.path.join(LOW, 'bk', 'predictions_per_lfunction.json')))
DB = np.abs(P['B']['D']); DR = np.abs(P['C']['D'])
mB = np.median(DB); mR = np.median(DR)
xs = [np.mean(np.log(DR[DR <= mR])), np.mean(np.log(DR[DR > mR])), np.mean(np.log(DB[DB <= mB])), np.mean(np.log(DB[DB > mB]))]
ys = np.array([s[1] for s in sub]); ws = np.array([1 / s[2] ** 2 for s in sub]); xs = np.array(xs)
xb = (ws * xs).sum() / ws.sum(); yb = (ws * ys).sum() / ws.sum(); Sxx = (ws * (xs - xb) ** 2).sum(); b = (ws * (xs - xb) * (ys - yb)).sum() / Sxx
chi_trend = float((ws * (ys - yb - b * (xs - xb)) ** 2).sum())
print('descriptive trend of gamma in mean ln|D| (4 half-samples): slope %.3f +- %.3f per unit ln|D|; chi2 about the line %.2f on 2 dof' % (b, 1 / math.sqrt(Sxx), chi_trend))
for nm in ('C', 'Dn', 'Bp'):
    x = float(np.mean(np.log(np.abs(out[nm]['D'])))); pr = yb + b * (x - xb); e = math.sqrt(1 / ws.sum() + (x - xb) ** 2 / Sxx)
    print('   trend predicts for %-2s (mean ln|D| %.3f): gamma %.3f +- %.3f (trend error), expected statistical se %.3f' % (nm, x, pr, e, REG[nm]['se']['gamma']))
f = os.path.join(HERE, 'predictions_alt.json'); json.dump(out, open(f, 'w'))
print('predictions_alt.json sha256 %s' % hashlib.sha256(open(f, 'rb').read()).hexdigest())
