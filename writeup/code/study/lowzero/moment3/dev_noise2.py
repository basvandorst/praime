import sys, os, json, math
import numpy as np
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dev_out.json')))
DEV = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
rows = d['triple']; Ds = [tuple(r['Ds']) for r in rows]
rp = {tuple(r['Ds']): r['meas'] / r['pred'] for r in d['pair']}
r3 = np.array([r['meas'] / r['pred']['full'] for r in rows])
c = np.array([rp[(a, b)] + rp[(a, cc)] + rp[(b, cc)] - 3 for a, b, cc in Ds])
X = np.array([[1.0] + [1.0 if ch in D else 0.0 for ch in DEV] for D in Ds])
def comps(v):
    b = np.linalg.lstsq(X, v, rcond=None)[0]; res = v - X @ b
    dof = len(v) - np.linalg.matrix_rank(X); se2 = (res ** 2).sum() / dof
    XtXi = np.linalg.pinv(X.T @ X); a = b[1:] - b[1:].mean()
    return se2, max(a.var(ddof=1) - se2 * np.mean(np.diag(XtXi)[1:]), 0.0)
def jk(stat):
    reps = np.array([stat([i for i, D in enumerate(Ds) if ch not in D]) for ch in DEV]); n = len(DEV)
    return math.sqrt((n - 1) / n * np.sum((reps - reps.mean()) ** 2))
bl = np.polyfit(c, r3 - 1, 1); print('linear CV: OLS slope of (r3 - 1) on c = %.3f (intercept %.4f); corr %.3f' % (bl[0], bl[1], np.corrcoef(c, r3)[0, 1]))
# slope through origin
b0 = np.sum(c * (r3 - 1)) / np.sum(c * c); print('   through origin %.3f' % b0)
for beta in [0.0, 1.0, 1.5, 2.0, bl[0]]:
    v = r3 - beta * c
    se2, sa2 = comps(v)
    print('beta %.2f: mean %.4f jk-se(dev) %.4f | sigma_e %.3f sigma_a %.3f | model se nc=56,N=2000: %.4f ; nc=200,N=4000: %.4f' % (
        beta, v.mean(), jk(lambda k: v[k].mean()), math.sqrt(se2), math.sqrt(sa2), math.sqrt(9 * sa2 / 56 + se2 / 2000), math.sqrt(9 * sa2 / 200 + se2 / 4000)))
# the mean of the pair ratios alone (control) and its variance components
rpv = np.array(list(rp.values())); Dp = list(rp.keys())
Xp = np.array([[1.0] + [1.0 if ch in D else 0.0 for ch in DEV] for D in Dp])
b = np.linalg.lstsq(Xp, rpv, rcond=None)[0]; res = rpv - Xp @ b; dof = len(rpv) - np.linalg.matrix_rank(Xp); se2 = (res ** 2).sum() / dof
a = b[1:] - b[1:].mean(); sa2 = max(a.var(ddof=1) - se2 * np.mean(np.diag(np.linalg.pinv(Xp.T @ Xp))[1:]), 0)
print('pairs: mean %.4f sigma_e %.3f sigma_a %.3f; model se nc=56 all 1540 pairs: %.4f' % (rpv.mean(), math.sqrt(se2), math.sqrt(sa2), math.sqrt(4 * sa2 / 56 + se2 / 1540)))
