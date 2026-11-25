import json, itertools, math, os, sys
import numpy as np
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dev_out.json')))
DEV = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
def jk(Ds_list, vals, stat):
    """jackknife over characters: delete every tuple containing character c."""
    full = stat(np.array(vals)); reps = []
    for c in DEV:
        keep = [v for D, v in zip(Ds_list, vals) if c not in D]
        reps.append(stat(np.array(keep)))
    reps = np.array(reps); n = len(DEV)
    return full, math.sqrt((n - 1) / n * np.sum((reps - reps.mean()) ** 2))
for kind in ['pair', 'triple']:
    rows = d[kind]
    Ds = [tuple(r['Ds']) for r in rows]
    if kind == 'pair':
        r = np.array([r['meas'] / r['pred'] for r in rows]); rw = np.array([r['meas_w'] / r['pred_w'] for r in rows])
    else:
        r = np.array([x['meas'] / x['pred']['full'] for x in rows]); rw = np.array([x['meas_w'] / x['pred_w']['full'] for x in rows])
    print('%s: n = %d  mean ratio %.4f  sd %.4f  sd(log) %.4f   tapered: mean %.4f sd %.4f' % (kind, len(r), r.mean(), r.std(ddof=1), np.log(r).std(ddof=1), rw.mean(), rw.std(ddof=1)))
    m, se = jk(Ds, r, np.mean); print('   jackknife (chars) mean ratio %.4f +- %.4f;  naive se %.4f' % (m, se, r.std(ddof=1) / math.sqrt(len(r))))
    m, se = jk(Ds, rw, np.mean); print('   tapered jackknife %.4f +- %.4f' % (m, se))
    # additive character effects on log ratio
    y = np.log(r); X = np.array([[1.0] + [1.0 if c in D else 0.0 for c in DEV] for D in Ds])
    b, res, rk, sv = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ b
    print('   log-ratio var %.4f, after per-character effects %.4f (dof %d); char effects sd %.3f' % (y.var(ddof=1), resid.var(ddof=X.shape[1] - 1) if False else (resid ** 2).sum() / (len(y) - np.linalg.matrix_rank(X)), len(y) - np.linalg.matrix_rank(X), np.std(b[1:])))
    W = np.array([x['win'] for x in rows])
    print('   window ratios (9 x 100): sd of log over windows (mean over tuples) %.3f; whole-range sd(log) %.3f; ratio %.2f (sqrt 9 = 3)' % (np.log(W).std(axis=1, ddof=1).mean(), np.log(r).std(ddof=1), np.log(W).std(axis=1, ddof=1).mean() / np.log(r).std(ddof=1)))
    # lag correlation of windows across tuples
    Wl = np.log(W) - np.log(W).mean(axis=1, keepdims=True)
    c1 = np.mean([np.corrcoef(w[:-1], w[1:])[0, 1] for w in Wl]); print('   adjacent-window corr %.3f' % c1)
rows = d['triple']; Ds = [tuple(r['Ds']) for r in rows]
print('\ncomparators on triples: mean of M/P_c (jackknife se), sd of log(M/P_c)')
for c in ['full', 'diag', 'kirk', 'no3', 'ind', 'lead']:
    r = np.array([x['meas'] / x['pred'][c] for x in rows])
    m, se = jk(Ds, r, np.mean)
    print('  %-5s mean %.4f +- %.4f   sd(log) %.3f   z(mean-1) %+.1f' % (c, m, se, np.log(r).std(ddof=1), (m - 1) / se))
# slope of log M on log P_full
lm = np.log([x['meas'] for x in rows]); lp = np.log([x['pred']['full'] for x in rows])
def slope(idx_rows):
    a = np.array(idx_rows)
    A = np.vstack([np.ones(len(a)), a[:, 1]]).T
    return np.linalg.lstsq(A, a[:, 0], rcond=None)[0][1]
s, se = jk(Ds, list(zip(lm, lp)), slope)
print('\nslope log M on log P_full: %.4f +- %.4f (range of log P: %.2f .. %.2f)' % (s, se, lp.min(), lp.max()))
# gamma amplitudes vs comparators
print('amplitude gamma_c = sum (y - x_c)(x_f - x_c) / sum (x_f - x_c)^2  (1 = full, 0 = comparator c):')
for c in ['diag', 'kirk', 'no3', 'ind', 'lead']:
    xc = np.log([x['pred'][c] for x in rows])
    def gam(arr):
        a = np.array(arr); y, xf, xc_ = a[:, 0], a[:, 1], a[:, 2]
        return np.sum((y - xc_) * (xf - xc_)) / np.sum((xf - xc_) ** 2)
    g, se = jk(Ds, list(zip(lm, lp, xc)), gam)
    print('  vs %-5s gamma %.3f +- %.3f  (sd of x_f - x_c: %.3f; mean %.3f)' % (c, g, se, np.std(lp - xc), np.mean(lp - xc)))
