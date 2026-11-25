import json, math, os
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); R = json.load(open(D + '/s1_rows.json'))
B = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/blind'
P = json.load(open(B + '/predictions_blind.json')); RB = json.load(open(B + '/results_blind.json'))
cx = lambda v: complex(v[0], v[1])
reg = [k for k in R if R[k]['status'] == 'COMPLETE']
def lines(keys):
    num = np.array([sum((cx(R[k]['Gl'][nm]) * cx(R[k]['lines_pred'][nm]).conjugate()).real for nm in R[k]['Gl']) for k in keys])
    den = np.array([sum(abs(cx(R[k]['lines_pred'][nm])) ** 2 for nm in R[k]['Gl']) for k in keys])
    var = np.array([sum(abs(cx(R[k]['lines_pred'][nm])) ** 2 * R[k]['sd_line'] ** 2 for nm in R[k]['Gl']) for k in keys])
    b = num.sum() / den.sum(); se = math.sqrt(var.sum()) / den.sum(); n = len(keys)
    loo = np.array([(num.sum() - num[i]) / (den.sum() - den[i]) for i in range(n)])
    return b, se, math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
for name, ks in (('ALL', reg), ('degree 1', [k for k in reg if R[k]['deg'] == 1]), ('degree 2', [k for k in reg if R[k]['deg'] == 2]),
                 ('degree 3', [k for k in reg if R[k]['deg'] == 3]), ('degree 4', [k for k in reg if R[k]['deg'] == 4]),
                 ('degree >= 2', [k for k in reg if R[k]['deg'] >= 2]), ('degree >= 3', [k for k in reg if R[k]['deg'] >= 3])):
    print('lines %-12s n=%3d  b = %.4f  model se %.4f  jackknife se %.4f' % ((name, len(ks)) + lines(ks)))
# products: residual scatter against registered sd
pr = RB['products']
for rd, meas in (('M1', 'M1'), ('M2', 'M2')):
    ks = [k for k, e in P['products'].items() if e['cls'] == 'T5']
    pred = np.array([P['products'][k]['readings'][rd]['pred'] for k in ks]); sd = np.array([P['products'][k]['readings'][rd]['sd'] for k in ks])
    m = np.array([pr[k][meas][0] for k in ks])
    b = np.sum(m * pred) / np.sum(pred ** 2); n = len(ks)
    loo = np.array([(np.sum(m * pred) - m[i] * pred[i]) / (np.sum(pred ** 2) - pred[i] ** 2) for i in range(n)])
    js = math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
    res = m - b * pred
    print('%s: slope %.3f, registered se %.3f, jackknife se %.3f, rms residual %.1f vs registered sd %.1f (median), rms measured %.1f' % (
        rd, b, math.sqrt(np.sum(pred ** 2 * sd ** 2)) / np.sum(pred ** 2), js, math.sqrt(np.mean(res ** 2)), np.median(sd), math.sqrt(np.mean(m ** 2))))
# M2 with an intercept: is it the mechanism or a constant offset?
ks = [k for k, e in P['products'].items() if e['cls'] == 'T5']
pred = np.array([P['products'][k]['readings']['M2']['pred'] for k in ks]); m = np.array([pr[k]['M2'][0] for k in ks])
X = np.stack([np.ones(len(ks)), pred], 1); bv, rss, *_ = np.linalg.lstsq(X, m, rcond=None)
s2 = rss[0] / (len(ks) - 2); cov = s2 * np.linalg.inv(X.T @ X)
print('M2 = a + b*mech: a = %.2f ± %.2f, b = %.3f ± %.3f (OLS, residual sd %.1f); mean M2 %.2f, mean mech %.2f' % (bv[0], math.sqrt(cov[0, 0]), bv[1], math.sqrt(cov[1, 1]), math.sqrt(s2), m.mean(), pred.mean()))
