import sys, os, json, math, itertools
import numpy as np
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dev_out.json')))
DEV = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
rows = d['triple']; Ds = [tuple(r['Ds']) for r in rows]
rp = {tuple(r['Ds']): r['meas'] / r['pred'] for r in d['pair']}
r3 = np.array([r['meas'] / r['pred']['full'] for r in rows])
y = np.log(r3)
z = np.array([math.log(rp[(a, b)]) + math.log(rp[(a, c)]) + math.log(rp[(b, c)]) for a, b, c in Ds])
zr = np.array([rp[(a, b)] * rp[(a, c)] * rp[(b, c)] for a, b, c in Ds])
X = np.array([[1.0] + [1.0 if c in D else 0.0 for c in DEV] for D in Ds])
def comps(v, X):
    b = np.linalg.lstsq(X, v, rcond=None)[0]; res = v - X @ b
    dof = len(v) - np.linalg.matrix_rank(X); se2 = (res ** 2).sum() / dof
    XtXi = np.linalg.pinv(X.T @ X)
    a = b[1:]; a = a - a.mean()
    var_hat = a.var(ddof=1); noise = se2 * np.mean(np.diag(XtXi)[1:])
    return se2, max(var_hat - noise, 0.0)
def jk(vals, stat):
    reps = []
    for c in DEV:
        keep = [i for i, D in enumerate(Ds) if c not in D]
        reps.append(stat(keep))
    reps = np.array(reps); n = len(DEV)
    return math.sqrt((n - 1) / n * np.sum((reps - reps.mean()) ** 2))
for beta in [0.0, 1.0, 1.25, 1.5]:
    v = y - beta * z
    se2, sa2 = comps(v, X)
    for nc, N in [(56, 2000), (56, 27720), (200, 2000), (200, 4000)]:
        pass
    pred = {(nc, N): math.sqrt(9 * sa2 / nc + se2 / N) for nc, N in [(56, 2000), (60, 2000), (150, 3000), (200, 4000)]}
    print('beta %.2f: sigma_e %.3f sigma_a %.3f | jackknife se of mean on dev %.3f | model se: ' % (beta, math.sqrt(se2), math.sqrt(sa2), jk(v, lambda k: v[k].mean())) +
          ' '.join('nc=%d,N=%d: %.4f' % (k[0], k[1], s) for k, s in pred.items()))
# ratio version of control variate: r3 / (pair ratio product)^beta
for beta in [0.0, 1.0, 1.5]:
    v = r3 / zr ** beta
    print('ratio CV beta %.1f: mean %.4f jackknife se %.4f' % (beta, v.mean(), jk(v, lambda k: v[k].mean())))
# gamma no3 / kirk with CV
lm = np.log([x['meas'] for x in rows]); lp = np.log([x['pred']['full'] for x in rows])
for c in ['no3', 'kirk', 'diag', 'ind']:
    xc = np.log([x['pred'][c] for x in rows])
    for beta in [0.0, 1.5]:
        yy = lm - beta * z
        def gam(k, centered=True):
            k = np.array(k); u = (xf := lp[k]) - xc[k]; w = yy[k] - xc[k]
            if centered: u = u - u.mean(); w = w - w.mean()
            return np.sum(w * u) / np.sum(u * u)
        g = gam(list(range(len(Ds))))
        g0 = gam(list(range(len(Ds))), False)
        print('%-5s beta %.1f: gamma centered %.3f +- %.3f   uncentered %.3f +- %.3f' % (c, beta, g, jk(None, gam), g0, jk(None, lambda k: gam(k, False))))
