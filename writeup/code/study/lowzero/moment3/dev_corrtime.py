"""dev_corrtime.py: per-tuple noise estimate of the integral over [100, 1000] from the variance and the correlation time
of the integrand (the 'random process' estimate), against the actual scatter of measured/recipe across tuples."""
import os, json, itertools, math
import numpy as np
DEV = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
V = {D: np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r149', 'dev_values', 'absL_%d.npy' % D)) ** 2 for D in DEV}
H = 0.05
def stats(f):
    m = f.mean(); v = f.var(); x = f - m
    n = len(f); maxlag = int(20 / H)
    ac = np.array([np.dot(x[:n - k], x[k:]) / ((n - k) * v) for k in range(0, maxlag)])
    # integrated autocorrelation time with cutoff at the first zero crossing (in t units), and the full sum to lag 20
    z = np.argmax(ac < 0) if np.any(ac < 0) else maxlag
    tau1 = H * (ac[:z].sum() - 0.5)
    tau20 = H * (ac.sum() - 0.5)
    se_rel = math.sqrt(v * 2 * tau1 * 900) / (m * 900)
    return m, v / m ** 2, tau1, tau20, se_rel, f.max() / m, (x ** 4).mean() / v ** 2
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dev_out.json')))
for k, rows in [(1, [(D,) for D in DEV]), (2, list(itertools.combinations(DEV, 2))), (3, list(itertools.combinations(DEV, 3)))]:
    S = []
    for t in rows:
        f = np.ones(18001)
        for D in t: f = f * V[D]
        S.append(stats(f))
    S = np.array(S)
    print('k = %d (%d tuples): rel var %.2f (median), max/mean %.0f, kurtosis %.0f, tau(first zero) %.3f, tau(sum to lag 20) %+.3f; '
          'naive rel se of the integral %.3f (median) [%.3f, %.3f]' % (k, len(rows), np.median(S[:, 1]), np.median(S[:, 5]), np.median(S[:, 6]),
          np.median(S[:, 2]), np.median(S[:, 3]), np.median(S[:, 4]), S[:, 4].min(), S[:, 4].max()))
print('actual sd of measured/recipe across tuples: singles %.3f, pairs %.3f, triples %.3f' % (
    np.std([r['meas'] / r['pred'] for r in d['single']], ddof=1), np.std([r['meas'] / r['pred'] for r in d['pair']], ddof=1),
    np.std([r['meas'] / r['pred']['full'] for r in d['triple']], ddof=1)))
