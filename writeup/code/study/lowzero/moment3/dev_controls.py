import sys, os, json, math, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import tmlib as T
DEV = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
t = 100.0 + 0.05 * np.arange(18001)
def trap(y): return 0.05 * (y.sum() - 0.5 * (y[0] + y[-1]))
V = {D: np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r149', 'dev_values', 'absL_%d.npy' % D)) ** 2 for D in DEV}
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dev_out.json')))
y4 = {}
print('fourth moment |L|^4 (Topacogullari Thm 1.1 polynomial): ratio measured / theorem; leading coef check')
for D in DEV:
    c = T.P4_coeffs(D)
    q = abs(D); fac = T.prime_factors(q); phi = q * np.prod([1 - 1 / p for p in fac])
    lead = 1 / (2 * math.pi ** 2) * (phi / q) ** 2 * np.prod([1 - 2 / (p + 1) for p in fac])
    f = lambda x: sum(c[k] * np.log(x) ** k for k in range(5))
    m = trap(V[D] ** 2); p = T.integrate(f, 100, 1000)
    y4[D] = math.log(m / p)
    # mid-range ratio of the 4th moment
    print('  %6d  ratio %.4f   c4 %.6e vs leading %.6e' % (D, m / p, c[4], lead))
ys = {r['D']: math.log(r['meas'] / r['pred']) for r in d['single']}
yp = {tuple(r['Ds']): math.log(r['meas'] / r['pred']) for r in d['pair']}
rows = d['triple']
y = np.array([math.log(r['meas'] / r['pred']['full']) for r in rows])
Ds = [tuple(r['Ds']) for r in rows]
z_pair = np.array([yp[(a, b)] + yp[(a, c)] + yp[(b, c)] for a, b, c in Ds])
z_four = np.array([y4[a] + y4[b] + y4[c] for a, b, c in Ds])
z_sing = np.array([ys[a] + ys[b] + ys[c] for a, b, c in Ds])
for name, z in [('sum of pair log-ratios', z_pair), ('sum of 4th-moment log-ratios', z_four), ('sum of single log-ratios', z_sing)]:
    cc = np.corrcoef(y, z)[0, 1]; b = np.polyfit(z, y, 1)[0]
    print('corr(triple log ratio, %s) = %.3f  slope %.2f  residual sd %.3f (raw %.3f)' % (name, cc, b, np.std(y - b * z), np.std(y)))
# character effects
X = np.array([[1.0] + [1.0 if c in D else 0.0 for c in DEV] for D in Ds])
b = np.linalg.lstsq(X, y, rcond=None)[0]
Xp = np.array([[1.0] + [1.0 if c in D else 0.0 for c in DEV] for D in yp])
bp = np.linalg.lstsq(Xp, np.array(list(yp.values())), rcond=None)[0]
print('\nchar   triple-effect  pair-effect  y4     y1     chi(2) chi(3) parity')
for i, c in enumerate(DEV):
    print('%6d  %+.3f        %+.3f       %+.3f  %+.3f  %+d     %+d     %s' % (c, b[i + 1], bp[i + 1], y4[c], ys[c], int(T.pari.kronecker(c, 2)), int(T.pari.kronecker(c, 3)), 'even' if c > 0 else 'odd'))
