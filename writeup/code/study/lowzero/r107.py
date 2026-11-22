"""Round 107 (preregister/round107.txt): lowest zero of every non-trivial character mod 1009 from S = 2 Re L'/L(1) +
ln(q/pi) + psi((1+a)/2), against a regression on Re chi(p), p <= 199. usage: r107.py"""
import json, math, os, cmath
import numpy as np, cypari2
from scipy.stats import spearmanr
from sympy import primerange, primitive_root
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
q = 1009; g0 = int(pari('znstar(1009,1).gen[1]'))
H = os.path.dirname(os.path.abspath(__file__))
ind = {}; v = 1
for e in range(q - 1): ind[v] = e; v = v * g0 % q
rows = []
for k in range(1, q - 1):
    L = pari.lfuncreate(pari('[znstar(%d,1),[%d]]' % (q, k)))
    T = 4; z = []
    while not z: z = [float(t) for t in pari.lfunzeros(L, pari('[%d,%d]' % (-T, T)))]; T *= 2
    gmin = min(abs(t) for t in z)
    L1 = complex(pari.lfun(L, 1)); dL1 = complex(pari.lfun(L, 1, 1))
    a = k % 2
    S = 2 * (dL1 / L1).real + math.log(q / math.pi) + float(pari.psi((1 + a) / 2.0))
    rows.append({'k': k, 'a': a, 'gmin': gmin, 'S': S})
S = np.array([r['S'] for r in rows]); gm = np.array([r['gmin'] for r in rows])
Smed = np.median(S); gmed = np.median(gm)
lhs = S - Smed + 1 / (0.25 + gmed ** 2)
pred = np.sqrt(np.clip(1 / np.maximum(lhs, 1e-9) - 0.25, 0, None))
rho = spearmanr(pred, gm)[0]
A1 = np.vstack([pred, np.ones_like(pred)]).T; c1 = np.linalg.lstsq(A1, gm, rcond=None)[0]
R2_pred = 1 - np.sum((gm - A1 @ c1) ** 2) / np.sum((gm - gm.mean()) ** 2)
P = list(primerange(2, 200))
X = np.array([[math.cos(2 * math.pi * r['k'] * ind[p % q] / (q - 1)) for p in P] + [r['a'], 1.0] for r in rows])
c2 = np.linalg.lstsq(X, gm, rcond=None)[0]
R2_small = 1 - np.sum((gm - X @ c2) ** 2) / np.sum((gm - gm.mean()) ** 2)
out = {'q': q, 'n': len(rows), 'spearman_pred_gmin': rho, 'R2_pred_line': R2_pred, 'R2_smallprimes_le199': R2_small,
       'median_abs_err_lowest100': float(np.median(np.abs(pred - gm)[np.argsort(gm)[:100]])), 'rows': rows}
print({k: v for k, v in out.items() if k != 'rows'})
json.dump(out, open(os.path.join(H, 'r107.json'), 'w'), indent=0)
