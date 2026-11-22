"""Round 118 (preregister/round118.txt): Sigma2(1) of round 117 regressed on 2|D, 3|D, 5|D and chi_D(2). usage: r118.py"""
import os, json
import numpy as np, cypari2
pari = cypari2.Pari()
H = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(H, 'r117.json')))['rows']
for r in rows: r['s2'] = int(pari.kronecker(r['D'], 2))
y = np.array([r['sigma2_1'] for r in rows]); X = np.array([[1, r['div2'], r['div3'], r['div5'], r['s2']] for r in rows], float)
c = np.linalg.lstsq(X, y, rcond=None)[0]; res = y - X @ c; s2 = res @ res / (len(y) - X.shape[1]); se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
names = ['const', '2|D', '3|D', '5|D', 'chi(2)']
for n, a, b in zip(names, c, se): print('%-7s %+.4f +- %.4f' % (n, a, b))
odd = [r for r in rows if not r['div2']]
print('odd D: chi(2)=+1 n=%d, chi(2)=-1 n=%d' % (sum(r['s2'] == 1 for r in odd), sum(r['s2'] == -1 for r in odd)))
out = {'coef': dict(zip(names, map(float, c))), 'se': dict(zip(names, map(float, se))), 'reading1': bool(abs(c[4]) < 2 * se[4])}
print('reading 1 (|coef chi(2)| < 2 se):', out['reading1'])
json.dump(out, open(os.path.join(H, 'r118.json'), 'w'), indent=1)
