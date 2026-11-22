"""Round 59 addendum (after the fifteenth review): a fairer baseline than the exponential law, the exponential times
the Hardy-Littlewood pair factor S({0,2,d,d+2})/S({0,2}) (round 60's discrete leading-order law), averaged over the
window with the twin density, on the same bands. usage: r59b.py"""
import os, json, math
import numpy as np
import common as C

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
r = json.load(open(os.path.join(OUT, 'r59.json')))
P3 = C.small_primes(20_000_000).astype(float); P3 = P3[P3 > 2]
C2 = math.exp(float(np.sum(np.log1p(-1 / (P3 - 1) ** 2))))
d = np.arange(6, 20001, 6); S = np.ones(len(d))
for p in C.small_primes(20001):
    p = int(p)
    if p < 5:
        continue
    nu = np.where(d % p == 0, 2, np.where((d % p == 2) | (d % p == p - 2), 3, 4))
    S *= (1 - nu / p) / (1 - 2 / p) / (1 - 1 / p) ** 2
u = np.linspace(math.log(1e9), math.log(2e9), 201); wt = np.exp(u) / u ** 2
P = np.zeros(len(d))
for ui, wi in zip(u, wt):
    rho = 2 * C2 / ui ** 2; x = S * np.exp(-rho * d); P += wi * x / x.sum()
P /= wt.sum()
law = dict(zip(d.tolist(), P.tolist()))
N = r['gaps']; data = {int(k): v for k, v in r['data'].items()}
chi = 0.0; dof = 0; rows = []
for row in r['compare']['rows']:
    lo, hi = row['band']; c = row['count']; e = N * sum(law[x] for x in range(lo, hi + 1, 6))
    rows.append({'band': [lo, hi], 'count': c, 'pair_factor_law': e, 'z': (c - e) / math.sqrt(e)})
    if e >= 50:
        chi += (c - e) ** 2 / e; dof += 1
out = {'chi2_pair_factor_baseline': chi, 'dof': dof, 'rows': rows}
print('chi2 pair-factor baseline %.1f on %d bands' % (chi, dof))
for x in rows[:4] + rows[-3:]:
    print(x['band'], x['count'], round(x['pair_factor_law'], 1), round(x['z'], 1))
json.dump(out, open(os.path.join(OUT, 'r59b.json'), 'w'), indent=1)
