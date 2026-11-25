#!/usr/bin/env python3
"""power.py -- expected sampling errors for the fresh moduli (scaled from the pilot) and the frozen registry file.

se_m(q, s) ~ se_m(pilot, s) * sqrt(n_pilot / n_q) * (St2(q, s)/St2(pilot, s))^(m/2),  n = (q-1)/2.
Reads pilot/q10007.npz (pilot statistic), pilot/pilot_pred.json, reg/predictions_mock.json; writes
reg/predictions_frozen.json (+ sha256 printed) and prints the power table.
"""
import hashlib, json, math, sys
import numpy as np
sys.path.insert(0, '.')
from mock_test import cumulants_se

z = np.load('pilot/q10007.npz')
qp = int(z['q']); npil = (qp - 1) / 2
pilot_pred = json.load(open('pilot/pilot_pred.json'))['cells']
pse, pS = {}, {}
for j, s in enumerate(z['sig']):
    m, k, se = cumulants_se(z['psums'][j, 0], npil)
    c = pilot_pred[f'{qp}_{float(s):.2f}']
    pse[round(float(s), 2)] = (se, c['St2'])
    st = z['diag'][j]
    pS[round(float(s), 2)] = {mm: (k[mm] - st[mm - 2]) - c.get(f'R{mm}', float('nan')) for mm in (2, 3, 4)}
reg = json.load(open('reg/predictions_mock.json'))
cells = reg['cells']
Lp = math.log(qp)
for key, c in cells.items():
    s = round(c['s'], 2); q = c['q']
    se_p, st2p = pse[s]
    fac = math.sqrt(npil / ((q - 1) / 2))
    for mm in (2, 3, 4):
        c[f'se{mm}_exp'] = se_p[mm] * fac * (c['St2'] / st2p) ** (mm / 2)
        U = c[f'U{mm}']
        Dp = pS[s][mm] - U                          # pilot deviation of the structured part from RMT
        c[f'Dpilot{mm}'] = Dp
        c[f'D1L{mm}'] = Dp * Lp / c['logq']         # calibrated 1/L law for the deviation (T3)
    c['ycut'] = 20.0 if q > 500000 else 50.0
out = {'family': 'all non-principal chi mod prime q', 'fhat': '(1-|u|/s)_+', 'pilot_q': qp,
       'cells': cells}
blob = json.dumps(cells, sort_keys=True).encode()
h = hashlib.sha256(blob).hexdigest()
out['sha256_cells'] = h
json.dump(out, open('reg/predictions_frozen.json', 'w'), indent=1, sort_keys=True)
print('sha256(cells) =', h)
# power table: pooled over the 4 moduli of each group
groups = {}
for c in cells.values():
    groups.setdefault((round(c['q'], -4) if c['q'] < 500000 else 1000000, round(c['s'], 2)), []).append(c)
print('group  s    | m=2: U2 se_pool z(U) | m=3: U3 se_pool z(U) | m=4: U4 se_pool z(U) | R4 (deficit) | Dpilot4')
for (gq, s), cl in sorted(groups.items()):
    line = f'{gq:<7d} {s:.1f} '
    for mm in (2, 3, 4):
        se = 1 / math.sqrt(sum(1 / c[f'se{mm}_exp'] ** 2 for c in cl))
        U = cl[0][f'U{mm}']
        line += f'| {U:+.5f} {se:.6f} {abs(U)/se:7.1f} '
    line += f"| {np.mean([c.get('R4', float('nan')) for c in cl]):+.5f} | {cl[0]['Dpilot4']:+.5f}"
    print(line)
