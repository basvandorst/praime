# Decisions M1-M3 of round 171 (sequential bootstrap).  Usage: decide171.py <eps> <X>
import sys, json, os, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass')
eps = int(sys.argv[1]); X = sys.argv[2]; par = 'odd' if eps else 'even'
key = [f for f in json.load(open(os.path.join(D, 'lmfdb_key.json'))) if f['parity'] == par and f['ap']]
rows = sorted((json.loads(l) for l in open(os.path.join(D, 'mixed_eps%d_X%s.jsonl' % (eps, X)))), key=lambda d: d['r'])
r = np.array([d['r'] for d in rows]); t = np.array([d['t_lp'] for d in rows]); tr = np.array([f['r'] for f in key if r[0] - 0.02 <= f['r'] <= r[-1]])
print('== %s, sequential, X = %s: %d grid points, r %.2f..%.2f; true eigenvalues there: %s' % (par, X, len(r), r[0], r[-1], ', '.join('%.4f' % x for x in tr)))
comps = []
for x in r[t >= -0.002]:
    if comps and x - comps[-1][1] < 0.0101: comps[-1][1] = x
    else: comps.append([x, x])
print('feasible grid points: %d in %d pieces: %s' % (int((t >= -0.002).sum()), len(comps), '; '.join('[%.2f, %.2f]' % tuple(c) for c in comps)))
m1 = all(np.abs(tr - c[0]).min() <= 0.02 + 1e-9 for c in comps) if len(comps) else True
far = [x for x in r[t >= -0.002] if not any(rj - 0.02 - 1e-9 <= x <= rj + 0.15 + 1e-9 for rj in tr)]
print('M1 (every piece begins within 0.02 of a true eigenvalue) -> %s' % ('PASS' if m1 else 'FAIL'))
print('M2 (no long tails): %d feasible points outside [r_j - 0.02, r_j + 0.15] -> %s' % (len(far), 'PASS' if not far else 'FAIL'))
m3 = True
for f in key:
    if not (r[0] <= f['r'] <= r[-1]): continue
    idx = [i for i in range(len(r)) if r[i] >= f['r'] - 0.005 and t[i] >= -0.002]
    if not idx: print('  r_j = %.4f: no feasible grid point at or after it' % f['r']); m3 = False; continue
    i = idx[0]; got = dict(zip(rows[i]['primes'], rows[i]['mean_ap'])); d2 = got[2] - f['ap']['2']; d3 = got[3] - f['ap']['3']; okk = abs(d2) <= 0.05 and abs(d3) <= 0.05; m3 &= okk
    print('  r_j = %.4f: first feasible grid point %.2f (t* = %+.5f, %d known forms): a_2 %+.3f (true %+.3f), a_3 %+.3f (true %+.3f) %s' % (f['r'], r[i], t[i], len(rows[i]['known']), got[2], f['ap']['2'], got[3], f['ap']['3'], 'ok' if okk else 'NO'))
print('M3 -> %s' % ('PASS' if m3 else 'FAIL'))
