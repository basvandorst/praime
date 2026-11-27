# Decisions H2-H4 of round 168 (odd Maass forms with the central zero).  Usage: decide168.py <X>
import sys, json, os, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass')
X = sys.argv[1]; key = [f for f in json.load(open(os.path.join(D, 'lmfdb_key.json'))) if f['parity'] == 'odd' and f['r'] < 20]
rows = sorted((json.loads(l) for l in open(os.path.join(D, 'scan_eps1_X%s_k2_J61_center.jsonl' % X))), key=lambda d: d['r'])
r = np.array([d['r'] for d in rows]); t = np.array([d['t_lp'] for d in rows]); tr = np.array([f['r'] for f in key])
print('== odd with central zero, X = %s: %d grid points, r %.2f..%.2f' % (X, len(r), r[0], r[-1]))
feas = r[t >= -0.002]; dist = np.array([np.abs(tr - x).min() for x in feas]) if len(feas) else np.array([])
false = feas[dist > 0.03] if len(feas) else feas
comps = []
for x in feas:
    if comps and x - comps[-1][1] < 0.0101: comps[-1][1] = x
    else: comps.append([x, x])
print('feasible grid points: %d in %d components: %s' % (len(feas), len(comps), '; '.join('[%.2f, %.2f]' % tuple(c) for c in comps)))
print('H2 (no false island): %d feasible points farther than 0.03 from every true odd eigenvalue%s -> %s' % (len(false), (' (r from %.2f to %.2f)' % (false.min(), false.max())) if len(false) else '', 'PASS' if len(false) == 0 else 'FAIL'))
ok3 = True; ok4 = True
for f in key:
    if f['r'] > r[-1] or f['r'] < r[0]: continue
    near = np.abs(r - f['r']) <= 0.02 + 1e-9; wide = np.abs(r - f['r']) <= 0.2 + 1e-9
    m = t[near].max(); p3 = m >= -0.02 and m >= t[wide].max() - 1e-12; ok3 &= p3
    i = int(np.argmin(np.abs(r - f['r']))); got = dict(zip(rows[i]['primes'], rows[i]['mean_ap']))
    d2, d3 = got[2] - f['ap']['2'], got[3] - f['ap']['3']; p4 = abs(d2) <= 0.05 and abs(d3) <= 0.05; ok4 &= p4
    print('  r_j = %.4f: largest t* within 0.02 = %+.5f, within 0.2 = %+.5f (%s); at r = %.2f: a_2 %+.3f (true %+.3f), a_3 %+.3f (true %+.3f) (%s)' % (
        f['r'], m, t[wide].max(), 'H3 ok' if p3 else 'H3 no', r[i], got[2], f['ap']['2'], got[3], f['ap']['3'], 'H4 ok' if p4 else 'H4 no'))
print('H3 -> %s;  H4 -> %s' % ('PASS' if ok3 else 'FAIL', 'PASS' if ok4 else 'FAIL'))
