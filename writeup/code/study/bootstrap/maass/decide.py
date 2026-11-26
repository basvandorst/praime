# Decisions D1-D5 of round 165 against the answer key.  Usage: decide.py <X> [tau]
import sys, json, os, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass')
X = sys.argv[1]; tau = float(sys.argv[2]) if len(sys.argv) > 2 else 0.002
key = json.load(open(os.path.join(D, 'lmfdb_key.json')))
for eps, par in [(1, 'odd'), (0, 'even')]:
    rows = sorted((json.loads(l) for l in open(os.path.join(D, 'scan_eps%d_X%s_k4_J61.jsonl' % (eps, X)))), key=lambda d: d['r'])
    r = np.array([d['r'] for d in rows]); t = np.array([d['t_lp'] for d in rows]); step = float(np.min(np.diff(r)))
    true = [f for f in key if f['parity'] == par and f['r'] <= 15]
    alltrue = [f['r'] for f in key if f['parity'] == par]
    print('== %s, X = %s, %d grid points (step %.2f), true eigenvalues <= 15: %s' % (par, X, len(r), step, ', '.join('%.4f' % f['r'] for f in true)))
    r1 = true[0]['r']; below = t[r < r1 - 0.3]
    print('D1 (nothing feasible below r1 - 0.3): max t* there = %.4f -> %s' % (below.max(), 'PASS' if below.max() < -tau else 'FAIL'))
    d2 = []
    for f in true:
        m = np.abs(r - f['r']) <= 0.05 + 1e-9; d2.append((f['r'], float(t[m].max()) if m.any() else np.nan))
    print('D2 (true eigenvalue feasible on the grid): ' + '; '.join('%.4f: best t* within 0.05 = %.4f' % x for x in d2) + ' -> %s' % ('PASS' if all(x[1] >= -tau for x in d2) else 'FAIL'))
    feas = t >= -tau; comps = []; i = 0
    while i < len(r):
        if feas[i]:
            j = i
            while j + 1 < len(r) and feas[j + 1]: j += 1
            comps.append((r[i], r[j])); i = j + 1
        else: i += 1
    comps15 = [c for c in comps if c[0] < 15]
    inside = [[x for x in alltrue if c[0] - step <= x <= c[1] + step] for c in comps15]
    ok3 = all(len(x) == 1 for x in inside) and all(any(c[0] - step <= f['r'] <= c[1] + step for c in comps15) for f in true)
    print('D3 (islands): components below 15: ' + ('; '.join('[%.2f, %.2f] holds %s' % (c[0], c[1], ['%.3f' % x for x in ins]) for c, ins in zip(comps15, inside)) or 'none') + ' -> %s' % ('PASS' if ok3 else 'FAIL'))
    mx = [r[i] for i in range(1, len(r) - 1) if t[i] > t[i - 1] and t[i] > t[i + 1]]
    d4 = [(f['r'], min((abs(m - f['r']) for m in mx), default=np.nan)) for f in true]
    print('D4 (strict local maximum within 0.1): ' + '; '.join('%.4f: nearest maximum at distance %.3f' % x for x in d4) + ' -> %s' % ('PASS' if all(x[1] <= 0.1 for x in d4) else 'FAIL'))
    ok5 = True
    for f in true:
        i = int(np.argmin(np.abs(r - f['r']))); row = rows[i]
        got = dict(zip(row['primes'], row['mean_ap'])); devs = [got[p] - f['ap'][str(p)] for p in (2, 3, 5, 7)]
        ok5 &= all(abs(d) <= 0.05 for d in devs)
        print('D5 at r = %.2f (true %.4f): found a_p %s, true %s, differences %s' % (r[i], f['r'], ' '.join('%+.3f' % got[p] for p in (2, 3, 5, 7)), ' '.join('%+.3f' % f['ap'][str(p)] for p in (2, 3, 5, 7)), ' '.join('%+.3f' % d for d in devs)))
    print('D5 -> %s' % ('PASS' if ok5 else 'FAIL'))
