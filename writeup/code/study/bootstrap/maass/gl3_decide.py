# Decisions G1-G4 of round 167 (SL(3,Z)).  Usage: gl3_decide.py <X>
import sys, json, os, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass')
X = sys.argv[1]; key = sorted(json.load(open(os.path.join(D, 'gl3_key.json')))['forms'], key=lambda f: f['laplace'])
loc = [json.loads(l) for l in open(os.path.join(D, 'gl3_local_X%s.jsonl' % X))]
g1 = g2 = g4 = True
for n, f in enumerate(key[:5]):
    l1, l2 = f['lambda']; e = [r for r in loc if r['kind'] == 'exact' and abs(r['l1'] - l1) < 1e-6 and abs(r['l2'] - l2) < 1e-6][0]
    grid = [r for r in loc if r['kind'] == 'grid' and abs(r['l1'] - round(l1, 1)) < 0.41 and abs(r['l2'] - round(l2, 1)) < 0.41]
    inf = [np.hypot(r['l1'] - l1, r['l2'] - l2) for r in grid if r['t_lp'] < -0.002]; dmin = min(inf) if inf else np.inf
    a2 = complex(*e['a_p'][0]); t2 = complex(*f['A_p1']['2']); a3 = complex(*e['a_p'][1]); t3 = complex(*f['A_p1']['3'])
    g1 &= e['t_lp'] >= -0.002; g2 &= dmin <= 0.2
    if n < 3: g4 &= abs(a2 - t2) <= 0.1
    print('form %d (%.4f, %.4f), Laplace %.2f: t* at the exact point %+.5f; local grid %d points, %d infeasible, nearest infeasible at %.3f; a_2 found %+.3f%+.3fi, true %+.3f%+.3fi, |difference| %.3f; a_3 |difference| %.3f' % (
        n + 1, l1, l2, f['laplace'], e['t_lp'], len(grid), len(inf), dmin, a2.real, a2.imag, t2.real, t2.imag, abs(a2 - t2), abs(a3 - t3)))
print('G1 (truth feasible) -> %s;  G2 (infeasible point within 0.2 of each of the five) -> %s;  G4 (a_2 of the three lowest within 0.1) -> %s' % tuple('PASS' if x else 'FAIL' for x in (g1, g2, g4)))
cf = os.path.join(D, 'gl3_coarse_X%s.jsonl' % X)
if os.path.exists(cf):
    co = [json.loads(l) for l in open(cf)]; q = lambda r: r['l1']**2 + r['l1']*r['l2'] + r['l2']**2
    low = [r for r in co if q(r) <= 180]; feas = [r for r in co if r['t_lp'] >= -0.002]
    print('coarse grid: %d of 861 points done; %d points with l1^2 + l1 l2 + l2^2 <= 180, of which feasible: %d -> G3 %s' % (len(co), len(low), sum(r['t_lp'] >= -0.002 for r in low), 'PASS' if all(r['t_lp'] < -0.002 for r in low) else 'FAIL'))
    if feas:
        m = min(feas, key=q); print('feasible coarse points: %d; smallest l1^2 + l1 l2 + l2^2 among them: %.2f at (%.1f, %.1f) (lowest form: %.2f)' % (len(feas), q(m), m['l1'], m['l2'], key[0]['laplace'] - 1))
        forms = [f for f in key if f['laplace'] < 400]
        inf = np.array([(r['l1'], r['l2']) for r in co if r['t_lp'] < -0.002])
        d = [np.hypot(inf[:, 0] - f['lambda'][0], inf[:, 1] - f['lambda'][1]).min() for f in forms]
        print('forms with Laplace eigenvalue < 400: %d; distance to the nearest infeasible coarse point: %s' % (len(forms), ' '.join('%.2f' % x for x in d)))
        print('feasible set by row (l2 downwards), x = feasible, . = infeasible, blank = not done:')
        done = {(round(r['l1'], 1), round(r['l2'], 1)): r['t_lp'] for r in co}
        for j in range(40, -1, -2):
            print('  l2=%4.1f ' % (0.5*j) + ''.join((' ' if (0.5*i, 0.5*j) not in done else 'x' if done[(0.5*i, 0.5*j)] >= -0.002 else '.') for i in range(41)))
