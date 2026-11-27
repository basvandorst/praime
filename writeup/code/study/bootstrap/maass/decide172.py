# Decisions N1-N3 of round 172 (the self-contained chain).  Usage: decide172.py <X>
import sys, json, os, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass')
X = sys.argv[1]; key = [f for f in json.load(open(os.path.join(D, 'lmfdb_key.json'))) if f['r'] < 20]
found = json.load(open(os.path.join(D, 'chain_X%s_found.json' % X))); rows = [json.loads(l) for l in open(os.path.join(D, 'chain_X%s.jsonl' % X))]
print('== self-contained chain, X = %s: %d grid points solved, %d forms recorded, %d true forms below 20' % (X, len(rows), len(found), len(key)))
hit = []; n2 = True; worst_r = 0; worst_a = 0
for g in found:
    c = [f for f in key if (f['parity'] == 'odd') == (g['eps'] == 1)]; f = min(c, key=lambda f: abs(f['r'] - g['r'])); hit.append(f['r'])
    d2 = g['ap']['2'] - f['ap']['2']; d3 = g['ap']['3'] - f['ap']['3']; ok = abs(d2) <= 0.02 and abs(d3) <= 0.02; n2 &= ok
    worst_r = max(worst_r, abs(g['r'] - f['r'])); worst_a = max(worst_a, abs(d2), abs(d3))
    print('  r^ = %.5f %s (true %.5f, difference %+.5f; %d fed-back forms): a_2 %+.4f (true %+.4f), a_3 %+.4f (true %+.4f)  max over p <= 7 of |difference| %.4f, over p = 11, 13 %.4f  %s' % (
        g['r'], 'odd ' if g['eps'] else 'even', f['r'], g['r'] - f['r'], g['known'], g['ap']['2'], f['ap']['2'], g['ap']['3'], f['ap']['3'],
        max(abs(g['ap'][p] - f['ap'][p]) for p in '2357'), max(abs(g['ap'][p] - f['ap'][p]) for p in ('11', '13')), 'ok' if ok else 'N2 no'))
n1 = len(found) == 10 and sorted(hit) == sorted(f['r'] for f in key) and worst_r <= 0.01
print('N1 (ten forms, each within 0.01 of a true eigenvalue of its parity, all ten hit once): largest difference %.5f -> %s' % (worst_r, 'PASS' if n1 else 'FAIL'))
print('N2 (a_2, a_3 within 0.02): largest difference %.4f -> %s' % (worst_a, 'PASS' if n2 else 'FAIL'))
far = [(d['r'], d['eps']) for d in rows if d['t_lp'] >= -0.002 and not any(f['r'] - 0.02 - 1e-9 <= d['r'] <= f['r'] + 0.15 + 1e-9 for f in key if (f['parity'] == 'odd') == (d['eps'] == 1))]
feas = sorted((d['eps'], d['r']) for d in rows if d['t_lp'] >= -0.002); pieces = []
for e, r in feas:
    if pieces and pieces[-1][0] == e and r - pieces[-1][2] < 0.0101: pieces[-1][2] = r
    else: pieces.append([e, r, r])
print('feasible grid points: %d, in %d pieces: %s' % (len(feas), len(pieces), '; '.join('%s [%.2f, %.2f]' % ('odd' if e else 'even', a, b) for e, a, b in sorted(pieces, key=lambda x: x[1]))))
print('N3 (no feasible grid point outside [r_j - 0.02, r_j + 0.15]): %d such points -> %s' % (len(far), 'PASS' if not far else 'FAIL'))
