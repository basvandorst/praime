# Comparison of a chain's recorded forms with the answer keys (rounds 173 and 174).  Usage: decide173.py <found.json> <grid.jsonl> [r_lo]
import sys, json, os, numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass')
key = [f for f in json.load(open(os.path.join(D, 'lmfdb_key.json'))) if f['r'] < 20] + json.load(open(os.path.join(D, 'key_above20.json')))
found = json.load(open(os.path.join(D, sys.argv[1]))); rows = [json.loads(l) for l in open(os.path.join(D, sys.argv[2]))]; rlo = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
rmax = max(d['r'] for d in rows)
print('== %s: %d forms recorded; grid up to %.2f; true forms up to 22.2: %d' % (sys.argv[1], len(found), rmax, len(key)))
for g in found:
    c = [f for f in key if (f['parity'] == 'odd') == (g['eps'] == 1)]; f = min(c, key=lambda f: abs(f['r'] - g['r']))
    a = {int(p): v for p, v in g['ap'].items()}; ta = {int(p): v for p, v in f['ap'].items()}
    print('  r^ = %.5f %s (true %.5f, difference %+.5f; %s; %d fed-back): a_2 %+.4f (%+.4f), a_3 %+.4f (%+.4f), a_5 %+.4f (%+.4f), a_7 %+.4f (%+.4f); largest |difference| p <= 3: %.4f, p <= 7: %.4f, p = 11, 13: %.4f' % (
        g['r'], 'odd ' if g['eps'] else 'even', f['r'], g['r'] - f['r'], g['how'], g['known'], a[2], ta[2], a[3], ta[3], a[5], ta[5], a[7], ta[7],
        max(abs(a[p] - ta[p]) for p in (2, 3)), max(abs(a[p] - ta[p]) for p in (2, 3, 5, 7)), max(abs(a[p] - ta[p]) for p in (11, 13))))
missed = [f for f in key if rlo < f['r'] < rmax - 0.1 and not any((f['parity'] == 'odd') == (g['eps'] == 1) and abs(g['r'] - f['r']) <= 0.02 for g in found)]
print('true eigenvalues in (%.2f, %.2f) with no recorded form within 0.02: %s' % (rlo, rmax - 0.1, ', '.join('%.4f %s' % (f['r'], f['parity']) for f in missed) or 'none'))
feas = sorted((d['eps'], d['r']) for d in rows if d['t_lp'] >= -0.002); pieces = []
for e, r in feas:
    if pieces and pieces[-1][0] == e and r - pieces[-1][2] < 0.0101: pieces[-1][2] = r
    else: pieces.append([e, r, r])
print('feasible pieces: %s' % '; '.join('%s [%.2f, %.2f]' % ('odd' if e else 'even', a, b) for e, a, b in sorted(pieces, key=lambda x: x[1])))
