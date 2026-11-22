"""Round 115 (preregister/round115.txt): spacing variance of the zeros of L(s, chi_D), 100 <= gamma <= 600, for every
third fundamental discriminant 100 <= |D| <= 400, against S (round 104). usage: r115.py PART NPARTS"""
import sys, os, json, math
import numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
part, nparts = int(sys.argv[1]), int(sys.argv[2])
rows = sorted([x for x in json.load(open(os.path.join(H, 'r104.json')))['rows'] if 100 <= abs(x['D']) <= 400], key=lambda x: (abs(x['D']), x['D']))[::3]
out = []
for r in rows[part::nparts]:
    D = r['D']; q = abs(D)
    g = np.array(sorted(float(t) for t in pari.lfunzeros(pari.lfuncreate(D), 600) if float(t) > 0))
    N = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e)); s = np.diff(N); mid = 0.5 * (g[1:] + g[:-1])
    x = s[mid >= 100]
    out.append({'D': D, 'nzeros': int(len(g)), 'nspacings': int(len(x)), 'var': float(np.var(x)), 'mean': float(np.mean(x))})
    print(D, len(g), round(float(np.var(x)), 4), flush=True)
json.dump(out, open(os.path.join(H, 'r115_part%d.json' % part), 'w'))
