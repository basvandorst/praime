"""Round 120: zeros (<= 600) of the last third (offset 2) of the fundamental discriminants 100 <= |D| <= 400.
usage: r120z.py PART NPARTS"""
import sys, os, json
import numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
part, nparts = int(sys.argv[1]), int(sys.argv[2])
rows = sorted([x for x in json.load(open(os.path.join(H, 'r104.json')))['rows'] if 100 <= abs(x['D']) <= 400], key=lambda x: (abs(x['D']), x['D']))[2::3]
Z = {}
for r in rows[part::nparts]:
    Z['D%d' % r['D']] = np.array(sorted(float(t) for t in pari.lfunzeros(pari.lfuncreate(r['D']), 600) if float(t) > 0))
np.savez_compressed(os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros', 'r120_part%d.npz' % part), **Z)
print(len(Z), 'done')
