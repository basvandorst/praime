"""Round 117: recompute and save the zeros (<= 600) of the 121 L-functions of rounds 115 and 116. usage: r117z.py PART NPARTS"""
import sys, os, json
import numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
part, nparts = int(sys.argv[1]), int(sys.argv[2])
Ds = [r['D'] for r in json.load(open(os.path.join(H, 'r115.json')))['rows']] + [r['D'] for r in json.load(open(os.path.join(H, 'r116.json')))['rows']]
Ds = sorted(Ds)[part::nparts]; Z = {}
for D in Ds:
    Z['D%d' % D] = np.array(sorted(float(t) for t in pari.lfunzeros(pari.lfuncreate(D), 600) if float(t) > 0))
np.savez_compressed(os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros', 'r117_part%d.npz' % part), **Z)
print(len(Ds), 'done')
