"""Round 122: zeros with 900 <= gamma <= 1500 of every third real fundamental discriminant 100 <= |D| <= 400 (offset 0),
one file per D in lzeros/r122/ (skips those already done). usage: r122z.py PART NPARTS"""
import sys, os, json
import numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros', 'r122')
part, nparts = int(sys.argv[1]), int(sys.argv[2])
rows = sorted([x for x in json.load(open(os.path.join(H, 'r104.json')))['rows'] if 100 <= abs(x['D']) <= 400], key=lambda x: (abs(x['D']), x['D']))[0::3]
for r in rows[part::nparts]:
    f = os.path.join(OUT, 'D%d.npy' % r['D'])
    if os.path.exists(f): continue
    np.save(f, np.array(sorted(float(t) for t in pari('lfunzeros(lfuncreate(%d),[900,1500])' % r['D']))))
    print(r['D'], flush=True)
