"""Round 128: zeros of the real characters with prime conductor q (D = q or -q) at heights where 7.6 <= ln(q t/2pi) <= 8.6.
One file per q in lzeros/r128/. usage: r128z.py PART NPARTS"""
import sys, os, math
import numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros', 'r128')
part, nparts = int(sys.argv[1]), int(sys.argv[2])
QS = [29, 31, 37, 41, 43, 47, 53, 59, 61, 67]
for q in QS[part::nparts]:
    D = q if q % 4 == 1 else -q
    f = os.path.join(OUT, 'D%d.npy' % D)
    if os.path.exists(f): continue
    a, b = 2 * math.pi * math.exp(7.6) / q, 2 * math.pi * math.exp(8.6) / q
    np.save(f, np.array(sorted(float(t) for t in pari('lfunzeros(lfuncreate(%d),[%f,%f])' % (D, a, b)))))
    print(D, flush=True)
