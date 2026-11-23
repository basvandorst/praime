"""Round 121: zeros (100 <= |gamma| <= 500) of complex primitive characters, conductors N = prime, 3p, 4p, 5p, 7p in
[100, 400]. usage: r121z.py PART NPARTS"""
import sys, os
import numpy as np, cypari2
from sympy import primerange
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
part, nparts = int(sys.argv[1]), int(sys.argv[2])
Ns = list(primerange(101, 398))[::2]
for c, bad in ((3, 3), (4, 2), (5, 5), (7, 7)):
    Ns += [c * p for p in primerange(3, 400) if p != bad and 100 <= c * p <= 400]
Ns = sorted(set(Ns))
Z = {}
for N in Ns[part::nparts]:
    G = pari('znstar(%d,1)' % N); m = None
    for a in range(2, N):
        if pari('gcd(%d,%d)' % (a, N)) != 1: continue
        if pari('znconreyconductor(znstar(%d,1),%d)' % (N, a)) == N and int(pari('charorder(znstar(%d,1),znconreylog(znstar(%d,1),%d))' % (N, N, a))) > 2:
            m = a; break
    z = np.array(sorted(float(t) for t in pari('lfunzeros(lfuncreate(Mod(%d,%d)),[-500,500])' % (m, N))))
    Z['N%d_m%d' % (N, m)] = z[np.abs(z) >= 100]
np.savez_compressed(os.path.join(H, '..', '..', '..', 'data', 'riemann', 'lzeros', 'r121_part%d.npz' % part), **Z)
print(len(Z), 'done')
