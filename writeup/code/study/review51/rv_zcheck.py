# Accuracy of stored zeros: Newton offset delta = Z(gamma)/Z'(gamma) with PARI lfunhardy, for a sample of zeros.
import os, sys, time
import numpy as np
import cypari2
pari = cypari2.Pari(); pari.allocatemem(800 * 10**6)
Z = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/study/r144/fresh_zeros.npz')
for s in sys.argv[1].split(','):
    t0 = time.time()
    a = [int(x.replace('m', '-')) for x in s.split('_')[1:]]
    pari('default(realprecision, 30)')
    pari('global(LL); LL = lfuninit(ellinit(%s), [1/2, 0, 405])' % str(a))
    z = np.sort(Z[s]); band = z[(z > 20) & (z <= 400)]
    sel = band[np.linspace(0, len(band) - 1, 25).astype(int)]
    d = []
    for g in sel:
        g = float(g)
        h = 1e-4
        zp = float(pari('lfunhardy(LL, %.12f)' % (g + h))); zm = float(pari('lfunhardy(LL, %.12f)' % (g - h))); z0 = float(pari('lfunhardy(LL, %.12f)' % g))
        d.append(z0 / ((zp - zm) / (2 * h)))
    d = np.array(d)
    print('%s: 25 zeros in band, Newton offsets |delta|: median %.2e max %.2e  (%.0fs)' % (s, np.median(np.abs(d)), np.abs(d).max(), time.time() - t0))
