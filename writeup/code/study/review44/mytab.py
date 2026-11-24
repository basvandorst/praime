"""Per-block transition tables (k = 1) for q = 3, 5, 7 and k = 2 for q = 10, H = 1e4, 1e5, 1e6 (reviewer's own code)."""
import sys, os, math
import numpy as np
SRC = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mywin.py')).read()
exec(SRC[SRC.index('def small_primes'):SRC.index('def block_stats')])
OUT = sys.argv[1]; WIN = [int(float(a)) for a in sys.argv[2:]]
SPAN = 10 ** 8; PAD = 4000
base = small_primes(int(math.isqrt(max(WIN) + SPAN + PAD)) + 2)
for x0 in WIN:
    P = sieve(x0 - PAD, x0 + SPAN + PAD, base); res = {}
    for q, k in ((3, 1), (5, 1), (7, 1), (10, 2)):
        r = red(q); C = len(r); lut = -np.ones(q, np.int64); lut[r] = np.arange(C); cls = lut[P % q]
        for H in (10 ** 4, 10 ** 5, 10 ** 6):
            NB = SPAN // H; blk = np.where(P >= x0, (P - x0) // H, -1)
            a = cls[:-k]; b = cls[k:]; bl = blk[:-k]; ok = (bl >= 0) & (bl < NB)
            res['X_%d_%d_%d' % (q, k, H)] = np.bincount((bl[ok] * C + a[ok]) * C + b[ok], minlength=NB * C * C).reshape(NB, C, C).astype(np.int32)
    np.savez_compressed(os.path.join(OUT, 't_%d.npz' % x0), **res)
    print('window %d' % x0, flush=True)
