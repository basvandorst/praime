# Discovery look (not a test): do the Artin coefficients a_K(p) = (#roots of f mod p) - 1 of S3 cubic fields, averaged over
# the fields with X <= |D| < 2X, oscillate as a function of y = p/|D|?  Compare with J_0(4 pi sqrt y) (weight-1 kernel, D < 0).
import sys, math, numpy as np, cypari2
from scipy.special import j0, y0, k0
pari = cypari2.Pari(); pari.allocatemem(600_000_000)
X = int(sys.argv[1]); YMAX = float(sys.argv[2])
for sgn, name in ((-1, 'complex (D < 0)'), (1, 'totally real (D > 0)')):
    fl = pari('nflist("S3", [%d, %d], %d)' % (X, 2 * X - 1, 1 if sgn < 0 else 0))
    nb = int(YMAX / 0.05); num = np.zeros(nb); cnt = np.zeros(nb); nf = 0
    pari('cnt(f, P) = apply(p -> #polrootsmod(f, p) - 1, P)')
    for f in fl:
        D = abs(int(pari.nfdisc(f))); P = pari.primes([2, int(YMAX * D)])
        a = np.array([int(v) for v in pari('cnt')(f, P)], float); p = np.array([int(v) for v in P], float)
        k = (D % p.astype(np.int64) != 0)
        b = np.minimum((p[k] / D / 0.05).astype(int), nb - 1)
        np.add.at(num, b, a[k] * math.sqrt(D)); np.add.at(cnt, b, 1.0); nf += 1
    m = num / cnt; err = math.sqrt(X * 1.5) / np.sqrt(cnt)        # a has variance about 1
    y = (np.arange(nb) + 0.5) * 0.05
    ker = j0(4 * math.pi * np.sqrt(y)) if sgn < 0 else (-y0(4 * math.pi * np.sqrt(y)) + (2 / math.pi) * k0(4 * math.pi * np.sqrt(y)))
    amp = float(np.sum(m * ker / err**2) / np.sum(ker**2 / err**2)); aerr = float(1 / math.sqrt(np.sum(ker**2 / err**2)))
    chi0 = float(np.sum((m / err) ** 2)); chi1 = float(np.sum(((m - amp * ker) / err) ** 2))
    print('%s: %d fields, %d bins; chi2 against zero %.0f, after the kernel fit %.0f; fitted amplitude %.3f +- %.3f' % (name, nf, nb, chi0, chi1, amp, aerr))
    print('   y   : ' + ' '.join('%5.2f' % v for v in y[:16]))
    print('   mean: ' + ' '.join('%5.2f' % v for v in m[:16]))
    print('   fit : ' + ' '.join('%5.2f' % v for v in (amp * ker)[:16]))
