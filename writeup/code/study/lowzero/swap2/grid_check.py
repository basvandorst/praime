"""u-grid of the J quadrature: step 0.05 to 6 (frozen) against step 0.025 to 8 and 0.1 to 6, Euler product to 3e4."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hlib import *
pari = pari_instance(); pr = primes_upto(30000)
_, X = real_table(pari, [-995, 997, 401, -403, 892, 321, -104, 104], 30000)
_, c1 = conrey_table(pari, 115, 2, 30000); _, c2 = conrey_table(pari, 388, 3, 30000)
cases = [('-995,997', X[0] * X[1]), ('401,-403', X[2] * X[3]), ('892,321', X[4] * X[5]), ('-104,104', X[6] * X[7]),
         ('115/2 x 388/3', c1 * np.conj(c2)), ('115/2 x conj', c1 * c1)]
for nm, psi in cases:
    psi = psi.astype(complex)
    a = J_pair_direct(psi, pr, (0.5, 1.0))
    b = J_pair_direct(psi, pr, (0.5, 1.0), ug=np.arange(0, 8.0 + 1e-9, 0.025), du=0.025)
    c = J_pair_direct(psi, pr, (0.5, 1.0), ug=np.arange(0, 6.0 + 1e-9, 0.1), du=0.1)
    print(nm, ' '.join('s=%.1f frozen %.7f%+.7fi | fine-long diff %.1e | coarse diff %.1e' % (s, a[s].real, a[s].imag, abs(a[s] - b[s]), abs(a[s] - c[s])) for s in (0.5, 1.0)))
