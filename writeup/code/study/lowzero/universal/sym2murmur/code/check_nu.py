"""Route C check of nu(m): average of pi H(D)/sqrt(D), D = m (4p - m N^2), over primes N in [N0,N1] and primes p with
m N^2/4 < p <= 3 N^2 (the class numbers of the N^2-terms of the trace formula), against nu(m) of model.py."""
import sys
sys.path.insert(0, sys.argv[1])
import numpy as np
from sympy import primerange
from hurwitz import hurwitz12
from model import nu
N0, N1 = 101, 199
Ns = list(primerange(N0, N1 + 1))
H = hurwitz12(16 * 12 * N1 * N1 + 10)
P = np.array(list(primerange(2, 3 * N1 * N1 + 1)), dtype=np.int64)
print(" m   route-C average (se)     nu(m) [Petersson lines]   ratio")
for m in range(1, 17):
    vals = []
    for N in Ns:
        pp = P[(4 * P > m * N * N) & (P <= 3 * N * N) & (P != N)]
        D = m * (4 * pp - m * N * N)
        vals.append(np.pi * H[D] / 12.0 / np.sqrt(D))
    v = np.concatenate(vals)
    print("%2d   %.4f (%.4f)   %.4f   %.4f" % (m, v.mean(), v.std() / np.sqrt(len(v)), nu(m), v.mean() / nu(m)))
