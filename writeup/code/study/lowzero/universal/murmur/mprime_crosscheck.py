import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from densities import M_prime, M_prime_bessel, C0, nu_prime
print("C0 =", C0, " nu'(1) =", nu_prime(1), " nu'(2) =", nu_prime(2), " nu'(6) =", nu_prime(6))
ys = np.array([0.1, 0.2, 0.3, 0.7, 1.3, 2.2, 3.7])
a = M_prime(ys)
for sm in (2000, 8000):
    b = M_prime_bessel(ys, smax=sm, dmax=200)
    print("smax", sm, " closed:", np.round(a, 4), "\n            bessel:", np.round(b, 4))
