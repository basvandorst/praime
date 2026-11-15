"""Independent check: Lucas-Lehmer test of 2^p - 1 for every prime p < X (gmpy2), prints the exponents found."""
import sys, gmpy2
from sympy import primerange
X = int(sys.argv[1])
out = [2]
for p in primerange(3, X):
    M = (gmpy2.mpz(1) << p) - 1; s = gmpy2.mpz(4)
    for _ in range(p - 2):
        s = (s * s - 2) % M
    if s == 0: out.append(p)
print(len(out), out)
