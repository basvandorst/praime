"""List the Gallagher intervals with 5 primes (sympy.isprime)."""
import math, sympy
n = 100_000
while True:
    L = round(math.log(n))
    if n + L - 1 > 500_000: break
    ps = [m for m in range(n, n + L) if sympy.isprime(m)]
    if len(ps) >= 5: print(n, L, ps, [x - ps[0] for x in ps])
    n += L
