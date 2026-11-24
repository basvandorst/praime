"""Independent check of sieve.py: (1) 2000 random listed numbers per stretch are prime (sympy.isprime, a proof below 2^64);
(2) completeness: in three sub-intervals of 2e5 per stretch every integer is tested with gmpy2.is_prime and the prime set
must equal the listed set."""
import numpy as np, sympy, gmpy2, random
random.seed(5)
for name, x0 in (('A', 150_000_000_000), ('B', 5_000_000_000)):
    P = np.load('primes_%s.npy' % name)
    samp = random.sample(range(len(P)), 2000)
    bad = [int(P[i]) for i in samp if not sympy.isprime(int(P[i]))]
    ok_c = True
    for lo in (x0, x0 + 50_000_000, x0 + 100_000_000 - 200_000):
        hi = lo + 200_000
        true = [n for n in range(lo, hi) if gmpy2.is_prime(n)]
        listed = P[(P >= lo) & (P < hi)].tolist()
        ok_c &= (true == listed)
    print(name, 'sample prime:', len(bad) == 0, ' completeness on 3 x 2e5:', ok_c)
