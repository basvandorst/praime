"""Second check of the validation sieve: gmpy2/sympy on random samples and on the gaps."""
import numpy as np, gmpy2, sympy, random
P = np.load('data/val_primes.npy'); X0 = 150_000_000_000
random.seed(1)
idx = random.sample(range(len(P) - 1), 3000)
bad = sum(1 for i in idx if not sympy.isprime(int(P[i])))
# no prime missed between consecutive listed primes (sample of 3000 gaps)
miss = 0
for i in idx:
    a, b = int(P[i]), int(P[i + 1])
    n = gmpy2.next_prime(a)
    if int(n) != b: miss += 1
print('sample of 3000: non-primes %d, gaps with a missed prime %d; first prime %d = next_prime(X0-1)? %s'
      % (bad, miss, P[0], int(gmpy2.next_prime(X0 - 1)) == int(P[0])))
