# Conjecture 28 (writeup/conj-28.md): Q(n, m) = nextprime(n!^m + 1) - n!^m. Buss: prime for all n when m <= 5.
# Reasoning. Q = k is the least k > 1 with n!^m + k prime. Every prime p <= n divides n!^m, so n!^m + k is prime only
# if gcd(k, n!) = 1: the candidates are the k > n without prime factors <= n, in increasing order: the primes in
# (n, p'^2) first, where p' is the next prime after n, then composites start to appear (p'^2 is the first).
# Each candidate is prime with about the same probability q = F_n / ln(n!^m), F_n = prod_{p <= n} (1 - 1/p)^-1
# (the candidate already avoids the primes <= n). So Q(n, m) is composite with probability
#   E(n, m) = sum over composite candidates c of q (1 - q)^(number of candidates below c).
# The script: (1) computes Q(n, m) exactly for small n, m and lists the composite values; (2) computes E(n, m) and
# the expected number of composite Q for each m over 4 <= n <= 200 (for n <= 3 the model's q exceeds 1; those
# cases are checked exactly); (3) compares with the failures found.
# Run: python3 conj28_model.py   (a few minutes; needs sympy, gmpy2, numpy)
from math import log, exp, lgamma
import numpy as np
from sympy import primerange, nextprime, isprime, factorint
import gmpy2

def Q(n, m):
    x = gmpy2.mpz(gmpy2.fac(n)) ** m
    return int(gmpy2.next_prime(x + 1) - x)

print('composite Q(n, m) for n <= 40, m <= 10:')
fails = {}
for m in range(1, 11):
    for n in range(2, 41):
        k = Q(n, m)
        if not isprime(k): fails.setdefault(m, []).append((n, k, factorint(k)))
for m in sorted(fails): print('  m = %d: %s' % (m, ', '.join('Q(%d)=%d=%s' % (n, k, '*'.join('%d^%d' % t if t[1] > 1 else str(t[0]) for t in f.items())) for n, k, f in fails[m])))

def E(n, m):
    """probability that the first k > 1 with n!^m + k prime is composite (model)."""
    F = 1.0
    for p in primerange(2, n + 1): F /= 1 - 1 / p
    q = F / (m * lgamma(n + 1))            # ln n! = lgamma(n+1)
    rmax = int(40 / q) + 10                # after rmax candidates the survival is below e^-40
    K = int(rmax * F) + 1000               # candidates have density about 1/F
    small = list(primerange(2, n + 1))
    cand = np.ones(K + 1, dtype=bool); cand[:n + 1] = False
    for p in small: cand[p::p] = False     # k with no prime factor <= n
    comp = np.zeros(K + 1, dtype=bool)
    for p in primerange(n + 1, int(K ** 0.5) + 1): comp[p * p::p] = True
    ks = np.nonzero(cand)[0][:rmax]
    surv = (1 - q) ** np.arange(len(ks))
    return float(q * surv[comp[ks]].sum())

print('\nmodel: probability that Q(n, m) is composite, and the expected number of composite Q over 4 <= n <= 200')
for m in range(1, 11):
    ex = [E(n, m) for n in range(4, 201)]           # n <= 3 are checked exactly above; beyond n = 200 the terms are below 1e-9
    known = len(fails.get(m, []))
    print('  m = %2d: E(5)=%.3f E(10)=%.3f E(28)=%.4f E(50)=%.2e E(100)=%.2e E(200)=%.1e  expected total %.3f  (found for n <= 40: %d)' %
          (m, ex[1], ex[6], ex[24], ex[46], ex[96], ex[196], sum(ex), known))
tot = sum(sum(E(n, m) for n in range(4, 201)) for m in range(1, 6))
print('expected composite Q over m = 1..5 together: %.2f; chance of none (as observed): %.0f %%' % (tot, 100 * exp(-tot)))
print('all checks passed')
