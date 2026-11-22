"""Independent check of round 45's hit(y) near 1e12: 200,000 consecutive primes from a random start (gmpy2 next_prime /
is_prime, gcd with the primorial of y instead of a sieve array). hit(y) = share of gaps whose composites all share a
factor with prod_{q<=y} q."""
import random, json, math
import gmpy2
random.seed(4512)
p = gmpy2.next_prime(10 ** 12 + random.randrange(10 ** 9))
YS = (2, 7, 97, 997)
prim = {}
for y in YS:
    P = gmpy2.mpz(1)
    for q in range(2, y + 1):
        if gmpy2.is_prime(q):
            P *= q
    prim[y] = P
hits = {y: 0 for y in YS}; n = 0
for _ in range(200000):
    q = gmpy2.next_prime(p)
    for y in YS:
        if all(gmpy2.gcd(m, prim[y]) > 1 for m in range(p + 1, q)):
            hits[y] += 1
    p = q; n += 1
out = {str(y): hits[y] / n for y in YS}
print(json.dumps({'start_near': '1e12', 'gaps': n, 'hit': out}))
