# Puzzle 1223 (writeup/1223.md): F(n, k) = p_n^(p_k - 1) + p_k - 1. Primes: (2, 2) -> 3, (3, 3) -> 11, (31, 31).
# Step 1 (proved): if p_n != p_k then p_k | F (Fermat) and F > p_k, so F is composite; so only n = k matters,
#   G(p) = p^(p-1) + p - 1.
# Step 2: G(p) is never divisible by 2 or 3 (p > 3): checked below and explained in the write-up.
# Step 3 (heuristic): for a prime q >= 5 the fraction f_q of primes p with q | G(p) follows from (p mod q, p mod q-1),
#   both uniform over units (CRT, gcd(q, q-1) = 1): f_q = #{(r, e): r^(e-1) = 1 - r mod q} / ((q-1) phi(q-1)).
#   P(G(p) prime) = C / ln G(p) with C = 2 * 3/2 * prod_{q >= 5} (1 - f_q)/(1 - 1/q); expected number of prime G(p)
#   over p > X is about C / ln X (sum of 1/((p-1) ln p)), finite.
# Step 4: calibration on the sieve (survivors of q <= 300 among primes p <= 20000) and the expected counts.
# Run: python3 p1223_model.py   (about a minute; needs sympy, gmpy2)
from math import log, gcd
from sympy import primerange, totient
import gmpy2

G = lambda p: p ** (p - 1) + p - 1
for p in primerange(5, 3000): assert G(p) % 2 and G(p) % 3
print('G(p) is odd and not divisible by 3 for 5 <= p < 3000')
pr = [p for p in primerange(2, 400) if gmpy2.is_prime(G(p))]
print('primes G(p) for p < 400:', pr)

def density(q):
    units_e = [e for e in range(q - 1) if gcd(e, q - 1) == 1]
    cnt = 0
    for r in range(1, q):
        for e in units_e:
            if pow(r, (e - 1) % (q - 1), q) == (1 - r) % q: cnt += 1
    return cnt / ((q - 1) * len(units_e))
# brute-force check of the density on actual primes for q < 60
for q in primerange(5, 60):
    ps = list(primerange(q + 1, 200000))
    emp = sum(1 for p in ps if (pow(p, p - 1, q) + p - 1) % q == 0) / len(ps)
    assert abs(emp - density(q)) < 0.02, (q, emp, density(q))
print('density formula agrees with the primes below 200000 for 5 <= q < 60')
C = 2 * 1.5; qs = list(primerange(5, 1500)); f = {}
for q in qs:
    f[q] = density(q); C *= (1 - f[q]) / (1 - 1 / q)
    if q in (97, 499, 997, 1499): print('  C (product to q = %d): %.4f' % (q, C))
# calibration: survivors of the sieve by 5 <= q <= 300 among primes 5 <= p <= 20000
ps = list(primerange(5, 20001)); small = [q for q in qs if q <= 300]
surv = sum(1 for p in ps if all((pow(p, p - 1, q) + p - 1) % q for q in small if q != p))
pred = len(ps)
for q in small: pred *= 1 - f[q]
print('calibration: %d of %d primes p <= 20000 have G(p) free of the primes 5..300; predicted %.1f (ratio %.3f)' % (surv, len(ps), pred, surv / pred))
def expected(lo, hi):
    return sum(C / ((p - 1) * log(p)) for p in primerange(lo, hi))
print('expected number of prime G(p): p in [5, 400): %.2f (found %d); p in [400, 27449]: %.3f (found 0)' % (expected(5, 400), len([p for p in pr if p >= 5]), expected(400, 27450)))
tail = C / log(27449)            # sum_{p > X} 1/((p-1) ln p) ~ 1/ln X
print('expected beyond p = 27449 (the 3000th prime, searched on the page): about %.3f; chance of no further prime: %.0f %%' % (tail, 100 * 2.718281828 ** -tail))
print('all checks passed')
