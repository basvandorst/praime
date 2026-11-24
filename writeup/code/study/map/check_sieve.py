"""independent check of map.py's sieve: gmpy2.is_prime on every odd integer of short ranges (no window statistic)"""
import sys, math, gmpy2
sys.path.insert(0, '.')
import map as M
for lo in (6 * 10 ** 9, 4 * 10 ** 10, 10 ** 12 - 10 ** 5):
    hi = lo + 2 * 10 ** 5
    P = M.sieve(lo, hi, M.small_primes(math.isqrt(hi) + 2)).tolist()
    G = [n for n in range(lo | 1, hi, 2) if gmpy2.is_prime(n)]
    print('[%d, %d): sieve %d primes, gmpy2 %d, identical: %s' % (lo, hi, len(P), len(G), P == G))
