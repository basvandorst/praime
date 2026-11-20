"""Independent check of euler.py: count n < N with n^2 + n + 41 prime by gmpy2.is_prime (no sieve), at N = 10^2..10^6.
usage: euler_check.py N out.txt"""
import sys, gmpy2
N = int(sys.argv[1]); out = open(sys.argv[2], 'w')
c = 0; marks = {10 ** e for e in range(2, 7)}
for n in range(N):
    c += gmpy2.is_prime(n * n + n + 41, 30)
    if n + 1 in marks:
        s = 'n < %d: %d primes' % (n + 1, c); print(s); out.write(s + '\n')
