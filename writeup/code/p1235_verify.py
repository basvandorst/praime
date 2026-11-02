# Puzzle 1235 (writeup/1235.md): checks on F(n) = 6^n + n^6.
# 1. F(n) is composite for 2 <= n <= 3000 (the page: to 150000), and F(1) = 7 is prime.
# 2. The small-prime classes: 2 | F(n) iff n even; 3 | F(n) iff 3 | n; 5 | F(n) iff n = +-2 (mod 5); 7 | F(n) iff
#    7 does not divide n and n is odd (i.e. the allowed odd n are 7 (mod 14)); the survivors mod 210 are
#    35, 49, 91, 119, 161, 175.
# 3. No finite set of primes divides every F(n) with n > 1: for any primes q_1..q_k >= 11, n = 35 (mod 210) and
#    n = 0 (mod q_1...q_k) gives F(n) = 6^n (mod q_i), not 0. Checked for the primes 11..97 by exhibiting such n.
# Run: python3 p1235_verify.py   (about 10 seconds; needs gmpy2, sympy)
import gmpy2
from sympy import primerange
from sympy.ntheory.modular import crt

F = lambda n: 6 ** n + n ** 6
assert gmpy2.is_prime(F(1))
for n in range(2, 3001): assert not gmpy2.is_prime(F(n)), n
print('F(1) = 7 is prime; F(n) is composite for 2 <= n <= 3000')
for n in range(1, 2000):
    assert (F(n) % 2 == 0) == (n % 2 == 0)
    assert (F(n) % 3 == 0) == (n % 3 == 0)
    assert (F(n) % 5 == 0) == (n % 5 in (2, 3))
    assert (F(n) % 7 == 0) == (n % 7 != 0 and n % 2 == 1)
surv = sorted({n % 210 for n in range(210, 420) if all(F(n) % q for q in (2, 3, 5, 7))})
assert surv == [35, 49, 91, 119, 161, 175]
print('small-prime rules checked for n < 2000; survivors mod 210:', surv)
qs = list(primerange(11, 100))
M = 1
for q in qs: M *= q
n = int(crt([210, M], [35, 0])[0])
assert all((pow(6, n, q) + pow(n, 6, q)) % q for q in (2, 3, 5, 7)) and all((pow(6, n, q) + pow(n, 6, q)) % q for q in qs)
print('n = %d is 35 (mod 210) and 0 (mod 11*13*...*97): F(n) has no prime factor below 100, so no covering set' % n)
print('all checks passed')
