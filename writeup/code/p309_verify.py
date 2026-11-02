"""Puzzle 309: check the linear identity c1*c2*(f(x) - x) = N - D*p on 300 random instances and the parity lemma."""
import random
from math import gcd
from sympy import totient as phi, divisor_sigma as sigma, isprime, nextprime
random.seed(1)
n = 0
while n < 300:
    Q = random.randrange(1, 300); c1 = random.randrange(2, 60, 2); c2 = random.randrange(2, 60, 2)
    r1 = nextprime(random.randrange(10**5, 10**7))
    p = c1 * r1 + 1
    if not isprime(p) or (p + 1) % c2 or not isprime((p + 1) // c2) or gcd(Q, p) > 1:
        continue
    r2 = (p + 1) // c2
    M1, M2 = int(phi(Q)) * c1, int(sigma(Q)) * c2
    if M1 % r1 == 0 or M2 % r2 == 0:
        continue
    x = Q * p
    f = int(sigma(phi(x))) - int(phi(sigma(x)))
    lin = int(sigma(M1)) * (r1 + 1) - int(phi(M2)) * (r2 - 1)
    assert f == lin
    D = c1 * c2 * Q + c1 * int(phi(M2)) - c2 * int(sigma(M1))
    Nn = c2 * int(sigma(M1)) * (c1 - 1) + c1 * int(phi(M2)) * (c2 - 1)
    assert c1 * c2 * (f - x) == Nn - D * p
    n += 1
print('identity holds on', n, 'random instances')
for x in range(3, 20000):
    f = int(sigma(phi(x))) - int(phi(sigma(x)))
    s = int(phi(x))
    square_like = any(s == k * k or s == 2 * k * k for k in range(1, int(s ** 0.5) + 2))
    assert (f % 2 == 1) == square_like
    if f == x:
        assert x == 5
print('parity lemma: f(x) odd iff phi(x) is a square or twice a square (x < 20000); only x = 5 solves it there')
