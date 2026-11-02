"""Puzzle 434: check the lemmas of 434.md on exact data."""
from sympy import divisor_sigma as sigma, isprime, factorint
from fractions import Fraction as Fr

A = lambda n: Fr(int(sigma(n)), n)
for k in range(2, 128):
    if isprime(2**k - 1):
        p = 2**k - 1
        assert sigma(sigma(p)) == sigma(p) + p
for a in range(1, 30):
    q_max = 2**(a + 1) - 1
    assert A(2**a) * Fr(q_max + 1, q_max) == 2
assert all(A(6 * k) >= 2 for k in range(1, 5000))
for k in range(1, 200):
    t = 2**k
    assert 2**k * (2**(k + 1) + 5) - 3 == (2 * t - 1) * (t + 3)
    q = 2**(k + 1) + 5
    if isprime(q):
        x = 2**k * q
        assert sigma(x) == 2 * x - 6
for x in range(1, 200000):
    y = int(sigma(x))
    if sigma(y) == y + x:
        assert A(y) == 1 + 1 / A(x) and A(y) < 2
        assert isprime(x) and x + 1 == 2**(x.bit_length())
x, p, r = 15633, 193, 97
assert x == 81 * p and p + 1 == 2 * r
assert int(sigma(sigma(x))) - int(sigma(x)) - x == -5 * r + 480
print('lemmas checked: Mersenne solutions, A(2^a (2^(a+1)-1)) = 2, 6k abundant, (2t-1)(t+3) factorisation, only Mersenne below 2e5, the 15633 family line')
