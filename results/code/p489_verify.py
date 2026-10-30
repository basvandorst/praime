"""Puzzle 489: checks for results/489.md (independent of the C search)."""
from sympy import divisor_sigma as s, totient as phi, isprime, factorint

f = lambda n: s(s(phi(n)))

# 1. Direct search for n = q + 1 <= 2*10^5 reproduces the page's list.
sols = [n for n in range(1, 200001) if f(n) == n]
assert sols == [1, 4, 8, 24, 32, 128, 1536, 8192, 131072], sols

# 2. Theorem on n = 2^a 3^b, checked for a <= 64, b <= 30:
#    solutions are exactly 1, 24, 1536 and 2^a with 2^a - 1 prime.
for a in range(0, 65):
    for b in range(0, 31):
        n = 2**a * 3**b
        expected = n in (1, 24, 1536) or (b == 0 and a >= 2 and isprime(2**a - 1))
        assert (f(n) == n) == expected, (a, b)

# 3. Key lemma: sigma(2^m - 1) = 3 * 2^(m-1) only for m = 4 and m = 10 (checked m <= 100).
assert [m for m in range(1, 101) if s(2**m - 1) == 3 * 2**(m - 1)] == [4, 10]

# 4. The inequality step: for a >= 1, b >= 2, M = A*B with A = 2^(a+1)-1, B = (3^b-1)/2,
#    n - M = B + (A+1)/2 while sigma(M) - M >= 1 + A + B.
for a in range(1, 40):
    for b in range(2, 25):
        A, B, n = 2**(a+1) - 1, (3**b - 1)//2, 2**a * 3**b
        assert A != B and n - A*B == B + (A + 1)//2 and s(A*B) - A*B >= 1 + A + B

# 5. Firoozbakht's family: 2^p - 1 prime and sigma(2^p + 1) = 3*2^(p-1); only p = 2, 5
#    among Mersenne exponents p <= 127 (the proof covers all p).
mers = [p for p in range(2, 128) if isprime(p) and isprime(2**p - 1)]
assert [p for p in mers if s(2**p + 1) == 3 * 2**(p - 1)] == [2, 5]
print("all checks pass")
