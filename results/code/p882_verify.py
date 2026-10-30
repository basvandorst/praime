"""Puzzle 882: check the run 95773717(6)_i1, i = 0..10, in results/882.md."""
from sympy import isprime, factorint

N = lambda A, Z, B, i: int(str(A) + str(Z) * i + str(B))

# Our run: 11 primes, all below 2^64, so sympy.isprime is a proof.
run = [N(95773717, 6, 1, i) for i in range(11)]
assert all(isprime(x) for x in run) and max(run) < 2**64
print(run[0], '...', run[-1])

# It stops at i = 11.
n11 = N(95773717, 6, 1, 11)
assert not isprime(n11)
print('i = 11:', n11, '=', factorint(n11))

# Taneja's example from the page: 291234(3)_i73, i = 0..9, and i = 10 composite.
assert all(isprime(N(291234, 3, 73, i)) for i in range(10)) and not isprime(N(291234, 3, 73, 10))

# Mod 3: N_i = digitsum(A) + digitsum(B) + i*Z (mod 3), so Z not = 0 mod 3 puts a multiple of 3 among i = 0, 1, 2.
for Z in range(10):
    for A in range(1, 300):
        for B in (1, 3, 7, 9):
            res = [N(A, Z, B, i) % 3 for i in range(3)]
            if Z % 3:
                assert 0 in res
print("all checks pass")
