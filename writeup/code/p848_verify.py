"""Puzzle 848: p = 95521311337, q = 7: p 0^n q prime for n = 0..9 (k = 9), with primality proofs."""
from sympy import isprime, factorint

def proven_prime(N):
    """Lucas (Pocklington) proof from the full factorisation of N - 1; recursive for factors above 2^64."""
    if N < 2**64:
        return isprime(N)                      # deterministic below 2^64
    fac = factorint(N - 1)
    for q in fac:
        if not proven_prime(q):
            return False
    for a in range(2, 200):
        if pow(a, N - 1, N) != 1:
            return False
        if all(pow(a, (N - 1) // q, N) != 1 for q in fac):
            return True                         # a has order N - 1, so N is prime
    raise RuntimeError("no witness found")

p, q = 95521311337, 7
chain = [int(str(p) + "0" * n + str(q)) for n in range(10)]
assert isprime(p) and all(proven_prime(N) for N in chain)
n10 = int(str(p) + "0" * 10 + str(q))
assert not isprime(n10)
print("k = 9:", chain[0], "...", chain[-1], "all proven; n = 10:", factorint(n10))

# Emmanuel Vantieghem's k = 8 record reproduces: 7 0^n 313292827.
assert all(isprime(int("7" + "0" * n + "313292827")) for n in range(9))

# Mod 7 (from the page): if neither p nor q is 7, one of n = 0..5 gives a multiple of 7.
for P in (11, 13, 17, 23, 29, 31):
    for Q in (3, 11, 13, 19, 23):
        assert any(int(str(P) + "0" * n + str(Q)) % 7 == 0 for n in range(6))

# Small x missed by the sieve (x <= 2741 is itself one of the sieving primes): brute force x < 3000.
small = [x for x in range(2, 3000) if isprime(x) and
         (all(isprime(int("7" + "0" * n + str(x))) for n in range(10)) or
          all(isprime(int(str(x) + "0" * n + "7")) for n in range(10)))]
assert small == [], small
print("all checks pass")
