"""Puzzle 788 Q2: the sequence is finite (at most 168 terms).

N_q = P*10^(L+d) + q^2*10^L + RP, L = digits of P, d = digits of q^2.
Mod 7, 7 | N_q  <=>  q^2 = -(x*10^d + y)  with x = P mod 7, y = RP*10^(-L) mod 7.
The check depends only on (x, y) and d mod 6, so 6*7 states cover every prime P != 7.
"""
from sympy import primerange, prime, isprime

def square_sets(bound):
    s = {d: set() for d in range(6)}
    for q in primerange(2, bound + 1):
        s[len(str(q * q)) % 6].add(q * q % 7)
    return s

def survivors(bound):
    s = square_sets(bound)
    return [(x, y) for x in range(1, 7) for y in range(7)
            if all((-(x * pow(10, d, 7) + y)) % 7 not in s[d] for d in range(6))]

assert prime(168) == 997 and prime(169) == 1009
assert survivors(997) == [(1, 3)], survivors(997)   # k = 168: one class left mod 7
assert survivors(1009) == []                        # k = 169: every P != 7 blocked
print("square sets for q <= 1009 by d mod 6:", {d: sorted(v) for d, v in square_sets(1009).items()})

# direct check on random primes: some q <= 1009 gives 7 | N_q
import random
random.seed(788)
for _ in range(2000):
    P = random.randrange(10, 10**random.randint(2, 30))
    if not isprime(P) or P == 7: continue
    RP = int(str(P)[::-1]); L = len(str(P))
    assert any((P * 10**(L + len(str(q*q))) + q*q * 10**L + RP) % 7 == 0 for q in primerange(2, 1010)), P
# P = 7: 7&4&7 = 747 = 3^2 * 83, fails already at k = 1
assert not isprime(747)
# the known terms still pass (the k = 6 example from the page)
P = 123638027; RP = 720836321
assert all(isprime(int(f"{P}{q*q}{RP}")) for q in (2, 3, 5, 7, 11, 13))
print("all checks passed: at most 168 terms")
