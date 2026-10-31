# Puzzle 61 Q7: last digit of N = a^p - a + 1 when N also equals b^2 - b + 1.
from sympy import isprime, primerange
# (1) the residue facts
assert {b * (b - 1) % 5 for b in range(5)} == {0, 1, 2}                          # square side: N mod 5 in {1,2,3}
for p in primerange(3, 200):
    vals = {(pow(a, p, 5) - a) % 5 for a in range(5)}
    assert vals == ({0} if p % 4 == 1 or p == 5 else {0, 1, 4}), p              # odd-p side
# last digit allowed: N odd, N mod 5 in {1,2,3} (square side) and in {1,2,0} (p = 3 mod 4) or {1} (p = 1 mod 4)
allowed = {d for d in range(10) if d % 2 == 1 and d % 5 in {1, 2}}
assert allowed == {1, 7}
# ending in 7 needs b = 3 mod 5, p = 3 mod 4 and a = 2 mod 5
assert {b % 5 for b in range(5) if (b * (b - 1) + 1) % 5 == 2} == {3}
assert {a for a in range(5) if (a**3 - a + 1) % 5 == 2} == {2}
# (2) the page's examples
pairs = [(7, (2, 3), (3, 2)), (31, (2, 5), (6, 2)), (211, (6, 3), (15, 2)), (241, (3, 5), (16, 2)),
         (8191, (2, 13), (91, 2)), (78121, (5, 7), (280, 2))]
for N, (a, p), (b, q) in pairs:
    assert a**p - a + 1 == N == b**q - b + 1 and isprime(N)
    reason = "p = 1 mod 4 or 5" if (p % 4 == 1 or p == 5) else f"p = 3 mod 4, a = {a % 5} mod 5"
    print(N, "ends in", N % 10, "|", reason, "| b mod 5 =", b % 5)
print("Q7 verified: a prime with a square expression ends in 1 or 7; 7 needs p = 3 (mod 4), a = 2 and b = 3 (mod 5)")
