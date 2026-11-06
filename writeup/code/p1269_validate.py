"""Puzzle 1269: triangles with prime sides on lattice points of Z^n.

Checks, with plain integers:
  1. parity: for a triangle with three odd prime sides the integer u.v = (p^2+q^2-r^2)/2 does not exist;
     the only valid triangles with a side 2 are (2,2,2), (2,2,3) and (2,p,p);
  2. brute force for small sides: which of these embed in Z^2, Z^3, Z^4 (vectors enumerated);
  3. the three prime lists: Z^3 (p^2-1 a sum of two squares), Z^4 only (Legendre), Z^5 only (p^2-1 = 4^k (8w+7)),
     each compared with the explicit characterisation through p-1 and p+1;
  4. explicit embeddings of (2,2,2), (2,3,3), (2,5,5), (2,31,31), (2,127,127) with their distances.
"""
from sympy import primerange, factorint, isprime
from itertools import product
from math import isqrt

# 1. parity and the shape of valid triangles -------------------------------
P = list(primerange(2, 60))
tri = [(p, q, r) for p in P for q in P for r in P if p <= q <= r and p + q > r]
odd3 = [t for t in tri if t[0] > 2]
assert all((p*p + q*q - r*r) % 2 == 1 for p, q, r in odd3)          # u.v would be a half-integer
with2 = [t for t in tri if t[0] == 2]
assert all(t in [(2, 2, 2), (2, 2, 3)] or t[1] == t[2] for t in with2)
assert (2*2 + 2*2 - 3*3) % 2 == 1                                       # (2,2,3) also fails parity

# 2. brute force for small sides ------------------------------------------
def vecs(n, N2):
    R = isqrt(N2)
    return [v for v in product(range(-R, R + 1), repeat=n) if sum(x*x for x in v) == N2]

def embeds(t, n):
    p, q, r = t
    if (p*p + q*q - r*r) % 2:
        return False
    d = (p*p + q*q - r*r) // 2
    U, V = vecs(n, p*p), vecs(n, q*q)
    return any(sum(a*b for a, b in zip(u, v)) == d for u in U if u[0] >= 0 for v in V)

least = {}
for t in [(2, 2, 2), (2, 2, 3)] + [(2, p, p) for p in primerange(3, 20)]:
    least[t] = next((n for n in (2, 3, 4) if (n < 4 or t[1] <= 13) and embeds(t, n)), None)
assert least[(2, 2, 2)] == 4 and least[(2, 2, 3)] is None
assert {p: least[(2, p, p)] for p in (3, 5, 7, 11, 13, 17, 19)} == {3: 3, 5: 4, 7: 4, 11: 4, 13: 4, 17: 3, 19: 3}
print("brute force, least n:", least)

# 3. the prime lists ---------------------------------------------------------
def two_squares(m):                 # m > 0 is a sum of two squares iff every prime 3 mod 4 has even exponent
    return all(e % 2 == 0 for q, e in factorint(m).items() if q % 4 == 3)

def three_squares(m):               # Legendre: not of the form 4^k (8w+7)
    while m % 4 == 0:
        m //= 4
    return m % 8 != 7

def v2(m):
    return (m & -m).bit_length() - 1

z3, z4, z5 = [], [], []
for p in primerange(3, 2000):
    m = p*p - 1
    if two_squares(m):
        z3.append(p)
    elif three_squares(m):
        z4.append(p)
    else:
        z5.append(p)
    # characterisations through p-1 and p+1 (gcd(p-1, p+1) = 2)
    assert two_squares(m) == (two_squares(p - 1) and two_squares(p + 1))
    s = v2(m)                                       # = 1 + v2(the one of p+-1 divisible by 4)
    assert (not three_squares(m)) == (s % 2 == 0 and (m >> s) % 8 == 7)
assert z3[:12] == [3, 17, 19, 73, 163, 179, 233, 467, 521, 577, 739, 809]   # the page prints 457: 456 = 8*3*19
assert not two_squares(457*457 - 1) and 457 in z4
assert z4[:10] == [5, 7, 11, 13, 23, 29, 37, 41, 43, 47]
assert z5[:12] == [31, 89, 103, 127, 167, 281, 359, 409, 487, 601, 743, 857]
print("Z^3:", z3[:12]); print("Z^4 only:", z4[:10]); print("Z^5 only:", z5[:12])

# brute force agrees with the lists for the small primes checked above
for p in (3, 17, 19):
    assert p in z3
for p in (5, 7, 11, 13):
    assert p in z4

# counts up to 10^6
c3 = c5 = tot = 0
for p in primerange(3, 10**6):
    tot += 1
    m = p*p - 1
    if two_squares(m):
        c3 += 1
    elif not three_squares(m):
        c5 += 1
print(f"odd primes below 10^6: {tot}; Z^3: {c3}; Z^4 only: {tot - c3 - c5}; Z^5 only: {c5}")

# 4. explicit embeddings -----------------------------------------------------
def dist2(a, b):
    return sum((x - y)**2 for x, y in zip(a, b))

examples = {
    (2, 2, 2): [(0, 0, 0, 0), (1, 1, 1, 1), (1, 1, 1, -1)],
    (2, 3, 3): [(0, 0, 0), (2, 0, 0), (1, 2, 2)],
    (2, 5, 5): [(0, 0, 0, 0), (2, 0, 0, 0), (1, 4, 2, 2)],
    (2, 31, 31): [(0,)*5, (2, 0, 0, 0, 0), (1, 28, 12, 4, 4)],
    (2, 127, 127): [(0,)*5, (2, 0, 0, 0, 0), (1, 112, 48, 32, 16)],
}
for (p, q, r), (A, B, C) in examples.items():
    assert sorted([dist2(A, B), dist2(A, C), dist2(B, C)]) == sorted([p*p, q*q, r*r])
    assert all(isprime(x) for x in (p, q, r))
# 960 and 16128 are not sums of three squares (so (2,31,31) and (2,127,127) need Z^5)
for m in (960, 16128):
    assert not any(a*a + b*b + c*c == m for a in range(isqrt(m) + 1) for b in range(a + 1) for c in range(b + 1))
    assert not three_squares(m)
print("explicit embeddings checked; 960 = 4^3*15 and 16128 = 4^4*63 are not sums of three squares")
