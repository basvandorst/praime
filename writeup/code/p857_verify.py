# Puzzle 857: a cascade of 24 primes 4k+1 under f(p) = 2ab + 1, p = a^2 + b^2, a < b (writeup/857.md).
# Run: python3 p857_verify.py   (a few seconds; needs cypari2 for primality proofs and factoring)
from itertools import combinations
from math import isqrt, gcd
import cypari2

pari = cypari2.Pari()
pari.allocatemem(10**9)
proven = lambda n: int(pari('isprime(%d)' % n)) == 1          # PARI isprime: a primality proof


def two_squares(p):
    """p prime, p = 1 (mod 4): x^2 = -1 mod p, then the Euclidean algorithm until the remainder is below sqrt p."""
    x = int(pari('lift(sqrt(Mod(-1, %d)))' % p))
    a, b, r = p, x, isqrt(p)
    while b > r:
        a, b = b, a % b
    c = isqrt(p - b * b)
    assert b * b + c * c == p
    return min(b, c), max(b, c)


def f(p):
    if p % 4 != 1 or not proven(p):
        return 0
    a, b = two_squares(p)
    return 2 * a * b + 1


# 1. The cascade from our 86-digit start P: 24 distinct proven primes, then a composite.
P = 36138098107659409165026482176946981669632220236460347529493984619966279471116709053877
T = 4065497829971192852885728611545617740636798200430307756885727107124653    # top of the page's 23-cascade
chain = [P]
while True:
    q = f(chain[-1])
    if q == 0 or q in chain or not proven(q):
        end = q
        break
    chain.append(q)
assert len(chain) == 24 and len(set(chain)) == 24 and chain[1] == T
assert all(proven(x) and x % 4 == 1 for x in chain)
assert end == 166381 and not proven(166381)              # 180541 -> 166381 = 379 * 439, composite: the chain ends
print('cascade of %d proven primes; start has %d digits; second term is the page\'s 23-cascade top' % (len(chain), len(str(P))))


# 2. Predecessors: f(p) = T with p = a^2 + b^2 means ab = (T - 1)/2, gcd(a, b) = 1, a, b of opposite parity, a < b.
#    So the predecessors come from the unitary splits of (T - 1)/2 (each prime power wholly in a or in b).
def predecessors(t):
    m = (t - 1) // 2
    fac = pari('factor(%d)' % m)
    pp = [int(fac[0][i]) ** int(fac[1][i]) for i in range(len(fac[0]))]
    assert all(proven(int(fac[0][i])) for i in range(len(fac[0])))       # complete factorisation, proven primes
    out, splits = [], 0
    for k in range(len(pp) + 1):
        for S in combinations(range(len(pp)), k):
            a = 1
            for i in S:
                a *= pp[i]
            b = m // a
            if a < b and gcd(a, b) == 1 and (a + b) % 2 == 1:
                splits += 1
                p = a * a + b * b
                if proven(p):
                    out.append(p)
    return out, splits


preT, sT = predecessors(T)
assert preT == [P]
preP, sP = predecessors(P)
assert preP == []
print('T has %d admissible splits and exactly one prime predecessor (P); P has %d splits and none' % (sT, sP))
print('all checks passed')
