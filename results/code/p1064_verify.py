# Puzzle 1064: checks every claim in results/1064.md (about 1 minute).
# G(p) = gcd(2^p + 1, 3^p + 1). Theorem: for an odd prime p, a prime q divides G(p) iff q = 1 + d p with -2 and -3 both
# d-th power residues mod q; then 6 | d, d != 4 (mod 8), and d = 2 (mod 4) forces p = d/2 (mod 4); so q = 1 or 19 (mod 24).
import gmpy2
from sympy import primerange, factorint, isprime

def G(p): return gmpy2.gcd(gmpy2.mpz(2)**p + 1, gmpy2.mpz(3)**p + 1)

# 1. every prime factor of G(p), all odd primes p < 10^5, exactly (full gcd, full factorisation)
nf = 0
for p in primerange(3, 10**5):
    for q in factorint(int(G(p))):
        nf += 1; d, r = divmod(q - 1, p)
        assert r == 0 and d % 6 == 0 and d % 8 != 4 and q % 24 in (1, 19)
        if d % 4 == 2: assert p % 4 == (d // 2) % 4
        assert pow(-2 % q, (q - 1) // d, q) == 1 and pow(q - 3, (q - 1) // d, q) == 1      # -2, -3 are d-th powers
assert nf == 156, nf                                                                       # 156 prime factors in all
assert G(2) == 5 and G(3) == 1                                                             # p = 2 is the one exception (5 != 1 mod 12)

# 2. the converse on a range: q = 1 + d p prime with -2, -3 d-th powers  =>  q | G(p)   (p < 2000, d <= 600)
for p in primerange(5, 2000):
    g = G(p)
    for d in range(2, 601, 2):
        q = 1 + d * p
        if isprime(q) and pow(q - 2, p, q) == 1 and pow(q - 3, p, q) == 1: assert g % q == 0

# 3. the composite G(p) listed on the page (Kalogeropoulos, Cavegn, Gupta, Vantieghem, van Delden, Volpatti), recomputed in full
page = {2243399: (53841577, 174985123), 2334547: (14007283, 70036411), 2743723: (16462339, 65849353),
        3932207: (23593243, 94372969), 4623107: (27738643, 249647779), 4716343: (28298059, 254682523),
        5482423: (32894539, 164472691), 5993411: (35960467, 179802331), 6490151: (38940907, 1752340771)}
for p, (a, b) in page.items():
    assert G(p) == a * b and isprime(a) and isprime(b) and isprime(p)
    da, db = (a - 1) // p, (b - 1) // p
    assert (a - 1) % p == 0 and (b - 1) % p == 0 and p % 4 == 3                           # all nine have p = 3 (mod 4)
    for d in (da, db): assert d % 6 == 0 and (d % 8 == 0 or (d % 4 == 2 and p % 4 == (d // 2) % 4))
print([((a - 1) // p, (b - 1) // p) for p, (a, b) in page.items()])
print("all claims verified")
