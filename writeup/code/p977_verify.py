# Puzzle 977 Q1: primes q, r, s with 4q, 3r, 2s consecutive integers, q of 566 digits (writeup/977.md).
# Run: python3 p977_verify.py          (under a second: LLR proofs; needs gmpy2 and sympy)
#      python3 p977_verify.py --pari   (also APRCL proofs with PARI through cypari2, a few minutes)
import sys
import gmpy2
from sympy import isprime


def lucas_V(h, P, N):
    """V_h(P, 1) mod N: V_(2m) = V_m^2 - 2, V_(2m+1) = V_m V_(m+1) - P."""
    a, b = gmpy2.mpz(2), gmpy2.mpz(P) % N
    for bit in bin(h)[2:]:
        if bit == '1':
            a, b = (a * b - P) % N, (b * b - 2) % N
        else:
            a, b = (a * a - 2) % N, (a * b - P) % N
    return a


def llr_prove(h, n):
    """Lucas-Lehmer-Riesel with Rodseth's start: for N = h*2^n - 1 (h odd, h < 2^n) returns True iff N is prime
    (None if no starting value P <= 2000 exists, which happens only for composite N)."""
    assert h % 2 == 1 and (1 << n) > h and n >= 3
    N = gmpy2.mpz(h) * gmpy2.mpz(2) ** n - 1
    P = 3
    while not (gmpy2.jacobi(P - 2, N) == 1 and gmpy2.jacobi(P + 2, N) == -1):
        P += 1
        if P > 2000:
            return None
    u = lucas_V(h, P, N)
    for _ in range(n - 2):
        u = (u * u - 2) % N
    return u == 0


# 1. The test itself, against sympy on all h*2^n - 1 with odd h < 400, 3 <= n < 40, h < 2^n.
for h in range(1, 400, 2):
    for n in range(3, 40):
        if (1 << n) > h:
            r = llr_prove(h, n)
            N = h * 2**n - 1
            assert r == isprime(N) or (r is None and not isprime(N))

# 2. The trios: m = k*2^1850, (q, r, s) = (3m - 1, 4m - 1, 6m - 1), so 4q, 3r, 2s = 12m - 4, 12m - 3, 12m - 2.
n = 1850
for k in (152439497, 139751997):
    m = k * 2**n
    q, r, s = 3 * m - 1, 4 * m - 1, 6 * m - 1
    assert 3 * r == 4 * q + 1 and 2 * s == 3 * r + 1
    assert llr_prove(3 * k, n) and llr_prove(k, n + 2) and llr_prove(3 * k, n + 1)
    assert all(gmpy2.is_prime(x, 25) for x in (q, r, s))
    print('k = %d: q, r, s proven prime (LLR); digits %d, %d, %d' % (k, len(str(q)), len(str(r)), len(str(s))))
    if '--pari' in sys.argv:
        import cypari2
        pari = cypari2.Pari()
        pari.allocatemem(2 * 10**9)
        assert all(int(pari('isprime(%d, 2)' % x)) == 1 for x in (q, r, s))
        print('   also proven by APRCL (PARI)')
print('all checks passed')
