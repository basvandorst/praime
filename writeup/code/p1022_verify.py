# Puzzle 1022: checks every claim of writeup/1022.md for x(n) = (n-1)^n + n.
# Run from anywhere: python3 p1022_verify.py          (about 15 seconds)
#                    python3 p1022_verify.py --full   (also proves x(n) composite for all 4 <= n <= 3000, about 4 minutes)
# Needs sympy and gmpy2.  The heuristic counts are recomputed here in plain Python, independently of p1022_model.c.
import sys
from math import gcd, exp, log
from sympy import isprime, primerange, n_order, factorint, legendre_symbol, symbols, factor_list, Poly
import gmpy2

x = lambda n: (n - 1)**n + n

# 1. x(n) is odd; x(n) = (-1)^n (mod n) and x(n) = 1 (mod n - 1), so no prime factor of n(n-1) divides x(n).
for n in range(3, 2001):
    X = x(n)
    assert X % 2 == 1 and X % n == (-1)**n % n and X % (n - 1) == 1 % (n - 1)

# 2. Algebraic factor: n = 2 (mod 3)  =>  (n^2 - n + 1) | x(n), and n^2 - n + 1 = 3 * (odd number not divisible by 3).
#    For n != 2 (mod 3): gcd(x(n), n^2 - n + 1) = 1 and 3 does not divide x(n).
for n in range(2, 2001):
    X, F = x(n), n*n - n + 1
    if n % 3 == 2:
        assert X % F == 0 and F % 3 == 0 and F % 9 != 0
    else:
        assert gcd(X, F) == 1 and X % 3 != 0

# 3. Selmer (1956): Y^n + Y + 1 is irreducible for n != 2 (mod 3) and (Y^2 + Y + 1) * irreducible for n = 2 (mod 3), n > 2.
Y = symbols('Y')
for n in range(2, 61):
    degs = sorted(Poly(f, Y).degree() for f, _ in factor_list(Y**n + Y + 1)[1])
    assert degs == ([n] if n % 3 != 2 else [2] if n == 2 else [2, n - 2])

# 4. Modulo n + 1: x(n) = (-2)^n - 1, so (n+1) | x(n) iff (-2)^n = 1 (mod n+1): always when n + 1 >= 5 is prime (Fermat),
#    and for composite n + 1 exactly at the base-2 Fermat pseudoprimes (OEIS A001567).
A001567 = [341, 561, 645, 1105, 1387, 1729, 1905, 2047, 2465, 2701, 2821, 3277, 4033, 4369, 4371, 4681,
           5461, 6601, 7957, 8321, 8481, 8911, 10261, 10585, 11305, 12801]
comp = []
for n in range(2, 13390):
    m = n + 1
    v = (pow(n - 1, n, m) + n) % m
    assert v == (pow(-2, n, m) - 1) % m
    assert gcd(v, m) == gcd((pow(-2, n, m) - 1) % m, m)
    if v == 0 and not isprime(m):
        comp.append(m)
    if isprime(m) and m >= 5:
        assert v == 0
assert comp == A001567

# 5. A prime r | n + 1 divides x(n) iff ord_r(2) | n (n even) or 2^n = -1 (mod r) (n odd).  Eric Snyder's
#    "almost all" for even n != 2 (mod 3) with n + 1 composite: 1045 of 2877 up to 13389 (36.3 %).
tot = hit = 0
for n in range(4, 13390):
    for r in factorint(n + 1):
        if r == 2:
            continue
        rule = (n % n_order(2, r) == 0) if n % 2 == 0 else (pow(2, n, r) == r - 1)
        assert ((pow(n - 1, n, r) + n) % r == 0) == rule
    if n % 2 == 0 and n % 3 != 2 and not isprime(n + 1):
        tot += 1
        hit += any((pow(n - 1, n, r) + n) % r == 0 for r in factorint(n + 1))
assert (tot, hit) == (2877, 1045)

# 6. Modulo 2n - 1: 2^n x(n) = n 2^n + (-1)^n (the Cullen number C_n for even n, the Woodall number W_n for odd n),
#    so (2n-1) | x(n) iff (-2)^(n-1) = 1 (mod 2n-1); for prime 2n - 1 that is (-2 / 2n-1) = 1, i.e. 2n - 1 = 1, 3 (mod 8).
comp2 = []
for n in range(2, 13390):
    M = 2*n - 1
    X = (pow(n - 1, n, M) + n) % M
    assert pow(2, n, M) * X % M == (n * pow(2, n, M) + (-1)**n) % M
    assert (X == 0) == (pow(-2, n - 1, M) == 1)
    if isprime(M) and M > 3:
        assert (X == 0) == (M % 8 in (1, 3)) == (legendre_symbol(-2 % M, M) == 1)
    elif X == 0 and not isprime(M):
        comp2.append(M)
assert comp2 == [341, 561, 1105, 1729, 1905, 2465, 4033, 4681, 5461, 6601, 8321, 8481, 10261, 10585, 12801,
                 15709, 15841, 16705, 18705, 25761]

# 7. No other linear family: for k r = j n - i (gcd(i, j) = 1) a prime r dividing x(n) divides
#    D = (i-j)^e j^j - (-i)^j j^e (e = k + i >= 0) or j^(j-e) - (-i)^j (i-j)^(-e) (e < 0), and D = 0 only for
#    (j, i, k) = (1, -1, 1) [r = n + 1] and (2, 1, 1) [r = 2n - 1].
def D(j, i, k):
    e = k + i
    return (i - j)**e * j**j - (-i)**j * j**e if e >= 0 else j**(j - e) - (-i)**j * (i - j)**(-e)
zero = []
for j in range(1, 7):
    for i in range(-40, 41):
        if gcd(i, j) != 1:
            continue
        for k in range(1, 13):
            d = D(j, i, k)
            if d == 0:
                zero.append((j, i, k))
            for n in range(2, 1501):
                t = j*n - i
                if t > 0 and t % k == 0 and isprime(t // k):
                    r = t // k
                    if (pow(n - 1, n, r) + n) % r == 0 and (j*i*(i - j)) % r:
                        assert d == 0 or d % r == 0
assert zero == [(1, -1, 1), (2, 1, 1)]

# 8. Perfect-power n: n = 4 b^4 gives a^4 + 4 b^4 (Sophie Germain), n = c^q with an odd prime q | c gives A^q + c^q.
a, b = 63**16, 2
assert x(4) == 85 == (9 - 6 + 2) * (9 + 6 + 2) and x(64) % (a*a - 2*a*b + 2*b*b) == 0
assert x(27) % (26**9 + 3) == 0 and x(216) % (215**72 + 6) == 0 and x(729) % (728**243 + 9) == 0

# 9. No covering set: n = 1 (mod r) gives x(n) = 0^n + 1 = 1 (mod r), and n = 0 (mod r) gives (-1)^n.
P = 1
for r in primerange(2, 200):
    P *= r
for t in range(1, 6):
    n = 1 + t*P
    assert (pow(n - 1, n, P) + n) % P == 1

# 10. Every prime factor r of x(n) has (-n / r) = 1 for even n (x = A^2 + n) and (-n(n-1) / r) = 1 for odd n
#     (x = (n-1) A^2 + n), with A = (n-1)^floor(n/2).
for n in range(4, 1200):
    Dn = -n if n % 2 == 0 else -n*(n - 1)
    for r in primerange(3, 5000):
        if (pow(n - 1, n, r) + n) % r == 0:
            assert legendre_symbol(Dn % r, r) == 1

# 11. The mod 4 law: 5 never divides x(n) for n = 2, 3 (mod 4); for n = 3 (mod 4) the factor n + 1 is impossible
#     (n + 1 even, x odd) and a prime 2n - 1 = 5 (mod 8) never divides.
for n in range(2, 4001):
    if n % 4 in (2, 3):
        assert x(n) % 5 != 0
    if n % 4 == 3 and isprime(2*n - 1):
        assert x(n) % (2*n - 1) != 0

# 12. Of the 13386 values 4 <= n <= 13389, the three proved divisors explain 7298 (54.5 %).
expl = sum(1 for n in range(4, 13390) if n % 3 == 2 or (pow(n - 1, n, n + 1) + n) % (n + 1) == 0
           or (pow(n - 1, n, 2*n - 1) + n) % (2*n - 1) == 0)
assert expl == 7298

# 13. Heuristic count (exact sieve to R, the two divisor families, then Mertens: P(prime) = e^gamma ln R / ln x).
def expected(fam, lo, hi, R=10000):
    ps = list(primerange(2, R + 1))
    eg = exp(0.5772156649015329) * log(R)
    E = 0.0
    surv = 0
    for n in range(lo, hi + 1):
        if fam == 0 and n % 3 == 2:
            continue
        dead = False
        for r in ps:
            e = (n - 1) % (r - 1) + 1
            v = (pow(n - 1, e, r) + n) if fam == 0 else (3*pow(n, e, r) + 2) if fam == 1 else ((n - 1)*pow(n, e, r) + 1)
            if v % r == 0:
                dead = True
                break
        if not dead and fam == 0:
            dead = gcd((pow(n - 1, n, n + 1) + n) % (n + 1), n + 1) > 1 or \
                   gcd((pow(n - 1, n, 2*n - 1) + n) % (2*n - 1), 2*n - 1) > 1
        if not dead:
            surv += 1
            lx = n*log(n - 1) if fam == 0 else n*log(n) + log(3) if fam == 1 else n*log(n) + log(n - 1)
            E += eg / lx
    return E, surv
E0, s0 = expected(0, 30, 13389)
E1, _ = expected(1, 30, 45000)
E2, _ = expected(2, 30, 2835)
print("x(n), 30 <= n <= 13389: %d survivors, expected primes %.3f, chance of none %.0f %%" % (s0, E0, 100*exp(-E0)))
print("calibration 3n^n+2 (A160360), 30 <= n <= 45000: expected %.2f, observed 2 (93, 605)" % E1)
print("calibration (n-1)n^n+1 (A108879), 30 <= n <= 2835: expected %.2f, observed 4 (93, 100, 1536, 2835)" % E2)
assert s0 == 1499 and abs(E0 - 1.528) < 0.001 and abs(E1 - 1.591) < 0.001 and abs(E2 - 2.540) < 0.001

# 14. x(n) is composite for 4 <= n <= LIM: every n is either in a class above or has a Miller-Rabin witness
#     (a witness proves compositeness).  Default LIM = 1000; --full gives LIM = 3000 (Martin Ruiz's range).
LIM = 3000 if '--full' in sys.argv else 1000
for n in range(4, LIM + 1):
    assert not gmpy2.is_prime(gmpy2.mpz(n - 1)**n + n, 1), n
print("x(n) composite for 4 <= n <= %d (Miller-Rabin witnesses)" % LIM)
print("all checks passed")
