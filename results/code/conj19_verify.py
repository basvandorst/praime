# Conjecture 19: an extreme Carmichael number C3 = p q r with q = p^2 + p - 1, r = (p^3 + p^2 - p + 1)/2 (Pinch-Nash bound).
# p = 25971283 * 2^560 + 1.  Proves p (Proth), q and r (Pocklington / Brillhart-Lehmer-Selfridge with N - 1 partly factored),
# checks Korselt's criterion, and optionally (--pari) re-proves q and r with PARI's APR-CL.
import sys, gmpy2
from math import gcd, isqrt
from sympy import factorint
k, n = 25971283, 560
p = k * 2**n + 1
q = p * p + p - 1
r = (p**3 + p * p - p + 1) // 2
assert (p**3 + p * p - p + 1) % 2 == 0 and q - 1 == (p - 1) * (p + 2) and 2 * (r - 1) == (p - 1) * (p + 1)**2

def proth(N, k, n):                                   # N = k 2^n + 1, k < 2^n: N prime iff a^((N-1)/2) = -1 for some a
    assert k < 2**n
    for a in range(3, 1000, 2):
        if gmpy2.jacobi(a, N) == -1: return pow(a, (N - 1) // 2, N) == N - 1
    return False

def pock_bls(N, F_factors):
    # F_factors: dict prime -> exponent with F | N - 1. Pocklington conditions for every prime of F, then
    # F > sqrt(N): prime;  N^(1/3) <= F <= sqrt(N): BLS Theorem 5, N = c2 F^2 + c1 F + 1, prime iff c1^2 - 4 c2 is not a square.
    F = 1
    for f, e in F_factors.items(): F *= f**e
    assert (N - 1) % F == 0
    for f in F_factors:
        for a in range(2, 200):
            if pow(a, N - 1, N) != 1: return False
            if gcd(pow(a, (N - 1) // f, N) - 1, N) == 1: break
        else: return False
    if F * F > N: return True
    assert F**3 >= N
    c2, rem = divmod(N - 1, F * F); c1 = rem // F
    assert c2 * F * F + c1 * F + 1 == N and 0 <= c1 < F
    d = c1 * c1 - 4 * c2
    return d < 0 or isqrt(d)**2 != d

kf = factorint(k)
assert proth(p, k, n)
# q - 1 = k 2^n (p + 2): add the small prime factors of p + 2 found by trial division
small = {}; m = p + 2
for f in range(3, 10**6, 2):
    if m % f == 0 and gmpy2.is_prime(f):
        e = 0
        while m % f == 0: m //= f; e += 1
        small[f] = e
Fq = dict(kf); Fq[2] = Fq.get(2, 0) + n
for f, e in small.items(): Fq[f] = Fq.get(f, 0) + e
assert pock_bls(q, Fq)
# r - 1 = k 2^(n+1) ((p+1)/2)^2
Fr = dict(kf); Fr[2] = Fr.get(2, 0) + n + 1
assert pock_bls(r, Fr)
N = p * q * r
assert all((N - 1) % (x - 1) == 0 for x in (p, q, r))           # Korselt: a Carmichael number
print("k factors:", kf, "| small factors of p+2:", small)
print("digits: p", len(str(p)), "q", len(str(q)), "r", len(str(r)), "C3", len(str(N)))
if "--pari" in sys.argv:
    import cypari2; pari = cypari2.Pari(); pari.allocatemem(2 * 10**9)
    assert pari.isprime(p) and pari.isprime(q) and pari.isprime(r)
    print("PARI isprime: p, q, r proven")
print("all claims verified")
