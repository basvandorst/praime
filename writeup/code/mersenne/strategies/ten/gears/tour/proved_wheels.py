"""Wheels that are proved, for comparison with the fast ones.
1. Wilson's wheel: one cell, multiplied by 1, 2, 3, ... turn after turn; after n - 1 turns it shows -1 exactly when n is prime (Wilson's theorem).
2. The test of Agrawal, Kayal and Saxena (2002), in the version of Lenstra's bound used in textbooks: the ring of r cells with
   offsets a = 1 .. floor(sqrt(phi(r)) log2 n), r the smallest number for which n has order above (log2 n)^2 mod r. Proved for every n.
Both are run on every n up to a small bound and compared with sympy; then the size of the proved machine is worked out for 2^4423 - 1."""
import math, sys
from sympy import isprime, totient, n_order, integer_nthroot
from ref import ring_test
W = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
A = int(sys.argv[2]) if len(sys.argv) > 2 else 400
bad = 0; cell = 1
for n in range(2, W + 1):
    f = 1
    for k in range(2, n): f = f*k % n
    if (f == n - 1) != isprime(n): bad += 1
print("Wilson's wheel, n = 2 ..", W, ': wrong answers', bad)
def aks_params(n):
    L = math.log2(n); r = 2
    while True:
        if math.gcd(r, n) == 1 and n_order(n, r) > L*L: break
        r += 1
    return r, int(math.isqrt(int(totient(r)))*L) if False else int(math.floor(math.sqrt(int(totient(r)))*L))
def aks(n):
    if n < 2: return False
    for b in range(2, n.bit_length() + 1):
        root, exact = integer_nthroot(n, b)
        if exact: return False
    r, amax = aks_params(n)
    for a in range(2, min(r, n - 1) + 1):
        if 1 < math.gcd(a, n) < n: return False
    if n <= r: return True
    for a in range(1, amax + 1):
        if not ring_test(n, r, a): return False
    return True
bad = 0; biggest = (0, 0, 0); ran = 0
for n in range(2, A + 1):
    v = aks(n)
    if v != isprime(n): bad += 1
    if isprime(n):
        r, amax = aks_params(n)
        if n > r: ran += 1; biggest = max(biggest, (r, amax, n))
print('the proved ring test, n = 2 ..', A, ': wrong answers', bad, '| primes that needed the rings:', ran, '| largest machine: ring of %d cells with %d offsets for n = %d' % biggest)
p = 4423; L = p
print('for n = 2^4423 - 1 (1,332 digits): the ring needs more than (log2 n)^2 = %s cells and about sqrt(r) log2 n = %s offsets; machine 10 uses 5 cells and one offset' % (format(p*p, ','), format(int(math.sqrt(p*p)*p), ',')))
from sympy import factorint
def lucas(n, amax=250):
    """Lucas (1876): n is prime if some a has a^(n-1) = 1 and a^((n-1)/q) != 1 (mod n) for every prime q dividing n - 1.
    Returns (proof found, step used, certainly composite)."""
    if n == 2: return True, 0, False
    if n % 2 == 0: return False, 0, True
    qs = list(factorint(n - 1))
    for a in range(2, min(amax, n - 1) + 1):
        if pow(a, n - 1, n) != 1: return False, a, True
        if all(pow(a, (n - 1)//q, n) != 1 for q in qs): return True, a, False
    return False, 0, False
LW = int(sys.argv[3]) if len(sys.argv) > 3 else 200000
bad = und = 0; amax = (0, 0)
for n in range(2, LW + 1):
    ok, a, comp = lucas(n)
    if ok != isprime(n): bad += 1
    if not ok and not comp: und += 1
    if ok: amax = max(amax, (a, n))
print("Lucas's wheel, n = 2 ..", LW, ': wrong answers', bad, '| not decided', und, '| largest step needed: %d (for n = %d)' % amax)
P = 3*2**2208 + 1; ok, a, comp = lucas(P)
print('3 * 2^2208 + 1 (%d digits): proof found %s with step %d' % (len(str(P)), ok, a))
from sympy import jacobi_symbol, primerange
from math import gcd
def dpow(a, c, e, n):
    x, y, A, B = 1 % n, 0, a % n, 1 % n
    while e:
        if e & 1: x, y = (x*A + c*y*B) % n, (x*B + y*A) % n
        A, B = (A*A + c*B*B) % n, 2*A*B % n
        e >>= 1
    return x, y
def nplus(n, amax=40):
    """Two cells with check turns (the N + 1 test of Lucas and Lehmer): z = a + sqrt(c), (c/n) = -1.
    n is prime if z^n = a - sqrt(c) and, for every prime q dividing n + 1, z^(2(n+1)/q) - N^((n+1)/q) has a norm
    without a factor in common with n (N = a^2 - c). Returns (proof found, c, a, certainly composite)."""
    if n == 2: return True, 0, 0, False
    if n % 2 == 0: return False, 0, 0, True
    c = 0
    for q in primerange(3, 200):
        if n % q == 0:
            if n == q: continue
            return False, 0, 0, True
        if jacobi_symbol(q, n) == -1: c = q; break
    if not c: return False, 0, 0, True
    qs = list(factorint(n + 1))
    for a in range(1, amax + 1):
        N = (a*a - c) % n
        if N == 0: continue
        if gcd(N, n) != 1: return False, c, a, True
        h = dpow(a, c, n, n)
        if h != (a % n, n - 1): return False, c, a, True
        good = True
        for q in qs:
            k = (n + 1)//q; w = dpow(a, c, 2*k, n); d0 = (w[0] - pow(N, k, n)) % n
            if gcd((d0*d0 - c*w[1]*w[1]) % n, n) != 1: good = False; break
        if good: return True, c, a, False
    return False, c, 0, False
NW = int(sys.argv[4]) if len(sys.argv) > 4 else 100000
bad = und = 0; amax = (0, 0)
for n in range(2, NW + 1):
    ok, c, a, comp = nplus(n)
    if ok != isprime(n): bad += 1
    if not ok and not comp: und += 1
    if ok: amax = max(amax, (a, n))
print('two cells with check turns, n = 2 ..', NW, ': wrong answers', bad, '| not decided', und, '| largest start needed: %d (for n = %d)' % amax)
for p in (521, 4423):
    ok, c, a, comp = nplus(2**p - 1)
    print('2^%d - 1: proof found %s with step %d and start %d' % (p, ok, c, a))
ok, c, a, comp = nplus(2**523 - 1); print('2^523 - 1: proof found %s, certainly not prime %s' % (ok, comp))
