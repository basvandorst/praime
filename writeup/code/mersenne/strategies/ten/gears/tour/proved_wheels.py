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
