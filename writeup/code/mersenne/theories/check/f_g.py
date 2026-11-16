# F: primes q whose order of 2 is an odd prime l  <=>  l | q-1, l odd prime, 2^l = 1 mod q
#    factor q-1 with a smallest-prime-factor sieve; check pow(2, l, q) for each odd prime l | q-1
# G: p = 1 mod 4, p <= 166666, 6p+1 prime; 6p+1 | 2^p - 1 ; 6p+1 = x^2 + 27 y^2 (brute force)
import numpy as np
from math import isqrt

N = 10**6
spf = np.zeros(N + 1, dtype=np.int32)
for i in range(2, N + 1):
    pass
# build spf with numpy
spf = np.arange(N + 1, dtype=np.int32)
for i in range(2, isqrt(N) + 1):
    if spf[i] == i:
        sl = spf[i * i::i]
        m = sl == np.arange(i * i, N + 1, i, dtype=np.int32)
        sl[m] = i
        spf[i * i::i] = sl
isp = (spf == np.arange(N + 1)) & (np.arange(N + 1) >= 2)
primes = np.nonzero(isp)[0]


def odd_prime_factors(n):
    s = set()
    while n > 1:
        r = int(spf[n])
        s.add(r)
        while n % r == 0:
            n //= r
    s.discard(2)
    return s


marks = [1000, 10**4, 10**5, 10**6]
count = {m: 0 for m in marks}
mod3_1 = 0
model = 0.0
for q in primes:
    q = int(q)
    if q == 2:
        continue
    ls = odd_prime_factors(q - 1)
    hit = False
    for l in ls:
        if pow(2, l, q) == 1:
            hit = True
    if hit:
        for m in marks:
            if q <= m:
                count[m] += 1
        if q % 3 == 1:
            mod3_1 += 1
    if q % 8 in (1, 7):
        for l in ls:
            model += 2 * (l - 1) / (q - 1)
print("F counts", count)
print("F of those <=1e6 with q=1 mod 3:", mod3_1)
print("F model sum", round(model, 2))

# G
tot = 0
div = set()
rep = set()
for p in range(5, 166667):
    if p % 4 != 1 or not isp[p]:
        continue
    q = 6 * p + 1
    if q > N:
        from sympy import isprime
        ok = isprime(q)
    else:
        ok = bool(isp[q])
    if not ok:
        continue
    tot += 1
    if pow(2, p, q) == 1:
        div.add(p)
    y = 1
    while 27 * y * y < q:
        r = q - 27 * y * y
        x = isqrt(r)
        if x * x == r:
            rep.add(p)
            break
        y += 1
print("G total", tot, "divides", len(div), "x2+27y2", len(rep), "same set", div == rep)
