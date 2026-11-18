"""Beyond the list: numbers above 2^64 that fool the doubling wheel, built on purpose (the two families of beyond64.py,
with larger members too), on the norm-2 wheel."""
import random
import gmpy2
from norm2 import home

B = 2**64
lines = []
F1 = []
k = 240000
while len(F1) < 4000 and k < 6000000:
    k += 1
    if gmpy2.is_prime(6 * k + 1) and gmpy2.is_prime(12 * k + 1) and gmpy2.is_prime(18 * k + 1):
        n = (6 * k + 1) * (12 * k + 1) * (18 * k + 1)
        if n > B:
            F1.append(n)
random.seed(64)
F2 = []
while len(F2) < 6000:
    p = int(gmpy2.next_prime(random.randrange(3 * 10**9, 10**12)))
    for kk in (2, 3, 4, 5, 6):
        q = kk * (p - 1) + 1
        if gmpy2.is_prime(q):
            n = p * q
            if n > B and pow(2, n - 1, n) == 1:
                F2.append(n)
F3 = []
while len(F3) < 2000:
    p = int(gmpy2.next_prime(random.randrange(10**19, 10**30)))
    for kk in (2, 3, 4, 5, 6):
        q = kk * (p - 1) + 1
        if gmpy2.is_prime(q):
            n = p * q
            if pow(2, n - 1, n) == 1:
                F3.append(n)
for name, F in (("Chernick's Carmichael numbers (6k+1)(12k+1)(18k+1)", F1), ('base-2 pseudoprimes p(k(p - 1) + 1), k = 2..6', F2), ('the same family with p of 20 to 30 digits', F3)):
    assert all(pow(2, n - 1, n) == 1 for n in F)
    h = sum(1 for n in F if home(n))
    s = '%s: %d numbers from %.2e to %.2e, all fool the doubling wheel; home on the norm-2 wheel: %d' % (name, len(F), min(F), max(F), h)
    print(s)
    lines.append(s)
open('norm2_beyond.txt', 'w').write('\n'.join(lines) + '\n')
