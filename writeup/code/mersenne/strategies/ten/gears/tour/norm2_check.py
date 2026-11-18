"""Checks of the norm-2 wheel in plain Python: every base-2 pseudoprime below 10^10, and the norm lemma on wheels that the number fits."""
import sys
from math import gcd
from sympy import isprime
from norm2 import home, start, power, jacobi, pairs

lines = []


def say(s):
    print(s, flush=True)
    lines.append(s)


ps = [int(x) for x in open('psp2_to1e10.txt').read().split()]
assert all(pow(2, n - 1, n) == 1 and not isprime(n) for n in ps[:2000])
came = [n for n in ps if home(n)]
amax = max((start(n)[0], n) for n in ps if start(n) and start(n)[0] != 'factor')
sq = [n for n in ps if start(n) is None]
say('base-2 pseudoprimes below 10^10: %d; came home on the norm-2 wheel: %d %s; perfect squares (no wheel): %s; largest a used: %d (n = %d)' % (len(ps), len(came), came[:5], sq, amax[0], amax[1]))

TOP = int(sys.argv[1]) if len(sys.argv) > 1 else 2000000
fit_home = 0
fit_not_psp = 0
first = []
for n in range(9, TOP, 2):
    if isprime(n):
        continue
    for a, c in pairs():
        if a > 6:
            break
        if jacobi(c, n) != 1 or gcd(2 * a * c, n) != 1:
            continue
        x, y = power(a, c, n, n)
        if x == a % n and (y == n - 1 or y == 1):
            fit_home += 1
            if len(first) < 6:
                first.append((n, a, c))
            if pow(2, n, n) != 2:
                fit_not_psp += 1
say('the norm lemma where it is not empty: odd composites below %d on wheels they fit ((c/n) = +1, a up to 6), counted when (a + sqrt c)^n = a - sqrt c or a + sqrt c: %d cases, first %s; of these not a base-2 pseudoprime: %d' % (TOP, fit_home, first, fit_not_psp))
idn = 0
for a in range(1, 31):
    for c, sg in ((a * a - 2, 1), (a * a + 2, -1)):
        if c < 2:
            continue
        e0, e1 = a * a - sg, a
        assert (a * a + c, 2 * a) == (2 * e0, 2 * e1) and e0 * e0 - c * e1 * e1 == 1 and a * a - c == 2 * sg
        idn += 1
say('the square of the start is twice a unit: (a + sqrt c)^2 = 2 ((a^2 -+ 1) + a sqrt c), and that unit has norm 1: exact for all %d starts with a up to 30' % idn)
ok = bad = 0
import random
from sympy import randprime
rng = random.Random(3)
for _ in range(300):
    n = randprime(10**6, 10**18)
    st = start(n)
    a, c = st
    sg = 1 if a * a - c == 2 else -1
    e = (a * a - sg, a)
    x, y = 1, 0
    A, B = e[0] % n, e[1] % n
    k = (n + 1) // 2
    while k:
        if k & 1:
            x, y = (x * A + c * y * B) % n, (x * B + y * A) % n
        A, B = (A * A + c * B * B) % n, 2 * A * B % n
        k >>= 1
    t = pow(2, (n - 1) // 2, n)
    if (t * x) % n == sg % n and (t * y) % n == 0:
        ok += 1
    else:
        bad += 1
say('for 300 random primes up to 10^18: 2^((n-1)/2) times the unit to the power (n+1)/2 equals the sign of the norm in %d, not in %d' % (ok, bad))
open('norm2_check.txt', 'w').write('\n'.join(lines) + '\n')
