"""The two-cell wheel from a start of norm 2: (a + sqrt c)^n = a - sqrt c (mod n), with c = a^2 - 2 or a^2 + 2 and (c/n) = -1.
The norm of the start is +2 or -2, so a composite that comes home has 2^n = 2 (mod n): it is on the list of base-2 pseudoprimes."""
import sys
from math import isqrt


def jacobi(a, n):
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def pairs():
    a = 1
    while True:
        if a >= 2:
            yield a, a * a - 2
        yield a, a * a + 2
        a += 1


def start(n, amax=10**6):
    """First (a, c) with (c/n) = -1; None for a perfect square; ('factor', g) when a small c shares a factor with n."""
    r = isqrt(n)
    if r * r == n:
        return None
    for a, c in pairs():
        j = jacobi(c, n)
        if j == -1:
            return a, c
        if j == 0 and c % n:
            return 'factor', c
        if a > amax:
            return None


def power(a, c, e, n):
    x, y = 1, 0
    A, B = a % n, 1
    while e:
        if e & 1:
            x, y = (x * A + c * y * B) % n, (x * B + y * A) % n
        A, B = (A * A + c * B * B) % n, 2 * A * B % n
        e >>= 1
    return x, y


def home(n):
    """True, False, or the reason why no wheel exists."""
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    s = start(n)
    if s is None:
        return False
    if s[0] == 'factor':
        from math import gcd
        return False if 1 < gcd(s[1], n) < n else None
    a, c = s
    x, y = power(a, c, n, n)
    return x == a % n and y == n - 1


if __name__ == '__main__':
    from sympy import isprime, primerange
    top = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    wrong = und = 0
    amax = (0, 0)
    for n in range(3, top, 2):
        h = home(n)
        if h is None:
            und += 1
            continue
        if h != isprime(n):
            wrong += 1
        s = start(n)
        if s and s[0] != 'factor' and s[0] > amax[0]:
            amax = (s[0], n)
    print('odd n from 3 to %d: wrong answers %d, undecided %d, largest a used %d (for n = %d)' % (top, wrong, und, amax[0], amax[1]))
