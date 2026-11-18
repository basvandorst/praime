"""The two-cell wheel on Mersenne numbers n = 2^p - 1 (p an odd prime) is the Lucas-Lehmer test.
Step: n = 1 (mod 3) and n = 3 (mod 4), so (3/n) = -1 and the wheel takes c = 3.
Home means (1 + s)^n = 1 - s with s^2 = 3, so (1 + s)^(n + 1) = -2. Since (1 + s)^2 = 2(2 + s) and n + 1 = 2^p,
(1 + s)^(n + 1) = 2^(2^(p - 1)) (2 + s)^(2^(p - 1)), and 2^(2^(p - 1)) = 2 (mod n) because 2^p = 1 (mod n) and 2^(p - 1) = 1 (mod p).
So home  <=>  (2 + s)^(2^(p - 1)) = -1 (mod n), which is the Lucas-Lehmer condition, and that is a proof of primality.
Checked here: wheel and Lucas-Lehmer agree for every odd prime p below LIM; then a small benchmark (one run each)."""
import sys, time, math
from gmpy2 import mpz, jacobi
from sympy import primerange
def fold(x, p, n):
    while x > n: x = (x & n) + (x >> p)
    return 0 if x == n else x
def ll(p):
    n = (mpz(1) << p) - 1; s = mpz(4)
    for _ in range(p - 2): s = fold(s*s - 2 + n, p, n)
    return s == 0
def step(n):
    for c in primerange(3, 200):
        if n % c == 0: return None
        if jacobi(c, n) == -1: return c
def duo(p):
    n = (mpz(1) << p) - 1; c = step(n); a, b = mpz(1), mpz(0)
    for _ in range(p):
        a, b = fold(a*a + c*b*b, p, n), fold(2*a*b, p, n)
        a, b = fold(a + c*b, p, n), fold(a + b, p, n)
    return c, (a == 1 and b == n - 1)
def ring(n):
    r = 5
    while True:
        m = int(n % r)
        if m and math.gcd(m, r) == 1 and m*m % r != 1: return r
        r += 1
def pas(p):
    n = (mpz(1) << p) - 1; r = ring(n); a = [mpz(0)]*r; a[0] = mpz(1)
    for t in range(p):
        s = [mpz(0)]*r
        for i in range(r):
            s[2*i % r] += a[i]*a[i]
            for j in range(i + 1, r): s[(i + j) % r] += 2*a[i]*a[j]
        a = [fold(v, p, n) for v in s]
        a = [fold(a[k] + a[k - 1], p, n) for k in range(r)]
    want = [mpz(0)]*r; want[0] += 1; want[int(n % r)] += 1
    return r, a == want
LIM = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
LP = int(sys.argv[2]) if len(sys.argv) > 2 else 700
ps = list(primerange(3, LIM)); bad = 0; home = []; steps = set(); badp = 0; np_ = 0
for p in ps:
    x = ll(p); c, y = duo(p); steps.add(c)
    if x != y: bad += 1
    if y: home.append(p)
    if p < LP:
        np_ += 1; r, z = pas(p)
        if z != x: badp += 1
print('odd prime exponents p below', LIM, ':', len(ps), '| steps used by the two-cell wheel:', sorted(steps), '| two-cell wheel differs from Lucas-Lehmer:', bad)
print('exponents whose wheel came home:', home)
print("Pascal's wheel, the", np_, 'exponents below', LP, ': differs from Lucas-Lehmer:', badp)
print('benchmark, seconds, one run (python with gmpy2):')
for p in [int(x) for x in sys.argv[3:]] or [4423, 9689, 11213]:
    t = time.time(); x = ll(p); t1 = time.time() - t
    t = time.time(); c, y = duo(p); t2 = time.time() - t
    t = time.time(); r, z = pas(p); t3 = time.time() - t
    print('  p = %d (%d digits): Lucas-Lehmer %s %.2f | two-cell wheel %s %.2f | Pascal on %d cells %s %.2f' % (p, int(p*math.log10(2)) + 1, x, t1, y, t2, r, z, t3))
