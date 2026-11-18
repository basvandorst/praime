"""The curve wheel (Goldwasser and Kilian, 1986), as on machine 16 of the page, in plain Python.
For n (not divisible by 2 or 3) take a curve y^2 = x^3 + a x + b mod n through the point P = (a + 1, 1), a = 1 .. 12
(for n below 60, where few curves have a suitable size, after that every a and every start (x0, 1)).
Count its points one by one as if n were prime: m = 1 + number of (x, y) with y^2 = x^3 + a x + b. If mP is not the point at infinity,
n is not prime. If it is, let q be the largest prime factor of m; when q > (n^(1/4) + 1)^2 and (m/q)P is not the point at infinity,
n is prime (proved: mod a prime factor p <= sqrt(n) the curve would have at most (sqrt(p) + 1)^2 < q points). Otherwise try the next curve.
A denominator that shares a factor with n also shows that n is not prime."""
import sys, math, json
from sympy import isprime, primerange
class Factor(Exception): pass
def inv(v, n):
    g = math.gcd(v % n, n)
    if g != 1: raise Factor(g)
    return pow(v, -1, n)
def add(P, Q, a, n):
    if P is None: return Q
    if Q is None: return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % n == 0: return None
        if P[1] != Q[1]: raise Factor(math.gcd(P[1] - Q[1], n))
        l = (3*P[0]*P[0] + a) * inv(2*P[1], n) % n
    else:
        l = (Q[1] - P[1]) * inv(Q[0] - P[0], n) % n
    x = (l*l - P[0] - Q[0]) % n
    return (x, (l*(P[0] - x) - P[1]) % n)
def mulp(k, P, a, n):
    R = None; B = P
    while k:
        if k & 1: R = add(R, B, a, n)
        k >>= 1
        if k:
            B = add(B, B, a, n)
            if B is None: return R
    return R
def largest_prime_factor(m):
    q = 1; d = 2
    while d*d <= m:
        while m % d == 0: q = d; m //= d
        d += 1
    return max(q, m) if m > 1 else q
def curves(n):
    """the curves tried: a = 1 .. 12 through (a + 1, 1); for n below 60 after that every a and every x0"""
    for a in range(1, 13): yield a, (a + 1) % n
    if n < 60:
        for a in range(n):
            for x0 in range(n): yield a, x0
def curve_wheel(n):
    """returns (verdict, curves tried, m, q): verdict 'prime' (proved), 'not prime' (certain) or 'undecided'"""
    if n in (2, 3): return 'prime', 0, 0, 0
    if n % 2 == 0 or n % 3 == 0: return 'not prime', 0, 0, 0
    cnt = [0]*n
    for y in range(n): cnt[y*y % n] += 1
    bound = (n**0.25 + 1)**2
    for a, x0 in curves(n):
        b = (1 - x0**3 - a*x0) % n; P = (x0, 1 % n)
        g = math.gcd((4*a**3 + 27*b*b) % n, n)
        if g == n: continue
        if g > 1: return 'not prime', a, 0, 0
        m = 1 + sum(cnt[(x*x*x + a*x + b) % n] for x in range(n))
        try:
            if mulp(m, P, a, n) is not None: return 'not prime', a, m, 0
            q = largest_prime_factor(m)
            if q <= bound: continue
            if mulp(m//q, P, a, n) is None: continue
        except Factor:
            return 'not prime', a, m, 0
        return 'prime', a, m, q
    return 'undecided', 12, 0, 0
if __name__ == "__main__":
    LIM = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    bad = und = 0; tries = {}; lit = []
    for n in range(2, LIM + 1):
        v, a, m, q = curve_wheel(n)
        if v == 'undecided': und += 1
        if (v == 'prime') != isprime(n) and v != 'undecided': bad += 1
        if v == 'prime': tries[a] = tries.get(a, 0) + 1; lit.append(n)
    print('curve wheel, n = 2 ..', LIM, ': wrong answers', bad, '| not decided', und, '| primes proved', len(lit), '| curve that gave the proof (a: primes):', dict(sorted(tries.items())))
    for n in (241, 4999, 65537, 999983, 1000003, 999985, 1000001):
        v, a, m, q = curve_wheel(n)
        print('  n = %d: %s%s' % (n, v, (' (curve a = %d, %d positions, check prime %d)' % (a, m, q)) if v == 'prime' else ''))
    MX = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127}; lit2 = [n for n in lit if n <= 250]; best = {'p': (0, 0), 'm': (0, 0), 'w': (0, 0)}; cur = []
    for t in range(2, 251):
        if t in lit2: cur.append(t)
        p = cur; mm = [k for k in p if k in MX]; w = [k for k in p if isprime(k - 2) or isprime(k + 2)]
        for key, arr in (('p', p), ('m', mm), ('w', w)):
            if len(arr) > best[key][0]: best[key] = (len(arr), t)
    row = list(best['p'] + best['m'] + best['w'])
    print('curve wheel, 250 wheels: primes %d (tick %d), Mersenne exponents %d (tick %d), twin primes %d (tick %d)' % tuple(row))
    json.dump({'row250': row}, open('curve_wheel.json', 'w'))
