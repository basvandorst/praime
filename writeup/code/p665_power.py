# Puzzle 665: (p + q)(q + r) = m^n for consecutive primes p < q < r.
# Finds every pair (A, g) with A(A + g) a perfect power, g even, g <= G, A <= AMAX, and keeps those with
# A = p + q, A + g = q + r.  usage: p665_power.py G AMAX [BRUTE]     e.g.  p665_power.py 3000 1e100
# With BRUTE the result is compared with a direct test of every A <= BRUTE.
# d = gcd(A, A+g) | g, A = d*u, A+g = d*v, gcd(u,v) = 1, d^2 u v = m^n.  A prime l not dividing d divides at most one
# of u, v, to a power divisible by n.  So u = alpha x^n, v = beta y^n with alpha, beta n-th-power-free, coprime,
# built from primes of d, and v_l(alpha) + v_l(beta) = -2 v_l(d) mod n.  Equation: beta y^n - alpha x^n = c = g/d.
import sys, math, gmpy2
from sympy import factorint, primerange
from decimal import Decimal
G = int(sys.argv[1]); AMAX = int(Decimal(sys.argv[2]))
NMAX = 2 * AMAX.bit_length() + 2
K = 2 * len(str(AMAX)) + 40                                  # decimal digits for the continued fraction
sols = set(); neq = 0; ncf = 0
def record(A, g, n):
    if A >= 1 and A <= AMAX and gmpy2.iroot(A * (A + g), n)[1]: sols.add((A, g, n))
def solve(alpha, beta, c, n, d, g):
    """all x >= 1 with beta y^n - alpha x^n = c, d*alpha*x^n <= AMAX"""
    global ncf
    xmax = int(gmpy2.iroot(AMAX // (d * alpha), n)[0])
    if xmax < 1: return
    if alpha == beta:                                        # = 1: y^n - x^n = c >= n x^(n-1)
        lim = min(xmax, int(gmpy2.iroot(c // n + 1, n - 1)[0]) + 1)
    else:                                                    # beyond this x, y/x is a convergent of (alpha/beta)^(1/n)
        lim = min(xmax, max(int(gmpy2.iroot(2 * c, n - 2)[0]) + 2 if n > 2 else xmax, 50))
    for x in range(1, lim + 1):
        t = alpha * x ** n + c
        if t % beta == 0 and gmpy2.iroot(t // beta, n)[1]: record(d * alpha * x ** n, g, n)
    if alpha == beta or lim >= xmax: return
    ncf += 1
    num = int(gmpy2.iroot(alpha * beta ** (n - 1) * 10 ** (n * K), n)[0]); den = beta * 10 ** K   # theta = num/den
    h0, h1, k0, k1 = 0, 1, 1, 0
    while den and k1 <= xmax:
        a, r = divmod(num, den); num, den = den, r
        h0, h1, k0, k1 = h1, a * h1 + h0, k1, a * k1 + k0
        if k1 < 1 or k1 > xmax: continue
        val = beta * h1 ** n - alpha * k1 ** n
        if val > 0 and c % val == 0:
            t, ex = gmpy2.iroot(c // val, n)
            if ex: record(d * alpha * (int(t) * k1) ** n, g, n)
primes_n = list(primerange(2, NMAX + 1))
for g in range(2, G + 1, 2):
    fg = factorint(g)
    divs = [1]
    for l, e in fg.items(): divs = [x * l ** j for x in divs for j in range(e + 1)]
    for d in divs:
        c = g // d; fd = factorint(d)
        for n in primes_n:
            if d * 2 ** 0 > AMAX: break
            parts = [(1, 1)]
            for l, v in fd.items():
                e = (-2 * v) % n
                parts = [(a * l ** e, b) for a, b in parts] + ([(a, b * l ** e) for a, b in parts] if e else [])
            for alpha, beta in parts:
                neq += 1; solve(alpha, beta, c, n, d, g)
print(f"G = {G}, AMAX = {AMAX}: {neq} equations, {ncf} continued fractions")
print(f"  {len(sols)} solutions (A, g, n); the largest A:")
for A, g, n in sorted(sols)[-8:]: print("    A =", A, "g =", g, "n =", n)
from sympy import isprime, prevprime, nextprime
found = []
for A, g, n in sorted(sols):                                 # A = p + q with consecutive primes p < q, and q + r = A + g
    if A < 5: continue
    p = prevprime(A // 2 + 1) if A // 2 + 1 > 2 else None    # p <= A/2
    if p is None: continue
    q = nextprime(p); r = nextprime(q)
    if p + q == A and q + r == A + g: found.append((p, q, r, n))
print("  from three consecutive primes:", found)
if len(sys.argv) > 3:                                         # compare with a direct test of every A <= BRUTE
    BR = int(sys.argv[3]); direct = set()
    for g in range(2, G + 1, 2):
        for A in range(1, BR + 1):
            N = A * (A + g)
            for n in primes_n:
                if 2 ** n > N: break
                if gmpy2.iroot(N, n)[1]: direct.add((A, g, n))
    mine = {t for t in sols if t[0] <= BR}
    assert mine == direct, (sorted(mine ^ direct)[:10])
    print(f"  direct test of every A <= {BR}: the same {len(direct)} solutions")
