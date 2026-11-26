"""Referee checks U1 (Hurwitz = sum_m m^-1 L(1, Kronecker)), U4 (factorisation, periodicity, Jacobi coefficients),
U6 (local sums and E(r)). Own code; plain Python integers."""
import math, random, cmath
from fractions import Fraction
from sympy import factorint, isprime, primerange

random.seed(7)


def jacobi(a, n):
    assert n > 0 and n % 2 == 1
    a %= n
    t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    return t if n == 1 else 0


def kron(D, n):
    """Kronecker symbol (D/n), n >= 1."""
    r = 1
    while n % 2 == 0:
        n //= 2
        if D % 2 == 0:
            return 0
        r *= 1 if D % 8 in (1, 7) else -1
    return r * jacobi(D, n)


def fund(Delta):
    """Delta < 0 discriminant -> (Delta0, f) with Delta = Delta0 f^2, Delta0 fundamental."""
    a = -Delta
    fac = factorint(a)
    core = 1; f = 1
    for l, e in fac.items():
        core *= l ** (e % 2); f *= l ** (e // 2)
    D0 = -core
    if D0 % 4 != 1:
        D0 *= 4; assert f % 2 == 0; f //= 2
    assert D0 * f * f == Delta and D0 % 4 in (0, 1)
    return D0, f


def L1_prim_odd(D0):
    """L(1, chi_D0) for fundamental D0 < 0 by the finite formula -pi/|D0|^(3/2) sum a chi(a)."""
    k = -D0
    s = sum(a * kron(D0, a) for a in range(1, k + 1))
    return -math.pi * s / k ** 1.5


def L1_kron(Delta):
    D0, f = fund(Delta)
    v = L1_prim_odd(D0)
    for l in factorint(f):
        v *= 1 - kron(D0, l) / l
    return v


def hurwitz_forms(D):
    """H(D): classes of all positive definite forms of discriminant -D, weights 1/2 (ax^2+ay^2-type) and 1/3."""
    H = Fraction(0)
    a = 1
    while 3 * a * a <= D:
        for b in range(-a + 1, a + 1):
            if (b * b + D) % (4 * a):
                continue
            c = (b * b + D) // (4 * a)
            if c < a or (c == a and b < 0):
                continue
            w = Fraction(1)
            if a == b == c:
                w = Fraction(1, 3)
            elif a == c and b == 0:
                w = Fraction(1, 2)
            H += w
        a += 1
    return H


def Mset(D):
    out = []
    m = 1
    while m * m <= D:
        if D % (m * m) == 0 and (-D // (m * m)) % 4 in (0, 1):
            out.append(m)
        m += 1
    return out


# ---------------- U1 on 300 values of D
Ds = set()
while len(Ds) < 120:
    D = random.randint(3, 30000)
    if D % 4 in (0, 3):
        Ds.add(D)
for k in range(1, 61):  # big square parts
    for base in (3, 4, 7, 8, 11, 15, 20, 23):
        D = base * k * k
        if D <= 60000 and D % 4 in (0, 3):
            Ds.add(D)
for N in (13, 17, 29, 37):
    for p in primerange(N + 1, 6 * N):
        for r in range(0, 6):
            D = N * (4 * p - r * r * N)
            if D > 0 and len(Ds) < 300:
                Ds.add(D)
Ds = sorted(Ds)[:300] if len(Ds) > 300 else sorted(Ds)
worst = 0.0; nonfund = 0
for D in Ds:
    S = sum(L1_kron(-D // (m * m)) / m for m in Mset(D))
    H = hurwitz_forms(D)
    worst = max(worst, abs(math.sqrt(D) / math.pi * S - float(H)))
    nonfund += len(Mset(D)) > 1 or fund(-D)[1] > 1
print(f"U1: {len(Ds)} D (max {max(Ds)}), {nonfund} with non-fundamental -D or several m; max |(sqrtD/pi) L* - H| = {worst:.2e}")
assert worst < 1e-8

# ---------------- U4: factorisation T(m,n;p) = (-N/q) T(m,n'';p) chi_q(4p - r^2 N), periodicity mod K = 2 n'' m^2


def T(m, n, p, r, N):
    D = N * (4 * p - r * r * N)
    if D <= 0 or D % (m * m):
        return 0
    Dl = -D // (m * m)
    if Dl % 4 not in (0, 1):
        return 0
    return kron(Dl, n)


def split(m, n, r, N):
    f = factorint(n)
    bad = set(factorint(2 * m * max(r, 1) * N)) if r >= 1 else set(factorint(2 * N))
    q = 1; qn = 1
    for l, e in f.items():
        if e % 2 == 1 and l not in bad:
            q *= l; qn *= l ** e
    return q, n // qn


Ns = [p for p in primerange(11, 400)]
cnt = 0; bad = 0; badper = 0; badcop = 0
for trial in range(6000):
    N = random.choice(Ns); r = random.randint(0, 6)
    m = random.choice([1, 1, 1, 2, 3, 4, 5, 6, 9, 10, 15])
    if r == 0:
        m = random.choice([1, 2])
    n = random.randint(1, 3000)
    q, n2 = split(m, n, r, N)
    K = 2 * n2 * m * m
    if math.gcd(K, q) != 1:
        badcop += 1
    lo = max(1, (r * r * N) // 4 + 1)
    for _ in range(4):
        p = random.randint(lo, lo + 50 * K)
        if p % N == 0:
            continue
        lhs = T(m, n, p, r, N)
        u = 4 * p - r * r * N
        rhs = jacobi(-N, q) * T(m, n2, p, r, N) * jacobi(u, q)
        cnt += 1
        bad += lhs != rhs
        p2 = p + K * random.randint(1, 5)
        badper += T(m, n2, p, r, N) != T(m, n2, p2, r, N)
print(f"U4: {cnt} (m,n,p,r,N) tests: factorisation failures {bad}, period-K failures {badper}, (K,q)>1 cases {badcop}")
assert bad == 0 and badper == 0 and badcop == 0


# ---------------- U4(d): coefficients of x -> chi_q(4x - c), q odd squarefree, (q, c) = 1
def char_coeffs_prime(l, c):
    g = next(g for g in range(2, l) if all(pow(g, (l - 1) // t, l) != 1 for t in factorint(l - 1)))
    log = {}
    x = 1
    for e in range(l - 1):
        log[x] = e; x = x * g % l
    out = []
    for k in range(l - 1):  # psi_k(g^e) = exp(2 pi i k e/(l-1))
        s = sum(jacobi(4 * x - c, l) * cmath.exp(-2j * math.pi * k * log[x] / (l - 1)) for x in range(1, l))
        out.append(s / (l - 1))
    return out


ok = True; info = []
for l in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]:
    for c in [1, 2, 5, 7 * 7 * 13, 4 * 101]:
        if c % l == 0:
            continue
        co = char_coeffs_prime(l, c)
        ok &= abs(co[0] + jacobi(-c, l) / (l - 1)) < 1e-12
        ok &= max(abs(z) for z in co[1:]) <= math.sqrt(l) / (l - 1) + 1e-12
        # exactly one coefficient of modulus 1/(l-1) among non-principal: the Legendre symbol k = (l-1)/2
        ok &= abs(abs(co[(l - 1) // 2]) - 1 / (l - 1)) < 1e-12
        ok &= all(abs(abs(z) - math.sqrt(l) / (l - 1)) < 1e-9 for k, z in enumerate(co) if k not in (0, (l - 1) // 2))
print("U4(d) prime moduli: ahat(psi0) = -(-c/l)/(l-1), |ahat| = sqrt(l)/(l-1) except Legendre (1/(l-1)):", ok)
assert ok


# ---------------- U6: e(m,n) by direct enumeration over (Z/P)^x, multiplicativity, local sums vs E_l(r)
def e_direct(m, n, r, N, P):
    s = 0; c = 0
    for p in range(P):
        if math.gcd(p, P) != 1:
            continue
        # T depends on p mod P only through the residue; use p + P*t with t large enough that D > 0
        pp = p + P * (r * r * N // (4 * P) + 2)
        s += T(m, n, pp, r, N); c += 1
    return Fraction(s, c)


N = 13
okm = True; tests = 0
for r in range(0, 4):
    for (m, n) in [(1, 3), (1, 5), (1, 15), (1, 9), (1, 2), (1, 6), (1, 10), (1, 4), (2, 3), (3, 2), (3, 5), (2, 15), (1, 45),
                   (6, 5), (1, 21), (3, 3), (1, 30)]:
        if r == 0 and m > 2:
            continue
        P = 8 * n * m * m * (1 if N % 1 else 1)
        e = e_direct(m, n, r, N, P)
        # product of local factors
        prod = Fraction(1)
        for l in set(factorint(2 * m * n)):
            ml = l ** (factorint(m).get(l, 0)); nl = l ** (factorint(n).get(l, 0))
            Pl = 8 * l ** (2 * factorint(m).get(l, 0) + factorint(n).get(l, 0) + 2)
            prod *= e_direct(ml, nl, r, N, Pl)
        okm &= e == prod; tests += 1
print(f"U6 multiplicativity e(m,n) = prod_l e_l (N=13, r=0..3): {tests} cases ok={okm}")
assert okm


def local_sum(l, r, N, I):
    """sum_{i<=I} l^-i mean_p 1[l^2i | D, disc] (1 - kappa/l)^-1 ; exact Fraction."""
    tot = Fraction(0)
    for i in range(I + 1):
        mod = l ** (2 * i + 1) if l != 2 else 2 ** (2 * i + 4)
        s = Fraction(0); c = 0
        for p in range(mod):
            if p % l == 0:
                continue
            pp = p + mod * (r * r * N // (4 * mod) + 2)
            D = N * (4 * pp - r * r * N)
            c += 1
            if D % l ** (2 * i):
                continue
            Dl = -D // l ** (2 * i)
            if l == 2 and Dl % 4 not in (0, 1):
                continue
            kap = kron(Dl, l)
            s += 1 / (1 - Fraction(kap, l))
        tot += Fraction(1, l ** i) * s / c
    return tot


def E_closed(l, r, N):
    if l == N:
        return Fraction(1)
    if l == 2:
        return Fraction(2, 3) if r % 2 else Fraction(4, 3)
    if r % l == 0:
        return Fraction(l * l, l * l - 1)
    return 1 - Fraction(1, (l - 1) ** 2 * (l + 1))


for N in (10007, 20011):
    worst = 0
    for r in range(0, 4):
        for l, I in [(2, 5), (3, 4), (5, 3), (7, 2), (11, 2), (13, 1)]:
            v = local_sum(l, r, N, I)
            tail = 4 * Fraction(1, l ** (3 * (I + 1)))  # sum_{i>I} l^-i P(l^2i | u) (1-1/l)^-1 crude: P <= 2 l^(1-2i)
            d = abs(float(v - E_closed(l, r, N)))
            worst = max(worst, d)
            assert d <= float(tail) * l + 1e-15, (N, r, l, d)
    print(f"U6 local sums = E_l(r) closed forms, N={N}, r=0..3, l<=13: max diff {worst:.2e} (within the i-tail)")

# global E(r): nu'(r) and zeta(2)(1 - N^-2)
N = 20011
for r in range(0, 4):
    E = 1.0; nu = 1.0
    for l in primerange(2, 3 * 10 ** 6):
        E *= float(E_closed(l, r, N))
        nu *= (l * l / (l * l - 1)) if (r % l == 0) else (1 - 1 / ((l - 1) ** 2 * (l + 1)))
    if r == 0:
        print(f"E(0) = {E:.6f}, zeta(2)(1-N^-2) = {math.pi**2/6*(1-N**-2):.6f}")
    else:
        print(f"E({r}) = {E:.6f}, nu'({r}) = {nu:.6f}")
