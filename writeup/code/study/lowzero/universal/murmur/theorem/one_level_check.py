"""one_level_check.py: checks of the numeric claims of murmurations_one_level.md (stand-alone, own code).

Run:  (ulimit -v 1200000; timeout 600 nice python -I -B one_level_check.py)      [sections: argv[1], default ABCDE]

 A  class numbers and the discriminants D = N(4p - r^2 N) of the family, N squarefree (steps 2-3 and U2-U4 of the note):
    A1 brute-force Hurwitz numbers = PARI; A2 H(D) = sum_{m in M(D)} h'(-D/m^2) with own reduced-form counts;
    A3 h'(Delta) = sqrt|Delta| L(1, xi_Delta)/pi (digamma formula); A4 N | Delta_0, (m, N) = 1, v_2(Delta_0) in {2, 3}
    for even N, F_l = 1 at l | N; A5 T(m, n; p) = 0 for (n, N) > 1, period n m^2, the factorisation of U4;
    A6 coefficients of chi_q(4x - c) in characters of x; A7 the weight bound w_d of U5(iii).
 B  local means by exact enumeration: B1 E_l^(N)(r) for l = 2, 3, 5, 7 at levels prime, 2q, 3q, 6q, 5q, 2310;
    B2 the local coefficients e_l(i, j) of U6 and their sums.
 C  the trace formula on the new space at four small squarefree levels (prime, 2q, 3q, 6q) against PARI, k = 2, 4, 6;
    old forms give trace 0; root numbers; dimensions.
 D  the densities and Euler products: C0, nu', E^(N)(0); special cases N prime, N = 2q; the average over squarefree N
    (Zubrilina's constants exactly, prime by prime, and by brute force); the repaired comparison M^(N) vs M^(S).
 E  one level of each type below 10^4 not used in the earlier notes (9973, 9998, 9987, 9978), weight 2, all p <= 4N:
    window means against the level's own density and every other density; class-number means; the truncation of U2.
    Decision rule fixed before the run: own chi2 <= 32 on 16 windows (errors: single-value standard errors, which
    overstate the error of a window mean) and every other density (the three other types and Zubrilina's average)
    at least 9 above the own chi2; class-number means |z| < 4.
No level in 30030..30530 is touched (blind test of round 163).
"""
import sys, math, time, random
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import gmpy2
import cypari2
from scipy.special import digamma

pari = cypari2.Pari()
pari.allocatemem(300 * 10 ** 6)
T0 = time.time()
NCHK = NFAIL = 0
FORBIDDEN = range(30030, 30531)


def check(name, ok, info=''):
    global NCHK, NFAIL
    NCHK += 1
    NFAIL += (not ok)
    print(('PASS ' if ok else 'FAIL ') + name + ('  ' + info if info else ''), flush=True)


def sqfree(n): return n > 1 and all(e == 1 for e in sp.factorint(n).values())
def primes_of(n): return tuple(sorted(sp.factorint(n)))
def omega(n): return len(sp.factorint(n))


# ============================================================ A
def hurwitz_bf(D):
    """H(D): classes of positive definite forms of discriminant -D, (a,0,a) weighted 1/2, (a,a,a) weighted 1/3."""
    s = Fr(0); a = 1
    while 3 * a * a <= D:
        for b in range(-a + 1, a + 1):
            if (b * b + D) % (4 * a):
                continue
            c = (b * b + D) // (4 * a)
            if c < a or (b < 0 and c == a):
                continue
            s += Fr(1, 2) if (b == 0 and a == c) else Fr(1, 3) if (b == a == c) else 1
        a += 1
    return s


def hprime_bf(Delta):
    """h'(Delta) = 2 h(Delta)/w(Delta), h = number of reduced primitive forms of discriminant Delta < 0."""
    D = -Delta; h = 0; a = 1
    while 3 * a * a <= D:
        for b in range(-a + 1, a + 1):
            if (b * b + D) % (4 * a):
                continue
            c = (b * b + D) // (4 * a)
            if c < a or (b < 0 and c == a) or math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            h += 1
        a += 1
    w = 6 if Delta == -3 else 4 if Delta == -4 else 2
    return Fr(2 * h, w)


def Mset(D):
    """M(D) = {m >= 1 : m^2 | D, -D/m^2 = 0, 1 mod 4}; m^2 | D iff m | f, D = c f^2 with c squarefree (PARI core)."""
    f = int(pari('core(%d,1)[2]' % D))
    return sorted(m for m in sp.divisors(f) if (-D // (m * m)) % 4 in (0, 1))


def H(D): return Fr(str(pari('qfbhclassno(%d)' % D)))


def T_fn(m, n, p, r, N):
    """T(m, n; p) = 1[m in M(D)] (Delta_m / n), D = N(4p - r^2 N), Delta_m = -D/m^2 (Kronecker symbol)."""
    u = 4 * p - r * r * N
    if u <= 0:
        return None
    D = N * u
    if D % (m * m):
        return 0
    Dm = -D // (m * m)
    if Dm % 4 not in (0, 1):
        return 0
    return gmpy2.kronecker(Dm, n)


def split_n(n, m, r, N):
    """q = prod{l odd : v_l(n) odd, l !| m r N} (r >= 1) or {l odd: v_l(n) odd, l !| m N} (r = 0); n'' = n / q-part."""
    q = 1; qn = 1
    for l, e in sp.factorint(n).items():
        if l == 2 or e % 2 == 0 or m % l == 0 or N % l == 0 or (r > 0 and r % l == 0):
            continue
        q *= l; qn *= l ** e
    return q, n // qn


def sec_A():
    bad = sum(hurwitz_bf(D) != H(D) for D in range(3, 801) if D % 4 in (0, 3))
    check('A1 brute-force Hurwitz numbers (reduced forms, weights 1/2, 1/3) = PARI qfbhclassno, D <= 800', bad == 0,
          '%d values' % len([D for D in range(3, 801) if D % 4 in (0, 3)]))
    bad = n = 0
    for D in range(3, 3001):
        if D % 4 not in (0, 3):
            continue
        n += 1
        bad += sum(hprime_bf(-D // (m * m)) for m in Mset(D)) != H(D)
    check('A2 U1: H(D) = sum_{m in M(D)} h\'(-D/m^2), h\' = 2h/w from own counts of primitive reduced forms, D <= 3000',
          bad == 0, '%d discriminants, exact' % n)
    worst = 0.0; n = 0
    for D in range(3, 1201):
        if D % 4 not in (0, 3):
            continue
        for m in Mset(D):
            Dl = -D // (m * m); q = -Dl
            a = np.arange(1, q)
            chi = np.array([gmpy2.kronecker(Dl, int(x)) for x in a], dtype=float)
            L1 = -float(np.sum(chi * digamma(a / q))) / q
            worst = max(worst, abs(math.sqrt(q) * L1 / math.pi - float(hprime_bf(Dl)))); n += 1
    check('A3 class number formula for orders h\'(Delta) = sqrt|Delta| L(1, xi_Delta)/pi, xi = Kronecker (Delta/.), '
          'L(1) = -(1/q) sum xi(a) psi(a/q), every Delta = -D/m^2, D <= 1200', worst < 1e-9, '%d discriminants, max |diff| %.1e' % (n, worst))
    # A4: the discriminants of the family at squarefree N
    levels = [N for N in range(2, 131) if sqfree(N)] + [210, 330, 390, 462, 510, 546, 1155, 2310, 4290, 39270]
    n = bad = 0; ev = set()
    for N in levels:
        assert N not in FORBIDDEN
        for p in sp.primerange(2, 400):
            if N % p == 0:
                continue
            r = 0
            while r * r * N < 4 * p:
                D = 4 * p * N - r * r * N * N
                d0 = int(pari('coredisc(%d)' % (-D)))
                f0 = math.isqrt(-D // d0); assert d0 * f0 * f0 == -D
                ok = d0 % N == 0 and math.gcd(f0, N) == 1 and all(math.gcd(m, N) == 1 for m in Mset(D))
                if N % 2 == 0:
                    v2 = (d0 & -d0).bit_length() - 1
                    ok = ok and v2 in (2, 3); ev.add(v2)
                ok = ok and all(gmpy2.kronecker(d0, l) == 0 for l in primes_of(N))      # chi(l) = 0, l !| f0: F_l = 1
                bad += not ok; n += 1; r += 1
    check('A4 squarefree N (all N <= 130 and ten with 3-6 prime factors, up to 39270), all p < 400, p !| N, all r: N divides '
          'the fundamental discriminant Delta_0 of -D, l !| f_0 and F_l(D) = 1 for l | N, every m in M(D) prime to N, '
          'v_2(Delta_0) in {2, 3} for even N', bad == 0, '%d cases, v_2 seen %s' % (n, sorted(ev)))
    # A5: T = 0 for (n, N) > 1; period n m^2; factorisation
    rnd = random.Random(20261126)
    Ns = [N for N in range(2, 400) if sqfree(N)] + [2310, 4290, 39270, 9978, 9998]
    nz = bz = npd = bpd = nf = bf = bK = 0
    for _ in range(40000):
        N = rnd.choice(Ns); r = rnd.randrange(0, 5); m = rnd.choice([1, 1, 1, 2, 3, 4, 5, 6, 9, 15]); n = rnd.randrange(1, 400)
        p = int(sp.nextprime(r * r * N // 4 + rnd.randrange(1, 50 * N)))
        if N % p == 0:
            continue
        t = T_fn(m, n, p, r, N)
        if math.gcd(n, N) > 1:
            bz += t != 0; nz += 1
        # period n m^2 (p any integer with u > 0)
        P = n * m * m; x = p + rnd.randrange(0, 3 * P); kk = rnd.randrange(1, 4)
        if T_fn(m, n, x, r, N) is not None:
            bpd += T_fn(m, n, x, r, N) != T_fn(m, n, x + kk * P, r, N); npd += 1
        q, n2 = split_n(n, m, r, N)
        if r >= 1:
            fac = gmpy2.jacobi((-N) % q, q) * gmpy2.jacobi((4 * p - r * r * N) % q, q) if q > 1 else 1
        else:
            fac = gmpy2.jacobi((-N) % q, q) * gmpy2.jacobi(p % q, q) if q > 1 else 1
        bf += t != fac * T_fn(m, n2, p, r, N); nf += 1
        bK += math.gcd(n2 * m * m, q) != 1
    check('A5 T(m,n;p) = 0 when (n, N) > 1', bz == 0, '%d cases' % nz)
    check('A5 T(m,n;.) periodic modulo n m^2 (the period 2 n m^2 of the earlier note is not needed)', bpd == 0, '%d tests' % npd)
    check('A5 factorisation T(m,n;p) = (-N/q) chi_q(4p - r^2 N) T(m,n\'\';p) (r = 0: chi_q(p)), (n\'\' m^2, q) = 1, '
          'random squarefree N < 400 and 2310, 4290, 39270, 9978, 9998', bf == 0 and bK == 0, '%d tests' % nf)
    # A6: chi_l(4x - c) = sum_psi ahat(psi) psi(x) on units, l prime
    bad = n = 0; worst = 0.0
    for l in sp.primerange(3, 62):
        g = int(sp.primitive_root(l)); lg = {pow(g, j, l): j for j in range(l - 1)}
        xs = np.arange(1, l); js = np.array([lg[int(x)] for x in xs])
        for N in (7, 15, 26, 42, 105, 2310):
            for r in (0, 1, 2, 3, 5):
                c = r * r * N
                f = np.array([int(sp.legendre_symbol((4 * int(x) - c) % l, l)) if (4 * int(x) - c) % l else 0 for x in xs], float)
                ah = [np.sum(f * np.exp(-2j * np.pi * t * js / (l - 1))) / (l - 1) for t in range(l - 1)]
                rec = sum(ah[t] * np.exp(2j * np.pi * t * js / (l - 1)) for t in range(l - 1))
                ok = np.max(np.abs(rec - f)) < 1e-9
                if c % l:
                    ok &= abs(ah[0] - (-int(sp.legendre_symbol((-c) % l, l)) / (l - 1))) < 1e-12
                    ok &= max(abs(a) for a in ah) <= math.sqrt(l) / (l - 1) + 1e-12
                    worst = max(worst, max(abs(a) for a in ah[1:]) * (l - 1) / math.sqrt(l))
                else:                       # chi_l(4x) = chi_l(x): ahat = [psi = Legendre], index (l-1)/2
                    ok &= all(abs(ah[t] - (1.0 if t == (l - 1) // 2 else 0.0)) < 1e-9 for t in range(l - 1))
                bad += not ok; n += 1
    # composite q: mean over (Z/q)^x of jacobi(4x - c, q) = prod_l ahat_l(psi_0)
    for q in (105, 1155, 255255):
        for N, r in ((38, 1), (46, 2), (19, 1), (2 * 19 * 23, 1)):
            c = r * r * N
            if math.gcd(c, q) > 1:
                continue
            xs = [x for x in range(1, q) if math.gcd(x, q) == 1]
            mean = Fr(sum(gmpy2.jacobi((4 * x - c) % q, q) for x in xs), len(xs))
            pred = Fr(1)
            for l in primes_of(q):
                pred *= Fr(-int(sp.legendre_symbol((-c) % l, l)), l - 1)
            bad += mean != pred; n += 1
    check('A6 chi_l(4x - c) = sum_psi ahat(psi) psi(x) on units; ahat(psi_0) = -(-c/l)/(l-1), |ahat(psi)| <= sqrt(l)/(l-1) '
          '(l !| c); ahat = [psi = Legendre] (l | c); composite q: mean = prod of local ahat(psi_0) (exact)', bad == 0,
          '%d cases (l < 62, c = r^2 N, N up to 2310), max |ahat|(l-1)/sqrt(l) over psi != psi_0: %.3f' % (n, worst))
    # A7: w_d = sum_{q in (Q,2Q], d | q} q^(-1/2)/phi(q) <= c0 Q^(-1/2)/phi(d), c0 = 2 zeta(2) zeta(3)/zeta(6)
    zz = float(sp.zeta(2) * sp.zeta(3) / sp.zeta(6))
    X = 2 ** 15
    ph = np.arange(2 * X + 1, dtype=float)                  # Euler phi by a sieve
    for l in sp.primerange(2, 2 * X + 1):
        ph[l::l] *= 1 - 1.0 / l
    ph = np.rint(ph)
    xs = np.arange(1, 2 * X + 1)
    cs = np.cumsum(xs / ph[1:])
    ok1 = bool(np.all(cs <= zz * xs + 1e-9))
    worst = 0.0
    for e in range(1, 15):
        Q = 2 ** e
        for d in range(1, 2 * Q + 1):
            qs = np.arange(d * (Q // d + 1), 2 * Q + 1, d)
            w = float(np.sum(qs ** -0.5 / ph[qs]))
            worst = max(worst, w * ph[d] * math.sqrt(Q))
    check('A7 sum_{n<=X} n/phi(n) <= (zeta(2)zeta(3)/zeta(6)) X (X <= 65536) and phi(d) Q^(1/2) w_d <= 2 zeta(2)zeta(3)/zeta(6) '
          '= %.4f (Q = 2..2^14, all d <= 2Q)' % (2 * zz), ok1 and worst <= 2 * zz, 'max phi(d) Q^(1/2) w_d = %.4f' % worst)


# ============================================================ densities (used in B, D, E)
PR = [int(p) for p in sp.primerange(2, 10 ** 6)]


def u_(l): return 1.0 / ((l - 1) ** 2 * (l + 1))


C0 = math.exp(sum(math.log1p(-u_(l)) for l in PR))       # tail beyond 1e6: sum u_l < 1e-12
Z2 = math.pi ** 2 / 6


def e_loc(l, r): return l * l / (l * l - 1.0) if r % l == 0 else 1.0 - u_(l)


def e_loc_exact(l, r):
    return Fr(l * l, l * l - 1) if r % l == 0 else 1 - Fr(1, (l - 1) ** 2 * (l + 1))


def nu_prime(r):
    v = C0
    for l in sp.factorint(r):
        v *= (l * l / (l * l - 1.0)) / (1.0 - u_(l))
    return v


def nu_N(r, S):
    v = nu_prime(r)
    for l in S:
        v /= e_loc(l, r)
    return v


def kappa(S): return float(np.prod([l / (l - 1.0) for l in S])) if S else 1.0
def piN(S): return float(np.prod([1 + 1.0 / l for l in S])) if S else 1.0
def E0(S): return Z2 * (float(np.prod([1 - 1.0 / l ** 2 for l in S])) if S else 1.0)


def Ucheb(n, t):
    if n == 0:
        return 1.0
    a, b = 1.0, 2 * t
    for _ in range(n - 1):
        a, b = b, 2 * t * b - a
    return b


def M_level(y, S, k=2):
    """M_k^(N)(y) = (1/(k-1)) [2 pi pi(N) sqrt y + (-1)^(k/2-1) (12/pi) kappa(N) sum_r nu_N(r) U_{k-2}(r/2sqrt y) sqrt(4y-r^2)]
    - 12 kappa(N) y [k = 2]; S = set of primes of N (S = () gives M_k of the prime-level theorem, q -> infinity)."""
    sg = (-1) ** (k // 2 - 1)
    out = 2 * math.pi * piN(S) * math.sqrt(y)
    r = 1
    while r * r < 4 * y:
        out += sg * (12 / math.pi) * kappa(S) * nu_N(r, S) * Ucheb(k - 2, r / (2 * math.sqrt(y))) * math.sqrt(4 * y - r * r)
        r += 1
    out /= (k - 1)
    if k == 2:
        out -= 12 * kappa(S) * y
    return out


def _ep(fn): return math.exp(sum(math.log(fn(float(p))) for p in PR))


Z_ALPHA = 2 * math.pi * _ep(lambda l: (l ** 4 - 2 * l ** 2 - l + 1) / (l ** 4 - 2 * l ** 2 + l))
Z_BETA = 2 * math.pi * _ep(lambda l: (l ** 3 + l ** 2 - 1) / (l * (l ** 2 + l - 1)))
Z_GAMMA = 12 * _ep(lambda l: l * (l + 1) / (l ** 2 + l - 1))


def nu_zub(r):
    v = 1.0
    for l in sp.factorint(r):
        v *= 1 + l * l / (l ** 4 - 2 * l * l - l + 1)
    return v


def M_zub(y, k=2):
    out = Z_BETA * math.sqrt(y) / (k - 1) - (Z_GAMMA * y if k == 2 else 0.0)
    r = 1
    while r * r < 4 * y:
        out += Z_ALPHA * (-1) ** (k // 2 - 1) / (k - 1) * nu_zub(r) * Ucheb(k - 2, r / (2 * math.sqrt(y))) * math.sqrt(4 * y - r * r)
        r += 1
    return out


# ============================================================ B
def S_(l, a): return (1 - Fr(1, l ** a)) * Fr(l, l - 1)


def local_mean_enclosure(l, J, N, r):
    """exact enclosure [lo, hi] of the Haar mean over p in Z_l^x of F_l(D), D = 4pN - r^2 N^2, from p mod l^J."""
    mod = l ** J
    p = np.arange(mod, dtype=np.int64); p = p[p % l != 0]
    X = (r * r * (N * N % mod) - 4 * p * (N % mod)) % mod              # -D mod l^J
    v = np.zeros(len(X), dtype=np.int64); Y = X.copy(); z = Y == 0
    for _ in range(J):
        sel = (Y % l == 0) & ~z
        v += sel; Y = np.where(sel, Y // l, Y)
    v = np.where(z, J, v)
    if l > 2:
        det = v <= J - 1
        a = v // 2
        leg = np.array([0] + [int(sp.legendre_symbol(t, l)) for t in range(1, l)])
        chi = np.where(v % 2 == 1, 0, leg[Y % l])
        a0 = J // 2
    else:
        det = v <= J - 3
        w8 = Y % 8
        a = np.where(v % 2 == 1, (v - 3) // 2, np.where(w8 % 4 == 1, v // 2, v // 2 - 1))
        chi = np.where((v % 2 == 1) | (w8 % 4 == 3), 0, np.where(w8 == 1, 1, -1))
        a0 = max(0, (J - 4) // 2)
    tot = Fr(0)
    keys, cnt = np.unique(np.stack([a[det], chi[det]]), axis=1, return_counts=True)
    for (aa, cc), c in zip(keys.T, cnt):
        assert aa >= 0
        tot += int(c) * (S_(l, int(aa)) + Fr(1, l ** int(aa)) / (1 - Fr(int(cc), l)))
    nund = int(np.sum(~det)); nall = len(X)
    lo = (tot + nund * S_(l, a0)) / nall
    hi = (tot + nund * Fr(l, l - 1)) / nall
    return lo, hi


LEVELS_B = (19, 23, 9973, 10009, 26, 9998, 21, 9987, 42, 9978, 35, 2310)    # primes in all classes mod 8; 2q; 3q; 6q; 5|N


def sec_B():
    Jl = {2: 18, 3: 11, 5: 7, 7: 6}
    n = 0; bad = []; width = Fr(0)
    for l in (2, 3, 5, 7):
        for N in LEVELS_B:
            for r in range(0, 11):
                lo, hi = local_mean_enclosure(l, Jl[l], N, r)
                pred = Fr(1) if N % l == 0 else e_loc_exact(l, r)
                ok = lo <= pred <= hi
                width = max(width, hi - lo); n += 1
                if not ok:
                    bad.append((l, N, r, float(lo), float(hi), float(pred)))
    check('B1 local means E_l^(N)(r) over p in Z_l^x by exact enumeration of p mod l^J (rigorous enclosure): = 1 for l | N '
          '(l = 2 included), = l^2/(l^2-1) (l | r, r = 0 included) or 1 - u_l (l !| r) for l !| N', not bad,
          '%d cases (l = 2, 3, 5, 7; levels %s; r = 0..10), enclosure width <= %.1e %s' % (n, LEVELS_B, float(width), bad[:2]))
    # B2: the coefficients e_l(i, j) of U6 (l odd) and their sums
    Jj = 10; n = 0; bad = []
    for l in (3, 5, 7):
        I = 3 if l < 7 else 2
        mod = l ** (2 * I + 1)
        p = np.arange(mod, dtype=np.int64); p = p[p % l != 0]; nall = len(p)
        leg = np.array([0] + [int(sp.legendre_symbol(t, l)) for t in range(1, l)])
        for N in (19, 9973, 26, 21, 42, 35, 2310):
            for r in range(0, 7):
                X = (r * r * (N * N % mod) - 4 * p * (N % mod)) % mod          # -D mod l^(2I+1)
                e = {}
                for i in range(I + 1):
                    div = X % l ** (2 * i) == 0
                    xi = leg[(X // l ** (2 * i)) % l]
                    for j in range(Jj + 1):
                        e[i, j] = Fr(int(np.sum(np.where(div, xi ** j if j else 1, 0))), nall)
                ok = True
                if N % l == 0:
                    ok = all(e[i, j] == (1 if i == j == 0 else 0) for i in range(I + 1) for j in range(Jj + 1))
                elif r % l == 0:
                    ok = all(e[0, j] == (1 if j % 2 == 0 else 0) for j in range(Jj + 1)) and \
                        all(e[i, j] == 0 for i in range(1, I + 1) for j in range(Jj + 1))
                else:
                    ok = e[0, 0] == 1 and all(e[0, j] == (Fr(-1, l - 1) if j % 2 else Fr(l - 2, l - 1)) for j in range(1, Jj + 1)) \
                        and all(e[i, 0] == Fr(l, l ** (2 * i) * (l - 1)) and all(abs(e[i, j]) <= e[i, 0] for j in range(Jj + 1))
                                for i in range(1, I + 1))
                s = sum(e[i, j] / l ** (i + j) for i in range(I + 1) for j in range(Jj + 1))
                tail = sum(e[i, 0] for i in range(I + 1)) * Fr(1, l ** (Jj + 1)) * Fr(l, l - 1) \
                    + e[I, 0] * Fr(1, l ** (I + 1)) * Fr(l, l - 1) * Fr(l, l - 1)
                pred = Fr(1) if N % l == 0 else e_loc_exact(l, r)
                ok = ok and abs(s - pred) <= tail
                n += 1
                if not ok:
                    bad.append((l, N, r, float(s), float(pred), float(tail)))
    check('B2 U6 local coefficients e_l(i,j) (exact, l = 3, 5, 7): [i=j=0] for l | N; e(0,odd) = -1/(l-1), e(0,even>=2) = '
          '(l-2)/(l-1), e(i,0) = l/(l^(2i)(l-1)) for l !| 2rN; and sum_{i,j} e_l(i,j) l^(-i-j) = E_l^(N)(r) within a rigorous tail',
          not bad, '%d cases (levels 19, 9973, 26, 21, 42, 35, 2310; r = 0..6) %s' % (n, bad[:2]))


# ============================================================ C
def U_t2(n, t2):
    return {0: Fr(1), 2: 4 * t2 - 1, 4: 16 * t2 * t2 - 12 * t2 + 1}[n]


def trace_rhs(p, N, k):
    s = H(4 * p * N) / 2
    r = 1
    while r * r * N < 4 * p:
        s += (-1) ** (k // 2 - 1) * U_t2(k - 2, Fr(r * r * N, 4 * p)) * H(4 * p * N - r * r * N * N)
        r += 1
    return s - (p + 1 if k == 2 else 0)


LEVELS_C = (23, 46, 69, 78)      # prime, 2q, 3q, 6q


def sec_C():
    for k, pmax in ((2, 500), (4, 250), (6, 120)):
        bad = n = 0
        for N in LEVELS_C:
            c = pari('mfcoefs(mftraceform([%d, %d], 0), %d)' % (N, k, N * pmax))
            sgn = (-1) ** (k // 2 + omega(N))
            for p in sp.primerange(2, pmax + 1):
                if N % p == 0:
                    continue
                lhs = sgn * Fr(str(c[N * p])) / Fr(N * p) ** (k // 2 - 1)      # = sum_f eps_f sqrt(p) lambda_f(p)
                bad += lhs != trace_rhs(p, N, k); n += 1
        check('C1 trace formula sum_f eps_f sqrt(p) lambda_f(p) = H(4pN)/2 + (-1)^(k/2-1) sum_r U_{k-2} H(4pN - r^2N^2) - '
              'delta_{k,2}(p+1) on S_%d^new(N), N = 23, 46, 69, 78, all primes p <= %d, p !| N' % (k, pmax), bad == 0,
              '%d pairs, exact rationals (PARI trace forms of the new space)' % n)
    bad = n = 0
    for N, k in ((46, 2), (78, 2), (69, 2), (46, 4), (69, 4)):
        trs = []
        for fl in (0, 1):
            pari('mf=mfinit([%d,%d],%d); A=mfatkininit(mf,%d)' % (N, k, fl, N))
            trs.append([pari('trace(mfheckemat(mf,%d)*A[2]/A[3])' % p) for p in sp.primerange(5, 80) if N % p])
        bad += trs[0] != trs[1]; n += len(trs[0])
    check('C2 old forms contribute 0: Tr(T_p W_N) on S_k(Gamma_0(N)) = on S_k^new(N), N = 46, 69, 78 (k = 2), 46, 69 (k = 4)',
          bad == 0, '%d traces' % n)
    bad = n = 0
    for N, k in ((46, 2), (69, 2), (78, 2), (46, 4), (78, 4)):
        pari('mf=mfinit([%d,%d],0); B=mfeigenbasis(mf)' % (N, k))
        for i in range(1, int(pari('#B')) + 1):
            if int(pari('poldegree(mffields(mf)[%d])' % i)) != 1:
                continue
            eps = int(pari('lfunrootres(lfunmf(mf,B[%d]))[3]' % i))
            wN = int(pari('mfatkineigenvalues(mf,%d)[%d][1]' % (N, i)))
            bad += eps != (-1) ** (k // 2) * wN
            for q in primes_of(N):
                wq = int(pari('mfatkineigenvalues(mf,%d)[%d][1]' % (q, i)))
                bad += int(pari('mfcoef(B[%d],%d)' % (i, q))) != -wq * q ** (k // 2 - 1)
            n += 1
    check('C3 root number (lfunrootres) = (-1)^(k/2) w_N(f), a_f(q) = -q^(k/2-1) w_q(f) for q | N (rational newforms)',
          bad == 0 and n > 0, '%d newforms' % n)
    zero = [N for N in (6, 10, 22) if int(pari('mfdim([%d,2],0)' % N)) == 0]
    worst = {2: 0.0, 4: 0.0}; bad = 0
    for N in range(2, 3001):
        if not sqfree(N):
            continue
        ph = int(sp.totient(N)); w = omega(N)
        for k in (2, 4):
            dev = abs(int(pari('mfdim([%d,%d],0)' % (N, k))) - Fr((k - 1) * ph, 12))
            worst[k] = max(worst[k], float(dev) / 2 ** w); bad += dev > 2 ** (w + 1)
    dev7 = abs(int(pari('mfdim([701569,2],0)')) - Fr(int(sp.totient(701569)), 12))
    check('C4 dimensions: S_2^new(N) = 0 for N = 6, 10, 22; |dim S_k^new(N) - (k-1) phi(N)/12| <= 2^(omega(N)+1) for squarefree '
          'N <= 3000, k = 2, 4; at N = 701569 = 11*23*47*59 the deviation is 25/3 > 2^(omega-1) = 8',
          zero == [6, 10, 22] and bad == 0 and dev7 == Fr(25, 3),
          'max |dev|/2^omega (N <= 3000): k=2 %.3f, k=4 %.3f; N = 701569: %s' % (worst[2], worst[4], dev7))


# ============================================================ D
TYPES = {'q': (), '2q': (2,), '3q': (3,), '6q': (2, 3)}


def sec_D():
    tail = 1.0 / (10 ** 6) ** 2          # sum_{l > 1e6} u_l <= sum_{n > 1e6} 2/n^3 <= 1e-12 (u_l <= 2/l^3 for l >= 3)
    check('D1 C0 = prod_l (1 - u_l) = nu\'(1) = %.9f (primes < 1e6; tail < %.0e); nu\'(r) = C0 prod_{l|r} (l^2/(l^2-1))/(1-u_l)'
          % (C0, tail), abs(C0 - 0.615132657) < 1e-8 and abs(nu_prime(1) - C0) < 1e-15)
    ok = all(abs(12 * kappa(S) * E0(S) / math.pi - 2 * math.pi * piN(S)) < 1e-12 for S in
             [(), (2,), (3,), (2, 3), (5,), (2, 3, 5, 7, 11)])
    check('D2 (12/pi) kappa(N) E^(N)(0) = 2 pi prod_{q|N}(1 + 1/q), E^(N)(0) = zeta(2) prod_{q|N}(1 - q^-2)', ok)
    q = 1000003
    ys = np.linspace(0.01, 6, 300)
    dev = max(abs(M_level(y, (q,), k) - M_level(y, (), k)) for k in (2, 4) for y in ys)
    check('D3 special case N = q prime: max_{y<=6} |M^(q) - M_k| <= 10/q at q = 1000003 (k = 2, 4)', dev <= 10.0 / q, '%.2e' % dev)

    def closed_2q(y):
        c2 = lambda r: 1.5 if r % 2 else 0.75
        s = 3 * math.pi * math.sqrt(y) - 24 * y
        r = 1
        while r * r < 4 * y:
            s += (24 / math.pi) * nu_prime(r) * c2(r) * math.sqrt(4 * y - r * r); r += 1
        return s
    dev = max(abs(M_level(y, (2, q)) - closed_2q(y)) for y in ys)
    check('D4 special case N = 2q: M^(2q) = 3 pi sqrt y - 24 y + (24/pi) sum_r nu\'(r) c_2(r) sqrt(4y - r^2) + O(1/q), c_2 = 3/2 '
          '(r odd), 3/4 (r even); max dev at q = 1000003, y <= 6', dev <= 20.0 / q, '%.2e' % dev)
    yt = (0.1, 0.2, 0.5, 1.0, 2.0)
    print('     table k = 2 (q -> infinity):  y = ' + ' '.join('%8.2f' % y for y in yt))
    for t, S in TYPES.items():
        print('     %-4s ' % t + ' '.join('%8.4f' % M_level(y, S) for y in yt))
    print('     zub  ' + ' '.join('%8.4f' % M_zub(y) for y in yt))
    print('     table k = 4:')
    for t, S in TYPES.items():
        print('     %-4s ' % t + ' '.join('%8.4f' % M_level(y, S, 4) for y in yt))
    print('     zub  ' + ' '.join('%8.4f' % M_zub(y, 4) for y in yt))
    for t, S in TYPES.items():
        print('     %-3s kappa %.4f  prod(1+1/q) %.4f  kappa*nu_N(1..4) %s' % (t, kappa(S), piN(S),
              ' '.join('%.5f' % (kappa(S) * nu_N(r, S)) for r in (1, 2, 3, 4))))
    # Proposition: exact Euler factors, prime by prime
    bad = 0
    for l in sp.primerange(2, 1000):
        L = Fr(l)
        avg = lambda e: (L * e + 1) / (L + 1)            # P(l | N) = 1/(l+1) among squarefree N; factor 1 at l | N
        e0 = e_loc_exact(l, 0); e1 = e_loc_exact(l, 1)
        phi_avg = (L + 1 - 1 / L) / (L + 1); e0_avg = (L + 1 - 1 / L ** 2) / (L + 1)
        bad += avg(e0) != 1 + L / ((L + 1) ** 2 * (L - 1))
        bad += avg(e1) != (L ** 4 - 2 * L ** 2 - L + 1) / (L ** 2 - 1) ** 2
        bad += avg(e0) / avg(e1) != 1 + L ** 2 / (L ** 4 - 2 * L ** 2 - L + 1)
        bad += e0_avg / phi_avg != (L ** 3 + L ** 2 - 1) / (L * (L ** 2 + L - 1))
        bad += 1 / phi_avg != L * (L + 1) / (L ** 2 + L - 1)
        bad += (1 - L ** -2) * e1 * (L + 1 / e1) / (L + 1) / phi_avg != (L ** 4 - 2 * L ** 2 - L + 1) / (L ** 4 - 2 * L ** 2 + L)
    check('D5 Proposition (average): per prime l < 1000, exactly: mean factors (l e_l + 1)/(l + 1) = her A, B; ratio = her nu(r); '
          'beta/2pi, gamma/12, alpha/2pi factors', bad == 0,
          'alpha = %.6f, beta = %.6f, gamma = %.6f' % (Z_ALPHA, Z_BETA, Z_GAMMA))
    X = 300000
    sq = np.ones(X + 1, bool); sq[0] = False
    for l in sp.primerange(2, math.isqrt(X) + 1):
        sq[l * l::l * l] = False
    fphi = np.ones(X + 1); f0 = np.ones(X + 1); fr = {r: np.ones(X + 1) for r in (1, 2, 3)}
    for l in sp.primerange(2, X + 1):
        fphi[l::l] *= 1 - 1.0 / l; f0[l::l] *= 1 - 1.0 / l ** 2
        for r in fr:
            fr[r][l::l] /= e_loc(l, r)
    Ns = np.arange(X + 1, dtype=float)
    phiN = (Ns * fphi)[sq]
    worst = 0.0; row = []
    for y in (0.1, 0.5, 1.0, 2.0, 3.3):
        tot = Ns[sq] * (2 * math.pi * math.sqrt(y) * f0[sq] - 12 * y)
        r = 1
        while r * r < 4 * y:
            tot = tot + Ns[sq] * (12 / math.pi) * nu_prime(r) * fr[r][sq] * math.sqrt(4 * y - r * r); r += 1
        avg = tot.sum() / phiN.sum()
        worst = max(worst, abs(avg - M_zub(y))); row.append('%.1f: %.4f/%.4f' % (y, avg, M_zub(y)))
    check('D5 sum_{N<=3e5 squarefree} phi(N) M^(N)(y) / sum phi(N) = M_zub(y), y = 0.1, 0.5, 1, 2, 3.3', worst < 5e-3, '; '.join(row))
    # repaired comparison: S = primes of N up to 2 sqrt(Y+1) (Y = 3: S subset {2, 3}); eps = sum_{q | N, q not in S} 1/q
    bad = 0; emax = 0.0
    spf = np.zeros(200001, dtype=np.int64)
    for l in sp.primerange(2, 200001):
        blk = spf[l::l]; blk[blk == 0] = l
    for N in range(2, 200001):
        if not sq[N]:
            continue
        P = []; x = N
        while x > 1:
            P.append(int(spf[x])); x //= int(spf[x])
        P = tuple(P); S = tuple(l for l in P if l <= 4); R = [l for l in P if l > 4]
        eps = sum(1.0 / l for l in R); emax = max(emax, eps)
        rk = kappa(P) / kappa(S); rp = piN(P) / piN(S)
        bad += not (1 <= rk <= math.exp(1.5 * eps) + 1e-12 and 1 <= rp <= math.exp(eps) + 1e-12)
        for r in (1, 2, 3):
            rn = nu_N(r, P) / nu_N(r, S)
            bad += not (1 - 1e-12 <= rn <= math.exp(eps) + 1e-12)
    epsP = sum(1.0 / l for l in sp.primerange(5, 128))
    check('D6 (repair "eps <= 1") kappa(N)/kappa(S) in [1, e^(1.5 eps)], prod(1+1/q) ratio in [1, e^eps], nu_N/nu_S in [1, e^eps], '
          'squarefree N <= 2e5, S = prime factors <= 4', bad == 0,
          'max eps here %.3f; eps > 1 occurs, e.g. N = prod of primes 5..127: eps = %.3f' % (emax, epsP))


# ============================================================ E
GP = ('lvl1(N,Y)=my(v=List());forprime(p=2,Y*N,if(N%p,my(h=vector(4),r=0,s=0);while(r^2*N<4*p,'
      'my(c=qfbhclassno(4*p*N-r^2*N^2));if(r<4,h[r+1]=c);s+=if(r,c,c/2);r++);'
      'listput(v,[p,h[1],h[2],h[3],h[4],s-(p+1)])));Vec(v)')
EDGES = np.arange(0, 4.0 + 1e-9, 0.25)
FRESH = {'q': 9973, '2q': 2 * 4999, '3q': 3 * 3329, '6q': 6 * 1663}
EARLIER = set([10007, 10009, 10037, 10039, 10061] + [2 * q for q in (5003, 5009, 5011, 5021, 5023)] +
              [3 * q for q in (3331, 3343, 3347, 3359, 3361)] + [6 * q for q in (1669, 1693, 1697, 1699, 1709)] +
              [30011, 30029, 20014, 3 * 6673, 6 * 3343, 5 * 4003, 20011])


def sec_E():
    pari(GP)
    okf = all(int(pari('isprime(%d)' % (N // m))) for N, m in zip(FRESH.values(), (1, 2, 3, 6))) and \
        not set(FRESH.values()) & EARLIER and not any(N in FORBIDDEN for N in FRESH.values())
    check('E0 levels 9973, 2*4999, 3*3329, 6*1663: of the four types, not used in the earlier notes or by their referees, '
          'outside 30030..30530', okf)
    res = {}; zmax = 0.0
    for t, N in FRESH.items():
        S = primes_of(N)
        rows = np.array([[float(x) for x in row] for row in pari('lvl1(%d,4)' % N)])
        ps = rows[:, 0]; y = ps / N; w = np.log(ps)
        dim = int(pari('mfdim([%d,2],0)' % N))
        E = rows[:, 5] / dim

        def own(yy, p):
            s = E0(S) * math.sqrt(yy); r = 1
            while r * r < 4 * yy:
                s += nu_N(r, S) * math.sqrt(4 * yy - r * r); r += 1
            return (N / math.pi * s - (p + 1)) / dim
        mod = {'own': np.array([own(a, b) for a, b in zip(y, ps)]), 'zub': np.array([M_zub(a) for a in y])}
        for t2, S2 in TYPES.items():
            if t2 != t:
                mod['type_' + t2] = np.array([M_level(a, S2) for a in y])
        chi = {k_: 0.0 for k_ in mod}; data = []
        for a, b in zip(EDGES[:-1], EDGES[1:]):
            sel = (y > a) & (y <= b); ww = w[sel]
            dm = np.sum(ww * E[sel]) / np.sum(ww)
            resid = E[sel] - mod['own'][sel]; rm = np.sum(ww * resid) / np.sum(ww)
            se = math.sqrt(np.sum(ww ** 2 * (resid - rm) ** 2)) / np.sum(ww)
            data.append(dm)
            for k_, v in mod.items():
                chi[k_] += ((dm - np.sum(ww * v[sel]) / np.sum(ww)) / se) ** 2
        res[t] = chi
        alt = {k_: v for k_, v in chi.items() if k_ != 'own'}
        ok = chi['own'] <= 32 and min(alt.values()) >= chi['own'] + 9
        check('E1 level %d (%s, dim %d, %d primes p <= 4N): chi2 on 16 windows (single-value s.e.) own %.1f; others %s'
              % (N, t, dim, len(ps), chi['own'], ', '.join('%s %.1f' % kv for kv in alt.items())), ok)
        print('     data ' + ' '.join('%.3f' % v for v in data))
        # class-number means = local predictions (key lemma main term at this level)
        line = []
        for r in range(4):
            Dr = 4 * ps * N - r * r * N * N; sel = Dr > 0
            Ls = math.pi * rows[sel, 1 + r] / np.sqrt(Dr[sel]); ww = w[sel]
            m = np.sum(ww * Ls) / np.sum(ww); se = math.sqrt(np.sum(ww ** 2 * (Ls - m) ** 2)) / np.sum(ww)
            pr = E0(S) if r == 0 else nu_N(r, S)
            z = (m - pr) / se; zmax = max(zmax, abs(z)); line.append('r=%d %.4f/%.4f (z %+.2f)' % (r, m, pr, z))
        print('     %s mean L*(D_r) / prediction: ' % t + '; '.join(line))
    check('E2 log-weighted means over p <= 4N of L*(D_r(p)) = pi H(D_r)/sqrt(D_r) equal E^(N)(0), nu_N(1..3) (|z| < 4)', zmax < 4,
          'max |z| %.2f' % zmax)
    # E3: the truncation of U2 at a composite level: per p, L(1, xi) cut at n <= M
    for N in (FRESH['2q'], FRESH['6q']):
        S = primes_of(N)
        Mc = {'N^(1/2)': int(N ** 0.5), 'N^(3/4)': int(N ** 0.75), 'N^(7/8)': int(N ** 0.875)}
        Mmax = max(Mc.values())
        ps = [int(p) for p in sp.primerange(N + 1, 2 * N)]
        msg = []
        okN = True
        for r in (0, 1):
            ex = []; tr = {k_: [] for k_ in Mc}
            for p in ps:
                D = 4 * p * N - r * r * N * N
                ex.append(math.pi * float(H(D)) / math.sqrt(D))
                acc = {k_: 0.0 for k_ in Mc}
                for m in Mset(D):
                    Dm = -D // (m * m)
                    v = np.cumsum([gmpy2.kronecker(Dm, n) / n for n in range(1, Mmax + 1)])
                    for k_, Mv in Mc.items():
                        acc[k_] += v[Mv - 1] / m
                for k_ in Mc:
                    tr[k_].append(acc[k_])
            ex = np.array(ex); rms = {k_: float(np.sqrt(np.mean((np.array(tr[k_]) - ex) ** 2))) for k_ in Mc}
            pr = E0(S) if r == 0 else nu_N(r, S)
            t78 = np.array(tr['N^(7/8)'])
            z = (t78.mean() - pr) / (t78.std() / math.sqrt(len(t78)))
            okN &= rms['N^(1/2)'] > rms['N^(3/4)'] > rms['N^(7/8)'] and abs(z) < 4
            msg.append('r=%d: rms error %s; mean of the series cut at N^(7/8) %.4f vs %.4f (z %+.2f)'
                       % (r, ', '.join('%s %.4f' % kv for kv in rms.items()), t78.mean(), pr, z))
        check('E3 U2 at the level %d: truncation error of sum_m m^-1 sum_{n<=M} xi(n)/n falls with M, and the mean over the %d primes '
              'in (N, 2N] of the series cut at N^(7/8) matches the local prediction' % (N, len(ps)), okN, ' | '.join(msg))


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'ABCDE'
    for ch in which:
        t = time.time()
        globals()['sec_' + ch]()
        print('     [section %s: %.0f s]' % (ch, time.time() - t), flush=True)
    print('run time %.0f s' % (time.time() - T0))
    print('%d checks, %d failed' % (NCHK, NFAIL))
