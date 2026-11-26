"""Checks for mprime_theorem.md (own code; PARI via cypari2 only for eigenform traces, dimensions and, after a
brute-force cross-check, Hurwitz class numbers of large arguments).

Sections (labels as in the note):
 A  Kloosterman sums, brute force: S(Nn,1;Nm) = -S(n,Nbar;m); S(Na,1;N^2 m') = 0; sum_a S(la,b;m) = mu(m) c_m(l).
 B  main-term algebra of step (c): mu(dg) c_dg(s^2 g^2)/phi(dg) = mu(g)/phi(d); A_nat(s,d) = (6/pi^2) s^-2 w(d).
 C  Euler-product identities behind nu'(r) and zeta(2).
 D  Poisson step (d): sum_s J_{k-1}(xs)/s in closed form (k = 2, 4); Moebius regrouping; closed form = Bessel form.
 E  class numbers: brute-force H(D) = PARI; H(D) = h'(D0) f0 prod Q_l; local means E_l(r) (exact enumeration).
 F  trace formula of T_p W_N at prime level against PARI trace forms (k = 2 and k = 4); dimensions.
 G  one prime level, no average over N: window means of eps_f a_f(p) against M_k (k = 2, 4); class-number means.
 H  the blind data file of round 160 against this note's M_2.
"""
import sys, math, cmath
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
from scipy.special import jv
import cypari2

pari = cypari2.Pari()
pari.allocatemem(300 * 10 ** 6)
DATA = '/home/bas/praime/writeup/data/study/r160/data_prime_hold.npz'
NCHK = 0
NFAIL = 0


def check(name, ok, info=''):
    global NCHK, NFAIL
    NCHK += 1
    NFAIL += (not ok)
    print(('PASS ' if ok else 'FAIL ') + name + ('  ' + info if info else ''), flush=True)


# ---------------------------------------------------------------- arithmetic helpers
def mu(n): return int(sp.mobius(n))
def phi(n): return int(sp.totient(n))
def ramanujan(m, n): return sum(mu(m // d) * d for d in sp.divisors(math.gcd(m, n)))
def sqfree(n): return all(e == 1 for e in sp.factorint(n).values())


def kloost(a, b, c):
    s = 0j
    for x in range(1, c + 1):
        if math.gcd(x, c) == 1:
            s += cmath.exp(2j * math.pi * (a * x + b * pow(x, -1, c)) / c)
    return s


# ---------------------------------------------------------------- A
def sec_A():
    worst = 0.0
    for N in (5, 7, 11):
        for m in range(1, 26):
            if m % N == 0:
                continue
            Nb = pow(N, -1, m) if m > 1 else 0
            for n in range(1, 13):
                worst = max(worst, abs(kloost(N * n, 1, N * m) + kloost(n, Nb, m)))
    check('A1 S(Nn,1;Nm) = -S(n,Nbar;m), N in {5,7,11}, m <= 25, n <= 12', worst < 1e-8, 'max dev %.1e' % worst)
    worst = 0.0
    for N in (5, 7):
        for mp in (1, 2, 3, 4, 6):
            for a in range(1, 7):
                if a % N:
                    worst = max(worst, abs(kloost(N * a, 1, N * N * mp)))
    check('A2 S(Na,1;N^2 m\') = 0 for N not dividing a', worst < 1e-8, 'max |S| %.1e' % worst)
    worst = 0.0
    for m in range(1, 31):
        for l in range(1, 21):
            for b in (1, 2, 5, 7):
                if math.gcd(b, m) != 1:
                    continue
                s = sum(kloost(l * a, b, m) for a in range(1, m + 1) if math.gcd(a, m) == 1)
                worst = max(worst, abs(s - mu(m) * ramanujan(m, l)))
    check('A3 sum_{a mod m, (a,m)=1} S(la,b;m) = mu(m) c_m(l), m <= 30, l <= 20', worst < 1e-7, 'max dev %.1e' % worst)
    worst = 0.0
    for c in range(1, 121):
        tc = len(sp.divisors(c))
        for x in range(1, 9):
            for y in range(1, 9):
                g = math.gcd(math.gcd(x, y), c)
                worst = max(worst, abs(kloost(x, y, c)) / (tc * math.sqrt(g * c)))
    check('A4 Weil bound |S(a,b;c)| <= tau(c) (a,b,c)^(1/2) c^(1/2), c <= 120, a, b <= 8', worst <= 1 + 1e-9, 'max ratio %.3f' % worst)


# ---------------------------------------------------------------- B
def w_sq(d):
    v = Fr(1)
    for p in sp.factorint(d):
        v *= Fr(p * p, (p - 1) ** 2 * (p + 1))
    return v


def sec_B():
    bad = 0; n = 0
    for s in range(1, 13):
        for d in range(1, 31):
            if not sqfree(d) or math.gcd(s, d) != 1:
                continue
            for g in range(1, 31):
                if not sqfree(g) or math.gcd(g, d) != 1:
                    continue
                m = d * g
                lhs = Fr(mu(m) * ramanujan(m, s * s * g * g), phi(m))
                bad += (lhs != Fr(mu(g), phi(d))); n += 1
    check('B1 mu(dg) c_dg(s^2 g^2)/phi(dg) = mu(g)/phi(d) (exact)', bad == 0, '%d cases' % n)
    worst = 0.0
    G = 3000
    for s, d in ((1, 1), (2, 1), (3, 1), (1, 2), (1, 3), (1, 6), (2, 3), (3, 2), (5, 6), (4, 1), (1, 5), (7, 10), (1, 30)):
        tot = 0.0
        for g in range(1, G + 1):
            t, m = s * g, d * g
            mm = mu(m)
            if mm:
                tot += mm * ramanujan(m, t * t) / (phi(m) * t * t)
        pred = 6 / math.pi ** 2 / s ** 2 * float(w_sq(d))
        worst = max(worst, abs(tot / pred - 1))
    check('B2 A_nat(s,d) = sum_{t/m=s/d} mu(m)c_m(t^2)/(phi(m)t^2) = (6/pi^2) s^-2 w(d), 13 pairs, g <= 3000',
          worst < 2e-3, 'max rel dev %.1e (tail ~ 1/G)' % worst)
    # d not squarefree -> 0
    tot = sum(mu(4 * g) * ramanujan(4 * g, g * g) / (phi(4 * g) * g * g) for g in range(1, 500) if mu(4 * g))
    check('B3 A_nat(1,4) = 0 (d not squarefree)', abs(tot) < 1e-15)


# ---------------------------------------------------------------- C
PR6 = list(sp.primerange(2, 10 ** 6))


def u_(p): return 1.0 / ((p - 1) ** 2 * (p + 1))


C0 = math.exp(sum(math.log1p(-u_(p)) for p in PR6))   # tail p > 1e6: < 1e-12


def nu_prime_def(r):
    """nu'(r) = sum_{delta | r, delta sqfree} c(delta), c(delta) = (C0/delta) prod_{p|delta} w(p)/(1-u_p)."""
    fr = list(sp.factorint(r))
    tot = 0.0
    for k in range(1 << len(fr)):
        sub = [q for i, q in enumerate(fr) if k >> i & 1]
        v = C0
        for q in sub:
            v *= float(w_sq(q)) / (1 - u_(q)) / q
        tot += v
    return tot


def nu_prime_prod(r):
    """Euler product form: prod_{l not| r} (1 - u_l) prod_{l | r} l^2/(l^2 - 1)."""
    v = C0
    for q in sp.factorint(r):
        v *= (q * q / (q * q - 1)) / (1 - u_(q))
    return v


def sec_C():
    worst = max(abs(nu_prime_def(r) / nu_prime_prod(r) - 1) for r in range(1, 61))
    check('C1 nu\'(r) (divisor-sum definition) = prod_{l!|r}(1-u_l) prod_{l|r} l^2/(l^2-1), r <= 60', worst < 1e-12,
          'max rel dev %.1e; C0 = %.9f, nu\'(1) = %.9f' % (worst, C0, nu_prime_prod(1)))
    bad = 0
    for p in sp.primerange(2, 1000):
        P = sp.Integer(p)
        wp = P ** 2 / ((P - 1) ** 2 * (P + 1)); up = 1 / ((P - 1) ** 2 * (P + 1))
        bad += sp.simplify(1 + wp * (P - 1) / P ** 2 - P ** 2 / (P ** 2 - 1)) != 0      # sum_d w(d)phi(d)/d^2 = zeta(2)
        bad += sp.simplify(1 - wp / P ** 2 - (1 - up)) != 0                           # sum_{(e,delta)=1} mu(e)w(e)/e^2
        bad += sp.simplify(1 - up + wp / P - P ** 2 / (P ** 2 - 1)) != 0                # l | r factor of nu'
        bad += sp.simplify(P * (P ** 2 - P - 1) / ((P - 1) ** 2 * (P + 1)) - (1 - up)) != 0  # local mean, l !| r
    check('C2 local identities (zeta(2), c(delta), nu\' factors, local mean) for all p < 1000 (exact)', bad == 0)
    z2 = math.exp(sum(-math.log1p(-1.0 / p ** 2) for p in PR6))
    check('C3 prod_p p^2/(p^2-1) = zeta(2) = pi^2/6', abs(z2 - math.pi ** 2 / 6) < 1e-6, '%.9f' % z2)


# ---------------------------------------------------------------- D
def U(n, t):
    if n == 0: return 1.0 + 0 * t
    a, b = 1.0 + 0 * t, 2 * t
    for _ in range(n - 1):
        a, b = b, 2 * t * b - a
    return b


def Fk_closed(x, k):
    """sum_{s>=1} J_{k-1}(x s)/s by Poisson (k even)."""
    nu = k - 1
    out = 1.0 / nu
    n = 1
    while 2 * math.pi * n < x:
        t = 2 * math.pi * n / x
        out += (2.0 / nu) * (-1) ** (k // 2 - 1) * math.sqrt(1 - t * t) * U(k - 2, t)
        n += 1
    if k == 2:
        out -= x / 4
    return out


def Fk_direct(x, k, S=400000):
    s = np.arange(1, S + 1, dtype=float)
    return float(np.sum(jv(k - 1, x * s) / s))


def M_closed(y, k=2):
    y = float(y)
    out = 2 * math.pi * math.sqrt(y)
    r = 1
    while r * r < 4 * y:
        out += (-1) ** (k // 2 - 1) * (12 / math.pi) * nu_prime_prod(r) * U(k - 2, r / (2 * math.sqrt(y))) * math.sqrt(4 * y - r * r)
        r += 1
    out /= (k - 1)
    if k == 2:
        out -= 12 * y
    return out


def _ep(fn):
    return math.exp(sum(math.log(fn(float(p))) for p in PR6))


ZA = 2 * math.pi * _ep(lambda p: (1 - p - 2 * p ** 2 + p ** 4) / (p ** 4 - 2 * p ** 2 + p))
ZB = 2 * math.pi * _ep(lambda p: (-1 + p ** 2 + p ** 3) / (p * (-1 + p + p ** 2)))
ZG = 12 * _ep(lambda p: p * (1 + p) / (-1 + p + p ** 2))


def M_zub(y, k=2):
    """Zubrilina, arXiv 2310.07681, Theorem 1 (squarefree levels averaged over N), as an alternative."""
    out = ZB * math.sqrt(y) / (k - 1) - (ZG * y if k == 2 else 0.0)
    r = 1
    while r * r < 4 * y:
        nu = 1.0
        for q in sp.factorint(r):
            nu *= 1 + q * q / (q ** 4 - 2 * q * q - q + 1)
        out += ZA * (-1) ** (k // 2 - 1) / (k - 1) * nu * math.sqrt(4 * y - r * r) * U(k - 2, r / (2 * math.sqrt(y)))
        r += 1
    return out



def sec_D():
    x = np.linspace(1e-6, 300, 600001)
    worst = 0.0
    for nu in (1, 3, 5, 9, 13):
        J = jv(nu, x); dJ = 0.5 * (jv(nu - 1, x) - jv(nu + 1, x))
        bound = np.minimum(1.0, (x / 2) ** nu / math.factorial(nu))
        worst = max(worst, float(np.max(np.abs(J) / bound)), float(np.max(np.abs(J) / (x / 2))), float(np.max(np.abs(dJ))))
    check('D0 |J_nu(x)| <= min(1, (x/2)^nu/nu!) <= x/2 and |J_nu\'(x)| <= 1 (nu = 1,3,5,9,13; 0 < x <= 300)', worst <= 1 + 1e-12,
          'max ratio %.4f' % worst)
    worst = 0.0
    for k in (2, 4):
        for x in (0.5, 3.0, 6.0, 7.0, 13.3, 25.1, 40.7):
            worst = max(worst, abs(Fk_direct(x, k) - Fk_closed(x, k)))
    check('D1 sum_s J_{k-1}(xs)/s = Poisson closed form (k = 2, 4; 7 values of x; s <= 4e5)', worst < 2e-5, 'max dev %.1e' % worst)
    # Moebius + Poisson per d: sum_{(s,d)=1} J_1(cs)/s = sum_{e|d} mu(e)/e F(ce)
    worst = 0.0
    for y in (0.3, 1.3, 3.7):
        for d in (1, 2, 3, 5, 6, 10, 15, 30):
            c = 4 * math.pi * math.sqrt(y) / d
            s = np.arange(1, 400001); s = s[np.gcd(s, d) == 1].astype(float)
            direct = float(np.sum(jv(1, c * s) / s))
            reg = sum(mu(e) / e * Fk_closed(c * e, 2) for e in sp.divisors(d))
            worst = max(worst, abs(direct - reg))
    check('D2 sum_{(s,d)=1} J_1(cs)/s = sum_{e|d} mu(e)/e F_0(ce) (Moebius + Poisson), 24 cases', worst < 2e-5, 'max dev %.1e' % worst)
    # full regrouping: Bessel form (12/pi) sqrt(y) sum_d w(d)/d sum_{e|d} mu(e)/e F_0(4 pi sqrt(y) e/d) = closed form
    D = 20000
    sq = [d for d in range(1, D + 1) if sqfree(d)]
    wd = {d: float(w_sq(d)) for d in sq}
    divs = {d: sp.divisors(d) for d in sq}
    worst = 0.0
    for y in (0.3, 1.3, 3.7, 7.9):
        tot = 0.0
        for d in sq:
            tot += wd[d] / d * sum(mu(e) / e * Fk_closed(4 * math.pi * math.sqrt(y) * e / d, 2) for e in divs[d])
        bes = 12 / math.pi * math.sqrt(y) * tot
        worst = max(worst, abs(bes - M_closed(y)))
    check('D3 Bessel form (d <= 2e4, s-sums exact) = closed form M_prime at y = 0.3, 1.3, 3.7, 7.9', worst < 3e-3, 'max dev %.1e' % worst)


# ---------------------------------------------------------------- E
def hurwitz_brute(D):
    """Hurwitz class number H(D): reduced forms of discriminant -D, weights 1/2 (a(x^2+y^2)) and 1/3 (a(x^2+xy+y^2))."""
    if D % 4 not in (0, 3):
        return Fr(0)
    tot = Fr(0)
    a = 1
    while 3 * a * a <= D:
        for b in range(-a + 1, a + 1):
            if (b * b + D) % (4 * a):
                continue
            c = (b * b + D) // (4 * a)
            if c < a or (c == a and b < 0):
                continue
            if a == c and b == 0:
                tot += Fr(1, 2)
            elif a == b == c:
                tot += Fr(1, 3)
            else:
                tot += 1
        a += 1
    return tot


def is_fund(Dl):  # Dl < 0
    if Dl % 4 == 1:
        return sqfree(-Dl)
    if Dl % 4 == 0:
        m = Dl // 4
        return m % 4 in (2, 3) and sqfree(-m)
    return False


def fund_part(D):
    """-D = D0 f0^2 with D0 fundamental."""
    for f in range(math.isqrt(D), 0, -1):
        if D % (f * f) == 0 and is_fund(-D // (f * f)):
            return -D // (f * f), f
    raise ValueError(D)


def kron(D0, l):
    if l == 2:
        if D0 % 2 == 0: return 0
        return 1 if D0 % 8 == 1 else -1
    return int(sp.jacobi_symbol(D0 % l, l)) if D0 % l else 0


def sec_E():
    bad = 0; n = 0
    for D in range(3, 3001):
        if D % 4 in (0, 3):
            bad += hurwitz_brute(D) != Fr(str(pari('qfbhclassno(%d)' % D))); n += 1
    check('E1 brute-force Hurwitz H(D) = PARI qfbhclassno, D <= 3000', bad == 0, '%d values' % n)
    bad = 0; n = 0
    for D in range(3, 3001):
        if D % 4 not in (0, 3):
            continue
        D0, f0 = fund_part(D)
        hp = hurwitz_brute(-D0)        # = h'(D0) for D0 fundamental (only primitive forms exist)
        pred = hp * f0
        for l, a in sp.factorint(f0).items():
            chi = kron(D0, l)
            pred *= (1 - Fr(chi, l)) * sum(Fr(1, l ** j) for j in range(a)) + Fr(1, l ** a)
        bad += pred != hurwitz_brute(D); n += 1
    check('E2 H(D) = h\'(D0) f0 prod_{l^a||f0} [(1 - chi(l)/l) sum_{j<a} l^-j + l^-a] (exact)', bad == 0, '%d values' % n)
    # local means E_l(r) by exact enumeration over p mod l^J (units), D = 4pN - r^2 N^2
    J = {2: 18, 3: 11, 5: 7, 7: 6}
    worst = 0.0; nc = 0; okall = True
    for l in (2, 3, 5, 7):
        q = l ** J[l]
        p = np.arange(q, dtype=np.int64); p = p[p % l != 0]
        for N in (10007, 100003, 13, l if l > 2 else 3):
            for r in range(0, 13):
                Dl = (r * r * N * N - 4 * p * N) % q        # Delta = -D mod l^J
                v = np.zeros(len(p), dtype=np.int64)
                for kk in range(1, J[l] + 1):
                    v += (Dl % l ** kk == 0)
                lo, hi = local_F(Dl, v, l, J[l])
                Elo, Ehi = lo.mean(), hi.mean()
                if l == N:
                    pred = 1.0
                elif l == 2:
                    pred = 2 / 3 if r % 2 else 4 / 3
                elif r % l == 0:
                    pred = l * l / (l * l - 1)
                else:
                    pred = 1 - 1 / ((l - 1) ** 2 * (l + 1))
                ok = (Elo - 1e-12 <= pred <= Ehi + 1e-12) and (Ehi - Elo < 2e-3)
                okall &= ok; nc += 1
                worst = max(worst, Ehi - Elo)
    check('E3 local means E_l(r) over p in Z_l^x: 2/3 | 4/3 (l=2, r odd|even), l^2/(l^2-1) (l|r), 1-u_l (l!|r), 1 (l=N)',
          okall, '%d cases (l = 2,3,5,7; 4 levels; r = 0..12), enclosure width <= %.1e' % (nc, worst))


def local_F(Dl, v, l, J):
    """lower/upper value of F_l(D) given Delta = -D mod l^J (v = valuation, capped at J), vectorized.
    F_l = sum_{j<a} l^-j + l^-a/(1 - chi(l)/l), -D = D0 f0^2, a = v_l(f0), chi = Kronecker symbol of D0 at l."""
    x = 1.0 / l
    S = lambda a: (1 - x ** a) / (1 - x)
    vv = np.minimum(v, J)
    pw = np.array([l ** j for j in range(J + 1)], dtype=np.int64)
    unit = (Dl // pw[vv]) if l > 2 else (Dl >> vv)
    if l > 2:
        und = v >= J
        leg = np.array([0] + [int(sp.legendre_symbol(u, l)) for u in range(1, l)])
        a = vv // 2
        chi = np.where(vv % 2 == 1, 0, leg[unit % l])
        a0 = J // 2
    else:
        und = v >= J - 3
        u = unit % 8
        odd_v = vv % 2 == 1
        a = np.where(odd_v, (vv - 3) // 2, np.where(u % 4 == 1, vv // 2, vv // 2 - 1))
        chi = np.where(odd_v, 0, np.where(u % 4 == 1, np.where(u == 1, 1, -1), 0))
        a0 = max((J - 6) // 2, 0)
    a = np.maximum(a, 0)
    F = S(a) + x ** a / (1 - chi * x)
    lo = np.where(und, S(a0) - x ** a0 / (1 - x), F)
    hi = np.where(und, S(a0) + x ** a0 / (1 - x), F)
    return lo, hi


# ---------------------------------------------------------------- F
def H(D): return Fr(str(pari('qfbhclassno(%d)' % D)))


def trace_rhs(p, N, k):
    s = H(4 * p * N) / 2
    r = 1
    while r * r * N < 4 * p:
        Ur = 1 if k == 2 else None
        t = Fr(r * r * N, 4 * p)          # t^2 with t = r sqrt(N)/(2 sqrt(p))
        if k == 4:
            Ur = 4 * t - 1                # U_2(t) = 4t^2 - 1
        s += (-1) ** (k // 2 - 1) * Ur * H(4 * p * N - r * r * N * N)
        r += 1
    if k == 2:
        s -= p + 1
    return s


def sec_F():
    for k, levels, pmax in ((2, (11, 37, 43, 61, 101, 163), 900), (4, (11, 13, 23, 61), 500)):
        bad = n = 0
        for N in levels:
            c = pari('mfcoefs(mftraceform([%d, %d], 0), %d)' % (N, k, N * pmax))
            for p in sp.primerange(2, pmax + 1):
                if p == N:
                    continue
                cNp = Fr(str(c[N * p]))
                lhs = cNp if k == 2 else -cNp / (N * p)          # sum_f eps_f sqrt(p) lambda_f(p)
                bad += lhs != trace_rhs(p, N, k); n += 1
        check('F%d trace formula (k = %d) = sum_f eps_f sqrt(p) lambda_f(p) from PARI trace forms, levels %s' % (k // 2, k, levels),
              bad == 0, '%d pairs (N, p)' % n)
    worst = 0.0; bad = 0
    for N in sp.primerange(5, 2000):
        dm = int(pari('mfdim([%d, 2], 0)' % N))
        nu2 = 1 + int(sp.jacobi_symbol(-1 % N, N)); nu3 = 1 + int(sp.jacobi_symbol(-3 % N, N))
        bad += Fr(N + 1, 12) - Fr(nu2, 4) - Fr(nu3, 3) != dm
        worst = max(worst, abs(12 * dm / N - 1) * N)
    check('F3 dim S_2(Gamma_0(N)) = (N+1)/12 - nu_2/4 - nu_3/3 for primes 5 <= N < 2000; |12 dim/N - 1| <= 14/N', bad == 0 and worst <= 14,
          'max N|12dim/N-1| = %.2f' % worst)


# ---------------------------------------------------------------- G
def sec_G():
    for N, ymax, h in ((10007, 12.0, 0.5), (100003, 4.0, 0.25)):
        ps = [int(q) for q in pari('primes([3, %d])' % int(ymax * N))]
        ps = [p for p in ps if p != N]
        T2 = []; T4 = []; Ls = {r: [] for r in range(0, 4)}
        for p in ps:
            s2 = Fr(0); s4 = Fr(0)
            h0 = H(4 * p * N)
            s2 += h0 / 2; s4 += h0 / 2
            Ls[0].append(math.pi * float(h0) / math.sqrt(4 * p * N))
            r = 1
            while r * r * N < 4 * p:
                D = 4 * p * N - r * r * N * N
                hv = H(D)
                s2 += hv; s4 -= (Fr(r * r * N, p) - 1) * hv
                if r < 4:
                    Ls[r].append(math.pi * float(hv) / math.sqrt(D))
                r += 1
            s2 -= p + 1
            T2.append(float(s2)); T4.append(float(s4))
        ps = np.array(ps, dtype=float); y = ps / N; wl = np.log(ps)
        for k, T in ((2, np.array(T2)), (4, np.array(T4))):
            dim = int(pari('mfdim([%d, %d], 0)' % (N, k)))
            E = T / dim
            res = {nm: E - np.array([fn(t, k) for t in y]) for nm, fn in (('M', M_closed), ('Zub', M_zub))}
            edges = np.arange(0, ymax + 1e-9, h)
            chi2 = {'M': 0.0, 'Zub': 0.0}; nb = 0; dev = []
            for a, b in zip(edges[:-1], edges[1:]):
                sel = (y > a) & (y <= b); w = wl[sel]
                for nm in chi2:
                    x = res[nm][sel]
                    m = np.sum(w * x) / np.sum(w)
                    se = math.sqrt(np.sum(w ** 2 * (x - m) ** 2)) / np.sum(w)
                    chi2[nm] += (m / se) ** 2
                    if nm == 'M':
                        dev.append(m)
                nb += 1
            rms = float(np.sqrt(np.mean(np.square(dev))))
            check('G%d one level N = %d, k = %d: window means (h = %.2f, y <= %g) of eps_f sqrt(p) lambda_f(p) vs M_k' % (k // 2, N, k, h, ymax),
                  chi2['M'] / nb < 2.0 and chi2['Zub'] > 3 * chi2['M'],
                  'chi2 %.1f on %d windows (Zubrilina squarefree-level density: %.1f), rms %.3f, %d primes'
                  % (chi2['M'], nb, chi2['Zub'], rms, len(ps)))
        zs = []
        for r in range(0, 4):
            a = np.array(Ls[r]); pred = math.pi ** 2 / 6 if r == 0 else nu_prime_prod(r)
            zs.append((a.mean() - pred) / (a.std() / math.sqrt(len(a))))
        check('G3 N = %d: mean over primes p of pi H(D_r)/sqrt(D_r) = zeta(2), nu\'(1), nu\'(2), nu\'(3) (r = 0..3)' % N,
              max(abs(z) for z in zs) < 4, 'z-scores ' + ', '.join('%.2f' % z for z in zs))


# ---------------------------------------------------------------- H
def sec_H():
    z = np.load(DATA)
    edges, T, dim, T2, cnt = z['edges'], z['T'], z['dim'], z['T2'], z['cnt']
    mean = T / dim
    err = np.sqrt(np.maximum(T2 / dim - mean ** 2, 0) / np.maximum(cnt, 1))
    model = []
    for a, b in zip(edges[:-1], edges[1:]):
        g = np.linspace(a + 1e-9, b, 401)
        model.append(float(np.trapezoid([M_closed(t) for t in g], g) / (b - a)))
    model = np.array(model)
    rms = float(np.sqrt(np.mean((mean - model) ** 2)))
    samp = z['prime'] / dim
    d2 = float(np.max(np.abs(samp - model)))
    check('H1 round-160 blind data (110 prime levels 2300-3200): rms(data - M_2 bin integral) <= 0.10', rms <= 0.10,
          'rms %.3f, chi2 %.1f on %d bins' % (rms, float(np.sum(((mean - model) / err) ** 2)), len(model)))
    check('H2 this note\'s M_2 (bin integral) = the registered model sampled at the primes, per bin', d2 < 0.02, 'max dev %.4f' % d2)


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'ABCDEFGH'
    for ch in which:
        globals()['sec_' + ch]()
    print('%d checks, %d failed' % (NCHK, NFAIL))
