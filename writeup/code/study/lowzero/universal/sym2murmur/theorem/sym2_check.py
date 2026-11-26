"""sym2_check.py -- checks every numeric claim of sym2_theorem.md (own code; nothing imported from the repository).

Sections (run all by default, or name some on the command line):
  A  trace identity: own Eichler-Selberg (Kaplan-Petrow Thm ESTFKL form) against PARI; the decomposition of Tr T_{p^2}
     into identity + hyperbolic + A2chi + A2N2; the reduction t = +-(2p - m N^2); decomposition = independent Gram route
  B  Hurwitz class numbers (own table) against PARI
  C  local means by exact enumeration: E_l(m) = nu_l(v_l(m)) (Lemma L), N-adic means, the twisted N-adic mean
     -2/(N^2-1), the untwisted means E'_l = 1 (symbolic + enumeration), nu(m) values, M_k values
  F  Fourier/Bessel identity, Weil bound, N-adic dual coefficients (used in Proposition H)
  D  class-number key lemma at one large level N = 10007 (PARI class numbers, random primes in a window)
  E  one-level test of the full statement: N = 211 (k = 2), N = 139 (k = 4), all primes p <= 3N^2, via the Gram route;
     exact split into identity + hyperbolic + class-number + twisted part R; per-prime fluctuation of R
Final line: "<n> checks, <f> failed".
"""
import sys, time, math, random
from math import isqrt, gcd, pi, sqrt, log
from fractions import Fraction
import numpy as np

T0 = time.time()
SECTIONS = set(sys.argv[1:]) or set("ABCFDE")
CHECKS = []


def check(name, ok, info=""):
    CHECKS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + ("  " + info if info else ""), flush=True)


_pari = None


def gp():
    global _pari
    if _pari is None:
        import cypari2
        _pari = cypari2.Pari()
        try:
            _pari.allocatemem(300 * 10 ** 6, silent=True)
        except TypeError:
            _pari.allocatemem(300 * 10 ** 6)
    return _pari


# ---------------------------------------------------------------- arithmetic helpers
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def factor_small(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def divisors(n):
    ds = [1]
    for q, e in factor_small(n).items():
        ds = [d * q ** j for d in ds for j in range(e + 1)]
    return sorted(ds)


def sigma(n, s=1):
    return sum(d ** s for d in divisors(n))


def legendre_table(N):
    leg = -np.ones(N, np.int64); leg[0] = 0
    leg[(np.arange(1, N, dtype=np.int64) ** 2) % N] = 1
    return leg


def hurwitz12_table(X):
    """12*H(D) for 0 <= D <= X (H(0) := 0): weighted count of reduced forms (a,b,c), 4ac - b^2 = D, |b| <= a <= c,
    b >= 0 if |b| = a or a = c; (a,0,a) has weight 1/2, (a,a,a) weight 1/3."""
    H = np.zeros(X + 1, np.int32)
    a = 1
    while 3 * a * a <= X:
        for b in range(a + 1):                       # c = a
            D = 4 * a * a - b * b
            if D <= X:
                H[D] += 6 if b == 0 else (4 if b == a else 12)
        st, base = 4 * a, 4 * a * (a + 1)
        for b in range(a + 1):                       # c > a
            s = base - b * b
            if s <= X:
                H[s::st] += 12 if (b == 0 or b == a) else 24
        a += 1
    return H


def Ucheb_int(t, n, k):
    """normalized Chebyshev U_{k-2}(t, n) = (al^{k-1} - alb^{k-1})/(al - alb), al + alb = t, al*alb = n (int64 arrays)."""
    if k == 2:
        return np.ones_like(t)
    if k == 4:
        return t * t - n
    raise ValueError


def Ucheb(x, k):
    """U_{k-2}(x), the Chebyshev polynomial of the second kind."""
    return np.ones_like(x) if k == 2 else 4 * x * x - 1


# ---------------------------------------------------------------- the trace formula (Kaplan-Petrow, Thm ESTFKL; chi = 1, N prime)
def trace_T(n, N, k, H12, leg):
    """Tr T_n on S_k(Gamma0(N)), N prime, gcd(n, N) = 1, k in {2, 4}: exact integer."""
    assert n % N != 0
    T = isqrt(4 * n)
    if T * T == 4 * n:
        T -= 1
    t = np.arange(-T, T + 1, dtype=np.int64)
    D = 4 * n - t * t
    term = H12[D].astype(np.int64) * (1 + leg[(-D) % N])
    msk = (D % (N * N) == 0)
    if msk.any():
        term[msk] += N * H12[D[msk] // (N * N)].astype(np.int64)
    S2 = int((Ucheb_int(t, n, k) * term).sum())          # 12 * elliptic sum
    s = isqrt(n)
    A1x24 = 2 * (k - 1) * (N + 1) * s ** (k - 2) if s * s == n else 0
    A3 = -sum(min(b, n // b) ** (k - 1) for b in divisors(n))    # tau in {1, N}, each phi(1) = 1: factor 2 * (1/2)
    A4 = sigma(n) if k == 2 else 0
    tot = A1x24 - S2 + 24 * (A3 + A4)
    assert tot % 24 == 0, (n, N, k, tot)
    return tot // 24


def trace_p2_decomposed(p, N, k, H12, leg):
    """Tr T_{p^2} | S_k(Gamma0(N)) (k with S_k(SL2Z) = 0) as identity + hyperbolic + A2chi + A2N2; returns the parts x24."""
    t = np.arange(-(2 * p - 1), 2 * p, dtype=np.int64)
    D = 4 * p * p - t * t
    U = Ucheb_int(t, p * p, k)
    chi = leg[(-D) % N]                                    # chi_N(t^2 - 4p^2), 0 when N | D
    A2chi_x24 = -int((U * H12[D].astype(np.int64) * chi).sum())
    ms = np.arange(1, (4 * p - 1) // (N * N) + 1, dtype=np.int64)
    ms = ms[ms * N * N < 4 * p]
    tm = 2 * p - ms * N * N
    Dm = ms * (4 * p - ms * N * N)
    A2N2_x24 = -2 * N * int((Ucheb_int(tm, p * p, k) * H12[Dm].astype(np.int64)).sum())
    I_x24 = 2 * (k - 1) * N * p ** (k - 2)
    Hyp_x24 = -12 * (2 + p ** (k - 1))
    return I_x24, Hyp_x24, A2chi_x24, A2N2_x24


# ---------------------------------------------------------------- Gram route (independent of the decomposition)
class Gram:
    """sum_f lambda_f(p)^2 over the Hecke eigenbasis of S_k(Gamma0(N)) from traces of T_n with n = m p, m in S (rank = dim)."""

    def __init__(self, N, k, H12, leg):
        self.N, self.k, self.H12, self.leg = N, k, H12, leg
        self.cache = {}
        self.dim = self.tr(1)
        S, m = [], 0
        while len(S) < self.dim:
            m += 1
            if m % N == 0:
                continue
            S2 = S + [m]
            G = self.Gn(S2)
            if np.linalg.matrix_rank(G, tol=1e-9 * np.abs(G).max()) == len(S2):
                S = S2
        self.S = S
        self.Gint = [[self.pair(a, b) for b in S] for a in S]
        self.G = self.Gn(S)
        self.Ginv = np.linalg.inv(self.G)

    def tr(self, n):
        if n not in self.cache:
            self.cache[n] = trace_T(n, self.N, self.k, self.H12, self.leg)
        return self.cache[n]

    def pair(self, a, b):
        g = gcd(a, b)
        return sum(d ** (self.k - 1) * self.tr(a * b // (d * d)) for d in divisors(g))

    def Gn(self, S):
        nr = np.array(S, float) ** ((self.k - 1) / 2)
        return np.array([[self.pair(a, b) for b in S] for a in S], float) / np.outer(nr, nr)

    def wvec(self, p):
        assert p > max(self.S)
        return [trace_T(m * p, self.N, self.k, self.H12, self.leg) for m in self.S]

    def s_of_p(self, p, w=None):
        if w is None:
            w = self.wvec(p)
        nr = (np.array(self.S, float) * p) ** ((self.k - 1) / 2)
        wn = np.array(w, float) / nr
        return float(wn @ self.Ginv @ wn)

    def s_exact(self, p):
        """exact rational sum_f a_f(p)^2 = w^T G^{-1} w (unnormalized, integer G)."""
        import sympy
        w = sympy.Matrix(self.wvec(p))
        G = sympy.Matrix(self.Gint)
        return (w.T * G.LUsolve(w))[0, 0]


# ---------------------------------------------------------------- local factors F_l(D) and nu
def F_local(D, l):
    """exact local factor F_l(D) = Q_l/(1 - chi(l)/l) as a Fraction (T2 of the template), D > 0, -D = 0, 1 mod 4."""
    v, u = 0, D
    while u % l == 0:
        u //= l; v += 1
    x = Fraction(1, l)
    if l != 2:
        a = v // 2
        chi = 0 if v % 2 else (1 if pow((-u) % l, (l - 1) // 2, l) == 1 else -1)
    else:
        w = (-u) % 8                                      # -D = 2^v * w
        if v % 2 == 0:
            if w % 4 == 1:
                a, chi = v // 2, (1 if w == 1 else -1)
            else:
                a, chi = v // 2 - 1, 0
        else:
            a, chi = (v - 3) // 2, 0
    S = sum(x ** j for j in range(a))
    return S + x ** a / (1 - chi * x)


def F_local_class(Dmod, l, J):
    """F_l for a residue Dmod mod l^J (Dmod in [0, l^J)): exact value if determined, else an enclosure (lo, hi)."""
    q = l ** J
    if Dmod == 0:
        v = J
    else:
        v, u = 0, Dmod
        while u % l == 0:
            u //= l; v += 1
    need = J - (3 if l == 2 else 1)                       # v and the unit part mod l (mod 8) are known if v <= need
    if v <= need:
        rep = Dmod
        # a positive representative with the same l-adic data; -D = 0,1 mod 4 is automatic for the residues used
        while rep <= 0 or (l != 2 and (-rep) % 4 not in (0, 1)):
            rep += q
        if l != 2:
            return F_local(rep, l), F_local(rep, l)
        return F_local(rep, l), F_local(rep, l)
    x = Fraction(1, l)
    amin = max(0, (need + 1) // 2 - 2)
    S = sum(x ** j for j in range(amin))
    return S + x ** amin / (1 + x), 1 / (1 - x)


def F_vec(Dm, l, J):
    """vectorized F_l on residues D mod l^J (float): exact where v_l(D) and the unit class are determined, else an enclosure."""
    q = l ** J
    u = np.asarray(Dm, dtype=np.int64) % q
    und = (u == 0)
    u = np.where(und, 1, u)
    v = np.zeros(len(u), np.int64)
    for _ in range(J):
        msk = (u % l == 0)
        if not msk.any():
            break
        u = np.where(msk, u // l, u); v += msk
    need = J - (3 if l == 2 else 1)
    und |= v > need
    x = 1.0 / l
    if l != 2:
        sq = np.zeros(l, bool); sq[(np.arange(1, l) ** 2) % l] = True
        a = v // 2
        chi = np.where(v % 2 == 1, 0, np.where(sq[(-u) % l], 1, -1))
    else:
        w = (-u) % 8
        even = (v % 2 == 0)
        a = np.where(even, np.where(w % 4 == 1, v // 2, v // 2 - 1), (v - 3) // 2)
        chi = np.where(even & (w % 4 == 1), np.where(w == 1, 1, -1), 0)
    assert (a[~und] >= 0).all()
    a = np.maximum(a, 0)
    F = (1 - x ** a) / (1 - x) + x ** a / (1 - chi * x)
    amin = max(0, (need + 1) // 2 - 2)
    lo = np.where(und, (1 - x ** amin) / (1 - x) + x ** amin / (1 + x), F)
    hi = np.where(und, 1 / (1 - x), F)
    return lo, hi


def nu_l(e, l):
    """closed form of Lemma L: nu_l(0) = 1 - u_l, nu_l(e) = S_a + x^a/(1-x^2) (e = 2a), S_a + x^a (e = 2a+1)."""
    x = Fraction(1, l)
    if e == 0:
        return 1 - Fraction(1, (l - 1) ** 2 * (l + 1))
    a = e // 2
    S = sum(x ** j for j in range(a))
    return S + (x ** a / (1 - x * x) if e % 2 == 0 else x ** a)


def nu_l_derivation(e, l):
    """the closed form printed in derivation.md, sect. 0."""
    if e == 0:
        return 1 - Fraction(1, (l - 1) * (l * l - 1))
    s = 1 + Fraction(1, l * l - 1) + Fraction(l, l + 1) * sum(Fraction(1, l ** j) for j in range(1, e // 2 + 1))
    if e % 2:
        s -= Fraction(1, l ** ((e - 1) // 2)) / (l * l - 1)
    return s


_C0 = None


def C0():
    global _C0
    if _C0 is None:
        ls = primes_upto(2 * 10 ** 6).astype(float)
        _C0 = float(np.exp(np.sum(np.log1p(-1.0 / ((ls - 1) ** 2 * (ls + 1))))))
    return _C0


def nu(m):
    r = C0()
    for l, e in factor_small(m).items():
        r *= float(nu_l(e, l) / nu_l(0, l))
    return r


def M_k(y, k):
    y = float(y)
    s = 0.0
    m = 1
    while m < 4 * y:
        s += nu(m) * sqrt(m * (4 * y - m)) * float(Ucheb(np.array(1 - m / (2 * y)), k))
        m += 1
    return 24 * (k == 2) * y ** 1.5 - 6 * sqrt(y) / (k - 1) - (12 / (pi * (k - 1) * sqrt(y))) * s if y > 0 else 0.0


# =====================================================================================================================
H12 = None


def need_H(X):
    global H12
    if H12 is None or len(H12) <= X:
        t = time.time()
        H12 = hurwitz12_table(X)
        print("  [H table to %d in %.1fs]" % (X, time.time() - t), flush=True)
    return H12


if "A" in SECTIONS or "B" in SECTIONS:
    need_H(9200000)

# ---------------------------------------------------------------- B: class numbers
if "B" in SECTIONS:
    P = gp()
    ok = all(int(H12[D]) == int(12 * P.qfbhclassno(D)) for D in range(1, 5001))
    check("B1 own 12H(D) = 12*PARI qfbhclassno(D) for all 1 <= D <= 5000", ok)
    rng = random.Random(1)
    Ds = [rng.randrange(1, len(H12)) for _ in range(1500)]
    ok = all(int(H12[D]) == int(12 * P.qfbhclassno(D)) for D in Ds)
    check("B2 own table = PARI at 1500 random D <= %d" % (len(H12) - 1), ok)

# ---------------------------------------------------------------- A: trace identity
if "A" in SECTIONS:
    P = gp()
    for N in (37, 43):
        for k in (2, 4):
            leg = legendre_table(N)
            co = [int(c) for c in P("mfcoefs(mftraceform([%d,%d],1),600)" % (N, k))]
            bad = [n for n in range(1, 601) if n % N and trace_T(n, N, k, H12, leg) != co[n]]
            check("A1 own Eichler-Selberg Tr T_n = PARI trace form, N=%d k=%d, all n <= 600 prime to N" % (N, k),
                  not bad, "mismatches %s" % bad[:5])
    for N in (11, 23):
        for k in (2, 4):
            leg = legendre_table(N)
            ps = [int(p) for p in primes_upto(150) if p != N]
            co = [int(c) for c in P("mfcoefs(mftraceform([%d,%d],1),%d)" % (N, k, 150 ** 2))]
            bad, nm = [], 0
            for p in ps:
                parts = trace_p2_decomposed(p, N, k, H12, leg)
                nm += parts[3] != 0
                if sum(parts) % 24 or sum(parts) // 24 != co[p * p] or trace_T(p * p, N, k, H12, leg) != co[p * p]:
                    bad.append(p)
            check("A2 decomposition I + Hyp + A2chi + A2N2 = PARI Tr T_{p^2}, N=%d k=%d, all primes p <= 150" % (N, k),
                  not bad, "(%d primes with m-terms) mismatches %s" % (nm, bad[:5]))
    ok, cnt = True, 0
    for N in (11, 23, 37, 211):
        for p in primes_upto(3 * N * N):
            p = int(p)
            if p == N:
                continue
            ts = [t for t in range(-(2 * p - 1), 2 * p) if (4 * p * p - t * t) % (N * N) == 0] if p < 400 else None
            ms = [m for m in range(1, 4 * p) if m * N * N < 4 * p]
            pred = sorted(set([2 * p - m * N * N for m in ms] + [-(2 * p - m * N * N) for m in ms]))
            if ts is not None:
                ok &= (sorted(ts) == pred); cnt += 1
            ok &= all((4 * p * p - t * t) == N * N * m * (4 * p - m * N * N) for m, t in zip(ms, [2 * p - m * N * N for m in ms]))
            ok &= len(pred) == 2 * len(ms)
    check("A3 {|t|<2p : N^2 | 4p^2-t^2} = {+-(2p - mN^2), 1 <= m < 4p/N^2}, quotient m(4p - mN^2), distinct",
          ok, "(%d primes enumerated in full)" % cnt)
    for N, k in ((37, 2), (37, 4), (43, 2), (43, 4), (211, 2), (139, 4)):
        leg = legendre_table(N)
        g = Gram(N, k, H12, leg)
        bad, worst = [], 0.0
        for p in primes_upto(1500):
            p = int(p)
            if p == N or p <= max(g.S):
                continue
            tr2 = sum(trace_p2_decomposed(p, N, k, H12, leg)) // 24
            sa = g.s_of_p(p) * p ** (k - 1)                # sum_f a_f(p)^2
            worst = max(worst, abs(sa - round(sa)) / max(1.0, abs(sa)))
            if round(sa) - p ** (k - 1) * g.dim != tr2:
                bad.append(p)
        check("A4 decomposition = Gram route (sum a_f(p)^2 - p^{k-1} dim), N=%d k=%d dim=%d |S|=%d, primes %d<p<=1500"
              % (N, k, g.dim, len(g.S), max(g.S), ), not bad and worst < 1e-9, "max rel. dev. from integer %.1e" % worst)
        if N in (37, 211):
            ps = [p for p in [1009, 1013, 1019] if p != N]
            ok = all(g.s_exact(p) == round(g.s_of_p(p) * p ** (k - 1)) for p in ps)
            check("A5 Gram route exact over Q (sympy) = float route, N=%d k=%d, p in %s" % (N, k, ps), ok)
    for N in (37, 43, 139, 211):
        for k in (2, 4):
            d = trace_T(1, N, k, H12, legendre_table(N))
            check("A6 dim = Tr T_1 = PARI mfdim, N=%d k=%d" % (N, k), d == int(gp()("mfdim([%d,%d],1)" % (N, k))), "dim %d" % d)


# ---------------------------------------------------------------- C: local means and nu
if "C" in SECTIONS:
    import sympy
    # C1 Lemma L: E_l(m) = mean over p in Z_l^x of F_l(m(4p - mN^2)) equals nu_l(v_l(m))
    worst, ok, cases = 0.0, True, 0
    for l, J in ((2, 18), (3, 11), (5, 8), (7, 7)):
        q = l ** J
        ps = np.array([p for p in range(1, q) if p % l], dtype=np.int64)
        for e in range(0, 6):
            for mp in (1, 5 if l != 5 else 3, 7 if l != 7 else 11):
                for N in (10007, 13, 101):
                    m = l ** e * mp
                    lo, hi = F_vec((m % q) * ((4 * ps - (m * N * N) % q) % q), l, J)
                    lo, hi = lo.mean(), hi.mean()
                    tgt = float(nu_l(e, l))
                    cases += 1
                    ok &= (lo - 1e-12 <= tgt <= hi + 1e-12)
                    worst = max(worst, hi - lo)
    check("C1 Lemma L by exact enumeration: nu_l(e) in [lo, hi] for l = 2,3,5,7, e = 0..5, 3 unit parts, 3 levels",
          ok and worst < 2e-3, "(%d cases, max enclosure width %.1e)" % (cases, worst))
    # C2 l = N: mean of F_N(D_m) over p mod N (N does not divide D_m) is N^2/(N^2-1)
    ok = True
    for N in (13, 101):
        for m in range(1, 12):
            vals = [F_local(((m * (4 * p - m * N * N)) % N) or N, N) if (m * (4 * p - m * N * N)) % N else None for p in range(1, N)]
            ok &= all(v is not None for v in vals) and sum(vals) / (N - 1) == Fraction(N * N, N * N - 1)
    check("C2 N-adic mean of F_N(m(4p - mN^2)) over p in (Z/N)^x = N^2/(N^2-1) (N = 13, 101; m < 12)", ok)
    # C3 closed form = derivation's formula
    ok = all(nu_l(e, l) == nu_l_derivation(e, l) for l in primes_upto(60) for e in range(0, 13))
    check("C3 closed form of Lemma L = the formula of derivation.md sect. 0 (exact, l < 60, e <= 12)", ok)
    # C4 twisted N-adic mean: E_t[chi_N(t^2 - 4p^2) F_N(4p^2 - t^2)] = -2/(N^2-1) for every p in (Z/N)^x
    ok = True
    for N in (5, 7, 11, 13, 101):
        for p in range(1, N):
            s = Fraction(0)
            for t in range(N):
                c = (t * t - 4 * p * p) % N
                if c:
                    chi = 1 if pow(c, (N - 1) // 2, N) == 1 else -1
                    s += chi / (1 - Fraction(chi, N))
            ok &= (s / N == Fraction(-2, N * N - 1))
    check("C4 twisted N-adic mean = -2/(N^2-1) exactly, N = 5,7,11,13,101, every p mod N", ok)
    # C5 untwisted means E'_l = E_t[F_l(4p^2 - t^2)] = 1 for l not dividing 2pN: symbolic (odd l) + enumeration
    x, L_ = sympy.symbols('x L', positive=True)
    S = lambda a: (1 - x ** a) / (1 - x)
    A_ = (1 / (1 - x)) * (1 / (1 - x ** 2) - 1 / (1 - x ** 3))          # sum_{a>=0} x^{2a} S_a
    B_ = 1 / (1 - x ** 3)                                                 # sum_{a>=0} x^{3a}
    Cc = (1 / (1 - x)) * (x / (1 - x ** 2) - x ** 2 / (1 - x ** 3))       # sum_{a>=1} x^{2a-1} S_a
    Dd = x ** 2 / ((1 - x ** 2) * (1 - x ** 3))                           # sum_{a>=1} x^{3a-1}/(1-x^2)
    l_ = 1 / x
    Eprime = x * ((l_ - 3) / 2 / (1 - x) + (l_ - 1) / 2 / (1 + x)) + 2 * x * (1 - x) * (A_ + B_ + Cc + Dd)
    sym_ok = sympy.simplify(Eprime - 1) == 0
    # numeric check of the four geometric sums
    xs = sympy.Rational(1, 7)
    num = sum(xs ** (2 * a) * S(a).subs(x, xs) for a in range(60)) - A_.subs(x, xs)
    sym_ok &= abs(float(num)) < 1e-30
    worst, ok = 0.0, True
    for l, J in ((2, 18), (3, 11), (5, 8), (7, 7)):
        q = l ** J
        ts = np.arange(q, dtype=np.int64)
        for p in ((1, 3, 5, 7) if l == 2 else (1, 2, 3)):
            if p % l == 0:
                continue
            lo, hi = F_vec((4 * p * p - (ts * ts) % q) % q, l, J)
            lo, hi = lo.mean(), hi.mean()
            ok &= lo - 1e-12 <= 1 <= hi + 1e-12
            worst = max(worst, hi - lo)
    check("C5 untwisted local means E'_l = 1: rational identity in x = 1/l (sympy) and enumeration l = 2,3,5,7",
          sym_ok and ok and worst < 2e-3, "(max enclosure width %.1e)" % worst)
    # C6 nu values printed in derivation.md
    ref = [0.6151, 0.9227, 0.6561, 1.5378, 0.6216, 0.9842, 0.6173, 1.3840, 0.9022, 0.9324, 0.6156, 1.6404, 0.6154, 0.9259,
           0.6630, 1.6916]
    dev = max(abs(nu(m) - r) for m, r in zip(range(1, 17), ref))
    check("C6 nu(1..16) and C0 = prod(1-u_l) agree with derivation.md", dev < 6e-5 and abs(C0() - 0.615133) < 1e-6,
          "C0 = %.6f, max dev %.1e" % (C0(), dev))
    ref2 = {2: [-1.138, 0.920, 2.540, 2.768, 3.523], 4: [-0.632, -0.307, 0.350, 2.382, 2.338]}
    dev = max(abs(M_k(y, k) - r) for k in (2, 4) for y, r in zip((0.1, 0.5, 1, 2, 3), ref2[k]))
    check("C7 M_k(y) at y = 0.1, 0.5, 1, 2, 3 (k = 2, 4) agree with derivation.md", dev < 1.5e-3, "max dev %.1e" % dev)

# ---------------------------------------------------------------- F: analytic inputs of Proposition H
if "F" in SECTIONS:
    from scipy.integrate import quad
    from scipy.special import jv
    dev = 0.0
    for k in (2, 4):
        for lam in (0.5, 3.0, 11.0, 40.0):
            val = quad(lambda th: math.sin((k - 1) * th) * math.sin(th) * math.cos(lam * math.cos(th)), 0, pi,
                       limit=400, epsabs=1e-13)[0]
            ref = pi * (k - 1) * (-1) ** ((k - 2) // 2) * jv(k - 1, lam) / lam
            dev = max(dev, abs(val - ref))
    check("F1 int_{-1}^{1} U_{k-2}(u) sqrt(1-u^2) e^{-i lam u} du = pi(k-1)(-i)^{k-2} J_{k-1}(lam)/lam (k = 2, 4)",
          dev < 1e-9, "max dev %.1e" % dev)
    N = 101
    leg = legendre_table(N)
    xs = np.arange(N)
    ratio, zero_ok = 0.0, True
    for c in range(1, N):
        f = leg[(xs * xs - c) % N].astype(float)
        F = np.fft.fft(f)
        zero_ok &= abs(F[0] + 1) < 1e-9
        ratio = max(ratio, np.abs(F[1:]).max() / sqrt(N))
    check("F2 Weil: |sum_x chi_N(x^2-c) e(ix/N)| <= 2 sqrt N for i != 0, = -1 for i = 0 (N = 101, all c != 0)",
          zero_ok and ratio <= 2 + 1e-9, "max ratio %.3f" % ratio)
    phimax, avgmax = 0.0, 0.0
    hats = {}
    for b in range(1, N):
        c = (xs * xs - 4 * b * b) % N
        chi = leg[c].astype(float)
        phi = np.where(c == 0, 0.0, chi / (1 - chi / N))
        hats[b] = np.fft.fft(phi) / N                       # hat(i) = (1/N) sum_x phi(x) e(-ix/N)
        phimax = max(phimax, np.abs(hats[b]).max() * sqrt(N))
    Hm = np.array([hats[b] for b in range(1, N)])               # rows b = 1..N-1, columns i
    bb = np.arange(1, N)[:, None]; ii = np.arange(N)[None, :]
    for sg in (1, -1):
        av = (Hm * np.exp(2j * pi * sg * 2 * bb * ii / N)).mean(axis=0)
        avgmax = max(avgmax, np.abs(av[1:]).max() * N ** 3)
    check("F3 |hat phi_N(i)| <= 3/sqrt(N) for the twisted N-adic factor (N = 101, every p mod N, every i)",
          phimax <= 3, "max sqrt(N)|hat| = %.3f" % phimax)
    check("F4 averaged dual coefficient |avg_b hat phi_N(i;b) e(+-2bi/N)| <= 3/N^3 for i != 0 (N = 101)",
          avgmax <= 3, "max N^3|avg| = %.3f" % avgmax)

# ---------------------------------------------------------------- D: class-number key lemma at N = 10007
if "D" in SECTIONS:
    import gmpy2
    P = gp()
    N = 10007
    rng = random.Random(7)
    ps = []
    while len(ps) < 2500:
        c = rng.randrange(2 * N * N + 1, 3 * N * N)
        if gmpy2.is_prime(c):
            ps.append(c)
    zs = []
    uN = 1 / ((N - 1) ** 2 * (N + 1))
    for m in range(1, 9):
        X = np.array([pi * float(P.qfbhclassno(m * (4 * p - m * N * N))) / sqrt(m * (4 * p - m * N * N)) for p in ps])
        tgt = nu(m) * N * N / ((N * N - 1) * (1 - uN))
        zs.append((X.mean() - tgt) / (X.std(ddof=1) / sqrt(len(X))))
        print("    m=%d mean %.4f +- %.4f  E(m) %.4f" % (m, X.mean(), X.std(ddof=1) / sqrt(len(X)), tgt))
    check("D1 one level N = 10007: mean of pi H(D_m)/sqrt(D_m) over 2500 random primes in (2N^2, 3N^2] = E(m), m = 1..8",
          max(abs(z) for z in zs) < 3.5, "z = " + " ".join("%.2f" % z for z in zs))

# ---------------------------------------------------------------- E: one-level test of the full statement
def trace_mp(m, p, N, k, H12, leg):
    """Tr T_{mp}, p prime > m, p != N (then A1 = 0, A3 = -2 sigma_{k-1}(m), A4 = sigma(m)(p+1))."""
    n = m * p
    T = isqrt(4 * n)
    if T * T == 4 * n:
        T -= 1
    t = np.arange(-T, T + 1, dtype=np.int64)
    D = 4 * n - t * t
    term = H12[D].astype(np.int64) * (1 + leg[(-D) % N])
    msk = (D % (N * N) == 0)
    if msk.any():
        term[msk] += N * H12[D[msk] // (N * N)].astype(np.int64)
    S2 = int((Ucheb_int(t, n, k) * term).sum())
    tot = -S2 - 48 * sigma(m, k - 1) + (24 * sigma(m) * (p + 1) if k == 2 else 0)
    assert tot % 24 == 0
    return tot // 24


def wmean(x, w):
    W = w.sum()
    mu = (w * x).sum() / W
    return mu, sqrt((w * w * (x - mu) ** 2).sum()) / W


if "E" in SECTIONS:
    need_H(9200000)
    for N, k in ((211, 2), (139, 4)):
        t1 = time.time()
        leg = legendre_table(N)
        g = Gram(N, k, H12, leg)
        d = g.dim
        sig = {m: (sigma(m, k - 1), sigma(m)) for m in g.S}
        ps = np.array([int(p) for p in primes_upto(3 * N * N) if p != N and p > max(g.S)])
        s = np.empty(len(ps))
        for i, p in enumerate(ps):
            s[i] = g.s_of_p(int(p), [trace_mp(m, int(p), N, k, H12, leg) for m in g.S])
        ok_range = bool((s > -1e-6).all() and (s < 4 * d + 1e-6).all())
        pf = ps.astype(float)
        SN = np.sqrt(pf) * (s - d) / d                      # (1/dim) sum_f sqrt(p) lambda_f(p^2)
        y = pf / N ** 2
        rho = (k - 1) * N / (12 * d)
        uN = 1 / ((N - 1) ** 2 * (N + 1))
        Em = {m: nu(m) * N * N / ((N * N - 1) * (1 - uN)) for m in range(1, 13)}
        planch = rho / np.sqrt(pf)
        hyp = -(np.sqrt(pf) / (2 * d)) * (1 + 2 * pf ** (1 - k))
        CNmod = np.zeros(len(ps)); CNex = np.zeros(len(ps)); CNnu1 = np.zeros(len(ps))
        for m in range(1, 13):
            sel = m * N * N < 4 * ps
            if not sel.any():
                continue
            Dm = m * (4 * ps[sel] - m * N * N)
            U = Ucheb(1 - m / (2 * y[sel]), k)
            CNex[sel] += -(N / (np.sqrt(pf[sel]) * d)) * U * H12[Dm] / 12.0
            base = -(N * N / (pi * np.sqrt(pf[sel]) * d)) * U * np.sqrt(m * (4 * y[sel] - m))
            CNmod[sel] += base * Em[m]
            CNnu1[sel] += base
        Rmain = (2 * pf ** 1.5 / ((N * N - 1) * d)) if k == 2 else np.zeros(len(ps))
        pred = planch + hyp + CNmod + Rmain
        Rex = SN - planch - hyp - CNex
        samp = range(0, len(ps), 40)
        dev0 = max(abs(pred[i] - rho * (1 / sqrt(pf[i]) + M_k(y[i], k)) + sqrt(pf[i]) * pf[i] ** (1 - k) / d) for i in samp)
        check("E0 N=%d k=%d: per-prime model = rho_N (1/sqrt p + M_k(p/N^2)) - sqrt(p) p^{1-k}/dim, up to factors N^2/(N^2-1) (%d primes)"
              % (N, k, len(samp)), dev0 < 5e-3, "max dev %.1e" % dev0)
        w = np.log(pf)
        edges = np.linspace(0, 3, 13)
        rows = []
        for a, b in zip(edges[:-1], edges[1:]):
            sel = (y > a) & (y <= b)
            o, se = wmean(SN[sel], w[sel])
            rows.append((a, b, o, se, wmean(pred[sel], w[sel])[0], wmean(hyp[sel], w[sel])[0],
                         wmean(CNmod[sel], w[sel])[0], wmean(CNnu1[sel], w[sel])[0], wmean(Rmain[sel], w[sel])[0],
                         wmean(Rex[sel], w[sel])[0], wmean(Rex[sel], w[sel])[1], wmean(CNex[sel], w[sel])[0], sel.sum()))
        R = np.array([r[2:12] for r in rows])
        o, se, pr, hy, cm, c1, rm, rx, rxse, cx = R.T
        chi = lambda model: float((((o - model) / se) ** 2).sum())
        X_model, X_zero, X_nohyp = chi(pr), chi(0 * pr), chi(pr - hy)
        X_nu1 = chi(pr - cm + c1)
        print("  N=%d k=%d dim=%d primes=%d (%.0fs)" % (N, k, d, len(ps), time.time() - t1))
        print("  window      obs     +-     M_k-model  R_obs   R_main")
        for r, rr in zip(rows, R):
            print("  (%.2f,%.2f] %7.3f %6.3f %7.3f   %7.3f %7.3f" % (r[0], r[1], rr[0], rr[1], rr[2], rr[7], rr[6]))
        print("  chi2 on 12 windows: model %.1f | zero %.1f | no hyperbolic %.1f | nu=1 %.1f" % (X_model, X_zero, X_nohyp, X_nu1)
              + (" | no R main term %.1f" % chi(pr - rm) if k == 2 else ""))
        check("E1 N=%d k=%d: 0 <= sum_f lambda_f(p)^2 <= 4 dim for all %d primes (Deligne; Gram route sane)" % (N, k, len(ps)),
              ok_range)
        check("E2 N=%d k=%d: window means of the full statistic vs rho_N(1/sqrt p + M_k): chi2 %.1f on 12 windows (< 26.2, p > 0.01)"
              % (N, k, X_model), X_model < 26.2,
              "alternatives: zero %.1f, no hyperbolic %.1f, nu = 1 %.1f" % (X_zero, X_nohyp, X_nu1))
        Xr = float((((rx - rm) / rxse) ** 2).sum())
        if k == 2:
            c = float(((rx * rm) / rxse ** 2).sum() / ((rm * rm) / rxse ** 2).sum()); sc = float(1 / sqrt(((rm * rm) / rxse ** 2).sum()))
            check("E3 N=%d k=2: twisted part R (exact split) vs 2p^{3/2}/((N^2-1)dim): chi2 %.1f on 12, coefficient %.3f +- %.3f"
                  % (N, Xr, c, sc), Xr < 26.2 and abs(c - 1) < 3 * sc)
        else:
            check("E3 N=%d k=4: twisted part R (exact split) vs 0: chi2 %.1f on 12 windows" % (N, Xr), Xr < 26.2)
        sel = CNmod != 0
        ratio = float((w[sel] * CNex[sel]).sum() / (w[sel] * CNmod[sel]).sum())
        jk = []                                             # jackknife over 20 blocks of primes
        idx = np.nonzero(sel)[0]
        for blk in np.array_split(idx, 20):
            keep = np.setdiff1d(idx, blk)
            jk.append((w[keep] * CNex[keep]).sum() / (w[keep] * CNmod[keep]).sum())
        jk = np.array(jk); sej = sqrt(19 / 20 * ((jk - jk.mean()) ** 2).sum())
        check("E4 N=%d k=%d: exact class-number term / its nu-model (pooled, y > 1/4) = %.4f +- %.4f" % (N, k, ratio, sej),
              abs(ratio - 1) < 4 * sej + 1e-3)
        if N == 211:
            # per-prime fluctuation: Sato-Tate size sqrt(p/dim) versus the random-sign model of the t-sum
            sel = (ps > 700) & (ps <= 1500)
            z = (SN[sel] - pred[sel]) * np.sqrt(d / pf[sel])
            vrs = []
            for p in ps[sel]:
                p = int(p)
                t = np.arange(-(2 * p - 1), 2 * p, dtype=np.int64)
                D = 4 * p * p - t * t
                Uv = np.ones(len(t)) if k == 2 else (t * t / (p * p) - 1.0)
                cnt = (leg[(-D) % N] != 0)
                vrs.append(((Uv * H12[D] / 12.0) ** 2 * cnt).sum() / (4 * p * d * d))
            vrs = np.array(vrs)
            r_obs = float(np.var(SN[sel] - pred[sel]) / vrs.mean())
            check("E5 N=211 k=2, 700 < p <= 1500 (%d primes): Var[(S_N - model) sqrt(dim/p)] = %.3f (Sato-Tate: 1); "
                  "observed variance / random-sign variance of the t-sum = %.4f" % (sel.sum(), float(np.var(z)), r_obs),
                  0.75 < np.var(z) < 1.3 and r_obs < 0.1)


print("[time %.0fs]" % (time.time() - T0), flush=True)
nf = sum(1 for _, o in CHECKS if not o)
print("%d checks, %d failed" % (len(CHECKS), nf))

