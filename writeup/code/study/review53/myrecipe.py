"""My own implementation of the CFKRS/Heap recipe for k = 1, 2, 3 distinct real primitive characters.
Independent of tmlib.py: own Kronecker symbol, own field discriminant of D_i D_j, own local factors, own Euler product
A (complete homogeneous polynomials by convolution of geometric series, not the e1/e2/e3 recurrence), multilinear
coefficients by finite differences of log h (Richardson), and a direct evaluation of the 2^k swap terms at small shifts.
PARI is used only for L(s, chi*) at real s (cross-checked at s = 1 by the class number formula)."""
import math, functools
import numpy as np
import mpmath as mp
import cypari2

pari = cypari2.Pari(); pari.allocatemem(200 * 10**6); pari.set_real_precision(30)
EG = float(mp.euler)


def factor(n):
    n = abs(int(n)); f = {}; p = 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def kron(D, p):
    """Kronecker symbol (D/p) for a prime p."""
    if p == 2:
        if D % 2 == 0:
            return 0
        return 1 if D % 8 in (1, 7) else -1
    a = D % p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def field_disc(n):
    """discriminant of Q(sqrt(n)) for a non-square integer n."""
    f = factor(n); sf = (1 if n > 0 else -1)
    for p, e in f.items():
        if e % 2:
            sf *= p
    return sf if sf % 4 == 1 else 4 * sf


@functools.lru_cache(maxsize=None)
def primes(N):
    s = np.ones(N + 1, bool); s[:2] = False
    for p in range(2, int(N ** .5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0]


@functools.lru_cache(maxsize=None)
def chivec(D, N):
    return np.array([kron(D, int(p)) for p in primes(N)], dtype=float)


@functools.lru_cache(maxsize=None)
def Ldata(Ds):
    return pari.lfuncreate(Ds)


@functools.lru_cache(maxsize=None)
def logLam(Di, Dj, u):
    """log Lambda_ij(1+u) = log sum_n chi_i(n) chi_j(n) n^{-1-u}, u real."""
    Ds = field_disc(Di * Dj)
    s = 1 + u
    val = math.log(float(pari.lfun(Ldata(Ds), pari(repr(s)))))
    for p in factor(Di * Dj):
        c = kron(Ds, p)
        if c != 0:
            val += math.log(1 - c * p ** (-s))
    return val


def logc(D, x):
    return sum(math.log(1 - p ** (-1 - x)) for p in factor(D))


def hn_series(xs, N):
    """complete homogeneous h_0..h_N of the variables xs (list of arrays over primes) by convolution."""
    n = np.arange(N + 1)
    acc = None
    for x in xs:
        pw = x[:, None] ** n[None, :]
        if acc is None:
            acc = pw
        else:
            new = np.zeros_like(acc)
            for a in range(N + 1):
                new[:, a:] += acc[:, :N + 1 - a] * pw[:, a:a + 1]
            acc = new
    return acc


PMAX = 10**6


@functools.lru_cache(maxsize=None)
def tailsum(a):
    """sum_{p > PMAX} p^{-2-a} ~ int_PMAX^inf x^{-2-a} / ln x dx."""
    return float(mp.quad(lambda x: x ** (-2 - a) / mp.log(x), [PMAX, 10 * PMAX, mp.inf]))


def logA(Ds, z):
    """log of A(z) = prod_p R_p(x, x) prod_{i,j} (1 - x_i x_j), x_j = chi_j(p) p^{-1/2-z_j}; slice alpha = beta = z."""
    P = primes(PMAX).astype(float)
    tot = 0.0
    for lo, hi, N in ((0, 50, 160), (50, 2000, 60), (2000, 10**5, 18), (10**5, PMAX + 1, 9)):
        sel = (P >= lo) & (P < hi)
        pp = P[sel]
        xs = [chivec(D, PMAX)[sel] * pp ** (-0.5 - zz) for D, zz in zip(Ds, z)]
        H = hn_series(xs, N)
        Rp = (H * H).sum(axis=1)
        la = np.log(Rp)
        for i in range(len(Ds)):
            for j in range(len(Ds)):
                la += np.log(1 - xs[i] * xs[j])
        tot += la.sum()
    k = len(Ds)
    for i in range(k):
        for j in range(i + 1, k):
            tot -= tailsum(2 * (z[i] + z[j]))
    return tot


def logA_pairs(Ds, z):
    """sum over pairs of log A_ij(z_i, z_j) (the no-three-body replacement), brute force R_p for two variables."""
    k = len(Ds); tot = 0.0
    for i in range(k):
        for j in range(i + 1, k):
            tot += logA((Ds[i], Ds[j]), (z[i], z[j]))
    return tot


def logh(Ds, z, no3=False):
    k = len(Ds)
    v = sum(logc(Ds[j], 2 * z[j]) for j in range(k))
    for i in range(k):
        for j in range(i + 1, k):
            v += 2 * logLam(Ds[i], Ds[j], round(z[i] + z[j], 15))
    if k >= 2:
        v += logA_pairs(Ds, z) if (no3 and k == 3) else logA(Ds, z)
    return v


def multilinear(Ds, no3=False, eps=(0.004, 0.002, 0.001)):
    """h[m] for all subsets m, from mixed partials of log h by central differences + Richardson (O(eps^2) -> O(eps^6))."""
    k = len(Ds)
    import itertools
    subsets = [m for r in range(k + 1) for m in itertools.combinations(range(k), r)]
    est = {}
    for e in eps:
        cache = {}
        def f(pt):
            if pt not in cache:
                cache[pt] = logh(Ds, [x * e for x in pt], no3)
            return cache[pt]
        d = {}
        for m in subsets:
            if len(m) == 0:
                d[m] = f(tuple([0] * k)); continue
            acc = 0.0
            for sg in itertools.product((-1, 1), repeat=len(m)):
                pt = [0] * k
                for jj, s in zip(m, sg):
                    pt[jj] = s
                acc += np.prod(sg) * f(tuple(pt))
            d[m] = acc / (2 * e) ** len(m)
        est[e] = d
    # Richardson: errors in powers of eps^2 (eps halves): two levels
    e0, e1, e2 = eps
    l = {}
    for m in subsets:
        a, b, c = est[e0][m], est[e1][m], est[e2][m]
        r1 = (4 * b - a) / 3; r2 = (4 * c - b) / 3
        l[m] = (16 * r2 - r1) / 15
    # exp: multilinear coefficients of exp(l) (only multilinear monomials)
    h0 = math.exp(l[()])
    hm = {(): h0}
    # general formula: h[m] = h0 * sum over set partitions of m of prod l[block]
    def partitions(s):
        if not s:
            yield []
            return
        first = s[0]
        rest = s[1:]
        for p in partitions(rest):
            yield [[first]] + p
            for i in range(len(p)):
                yield p[:i] + [[first] + p[i]] + p[i + 1:]
    for m in subsets[1:]:
        hm[m] = h0 * sum(np.prod([l[tuple(sorted(b))] for b in part]) for part in partitions(list(m)))
    return hm, l


def P_at(Ds, hm, t):
    q = [abs(D) for D in Ds]; k = len(Ds)
    L = [math.log(qq * t / (2 * math.pi)) + 2 * EG for qq in q]
    out = 0.0
    for m, c in hm.items():
        term = c
        for j in range(k):
            if j not in m:
                term *= L[j]
        out += term
    return out


def P_poly(Ds, hm):
    """P(t) as polynomial in u = ln t: coefficients (numpy poly1d)."""
    q = [abs(D) for D in Ds]; k = len(Ds)
    out = np.poly1d([0.0])
    for m, c in hm.items():
        term = np.poly1d([c])
        for j in range(k):
            if j not in m:
                term = term * np.poly1d([1.0, math.log(q[j] / (2 * math.pi)) + 2 * EG])
        out = out + term
    return out


def int_poly_lnt(poly, T1=100.0, T2=1000.0):
    """exact int_T1^T2 poly(ln t) dt: int u^k e^u du."""
    def F(T):
        u = math.log(T); tot = 0.0
        for k, a in enumerate(poly.coeffs[::-1]):
            # int u^k e^u du = e^u sum_{j=0}^k (-1)^{k-j} k!/j! u^j
            s = sum((-1) ** (k - j) * math.factorial(k) / math.factorial(j) * u ** j for j in range(k + 1))
            tot += a * T * s
        return tot
    return F(T2) - F(T1)


def swap_sum(Ds, t, e, dirs=(1.0, 1.31, 1.73)):
    """direct sum of the 2^k swap terms at alpha_j = beta_j = z_j = e * dirs_j (mpmath zeta, own h)."""
    import itertools
    k = len(Ds); z = [e * d for d in dirs[:k]]
    X = [abs(D) * t / (2 * math.pi) for D in Ds]
    tot = mp.mpf(0)
    for s in itertools.product((1, -1), repeat=k):
        sz = [si * zi for si, zi in zip(s, z)]
        term = mp.e ** mp.mpf(logh(Ds, sz))
        for j in range(k):
            term *= mp.zeta(1 + 2 * sz[j])
            if s[j] == -1:
                term *= mp.mpf(X[j]) ** (-2 * z[j])
        tot += term
    return float(tot)


def diag_coeffs(Ds, eps=(0.004, 0.002, 0.001)):
    """g(sigma) = h(sigma/2,...) Taylor coefficients g_0..g_k by central differences (Richardson)."""
    k = len(Ds)
    res = []
    for e in eps:
        f = {i: math.exp(logh(Ds, [i * e / 2] * k)) for i in range(-3, 4)}
        g0 = f[0]
        g1 = (f[1] - f[-1]) / (2 * e)
        g2 = (f[1] - 2 * f[0] + f[-1]) / e ** 2 / 2
        g3 = (f[2] - 2 * f[1] + 2 * f[-1] - f[-2]) / (2 * e ** 3) / 6
        res.append([g0, g1, g2, g3][:k + 1])
    a, b, c = [np.array(r) for r in res]
    r1 = (4 * b - a) / 3; r2 = (4 * c - b) / 3
    return (16 * r2 - r1) / 15


def stieltjes_szeta(k):
    """coefficients of (sigma zeta(1+sigma))^k up to sigma^k."""
    st = [float(mp.stieltjes(n)) for n in range(k + 1)]
    zs = [1.0] + [(-1) ** n * st[n] / math.factorial(n) for n in range(k)]
    zk = [1.0] + [0.0] * k
    for _ in range(k):
        new = [0.0] * (k + 1)
        for a in range(k + 1):
            for b in range(k + 1 - a):
                new[a + b] += zk[a] * zs[b]
        zk = new
    return zk


def diag_density(Ds, g, t, W=None):
    """2 Res_{s=0} zeta(1+s)^k g(s) C^{s/2} Wt(s), Wt(s) = 1/s + w0 + w1 s + w2 s^2 (+...); W = [w0, w1, w2] (sharp: 0)."""
    k = len(Ds)
    zk = stieltjes_szeta(k)
    c = [sum(zk[a] * g[n - a] for a in range(n + 1)) for n in range(k + 1)]
    ell = sum(math.log(abs(D) * t / (2 * math.pi)) for D in Ds) / 2
    Wc = [1.0] + list(W if W is not None else [0.0] * k)
    # coefficient of s^k in c(s) * exp(ell s) * (1 + w0 s + w1 s^2 + ...)
    E = [ell ** n / math.factorial(n) for n in range(k + 1)]
    tot = 0.0
    for a in range(k + 1):
        for b in range(k + 1 - a):
            cidx = k - a - b
            if cidx < len(Wc):
                tot += c[a] * E[b] * Wc[cidx]
    return 2 * tot


def classno_check(D):
    """L(1, chi_D) by PARI lfun vs class number formula."""
    L1 = float(pari.lfun(Ldata(D), 1))
    if D < 0:
        h = int(pari.qfbclassno(D)); w = 6 if D == -3 else (4 if D == -4 else 2)
        ref = 2 * math.pi * h / (w * math.sqrt(-D))
    else:
        q = pari.quadclassunit(D); h = int(q[0]); reg = float(q[3])
        ref = 2 * h * reg / math.sqrt(D)
    return L1, ref
