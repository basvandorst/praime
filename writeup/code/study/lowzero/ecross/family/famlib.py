"""famlib: Kloosterman sums, Petersson series, phases, Sato-Tate local algebra for the family one-point note.

Conventions: f a holomorphic newform of weight k, level N (prime or 1), trivial character, Hecke eigenvalues lambda_f(n)
(Deligne-normalised), L(s,f) = sum lambda_f(n) n^-s, Lambda(s,f) = N^{s/2} Gamma_C(s+(k-1)/2) L(s,f) = eps_f Lambda(1-s,f).
theta_f(t) = (t/2) ln N - t ln 2pi + Im lnGamma(k/2 + i t)  (common to the family), X(1/2+it)^{-1} = e^{2 i theta(t)}.
Harmonic weights omega_f = Gamma(k-1)/((4pi)^{k-1} <f,f>_{Gamma_0(N)\\H}) (Petersson: sum_f omega_f lam(m) lam(n) =
delta(m,n) + 2 pi i^{-k} sum_{c = 0 mod N} S(m,n;c)/c J_{k-1}(4 pi sqrt(mn)/c)).
"""
import math
import numpy as np
from math import gcd
from scipy.special import jv, loggamma

TWO_PI = 2 * math.pi


def theta(k, N, t):
    t = np.asarray(t, float)
    return t / 2 * math.log(N) - t * math.log(TWO_PI) + np.imag(loggamma(k / 2 + 1j * t))


def _powmod_vec(x, e, c):
    """x**e mod c elementwise, int64 (c < 3e9 fine since products < 9.2e18 need c < 3.03e9)."""
    r = np.ones_like(x)
    b = x % c
    while e:
        if e & 1:
            r = (r * b) % c
        b = (b * b) % c
        e >>= 1
    return r


def _phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


_unit_cache = {}


def units_and_inverses(c):
    if c in _unit_cache:
        return _unit_cache[c]
    if c == 1:
        res = (np.array([0], np.int64), np.array([0], np.int64))
    else:
        x = np.arange(1, c, dtype=np.int64)
        g = np.gcd(x, c)
        x = x[g == 1]
        inv = _powmod_vec(x, _phi(c) - 1, c)
        assert np.all((x * inv) % c == 1)
        res = (x, inv)
    if c < 20000:
        _unit_cache[c] = res
    return res


def kloosterman(m, n, c):
    """S(m,n;c) = sum_{x mod c, (x,c)=1} e((m x + n xbar)/c) (real)."""
    if c == 1:
        return 1.0
    x, inv = units_and_inverses(c)
    ph = ((m % c) * x + (n % c) * inv) % c
    return float(np.sum(np.cos(TWO_PI * ph / c)))


def kloosterman_vec_m(mvals, n, c):
    """S(m,n;c) for an array of m (same c)."""
    if c == 1:
        return np.ones(len(mvals))
    x, inv = units_and_inverses(c)
    mv = np.asarray(mvals, np.int64) % c
    ph = (np.outer(mv, x) + (n % c) * inv[None, :]) % c
    return np.cos(TWO_PI * ph / c).sum(axis=1)


def ramanujan(c, n):
    return sum(d * mobius(c // d) for d in range(1, c + 1) if c % d == 0 and n % d == 0)


def mobius(n):
    if n == 1:
        return 1
    r, m, p = 1, n, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            r = -r
        p += 1
    if m > 1:
        r = -r
    return r


def petersson_delta(k, N, m, n, Cmax):
    """Delta_{k,N}(m,n) = delta(m,n) + 2 pi i^{-k} sum_{c<=Cmax, N|c} S(m,n;c)/c J_{k-1}(4 pi sqrt(mn)/c)."""
    s = 0.0
    sgn = (-1) ** (k // 2)  # i^{-k} for k even
    for c in range(N, Cmax + 1, N):
        s += kloosterman(m, n, c) / c * jv(k - 1, 4 * math.pi * math.sqrt(m * n) / c)
    return (1.0 if m == n else 0.0) + 2 * math.pi * sgn * s


def rootnumber_petersson(k, N, n, Mmax):
    """sum_f omega_f eps_f lambda_f(n) via eps_f = i^k mu(N) lambda_f(N) sqrt N and lambda_f(N)lambda_f(n) = lambda_f(Nn):
    = i^k mu(N) sqrt N * Delta_{k,N}(Nn,1) = -2 pi N^{-1/2} sum_m S(Nn,1;Nm)/m J_{k-1}(4 pi sqrt(n/N)/m)  (N prime)."""
    s = 0.0
    for m in range(1, Mmax + 1):
        s += kloosterman(N * n, 1, N * m) / m * jv(k - 1, 4 * math.pi * math.sqrt(n / N) / m)
    return -2 * math.pi / math.sqrt(N) * s


def primes_upto(n):
    s = np.ones(n + 1, bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


# ------------------------------------------------------------------ PARI family data
def family_data(pari, N, k, P, Msym=2500):
    """All newforms of S_k(Gamma_0(N)) (N prime with S_k(1)=0, or N=1): list of dicts with
    lam[n] (n<=P, Deligne-normalised, real), eps (i^k mu(N) lam(N) sqrt N, N prime; i^k for N=1), LS1 = L(1,Sym^2 f)
    (PARI lfuncreate: motivic Vga=[2-k,0,1], weight 2k-1, conductor N^2, bad factor (1-N^{k-2}X)^{-1}), omega = 2pi^2/((k-1) N LS1)."""
    pari('mf=mfinit([%d,%d],0); B=mfeigenbasis(mf);' % (N, k))
    nb = int(pari('#B'))
    out = []
    Mc = max(P, Msym)
    for j in range(1, nb + 1):
        E = pari('mfembed(B[%d], mfcoefs(B[%d], %d))' % (j, j, Mc))
        if str(pari('type(mfembed(B[%d], mfcoefs(B[%d], 2))[1])' % (j, j))) != 't_VEC':
            embs = [E]
        else:
            embs = list(E)
        for e in embs:
            a = np.array([float(pari('real(%s)' % str(x))) if False else float(complex(x).real) for x in e])
            n = np.arange(len(a), dtype=float)
            lam = np.zeros_like(a)
            lam[1:] = a[1:] / n[1:] ** ((k - 1) / 2)
            if N == 1:
                eps = (-1) ** (k // 2)
            else:
                eps = int(round((-1) ** (k // 2) * (-1) * lam[N] * math.sqrt(N)))
            # Sym^2 coefficients (motivic) via Euler product, numeric
            M = Msym
            co = np.ones(M + 1)
            pr = primes_upto(M)
            for p in pr:
                p = int(p)
                if p == N:
                    # (1 - N^{k-2} X)^{-1}
                    c = [float(N) ** ((k - 2) * e_) for e_ in range(0, 40)]
                else:
                    ap = a[p]
                    q = float(p) ** (k - 1)
                    poly = [1.0, -(ap * ap - q), q * (ap * ap - q), -q ** 3]
                    # series of 1/poly
                    L = int(math.log(M) / math.log(p)) + 2
                    c = [1.0] + [0.0] * L
                    for i_ in range(1, L + 1):
                        s = 0.0
                        for d_ in range(1, 4):
                            if i_ - d_ >= 0:
                                s -= poly[d_] * c[i_ - d_]
                        c[i_] = s
                pk = p
                e_ = 1
                while pk <= M:
                    # multiply co[m] for m with v_p(m) = e_ by c[e_]: do it by stepping
                    pk *= p
                    e_ += 1
                # apply: for each m, v_p(m)
                mm = np.arange(p, M + 1, p)
                v = np.ones(len(mm), int)
                t_ = mm // p
                while True:
                    sel = (t_ % p == 0)
                    if not sel.any():
                        break
                    v[sel] += 1
                    t_[sel] //= p
                co[mm] *= np.array([c[x] for x in v])
            pari.set_real_precision(30)
            vec = pari('[' + ','.join('%.17g' % x for x in co[1:]) + ']')
            Lsym = pari.lfuncreate([vec, 0, [2 - k, 0, 1], 2 * k - 1, N * N, 1])
            feq = int(pari.lfuncheckfeq(Lsym))
            LS1 = float(pari.lfun(Lsym, k).real())
            omega = 2 * math.pi ** 2 / ((k - 1) * N * LS1)
            out.append(dict(lam=lam[:P + 1], eps=eps, LS1=LS1, omega=omega, feq=feq, a=a[:P + 1]))
    return out


def c_f(lam, N, P):
    """c_f = B_f / L(1,Sym^2 f): returns B_f (truncated Euler product over good p <= P; bad factor (1-1/N)^{-1})."""
    pr = primes_upto(P)
    pr = pr[pr != N]
    lp = lam[pr]
    e1 = lp ** 2 - 1
    p = pr.astype(float)
    num = 1 - e1 / p                      # may be <= 0 at p = 2 (lambda_2^2 > 3): c_f can be negative
    den = 1 - e1 / p + e1 / p ** 2 - 1 / p ** 3   # = (1 - 1/p)|1 - alpha^2/p|^2 > 0
    sgn = np.prod(np.sign(num))
    lnB = np.sum(np.log(np.abs(num)) - np.log(den))
    if N > 1:
        lnB += -math.log1p(-1 / N)
    return float(sgn * math.exp(lnB))


def stat_weights(k, N, n, Hfun):
    """Exact-point leading stationary phase weight of Proposition 1's n-term (without Lambda(n)):
    -(1/2pi) n^{-1/2} (2pi/phi'')^{1/2} e^{i(phi(t_n)+pi/4)} H(t_n), phi(t) = 2theta(t) - t ln n, 2theta'(t_n) = ln n."""
    import mpmath as mp
    n = np.asarray(n, float)
    out = np.zeros(len(n), complex)
    for i, nv in enumerate(n):
        ln_n = math.log(nv)
        t = 2 * math.pi * math.sqrt(nv / N)
        for _ in range(30):
            f = math.log(N) - 2 * math.log(2 * math.pi) + 2 * float(mp.re(mp.digamma(k / 2 + 1j * t))) - ln_n
            fp = -2 * float(mp.im(mp.psi(1, k / 2 + 1j * t)))
            dt = f / fp
            t -= dt
            if abs(dt) < 1e-12:
                break
        th = t / 2 * math.log(N) - t * math.log(2 * math.pi) + float(mp.im(mp.loggamma(k / 2 + 1j * t)))
        ph = 2 * th - t * ln_n
        php = -2 * float(mp.im(mp.psi(1, k / 2 + 1j * t)))
        out[i] = -(1 / (2 * math.pi)) * nv ** -0.5 * math.sqrt(2 * math.pi / php) * np.exp(1j * (ph + math.pi / 4)) * Hfun(np.array([t]))[0]
    return out
