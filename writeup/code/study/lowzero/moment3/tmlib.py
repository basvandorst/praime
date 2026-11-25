"""tmlib.py: the CFKRS/Heap recipe for the mean square of a product of k = 1, 2, 3 distinct real primitive Dirichlet
L-functions, and the comparators.  Pure Python + numpy + cypari2 (L-values at s = 1 only).

Notation (derivation.md):
  chi_j = Kronecker symbol of the fundamental discriminant D_j, q_j = |D_j|, X_j(t) = q_j t / (2 pi), L_j = ln X_j.
  Recipe (Heap's modification: swaps per L-function): with alpha_j = beta_j = z_j,
     P(t) = CT_z [ prod_j (2 zeta(1 + 2 z_j) X_j^{z_j}) h(z) ] = sum_{m subset of {1..k}} h[m] prod_{j not in m} (L_j + 2 gamma)
  where h(z) = prod_j c_j(2 z_j) * prod_{i<j} Lam_ij(1 + z_i + z_j)^2 * A(z),
     c_j(x) = prod_{p | q_j} (1 - p^{-1-x}),
     Lam_ij(s) = sum_n chi_i(n) chi_j(n) n^{-s}  (imprimitive if (q_i, q_j) > 1),
     A(z) = prod_p A_p,  A_p = R_p(x) prod_{i,j} (1 - x_i x_j),  x_j = chi_j(p) p^{-1/2 - z_j},  R_p = sum_n h_n(x)^2,
  and h[m] = the coefficient of prod_{j in m} z_j in h (multilinear jet).
"""
import math, itertools, functools
import numpy as np
import cypari2

pari = cypari2.Pari()
pari.allocatemem(300 * 10**6)
pari.set_real_precision(38)

EULER_GAMMA = 0.57721566490153286061
PMAX = 10**6                      # explicit Euler product of A up to PMAX, quartic principal tail beyond (derivation.md 2.4)


# ----------------------------------------------------------------------------------------------- truncated algebras
class Alg:
    """Truncated polynomial algebra over R in nv variables; monomials given by a list of exponent tuples
    (closed under division).  Elements: numpy arrays of shape (..., M)."""

    def __init__(self, monos):
        self.monos = [tuple(m) for m in monos]
        self.M = len(self.monos)
        self.idx = {m: i for i, m in enumerate(self.monos)}
        self.nv = len(self.monos[0])
        self.table = []
        for i, a in enumerate(self.monos):
            for j, b in enumerate(self.monos):
                c = tuple(x + y for x, y in zip(a, b))
                if c in self.idx:
                    self.table.append((i, j, self.idx[c]))
        self.deg = max(sum(m) for m in self.monos)       # nilpotency: n^(deg+1) = 0 for n with zero constant term

    def const(self, c, shape=()):
        out = np.zeros(shape + (self.M,))
        out[..., 0] = c
        return out

    def var(self, k, shape=()):
        out = np.zeros(shape + (self.M,))
        e = tuple(1 if i == k else 0 for i in range(self.nv))
        out[..., self.idx[e]] = 1.0
        return out

    def mul(self, a, b):
        a, b = np.broadcast_arrays(a, b)
        out = np.zeros(a.shape)
        for i, j, k in self.table:
            out[..., k] += a[..., i] * b[..., j]
        return out

    def series(self, a, coefs):
        """sum_n coefs[n] * (a - a0)^n, coefs a list of arrays broadcastable to a[..., 0]."""
        n = a.copy(); n[..., 0] = 0.0
        out = np.zeros(a.shape); out[..., 0] = coefs[0]
        pw = None
        for k in range(1, min(len(coefs), self.deg + 1)):
            pw = n if pw is None else self.mul(pw, n)
            out += np.asarray(coefs[k])[..., None] * pw
        return out

    def log(self, a):
        a0 = a[..., 0]
        u = a / a0[..., None]
        co = [np.log(a0)] + [((-1.0) ** (k + 1)) / k * np.ones_like(a0) for k in range(1, self.deg + 1)]
        return self.series(u, co)

    def exp(self, a):
        a0 = a[..., 0]
        co = [np.exp(a0) / math.factorial(k) for k in range(self.deg + 1)]
        return self.series(a, co)

    def exp_scaled(self, a, scale):
        """exp(scale * a) for array scale broadcast over leading dims."""
        return self.exp(a * np.asarray(scale)[..., None])

    def compose_taylor(self, s, coefs):
        """f(s) = sum coefs[k] s^k for a jet s with zero constant term (coefs scalars)."""
        out = self.const(coefs[0], s.shape[:-1])
        pw = None
        for k in range(1, min(len(coefs), self.deg + 1)):
            pw = s if pw is None else self.mul(pw, s)
            out = out + coefs[k] * pw
        return out


def multilinear(k):
    return Alg(list(itertools.product((0, 1), repeat=k)))


def univariate(deg):
    return Alg([(d,) for d in range(deg + 1)])


# ----------------------------------------------------------------------------------------------- characters
def is_fundamental(D):
    return D != 1 and bool(pari.isfundamental(D))


@functools.lru_cache(maxsize=None)
def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0].astype(np.int64)


def powmod_vec(b, e, m):
    """b^e mod m elementwise (int64; m < 3e9 so products < 2^63)."""
    b = b % m; r = np.ones_like(b); e = e.copy()
    while np.any(e > 0):
        odd = (e & 1) == 1
        r = np.where(odd, (r * b) % m, r)
        b = (b * b) % m
        e >>= 1
    return r


@functools.lru_cache(maxsize=None)
def kron_primes(D, n=PMAX):
    """chi_D(p) for all primes p <= n (Kronecker symbol), int8 array aligned with primes_upto(n)."""
    P = primes_upto(n)
    out = np.zeros(len(P), dtype=np.int8)
    # p = 2
    if D % 2 == 0:
        out[0] = 0
    else:
        out[0] = 1 if D % 8 in (1, 7) else -1
    po = P[1:]
    Dm = np.int64(D) % po
    e = (po - 1) // 2
    leg = powmod_vec(Dm, e, po)
    v = np.where(Dm == 0, 0, np.where(leg == 1, 1, -1))
    out[1:] = v
    return out


def prime_factors(n):
    n = abs(int(n)); f = []; p = 2
    while p * p <= n:
        if n % p == 0:
            f.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        f.append(n)
    return f


@functools.lru_cache(maxsize=None)
def lam_taylor(Di, Dj, order=3):
    """Taylor coefficients l_0..l_order of log Lam(1 + s), Lam(s) = sum_n chi_Di(n) chi_Dj(n) n^{-s}."""
    Ds = int(pari.coredisc(Di * Dj))
    L = pari.lfuncreate(Ds)
    der = [pari.lfun(L, 1, k) for k in range(order + 1)]          # PARI: k-th derivative at s = 1
    der = [float(pari.real(d)) for d in der]
    # log L(1+s) Taylor from derivatives: L(1+s) = sum der[k] s^k / k!
    a = [der[k] / math.factorial(k) for k in range(order + 1)]
    lg = log_series(a, order)
    # Euler factors at p | Di Dj, p not dividing Ds: (1 - chi_Ds(p) p^{-1-s}); at p | Ds factor is 1.
    for p in prime_factors(Di * Dj):
        c = int(pari.kronecker(Ds, p))
        if c == 0:
            continue
        # log(1 - c p^{-1} e^{-s ln p}) Taylor in s
        lp = math.log(p)
        ser = [0.0] * (order + 1)
        # -sum_{m>=1} c^m p^{-m} e^{-m s ln p} / m ; coefficient of s^k: -sum_m c^m p^{-m} (-m lp)^k / (m k!)
        for k in range(order + 1):
            ser[k] = -sum((c ** m) * p ** (-m) * (-m * lp) ** k / (m * math.factorial(k)) for m in range(1, 200))
        lg = [x + y for x, y in zip(lg, ser)]
    return lg, Ds


def log_series(a, order):
    """Taylor coefficients of log(sum a_k s^k) up to order (a_0 > 0)."""
    b = [x / a[0] for x in a]
    # log(1 + u), u = sum_{k>=1} b_k s^k
    out = [math.log(a[0])] + [0.0] * order
    upow = [1.0] + [0.0] * order
    for m in range(1, order + 1):
        new = [0.0] * (order + 1)
        for i in range(order + 1):
            if upow[i] == 0.0:
                continue
            for k in range(1, order + 1 - i):
                new[i + k] += upow[i] * b[k]
        upow = new
        for i in range(order + 1):
            out[i] += ((-1) ** (m + 1)) / m * upow[i]
    return out


# ----------------------------------------------------------------------------------------------- h(z) as a jet
def jet_c(alg, D, xjet):
    """log c(x) = sum_{p | D} log(1 - p^{-1-x}) for a jet x (zero constant)."""
    out = alg.const(0.0)
    for p in prime_factors(D):
        u = alg.exp(xjet * (-math.log(p))) * (1.0 / p)
        out = out + alg.log(alg.const(1.0) - u)
    return out


def jet_logA(alg, Ds, zjets, pmax=PMAX, chunk=20000):
    """log A(z) = sum_p log A_p as a jet; zjets: list of k jets (the shifts z_j as jets)."""
    k = len(Ds)
    P = primes_upto(pmax)
    chis = [kron_primes(D, pmax).astype(float) for D in Ds]
    total = alg.const(0.0)
    bounds = sorted(set([0] + [int(np.searchsorted(P, thr)) for thr in (50, 2000, 10**5)]
                        + list(range(int(np.searchsorted(P, 10**5)), len(P), chunk)) + [len(P)]))
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        pp = P[lo:hi].astype(float)
        lp = np.log(pp)
        x = []
        for j in range(k):
            xj = alg.exp(zjets[j][None, :] * (-lp)[:, None]) * (chis[j][lo:hi] / np.sqrt(pp))[:, None]
            x.append(xj)
        # elementary symmetric jets
        if k == 1:
            e = [x[0]]
        elif k == 2:
            e = [x[0] + x[1], alg.mul(x[0], x[1])]
        else:
            e = [x[0] + x[1] + x[2], alg.mul(x[0], x[1]) + alg.mul(x[0], x[2]) + alg.mul(x[1], x[2]),
                 alg.mul(alg.mul(x[0], x[1]), x[2])]
        pmin = pp[0]
        nmax = 140 if pmin < 50 else (40 if pmin < 2000 else (16 if pmin < 10**5 else 10))
        hs = [alg.const(1.0, (len(pp),))]
        R = alg.const(1.0, (len(pp),))
        for n in range(1, nmax + 1):
            hn = alg.mul(e[0], hs[-1])
            if k >= 2 and n >= 2:
                hn = hn - alg.mul(e[1], hs[-2])
            if k >= 3 and n >= 3:
                hn = hn + alg.mul(e[2], hs[-3])
            hs.append(hn)
            if len(hs) > 4:
                hs.pop(0)
            R = R + alg.mul(hn, hn)
        lA = alg.log(R)
        for i in range(k):
            for j in range(k):
                lA = lA + alg.log(alg.const(1.0) - alg.mul(x[i], x[j]))
        total = total + lA.sum(axis=0)
    # principal quartic tail beyond pmax: log A_p ~ -sum_{i<j} p^{-2-2 z_i - 2 z_j}  (all p > pmax are good)
    if k >= 2:
        for i in range(k):
            for j in range(i + 1, k):
                s = (zjets[i] + zjets[j]) * 2.0           # jet of 2(z_i + z_j), zero constant
                # sum_{p > pmax} p^{-2} e^{-s ln p} = sum_n (-s)^n/n! S_n, S_n = sum_{p>pmax} (ln p)^n p^{-2}
                co = [-tail_S(nn, pmax) * ((-1.0) ** nn) / math.factorial(nn) for nn in range(alg.deg + 1)]
                total = total + alg.compose_taylor(s, co)
    return total


@functools.lru_cache(maxsize=None)
def tail_S(n, pmax):
    """S_n = sum_{p > pmax} (ln p)^n / p^2 ~ int_pmax^inf (ln x)^(n-1) x^-2 dx (prime number theorem; rel. error ~ 1/ln^2)."""
    import mpmath
    return float(mpmath.quad(lambda x: mpmath.log(x) ** (n - 1) / x ** 2, [pmax, 10 * pmax, mpmath.inf]))


def h_jet(alg, Ds, zjets, pmax=PMAX):
    """Jet of h(z) = prod_j c_j(2 z_j) prod_{i<j} Lam_ij(1 + z_i + z_j)^2 A(z)."""
    k = len(Ds)
    lh = alg.const(0.0)
    for j in range(k):
        lh = lh + jet_c(alg, Ds[j], zjets[j] * 2.0)
    for i in range(k):
        for j in range(i + 1, k):
            lt, _ = lam_taylor(Ds[i], Ds[j], 3)
            s = zjets[i] + zjets[j]
            lh = lh + 2.0 * alg.compose_taylor(s, lt)
    if k >= 2:
        lh = lh + jet_logA(alg, Ds, zjets, pmax)
    return alg.exp(lh)


# ----------------------------------------------------------------------------------------------- predictions
class Triple:
    """Recipe and comparators for a tuple of distinct fundamental discriminants (k = 1, 2, 3)."""

    def __init__(self, Ds, pmax=PMAX):
        self.Ds = tuple(int(D) for D in Ds)
        self.k = len(self.Ds)
        self.q = np.array([abs(D) for D in self.Ds], dtype=float)
        ml = multilinear(self.k)
        z = [ml.var(j) for j in range(self.k)]
        self.hm = h_jet(ml, self.Ds, z, pmax)             # multilinear coefficients h[m]
        self.ml = ml
        un = univariate(self.k)
        s = un.var(0)
        self.g = h_jet(un, self.Ds, [s * 0.5] * self.k, pmax)   # g(sigma) = h(sigma/2, ..., sigma/2), coefficients g_0..g_k

    def logs(self, t):
        t = np.asarray(t, dtype=float)
        return [np.log(qq * t / (2 * math.pi)) for qq in self.q]

    def P(self, t):
        """recipe density P(t)."""
        Ls = self.logs(t)
        out = 0.0
        for i, m in enumerate(self.ml.monos):
            term = self.hm[i]
            for j in range(self.k):
                if m[j] == 0:
                    term = term * (Ls[j] + 2 * EULER_GAMMA)
            out = out + term
        return out

    def P_terms(self, t):
        """the recipe split by the multilinear index m (which shifts were differentiated)."""
        Ls = self.logs(t)
        out = {}
        for i, m in enumerate(self.ml.monos):
            term = self.hm[i] * np.ones_like(Ls[0])
            for j in range(self.k):
                if m[j] == 0:
                    term = term * (Ls[j] + 2 * EULER_GAMMA)
            out[m] = term
        return out

    def lead(self, t):
        Ls = self.logs(t)
        out = self.hm[0]
        for L in Ls:
            out = out * L
        return out

    def diag(self, t):
        """AFE diagonal (sharp cutoff at sqrt(C), C = prod X_j): 2 Res_{s=0} zeta(1+s)^k g(s) C^{s/2} / s."""
        Ls = self.logs(t)
        lnC = sum(Ls)
        k = self.k
        # sigma*zeta(1+sigma) = 1 + g0' sigma + ... ; Stieltjes: zeta(1+s) = 1/s + gamma0 - gamma1 s + gamma2 s^2/2 - ...
        st = [EULER_GAMMA, -0.0728158454836767248605863758749547, -0.00969036319287231848453038603521]
        zs = [1.0, st[0], -st[1], st[2] / 2.0]                     # coefficients of sigma*zeta(1+sigma)
        # series of (sigma zeta)^k
        zk = [1.0] + [0.0] * k
        for _ in range(k):
            new = [0.0] * (k + 1)
            for a in range(k + 1):
                for b in range(k + 1 - a):
                    new[a + b] += zk[a] * zs[b]
            zk = new
        # (sigma zeta)^k * g(sigma): coefficients c_0..c_k
        c = [sum(zk[a] * self.g[n - a] for a in range(n + 1)) for n in range(k + 1)]
        # 2 * coefficient of sigma^k in c(sigma) * exp(sigma lnC / 2)
        out = 0.0
        for n in range(k + 1):
            out = out + c[n] * (lnC / 2.0) ** (k - n) / math.factorial(k - n)
        return 2.0 * out


def integrate(f, T1, T2, n=64):
    """Gauss-Legendre in u = ln t of f(t) dt (P is a polynomial in ln t)."""
    x, w = np.polynomial.legendre.leggauss(n)
    a, b = math.log(T1), math.log(T2)
    u = 0.5 * (b - a) * x + 0.5 * (b + a)
    t = np.exp(u)
    return float(np.sum(w * f(t) * t) * 0.5 * (b - a))


# ----------------------------------------------------------------------------------------------- comparators
def h_jet_no3(alg, Ds, zjets, pmax=PMAX):
    """h with the three-body Euler product replaced by the product of the pair factors A_ij(z_i + z_j)."""
    k = len(Ds)
    lh = alg.const(0.0)
    for j in range(k):
        lh = lh + jet_c(alg, Ds[j], zjets[j] * 2.0)
    for i in range(k):
        for j in range(i + 1, k):
            lt, _ = lam_taylor(Ds[i], Ds[j], 3)
            lh = lh + 2.0 * alg.compose_taylor(zjets[i] + zjets[j], lt)
            lh = lh + jet_logA(alg, (Ds[i], Ds[j]), [zjets[i], zjets[j]], pmax)
    return alg.exp(lh)


class Family:
    """All predictions for one k-tuple (k = 3): recipe, leading, AFE diagonal, Kirkwood (pairs/singles),
    independence (product of singles), no-3-body (A_3 -> prod A_ij)."""

    def __init__(self, Ds, pmax=10**5):
        self.Ds = tuple(Ds); self.k = len(Ds)
        self.full = Triple(Ds, pmax)
        self.singles = [Triple((D,), pmax) for D in Ds]
        self.pairs = {}
        if self.k == 3:
            for i in range(3):
                for j in range(i + 1, 3):
                    self.pairs[(i, j)] = Triple((Ds[i], Ds[j]), pmax)
            ml = self.full.ml
            z = [ml.var(j) for j in range(3)]
            self.hm_no3 = h_jet_no3(ml, self.Ds, z, pmax)

    def P_no3(self, t):
        Ls = self.full.logs(t); out = 0.0
        for i, m in enumerate(self.full.ml.monos):
            term = self.hm_no3[i]
            for j in range(self.k):
                if m[j] == 0:
                    term = term * (Ls[j] + 2 * EULER_GAMMA)
            out = out + term
        return out

    def P_ind(self, t):
        out = 1.0
        for s in self.singles:
            out = out * s.P(t)
        return out

    def P_kirk(self, t):
        out = 1.0
        for p in self.pairs.values():
            out = out * p.P(t)
        for s in self.singles:
            out = out / s.P(t)
        return out

    def predictions(self, T1=100.0, T2=1000.0, w=None):
        """integrals over [T1, T2] (weight w(t) optional) of every comparator density."""
        f = {'full': self.full.P, 'lead': self.full.lead, 'diag': self.full.diag}
        if self.k == 3:
            f.update({'kirk': self.P_kirk, 'ind': self.P_ind, 'no3': self.P_no3})
        out = {}
        for name, fn in f.items():
            g = fn if w is None else (lambda t, fn=fn: fn(t) * w(t))
            out[name] = integrate(g, T1, T2, n=96) if w is None else integrate_w(g, T1, T2)
        return out


def integrate_w(f, T1, T2, n=4000):
    """composite Gauss-Legendre (weights with kinks need panels)."""
    edges = np.linspace(T1, T2, 41)
    x, wq = np.polynomial.legendre.leggauss(24)
    tot = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        t = 0.5 * (b - a) * x + 0.5 * (b + a)
        tot += float(np.sum(wq * f(t)) * 0.5 * (b - a))
    return tot


# ----------------------------------------------------------------------------------------------- fourth moment (theorem)
def P4_coeffs(D, N=64, r=0.2):
    """Coefficients c_0..c_4 of Topacogullari's (arXiv 1909.11517, eq. (5508)) polynomial for |L(1/2+it, chi_D)|^4:
    P(log t) = Res_{z=0} q^z psi_z(q)^6 / (psi_0 psi_2z psi_{1+2z}) zeta(1+z)^6 / ((2pi)^z zeta(1+2z) zeta(2+2z)) t^z
    = sum_k c_k (log t)^k  with c_k = Res [F(z) z^k / k!]."""
    import mpmath as mp
    mp.mp.dps = 30
    q = abs(D); fac = prime_factors(q)
    def psi(z):
        out = mp.mpf(1)
        for p in fac:
            out *= (1 - mp.mpf(p) ** (-1 - z))
        return out
    c = [mp.mpf(0)] * 5
    for kk in range(N):
        z = r * mp.exp(2j * mp.pi * (kk + 0.5) / N)
        F = (mp.mpf(q) ** z * psi(z) ** 6 / (psi(0) * psi(2 * z) * psi(1 + 2 * z))
             * mp.zeta(1 + z) ** 6 / ((2 * mp.pi) ** z * mp.zeta(1 + 2 * z) * mp.zeta(2 + 2 * z)))
        for k in range(5):
            c[k] += F * z ** k / mp.factorial(k) * z / N
    return [float(mp.re(x)) for x in c]
