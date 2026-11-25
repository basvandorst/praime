"""smlib: second moment of the one-point sum over the Petersson family (N prime, k in {2,4,6,8,10,14}).

Z_f(H) = sum_m beta_m lambda_f(m) + E_0 (Proposition 1 of the note), with
  w_{p^j} := -(1/2pi) (p^j)^{-1/2} I_{p^j} ln p,   I_u := int H(t) exp(i(2 theta(t) - t ln u)) dt,
  beta_{p^j} = w_{p^j} - w_{p^{j+2}}  (p != N),  beta_{N^j} = w_{N^j},  beta_1 = -sum_{p != N} w_{p^2}
(Lambda_f(p^j) = (lambda_f(p^j) - lambda_f(p^{j-2})) ln p for p != N, lambda_f(N)^j ln N for p = N).
Family sums of quadratic forms in lambda_f by Petersson:
  sum^h lambda_f(m) lambda_f(n) = delta_{mn} + 2 pi i^{-k} sum_{c = 0 (N)} S(m,n;c)/c J_{k-1}(4 pi sqrt(mn)/c).
Off-diagonal modulus c = N*mm, (mm, N) = 1: bilinear form sum_{m,n} x_m y_n S(m,n;c) J(...) evaluated as
  sum_r lam_r sum_{x mod c}^* Xhat_r(x) Yhat_r(xbar),  Xhat_r(x) = sum_n x_n F_r(n) e(n x/c),
with a low-rank (Chebyshev + eigen) factorisation J_{k-1}(s t/(pi mm)) = sum_r lam_r F_r(s) F_r(t), s = 2 pi sqrt(n/N).
"""
import math
import numpy as np
from scipy.special import jv, loggamma

TWO_PI = 2 * math.pi


def primes_upto(n):
    s = np.ones(n + 1, bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def theta(k, N, t):
    t = np.asarray(t, float)
    return t / 2 * math.log(N) - t * math.log(TWO_PI) + np.imag(loggamma(k / 2 + 1j * t))


def Wfun(u):
    u = np.asarray(u, float)
    z = 1 - (2 * u - 3) ** 2
    out = np.zeros_like(u)
    sel = z > 0
    out[sel] = np.exp(1 - 1 / z[sel])
    return out


class Window:
    """H = W(./T) * g_1 on a grid; I_u by trapezoid quadrature (H has Gaussian tails, integrand analytic)."""

    def __init__(self, N, k, T, h=0.02, pad=8.0):
        self.N, self.k, self.T, self.h = N, k, T, h
        self.tt = np.arange(T - pad, 2 * T + pad + h / 2, h)
        ss = np.arange(T, 2 * T, 0.002)
        Ws = Wfun(ss / T)
        # H(t) = int W(s/T) e^{-(t-s)^2}/sqrt(pi) ds
        H = np.empty(len(self.tt))
        for i in range(0, len(self.tt), 500):
            d = self.tt[i:i + 500, None] - ss[None, :]
            H[i:i + 500] = (np.exp(-d * d) * Ws[None, :]).sum(axis=1) * 0.002 / math.sqrt(math.pi)
        self.H = H
        self.intH = float(np.sum(H) * h)
        self.intH2 = float(np.sum(H * H) * h)
        self.L = self.intH / TWO_PI
        self.A0 = H * np.exp(2j * theta(k, N, self.tt)) * h

    def I(self, u):
        u = np.asarray(u, float)
        out = np.empty(len(u), complex)
        lu = np.log(u)
        for i in range(0, len(u), 1000):
            out[i:i + 1000] = np.exp(-1j * np.outer(lu[i:i + 1000], self.tt)) @ self.A0
        return out


def prime_power_coeffs(win, nlo, nhi, model='primes'):
    """beta as dict-free arrays: (n, beta_n) for the index set; n = 1 included (index 0) for model 'primes'.
    model 'integers': Y_f = sum_n -(1/2pi) n^{-1/2} I_n lambda_f(n) over all n in [nlo, nhi] (control model)."""
    N = win.N
    if model == 'integers':
        n = np.arange(max(2, nlo), nhi + 1, dtype=np.int64)
        b = -(1 / TWO_PI) * n ** -0.5 * win.I(n.astype(float))
        return n, b
    pr = primes_upto(nhi)
    ns, ps, js = [], [], []
    for p in pr:
        p = int(p)
        q, j = p, 1
        while q <= nhi:
            if q >= nlo or (j >= 3 and q * p * p >= nlo):
                ns.append(q); ps.append(p); js.append(j)
            q *= p
            j += 1
    ns = np.array(ns, np.int64); ps = np.array(ps, np.int64); js = np.array(js)
    w = -(1 / TWO_PI) * ns ** -0.5 * win.I(ns.astype(float)) * np.log(ps)
    w[(ns < nlo)] = 0.0
    wd = dict(zip(ns.tolist(), w.tolist()))
    beta = w.copy()
    b1 = 0.0 + 0.0j
    for i in range(len(ns)):
        p, j = int(ps[i]), int(js[i])
        if p == N:
            continue
        beta[i] -= wd.get(int(ns[i]) * p * p, 0.0)
        if j == 2:
            b1 -= w[i]
    n = np.concatenate([[1], ns])
    b = np.concatenate([[b1], beta])
    keep = (b != 0)
    return n[keep], b[keep]


# ---------------------------------------------------------------------------- Kloosterman / Petersson tools
def inv_table(c):
    """inverse mod c for all residues (0 where not a unit)."""
    x = np.arange(c, dtype=np.int64)
    inv = np.zeros(c, np.int64)
    if c == 1:
        return inv
    g = np.gcd(x, c)
    u = x[g == 1]
    # Euler: u^(phi(c)-1)
    ph = int(np.count_nonzero(g == 1))
    r = np.ones_like(u)
    b = u.copy()
    e = ph - 1
    while e:
        if e & 1:
            r = (r * b) % c
        b = (b * b) % c
        e >>= 1
    inv[u] = r
    return inv


class Moduli:
    def __init__(self, N):
        self.N = N
        self.invN = inv_table(N)
        self._invm = {}

    def invm(self, m):
        if m not in self._invm:
            self._invm[m] = inv_table(m)
        return self._invm[m]

    def units_and_bar(self, m):
        """for c = N m with (m,N)=1: arrays x (units mod c) and xbar."""
        N = self.N
        c = N * m
        if m % N == 0:
            inv = inv_table(c)
            x = np.nonzero(inv)[0].astype(np.int64)
            return x, inv[x]
        x = np.arange(c, dtype=np.int64)
        xN = x % N
        if m == 1:
            sel = xN != 0
            x = x[sel]
            return x, self.invN[x]
        xm = x % m
        im = self.invm(m)
        sel = (xN != 0) & (np.gcd(xm, m) == 1)
        x = x[sel]
        xN = xN[sel]; xm = xm[sel]
        # CRT: xbar = a mod N, b mod m -> a*m*(m^-1 mod N) + b*N*(N^-1 mod m)
        eN = m * int(pow(m, -1, N)) % c
        em = N * int(pow(N, -1, m)) % c
        xbar = (self.invN[xN] * eN + im[xm] * em) % c
        return x, xbar


def kloost(a, b, c):
    if c == 1:
        return 1.0
    inv = inv_table(c)
    x = np.arange(c, dtype=np.int64)
    sel = inv > 0
    return float(np.cos(TWO_PI * ((a * x[sel] + b * inv[sel]) % c) / c).sum())


def cheb_nodes(a, b, n):
    j = np.arange(n)
    return 0.5 * (a + b) + 0.5 * (b - a) * np.cos(math.pi * (j + 0.5) / n)


def bary_matrix(nodes, s):
    """Lagrange basis at Chebyshev (1st kind) nodes, barycentric form: returns len(s) x len(nodes)."""
    n = len(nodes)
    j = np.arange(n)
    wts = (-1.0) ** j * np.sin(math.pi * (j + 0.5) / n)
    d = s[:, None] - nodes[None, :]
    hit = np.abs(d) < 1e-14
    d[hit] = 1e-14
    M = wts[None, :] / d
    M /= M.sum(axis=1, keepdims=True)
    return M


def kernel_factor(k, mm, nodes, tol=1e-12):
    """eigen-factorisation of K(s,t) = J_{k-1}(s t/(pi mm)) on the node grid; returns (lam, E) with |lam| > tol*max."""
    K = jv(k - 1, np.outer(nodes, nodes) / (math.pi * mm))
    lam, E = np.linalg.eigh(K)
    keep = np.abs(lam) > tol * np.abs(lam).max()
    return lam[keep], E[:, keep]


def offdiag_modulus(mod, mm, n, xw, yw, Fr, lam, conj_second=False):
    """sum_{n,n'} xw_n yw_n' S(n,n';c) sum_r lam_r F_r(n) F_r(n') with c = N mm, (mm,N)=1.
    If conj_second: the second factor uses yw with e(n' xbar/c) (same) -- caller passes yw already conjugated."""
    c = mod.N * mm
    x, xbar = mod.units_and_bar(mm)
    res = n % c
    tot = 0.0 + 0.0j
    for r in range(len(lam)):
        A = np.bincount(res, weights=(xw * Fr[:, r]).real, minlength=c) + 1j * np.bincount(res, weights=(xw * Fr[:, r]).imag, minlength=c)
        Ah = np.fft.ifft(A) * c          # sum_b A(b) e(b x/c)
        if yw is xw:
            Bh = Ah
        else:
            B = np.bincount(res, weights=(yw * Fr[:, r]).real, minlength=c) + 1j * np.bincount(res, weights=(yw * Fr[:, r]).imag, minlength=c)
            Bh = np.fft.ifft(B) * c
        tot += lam[r] * np.sum(Ah[x] * Bh[xbar])
    return tot
