"""ecross library: the swapped (Gram-lattice) cross term between the zeros of two elliptic-curve L-functions.

Normalisation (as round 144): L(s,E) = sum lambda(n) n^-s, lambda(n) = a_n/sqrt n; Lambda(s) = N^{s/2} Gamma_C(s+1/2) L(s,E)
= w Lambda(1-s); zeros 1/2 + i gamma; theta_E(t) = (t/2) ln N - t ln 2 pi + Im lnGamma(1 + i t), so that
L(1/2+it) = w e^{-2 i theta_E(t)} L(1/2-it).

One-point law (recipe):  sum_gamma e^{2 i theta_E(gamma)} ~ -w c_E (T2 - T1)/2pi,
   c_E = lim_{x->0} F(1+x)/x,  F(s) = sum lambda(n) mu_E(n) n^-s  (mu_E: coefficients of 1/L(s,E))
       = B_E / L(1, Sym^2 E),  B_E = prod_good B_p * prod_{p|N} (1-1/p)^-1,
   B_p = (1 - e1/p)/(1 - e1/p + e1/p^2 - 1/p^3),  e1 = lambda_p^2 - 1   (absolutely convergent).
Lines of the one-point sum (recipe): sum_gamma e^{2 i theta(gamma)} r^{-i gamma} ~ -w c_E rho(r)/sqrt(ab) (T2-T1)/2pi for
   r = a/b in lowest terms, rho multiplicative, rho(p^v) = p^{|v|/2} S_p(v)/F_p(1/p),
   S_p(v) = sum_{m-h=v} lambda(p^m) mu(p^h) p^{-(m+h)/2},  F_p(X) = sum_j lambda(p^j) mu(p^j) X^j.
Two-curve law (recipe, double swap):  E X_12 = (1/4 pi^2) w1 w2 c1 c2 J_12,  J_12 = int H_12(u) K(u) du,
   H_12(u) = prod_p h_p(iu),  h_p(x) = Y_p(x)/(F_1p(1/p) F_2p(1/p)) = sum_v rho_1(p^v) rho_2(p^v) p^{-|v|} p^{x v},
   Y_p(x) = sum_{m+k=n+h} lambda_1(p^m) lambda_2(p^n) mu_1(p^h) mu_2(p^k) p^{-(m+n+h+k)/2} p^{x(m-h)}
          = 1 + lambda_1 lambda_2 (p^{x-1} + p^{-x-1}) - (lambda_1^2 + lambda_2^2)/p + O(p^-2),
   so H_12(u) = |L(1+iu, E1 x E2)|^2 A_12(iu), A_12 absolutely convergent.
"""
import math
import numpy as np
from scipy.special import loggamma

TWO_PI = 2 * math.pi


def theta_E(N, t):
    t = np.asarray(t, float)
    return t / 2 * math.log(N) - t * math.log(TWO_PI) + np.imag(loggamma(1 + 1j * t))


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]


# ------------------------------------------------------------------ curves (PARI)
def curve_info(pari, ainv, P=10 ** 5):
    """N, w, bad primes, a_p for all p <= P (dict p -> a_p as array aligned with primes_upto(P)), L(1, Sym^2 E)."""
    E = pari.ellinit(ainv)
    N = int(pari.ellglobalred(E)[0]); fa = pari.factor(N)
    bad = [int(p) for p in fa[0]]
    assert all(int(e) == 1 for e in fa[1]), 'semistable only'
    w = int(pari.ellrootno(E))
    an = pari.ellan(E, P); pr = primes_upto(P)
    ap = np.array([int(an[int(p) - 1]) for p in pr], dtype=np.int64)
    LS1 = complex(pari.lfun(pari.lfunsympow(E, 2), 2))      # analytic L(1, Sym^2 E) = PARI value at s = 2
    assert abs(LS1.imag) < 1e-12
    return dict(N=N, w=w, bad=bad, primes=pr, ap=ap, LS1=LS1.real, ainv=list(ainv))


def c_E(info, P=None):
    """one-point constant c_E = B_E/L(1,Sym^2 E); also the direct (conditionally convergent) partial product
    prod_{good p<=P} (1 - lambda(p^2)/p) prod_{p|N} (1 + 1/p) as a check."""
    pr = info['primes'].astype(float); ap = info['ap'].astype(float)
    if P is not None:
        m = pr <= P; pr, ap = pr[m], ap[m]
    good = ~np.isin(pr.astype(np.int64), info['bad'])
    lam2 = ap[good] ** 2 / pr[good]; e1 = lam2 - 1; p = pr[good]
    lnB = np.sum(np.log1p(-e1 / p) - np.log1p(-e1 / p + e1 / p ** 2 - 1 / p ** 3))
    lnB += sum(-math.log1p(-1 / q) for q in info['bad'])
    B = math.exp(lnB)
    direct = math.exp(np.sum(np.log1p(-e1 / p)) + sum(math.log1p(1 / q) for q in info['bad']))
    return B / info['LS1'], B, direct


# ------------------------------------------------------------------ local factors
def lam_pows(p, a, bad, jmax):
    """lambda(p^j), j = 0..jmax (array over curves: a and bad are arrays), mu(p^h) h = 0..2, F_p(1/p)."""
    a = np.asarray(a, float); bad = np.asarray(bad, bool)
    lam = a / math.sqrt(p)
    L = np.zeros(a.shape + (jmax + 1,))
    L[..., 0] = 1.0
    if jmax >= 1: L[..., 1] = lam
    for j in range(2, jmax + 1):
        L[..., j] = np.where(bad, lam * L[..., j - 1], lam * L[..., j - 1] - L[..., j - 2])
    mu = np.stack([np.ones_like(lam), -lam, np.where(bad, 0.0, 1.0)], axis=-1)
    F = np.where(bad, 1 - lam ** 2 / p, 1 - lam ** 2 / p + (lam ** 2 - 1) / p ** 2)
    return L, mu, F


def jmax_for(p, tol=1e-17):
    j = 1
    while p ** (-j) * (j + 3) ** 2 > tol: j += 1
    return j + 2


def h_table(p, a1, b1, a2, b2, x):
    """h_p(x) = Y_p(x)/(F_1 F_2) for every combination of the value lists (a1, b1) (curve 1: a_p, bad flag) and
    (a2, b2), at the points x (complex array). Returns array [len(a1), len(a2), len(x)]."""
    J = jmax_for(p)
    L1, mu1, F1 = lam_pows(p, a1, b1, J + 2)
    L2, mu2, F2 = lam_pows(p, a2, b2, J + 2)
    x = np.asarray(x, complex)
    Y = np.zeros((len(a1), len(a2), len(x)), complex)
    n = np.arange(J + 1)
    Xn = np.exp(-np.outer(1 - x, n) * math.log(p))                  # [nx, n]: p^{-n(1-x)}
    for h in range(3):
        for k in range(3):
            d = h - k
            # lambda_1(p^{n+d}) for n >= max(0, -d)
            idx = n + d; ok = idx >= 0
            A = np.where(ok[None, :], L1[:, np.clip(idx, 0, None)], 0.0)       # [m1, n]
            B = L2[:, :J + 1]                                                   # [m2, n]
            S = np.einsum('in,jn,un->iju', A, B, Xn)
            Y += (mu1[:, h][:, None, None] * mu2[:, k][None, :, None] * p ** (-h)) * np.exp(-k * x * math.log(p))[None, None, :] * S
    return Y / (F1[:, None, None] * F2[None, :, None])


def S_line(p, a, bad, v):
    """S_p(v) = sum_{m-h=v} lambda(p^m) mu(p^h) p^{-(m+h)/2} for arrays a, bad (curves)."""
    J = jmax_for(p) + abs(v) + 2
    L, mu, F = lam_pows(p, a, bad, J)
    out = np.zeros(np.shape(a))
    for h in range(3):
        m = h + v
        if m < 0: continue
        # sum over m' = m + j, h' = h + j is not the structure: mu has support h <= 2, so S(v) = sum_h lambda(p^{h+v}) mu(p^h) p^{-(2h+v)/2}
        out = out + L[..., m] * mu[..., h] * p ** (-(2 * h + v) / 2)
    return out, F


def rho_pv(p, a, bad, v):
    """rho(p^v) = p^{|v|/2} S_p(v)/F_p(1/p)."""
    S, F = S_line(p, a, bad, v)
    return p ** (abs(v) / 2) * S / F


# ------------------------------------------------------------------ H for all pairs of a set
def S_vec(p, a, bad, tol=1e-15):
    """S_p(v) for v = -2..V (arrays over curves) and F_p(1/p); V chosen so that p^{-v/2}(v+3) < tol."""
    V = 2
    while p ** (-V / 2) * (V + 3) > tol: V += 1
    S = np.stack([S_line(p, a, bad, v)[0] for v in range(-2, V + 1)], axis=-1)
    F = S_line(p, a, bad, 0)[1]
    return S, F, np.arange(-2, V + 1)


def h_lines(p, a1, b1, a2, b2, x):
    """h_p(x) = sum_v S_1(v) S_2(v) p^{x v}/(F_1 F_2) (the lines form; must equal h_table)."""
    S1, F1, v = S_vec(p, a1, b1); S2, F2, _ = S_vec(p, a2, b2)
    E = np.exp(np.outer(np.asarray(x, complex), v) * math.log(p))          # [nx, nv]
    return np.einsum('iv,jv,uv->iju', S1, S2, E) / (F1[:, None, None] * F2[None, :, None])


def lnH_matrix(curves, u, P0=3000, P=10 ** 5, curves2=None, reverse=False, second=True):
    """ln H_ij(u) for curves i (list curves) and j (list curves2, default the same): exact local factors (lines form)
    for p <= P0; for P0 < p <= P (good primes; bad primes are always exact) the first-order lines
    2 lambda_i lambda_j cos(u ln p)/p (reverse: with the opposite sign) plus, if second, the non-oscillating second-order
    term -lambda_i^2 lambda_j^2/p^2 (the oscillating O(p^-2) terms beyond P0 are dropped: they sit at frequencies
    >= ln P0 where the Gaussian kernel is negligible).  reverse=True pairs the line r of curve i with the line 1/r of
    curve j (the 'reversed direction' alternative).  Returns [n1, n2, len(u)]."""
    c2 = curves if curves2 is None else curves2
    x = 1j * np.asarray(u, float)
    pr = curves[0]['primes']
    A1 = np.array([c['ap'] for c in curves], dtype=float); A2 = np.array([c['ap'] for c in c2], dtype=float)
    allbad = sorted(set(b for c in list(curves) + list(c2) for b in c['bad']))
    lnH = np.zeros((len(curves), len(c2), len(u)), complex)
    k0 = int(np.searchsorted(pr, P0, 'right')); k1 = int(np.searchsorted(pr, P, 'right'))
    exact = list(range(k0)) + [int(np.searchsorted(pr, b)) for b in allbad if b > P0]   # bad primes > P0: exact too
    for k in exact:
        p = int(pr[k])
        b1 = np.array([p in c['bad'] for c in curves]); b2 = np.array([p in c['bad'] for c in c2])
        code1 = (A1[:, k] * 2 + b1).astype(np.int64); code2 = (A2[:, k] * 2 + b2).astype(np.int64)
        allv = np.unique(np.concatenate([code1, code2]))
        av = (allv // 2).astype(float); bv = (allv % 2).astype(bool)
        if reverse:
            S, F, v = S_vec(p, av, bv)
            Sr = np.zeros_like(S); nv = len(v)                          # S(-v): defined for -v >= -2, i.e. v <= 2
            for a_, vv in enumerate(v):
                if -vv in v: Sr[:, a_] = S[:, list(v).index(-vv)]
            E = np.exp(np.outer(x, v) * math.log(p))
            tab = np.log(np.einsum('iv,jv,uv->iju', S, Sr, E) / (F[:, None, None] * F[None, :, None]))
        else:
            tab = np.log(h_lines(p, av, bv, av, bv, x))
        i1 = np.searchsorted(allv, code1); i2 = np.searchsorted(allv, code2)
        lnH += tab[i1][:, i2]
    if k1 > k0:
        cols = np.array([k for k in range(k0, k1) if int(pr[k]) not in allbad], dtype=np.int64)
        p = pr[cols].astype(float); lp = np.log(p)
        l1 = A1[:, cols] / np.sqrt(p)[None, :]; l2 = A2[:, cols] / np.sqrt(p)[None, :]
        sg = -1.0 if reverse else 1.0
        for iu, uu in enumerate(u):
            wv = sg * 2 * np.cos(uu * lp) / p
            lnH[:, :, iu] += (l1 * wv[None, :]) @ l2.T
        if second:
            lnH += (-((l1 ** 2) / p[None, :] ** 2) @ (l2 ** 2).T)[:, :, None]
    return lnH


def J_from_H(H, u, sigma):
    """J = int H(u) K(u) du over u in (-inf, inf), K = exp(-u^2/2 sigma^2), using H(-u) = conj H(u) and the grid u >= 0
    (uniform step, u[0] = 0); returns the real J (trapezoid)."""
    u = np.asarray(u, float); hstep = u[1] - u[0]
    K = np.exp(-u ** 2 / (2 * sigma ** 2)); wts = np.full(len(u), hstep); wts[0] = hstep / 2
    return 2 * np.real(H @ (wts * K))


def stationary_R(N, w_taper, n, Lam, t_lo, t_hi):
    """prime side of G = sum_gamma e^{2 i theta(gamma)} w_taper(gamma): -(1/2pi) sum_n Lam(n) n^{-1/2} hhat(ln n), with
    hhat(ln n) = int e^{i(2 theta(t) - t ln n)} w(t) dt by stationary phase at 2 theta'(t_n) = ln n (exact theta)."""
    n = np.asarray(n, float); ln_n = np.log(n)
    t = TWO_PI * np.sqrt(n / N)                              # first guess: 2 theta' ~ ln(N t^2/4 pi^2)
    from scipy.special import digamma
    def trigamma(z):                                     # psi'(z), z complex: shift by 10, then the asymptotic series
        acc = 0
        for j in range(10): acc = acc + 1 / (z + j) ** 2
        y = z + 10
        return acc + 1 / y + 1 / (2 * y ** 2) + 1 / (6 * y ** 3) - 1 / (30 * y ** 5) + 1 / (42 * y ** 7) - 1 / (30 * y ** 9)
    for _ in range(4):
        z = 1 + 1j * t
        f = math.log(N) - 2 * math.log(TWO_PI) + 2 * np.real(digamma(z)) - ln_n        # 2 theta'(t) - ln n
        fp = -2 * np.imag(trigamma(z))                                              # d/dt 2 Re psi(1+it) = -2 Im psi'(1+it)
        t = t - f / fp
    z = 1 + 1j * t
    phi2 = -2 * np.imag(trigamma(z))                                                 # 2 theta''(t) > 0
    phi = 2 * theta_E(N, t) - t * ln_n
    m = (t > t_lo) & (t < t_hi)
    amp = w_taper(t[m]) * np.sqrt(TWO_PI / phi2[m]) * np.exp(1j * (phi[m] + math.pi / 4))
    return -(1 / TWO_PI) * np.sum(Lam[m] * n[m] ** -0.5 * amp)
