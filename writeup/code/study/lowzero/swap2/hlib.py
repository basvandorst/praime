"""swap2 prediction library (own code; no repo imports).

Swapped term of the cross two-point function of the zeros of two primitive Dirichlet L-functions (ratios recipe):
  R12_swap(t1, t2) = (1/4 pi^2) { k1 k2 eps1 conj(eps2) e^{2i(theta2(t2) - theta1(t1))} H12(t2 - t1) + c.c. },
  H12(u) = L(1 - iu, psi) L(1 + iu, conj psi) A12(-iu),  psi = chi1 conj(chi2) (mod lcm, psi(p) = 0 for p | q1 q2),
  per prime p not dividing q1 q2:  h_p(u) = (1 - 2/p + y) / ((1 - conj y)(1 - 1/p)^2),  y = psi(p) p^{iu} / p,
  h_p = 1 for p | q1 q2;  k_j = q_j/phi(q_j).
Statistic  X12 = (1/T) sum_{g1, g2} e^{2i theta1(g1) - 2i theta2(g2)} K(g2 - g1),  K(u) = exp(-u^2 / 2 sigma^2):
  E X12 = (1/4 pi^2) k1 k2 eps1 conj(eps2) J12(sigma),   J12 = int H12(u) K(u) du.
One-point law: sum_{T1 < g <= T2} e^{2i theta_chi(g)} = -eps_chi k_chi (T2 - T1)/2pi + o(T).
"""
import math
import numpy as np

UMAX, DU = 6.0, 0.05
UG = np.arange(0.0, UMAX + DU / 2, DU)          # u >= 0 grid; K(6 sigma) = e^-18


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]


def phi_ratio(k):
    r = 1.0; n = k; p = 2
    while p * p <= n:
        if n % p == 0:
            r *= p / (p - 1)
            while n % p == 0: n //= p
        p += 1
    if n > 1: r *= n / (n - 1)
    return r


def quad_weights(sig, ug=UG, du=DU):
    """trapezoid weights on the symmetric grid folded to u >= 0: J = sum_u w_u [H(u) + H(-u)] (u > 0) + w_0 H(0)."""
    w = np.exp(-ug ** 2 / (2 * sig * sig)) * du
    w[0] *= 0.5        # folded: the u = 0 node is counted once in H(u) + H(-u) at u = 0, i.e. 2 * w0/2
    return w


# ---------------------------------------------------------------- character tables (PARI)
def pari_instance():
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
    return pari


def real_table(pari, Ds, P):
    """chi_D(p) for p <= P as int8 (Kronecker symbol), rows = Ds."""
    pr = primes_upto(P)
    X = np.zeros((len(Ds), len(pr)), np.int8)
    for i, D in enumerate(Ds):
        v = pari('apply(p->kronecker(%d,p), primes([2,%d]))' % (D, P))
        X[i] = np.array([int(x) for x in v], np.int8)
    return pr, X


def conrey_table(pari, N, m, P):
    """chi(p) (complex, 0 if p | N) for the Conrey character chi_N(m, .), p <= P."""
    pr = primes_upto(P)
    v = pari('my(G=znstar(%d,1),c=znconreylog(G,%d)); apply(p->if(%d%%p,chareval(G,c,p),-1), primes([2,%d]))' % (N, m, N, P))
    out = np.zeros(len(pr), complex)
    for i, x in enumerate(v):
        fr = float(x)
        if fr < 0:
            continue
        out[i] = np.exp(2j * math.pi * fr)
    return pr, out


def conrey_info(pari, N, m):
    """(conductor, parity a, root number eps = tau(chi)/(i^a sqrt q) from PARI lfunrootres; also tau computed directly)."""
    pari('G=znstar(%d,1); c=znconreylog(G,%d)' % (N, m))
    q = int(pari('znconreyconductor(G,%d)' % m)); a = int(pari('zncharisodd(G,c)'))
    eps = complex(pari('lfunrootres(lfuncreate(Mod(%d,%d)))[3]' % (m, N)))
    tau = complex(pari('sum(n=1,%d,if(gcd(n,%d)==1,exp(2*Pi*I*chareval(G,c,n))*exp(2*Pi*I*n/%d),0))' % (N, N, N)))
    eps_tau = tau / (1j ** a * math.sqrt(q))
    return q, a, eps, eps_tau


# ---------------------------------------------------------------- H on the u-grid
def h_factor(psi, p, x):
    """Euler factor h_p for psi values (array, broadcast) and x = p^{iu}."""
    y = psi * x / p
    return (1 - 2.0 / p + y) / ((1 - np.conj(y)) * (1 - 1.0 / p) ** 2)


def J_real_all(X, pr, sigmas, ug=UG, du=DU, verbose=False):
    """All pairs of real characters at once (X: n x np int8 of Kronecker symbols at primes pr).
    log H_ij(u) = sum_{p>=5} [chi_i chi_j (a_p - b_p)/2 + chi_i^2 chi_j^2 (a_p + b_p)/2],  a_p = ln h_p(+1), b_p = ln h_p(-1);
    p = 2, 3 multiplied directly (h_3 vanishes at u = 0 for psi(3) = -1).  Returns {sigma: n x n real J} and H(0) matrix."""
    n = X.shape[0]
    big = pr >= 5
    Xb = X[:, big].astype(float); Pb = pr[big].astype(float); X2b = Xb * Xb
    lp = np.log(Pb)
    x2 = X[:, 0].astype(float); x3 = X[:, 1].astype(float)       # chi(2), chi(3)
    assert pr[0] == 2 and pr[1] == 3
    W = {s: quad_weights(s, ug, du) for s in sigmas}
    J = {s: np.zeros((n, n)) for s in sigmas}
    H0 = None
    A = np.hstack([Xb, X2b])
    for k, u in enumerate(ug):
        xp = np.exp(1j * u * lp)
        a = np.log(h_factor(1.0, Pb, xp)); b = np.log(h_factor(-1.0, Pb, xp))
        c = np.concatenate([(a - b) / 2, (a + b) / 2])
        L = (A * c.real) @ A.T + 1j * ((A * c.imag) @ A.T)
        H = np.exp(L)
        for p, xx in ((2, x2), (3, x3)):
            psi = np.outer(xx, xx)                       # real psi(p) in {-1, 0, 1}
            hp = np.where(psi != 0, h_factor(psi, p, np.exp(1j * u * math.log(p))), 1.0)
            H = H * hp
        if k == 0: H0 = H.real.copy()
        for s in sigmas:
            J[s] += 2 * W[s][k] * H.real                 # H(u) + H(-u) = 2 Re H(u) for real psi
        if verbose and k % 20 == 0: print('  u = %.2f' % u, flush=True)
    return J, H0


def J_general_all(C, pr, sigmas, ug=UG, du=DU, P0=100, K=4, verbose=False):
    """All ordered pairs (i, j) of characters given by complex tables C (n x np; 0 where p | q), psi = chi_i conj chi_j.
    p <= P0 multiplied directly; p > P0 via the series of ln h_p in y to order K.
    Returns {sigma: n x n complex J_ij} (J_ji = conj J_ij is a check)."""
    n = C.shape[0]
    small = pr <= P0; Ps = pr[small].astype(float); Cs = C[:, small]
    Pb = pr[~small].astype(float); Cb = C[:, ~small]; lpb = np.log(Pb)
    nz = (np.abs(Cb) > 0).astype(float)
    c0 = (np.log(1 - 2 / Pb) - 2 * np.log(1 - 1 / Pb))
    L0 = (nz * c0) @ nz.T                                    # u-independent part
    Mk = [Cb ** k for k in range(1, K + 1)]
    W = {s: quad_weights(s, ug, du) for s in sigmas}
    J = {s: np.zeros((n, n), complex) for s in sigmas}
    PSI = Cs[:, None, :] * np.conj(Cs[None, :, :])          # n x n x ns
    for kk, u in enumerate(ug):
        Hs = []
        for uu in ((u,) if kk == 0 else (u, -u)):
            x = np.exp(1j * uu * lpb)
            L = L0.astype(complex)
            for k in range(1, K + 1):
                al = ((-1) ** (k + 1) / k) * (x / Pb) ** k / (1 - 2 / Pb) ** k
                be = (1.0 / k) * (np.conj(x) / Pb) ** k
                M = Mk[k - 1]
                L = L + (M * al) @ np.conj(M).T + (np.conj(M) * be) @ M.T
            xs = np.exp(1j * uu * np.log(Ps))
            hs = np.where(PSI != 0, h_factor(PSI, Ps, xs), 1.0)
            Hs.append(np.exp(L) * np.prod(hs, axis=2))
        Hsum = Hs[0] if kk == 0 else Hs[0] + Hs[1]
        for s in sigmas:
            J[s] += W[s][kk] * (Hsum if kk else 2 * Hsum)   # fold: node 0 weight already halved
        if verbose and kk % 20 == 0: print('  u = %.2f' % u, flush=True)
    return J


def J_pair_direct(psi, pr, sigmas, ug=UG, du=DU):
    """one pair, plain Euler product over all given primes (psi array over pr), both signs of u."""
    W = {s: quad_weights(s, ug, du) for s in sigmas}
    keep = psi != 0; P = pr[keep].astype(float); s_ = psi[keep]; lp = np.log(P)
    out = {s: 0j for s in sigmas}
    for k, u in enumerate(ug):
        tot = 0j
        for uu in ((u,) if k == 0 else (u, -u)):
            tot += np.prod(h_factor(s_, P, np.exp(1j * uu * lp)))
        for s in sigmas:
            out[s] += W[s][k] * (tot if k else 2 * tot)
    return out


def predict(kap_i, kap_j, eps_i, eps_j, J):
    return kap_i * kap_j * eps_i * np.conj(eps_j) * J / (4 * math.pi ** 2)
