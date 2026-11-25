"""H2: the local factor of the constant in degree 4, from scratch, on random Satake parameters.
(a) USp(4) Haar (Weyl density (cos t1 - cos t2)^2 sin^2 t1 sin^2 t2), Sym^3 of SU(2) (Sato-Tate), and E1 x E2 (two independent
SU(2)): F_p(X) = sum_k conj(h_k(alpha)) (-1)^k e_k(alpha) X^k; F(1) = 0; G = F/(1 - X) of degree 3 with
g1 = 1 - |lambda|^2, top -conj h_4 (central character 1); explicit g2 = 1 - e1^2 + e1^2 e2 - e2^2 for USp(4).
(b) |lambda|^2 - 1 against characters: USp(4): Sym^2 (10) + Lambda^2_0 (5); Sym^3: Sym^2 + Sym^4 + Sym^6; RS: S2(g1) + S2(g2) + S2 S2.
Haar means, norms <chi, chi>, cross terms, and Frobenius-Schur indicators E[chi(g^2)] (the prime-square bias nu).
(c) B_p(X) = F_p(X) / prod_{i,j}(1 - alpha_i conj(alpha_j) X) = 1 + O(X^2) (so prod G_p(1/p) = B(1)/L(1, Ad) with Ad the
15-dim adjoint L x Ltilde / zeta), and prod_{i,j}(1 - a_i conj a_j X) = (1 - X) det(1 - X Sym^2) det(1 - X Lambda^2_0) for USp(4).
(d) E over the group of log G(1/p) to second order, i.e. the size of the O(1/p^2) terms the first-order tail drops."""
import numpy as np, math
rng = np.random.default_rng(7)
def sample_usp4(n):
    out = []
    while len(out) < n:
        t = rng.uniform(0, math.pi, (4 * n, 2)); w = (np.cos(t[:, 0]) - np.cos(t[:, 1])) ** 2 * np.sin(t[:, 0]) ** 2 * np.sin(t[:, 1]) ** 2
        acc = rng.uniform(0, 4, len(w)) < w; out.extend(t[acc].tolist())
    t = np.array(out[:n]); return np.exp(1j * np.stack([t[:, 0], -t[:, 0], t[:, 1], -t[:, 1]], 1)), t
def sample_su2(n):
    out = []
    while len(out) < n:
        t = rng.uniform(0, math.pi, 4 * n); acc = rng.uniform(0, 1, len(t)) < np.sin(t) ** 2; out.extend(t[acc].tolist())
    return np.array(out[:n])
def ek(a):
    return np.array([np.poly(x)[::1] * np.array([(-1) ** k for k in range(len(x) + 1)]) for x in a])   # e_0..e_4
def hk(e, K):
    """complete homogeneous h_0..h_K from e (Newton: sum_j (-1)^j e_j h_{k-j} = 0)."""
    n = e.shape[0]; h = np.zeros((n, K + 1), complex); h[:, 0] = 1
    for k in range(1, K + 1):
        s = 0
        for j in range(1, min(k, e.shape[1] - 1) + 1): s = s + (-1) ** (j + 1) * e[:, j] * h[:, k - j]
        h[:, k] = s
    return h
def check(name, a):
    e = ek(a); h = hk(e, 6); hb = hk(ek(a.conj()), 6)
    F = np.stack([hb[:, k] * (-1) ** k * e[:, k] for k in range(5)], 1)              # F_p coefficients X^0..X^4
    F1 = np.abs(F.sum(1)).max()
    G = np.cumsum(F, 1)[:, :4]                                                       # G = F/(1 - X): partial sums; remainder F(1)
    lam2 = np.abs(e[:, 1]) ** 2
    print('%-8s max|F(1)| %.1e | max|g1 - (1 - |lam|^2)| %.1e | max|g3 + conj h4| %.1e' % (name, F1, np.abs(G[:, 1] - (1 - lam2)).max(), np.abs(G[:, 3] + hb[:, 4]).max()))
    # (c) B_p = F / prod(1 - a_i conj a_j X): coefficient of X
    prodc = np.ones((len(a), 1), complex)
    for i in range(4):
        for j in range(4):
            prodc = np.array([np.convolve(pc, [1, -a[m, i] * a[m, j].conj()]) for m, pc in enumerate(prodc)])
    inv = np.zeros((len(a), 6), complex); inv[:, 0] = 1                              # 1/prodc to X^5
    for k in range(1, 6):
        inv[:, k] = -sum(prodc[:, j] * inv[:, k - j] for j in range(1, min(k, prodc.shape[1] - 1) + 1))
    B = np.array([np.convolve(F[m], inv[m])[:4] for m in range(len(a))])
    print('         B_p = F_p/prod(1 - a_i conj a_j X): max|coef X^1| %.1e, mean coef X^2 %+.3f, max|coef X^2| %.2f' % (np.abs(B[:, 1]).max(), B[:, 2].real.mean(), np.abs(B[:, 2]).max()))
    return e, h, G, lam2
n = 20000
a, t = sample_usp4(n); e, h, G, lam2 = check('USp(4)', a)
e1 = e[:, 1].real; e2 = e[:, 2].real
print('         max|g2 - (1 - e1^2 + e1^2 e2 - e2^2)| %.1e' % np.abs(G[:, 2] - (1 - e1 ** 2 + e1 ** 2 * e2 - e2 ** 2)).max())
chi10 = h[:, 2].real; chi5 = e2 - 1
print('         |lam|^2 - 1 - (chi10 + chi5): max %.1e' % np.abs(lam2 - 1 - chi10 - chi5).max())
a2 = a ** 2; e_sq = ek(a2); h_sq = hk(e_sq, 2)
fs10 = h_sq[:, 2].real.mean(); fs5 = (e_sq[:, 2].real - 1).mean()
print('         Haar: E chi10 %+.3f E chi5 %+.3f | E chi10^2 %.3f E chi5^2 %.3f E chi10 chi5 %+.3f | FS: E chi10(g^2) %+.3f E chi5(g^2) %+.3f -> nu = %.2f' % (
      chi10.mean(), chi5.mean(), (chi10 ** 2).mean(), (chi5 ** 2).mean(), (chi10 * chi5).mean(), fs10, fs5, fs10 + fs5))
print('         moments of e1: M2 %.3f M4 %.3f M6 %.3f M8 %.2f ; E[g2] %+.3f, E[e2 h2] %+.3f' % tuple([(e1 ** k).mean() for k in (2, 4, 6, 8)] + [G[:, 2].real.mean(), (e2 * chi10).mean()]))
# factorisation of prod(1 - a_i conj a_j X) for USp(4)
th1, th2 = t[:, 0], t[:, 1]; X = 0.37
lhs = np.prod([[1 - a[:, i] * a[:, j].conj() * X for j in range(4)] for i in range(4)], axis=(0, 1))
s2 = (1 - X) ** 2 * np.prod([1 - 2 * np.cos(u) * X + X * X for u in (2 * th1, 2 * th2, th1 + th2, th1 - th2)], 0)
l5 = (1 - X) * np.prod([1 - 2 * np.cos(u) * X + X * X for u in (th1 + th2, th1 - th2)], 0)
print('         prod(1 - a_i conj a_j X) vs (1 - X) det(1 - X Sym2) det(1 - X L2_0) at X = 0.37: max rel diff %.1e' % np.abs(lhs / ((1 - X) * s2 * l5) - 1).max())
# Sym^3 of SU(2)
ph = sample_su2(n); a3 = np.exp(1j * np.outer(ph, [3, 1, -1, -3])); e, h, G, lam2 = check('Sym^3', a3)
S = lambda k, x: np.sin((k + 1) * x) / np.sin(x)
print('         |lam|^2 - 1 - (S2 + S4 + S6): max %.1e ; FS nu = E[S2+S4+S6](g^2) = %.2f' % (np.abs(lam2 - 1 - S(2, ph) - S(4, ph) - S(6, ph)).max(), (S(2, 2 * ph) + S(4, 2 * ph) + S(6, 2 * ph)).mean()))
# RS
p1, p2 = sample_su2(n), sample_su2(n); ar = np.stack([np.exp(1j * (p1 + p2)), np.exp(1j * (p1 - p2)), np.exp(1j * (-p1 + p2)), np.exp(-1j * (p1 + p2))], 1)
e, h, G, lam2 = check('E1xE2', ar)
print('         |lam|^2 - 1 - (S2(1) + S2(2) + S2 S2): max %.1e ; FS nu = %.2f' % (np.abs(lam2 - 1 - S(2, p1) - S(2, p2) - S(2, p1) * S(2, p2)).max(), (S(2, 2 * p1) + S(2, 2 * p2) + S(2, 2 * p1) * S(2, 2 * p2)).mean()))
# the trivial representation inside L x Ltilde: E|lam|^2 (pole order) and E lam^2 (sign of the pairing)
for nm, aa in (('USp(4)', a), ('Sym^3', a3), ('E1xE2', ar)):
    l = aa.sum(1); print('  %-7s E|lam|^2 = %.3f (pole of L x Ltilde), E[lam(g)^2 - lam(g^2)]/2 = E tr Lambda^2 = %.3f, E tr g^2 = %+.3f' % (nm, (np.abs(l) ** 2).mean(), ((l ** 2 - (aa ** 2).sum(1)) / 2).real.mean(), (aa ** 2).sum(1).real.mean()))
