"""Recipe (ratios conjecture; Conrey-Snaith 2007 sect. 7.1) for the density in t of sum_gamma |L'(rho,chi)|^2, chi real
primitive mod q: zeta -> zeta_q(s) = zeta(s) prod_{p|q} (1 - p^-s), t/2pi -> q t/2pi, L = ln(q t/2pi).
  m1(L) = (1/2pi) [ 2 d_beta d_gamma d_delta S + L f''(0) ],
  S = zq(1+b+d) zq(1+g+d)/zq(1+d) + e^{-(b+d)L} zq(1-d-b) zq(1+g-b)/zq(1-b) + e^{-(g+d)L} zq(1+b-g) zq(1-d-g)/zq(1-g),
  f(s) = zq(1+s) + e^{-sL} zq(1-s).
Mixed derivatives by Cauchy's formula on the torus |b| = r, |g| = 2r, |d| = 3r (no denominator vanishes there; S is analytic
inside). Checked at q = 1 against Conrey-Snaith's polynomial; then the bad-prime ratio m1(q=p)/m1(q=1) at equal L is
compared with the local-density-modulation (LDM) model."""
import mpmath as mp, numpy as np, itertools
mp.mp.dps = 40
def zq(s, P): return mp.zeta(s) * mp.fprod([1 - mp.mpf(p) ** (-s) for p in P])
def S(b, g, d, L, P):
    return (zq(1 + b + d, P) * zq(1 + g + d, P) / zq(1 + d, P) + mp.exp(-(b + d) * L) * zq(1 - d - b, P) * zq(1 + g - b, P) / zq(1 - b, P)
            + mp.exp(-(g + d) * L) * zq(1 + b - g, P) * zq(1 - d - g, P) / zq(1 - g, P))
def m1(L, P, n=8, r=mp.mpf('0.02')):
    L = mp.mpf(L); pts = [mp.exp(2j * mp.pi * (k + mp.mpf(1) / 3) / n) for k in range(n)]
    acc = 0
    for u, v, w in itertools.product(pts, pts, pts):
        b, g, d = r * u, 2 * r * v, 3 * r * w
        acc += S(b, g, d, L, P) / (b * g * d)
    c3 = acc / n ** 3                      # coefficient of b g d
    f = lambda s: zq(1 + s, P) + mp.exp(-s * L) * zq(1 - s, P)
    c2 = sum(f(r * u) / (r * u) ** 2 for u in pts) / n   # coefficient of s^2
    return (2 * mp.re(c3) + L * 2 * mp.re(c2)) / (2 * mp.pi)
g = [mp.stieltjes(k) for k in range(4)]
def cs(L):
    return (L ** 4 / 24 + g[0] / 3 * L ** 3 + (g[0] ** 2 / 2 - g[1]) * L ** 2 - (g[0] ** 3 + 5 * g[0] * g[1] + g[2] / 2) * L
            + g[0] ** 4 + 6 * g[0] ** 2 * g[1] + 7 * g[1] ** 2 + 4 * g[0] * g[2] + mp.mpf(5) / 3 * g[3]) / mp.pi
for L in (5, 10):
    print('q=1 check L=%d: Cauchy %.12f  Conrey-Snaith %.12f' % (L, float(m1(L, [])), float(cs(L))), flush=True)
z0 = np.exp(1j * 2 * np.pi * (np.arange(1 << 14) + 0.5) / (1 << 14))
print('bad-prime ratio of the k = 1 moment at equal L: recipe | LDM first order | LDM resummed (1+eps)^4')
for p in (2, 3, 5, 7):
    wv = p ** -0.5 * z0; R = np.real(wv / (1 - wv)); W = np.abs(1 - wv) ** -2.0
    for L in (7.5, 9, 10, 11, 12, 20, 50):
        rec = float(m1(L, [p]) / m1(L, [])); eps = -2 * np.log(p) / L * R
        first = (1 - 1 / p) * (1 + 8 * np.log(p) / ((p - 1) * L)); res = (1 - 1 / p) * np.mean(W) / np.mean(W * (1 + eps) ** 4)
        print('  p=%d L=%5.1f  %.5f | %.5f | %.5f   (recipe-lead)*L/(ln p (1-1/p)) %.3f' % (p, L, rec, first, res, (rec - (1 - 1 / p)) * L / (np.log(p) * (1 - 1 / p))), flush=True)
