"""Exact (recipe / ratios conjecture, Conrey-Snaith 2007 sect. 7.1) lower-order terms of sum_gamma |L'(rho, chi)|^2 for a
real primitive chi mod q: zeta(1+s) -> zeta_q(1+s) = zeta(1+s) prod_{p|q}(1 - p^{-1-s}), t/2pi -> q t/2pi.
Density per unit t: m1(t) = (1/2pi) [2 c_{beta gamma delta} S + L f''(0)], S = the three recipe terms, f(s) = Z(s) + e^{-sL} Z(-s).
Checks the q = 1 polynomial against Conrey-Snaith's theorem, then compares the bad-prime ratio with the LDM model."""
import sympy as sp, mpmath as mp, numpy as np
mp.mp.dps = 30
L, e, x, y, w = sp.symbols('L e x y w')
c = sp.symbols('c0:6')
def Z(s): return 1 / s + sum(c[j] * s ** j for j in range(6))
b_, g_, d_ = e * x, e * y, e * w
S = (Z(b_ + d_) * Z(g_ + d_) / Z(d_) + sp.exp(-(b_ + d_) * L) * Z(-d_ - b_) * Z(g_ - b_) / Z(-b_)
     + sp.exp(-(g_ + d_) * L) * Z(b_ - g_) * Z(-d_ - g_) / Z(-g_))
ser = sp.series(S, e, 0, 4).removeO()
c3 = sp.expand(ser.coeff(e, 3)); cxyw = sp.expand(c3).coeff(x, 1).coeff(y, 1).coeff(w, 1)
s = sp.symbols('s'); f = Z(s) + sp.exp(-s * L) * Z(-s)
f2 = sp.limit(sp.diff(sp.series(f, s, 0, 4).removeO(), s, 2), s, 0)
poly = sp.expand(2 * cxyw + L * f2)   # times (1/2pi) c_q
print('m1 * 2pi / c_q =', sp.collect(poly, L))
# q = 1 check: zeta(1+s) = 1/s + sum (-1)^n gamma_n s^n / n!
gam = [mp.stieltjes(n) for n in range(5)]
zc = {c[n]: sp.Float(str((-1) ** n * gam[n] / mp.factorial(n)), 30) for n in range(5)}; zc[c[5]] = 0
pz = sp.Poly(sp.expand(poly.subs(zc)), L)
g0, g1, g2, g3 = [sp.Float(str(gam[n]), 30) for n in range(4)]
cs = {4: sp.Rational(1, 24) * 2, 3: g0 / 3 * 2, 2: (g0 ** 2 / 2 - g1) * 2, 1: -(g0 ** 3 + 5 * g0 * g1 + g2 / 2) * 2,
      0: (g0 ** 4 + 6 * g0 ** 2 * g1 + 7 * g1 ** 2 + 4 * g0 * g2 + sp.Rational(5, 3) * g3) * 2}
print('q=1 check against Conrey-Snaith (coefficients of L^j, ours vs CS*2pi/pi... ):')
for j in range(4, -1, -1): print('  L^%d  %.10f  %.10f' % (j, float(pz.coeff_monomial(L ** j)), float(cs[j]) / 2 * 2))
# zeta_q Laurent coefficients
def coeffs(primes):
    def zq(t): return mp.zeta(1 + t) * mp.fprod([(1 - mp.mpf(p) ** (-1 - t)) for p in primes]) / mp.fprod([(1 - mp.mpf(1) / p) for p in primes])
    # c_j of Z(s) = zq(s) (normalised so the pole is 1/s): Taylor of s*zq(s)
    tay = mp.taylor(lambda t: t * zq(t) if t != 0 else mp.mpf(1), 0, 6)
    return {c[j]: sp.Float(str(tay[j + 1]), 25) for j in range(6)}
P1 = sp.lambdify(L, poly.subs(coeffs([])), 'mpmath')
print('\nbad-prime ratio at k = 1: recipe vs LDM model (first order and resummed)')
z0 = np.exp(1j * 2 * np.pi * (np.arange(1 << 14) + 0.5) / (1 << 14))
for p in (2, 3, 5):
    Pp = sp.lambdify(L, poly.subs(coeffs([p])), 'mpmath')
    wv = p ** -0.5 * z0; R = np.real(wv / (1 - wv)); W = np.abs(1 - wv) ** -2.0
    for Lv in (7.5, 9.0, 10.0, 11.0, 12.0, 20.0, 50.0):
        rec = (1 - 1 / p) * float(Pp(Lv) / P1(Lv)); eps = -2 * np.log(p) / Lv * R
        first = (1 - 1 / p) * (1 + 8 * np.log(p) / ((p - 1) * Lv)); res = (1 - 1 / p) * np.mean(W) / np.mean(W * (1 + eps) ** 4)
        lin = (1 - 1 / p) / (1 + 4 * np.mean(W * eps) / np.mean(W))
        print('  p=%d L=%5.1f  recipe %.5f   LDM first-order %.5f   LDM 1/(1+4<eps>) %.5f   LDM resummed (1+eps)^4 %.5f' % (p, Lv, rec, first, lin, res))
