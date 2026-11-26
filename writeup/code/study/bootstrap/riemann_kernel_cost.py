# Cost of the certificate with transform Xi(t)^2 (g = Phi*Phi, Phi Riemann's kernel), by the cost identity
#   W(g) = ghat(i/2) + ghat(-i/2) + (1/2pi) int ghat(t) Omega(t) dt,  Omega = Re psi(1/4 + it/2) - log pi,
# with ghat(+-i/2) = xi(0)^2 = 1/4.  Also the constant of the truncation lemma.
from mpmath import mp, mpf, quad, zeta, gamma, pi, digamma, log, re, inf, catalan
mp.dps = 25
xi = lambda s: s*(s-1)/2*pi**(-s/2)*gamma(s/2)*zeta(s)
Xi2 = lambda t: re(xi(mpf(1)/2 + 1j*t))**2
Om = lambda t: re(digamma(mpf(1)/4 + 1j*t/2)) - log(pi)
pts = [0, 5, 10, 14.134725, 21.02204, 25.010858, 30.424876, 32.935062, 37.586178, 40.918719, 43.327073, 48.005151, 49.773832, 60, 80, 120]
g0 = 2*quad(Xi2, pts)/(2*pi); arch = 2*quad(lambda t: Xi2(t)*Om(t), pts)/(2*pi)
print('g(0) = %s   pole part = 1/2   archimedean part = %s' % (mp.nstr(g0, 10), mp.nstr(arch, 12)))
print('W(g)/g(0) = %s' % mp.nstr((mpf(1)/2 + arch)/g0, 10))
c1 = 7*zeta(3) + pi**3/4
print('c1 = 7 zeta(3) + pi^3/4 = %s   quadrature: %s   12 c1 = %s' % (mp.nstr(c1, 10), mp.nstr(quad(lambda u: u**2*mp.exp(-u/2)/(1-mp.exp(-2*u)), [0, 1, 10, inf]), 10), mp.nstr(12*c1, 6)))
