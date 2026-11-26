"""Sanity checks of core.py: Bessel closed form vs direct convolution of Phi, evenness, Xi moments."""
import mpmath as mp
from core import Basis, Phi, Xi, xi_moments
import time

mp.mp.dps = 40
# 1. int Phi = Xi(0) = xi(1/2)
h = mp.mpf(1) / 100
s = Phi(0) + 2 * sum(Phi(i * h) for i in range(1, 300))
print("int Phi (trapezoid)        =", mp.nstr(s * h, 25))
print("Xi(0)                      =", mp.nstr(Xi(0), 25))

# 2. Phi*Phi at u=0, 0.7 directly by trapezoid vs Bessel closed form
bas = Basis(3, 1)
for u in [mp.mpf(0), mp.log(2), mp.mpf('1.3')]:
    f = lambda w: Phi(w) * Phi(u - w)
    # centre at u/2
    c = u / 2
    tr = f(c) + sum(f(c + i * h) + f(c - i * h) for i in range(1, 300))
    v = bas.eval(u)
    print("u=%s direct=%s bessel=%s" % (mp.nstr(u, 6), mp.nstr(tr * h, 25), mp.nstr(v[0][0], 25)))

# 3. evenness of the closed form (needs many N at negative u)
for u in ['0.2', '0.05']:
    a = bas.eval(mp.mpf(u))
    b = bas.eval(-mp.mpf(u))
    print("even check u=%s: %s %s | g1: %s %s" % (u, mp.nstr(a[0][0], 20), mp.nstr(b[0][0], 20),
                                                   mp.nstr(a[1][0], 20), mp.nstr(b[1][0], 20)))

# 4. moments m_j vs g_j(0)
t0 = time.time()
mp.mp.dps = 30
m, a = xi_moments(3, T=120)
print("xi moments time", time.time() - t0)
mp.mp.dps = 40
v = bas.eval(0)
for j in range(4):
    print("j=%d  g_j(0)=%s  m_j=%s" % (j, mp.nstr(v[j][0], 20), mp.nstr(m[j], 20)))
