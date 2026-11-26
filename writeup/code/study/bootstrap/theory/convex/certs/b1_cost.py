"""Deliverable 1: cost/g(0) for B = 1 (g = Phi*Phi), by the prime sum and by the archimedean formula."""
import mpmath as mp
from core import Basis, xi_moments, von_mangoldt

mp.mp.dps = 40
bas = Basis(0, 0)
g0 = bas.eval(0)[0][0]
prime = mp.mpf(0)
terms = []
for n in range(2, 80):
    L = von_mangoldt(n)
    if L == 0:
        continue
    t = 2 * L / mp.sqrt(n) * bas.eval(mp.log(n))[0][0]
    prime += t
    terms.append((n, t))
print("g(0) = (Phi*Phi)(0)            =", mp.nstr(g0, 30))
print("prime-side cost 2 sum L(n) n^-1/2 g(log n), n<80 =", mp.nstr(prime, 30))
for n, t in terms[:8]:
    print("   n=%2d  2 Lambda(n) n^-1/2 g(log n) = %s   (share %.3e)" % (n, mp.nstr(t, 15), float(t / prime)))

# archimedean side: 2 int g cosh(u/2) + (1/2pi) int ghat (Re psi(1/4+it/2) - log pi)
mp.mp.dps = 30
m, a = xi_moments(0, T=130, panel=1.0, nodes=24)
mp.mp.dps = 40
# direct trapezoid of 2 int g cosh(u/2) (should be B(i/2)/2 = 1/2)
h = mp.mpf(1) / 200
s = bas.eval(0)[0][0] + 2 * sum(bas.eval(i * h)[0][0] * mp.cosh(i * h / 2) for i in range(1, 900))
print("2 int g cosh(u/2) (trapezoid)  =", mp.nstr(2 * s * h, 25), " (exact 1/2 = 2 Xi(i/2)^2 B(i/2))")
arch = mp.mpf(1) / 2 + a[0]
print("arch-side cost 1/2 + (1/2pi) int Xi^2 (Re psi - log pi) =", mp.nstr(arch, 25))
print("difference prime - arch        =", mp.nstr(prime - arch, 5))
print("cost/g(0) (prime side)         =", mp.nstr(prime / g0, 20))
print("cost/g(0) (arch side)          =", mp.nstr(arch / g0, 20))
