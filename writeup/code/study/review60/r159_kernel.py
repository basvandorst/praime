"""Own check of chi(-2w) Ktilde21(w) = G_3(t = i w), mu = (1, k-1, k); Ktilde21(w) = int_0^oo K21(z) z^{w-1} dz,
K21(z) = 4 pi (-1)^{k/2} sqrt z cos(4 pi z) J_{k-1}(4 pi z). Ktilde21 = (-1)^{k/2} (4 pi)^{1/2-w} I(w), I(w) = int x^{w-1/2} J_{k-1}(x) cos x dx,
I by DLMF 10.22.57 with cos x = sqrt(pi x/2) J_{-1/2}(x); I checked once by direct quadrature."""
import mpmath as mp
mp.mp.dps = 30
def I_closed(w, k):
    m, n, lam = k - 1, mp.mpf(-1) / 2, -w
    return mp.sqrt(mp.pi / 2) * mp.gamma(lam) * mp.gamma((n + m - lam + 1) / 2) / (
        2 ** lam * mp.gamma((-n + m + lam + 1) / 2) * mp.gamma((n + m + lam + 1) / 2) * mp.gamma((n - m + lam + 1) / 2))
def I_quad(w, k):
    f = lambda x: x ** (w - mp.mpf(1) / 2) * mp.besselj(k - 1, x) * mp.cos(x)
    return mp.quadosc(f, [0, mp.inf], period=mp.pi)
GR = lambda s: mp.pi ** (-s / 2) * mp.gamma(s / 2)
chi = lambda s: mp.zeta(s) / mp.zeta(1 - s)
def G3(w, k):
    return mp.fprod(GR(mp.mpf(1) / 2 + w + mu) / GR(mp.mpf(1) / 2 - w + mu) for mu in (1, k - 1, k))
def Kt(w, k):
    return (-1) ** (k // 2) * (4 * mp.pi) ** (mp.mpf(1) / 2 - w) * I_closed(w, k)
mp.mp.dps = 20
for w, k in ((mp.mpf(-0.6), 2), (mp.mpf(-1.3), 4)):
    a, b = I_closed(w, k), I_quad(w, k)
    print("quadrature check I(%s), k=%d: closed %s quad %s rel %.1e" % (w, k, mp.nstr(a, 12), mp.nstr(b, 12), abs(a / b - 1)))
mp.mp.dps = 30
for w in (mp.mpc(-0.25, 3.0), mp.mpc(-0.4, 10.5), mp.mpc(0.3, -25.0)):
    for k in (2, 4, 6):
        r = chi(-2 * w) * Kt(w, k) / G3(w, k)
        print("w=%s k=%d: chi(-2w)Kt/G3 - 1 = %.1e" % (mp.nstr(w, 4), k, abs(r - 1)))
import numpy as np
P = [p for p in range(2, 2000000) if all(p % q for q in range(2, int(p ** .5) + 1))] if False else None
