"""Is chi(2it) * Ktilde21(-it) a constant multiple of the Sym^2 gamma ratio G3(t)?  (Poisson in q of the Petersson lines.)
Ktilde21(w) = int_0^inf z^{w-1} K21(z) dz in closed form (Weber-Schafheitlin):  K21 = 4 pi (-1)^{k/2} sqrt z cos(4 pi z) J_nu(4 pi z)."""
import mpmath as mp
mp.mp.dps = 30

def MJJ(s, mu, nu):   # int_0^inf x^{s-1} J_mu(x) J_nu(x) dx
    return 2 ** (s - 1) * mp.gamma(1 - s) * mp.gamma((mu + nu + s) / 2) / (
        mp.gamma(1 + (nu - mu - s) / 2) * mp.gamma(1 + (mu + nu - s) / 2) * mp.gamma(1 + (mu - nu - s) / 2))

def K21t(w, k):
    nu = k - 1
    return 4 * mp.pi * (-1) ** (k // 2) * (4 * mp.pi) ** (-w - 0.5) * mp.sqrt(mp.pi / 2) * MJJ(w + 1, -0.5, nu)

def GR(s):
    return mp.pi ** (-s / 2) * mp.gamma(s / 2)

def G3(t, k):
    mus = [1, k - 1, k]
    v = mp.mpf(1)
    for m in mus:
        v *= GR(0.5 - 1j * t + m) / GR(0.5 + 1j * t + m)
    return v

def chi(s):
    return mp.pi ** (s - 0.5) * mp.gamma((1 - s) / 2) / mp.gamma(s / 2)

for k in (2, 4, 6):
    print("k =", k)
    for t in (0.7, 3.1, 9.4, 25.0, 61.0):
        w = -1j * t
        r = chi(2j * t) * K21t(w, k) / G3(t, k)
        print("   t=%5.1f  ratio = %s" % (t, mp.nstr(r, 12)))
