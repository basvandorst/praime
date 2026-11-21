"""Reference check of the Guinand-Weil explicit formula for zeta with phi(u) = cos^2(pi u / 2a) cos(w u) on |u| <= a.

    sum_rho h(gamma) = h(i/2) + h(-i/2) - g(0) log pi + (1/2pi) int h(r) Re psi(1/4 + i r/2) dr
                       - 2 sum_n Lambda(n) n^(-1/2) g(log n)
with g = phi * phi and h(r) = phi_hat(r)^2, phi_hat(r) = int phi(u) e^{iru} du.
Zero side: zeros from mpmath.zetazero, tail estimated by the zero density.
Usage: python weil_ref.py  -> weil_ref.json
"""
import json, math
import mpmath as mp

mp.mp.dps = 20
NZ = 400
ZEROS = [float(mp.zetazero(k).imag) for k in range(1, NZ + 1)]


def B(a, z):
    c = (mp.pi / a) ** 2
    if abs(z) < 1e-12:
        return mp.mpf(a)
    if abs(z * z - c) < 1e-12:
        return mp.mpf(a) / 2
    return c * mp.sin(a * z) / (z * (c - z * z))


def phihat(a, w, z):
    return (B(a, z - w) + B(a, z + w)) / 2


def h(a, w, z):
    return phihat(a, w, z) ** 2


def phi(a, w, u):
    return mp.cos(mp.pi * u / (2 * a)) ** 2 * mp.cos(w * u) if abs(u) <= a else mp.mpf(0)


def g(a, w, x):
    x = abs(x)
    if x >= 2 * a:
        return mp.mpf(0)
    return mp.quad(lambda t: phi(a, w, t) * phi(a, w, x - t), [x - a, (x - a + a) / 2, a])


def mangoldt(n):
    for p in range(2, int(n ** 0.5) + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return math.log(p) if m == 1 else 0.0
    return math.log(n)


def sides(a, w):
    pole = 2 * h(a, w, mp.mpc(0, 0.5)).real
    g0 = g(a, w, 0) * mp.log(mp.pi)
    arch = mp.quad(lambda r: h(a, w, r) * mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * r / 2)),
                   [0, w, w + 10, w + 40, w + 200, mp.inf]) / mp.pi
    N = int(math.exp(2 * a))
    primes = 0
    for n in range(2, N + 1):
        L = mangoldt(n)
        if L:
            primes += L / math.sqrt(n) * g(a, w, math.log(n))
    primes *= 2
    rhs = pole - g0 + arch - primes
    zs = 2 * sum(h(a, w, t) for t in ZEROS)
    T = ZEROS[-1]
    tail = 2 * mp.quad(lambda r: h(a, w, r) * mp.log(r / (2 * mp.pi)) / (2 * mp.pi), [T, T + 100, mp.inf])
    return {"a": a, "w": w, "pole": float(pole), "g0logpi": float(g0), "arch": float(arch), "primes": float(primes),
            "rhs": float(rhs), "zeros": float(zs), "tail_est": float(tail), "diff": float(rhs - zs - tail)}


if __name__ == "__main__":
    out = []
    for a, w in [(1.0, 0.0), (1.5, 0.0), (2.0, 14.13), (2.0, 17.0), (3.0, 21.02), (3.0, 0.0), (3.0, 40.0)]:
        r = sides(a, w)
        out.append(r)
        print(r)
    json.dump({"zeros": ZEROS, "rows": out}, open("weil_ref.json", "w"))
