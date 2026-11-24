"""Check of the arithmetic factor B1 of the triple correlation from first principles: the local factor of the ratios
recipe for the swapped term (A' = {-b1}, B' = {-a, b2}, C = {a}, D = {b1, b2}) is computed as a theta-integral
(trapezoid, exponentially convergent), its log-derivative in b2 (B' only) by complex-step-free central differences,
and compared with the closed forms of Conrey-Snaith 2008 (B1) and 2007 (P). usage: check_B1.py -> check_B1_out.txt"""
import numpy as np
M = 4096; th = np.arange(M) / M; et = np.exp(2j * np.pi * th)
def zpt(x, p, sgn):   # z_{p, sgn*theta}(x) = (1 - e(sgn theta) p^{-x})^{-1}
    return 1 / (1 - (et if sgn > 0 else 1 / et) * p ** (-x))
def Zp(A, B, p):
    out = 1.0 + 0j
    for a in A:
        for b in B: out *= (1 - p ** (-1 - a - b))
    return out
def logA(p, a, b1, b2v, b2):
    Ap = [-b1]; Bp = [-a, b2v]; C = [a]; D = [b1, b2]
    Z = Zp(Ap, Bp, p) * Zp(C, D, p) / (Zp(Ap, D, p) * Zp(Bp, C, p))
    integrand = np.ones(M, complex)
    for x in Ap: integrand *= zpt(0.5 + x, p, -1)
    for x in Bp: integrand *= zpt(0.5 + x, p, +1)
    for x in C: integrand /= zpt(0.5 + x, p, -1)
    for x in D: integrand /= zpt(0.5 + x, p, +1)
    return np.log(Z * integrand.mean())
def b1_cs08(p, x, y):
    lp = np.log(p)
    return (1 - p ** -x) * (1 - p ** -x - p ** -y + p ** (-1 - y)) * lp / ((1 - p ** (-1 + x - y)) * (1 - p ** (-1 - y)) * (1 - 2 / p + p ** (-1 - x)) * p ** (2 - x + y))
def b1_cs07(p, x, y):   # the P(x,y) summand of the 2007 paper, without the A-product and with log(1/p) -> sign
    lp = np.log(p)
    return -(1 - p ** -x) * (1 - p ** -x - p ** -y + p ** (-1 - y)) * lp / ((1 / (1 - p ** (1 - x + y))) * (1 - p ** (-1 - y)) * (1 - 2 / p + p ** (-1 - x)) * p ** (2 - x + y))
rng = np.random.default_rng(1)
print('p   a        b1       b2       | -dlogA/db2 (recipe)        | CS2008 B1 summand           | CS2007 P summand/A')
worst = 0
for p in (2, 3, 5, 7, 11, 101):
    for _ in range(3):
        a, b1, b2 = 1j * rng.uniform(-3, 3, 3); h = 1e-6
        num = -(logA(p, a, b1, b2 + h, b2) - logA(p, a, b1, b2 - h, b2)) / (2 * h)
        x, y = a + b1, a + b2
        c8 = b1_cs08(p, x, y); c7 = b1_cs07(p, x, y)
        worst = max(worst, abs(num - c8))
        print('%-3d %+.3fj %+.3fj %+.3fj | %+.10f%+.10fj | %+.10f%+.10fj | %+.6f%+.6fj' % (p, a.imag, b1.imag, b2.imag, num.real, num.imag, c8.real, c8.imag, c7.real, c7.imag))
print('largest |recipe - CS2008| = %.2e' % worst)
