"""I_inf: closed form vs direct enumeration of local residue averages; and the MI of the discrete leading-order law."""
import math, itertools
import numpy as np
from sympy import primerange

def Ip_enum(p):
    tot = 0.0; wsum = 0.0
    for a, b in itertools.product(range(p), repeat=2):
        nu3 = len({0, a % p, (a + b) % p}); nua = len({0, a % p}); nub = len({0, b % p})
        w = (1 - nu3 / p) / (1 - 1 / p) ** 3 / p ** 2
        wsum += w
        if w > 0:
            R = (1 - nu3 / p) * (1 - 1 / p) / ((1 - nua / p) * (1 - nub / p))
            tot += w * math.log(R)
    return tot, wsum

def Ip_closed(p):
    fA = (1 - 2 / p) / (1 - 1 / p) ** 3; fB = (1 - 3 / p) / (1 - 1 / p) ** 3
    v = (p - 1) / p ** 2 * fA * math.log((1 - 1 / p) / (1 - 2 / p))
    if p > 3:
        v += (p - 1) * (p - 2) / p ** 2 * fB * math.log((1 - 3 / p) * (1 - 1 / p) / (1 - 2 / p) ** 2)
    return v

print('p=2 enum', Ip_enum(2))
mx = 0
for p in primerange(3, 62):
    e, ws = Ip_enum(p); c = Ip_closed(p); mx = max(mx, abs(e - c)); assert abs(ws - 1) < 1e-12
print('max |enum-closed| p<=61: %.2e' % mx)
P = np.array(list(primerange(3, 10 ** 7)), float)
fA = (1 - 2 / P) / (1 - 1 / P) ** 3; fB = (1 - 3 / P) / (1 - 1 / P) ** 3
rB = np.where(P > 3, (1 - 3 / P) * (1 - 1 / P) / (1 - 2 / P) ** 2, 1.0)
Ip = (P - 1) / P ** 2 * fA * np.log((1 - 1 / P) / (1 - 2 / P)) + (P - 1) * (P - 2) / P ** 2 * fB * np.log(rB)
# tail: I_p ~ c/p^3; estimate c from the last primes and sum p^-3 over p > 1e7 ~ 1/(2 X^2 ln X)
c = float(np.median(Ip[-1000:] * P[-1000:] ** 3)); X = 1e7
tail = c / (2 * X ** 2 * math.log(X))
I = float(math.fsum(Ip)) + tail
print('c (I_p p^3) = %.4f, tail %.2e' % (c, tail))
print('I_inf nats %.7f bits %.7f; p=3 bits %.6f; p=5 bits %.6f' % (I, I / math.log(2), Ip[0] / math.log(2), Ip[1] / math.log(2)))
print('partial up to 7, 100: bits', Ip[:3].sum() / math.log(2), Ip[P < 100].sum() / math.log(2))

# Discrete leading-order law at finite lam: J(g1,g2) ~ S3({0,g1,g1+g2}) e^{-(g1+g2)/lam}, g even.
def mi_discrete(lam, G=700):
    g = np.arange(2, G + 1, 2); g1, g2 = np.meshgrid(g, g, indexing='ij')
    lnS = np.zeros(g1.shape)
    for p in primerange(3, 2 * G + 3):
        h = {0}
        a = g1 % p; b = (g1 + g2) % p
        nu = 1 + (a != 0) + ((b != 0) & (b != a))
        lnS += np.log(np.maximum(1 - nu / p, 1e-300)) - 3 * math.log(1 - 1 / p)
        lnS[(1 - nu / p) <= 0] = -np.inf
    J = np.exp(lnS - (g1 + g2) / lam); J /= J.sum()
    m1 = J.sum(1); m2 = J.sum(0)
    nz = J > 0
    return float(np.sum(J[nz] * np.log(J[nz] / np.outer(m1, m2)[nz]))) / math.log(2)
for lam in (17.8, 21.1, 23.0, 27.6, 34.5, 50, 100, 200):
    G = int(min(3000, 24 * lam)) // 2 * 2
    print('lam %.1f discrete leading-order MI %.4f bits (G=%d)' % (lam, mi_discrete(lam, G), G), flush=True)
