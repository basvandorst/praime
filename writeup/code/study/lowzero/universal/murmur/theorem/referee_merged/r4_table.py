# Referee: table of Section 2 (k = 2, q -> infinity) and Zubrilina's density from the Euler products, own code.
import math
from sympy import primerange
P = list(primerange(2, 10 ** 6))
def u(l): return 1.0 / ((l - 1) ** 2 * (l + 1))
C0 = 1.0
for l in P: C0 *= 1 - u(l)
def nup(r):  # nu'(r)
    v = C0
    for l in P:
        if l > r: break
        if r % l == 0: v *= (l * l / (l * l - 1.0)) / (1 - u(l))
    return v
def e(l, r): return l * l / (l * l - 1.0) if r % l == 0 else 1 - u(l)
def M2(y, S):  # S = small prime set of N = prod(S) * q, q -> infinity
    kap = 1.0; pi_ = 1.0
    for l in S: kap *= l / (l - 1.0); pi_ *= 1 + 1.0 / l
    s = 0.0; r = 1
    while r < 2 * math.sqrt(y):
        c = 1.0
        for l in S: c /= e(l, r)
        s += nup(r) * c * math.sqrt(4 * y - r * r); r += 1
    return 2 * math.pi * pi_ * math.sqrt(y) + (12 / math.pi) * kap * s - 12 * kap * y
al = 2 * math.pi; be = 2 * math.pi; ga = 12.0
for l in P:
    al *= (l ** 4 - 2 * l ** 2 - l + 1) / (l ** 4 - 2 * l ** 2 + l)
    be *= (l ** 3 + l ** 2 - 1) / (l * (l ** 2 + l - 1))
    ga *= l * (l + 1) / (l ** 2 + l - 1)
def nuz(r):
    v = 1.0
    for l in P:
        if l > r: break
        if r % l == 0: v *= 1 + l * l / (l ** 4 - 2 * l ** 2 - l + 1)
    return v
def Mz(y):
    s = 0.0; r = 1
    while r < 2 * math.sqrt(y):
        s += nuz(r) * math.sqrt(4 * y - r * r); r += 1
    return be * math.sqrt(y) + al * s - ga * y
ys = [0.1, 0.2, 0.5, 1, 2]
for name, S in [("q", []), ("2q", [2]), ("3q", [3]), ("6q", [2, 3])]:
    print("%-4s" % name, " ".join("%8.4f" % M2(y, S) for y in ys))
print("zub ", " ".join("%8.4f" % Mz(y) for y in ys))
# her nu(r) = ratio of the l | r and l !| r mixtures
for l in [2, 3, 5, 7]:
    A = (l * (l * l / (l * l - 1.0)) + 1) / (l + 1); B = (l * (1 - u(l)) + 1) / (l + 1)
    print("l=%d: A/B = %.12f, 1 + l^2/(l^4-2l^2-l+1) = %.12f" % (l, A / B, 1 + l * l / (l ** 4 - 2 * l * l - l + 1)))
print("kappa*nu_N(1) for q, 2q, 3q, 6q:", ["%.3f" % (k * C0 * c) for k, c in [(1, 1), (2, 1.5), (1.5, 16 / 15), (3, 1.5 * 16 / 15)]])
