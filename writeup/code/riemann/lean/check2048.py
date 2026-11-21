import math, numpy as np
from sympy import primerange
def sieve(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]: s[p * p::p] = False
    return np.nonzero(s)[0]
P = sieve(100_000_000)
# 2048: the prime sum G(y) = sum_{n<=y} Lambda(n)(y/n - 1) - y ln y + (1+gamma) y stays within (2 + gamma - ln 4pi) sqrt(y) of ln 2pi - T(y)
g2048 = 0.57721566490153286061; lim2048 = 2 + g2048 - math.log(4 * math.pi)
assert abs(lim2048 - 0.0461914179) < 1e-10
def triv2048(y): u = 1 / y; return -0.5 * math.log(1 - u * u) - math.atanh(u) / u + 1 if y > 1.0001 else 1 - math.log(2)
dev1 = (1 + g2048) - math.log(2 * math.pi) + (1 - math.log(2))
assert abs(dev1 - lim2048) < 1e-15                      # at y = 1 the curve is on the edge
P8 = P[P <= 100_000_000]
nn, ww = [P8.astype(np.float64)], [np.log(P8.astype(np.float64))]
for k in range(2, 27):
    pk = P8[P8 <= int(1e8 ** (1 / k)) + 1]; q = pk.astype(np.int64) ** k; m = q <= 100_000_000
    nn.append(q[m].astype(np.float64)); ww.append(np.log(pk[m].astype(np.float64)))
n2048 = np.concatenate(nn); w2048 = np.concatenate(ww); o = np.argsort(n2048); n2048 = n2048[o]; w2048 = w2048[o]
A = np.cumsum((w2048 / n2048).astype(np.longdouble)); S = np.cumsum(w2048.astype(np.longdouble)); N = n2048.astype(np.longdouble)
dev = ((N * A - S - N * np.log(N) + (1 + g2048) * N) - np.longdouble(math.log(2 * math.pi))) / np.sqrt(N)
ys = np.exp(A + g2048); nxt = np.append(N[1:], np.longdouble(1e8)); ins = (ys > N) & (ys < nxt)
devs = ((ys[ins] * A[ins] - S[ins] - ys[ins] * np.log(ys[ins]) + (1 + g2048) * ys[ins]) - np.longdouble(math.log(2 * math.pi))) / np.sqrt(ys[ins])
big, bigs = N >= 100, ys[ins] >= 100
lo = float(min(dev[big].min(), devs[bigs].min())); hi = float(max(dev[big].max(), devs[bigs].max()))
assert abs(lo + 0.022681) < 5e-7 and float(N[big][dev[big].argmin()]) == 319439
assert abs(hi - 0.024694) < 5e-7 and float(N[big][dev[big].argmax()]) == 34186499
small = [(int(n), math.log(p)) for p in primerange(2, 2001) for n in [p ** k for k in range(1, 12)] if n <= 2000]
worst = 0.0
for i in range(1, 100_000):
    y = 1 + i * 0.02
    G = sum(L * (y / n - 1) for n, L in small if n <= y) - y * math.log(y) + (1 + g2048) * y
    worst = max(worst, abs(G - math.log(2 * math.pi) + triv2048(y)) / math.sqrt(y))
assert worst < lim2048
print("2048 ok", lo, hi, worst)
