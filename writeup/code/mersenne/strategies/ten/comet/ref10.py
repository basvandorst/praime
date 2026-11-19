"""Independent reference for prime10.html (numpy counting per n, no Fourier product): for even n <= X the number of
ordered sums r(n); the n with exactly one unordered sum; the fewest unordered sums for n >= X/2; the band ratios for
3, 5, 7, 11, 13 (mean of r / (n / ln^2 n) for p | n against p not | n, n >= X/2); the width: sd of r / HL over even
n >= X/2 after dividing by a least-squares quadratic, HL = 2 C2 S(n) I(n) with I by the trapezoid rule in ln t (4001 points).
usage: ref10.py X out.json"""
import sys, json, math
import numpy as np
from sympy import primerange, isprime
X = int(sys.argv[1])
P = np.array(list(primerange(2, X + 1))); isP = np.zeros(X + 1, bool); isP[P] = True
r = np.zeros(X + 1, np.int64)
for n in range(4, X + 1, 2):
    q = P[P < n]; r[n] = int(isP[n - q].sum())
assert all(r[n] == sum(1 for p in range(2, n) if isprime(p) and isprime(n - p)) for n in (4, 6, 100, 998, X - X % 2))
G = lambda n: (r[n] + int(isP[n // 2])) // 2          # unordered sums p <= q
lo = X // 2
ones = [n for n in range(4, X + 1, 2) if G(n) == 1]
up = [n for n in range(lo + lo % 2, X + 1, 2)]
low = min((G(n), n) for n in up)
S = np.ones(X + 1)
for p in P[1:]: S[p::p] *= (p - 1) / (p - 2)
def I(n):
    v = np.linspace(np.log(2), np.log(n / 2), 4001); t = np.exp(v); f = t / (np.log(t) * np.log(n - t)); return 2 * float(np.sum((f[1:] + f[:-1]) / 2) * (v[1] - v[0]))   # trapezoid in ln t, symmetric halves (the plain trapezoid in t overestimated I by 0.3 %)
C2 = 0.6601618158468696
ratios = []
for p in (3, 5, 7, 11, 13):
    a = [r[n] / (n / math.log(n) ** 2) for n in up if n % p == 0]; b = [r[n] / (n / math.log(n) ** 2) for n in up if n % p]
    ratios.append(float(np.mean(a) / np.mean(b)))
ns = np.array(up, float); v = np.array([r[n] / (2 * C2 * S[n] * I(n)) for n in up])
fit = np.polyval(np.polyfit(ns, v, 2), ns); width = float((v / fit).std())
out = dict(ones=ones, low=[int(low[0]), int(low[1])], ratios=[round(x, 6) for x in ratios], width=round(width, 6), rsum=int(r.sum()))
json.dump(out, open(sys.argv[2], 'w')); print(out)
