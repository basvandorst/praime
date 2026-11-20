"""Two-square tiles: B(x) = #{n <= x : n = a^2 + b^2}.  n is a sum of two squares iff every prime q = 3 mod 4 divides n
to an even power (Fermat, Euler).  Landau (1908): B(x) ~ K x / sqrt(ln x), K = 0.7642236535892206 (Landau-Ramanujan).
Second order (Shanks 1964): B(x) ~ K x / sqrt(ln x) * (1 + C / ln x), C = 0.5819486593 (literature value).
Ramanujan's version K * integral_2^x dt / sqrt(ln t).
Prediction written before the run: the first-order formula is off by more than 5 % at 10^8, the second-order one
within 0.5 %.   usage: tiles.py X out.txt"""
import sys, math
import numpy as np
from scipy.integrate import quad
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
K, C = 0.7642236535892206, 0.5819486593
par = np.zeros(X + 1, np.uint8)
s = np.ones(int(X ** .5) + 2, bool); s[:2] = False
for i in range(2, int(len(s) ** .5) + 1):
    if s[i]: s[i * i::i] = False
small = np.nonzero(s)[0]
# primes q = 3 mod 4 up to X: segmented sieve to list them
seg = 1 << 22
qs = []
for lo in range(0, X + 1, seg):
    hi = min(X + 1, lo + seg); m = np.ones(hi - lo, bool)
    if lo == 0: m[:2] = False
    for p in small:
        if p * p >= hi: break
        st = max(p * p, (lo + p - 1) // p * p); m[st - lo::p] = False
    pp = np.nonzero(m)[0] + lo
    qs.append(pp[pp % 4 == 3])
qs = np.concatenate(qs)
for q in qs:                      # par[n] = number of q = 3 mod 4 dividing n to an odd power: +1 at q^1, -1 at q^2, +1 at q^3 ...
    q = int(q); qk = q; sg = 1
    while qk <= X:
        if sg == 1: par[qk::qk] += 1
        else: par[qk::qk] -= 1
        qk *= q; sg = -sg
ok = (par == 0); ok[0] = False
cs = np.cumsum(ok, dtype=np.int64)
pr('X = %d: primes q = 3 mod 4 used: %d; B(100) = %d, first: %s' % (X, len(qs), int(cs[100]), list(np.nonzero(ok[:30])[0])))
pr('      x        B(x)   K x/sqrt(ln x)  ratio   K x/sqrt(ln x)(1+C/ln x)  ratio   K int_2^x dt/sqrt(ln t)  ratio')
for e in range(2, int(math.log10(X)) + 1):
    x = 10 ** e; b = int(cs[x]); L = math.log(x)
    f1 = K * x / math.sqrt(L); f2 = f1 * (1 + C / L)
    ri = K * quad(lambda u: math.exp(u) / math.sqrt(u), math.log(2), L, limit=200)[0]
    pr('%8.0e %11d %15.1f %7.4f %22.1f %7.4f %22.1f %7.4f' % (x, b, f1, b / f1, f2, b / f2, ri, b / ri))
