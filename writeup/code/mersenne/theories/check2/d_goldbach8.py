# Item D: Goldbach pairs by residue mod 8, and S(n) = #{primes s>=3, s^2<n, n-s^2 prime}.
import numpy as np, math
from sympy import primerange, isprime

NMAX = 2_000_000
s = np.ones(NMAX + 1, bool); s[:2] = False
for i in range(2, int(NMAX**0.5) + 1):
    if s[i]: s[i*i::i] = False
ar = np.arange(NMAX + 1)
F = {a: (s & (ar % 8 == a)).astype(np.float64) for a in (1, 3, 5, 7)}
L = 1 << 22
FT = {a: np.fft.rfft(F[a], L) for a in F}
def conv(a, b):
    c = np.fft.irfft(FT[a] * FT[b], L)[:NMAX + 1]
    return np.rint(c).astype(np.int64)
pairs = {}
for (a, b) in [(1, 7), (7, 1), (3, 5), (5, 3), (1, 1), (5, 5), (1, 3), (3, 1), (5, 7), (7, 5), (1, 5), (5, 1), (3, 3), (7, 7)]:
    pairs[(a, b)] = conv(a, b)

n = np.arange(1_000_000, NMAX + 1, 2)
D = np.zeros(len(n), np.int64)
cls = n % 8
P = lambda a, b: pairs[(a, b)][n]
D = np.where(cls == 0, P(1, 7) + P(7, 1) - P(3, 5) - P(5, 3),
    np.where(cls == 2, P(1, 1) - P(5, 5),
    np.where(cls == 4, P(1, 3) + P(3, 1) - P(5, 7) - P(7, 5),
                       P(1, 5) + P(5, 1) - P(3, 3) - P(7, 7))))

# direct check on a few n
PR = list(primerange(3, NMAX))
def direct(nn):
    c = {}
    for p in PR:
        if p >= nn: break
        q = nn - p
        if isprime(q):
            c[(p % 8, q % 8)] = c.get((p % 8, q % 8), 0) + 1
    g = lambda a, b: c.get((a, b), 0)
    r = nn % 8
    if r == 0: return g(1, 7) + g(7, 1) - g(3, 5) - g(5, 3)
    if r == 2: return g(1, 1) - g(5, 5)
    if r == 4: return g(1, 3) + g(3, 1) - g(5, 7) - g(7, 5)
    return g(1, 5) + g(5, 1) - g(3, 3) - g(7, 7)
for nn in (1_000_000, 1_000_002, 1_234_564, 1_999_998, 2_000_000):
    i = (nn - 1_000_000) // 2
    assert direct(nn) == D[i], (nn, direct(nn), D[i])

# S(n)
S = np.zeros(len(n), np.int64)
for p in primerange(3, int(math.isqrt(NMAX)) + 2):
    sq = p * p
    ok = n > sq
    S[ok] += s[n[ok] - sq]
for nn in (1_000_000, 1_500_002, 2_000_000):
    i = (nn - 1_000_000) // 2
    assert S[i] == sum(1 for p in primerange(3, math.isqrt(nn) + 1) if p * p < nn and isprime(nn - p * p))

for r in (0, 2, 4, 6):
    sel = cls == r
    print('n = %d mod 8: count %d, mean Delta = %.2f, mean S = %.2f' % (r, sel.sum(), D[sel].mean(), S[sel].mean()))
