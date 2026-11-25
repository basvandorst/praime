"""Exploratory (not registered): the prime side of the one-point law for the Artin L-function of the cubic field of discriminant -23,
read on ordinary primes: a(p) = (number of roots of x^3 - x - 1 mod p) - 1 (2 if p is represented by x^2 + xy + 6y^2),
S(x) = sum_{p <= x} a(p) log p e(-2 sqrt(p/23)); law S(x) ~ (2 sqrt 2/3) c e^{-i pi/4} (x^3/23)^{1/4}, c = 1.853352.
usage: cubic23.py X"""
import numpy as np, math, sys
X = int(float(sys.argv[1]))
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** 0.5) + 1):
    if sv[i]: sv[i * i::i] = False
# primes represented by the principal form x^2 + xy + 6y^2 of discriminant -23: there x^3 - x - 1 has three roots mod p
rep = np.zeros(X + 1, bool); Y = int(math.sqrt(4 * X / 23)) + 1
for y in range(0, Y + 1):
    # x^2 + xy + 6y^2 <= X
    disc = y * y - 4 * (6 * y * y - X)
    if disc < 0: continue
    r = math.isqrt(disc); x = np.arange((-y - r) // 2 - 1, (-y + r) // 2 + 2, dtype=np.int64); n = x * x + x * y + 6 * y * y; n = n[(n >= 0) & (n <= X)]; rep[n] = True
P = np.nonzero(sv)[0]; P = P[P != 23]
# Legendre symbol (-23/p) = (p/23) for odd p; p = 2: -23 = 1 mod 8 -> +1
qr = np.zeros(23, bool); qr[(np.arange(1, 23) ** 2) % 23] = True
split = qr[P % 23]
a = np.where(rep[P], 2, np.where(split, -1, 0)).astype(float)
print('primes to %.0e: %d; three roots %d, no root %d, one root %d (expected 1/6, 1/3, 1/2)' % (X, len(P), (a == 2).sum(), (a == -1).sum(), (a == 0).sum()))
# check on small primes by counting roots directly
for p in (59, 101, 167, 173, 211, 223, 2, 3, 5, 7):
    roots = sum(1 for x in range(p) if (x ** 3 - x - 1) % p == 0); assert roots - 1 == a[np.searchsorted(P, p)], p
p = P.astype(float); lp = np.log(p)
S = np.cumsum(a * lp * np.exp(-4j * math.pi * np.sqrt(p / 23))); Sc = np.cumsum(a * lp * np.exp(-4j * math.pi * 1.37 * np.sqrt(p / 23)))
pred = (2 * math.sqrt(2) / 3) * 1.853352 / 23 ** 0.25 * np.exp(-1j * math.pi / 4)
for x in (10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, X):
    if x > X: continue
    i = np.searchsorted(P, x, 'right') - 1
    print('x = %.0e  S/x^(3/4) = %+.4f %+.4fi   law %+.4f %+.4fi   control frequency |.| = %.4f' % (x, (S[i] / x ** 0.75).real, (S[i] / x ** 0.75).imag, pred.real, pred.imag, abs(Sc[i]) / x ** 0.75))
