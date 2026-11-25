"""Exploratory (not registered): the prime side of the one-point law (fact 2055) for the CM forms of Q(i) with Hecke character
(z/|z|)^m, m = 4, 8, 12, read directly on the primes p = a^2 + b^2: S(x) = sum_{p <= x} 2 cos(m theta_p) log p e(-sqrt p),
theta_p = arctan(b/a); law S(x) ~ (2/3) c e^{-i pi/4} x^{3/4}, c the Euler product of fact 2057; a control frequency.
usage: gaussian_primes.py X"""
import numpy as np, math, sys
X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 2 * 10 ** 7
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** 0.5) + 1):
    if sv[i]: sv[i * i::i] = False
# Gaussian primes a + bi with a > b > 0, a^2 + b^2 = p prime (p = 1 mod 4): one per prime
A = int(X ** 0.5) + 1; ps, th = [], []
for a in range(1, A):
    b = np.arange(1, min(a, int(math.isqrt(max(X - a * a, 0))) + 1)); n = a * a + b * b; m = sv[n]
    ps.append(n[m]); th.append(np.arctan2(b[m], a))
p = np.concatenate(ps); th = np.concatenate(th); o = np.argsort(p); p = p[o].astype(float); th = th[o]
print('primes 1 mod 4 up to %.0e: %d' % (X, len(p)))
lp = np.log(p)
for m4 in (4, 8, 12):
    lam = 2 * np.cos(m4 * th)                                   # lambda(p) of the CM form with psi(z) = (z/|z|)^m4
    S = np.cumsum(lam * lp * np.exp(-2j * math.pi * np.sqrt(p)))    # N = 4: e(-2 sqrt(p/4)) = e(-sqrt p)
    # control: wrong frequency
    Sc = np.cumsum(lam * lp * np.exp(-2j * math.pi * 1.37 * np.sqrt(p)))
    # constant: prod over primes (1 - (lambda^2 - 1)/q), inert primes q = 3 mod 4: (1 + 1/q)
    q = np.nonzero(sv)[0]; q3 = q[q % 4 == 3].astype(float)
    logc = np.cumsum(np.log1p(-(lam ** 2 - 1) / p)); c1 = math.exp(logc[np.searchsorted(p, 10 ** 6)]) ; c3 = math.exp(np.sum(np.log1p(1 / q3[q3 <= 10 ** 6])))
    print('psi = (z/|z|)^%d:' % m4)
    for x in (10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, X):
        i = np.searchsorted(p, x, 'right') - 1
        print('   x = %.0e   S/x^(3/4) = %+.4f %+.4fi   |.| = %.4f   control frequency: |S_c|/x^(3/4) = %.4f   (|S|/sqrt(x) = %.1f)' % (x, (S[i] / x ** 0.75).real, (S[i] / x ** 0.75).imag, abs(S[i]) / x ** 0.75, abs(Sc[i]) / x ** 0.75, abs(S[i]) / x ** 0.5))
    print('   Euler product to 1e6 without the factor at 2: split part %.4f, inert part %.4f, product %.4f; (2/3) * product = %.4f' % (c1, c3, c1 * c3, 2 / 3 * c1 * c3))
