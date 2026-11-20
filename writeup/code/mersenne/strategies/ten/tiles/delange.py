"""Selberg-Delange coefficients for B(x) = #{n <= x : n = a^2 + b^2}.  The Dirichlet series of the indicator is
F(s) = zeta(s)^(1/2) L(s, chi_4)^(1/2) (1 - 2^-s)^(-1/2) P3(2s)^(1/2), with P3(t) = prod_{q = 3 mod 4} (1 - q^-t)^-1.
Since zeta(t)(1 - 2^-t) / L(t, chi_4) = P3(t)^2 / P3(2t), P3(t) = prod_{j >= 0} [zeta(2^j t)(1 - 2^-(2^j t)) / L(2^j t)]^(1/2^(j+1)).
With H(s) = ((s - 1) zeta(s))^(1/2) * (rest) / s:  B(x) = x / sqrt(ln x) * sum_k lambda_k / ln^k x,
lambda_k = H^(k)(1) / (k! Gamma(1/2 - k)).  Checks: lambda_0 = K = 0.76422365358922, lambda_1/lambda_0 = Shanks's C = 0.58194865932.
usage: delange.py out.txt"""
import sys
from mpmath import mp, zeta, dirichlet, gamma, diff, factorial, mpf
mp.dps = 40
L = lambda s: dirichlet(s, [0, 1, 0, -1])
def P3(t):
    r = mpf(1)
    for j in range(0, 12):
        u = 2 ** j * t
        r *= (zeta(u) * (1 - mpf(2) ** (-u)) / L(u)) ** (mpf(1) / 2 ** (j + 1))
    return r
def H(s):
    a = mpf(1) if s == 1 else (s - 1) * zeta(s)
    return (a * L(s) / (1 - mpf(2) ** (-s)) * P3(2 * s)) ** mpf(0.5) / s
out = open(sys.argv[1], 'w')
lam = [diff(H, 1, k) / (factorial(k) * gamma(mpf(0.5) - k)) for k in range(4)]
for k in range(4):
    s = 'lambda_%d = %s   c_%d = lambda_%d / lambda_0 = %s' % (k, mp.nstr(lam[k], 15), k, k, mp.nstr(lam[k] / lam[0], 12)); print(s); out.write(s + '\n')
