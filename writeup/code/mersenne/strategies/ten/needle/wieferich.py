"""Wieferich's needle: for odd primes p, the Fermat quotient q(p) = (2^(p-1) - 1)/p mod p.  p is a Wieferich prime when
q(p) = 0 (only 1093 and 3511 are known).  Heuristic: q(p)/p behaves like a uniform number in [0, 1), so
#{p <= x : q(p) = 0} ~ sum 1/p ~ ln ln x, and #{p <= x : min(q, p - q)/p < e} ~ 2 e pi(x).
Prediction written before the run: the 10 bins of q(p)/p for p < 10^7 pass a chi-square test (statistic below 21.7,
the 1 % point on 9 degrees of freedom), and the count with min(q, p - q)/p < 10^-3 is within 3 sd of 2 * 10^-3 * pi(x).
usage: wieferich.py X out.txt"""
import sys, math, gmpy2
from sympy import primerange
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
bins = [0] * 10; zero = []; near = 0; n = 0; e = 1e-3; best = []
for p in primerange(3, X + 1):
    q = (int(gmpy2.powmod(2, p - 1, p * p)) - 1) // p
    n += 1
    if q == 0: zero.append(p); bins[0] += 1; continue
    u = q / p; bins[min(9, int(u * 10))] += 1
    d = min(q, p - q) / p
    if d < e: near += 1
    best.append((d, p, q if q <= p - q else q - p))
best.sort()
pr('odd primes p < %d: %d; Wieferich primes: %s' % (X, n, zero))
exp_ = n / 10; chi = sum((b - exp_) ** 2 / exp_ for b in bins)
pr('bins of q(p)/p: %s; chi-square %.2f on 9 degrees of freedom (1 %% point 21.67)' % (bins, chi))
pr('min(q, p - q)/p < 1e-3: %d primes; expected 2e-3 * %d = %.1f, sd %.1f' % (near, n, 2 * e * n, math.sqrt(2 * e * n * (1 - 2 * e))))
s1 = sum(1.0 / p for p in primerange(3, X + 1))
pr('sum of 1/p over odd primes p < %d: %.4f (expected number of Wieferich primes under the heuristic, counting from p = 3)' % (X, s1))
pr('closest non-Wieferich primes (smallest |q|/p, signed q): ' + ', '.join('%d (q = %d, %.2e)' % (p, qq, d) for d, p, qq in best[:6]))
