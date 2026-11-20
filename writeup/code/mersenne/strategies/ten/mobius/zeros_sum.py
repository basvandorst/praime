"""Sum over the first N zeros rho = 1/2 + i gamma of zeta of 1/|rho zeta'(rho)|^2 (both signs of gamma: twice the sum over
gamma > 0).  Under RH, simple zeros and a growth condition, Ng (2004) shows the log-mean of (M(x)/sqrt(x))^2 equals
this sum over all zeros.  usage: zeros_sum.py N out.txt"""
import sys, mpmath
mpmath.mp.dps = 20
N = int(sys.argv[1]); out = open(sys.argv[2], 'w')
tot = mpmath.mpf(0); part = {}
for k in range(1, N + 1):
    rho = mpmath.zetazero(k); d = mpmath.zeta(rho, derivative=1)
    tot += 2 / (abs(rho) ** 2 * abs(d) ** 2)
    if k in (1, 2, 5, 10, 20, 50, 100, 200): part[k] = tot
for k, v in part.items():
    s = 'first %d zeros: %.6f' % (k, float(v)); print(s); out.write(s + '\n')
