"""Euclid's orchard: V(N) = #{(a, b) : 1 <= a, b <= N, gcd(a, b) = 1} = 2 sum_{n <= N} phi(n) - 1, against (6/pi^2) N^2.
The error D(N) = V(N) - 6 N^2 / pi^2: its size relative to N (Mertens: O(N log N)), the mean, the extremes and how often
it is positive, for N up to X.   usage: orchard.py X out.txt"""
import sys, math
import numpy as np
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
phi = np.arange(X + 1, dtype=np.int64)
for p in range(2, X + 1):
    if phi[p] == p: phi[p::p] -= phi[p::p] // p
V = 2 * np.cumsum(phi) - 1                     # V[N]; V[1] = 1
N = np.arange(X + 1, dtype=float); D = V - 6 / math.pi ** 2 * N * N
r = D[1:] / N[1:]
pr('X = %d: V(X) = %d, V/X^2 = %.8f (6/pi^2 = %.8f)' % (X, V[X], V[X] / X ** 2, 6 / math.pi ** 2))
for lo, hi in ((10, 10 ** 3), (10 ** 3, 10 ** 5), (10 ** 5, X)):
    s = r[lo - 1:hi]
    pr('N in [%d, %d]: D/N mean %.4f, sd %.4f, min %.3f, max %.3f, share positive %.3f' % (lo, hi, s.mean(), s.std(), s.min(), s.max(), np.mean(s > 0)))
k = int(np.argmax(np.abs(r))) + 1
pr('largest |D/N|: %.3f at N = %d; largest |D|/(N ln N) for N >= 100: %.4f' % (r[k - 1], k, float(np.max(np.abs(D[100:]) / (N[100:] * np.log(N[100:]))))))
