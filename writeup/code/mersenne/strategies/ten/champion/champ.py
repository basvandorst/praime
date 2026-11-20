"""Jumping champions: the most frequent gap between consecutive primes up to x.  Counts N(d, x) of gaps d (p, p+d
consecutive, p + d <= x) at x = 10^2 .. 10^8, the champion(s), and the ratios N(30)/N(6) and N(6)/N(2) against the
heuristic N(d) proportional to S(d) exp(-d / ln x), S(d) = prod_{p | d, p > 2} (p - 1)/(p - 2) (Hardy-Littlewood singular
series times the Cramer factor; Odlyzko, Rubinstein and Wolf 1999 predict 30 overtakes 6 near 1.7e35).
Prediction written before the run: 6 is the sole champion at every power of ten from 10^3 to 10^8, and the model ratio
N(30)/N(6) is within 20 % of the count at 10^8.   usage: champ.py X out.txt"""
import sys, math
import numpy as np
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(s): print(s); out.write(s + '\n')
s = np.ones(X + 1, bool); s[:2] = False
for i in range(2, int(X ** .5) + 1):
    if s[i]: s[i * i::i] = False
P = np.nonzero(s)[0]; del s
G = np.diff(P); Q = P[1:]
def S(d):
    r = 1.0; m = d; q = 3
    while m % 2 == 0: m //= 2
    while q <= m:
        if m % q == 0:
            r *= (q - 1) / (q - 2)
            while m % q == 0: m //= q
        q += 2
    return r
pr('    x   champion(s)   N(2)   N(4)   N(6)   N(30)   N(6)/N(2) meas  model   N(30)/N(6) meas  model')
for e in range(2, int(math.log10(X)) + 1):
    x = 10 ** e; g = G[Q <= x]; c = np.bincount(g, minlength=40)
    top = c.max(); ch = [d for d in range(len(c)) if c[d] == top]
    L = math.log(x)
    m62 = S(6) / S(2) * math.exp(-4 / L); m306 = S(30) / S(6) * math.exp(-24 / L)
    pr('%6.0e %12s %7d %7d %7d %7d %10.3f %8.3f %12.3f %8.3f' % (x, ','.join(map(str, ch)), c[2], c[4], c[6], c[30], c[6] / c[2], m62, c[30] / c[6], m306))
