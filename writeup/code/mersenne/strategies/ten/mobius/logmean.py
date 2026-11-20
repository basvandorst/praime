"""Log-mean of (M(x)/sqrt(x))^2 over x in [A, X]: (1/ln(X/A)) * integral of (M(x)/sqrt x)^2 dx/x, done exactly since M is
constant between integers: sum over n of M(n)^2 * (1/n - 1/(n+1)).  Compared with the sum over the zeta zeros of
1/|rho zeta'(rho)|^2 (Ng 2004 under RH and further hypotheses).  usage: logmean.py X out.txt"""
import sys, math
import numpy as np
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
spf = np.zeros(X + 1, np.int32)
for p in range(2, int(X ** .5) + 1):
    if spf[p] == 0:
        blk = spf[p * p::p]; blk[blk == 0] = p
mu = np.ones(X + 1, np.int8); mu[0] = 0
n = np.arange(X + 1, dtype=np.int64); m = n.copy(); sgn = np.ones(X + 1, np.int8)
# factor by trial with smallest prime factors (vectorised loop over at most ~24 rounds)
p = spf.astype(np.int64); p[p == 0] = m[p == 0]
for _ in range(30):
    act = m > 1
    if not act.any(): break
    q = spf[m].astype(np.int64); q[q == 0] = m[q == 0]
    q = np.where(act, q, 1)
    m2 = m // q
    sq = act & (m2 % q == 0)
    mu[sq] = 0
    sgn[act] *= -1
    m = m2
mu = np.where(mu == 0, 0, sgn).astype(np.int8); mu[0] = 0
M = np.cumsum(mu, dtype=np.int64)
def pr(s): print(s); out.write(s + '\n')
pr('M(10^k), k = 2..%d: %s' % (int(math.log10(X)), ', '.join(str(int(M[10 ** k])) for k in range(2, int(math.log10(X)) + 1))))
w = 1.0 / n[1:X] - 1.0 / (n[1:X] + 1)
for A in (100, 10 ** 4):
    for B in sorted(set([10 ** 5, 10 ** 6, X])):
        if B > X or B <= A: continue
        s = float(np.sum((M[A:B].astype(np.float64) ** 2) * w[A - 1:B - 1]))
        pr('log-mean of (M/sqrt x)^2 on [%d, %d]: %.4f (rms %.3f)' % (A, B, s / math.log(B / A), math.sqrt(s / math.log(B / A))))
