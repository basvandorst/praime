"""Continuous smooth model of the trains far beyond the reach of counting: q(m) = li^{-1}(m) (q = 2 below m = 1),
G(m) = int_0^m q, N_k(x) = the largest m0 with G(m0 + k) - G(m0) <= x, count(x) = sum_k N_k(x), mean = count / x.
k = 1..2000 summed exactly; larger k by the trapezoid rule on a geometric grid of 6000 points up to K(x) = the largest k
with G(k) <= x (N_k is smooth there).  Checked against the discrete smooth count (trains_smooth.txt) at 10^4 to 10^7 and
against the exact prime counts (trains_count.txt) up to 10^9.  Question: where does the mean stop falling?
usage: trains_far.py out.txt"""
import sys, math
import numpy as np
from scipy.special import expi
out = open(sys.argv[1], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
def li(T): return expi(np.log(T)) - expi(math.log(2))
def li_inv(m):
    m = np.asarray(m, float); T = np.maximum(3.0, m * np.log(np.maximum(m, 2)))
    for _ in range(60): T = np.maximum(2.0000001, T - (li(T) - m) * np.log(T))
    return T
def mean(x):
    M = float(li(np.array([x]))[0]) * 1.2 + 10
    m = np.concatenate([np.linspace(0, 1000, 20001), np.geomspace(1000, M, 400000)[1:]])
    q = np.where(m < 1, 2.0, li_inv(np.maximum(m, 1)))
    G = np.concatenate([[0], np.cumsum((q[1:] + q[:-1]) / 2 * np.diff(m))])
    Kmax = float(np.interp(x, G, m))
    def N(k):
        k = np.asarray(k, float); lo = np.zeros_like(k); hi = np.full_like(k, M)
        for _ in range(90):
            mid = (lo + hi) / 2; ok = np.interp(mid + k, m, G) - np.interp(mid, m, G) <= x
            lo = np.where(ok, mid, lo); hi = np.where(ok, hi, mid)
        return lo
    k1 = np.arange(1, min(2000, int(Kmax)) + 1, dtype=float); tot = float(np.sum(N(k1)))
    if Kmax > 2001:
        kg = np.geomspace(2000, Kmax, 6000); Ng = N(kg)
        tot += float(np.sum((Ng[1:] + Ng[:-1]) / 2 * np.diff(kg))) - Ng[0] / 2
    return tot / x
pr('        x   mean (continuous smooth model)   ln 2 - mean')
for e in range(4, 21):
    v = mean(10.0 ** e); pr('   10^%-3d %10.5f %20.5f' % (e, v, math.log(2) - v))
